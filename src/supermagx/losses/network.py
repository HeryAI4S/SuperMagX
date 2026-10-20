"""Coupled axisymmetric conducting-cell model for pulse heating."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Hashable

import numpy as np
from scipy.linalg import eigh

from supermagx.electromagnetics import AxisymmetricMagnetostatics
from supermagx.excitation import PiecewiseLinearCurrent
from supermagx.geometry import Magnet, RectTurn

from .eddy import _filtered_slope_integral


@dataclass(frozen=True)
class CoupledLossResult:
    transport_j: float
    eddy_j: float
    eddy_before_tail_j: float
    tail_j: float
    total_j: float
    transport_by_zone_j: dict[Hashable, float]
    resistance_by_zone_ohm: dict[Hashable, float]
    cells: int
    modes: int
    shortest_tau_s: float | None
    longest_tau_s: float | None
    model: str = "coupled axisymmetric azimuthal-cell RL Galerkin"


def _cells(magnet: Magnet, nr: int, nz: int) -> tuple[Magnet, np.ndarray]:
    cells = []
    parent = []
    for turn_index, turn in enumerate(magnet.turns):
        radii = np.linspace(turn.r_inner, turn.r_outer, nr+1)
        heights = np.linspace(turn.z_lower, turn.z_upper, nz+1)
        for ir in range(nr):
            for iz in range(nz):
                cells.append(RectTurn(float(radii[ir]), float(radii[ir+1]),
                                      float(heights[iz]), float(heights[iz+1]),
                                      len(cells), source_profile="dc_equilibrium"))
                parent.append(turn_index)
    return Magnet(tuple(cells)), np.array(parent, dtype=int)


def coupled_network_loss(
    magnet: Magnet,
    waveform: PiecewiseLinearCurrent,
    *,
    resistivity: float,
    radial_cells: int = 2,
    axial_cells: int = 2,
    mutual_rtol: float = 1e-4,
    mutual_max_order: int = 24,
) -> CoupledLossResult:
    """Transport and induced heat, including all post-waveform eddy decay.

    Cell currents obey ``L i' + R i = C.T v, C i = I``. Write ``i=G I+Q q``
    with the resistive equilibrium splitter G and zero-sum basis Q. Then
    ``Q.T R G=0`` and ``(Q.T L Q) q'+(Q.T R Q) q=-Q.T L G I'``. A symmetric
    generalized eigensolve reduces the equation to independent RL modes.
    Each mode's heat is integrated over the exact piecewise-linear waveform.
    The initial internal currents are at DC equilibrium; rho is constant.
    """
    if not np.isfinite(resistivity) or resistivity <= 0:
        raise ValueError("resistivity must be finite and positive")
    for name, value in (("radial_cells", radial_cells),
                        ("axial_cells", axial_cells)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    unknown = set(waveform.currents)-set(magnet.zones)
    if unknown:
        raise ValueError(f"unknown zone: {unknown}")

    cell_magnet, parent = _cells(magnet, radial_cells, axial_cells)
    resistances = np.array([c.resistance(resistivity)
                            for c in cell_magnet.turns])
    nturn = len(magnet.turns)
    ncell = len(cell_magnet.turns)
    weights = np.zeros(ncell)
    effective_resistance = np.empty(nturn)
    for t in range(nturn):
        at = np.flatnonzero(parent == t)
        conductance = 1/resistances[at]
        weights[at] = conductance/conductance.sum()
        effective_resistance[t] = 1/conductance.sum()

    times = np.asarray(waveform.times)
    durations = np.diff(times)
    currents = np.array([turn.direction*np.asarray(waveform.values(turn.zone))
                         for turn in magnet.turns])
    squared_integrals = np.sum(durations*(currents[:, :-1]**2
                         + currents[:, :-1]*currents[:, 1:]
                         + currents[:, 1:]**2)/3, axis=1)
    transport_by_zone = {zone: 0.0 for zone in magnet.zones}
    resistance_by_zone = {zone: 0.0 for zone in magnet.zones}
    for t, turn in enumerate(magnet.turns):
        transport_by_zone[turn.zone] += float(effective_resistance[t]
                                              *squared_integrals[t])
        resistance_by_zone[turn.zone] += float(effective_resistance[t])
    transport = sum(transport_by_zone.values())
    if ncell == nturn:
        return CoupledLossResult(transport, 0.0, 0.0, 0.0, transport,
                                 transport_by_zone, resistance_by_zone,
                                 ncell, 0, None, None)

    inductance = AxisymmetricMagnetostatics(cell_magnet).mutual_matrix(
        rtol=mutual_rtol, max_order=mutual_max_order)
    last = np.array([np.flatnonzero(parent == t)[-1]
                     for t in range(nturn)])
    free = np.flatnonzero(~np.isin(np.arange(ncell), last))
    last_for_free = last[parent[free]]
    lq = (inductance[np.ix_(free, free)]
          - inductance[np.ix_(free, last_for_free)]
          - inductance[np.ix_(last_for_free, free)]
          + inductance[np.ix_(last_for_free, last_for_free)])
    rq = np.diag(resistances[free])
    for t in range(nturn):
        at = np.flatnonzero(parent[free] == t)
        rq[np.ix_(at, at)] += resistances[last[t]]
    lg = np.empty((ncell, nturn))
    for t in range(nturn):
        at = np.flatnonzero(parent == t)
        lg[:, t] = inductance[:, at]@weights[at]
    coupling = lg[free, :] - lg[last_for_free, :]
    try:
        rates, eigenvectors = eigh(rq, lq, check_finite=True)
    except np.linalg.LinAlgError as exc:
        raise ArithmeticError("projected inductance is not positive definite; "
                              "refine mutual quadrature or cell mesh") from exc
    if np.any(~np.isfinite(rates)) or np.any(rates <= 0):
        raise ArithmeticError("nonpositive or nonfinite diffusion rate")
    modal_slopes = (eigenvectors.T@coupling)@(
        np.diff(currents, axis=1)/durations)
    before = tail = 0.0
    for rate, slopes in zip(rates, modal_slopes):
        tau = 1/float(rate)
        integrated, remaining = _filtered_slope_integral(-slopes*tau,
                                                          durations, tau)
        before += float(rate)*integrated
        tail += float(rate)*remaining
    if before < 0 or tail < 0:
        raise ArithmeticError("negative modal heat from time integration")
    induced = before+tail
    return CoupledLossResult(transport, induced, before, tail,
                             transport+induced, transport_by_zone,
                             resistance_by_zone, ncell, len(rates),
                             float(1/rates[-1]), float(1/rates[0]))
