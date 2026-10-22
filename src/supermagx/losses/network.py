"""Coupled axisymmetric conducting-cell model for pulse heating.

Compile a geometry once, then evaluate many resistivities and waveforms.
The same operator also yields the slow-ramp quadratic loss matrix.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Hashable

import numpy as np
from scipy.linalg import eigh

from supermagx.electromagnetics import AxisymmetricMagnetostatics
from supermagx.excitation import PiecewiseLinearCurrent
from supermagx.geometry import Magnet, RectTurn

from .eddy import _filtered_slope_integrals


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


@dataclass(frozen=True)
class LowFrequencyLossCoefficients:
    """P = rho*I.T*G*I + I_dot.T*K*I_dot/rho in the slow limit.

    Currents are the independent zone currents in ``zones`` order. G is
    diagonal for series turns; K includes all cross-zone induced coupling.
    K is the coefficient of the metallic eddy increment, not magnetic energy.
    """

    zones: tuple[Hashable, ...]
    transport_resistance_per_rho: np.ndarray
    eddy_kernel_times_rho: np.ndarray


@dataclass(frozen=True)
class CoupledLossOperator:
    """Geometry-only modal operator reusable over rho and current histories."""

    magnet: Magnet
    radial_cells: int
    axial_cells: int
    cell_distribution: str
    resistance_per_rho: np.ndarray
    rates_per_rho: np.ndarray
    forcing: np.ndarray

    def low_frequency_coefficients(self) -> LowFrequencyLossCoefficients:
        """Return the discretized slow-limit alpha*rho + beta/rho forms.

        For slow changes ``y_m ≈ -F_m I_dot/(rho*lambda_m^unit)``. Thus
        ``K = F_zone.T diag(1/lambda_unit) F_zone`` is symmetric positive
        semidefinite. The approximation requires every relevant ramp time
        to be long compared with the retained magnetic diffusion times.
        """
        zones = self.magnet.zones
        zone_index = {zone: i for i, zone in enumerate(zones)}
        winding = np.zeros((len(self.magnet.turns), len(zones)))
        gdiag = np.zeros(len(zones))
        for turn_index, turn in enumerate(self.magnet.turns):
            zi = zone_index[turn.zone]
            winding[turn_index, zi] = turn.direction
            gdiag[zi] += self.resistance_per_rho[turn_index]
        fzone = self.forcing@winding
        if len(self.rates_per_rho):
            kernel = fzone.T@(fzone/self.rates_per_rho[:, None])
            kernel = (kernel+kernel.T)/2
        else:
            kernel = np.zeros((len(zones), len(zones)))
        return LowFrequencyLossCoefficients(zones, np.diag(gdiag), kernel)

    def evaluate(self, waveform: PiecewiseLinearCurrent, *,
                 resistivity: float) -> CoupledLossResult:
        """Evaluate exact modal response to a piecewise-linear zone pulse.

        The initial internal currents are assumed at DC equilibrium. The
        terminal transport heat is integrated only up to the final knot;
        induced currents decay to infinity after it.
        """
        if not np.isfinite(resistivity) or resistivity <= 0:
            raise ValueError("resistivity must be finite and positive")
        unknown = set(waveform.currents)-set(self.magnet.zones)
        if unknown:
            raise ValueError(f"unknown zone: {unknown}")
        durations = np.diff(np.asarray(waveform.times))
        currents = np.array([
            turn.direction*np.asarray(waveform.values(turn.zone))
            for turn in self.magnet.turns])
        squared_integrals = np.sum(
            durations*(currents[:, :-1]**2
                       + currents[:, :-1]*currents[:, 1:]
                       + currents[:, 1:]**2)/3, axis=1)
        by_zone = {zone: 0.0 for zone in self.magnet.zones}
        resistance_by_zone = {zone: 0.0 for zone in self.magnet.zones}
        for t, turn in enumerate(self.magnet.turns):
            resistance = resistivity*self.resistance_per_rho[t]
            by_zone[turn.zone] += float(resistance*squared_integrals[t])
            resistance_by_zone[turn.zone] += float(resistance)
        transport = sum(by_zone.values())
        if not len(self.rates_per_rho):
            return CoupledLossResult(
                transport, 0.0, 0.0, 0.0, transport, by_zone,
                resistance_by_zone, len(self.magnet.turns), 0, None, None)

        rates = resistivity*self.rates_per_rho
        modal_slopes = self.forcing@(np.diff(currents, axis=1)/durations)
        taus = 1/rates
        integrals, remainders = _filtered_slope_integrals(
            -modal_slopes*taus[:, None], durations, taus)
        before = float(np.dot(rates, integrals))
        tail = float(np.dot(rates, remainders))
        if before < 0 or tail < 0:
            raise ArithmeticError("negative modal heat from time integration")
        induced = before+tail
        return CoupledLossResult(
            transport, induced, before, tail, transport+induced, by_zone,
            resistance_by_zone, len(self.magnet.turns)*self.radial_cells
            *self.axial_cells, len(rates), float(1/rates[-1]),
            float(1/rates[0]))


def _cells(magnet: Magnet, nr: int, nz: int,
           distribution: str) -> tuple[Magnet, np.ndarray]:
    cells = []
    parent = []
    radial_nodes = np.linspace(0.0, 1.0, nr+1)
    axial_nodes = np.linspace(0.0, 1.0, nz+1)
    if distribution == "edge_cosine":
        radial_nodes = (1-np.cos(np.pi*radial_nodes))/2
        axial_nodes = (1-np.cos(np.pi*axial_nodes))/2
    for turn_index, turn in enumerate(magnet.turns):
        radii = turn.r_inner+turn.width*radial_nodes
        heights = turn.z_lower+turn.height*axial_nodes
        for ir in range(nr):
            for iz in range(nz):
                cells.append(RectTurn(float(radii[ir]), float(radii[ir+1]),
                                      float(heights[iz]), float(heights[iz+1]),
                                      len(cells), source_profile="dc_equilibrium"))
                parent.append(turn_index)
    return Magnet(tuple(cells)), np.array(parent, dtype=int)


def compile_coupled_loss(
    magnet: Magnet,
    *,
    radial_cells: int = 2,
    axial_cells: int = 2,
    cell_distribution: str = "uniform",
    mutual_rtol: float = 1e-4,
    mutual_max_order: int = 24,
) -> CoupledLossOperator:
    """Compile finite-cell mutuals and diffusion modes once per geometry.

    With cell currents i, L i' + R i = C.T v and C i = I. Set i=G I+Q q,
    where G is the resistive equilibrium splitter and C Q=0. Projection
    eliminates v: (Q.T L Q) q' + (Q.T R Q) q = -Q.T L G I'. R is assembled
    at unit resistivity, so its generalized eigenvalues scale linearly with
    rho while its eigenvectors and forcing do not. All turn and zone
    directions remain arbitrary.
    """
    for name, value in (("radial_cells", radial_cells),
                        ("axial_cells", axial_cells)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    if not (0 < mutual_rtol < 1) or mutual_max_order < 8:
        raise ValueError("invalid mutual quadrature settings")
    if cell_distribution not in ("uniform", "edge_cosine"):
        raise ValueError("cell_distribution must be uniform or edge_cosine")
    if any(turn.source_profile != "dc_equilibrium" for turn in magnet.turns):
        raise ValueError("coupled loss requires dc_equilibrium turn profiles")
    cell_magnet, parent = _cells(magnet, radial_cells, axial_cells,
                                 cell_distribution)
    resistances = np.array([c.resistance(1.0) for c in cell_magnet.turns])
    nturn = len(magnet.turns)
    ncell = len(cell_magnet.turns)
    weights = np.zeros(ncell)
    effective_resistance = np.empty(nturn)
    for t in range(nturn):
        at = np.flatnonzero(parent == t)
        conductance = 1/resistances[at]
        weights[at] = conductance/conductance.sum()
        effective_resistance[t] = 1/conductance.sum()
    if ncell == nturn:
        return CoupledLossOperator(magnet, radial_cells, axial_cells,
                                   cell_distribution,
                                   effective_resistance, np.empty(0),
                                   np.empty((0, nturn)))

    inductance = AxisymmetricMagnetostatics(cell_magnet).mutual_matrix(
        rtol=mutual_rtol, max_order=mutual_max_order)
    # Each Q column has +1 on a free cell and -1 on its parent's final cell.
    last = np.array([np.flatnonzero(parent == t)[-1] for t in range(nturn)])
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
    forcing = eigenvectors.T@coupling
    return CoupledLossOperator(magnet, radial_cells, axial_cells,
                               cell_distribution,
                               effective_resistance, rates, forcing)


def coupled_network_loss(
    magnet: Magnet,
    waveform: PiecewiseLinearCurrent,
    *,
    resistivity: float,
    radial_cells: int = 2,
    axial_cells: int = 2,
    cell_distribution: str = "uniform",
    mutual_rtol: float = 1e-4,
    mutual_max_order: int = 24,
) -> CoupledLossResult:
    """Convenience wrapper; reuse ``compile_coupled_loss`` for pulse scans."""
    operator = compile_coupled_loss(
        magnet, radial_cells=radial_cells, axial_cells=axial_cells,
        cell_distribution=cell_distribution,
        mutual_rtol=mutual_rtol, mutual_max_order=mutual_max_order)
    return operator.evaluate(waveform, resistivity=resistivity)
