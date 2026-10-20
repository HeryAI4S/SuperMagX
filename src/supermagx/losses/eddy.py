"""Screening-level induced loss with one-pole or ideal-slab diffusion modes.

The turn-centre field is used for a local thin-rectangle approximation.
The one-pole time constant is calibrated; ideal-slab modes follow the 1D
diffusion equation. Neither is a resolved 2D eddy-current Maxwell solve.
Energy after the last waveform knot is integrated analytically.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import expm1, pi
from typing import Hashable

import numpy as np
from numpy.polynomial.legendre import leggauss

from supermagx.electromagnetics import AxisymmetricMagnetostatics
from supermagx.excitation import PiecewiseLinearCurrent
from supermagx.geometry import Magnet

MU0 = 4e-7*pi
_GAUSS_X, _GAUSS_W = leggauss(4)
_GAUSS_X = (_GAUSS_X+1)/2
_GAUSS_W = _GAUSS_W/2


@dataclass(frozen=True)
class EddyResult:
    total_j: float
    before_tail_j: float
    tail_j: float
    by_zone: dict[Hashable, float]
    model: str = "turn-centre thin-rectangle, calibrated one-pole diffusion"


def _filtered_slope_integral(slopes: np.ndarray, durations: np.ndarray,
                             tau: float) -> tuple[float, float]:
    """Integral of y², tau*y' + y = slope(t); includes the final decay."""
    y = 0.0
    integral = 0.0
    for drive, duration in zip(slopes, durations):
        if tau == 0:
            integral += float(drive*drive*duration)
            y = float(drive)
            continue
        x = duration/tau
        if x < .2:
            # The closed form subtracts O(duration) terms to obtain an
            # O(duration^3) result when starting from zero. Four-point
            # Gauss integrates the smooth exact exponential response without
            # catastrophic cancellation in this short-step regime.
            values = y+(drive-y)*(-np.expm1(-x*_GAUSS_X))
            integral += float(duration*np.dot(_GAUSS_W, values*values))
        else:
            delta = y-drive
            one_minus_decay = -expm1(-x)
            one_minus_decay2 = -expm1(-2*x)
            integral += (drive*drive*duration
                         + 2*drive*delta*tau*one_minus_decay
                         + delta*delta*tau/2*one_minus_decay2)
        y += (drive-y)*(-expm1(-x))
    tail = y*y*tau/2 if tau else 0.0
    return integral, tail


def eddy_energy(magnet: Magnet, waveform: PiecewiseLinearCurrent, *,
                resistivity: float, gamma: float,
                field_order: int = 6) -> EddyResult:
    """Induced heat in metal, including all post-waveform diffusion heat.

    ``gamma`` must be supplied explicitly because its calibration is specific
    to conductor physics. Heat is attributed to the turn's zone even when a
    different zone's changing field induced it. It is not terminal power.
    """
    if not np.isfinite(resistivity) or resistivity <= 0:
        raise ValueError("resistivity must be finite and positive")
    if not np.isfinite(gamma) or gamma < 0:
        raise ValueError("gamma must be finite and nonnegative")
    if field_order < 2 or field_order % 2:
        raise ValueError("field_order must be an even integer >= 2")
    unknown = set(waveform.currents)-set(magnet.zones)
    if unknown:
        raise ValueError(f"unknown zone: {unknown}")
    zones = magnet.zones
    times = np.asarray(waveform.times)
    durations = np.diff(times)
    currents = np.array([waveform.values(z) for z in zones])
    current_slopes = np.diff(currents, axis=1)/durations
    model = AxisymmetricMagnetostatics(magnet)
    by_zone = {zone: 0.0 for zone in zones}
    before_tail = tail_total = 0.0
    for turn in magnet.turns:
        br, bz = model.screening_coefficients(turn.r_center, turn.z_center,
                                               order=field_order)
        local_before = local_tail = 0.0
        for dimension, coefficients in ((turn.width, bz), (turn.height, br)):
            tau = gamma*MU0*dimension**2/(pi*pi*resistivity)
            slopes = coefficients@current_slopes
            integral, tail = _filtered_slope_integral(slopes, durations, tau)
            factor = turn.metal_volume*dimension**2/(12*resistivity)
            local_before += factor*integral
            local_tail += factor*tail
        by_zone[turn.zone] += local_before+local_tail
        before_tail += local_before
        tail_total += local_tail
    return EddyResult(before_tail+tail_total, before_tail, tail_total, by_zone)


def eddy_energy_slab_modes(magnet: Magnet, waveform: PiecewiseLinearCurrent,
                           *, resistivity: float, modes: int = 16,
                           field_order: int = 6) -> EddyResult:
    """Ideal one-dimensional slab diffusion, summed over odd spatial modes.

    The conductor-centre applied field remains a screening approximation.
    Each rectangular dimension is treated as an independent infinite slab:
    tau_n = mu0*d²/(rho*pi²*(2n+1)²), and its heat weight is
    8*V*d²/(rho*pi⁴*(2n+1)⁴). Thus the slow-ramp limit sums exactly to
    V*d²/(12*rho)*integral(dB/dt)², without a fitted time constant.
    """
    if not np.isfinite(resistivity) or resistivity <= 0:
        raise ValueError("resistivity must be finite and positive")
    if isinstance(modes, bool) or not isinstance(modes, int) or modes < 1:
        raise ValueError("modes must be a positive integer")
    if field_order < 2 or field_order % 2:
        raise ValueError("field_order must be an even integer >= 2")
    unknown = set(waveform.currents)-set(magnet.zones)
    if unknown:
        raise ValueError(f"unknown zone: {unknown}")
    zones = magnet.zones
    durations = np.diff(np.asarray(waveform.times))
    currents = np.array([waveform.values(zone) for zone in zones])
    current_slopes = np.diff(currents, axis=1)/durations
    field_model = AxisymmetricMagnetostatics(magnet)
    by_zone = {zone: 0.0 for zone in zones}
    before_tail = tail_total = 0.0
    for turn in magnet.turns:
        br, bz = field_model.screening_coefficients(
            turn.r_center, turn.z_center, order=field_order)
        local_before = local_tail = 0.0
        for dimension, coefficients in ((turn.width, bz), (turn.height, br)):
            slopes = coefficients@current_slopes
            for mode in range(modes):
                odd = 2*mode+1
                tau = MU0*dimension**2/(pi*pi*resistivity*odd*odd)
                factor = (8*turn.metal_volume*dimension**2
                          /(pi**4*resistivity*odd**4))
                integral, tail = _filtered_slope_integral(slopes, durations,
                                                          tau)
                local_before += factor*integral
                local_tail += factor*tail
        by_zone[turn.zone] += local_before+local_tail
        before_tail += local_before
        tail_total += local_tail
    return EddyResult(before_tail+tail_total, before_tail, tail_total, by_zone,
                      model=f"1D slab diffusion, {modes} odd modes per axis; "
                            "turn-centre applied field")
