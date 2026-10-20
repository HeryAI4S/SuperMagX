"""Screening-level induced loss with calibrated one-pole diffusion per axis.

The turn-centre field is used for a local thin-rectangle approximation.
The time constants are phenomenological; this is not a resolved eddy-current
Maxwell solve. Energy after the last waveform knot is integrated analytically.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, pi
from typing import Hashable

import numpy as np

from supermagx.electromagnetics import AxisymmetricMagnetostatics
from supermagx.excitation import PiecewiseLinearCurrent
from supermagx.geometry import Magnet

MU0 = 4e-7*pi


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
        decay = exp(-duration/tau)
        delta = y-drive
        integral += (drive*drive*duration
                     + 2*drive*delta*tau*(1-decay)
                     + delta*delta*tau/2*(1-decay*decay))
        y = drive+delta*decay
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
