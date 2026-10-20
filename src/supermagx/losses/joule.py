"""Exact integration of I²R for piecewise-linear zone currents and fixed rho."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Hashable

from supermagx.excitation import PiecewiseLinearCurrent
from supermagx.geometry import Magnet


@dataclass(frozen=True)
class JouleResult:
    total_j: float
    by_zone: dict[Hashable, float]
    resistances_ohm: dict[Hashable, float]


def joule_energy(magnet: Magnet, waveform: PiecewiseLinearCurrent,
                 *, resistivity: float) -> JouleResult:
    """Energy in joules; rho is held constant for the entire pulse.

    Zone turns are in series. Independently driven zones may have unrelated
    current time histories. Unknown waveform zones are rejected.
    """
    unknown = set(waveform.currents)-set(magnet.zones)
    if unknown:
        raise ValueError(f"unknown zone: {unknown}")
    resistances = magnet.resistances(resistivity)
    result = {}
    for zone, resistance in resistances.items():
        amps = waveform.values(zone)
        integral = sum(
            (t1-t0)*(i0*i0 + i0*i1 + i1*i1)/3
            for t0, t1, i0, i1 in zip(waveform.times,
                                      waveform.times[1:], amps, amps[1:])
        )
        result[zone] = resistance*integral
    return JouleResult(sum(result.values()), result, resistances)
