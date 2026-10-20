"""Independent zone currents, represented by a shared piecewise-linear time grid."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Hashable, Mapping


@dataclass(frozen=True)
class PiecewiseLinearCurrent:
    times: tuple[float, ...]
    currents: Mapping[Hashable, tuple[float, ...]]

    def __post_init__(self) -> None:
        times = tuple(float(t) for t in self.times)
        if len(times) < 2 or not all(isfinite(t) for t in times):
            raise ValueError("need at least two finite time points")
        if any(b <= a for a, b in zip(times, times[1:])):
            raise ValueError("times must increase strictly")
        currents = {zone: tuple(float(x) for x in values)
                    for zone, values in self.currents.items()}
        if any(len(values) != len(times) or not all(isfinite(x) for x in values)
               for values in currents.values()):
            raise ValueError("each zone needs one finite current at every time")
        object.__setattr__(self, "times", times)
        object.__setattr__(self, "currents", currents)

    def values(self, zone: Hashable) -> tuple[float, ...]:
        return self.currents.get(zone, (0.0,) * len(self.times))
