"""Geometry and winding-to-power-zone mapping. All dimensions are SI metres."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, log, pi
from typing import Hashable, Iterable


@dataclass(frozen=True)
class RectTurn:
    """One circular conductor with a rectangular r-z metal cross-section.

    ``dc_equilibrium`` gives J_phi proportional to 1/r. ``uniform`` prescribes
    constant J_phi as in a uniform external-current-density FEM source.
    A negative direction reverses its magnetic polarity, not its resistance.
    """

    r_inner: float
    r_outer: float
    z_lower: float
    z_upper: float
    zone: Hashable
    direction: int = 1
    source_profile: str = "dc_equilibrium"

    def __post_init__(self) -> None:
        values = (self.r_inner, self.r_outer, self.z_lower, self.z_upper)
        if not all(isfinite(x) for x in values):
            raise ValueError("turn coordinates must be finite")
        if self.r_inner <= 0 or self.r_outer <= self.r_inner:
            raise ValueError("require 0 < r_inner < r_outer")
        if self.z_upper <= self.z_lower:
            raise ValueError("require z_lower < z_upper")
        if self.direction not in (-1, 1):
            raise ValueError("direction must be -1 or +1")
        if self.zone is None:
            raise ValueError("zone cannot be None")
        if self.source_profile not in ("dc_equilibrium", "uniform"):
            raise ValueError("source_profile must be dc_equilibrium or uniform")

    @property
    def width(self) -> float:
        return self.r_outer - self.r_inner

    @property
    def height(self) -> float:
        return self.z_upper - self.z_lower

    @property
    def r_center(self) -> float:
        return (self.r_inner + self.r_outer) / 2

    @property
    def z_center(self) -> float:
        return (self.z_lower + self.z_upper) / 2

    @property
    def metal_volume(self) -> float:
        return pi * (self.r_outer**2 - self.r_inner**2) * self.height

    @property
    def h_integral(self) -> float:
        """Integral of 1/r over the metal cross-section, in metres."""
        return self.height * log(self.r_outer / self.r_inner)

    def resistance(self, resistivity: float) -> float:
        """Joule-equivalent ring resistance for the selected source profile."""
        if not isfinite(resistivity) or resistivity <= 0:
            raise ValueError("resistivity must be finite and positive")
        if self.source_profile == "uniform":
            return 2*pi*resistivity*self.r_center/(self.width*self.height)
        return 2*pi*resistivity/self.h_integral


@dataclass(frozen=True)
class Magnet:
    turns: tuple[RectTurn, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "turns", tuple(self.turns))
        if not self.turns:
            raise ValueError("magnet needs at least one turn")
        for i, a in enumerate(self.turns):
            for j in range(i+1, len(self.turns)):
                b = self.turns[j]
                if (max(a.r_inner, b.r_inner) < min(a.r_outer, b.r_outer)
                        and max(a.z_lower, b.z_lower) < min(a.z_upper, b.z_upper)):
                    raise ValueError(f"metal turns {i} and {j} overlap")

    @property
    def zones(self) -> tuple[Hashable, ...]:
        return tuple(dict.fromkeys(turn.zone for turn in self.turns))

    def resistances(self, resistivity: float) -> dict[Hashable, float]:
        out = {zone: 0.0 for zone in self.zones}
        for turn in self.turns:
            out[turn.zone] += turn.resistance(resistivity)
        return out


def stacked_layer(
    r_inner: float,
    width: float,
    height: float,
    axial_centers: Iterable[float],
    zone: Hashable,
    direction: int = 1,
    source_profile: str = "dc_equilibrium",
) -> tuple[RectTurn, ...]:
    """Generate one radial layer; each axial center represents a separate turn."""
    return tuple(
        RectTurn(r_inner, r_inner + width, z - height / 2, z + height / 2,
                 zone, direction, source_profile)
        for z in axial_centers
    )


@dataclass(frozen=True)
class RadialGroup:
    """One power-zone group of consecutive radial layers.

    ``gap_after`` is the total radial space to the next group; it can contain
    insulation, reinforcement, or a cooling channel. Those materials do not
    carry transport current in the phase-one electromagnetic model.
    """

    zone: Hashable
    layers: int
    turns_per_layer: int
    metal_width: float
    metal_height: float
    radial_pitch: float
    axial_pitch: float
    gap_after: float = 0.0
    direction: int = 1
    source_profile: str = "dc_equilibrium"

    def __post_init__(self) -> None:
        if self.layers < 1 or self.turns_per_layer < 1:
            raise ValueError("layers and turns_per_layer must be positive")
        lengths = (self.metal_width, self.metal_height, self.radial_pitch,
                   self.axial_pitch, self.gap_after)
        if not all(isfinite(x) for x in lengths):
            raise ValueError("group dimensions must be finite")
        if self.metal_width <= 0 or self.metal_height <= 0:
            raise ValueError("metal dimensions must be positive")
        if self.radial_pitch < self.metal_width or self.axial_pitch < self.metal_height:
            raise ValueError("pitches must be at least the metal dimensions")
        if self.gap_after < 0 or self.direction not in (-1, 1):
            raise ValueError("gap must be nonnegative and direction +/-1")
        if self.source_profile not in ("dc_equilibrium", "uniform"):
            raise ValueError("source_profile must be dc_equilibrium or uniform")


def radial_stack(r_start: float,
                 groups: Iterable[RadialGroup]) -> Magnet:
    """Build a magnet with any number and ordering of series-fed groups.

    ``r_start`` is the cell boundary before the first radial conductor layer.
    Units throughout are metres. Reusing a zone identifier connects groups
    in series in the electrical and magnetic model.
    """
    if not isfinite(r_start) or r_start <= 0:
        raise ValueError("r_start must be finite and positive")
    turns = []
    radius = r_start
    for group in groups:
        centers = ((j-(group.turns_per_layer-1)/2)*group.axial_pitch
                   for j in range(group.turns_per_layer))
        z_centers = tuple(centers)
        for layer in range(group.layers):
            metal_inner = (radius + (layer+.5)*group.radial_pitch
                           - group.metal_width/2)
            turns.extend(stacked_layer(metal_inner, group.metal_width,
                                       group.metal_height, z_centers,
                                       group.zone, group.direction,
                                       group.source_profile))
        radius += group.layers*group.radial_pitch + group.gap_after
    return Magnet(tuple(turns))
