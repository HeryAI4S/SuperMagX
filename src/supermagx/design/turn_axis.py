"""Exact on-axis field of finite rectangular turns and supply field feasibility.

For prescribed uniform or 1/r azimuthal current, the radial and axial
Biot-Savart integrals can be evaluated in closed form on the symmetry axis.
No structural, circuit-transient, or rounded-corner fidelity is claimed.
"""
from __future__ import annotations
from dataclasses import dataclass
from math import asinh, isfinite, log1p, pi
from typing import Hashable, Mapping
from supermagx.geometry import Magnet, RectTurn

MU0 = 4e-7*pi

def _uniform_primitive(c: float, a: float, b: float) -> float:
    if c == 0:
        return 0.
    return c*(asinh(b/abs(c))-asinh(a/abs(c)))

def _resistive_primitive(c: float, a: float, b: float) -> float:
    return asinh(c/a)-asinh(c/b)

def exact_turn_axis_coefficient(turn: RectTurn, z_m: float = 0.) -> float:
    """Signed Bz(0,z)/I [T/A], integrated over the FULL rectangular metal turn."""
    if not isfinite(z_m):
        raise ValueError("z_m must be finite")
    a, b, h = turn.r_inner, turn.r_outer, turn.height
    upper, lower = turn.z_upper-z_m, turn.z_lower-z_m
    if turn.source_profile == "uniform":
        kernel = (_uniform_primitive(upper,a,b)-
                  _uniform_primitive(lower,a,b))/(turn.width*h)
    else:
        kernel = (_resistive_primitive(upper,a,b)-
                  _resistive_primitive(lower,a,b))/(h*log1p(turn.width/a))
    return turn.direction*MU0*kernel/2

def zone_axis_coefficients(magnet: Magnet, z_m: float = 0.) -> dict[Hashable,float]:
    """Signed on-axis field coefficients [T/A] for all series-connected zones."""
    result = {z: 0. for z in magnet.zones}
    for t in magnet.turns:
        result[t.zone] += exact_turn_axis_coefficient(t,z_m)
    return result

def turn_current_density_ceiling(turn: RectTurn, j_limit_a_m2: float) -> float:
    """Peak local J cap converted to transport current [A], for RECTANGULAR metal.

    Uniform source: J=I/(width*height).
    DC equilibrium source: J(r)=I/[r*height*ln(b/a)]; max is at inner radius.
    """
    if not isfinite(j_limit_a_m2) or j_limit_a_m2 <= 0:
        raise ValueError("j_limit_a_m2 must be finite and positive")
    if turn.source_profile == "uniform":
        return j_limit_a_m2*turn.width*turn.height
    return j_limit_a_m2*turn.r_inner*turn.height*log1p(turn.width/turn.r_inner)

def zone_current_density_ceilings(magnet: Magnet, j_limit_a_m2: float,
                                  *, one_section_per_zone: bool=False
                                  ) -> dict[Hashable,float]:
    """The MIN turn current cap in each series circuit; optional section check."""
    ceilings = {z: float("inf") for z in magnet.zones}
    sections = {}
    for t in magnet.turns:
        section = (t.width,t.height)
        if one_section_per_zone and t.zone in sections:
            former=sections[t.zone]
            if any(abs(x-y)>1e-12*max(1.,abs(x),abs(y))
                   for x,y in zip(former,section)):
                raise ValueError(f"mixed conductor sections in zone {t.zone!r}")
        sections.setdefault(t.zone,section)
        ceilings[t.zone]=min(ceilings[t.zone],
                             turn_current_density_ceiling(t,j_limit_a_m2))
    return ceilings

@dataclass(frozen=True)
class AxisFieldEnvelope:
    """Positive-field supremum for nonnegative series currents, fixed geometry."""
    coefficients_t_per_a: Mapping[Hashable,float]
    current_ceilings_a: Mapping[Hashable,float]
    field_supremum_t: float
    target_field_t: float
    required_fraction_of_ceilings: float
    proportional_currents_a: Mapping[Hashable,float]

def axis_field_envelope(
    magnet: Magnet, target_field_t: float, j_limit_a_m2: float, *,
    z_m: float=0., supply_current_ceiling_a: float|None=None,
    expected_supply_count: int|None=None, one_section_per_zone: bool=True,
) -> AxisFieldEnvelope:
    """Best positive Bz with 0<=I_zone<=cap; NOT a structural optimum.

    Supply strict current ceilings yield a SUPREMUM, not necessarily an
    attainable maximum at equality. Negative field coefficients are not
    energized when seeking maximum positive field.
    """
    if not isfinite(target_field_t) or target_field_t < 0:
        raise ValueError("target_field_t must be finite and nonnegative")
    if expected_supply_count is not None and len(magnet.zones)!=expected_supply_count:
        raise ValueError("incorrect number of independent supply zones")
    if supply_current_ceiling_a is not None and (
        not isfinite(supply_current_ceiling_a) or supply_current_ceiling_a<=0
    ):
        raise ValueError("supply_current_ceiling_a must be positive")
    k=zone_axis_coefficients(magnet,z_m)
    caps=zone_current_density_ceilings(magnet,j_limit_a_m2,
                                       one_section_per_zone=one_section_per_zone)
    if supply_current_ceiling_a is not None:
        caps={z:min(c,supply_current_ceiling_a) for z,c in caps.items()}
    positive={z:(caps[z] if k[z]>0 else 0.) for z in magnet.zones}
    bmax=sum(k[z]*positive[z] for z in magnet.zones)
    ratio=(target_field_t/bmax if bmax>0
           else (0. if target_field_t==0 else float("inf")))
    currents={z:(positive[z]*ratio if ratio<=1 else 0.)
              for z in magnet.zones}
    return AxisFieldEnvelope(k,caps,bmax,target_field_t,ratio,currents)
