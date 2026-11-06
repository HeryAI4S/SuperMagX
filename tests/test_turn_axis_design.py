"""Independent Biot-Savart quadrature and M26-A 480-turn regression."""
import itertools
import math
import pytest
from scipy.integrate import dblquad
from supermagx.geometry import Magnet, RectTurn
from supermagx.design.turn_axis import (
    MU0, exact_turn_axis_coefficient, zone_axis_coefficients,
    turn_current_density_ceiling, zone_current_density_ceilings,
    axis_field_envelope,
)

@pytest.mark.parametrize("profile",["uniform","dc_equilibrium"])
@pytest.mark.parametrize("z0",[0., .02, .3, -.15])
def test_closed_form_against_independent_double_integral(profile,z0):
    t=RectTurn(.13,.141,-.11,.06,"PS1",source_profile=profile)
    a,b,h=t.r_inner,t.r_outer,t.height
    norm=t.width*h if profile=="uniform" else h*math.log(b/a)
    numeric=dblquad(
        lambda z,r: MU0/2*(1. if profile=="uniform" else 1/r)/norm*
        r*r/(r*r+(z0-z)**2)**1.5,
        a,b,lambda _: t.z_lower,lambda _: t.z_upper,
        epsabs=1e-13,epsrel=2e-12,
    )[0]
    assert exact_turn_axis_coefficient(t,z0)==pytest.approx(numeric,rel=2e-11)

def test_signed_zone_map_inner_edge_peak_j():
    a=RectTurn(.13,.134,-.004,.004,"PS1",source_profile="dc_equilibrium")
    b=RectTurn(.15,.154,-.004,.004,"PS2",source_profile="uniform")
    c=RectTurn(.18,.184,-.004,.004,"PS3",direction=-1,source_profile="uniform")
    m=Magnet((a,b,c))
    k=zone_axis_coefficients(m)
    assert k["PS1"]>0 and k["PS2"]>0 and k["PS3"]<0
    j=600e6
    assert turn_current_density_ceiling(b,j)==pytest.approx(j*.004*.008)
    assert turn_current_density_ceiling(a,j)==pytest.approx(
        j*.13*.008*math.log(.134/.13))
    assert turn_current_density_ceiling(a,j)<j*.004*.008
    e=axis_field_envelope(m,.001,j,expected_supply_count=3,
                           supply_current_ceiling_a=25000.)
    assert e.proportional_currents_a["PS3"]==0
    assert e.field_supremum_t==pytest.approx(
        k["PS1"]*e.current_ceilings_a["PS1"]+
        k["PS2"]*e.current_ceilings_a["PS2"])

def test_three_supply_envelope_hits_target():
    m=Magnet(tuple(RectTurn(.13+i*.02,.134+i*.02,-.04,.04,
                            f"PS{i+1}") for i in range(3)))
    e=axis_field_envelope(m,.3,600e6,expected_supply_count=3)
    assert e.required_fraction_of_ceilings<1
    assert sum(e.coefficients_t_per_a[k]*v
               for k,v in e.proportional_currents_a.items())==pytest.approx(.3,rel=2e-14)
    e2=axis_field_envelope(m,e.field_supremum_t*1.1,600e6)
    assert e2.required_fraction_of_ceilings>1
    assert all(v==0 for v in e2.proportional_currents_a.values())

def test_mixed_section_series_rejected():
    m=Magnet((RectTurn(.13,.134,-.004,.004,"PS1"),
              RectTurn(.15,.153,-.004,.004,"PS1"),
              RectTurn(.17,.174,-.004,.004,"PS2")))
    with pytest.raises(ValueError,match="mixed conductor sections"):
        zone_current_density_ceilings(m,600e6,one_section_per_zone=True)
    assert zone_current_density_ceilings(m,600e6)["PS1"]<(
        zone_current_density_ceilings(m,600e6)["PS2"])

def test_m26_a_480_turn_reference():
    radii=((132.,136.3),(150.6,154.9),(171.2,175.5),(187.8,192.1))
    turns=[]
    for mod,rs in enumerate(radii):
        zone=f"PS{min(mod+1,3)}"
        for r,j in itertools.product(rs,range(60)):
            z=-244.85+8.3*j
            turns.append(RectTurn((r-2)*1e-3,(r+2)*1e-3,
                                  (z-4)*1e-3,(z+4)*1e-3,
                                  zone,source_profile="uniform"))
    m=Magnet(tuple(turns))
    assert len(m.turns)==480 and len(m.zones)==3
    k=zone_axis_coefficients(m)
    nominal={"PS1":14600.785369366342,"PS2":17402.973451509286,
             "PS3":23740.635279040147}
    assert sum(k[z]*nominal[z] for z in m.zones)==pytest.approx(
        19.99944488697459,rel=3e-11)
    e=axis_field_envelope(m,20.,600e6,supply_current_ceiling_a=25000.,
                           expected_supply_count=3)
    assert e.field_supremum_t==pytest.approx(19.467738741260266,rel=3e-11)
    assert e.required_fraction_of_ceilings>1
    assert all(x==pytest.approx(19200.,rel=1e-12)
               for x in e.current_ceilings_a.values())

def test_input_guard_and_reversed_turn():
    t=RectTurn(.13,.134,-.004,.004,"z")
    back=RectTurn(.13,.134,-.004,.004,"z",direction=-1)
    assert exact_turn_axis_coefficient(t)==pytest.approx(-exact_turn_axis_coefficient(back))
    with pytest.raises(ValueError):
        exact_turn_axis_coefficient(t,float("nan"))
    with pytest.raises(ValueError):
        turn_current_density_ceiling(t,-1)
    with pytest.raises(ValueError):
        axis_field_envelope(Magnet((t,)),20.,600e6,expected_supply_count=3)
