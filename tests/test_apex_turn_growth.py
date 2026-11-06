"""M31 turn growth: exact-field consistency and genuine design assumptions."""
import math
import runpy
from pathlib import Path
import pytest
from supermagx.geometry import RectTurn
from supermagx.design.turn_axis import exact_turn_axis_coefficient, zone_axis_coefficients
from supermagx.geometry import Magnet

NAMESPACE=runpy.run_path(str(Path(__file__).resolve().parents[1]/
                             "studies/apex_turn_growth/scan.py"))
study_case=NAMESPACE["study_case"]

def test_exact_m26_geometry_regression():
    result=study_case((2,2,2,2),60)
    # For identical supply currents the field coefficient is the sum
    # of the three exact rectangular metal turn coefficients.
    assert result["required_equal_supply_current_kA"]==pytest.approx(19.725,rel=5e-4)
    assert result["turns_total"]==480
    assert result["outer_diameter_mm"]==pytest.approx(420.8,abs=1e-9)

def test_axial_and_radial_turn_increases_improve_field_per_peak_current():
    short=study_case((2,2,2,2),60)
    long=study_case((2,2,2,2),80)
    thick=study_case((3,3,3,3),80)
    assert long["required_equal_supply_current_kA"]<short["required_equal_supply_current_kA"]
    assert thick["required_equal_supply_current_kA"]<long["required_equal_supply_current_kA"]
    assert not short["under_600_A_mm2_area_screen"]
    assert long["under_600_A_mm2_area_screen"]
    assert thick["under_600_A_mm2_area_screen"]

def test_axis_turn_growth_saturation_not_linear_with_turns():
    a=study_case((2,2,2,2),60)
    b=study_case((2,2,2,2),80)
    c=study_case((2,2,2,2),120)
    assert a["required_equal_supply_current_kA"]-b["required_equal_supply_current_kA"] > (
        b["required_equal_supply_current_kA"]-c["required_equal_supply_current_kA"])
    assert c["axial_envelope_mm"]>b["axial_envelope_mm"]

def test_same_coefficient_as_existing_turn_axis_model():
    # 4 modules / 8 layers / one turn per layer: direct zone aggregation checks
    # the actual existing closed-form kernel rather than an independent formula.
    r=.13
    turns=[]
    for module in range(4):
        z=f"PS{min(module+1,3)}"
        for i in range(2):
            a=r+i*.0043
            turns.append(RectTurn(a,a+.004,-.004,.004,z,source_profile="uniform"))
        r+=.0086+((.010,.012,.008)[module] if module<3 else 0)
    coeff=zone_axis_coefficients(Magnet(tuple(turns)))
    computed=study_case((2,2,2,2),1)
    assert computed["required_equal_supply_current_kA"]*1000*sum(coeff.values()) == pytest.approx(20.,rel=1e-13)

def test_invalid_topology():
    with pytest.raises(ValueError):
        study_case((2,2,2),60)
    with pytest.raises(ValueError):
        study_case((2,2,2,0),60)
    with pytest.raises(ValueError):
        study_case((2,2,2,2),0)
