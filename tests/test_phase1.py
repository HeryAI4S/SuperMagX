import math

import numpy as np
import pytest
from scipy.integrate import dblquad

from supermagx.geometry import Magnet, RadialGroup, RectTurn, radial_stack, stacked_layer
from supermagx.excitation import PiecewiseLinearCurrent
from supermagx.electromagnetics import AxisymmetricMagnetostatics, loop_mutual
from supermagx.losses import joule_energy, eddy_energy
from supermagx.case import run_case


def test_arbitrary_zone_map_and_exact_ring_resistance():
    turns = (
        RectTurn(0.20, 0.21, -0.01, 0.01, "supply-B"),
        RectTurn(0.24, 0.25, -0.01, 0.01, "supply-A", direction=-1),
        RectTurn(0.27, 0.28, -0.01, 0.01, "supply-B"),
    )
    magnet = Magnet(turns)
    assert magnet.zones == ("supply-B", "supply-A")
    rho = 6e-9
    expected = rho * 2 * math.pi / (0.02 * math.log(0.21 / 0.20))
    assert turns[0].resistance(rho) == pytest.approx(expected, rel=1e-14)
    assert len(stacked_layer(0.2, 0.01, 0.012, [-0.012, 0, 0.012], "x")) == 3


def test_finite_cross_section_axis_field_against_independent_integral():
    turn = RectTurn(0.20, 0.24, -0.04, 0.02, "A")
    model = AxisymmetricMagnetostatics(Magnet((turn,)), rtol=1e-8)
    result = model.field(0.0, 0.0, {"A": 1234.0})
    h = turn.z_upper - turn.z_lower
    norm = h * math.log(turn.r_outer / turn.r_inner)
    expected = dblquad(
        lambda z, r: 1234.0 * (1 / r) / norm
        * (4e-7 * math.pi) * r * r / (2 * (r * r + z * z) ** 1.5),
        turn.r_inner, turn.r_outer,
        lambda _: turn.z_lower, lambda _: turn.z_upper,
        epsabs=1e-12,
    )[0]
    assert result.br == pytest.approx(0.0, abs=1e-14)
    assert result.bz == pytest.approx(expected, rel=2e-8)
    assert result.relative_change < 1e-8


def test_any_number_of_zones_superposes_and_rejects_missing_current():
    magnet = Magnet(tuple(RectTurn(0.2 + i * .02, 0.205 + i * .02,
                                   -.01, .01, f"z{i}") for i in range(5)))
    model = AxisymmetricMagnetostatics(magnet)
    amps = {f"z{i}": (i + 1) * 1000.0 for i in range(5)}
    whole = model.field(0.04, 0.03, amps)
    parts = [model.field(0.04, 0.03, {f"z{i}": amps[f"z{i}"]})
             for i in range(5)]
    assert whole.bz == pytest.approx(sum(p.bz for p in parts), rel=1e-11)
    assert whole.br == pytest.approx(sum(p.br for p in parts), rel=1e-11)
    with pytest.raises(ValueError, match="unknown zone"):
        model.field(0, 0, {"other": 1})


def test_mutual_matrix_is_reciprocal_and_finite_size_converges():
    magnet = Magnet((RectTurn(.20, .21, -.01, .01, "A"),
                     RectTurn(.25, .26, -.01, .01, "B")))
    model = AxisymmetricMagnetostatics(magnet)
    m = model.mutual_matrix(rtol=3e-4)
    assert m.shape == (2, 2)
    assert np.allclose(m, m.T, rtol=0, atol=1e-14)
    assert np.linalg.eigvalsh(m).min() > 0
    thin = loop_mutual(.205, .255, 0.0)
    assert m[0, 1] == pytest.approx(thin, rel=.02)


def test_joule_integrates_each_independent_zone_exactly():
    magnet = Magnet((RectTurn(.2, .21, -.01, .01, "A"),
                     RectTurn(.25, .26, -.01, .01, "B")))
    wave = PiecewiseLinearCurrent(
        times=(0.0, .01, .011, .021),
        currents={"A": (0, 30000, 30000, 0),
                  "B": (0, 20000, 20000, 0)},
    )
    result = joule_energy(magnet, wave, resistivity=6e-9)
    for zone, peak in (("A", 30000), ("B", 20000)):
        resistance = next(t for t in magnet.turns if t.zone == zone).resistance(6e-9)
        assert result.by_zone[zone] == pytest.approx(
            resistance * peak**2 * (.001 + 2 * .01 / 3), rel=1e-13)


def test_eddy_model_includes_tail_and_accepts_unequal_zone_waveforms():
    magnet = Magnet((RectTurn(.2, .208, -.006, .006, "A"),
                     RectTurn(.23, .238, -.006, .006, "B")))
    wave = PiecewiseLinearCurrent(
        times=(0, .003, .004, .008),
        currents={"A": (0, 20000, 20000, 0),
                  "B": (0, 10000, 8000, 0)},
    )
    result = eddy_energy(magnet, wave, resistivity=6e-9, gamma=.7521616277762179)
    assert result.total_j > result.before_tail_j > 0
    assert result.tail_j == pytest.approx(result.total_j - result.before_tail_j)
    assert set(result.by_zone) == {"A", "B"}


def test_magnetostatic_input_validation_and_internal_field_guard():
    with pytest.raises(ValueError):
        RectTurn(.2, .2, -.01, .01, "A")
    model = AxisymmetricMagnetostatics(Magnet((RectTurn(.2, .21, -.01, .01, "A"),)))
    with pytest.raises(ValueError, match="inside a conductor"):
        model.field(.205, 0.0, {"A": 1.0})


def test_parameterized_radial_stack_has_arbitrary_groups_and_turn_counts():
    magnet = radial_stack(.15, (
        RadialGroup("outer-supply", layers=2, turns_per_layer=3,
                    metal_width=.008, metal_height=.012,
                    radial_pitch=.009, axial_pitch=.013, gap_after=.018),
        RadialGroup("inner-supply", layers=1, turns_per_layer=5,
                    metal_width=.010, metal_height=.010,
                    radial_pitch=.011, axial_pitch=.011),
    ))
    assert len(magnet.turns) == 11
    assert magnet.zones == ("outer-supply", "inner-supply")
    assert magnet.turns[0].r_inner == pytest.approx(.1505)
    assert magnet.turns[-1].r_inner == pytest.approx(.1865)
    assert magnet.turns[-1].z_center == pytest.approx(.022)


def test_json_case_connects_geometry_em_and_losses_without_fixed_zones():
    case = {
        "geometry": {"kind": "radial_stack", "r_start_m": .15,
                     "groups": [
                         {"zone": "bus-2", "layers": 1, "turns_per_layer": 2,
                          "metal_width": .008, "metal_height": .012,
                          "radial_pitch": .009, "axial_pitch": .013,
                          "gap_after": .02},
                         {"zone": "bus-1", "layers": 1, "turns_per_layer": 1,
                          "metal_width": .008, "metal_height": .012,
                          "radial_pitch": .009, "axial_pitch": .013},
                     ]},
        "pulse": {"times_s": [0, .01, .011, .021],
                  "currents_a": {"bus-2": [0, 10000, 10000, 0],
                                 "bus-1": [0, 8000, 8000, 0]}},
        "resistivity_ohm_m": 6e-9,
        "probes_m": [[0, 0]],
        "eddy": {"gamma": .7521616277762179},
    }
    result = run_case(case)
    assert result["geometry"]["turns"] == 3
    assert result["geometry"]["zones"] == ["bus-2", "bus-1"]
    assert result["field_probes"][0]["bz_t"] > 0
    assert result["losses"]["joule_j"] > 0
    assert result["losses"]["eddy_j"] > result["losses"]["eddy_tail_j"] > 0


def test_overlapping_metal_turns_are_rejected():
    with pytest.raises(ValueError, match="overlap"):
        Magnet((RectTurn(.20, .21, -.01, .01, "A"),
                RectTurn(.205, .215, -.01, .01, "B")))
