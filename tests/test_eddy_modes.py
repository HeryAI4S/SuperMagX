"""Magnetic-diffusion limits and energy bookkeeping."""

import numpy as np
import pytest

from supermagx.case import run_case
from supermagx.electromagnetics import AxisymmetricMagnetostatics
from supermagx.excitation import PiecewiseLinearCurrent
from supermagx.geometry import Magnet, RectTurn
from supermagx.losses.eddy import _filtered_slope_integral


def test_short_ramp_filtered_energy_is_positive_and_has_cubic_limit():
    duration = 1e-8
    before, tail = _filtered_slope_integral(
        np.array([1.0]), np.array([duration]), 1.0)
    assert before == pytest.approx(duration**3/3, rel=2e-8, abs=0)
    assert tail == pytest.approx(duration**2/2, rel=2e-8, abs=0)


def test_slab_modes_converge_to_quasistatic_joule_formula_for_slow_ramp():
    from supermagx.losses import eddy_energy_slab_modes

    turn = RectTurn(.2, .208, -.006, .006, "A")
    magnet = Magnet((turn,))
    peak, ramp_s, rho = 1000.0, 100.0, 6e-9
    waveform = PiecewiseLinearCurrent(
        (0.0, ramp_s, 2*ramp_s), {"A": (0.0, peak, 0.0)})
    result = eddy_energy_slab_modes(magnet, waveform,
                                    resistivity=rho, modes=32)
    one_mode = eddy_energy_slab_modes(magnet, waveform,
                                      resistivity=rho, modes=1)
    br, bz = AxisymmetricMagnetostatics(magnet).screening_coefficients(
        turn.r_center, turn.z_center)
    field_slope = peak/ramp_s
    quasistatic = (turn.metal_volume/(12*rho)
                   * (turn.width**2*(bz[0]*field_slope)**2
                      + turn.height**2*(br[0]*field_slope)**2)
                   * 2*ramp_s)
    assert result.total_j == pytest.approx(quasistatic, rel=2e-4)
    assert 0 < one_mode.total_j < result.total_j
    assert result.tail_j > 0
    assert result.total_j == pytest.approx(result.before_tail_j+result.tail_j)


def test_json_case_selects_uncalibrated_slab_modes_without_gamma():
    case = {
        "geometry": {"kind": "explicit_turns", "turns": [
            {"r_inner": .2, "r_outer": .208,
             "z_lower": -.006, "z_upper": .006, "zone": "A"}]},
        "pulse": {"times_s": [0, .01, .011, .021],
                  "currents_a": {"A": [0, 10000, 10000, 0]}},
        "resistivity_ohm_m": 6e-9,
        "eddy": {"model": "slab_modes", "modes": 8},
    }
    result = run_case(case)
    assert result["losses"]["eddy_model"].startswith("1D slab diffusion")
    assert result["losses"]["eddy_tail_j"] > 0
