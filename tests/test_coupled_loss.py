"""Independent energy and circuit checks for the resolved loss option."""

import numpy as np
import pytest
from scipy.integrate import quad

from supermagx.excitation import PiecewiseLinearCurrent
from supermagx.geometry import Magnet, RectTurn
from supermagx.losses import compile_coupled_loss, coupled_network_loss
from supermagx.case import run_case


def _one_turn():
    return Magnet((RectTurn(.10, .12, -.01, .01, "a"),))


def _pulse():
    return PiecewiseLinearCurrent((0, .01, .011, .021),
                                  {"a": (0, 1000, 1000, 0)})


def test_one_cell_is_exact_dc_joule_and_has_no_resolved_eddy():
    turn = _one_turn().turns[0]
    result = coupled_network_loss(_one_turn(), _pulse(), resistivity=6e-9,
                                  radial_cells=1, axial_cells=1)
    integral_i2 = 1000**2*(.01/3 + .001 + .01/3)
    assert result.transport_j == pytest.approx(turn.resistance(6e-9)*integral_i2)
    assert result.eddy_j == 0
    assert result.tail_j == 0


def test_two_cells_match_independent_scalar_rl_solution():
    """For one turn and two cells, the projected circuit has exactly one mode."""
    from supermagx.electromagnetics import AxisymmetricMagnetostatics

    rho = 6e-9
    original = _one_turn().turns[0]
    mid = (original.r_inner+original.r_outer)/2
    cells = (RectTurn(original.r_inner, mid, -.01, .01, 0),
             RectTurn(mid, original.r_outer, -.01, .01, 1))
    r1, r2 = (cell.resistance(rho) for cell in cells)
    matrix = AxisymmetricMagnetostatics(Magnet(cells)).mutual_matrix(
        rtol=2e-4, max_order=24)
    a = matrix[0, 0]+matrix[1, 1]-2*matrix[0, 1]
    g1, g2 = r2/(r1+r2), r1/(r1+r2)
    f = ((matrix[0, 0]-matrix[1, 0])*g1
         + (matrix[0, 1]-matrix[1, 1])*g2)/a
    lam = (r1+r2)/a
    # Independent quadrature of the exact scalar state on all three segments.
    times = _pulse().times
    values = _pulse().values("a")
    q = 0.0
    heat = 0.0
    for t0, t1, i0, i1 in zip(times[:-1], times[1:], values[:-1], values[1:]):
        slope = (i1-i0)/(t1-t0)
        steady = -f*slope/lam
        heat += (r1+r2)*quad(
            lambda t: (steady+(q-steady)*np.exp(-lam*t))**2,
            0, t1-t0, epsabs=1e-11)[0]
        q = steady+(q-steady)*np.exp(-lam*(t1-t0))
    tail = (r1+r2)*q*q/(2*lam)
    result = coupled_network_loss(_one_turn(), _pulse(), resistivity=rho,
                                  radial_cells=2, axial_cells=1,
                                  mutual_rtol=2e-4, mutual_max_order=24)
    assert result.eddy_j == pytest.approx(heat+tail, rel=1e-8)
    assert result.tail_j == pytest.approx(tail, rel=1e-8)
    assert result.tail_j > 0


def test_cross_turn_coupling_and_zone_current_constraints():
    turns = Magnet((RectTurn(.10, .12, -.02, -.005, "a"),
                    RectTurn(.10, .12, .005, .02, "b", direction=-1)))
    waveform = PiecewiseLinearCurrent((0, .01, .02),
                                      {"a": (0, 1000, 0),
                                       "b": (0, 400, 0)})
    result = coupled_network_loss(turns, waveform, resistivity=6e-9,
                                  radial_cells=2, axial_cells=1,
                                  mutual_rtol=2e-4, mutual_max_order=24)
    assert result.eddy_j >= result.tail_j >= 0
    assert result.transport_j > 0
    assert result.total_j == pytest.approx(result.transport_j+result.eddy_j)
    assert set(result.transport_by_zone_j) == {"a", "b"}
    assert sum(result.transport_by_zone_j.values()) == pytest.approx(
        result.transport_j)


def test_json_case_reports_mesh_changes_without_claiming_fem_error():
    case = {
        "geometry": {"kind": "explicit_turns", "turns": [
            {"r_inner": .1, "r_outer": .12, "z_lower": -.01,
             "z_upper": .01, "zone": "a"}]},
        "pulse": {"times_s": [0, .01, .011, .021],
                  "currents_a": {"a": [0, 1000, 1000, 0]}},
        "resistivity_ohm_m": 6e-9,
        "eddy": {"model": "coupled_network", "radial_cells": 2,
                 "axial_cells": 2, "mutual_rtol": 2e-4,
                 "refinement_meshes": [[2, 2], [3, 3]]},
    }
    losses = run_case(case)["losses"]
    assert losses["eddy_j"] > losses["eddy_tail_j"] > 0
    assert losses["network_modes"] == 3
    assert "eddy_by_zone_j" not in losses
    assert losses["network_refinement"][1]["relative_eddy_change"] > 0


def test_compiled_operator_reuses_geometry_and_exposes_low_frequency_kernel():
    magnet = _one_turn()
    operator = compile_coupled_loss(magnet, radial_cells=2, axial_cells=1,
                                    mutual_rtol=2e-4, mutual_max_order=24)
    coefficients = operator.low_frequency_coefficients()
    assert coefficients.zones == ("a",)
    assert coefficients.transport_resistance_per_rho[0, 0] == pytest.approx(
        magnet.turns[0].resistance(1.0))
    assert coefficients.eddy_kernel_times_rho[0, 0] > 0
    assert operator.evaluate(_pulse(), resistivity=6e-9).total_j == pytest.approx(
        coupled_network_loss(magnet, _pulse(), resistivity=6e-9,
                             radial_cells=2, axial_cells=1,
                             mutual_rtol=2e-4, mutual_max_order=24).total_j)
    for rho in (3e-9, 1.2e-8):
        reusable = operator.evaluate(_pulse(), resistivity=rho)
        rebuilt = coupled_network_loss(magnet, _pulse(), resistivity=rho,
                                       radial_cells=2, axial_cells=1,
                                       mutual_rtol=2e-4,
                                       mutual_max_order=24)
        assert reusable.eddy_j == pytest.approx(rebuilt.eddy_j, rel=1e-12)
        assert reusable.transport_j == pytest.approx(rebuilt.transport_j,
                                                      rel=1e-12)
    slow = PiecewiseLinearCurrent((0, 1, 1.001, 2.001),
                                  {"a": (0, 1000, 1000, 0)})
    exact = operator.evaluate(slow, resistivity=6e-9).eddy_j
    limit = coefficients.eddy_kernel_times_rho[0, 0]/6e-9*2*1000**2
    assert exact == pytest.approx(limit, rel=.02)


def test_edge_graded_cells_preserve_dc_resistance_and_positive_eddy_heat():
    magnet = _one_turn()
    uniform = compile_coupled_loss(magnet, radial_cells=4, axial_cells=4)
    graded = compile_coupled_loss(magnet, radial_cells=4, axial_cells=4,
                                  cell_distribution="edge_cosine")
    assert graded.resistance_per_rho[0] == pytest.approx(
        uniform.resistance_per_rho[0], rel=1e-13)
    assert graded.evaluate(_pulse(), resistivity=6e-9).eddy_j > 0
    with pytest.raises(ValueError, match="cell_distribution"):
        compile_coupled_loss(magnet, cell_distribution="unknown")
    uniform_current = Magnet((RectTurn(.1, .12, -.01, .01, "a",
                                       source_profile="uniform"),))
    with pytest.raises(ValueError, match="dc_equilibrium"):
        compile_coupled_loss(uniform_current)
