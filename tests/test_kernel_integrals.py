"""Independent checks for the singular and reduced-dimensional kernels."""

from math import log, pi

import numpy as np
import pytest
from scipy.integrate import dblquad

import supermagx.electromagnetics.axisymmetric as kernels
from supermagx.electromagnetics.axisymmetric import (
    MU0,
    _pair_mutual,
    _self_mutual,
    loop_field,
    loop_mutual,
)
from supermagx.geometry import RectTurn


def test_distinct_nearly_coincident_filaments_remain_finite():
    gap = 1e-10
    mutual = loop_mutual(1.0, 1.0 + gap, 0.0)
    expected = MU0 * (log(8 / gap) - 2)
    assert np.isfinite(mutual)
    assert mutual == pytest.approx(expected, rel=2e-9)
    br, bz = loop_field(1.0, 0.0, 1.0 + gap)
    assert np.isfinite(br) and np.isfinite(bz)
    assert bz == pytest.approx(-MU0 / (2*pi*gap), rel=2e-7)


def test_near_axis_radial_field_uses_maxwell_limit():
    radius, dz, probe_r = 1.0, 0.4, 1e-10
    br, bz = loop_field(radius, dz, probe_r)
    expected_br = (3*MU0*radius**2*dz*probe_r
                   / (4*(radius**2+dz**2)**2.5))
    expected_bz = MU0*radius**2/(2*(radius**2+dz**2)**1.5)
    assert br == pytest.approx(expected_br, rel=2e-8, abs=0)
    assert bz == pytest.approx(expected_bz, rel=2e-12, abs=0)


@pytest.mark.parametrize("ha,hb,dz", [(.012, .012, .0),
                                           (.012, .020, .031),
                                           (.020, .012, -.013)])
def test_axial_convolution_matches_independent_double_integral(ha, hb, dz):
    source_r, target_r = .20, .28
    reference = dblquad(
        lambda zb, za: loop_mutual(source_r, target_r, dz+za-zb)
        / (ha*hb),
        -ha/2, ha/2, lambda _: -hb/2, lambda _: hb/2,
        epsabs=1e-12,
    )[0]
    reduced = kernels._axial_average(source_r, target_r, dz, ha, hb, 20)
    assert reduced == pytest.approx(reference, rel=2e-10)


def test_reduced_pair_integral_is_symmetric_and_handles_adjacent_turns():
    a = RectTurn(.20, .208, -.006, .006, "A", source_profile="uniform")
    b = RectTurn(.20, .208, .006, .018, "B", source_profile="dc_equilibrium")
    ab = _pair_mutual(a, b, 24)
    ba = _pair_mutual(b, a, 24)
    assert np.isfinite(ab) and ab > 0
    assert ab == pytest.approx(ba, rel=2e-6)


def test_self_quadrature_stays_finite_at_high_order():
    turn = RectTurn(.15, .158, -.006, .006, "A")
    assert np.isfinite(_self_mutual(turn, 64))


def test_narrow_conductor_dc_resistance_uses_stable_logarithm():
    turn = RectTurn(1.3, 1.3 + 1e-12, -.01, .01, "A")
    expected = 2*pi*6e-9/(.02*np.log1p(turn.width/turn.r_inner))
    assert turn.resistance(6e-9) == pytest.approx(expected, rel=2e-12)
