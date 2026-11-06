"""Biot-Savart and load-balance tests using independent numerical quadrature."""
import math
import pytest
from scipy.integrate import dblquad, quad
from supermagx.design import (
    MU0, ideal_annulus_axis_field, ideal_long_solenoid_hoop_resultant,
    infinite_solenoid_active_build, magnetic_pressure, rectangular_winding_fill,
    required_outer_radius, rounded_rectangle_area,
)

def test_axis_field_independent_quadrature():
    a, b, L, J = .13, .198, .64, 3.2e8
    for z0 in (0., .0025, .17, .50):
        numeric = dblquad(
            lambda z, r: MU0*J*r*r/(2*(r*r+(z0-z)**2)**1.5),
            a, b, lambda _: -L/2, lambda _: L/2,
            epsabs=1e-10, epsrel=1e-11,
        )[0]
        assert ideal_annulus_axis_field(a, b, L, J, z0) == pytest.approx(numeric, rel=2e-11)

def test_parity_sign_and_zero():
    args = (.13, .20, .5, 2.4e8)
    assert ideal_annulus_axis_field(*args, .15) == pytest.approx(
        ideal_annulus_axis_field(*args, -.15), rel=1e-13)
    assert ideal_annulus_axis_field(.13, .13, .5, 2.4e8) == 0
    assert ideal_annulus_axis_field(*args) == pytest.approx(
        -ideal_annulus_axis_field(.13, .20, .5, -2.4e8), rel=1e-14)

def test_inverse_and_long_coil_limit():
    eta = rectangular_winding_fill(.004, .008, .0043, .0083, .0005)
    assert eta == pytest.approx((32-(4-math.pi)*.25)/(4.3*8.3), rel=1e-14)
    builds = []
    for L in (.3, .5, .7, 1.0, 5.0):
        b = required_outer_radius(.13, L, 20., 600e6, eta)
        assert ideal_annulus_axis_field(.13, b, L, 600e6*eta) == pytest.approx(20., rel=1e-12)
        builds.append(b-.13)
    assert all(x > y for x, y in zip(builds, builds[1:]))
    limit = infinite_solenoid_active_build(20., 600e6, eta)
    assert builds[-1] > limit
    assert (builds[-1]-limit)/limit < .005

def test_hoop_identity_independent_force_quadrature():
    a, t, B = .13, .05, 20.
    J = B/(MU0*t)
    force = quad(lambda r: r*J*(MU0*J*(a+t-r)), a, a+t, epsrel=1e-12)[0]
    assert ideal_long_solenoid_hoop_resultant(a, t, B) == pytest.approx(force, rel=5e-14)
    assert magnetic_pressure(B) == pytest.approx(159154943.09189534, rel=1e-12)

def test_invalid_parameters():
    with pytest.raises(ValueError):
        required_outer_radius(.13, .5, 20., 600e6, 1.1)
    with pytest.raises(ValueError):
        required_outer_radius(.13, -.5, 20., 600e6)
    with pytest.raises(ValueError):
        rectangular_winding_fill(.004, .008, .003, .0083)
    with pytest.raises(ValueError):
        rounded_rectangle_area(.004, .008, .003)
    with pytest.raises(ValueError):
        infinite_solenoid_active_build(-20., 600e6)
    with pytest.raises(ValueError):
        ideal_annulus_axis_field(.13, .11, .5, 2e8)
