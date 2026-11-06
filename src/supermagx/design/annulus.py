"""First-principles ideal-annulus magnetic benchmarks; all dimensions SI.

Uniform equivalent azimuthal current density in an annular finite solenoid.
This is not a discrete-turn, circuit, contact, or structural solver.
"""
from __future__ import annotations
from math import asinh, isfinite, pi, sinh

MU0 = 4.0e-7 * pi

def _positive(name: str, value: float) -> None:
    if not isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")

def rounded_rectangle_area(width_m: float, height_m: float, corner_radius_m: float = 0.) -> float:
    """Metal cross-section of an externally rounded rectangular wire [m²]."""
    _positive("width_m", width_m)
    _positive("height_m", height_m)
    if not isfinite(corner_radius_m) or not 0 <= corner_radius_m <= min(width_m, height_m)/2:
        raise ValueError("invalid corner radius")
    return width_m*height_m - (4-pi)*corner_radius_m**2

def rectangular_winding_fill(width_m: float, height_m: float, radial_pitch_m: float,
                             axial_pitch_m: float, corner_radius_m: float = 0.) -> float:
    """Local 2D metal filling factor. EXCLUDES extra Zylon radial gaps."""
    _positive("radial_pitch_m", radial_pitch_m)
    _positive("axial_pitch_m", axial_pitch_m)
    if radial_pitch_m < width_m or axial_pitch_m < height_m:
        raise ValueError("pitch smaller than metal dimensions")
    return rounded_rectangle_area(width_m, height_m, corner_radius_m)/(radial_pitch_m*axial_pitch_m)

def ideal_annulus_axis_field(r_inner_m: float, r_outer_m: float, axial_length_m: float,
                             equivalent_j_a_m2: float, z_m: float = 0.) -> float:
    """Exact axial Bz [T] of uniformly filled annulus, a<=r<=b, |z|<=L/2.

    Positive or negative equivalent J is supported. Formula follows direct
    Biot-Savart integration; no FEM or fitted coefficient.
    """
    _positive("r_inner_m", r_inner_m)
    _positive("axial_length_m", axial_length_m)
    if not isfinite(r_outer_m) or r_outer_m < r_inner_m:
        raise ValueError("r_outer_m must be finite and >= r_inner_m")
    if not isfinite(equivalent_j_a_m2) or not isfinite(z_m):
        raise ValueError("J and z must be finite")
    h = axial_length_m/2
    def antiderivative(c: float) -> float:
        if c == 0:
            return 0.
        return c*(asinh(r_outer_m/abs(c))-asinh(r_inner_m/abs(c)))
    return MU0*equivalent_j_a_m2/2 * (antiderivative(z_m+h)-antiderivative(z_m-h))

def required_outer_radius(r_inner_m: float, axial_length_m: float,
                          target_center_field_t: float, max_metal_j_a_m2: float,
                          assumed_metal_fraction: float = 1.) -> float:
    """Exact inverse under SATURATED, CONSTANT homogenized J=eta*J_metal.

    NOT a universal geometry lower bound for arbitrary graded/discrete turns.
    Not a prediction of shell thickness or structural feasibility.
    """
    _positive("r_inner_m", r_inner_m)
    _positive("axial_length_m", axial_length_m)
    _positive("max_metal_j_a_m2", max_metal_j_a_m2)
    if not isfinite(target_center_field_t) or target_center_field_t < 0:
        raise ValueError("target B must be finite and nonnegative")
    if not isfinite(assumed_metal_fraction) or not 0 < assumed_metal_fraction <= 1:
        raise ValueError("metal fraction must be in (0, 1]")
    h = axial_length_m/2
    arg = asinh(r_inner_m/h)+target_center_field_t/(MU0*max_metal_j_a_m2*assumed_metal_fraction*h)
    try:
        b = h*sinh(arg)
    except OverflowError as exc:
        raise OverflowError("radius overflows for specified short coil") from exc
    if not isfinite(b):
        raise OverflowError("outer radius nonfinite")
    return b

def infinite_solenoid_active_build(target_field_t: float, max_metal_j_a_m2: float,
                                   assumed_metal_fraction: float = 1.) -> float:
    """Long-coil limit: B/(mu0 eta Jmetal), active radial build [m]."""
    _positive("max_metal_j_a_m2", max_metal_j_a_m2)
    if not isfinite(target_field_t) or target_field_t < 0:
        raise ValueError("target field must be finite and nonnegative")
    if not isfinite(assumed_metal_fraction) or not 0 < assumed_metal_fraction <= 1:
        raise ValueError("metal fraction must be in (0, 1]")
    return target_field_t/(MU0*max_metal_j_a_m2*assumed_metal_fraction)

def magnetic_pressure(field_t: float) -> float:
    """B²/(2mu0) [Pa], a characteristic energy density, not a wall traction."""
    if not isfinite(field_t):
        raise ValueError("field must be finite")
    return field_t**2/(2*MU0)

def ideal_long_solenoid_hoop_resultant(r_inner_m: float, active_build_m: float,
                                       bore_field_t: float) -> float:
    """Integrated hoop-load demand [N/m] for an ideal LONG annular solenoid.

    With uniform J, free inner/outer radial boundaries, and negligible axial
    shear gradient, integrating radial equilibrium gives:
    int(sigma_theta dr) = int(r J Bz dr) = p_mag*(a+t/3).
    This is NOT the Cu-Nb or Zylon peak stress or failure utilization.
    """
    _positive("r_inner_m", r_inner_m)
    _positive("active_build_m", active_build_m)
    return magnetic_pressure(bore_field_t)*(r_inner_m+active_build_m/3)
