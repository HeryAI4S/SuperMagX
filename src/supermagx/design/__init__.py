"""Transparent ideal-annulus physics checks, not a validated structural solver."""
from .annulus import (
    MU0, ideal_annulus_axis_field, ideal_long_solenoid_hoop_resultant,
    infinite_solenoid_active_build, magnetic_pressure, rectangular_winding_fill,
    required_outer_radius, rounded_rectangle_area,
)
__all__ = [
    "MU0", "ideal_annulus_axis_field", "ideal_long_solenoid_hoop_resultant",
    "infinite_solenoid_active_build", "magnetic_pressure",
    "rectangular_winding_fill", "required_outer_radius", "rounded_rectangle_area",
]
