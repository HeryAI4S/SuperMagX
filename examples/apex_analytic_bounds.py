"""Run: python examples/apex_analytic_bounds.py (after pip install -e .)"""
from supermagx.design import (
    ideal_annulus_axis_field, ideal_long_solenoid_hoop_resultant,
    infinite_solenoid_active_build, rectangular_winding_fill,
    required_outer_radius,
)
eta = rectangular_winding_fill(.004, .008, .0043, .0083, .0005)
j = 600e6
a = .13
print(f"Metal cell fill fraction: {eta:.6f}")
print(f"Infinite-coil active build: {infinite_solenoid_active_build(20., j, eta)*1e3:.3f} mm")
for L in (.3, .5, .7, 1.):
    b = required_outer_radius(a, L, 20., j, eta)
    at_center = ideal_annulus_axis_field(a, b, L, j*eta)
    at_end = ideal_annulus_axis_field(a, b, L, j*eta, .0025)
    print(f"L={L:.2f} m: active build={(b-a)*1e3:.2f} mm; B0={at_center:.6f} T; B(+2.5mm)={at_end:.6f} T")
for t in (.03, .05, .10):
    N = ideal_long_solenoid_hoop_resultant(a, t, 20.)
    print(f"LONG coil active build={t*1e3:.0f} mm: integrated hoop demand={N/1e6:.3f} MN/m; nominal all-material mean={N/t/1e6:.3f} MPa (NOT Cu-Nb peak)")
