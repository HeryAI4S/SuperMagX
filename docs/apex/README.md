# APEX current model specification

This directory integrates the **current APEX physics-model program** with the existing SuperMagX implementation. Existing algorithms, studies, and FEM comparison records remain unchanged.

- [Design constraints and provenance](design_constraints.md)
- [M1 electromagnetic specification and code mapping](M1_electromagnetic_model.md)
- [M1 verification and pending gaps](M1_verification.md)
- [Constraint snapshot](constraints.v1.json) — reference only, not an active optimization gate.

SuperMagX already implements the finite-cross-section axisymmetric Green-kernel model in [axisymmetric.py](../../src/supermagx/electromagnetics/axisymmetric.py); the [first-principles derivation](../first_principles.md), [accuracy notes](../accuracy.md), and [COMSOL loss comparison](../loss_theory_validation.md) are retained. M1 is an interface and verification specification, **not a new parallel solver** or a claim of validation of the current full APEX geometry.

The historical [OP5 screening study](../../studies/apex_cunb_op5/op5_report.md) used an old <40 kA limit and a 29.785 kA candidate. It remains reproducible but does **not** satisfy the active <25 kA limit.

## Updated research baseline (2026-11-06)

- [Current APEX priority/constraint manifest](constraints.current.json): field → structure → voltage → thermal; three supplies; no compactness/turn-count maximum.
- [M29 exact finite-annulus design identities](M29_design_bounds.md): finite-length on-axis field/inverse radius, metal cell fill, long-coil hoop resultant, verification and explicit limits.
- Python module: [supermagx.design.annulus](../../src/supermagx/design/annulus.py), example: [apex_analytic_bounds.py](../../examples/apex_analytic_bounds.py), [independent quadrature tests](../../tests/test_annulus_design.py).
- The historical [constraints.v1.json](constraints.v1.json) remains unchanged for provenance; **do not** treat its older 30 kV / ≥3 supply scenario as the latest design requirement.

**No new claimed structural optimum:** these exact idealized annulus equations are not a predictor of actual Cu–Nb maximum stress or Zylon transverse failure.


### M30 exact discrete-turn on-axis field feasibility (2026-11-06)

[Closed-form field and 3-supply limit](M30_discrete_field_envelope.md) adds independent analytic `Bz(0,z)` for `uniform` and `dc_equilibrium` source profiles and a local metal current-density ceiling. The historical 480-turn M26-A winding attains 19.99944 T at its original currents but only **19.46774 T** with `J≤600 A/mm²` using ideal 4×8 mm² rectangular conductor. This is a fixed-geometry magnetic infeasibility test, **not** a mechanical stress or rounded-corner result. Reproduce via [example](../../examples/apex_m26_field_envelope.py) and [tests](../../tests/test_turn_axis_design.py).


### M31: Increase turn count without assuming structural safety

[Deterministic axial/radial turn-growth calculation](../../studies/apex_turn_growth/README.md) uses the already-verified finite-turn closed-form field coefficient and three-supply mapping. Increasing 8-layer turns from 480 to 640 reduces the calculated equal-current 20 T load from 19.725 to 18.398 kA, under a rounded-area 600 A/mm² screening threshold (field geometry still rectangular). This does **not** establish a structural or voltage optimum.
