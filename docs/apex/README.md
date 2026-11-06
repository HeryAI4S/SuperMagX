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
