# M1 validation: existing evidence and remaining gates

**Status:** plan and gap register, not completed whole-APEX validation. Accuracy targets are provisional project gates, not certified uncertainty bounds.

| Test | Reference / invariant | Present evidence | Pending on active APEX case |
|---|---|---|---|
| Single loop / finite solenoid axial field | Biot–Savart analytical field | Existing finite-rectangle field tests | APEX-specific geometry and axis/target profile check |
| Coaxial loop mutual / finite-section self | Neumann elliptic kernel; self-singularity integrable | Existing adaptive quadrature and Duffy treatment, tested | Current conductor cross-section, independent convergence |
| Arbitrary-zone superposition | `B(I)=Σ I_j b_j` | Covered by existing unit tests | New zone map + polarity verification |
| Reciprocity / positivity | `Lij=Lji`, `L>=0` | Existing small-winding unit tests | Active full-zone matrix, conditioning and convergence |
| Magnetic field-energy | `W_B=IᵀLI/2` | Energy computed in existing model | Independent FEM volume-energy and finite-domain convergence |
| Source-profile consistency | `dc_equilibrium` versus `uniform` | Both implemented; transport resistance matched to source | Same source/geometry for all M1/M2/FEM comparisons |
| Force density | `fr=Jphi Bz`, `fz=-Jphi Br` | Load formulas specified | Full spatial force-field extraction and mesh convergence |
| Axial force virtual work | `Fz=I1 I2 dM/dh` | Analytic path available | FEM and finite-difference comparison on matched geometry |
| Radial expansion virtual work | `Fq=0.5 Iᵀ(dL/dq)I` | Analytic path available | Explicit displacement mode and volumetric force comparison |
| Air-domain/truncation | Refine `B,L,W,f` independently | Existing COMSOL transient small-ring domain sensitivity | M1 full-winding stationary study |

## Suggested initial numerical acceptance targets

- Center field grid-refinement change <0.1%.
- Full-zone inductance reciprocity normalized difference <0.1%, energy agreement <0.5% (with outer-domain boundary audit).
- Compare force and virtual-work outcomes with stated numerical tolerance only after mesh and finite-difference-step convergence; a preliminary 1% goal is **not** an uncertainty certificate.
- Do not use componentwise relative error for physically near-zero mutuals/forces; use an appropriate scale.

## Reproducibility and COMSOL conventions

Each result must store design revision/hash, exact metal shape and source profile, zone order/polarity, boundary/air settings, mesh/quadrature, `B` probe coordinates, `L` units, energy and force definitions, test tolerances and code/FEM version. Use the existing package's `python -m pytest -q` and existing [accuracy notes](../accuracy.md). The current independent COMSOL 6.4 reference cases in [loss_theory_validation.md](../loss_theory_validation.md) establish limited **transient loss** checks, **not** a full APEX M1 static force/FEM certification.

**Volume integration:** default COMSOL GUI axisymmetric Surface Integration with *Compute volume integral* applies the (2πr) factor. The repository's COMSOL Java `IntSurface` benchmark does **not**, per its explicit area audit; use `2*pi*r*expr` in that API. Only export `fr,fz` on actual Cu–Nb metal to M4 and keep output as N/m³, not pre-multiplied by (2πr).

**Exit criterion:** do not mark M1 `validated_for_declared_scope` until the active case has consistent source/geometry, multiple independent checks and convergence of field, inductance, energy and distributed force.
