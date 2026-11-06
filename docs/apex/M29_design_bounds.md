# APEX M29 — first-principles design bounds (2026-11-06)

**Status:** exact analytical formulas for an ideal uniform-current annulus, not a structure calculation, a globally optimized geometry, or an FEM result. The existing SuperMagX discrete finite-cross-section Green-kernel solver is preserved.

## Current design order

1. Center magnetic field **20 T** in **200 mm** warm clear bore; target axial zone **5 mm**, but no invented uniformity tolerance.
2. Demonstrated mechanical reliability of Cu–Nb, Zylon/epoxy, GFRP and insulation, including axial-force transfer and transverse composite failure. **750 MPa** Cu–Nb is a screening limit, **not** a qualified allowable.
3. **Exactly three** independent supplies, each **I < 25 kA**. Working **|V| ≤ 15 kV** with **8 kV** as lower-voltage comparison. A circuit can drive several mechanically separate series modules, all using the same conductor section.
4. Thermal loss optimization last, but material thermal/cyclic safety remains mandatory.

No imposed outside diameter, winding length, radial-layer or turn-count maximum. Prefer small, domestically manufacturable high-strength Cu–Nb. The **600 A/mm² metal J** limit is a current screening target, not an independently measured failure threshold. Metallic load-bearing shells are excluded.

## Exact finite-annulus axis field

Uniform equivalent azimuthal current density J occupies a ≤ r ≤ b and -h ≤ z ≤ h, with h = L/2. Direct integration of the circular-loop Biot–Savart kernel gives:

**Bz(0,z₀) = (μ₀J/2)[G(z₀+h) − G(z₀−h)]**,

where **G(c) = c[asinh(b/|c|) − asinh(a/|c|)]** for nonzero c, G(0) = 0.

At the center:

**B₀ = μ₀Jh[asinh(b/h) − asinh(a/h)]**.

For fixed positive J, B₀, a and h, invert exactly:

**b = h sinh[asinh(a/h) + B₀/(μ₀Jh)]**.

Under a *uniformly homogenized cell approximation* use **J = η J_metal** and η = A_metal/(pitch_r pitch_z), with A_metal = width×height − (4−π) corner_radius². η does not include additional Zylon shells, axial stagger, nonuniform modules or radial gaps. The inverse formula is a **conditional idealized benchmark, not a rigorous lower bound** for arbitrary graded/discrete windings. The long-coil limit is **b−a = B₀/(μ₀ηJ_metal)**.

### Illustrative calculation, not a recommended APEX winding

For a = 130 mm, Cu–Nb 4×8 mm², R0.5 mm, radial/axial pitch 4.3/8.3 mm and J_metal = 600 A/mm²:

| Axial length | Ideal equivalent active radial build |
|---|---:|
| Infinite | 29.78 mm |
| 300 mm | 42.24 mm |
| 500 mm | 34.57 mm |
| 700 mm | 32.28 mm |
| 1000 mm | 31.02 mm |

Local cell η = **0.890597**; rounded wire metal area = **31.7854 mm²**. These builds exclude the extra radial Zylon and cryostat requirements. Field variation across the 5 mm target zone can be computed directly but has no agreed acceptance threshold yet.

## Hoop loading from radial equilibrium, ideal long coil only

Inside a long uniform-J annulus, Bz(r) = B₀(b−r)/(b−a), Lorentz force density f_r = J Bz. Radial equilibrium gives dσr/dr + (σr−σθ)/r + f_r = 0. For traction-free radial boundaries and no axial shear-gradient contribution:

**Nθ = ∫σθ dr = ∫r J Bz dr = [B₀²/(2μ₀)](a + t/3)**, t=b−a.

At **20 T**, B₀²/(2μ₀) = **159.15 MPa** is a characteristic pressure/energy-density scale, not a universally applied wall pressure.

At fixed a and B₀, **Nθ increases with t**, while its average over the entire active thickness **Nθ/t decreases**. Neither is a Cu–Nb local von Mises stress or Zylon transverse-strength prediction. Radial contact, Zylon prestress, axial force, varying layer lengths, insulator shear and real electrical circuits require additional models. Do not rank any physical winding by this scalar alone.

## Reproducibility

Run **pytest -q tests/test_annulus_design.py** and **python examples/apex_analytic_bounds.py**. The tests independently integrate the Biot–Savart kernel and long-coil Lorentz force, verify the inverse formula, monotonic limit, symmetry, and invalid inputs. These checks support the stated analytical identities, **not** mechanical design certification.

Next gate: compare real multi-supply discrete-turn fields, radial/axial force distributions, conductor/manufacturing constraints and 2D structural load closure *before* launching structural Pareto ranking.
