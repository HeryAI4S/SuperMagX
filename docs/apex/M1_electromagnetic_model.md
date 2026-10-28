# M1 — APEX electromagnetic model, V1.0 (integrated with SuperMagX)

**M1 status:** governing equations and interface defined; existing Python free-space static solver already implemented. **No claim** that the current APEX complete winding has passed all M1/FEM tests. See [active constraints](design_constraints.md).

## 1. Current implementation, not a new solver

- [\`src/supermagx/electromagnetics/axisymmetric.py\`](../../src/supermagx/electromagnetics/axisymmetric.py): free-space circular-loop Green functions integrated over finite rectangular metal turns; adaptive field and mutual quadrature; near-axis and near-coincident stability; Duffy treatment of finite-section self terms.
- [\`src/supermagx/geometry.py\`](../../src/supermagx/geometry.py): physical turns, radial groups, polarity and zone identifiers; explicit turns for irregular windings.
- [\`docs/first_principles.md\`](../first_principles.md): actual \`1/r\` and uniform profile definitions, analytic reductions and mutual/energy conventions.
- [\`docs/accuracy.md\`](../accuracy.md): implemented tests, numerical convergence versus external validity.
- [\`docs/loss_theory_validation.md\`](../loss_theory_validation.md): existing constant-resistivity single/two-ring transient COMSOL benchmarks, **not whole-magnet static M1 force validation**.

The Green-kernel solver is an integral-equation realization of the same magnetostatics as the axisymmetric magnetic-vector-potential PDE. It is **not** a separate generic FEM PDE solver.

## 2. Assumptions and source normalization

Use \`(r,phi,z)\`, \`d/dphi=0\`, prescribed closed azimuthal currents and \`mu=mu0\`. For this class of current and boundary conditions a gauge with \`A=A_phi(r,z) e_phi\` is available (not entailed by axial symmetry alone).

Each turn occupies metal cross-section \(\Omega_k\), belongs to zone \(z(k)\), and has polarity \(s_k=\pm1\). With \(\int_{\Omega_k}p_k\,dr\,dz=1\),

\[
J_\varphi(r,z)=\sum_{j=1}^{N_s} I_j g_j(r,z), \qquad
g_j=\sum_{k:z(k)=j}s_k\, p_k(r,z)\mathbf1_{\Omega_k}(r,z).
\]

**Source profiles must match current code**:
- \`uniform\`: \(p_k=1/S_k\), \(S_k=\int_{\Omega_k}dr\,dz\), as for a prescribed FEM uniform external current density.
- \`dc_equilibrium\` (**repo default**): \(p_k=1/(rH_k)\), \(H_k=\int_{\Omega_k}(1/r)\,dr\,dz\). This is the constant-resistivity steady-state ring distribution, \(J_\varphi\propto1/r\).

These are not interchangeable. The actual transient current distribution is determined by magnetic diffusion/transport constraints (M2); the existing coupled cell loss solver requires the \`dc_equilibrium\` profile.

## 3. Maxwell problem, boundary conditions and flux

\[
\nabla\cdot\boldsymbol B=0,\quad \nabla\times\boldsymbol H=\boldsymbol J,\quad
\boldsymbol B=\mu_0\boldsymbol H,\quad \boldsymbol B=\nabla\times\boldsymbol A.
\]

\[
B_r=-\partial_z A_\varphi,\qquad
B_z=\frac1r\partial_r(rA_\varphi),
\]

\[
\boxed{
\partial_r^2 A_\varphi+\frac1r\partial_rA_\varphi
-\frac{A_\varphi}{r^2}+\partial_z^2 A_\varphi=-\mu_0J_\varphi.
}
\]

For \(\psi=r A_\varphi\), \(B_r=-r^{-1}\partial_z\psi\), \(B_z=r^{-1}\partial_r\psi\) and the circular flux \(\Phi(r,z)=2\pi\psi(r,z)\). Regularity at the axis: \(A_\varphi=O(r)\). For localized current, the ideal outer condition is decaying field at infinity. FEM implementations use finite-air truncation with convergence checking or infinite elements. With no free interface sheet current: \(n\cdot[B]=0\), \(n\times[H]=0\).

## 4. Unit-current response and inductance

For each named zone \(j\), drive 1 A into that zone, zero into other zones, and preserve turn polarity/profile. Write \(a_{\varphi,j},\boldsymbol b_j=(b_{r,j},b_{z,j})\) for the resulting unit-current response.

\[
A_\varphi=\sum_j I_ja_{\varphi,j},\quad
\boldsymbol B(r,z)=\sum_j I_j\boldsymbol b_j(r,z),\quad
B_z(0,0)=\boldsymbol k^\mathsf T\boldsymbol I,\quad k_j=b_{z,j}(0,0).
\]

\[
L_{ij}=2\pi\int_{\Omega_{\rm metal}}r\,g_i a_{\varphi,j}\,dr\,dz
=\frac{2\pi}{\mu_0}\int_{\Omega_{\rm all}}
r(b_{r,i}b_{r,j}+b_{z,i}b_{z,j})\,dr\,dz
\]
for an infinite-domain consistent linear model (finite-domain boundary terms must be handled or shown negligible).

\[
\boxed{\boldsymbol L=\boldsymbol L^\mathsf T,\quad
\boldsymbol L\succeq0,\quad
W_m=\tfrac12\boldsymbol I^\mathsf T\boldsymbol L\boldsymbol I
=\int_{\mathbb R^3}|\boldsymbol B|^2/(2\mu_0)\,dV .}
\]

The existing code computes the finite-section inductance using area-weighted circular-filament mutuals; a zero-section filament self-inductance **diverges** and cannot substitute for finite-turn self terms. Use the existing self-quadrature and compare against field-energy and FEM only at matched current profiles and boundary conditions.

## 5. Lorentz-force load for M4

In actual **Cu–Nb metal only**:

\[
\boxed{ f_r=J_\varphi B_z,\qquad f_z=-J_\varphi B_r }\quad[\mathrm{N/m^3}].
\]

Within zone \(i\), the self term scales as \(I_i^2g_i b_i\); cross-zone terms as \(I_iI_jg_i b_j\). These spatial force-density arrays, not total radial Cartesian force, are exported for mechanical coupling. The net Cartesian transverse force of a complete axisymmetric winding is zero by symmetry even though hoop expansion loads are large.

For a specified mechanical virtual displacement \(q\), under fixed currents:
\[
F_q=\tfrac12\boldsymbol I^\mathsf T
(\partial\boldsymbol L/\partial q)\boldsymbol I.
\]
Check it against the matching \( \int\boldsymbol f\cdot(\partial\boldsymbol u/\partial q)\,dV\). The current/geometry variation, conductor domain, and prescribed displacement pattern must be identical. The derived \(rJ_\varphi B_z\), \(B^2/(2\mu_0)\) are only stress-related **proxies**, never Cu–Nb von Mises stress.

## 6. Explicit limits and cross-module contract

- **M1→M2:** zone order/polarity, \(k_j\), \(\boldsymbol L\), unit fields, static source profile. M2 supplies currents, circuit voltages, resistive/induced redistribution, and dynamic magnetic field. A fixed-\(L\) model alone does not validate millisecond magnetic response.
- **M1→M4:** actual per-metal-domain \(f_r,f_z\) in N/m³, with spatial mapping and units. M4 supplies thermal shrinkage, contact, preload, material tensors and structure; magnetic force alone does **not** give stress.
- **M1→M3/M5:** geometry and field as inputs to subsequent heating/material models; no standalone hotspot, thermal cycle or structural allowable inferred.

**Axisymmetric integration audit:** GUI \`Derived Values > Surface Integration\` with *Compute volume integral* enabled already supplies \(2\pi r\); do not multiply again. The existing **COMSOL Java \`IntSurface\`** benchmark in [loss_theory_validation.md](../loss_theory_validation.md) was verified to integrate only the \(r-z\) cross-section, and explicitly requires \`2*pi*r*mf.Qh\`. Match the API and selection; never conflate the two.

See [M1 verification protocol](M1_verification.md) for conditions under which output may be declared validated.
