# CuNb pulse-loss derivation and independent COMSOL check

## Scope and conventions

Each insulated physical turn is a solid, axisymmetric conducting ring with
metal cross-section \(A_n\) in the \(r,z\) plane. Turns may belong to any
number of independently powered zones. The imposed signed turn current is
\(I_n=s_n I_{k(n)}\), where \(s_n\in\{-1,+1\}\). This calculation uses one
constant, isotropic resistivity \(\rho\) throughout all CuNb and free-space
permeability \(\mu_0\); the benchmark uses \(\rho=6\times10^{-9}\)
\(\Omega\,\mathrm m\) unless stated otherwise. Displacement current, iron,
conductive supports, leads, strand insulation, and heating-induced changes in
\(\rho\) are excluded. Each \(r,z\) metal cross-section is a **continuous
azimuthal conductor**; changing that connectivity changes the problem.

The independent COMSOL benchmark uses its [single-conductor Coil feature](https://doc.comsol.com/6.4/doc/com.comsol.help.acdc/acdc_ug_magnetic_fields.08.061.html),
which resolves skin and proximity current redistribution under a prescribed
current. It is a magnetic transient, without thermal feedback. COMSOL defines
the electromagnetic loss density as [`mf.Qh`](https://doc.comsol.com/6.4/doc/com.comsol.help.acdc/acdc_ug_modeling.06.23.html).

## Continuum equations and exact heat split

For azimuthal vector potential \(A_\phi\), quasistatic Maxwell and Ohm
equations give

\[
\nabla\times\!\left(\mu_0^{-1}\nabla\times\boldsymbol A\right)
=J_\phi\boldsymbol e_\phi,\qquad
\rho J_\phi=\frac{V_n}{2\pi r}-\partial_t A_\phi
\quad (r,z)\in A_n.
\]

The unknown \(V_n\) enforces \(\int_{A_n}J_\phi\,\mathrm dA=I_n\); it is the
effective voltage around that closed turn. Define

\[
H_n=\int_{A_n}\frac{\mathrm dA}{2\pi r},\quad
J_{0,n}=\frac{I_n}{2\pi rH_n},\quad
R_n=\frac{\rho}{H_n}.
\]

For a rectangle of height \(h_n\) and radii \(r_i,r_o\),
\(H_n=h_n\ln(r_o/r_i)/(2\pi)\). The actual current is
\(J_\phi=J_{0,n}+J_{e,n}\), with \(\int_{A_n}J_{e,n}\,\mathrm dA=0\). Thus,

\[
\int_{A_n}2\pi r\rho J_\phi^2\,\mathrm dA
=R_n I_n^2+\int_{A_n}2\pi r\rho J_{e,n}^2\,\mathrm dA.
\]

The cross term vanishes **exactly**, even with finite skin depth, because
\(2\pi rJ_{0,n}=I_n/H_n\) is constant on a turn. Therefore the CuNb heat has
a precise transport term and a nonnegative induced-current increment. The
all-turn transport power is \(\rho\boldsymbol I^\mathsf TG\boldsymbol I\),
where \(G_{kk}=\sum_{n\in k}1/H_n\). This is a sum of each zone's
\(I_k^2R_k\); it is not the resistance of one series supply across zones.

Eliminating \(V_n\) gives an implicit field-current relation. If
\(\Pi_n f=f-(\int_{A_n}f\,\mathrm dA)/(2\pi rH_n)\), then

\[
J_{e,n}=-\rho^{-1}\Pi_n\partial_t A_\phi.
\]

The projection removes precisely the net-turn-current component. Since
\(A_\phi\) includes the field of \(J_e\) itself, evaluating this expression
requires magnetic diffusion; inserting only the transport field is a
slow-change approximation.

## Derivation of the \(\alpha\rho+\beta/\rho\) limit

Let \(a_k(r,z)\) be the magnetostatic vector potential per ampere from
zone \(k\)'s relaxed transport current distribution. For ramps slow compared
with magnetic diffusion, \(\partial_t A_\phi\simeq\sum_k a_k\dot I_k\).
Then

\[
P_{e,\mathrm{slow}}=\frac{1}{\rho}\dot{\boldsymbol I}^{\mathsf T}
C\dot{\boldsymbol I},\quad
C_{kl}=\sum_n\left[
\int_{A_n}2\pi r a_k a_l\,\mathrm dA
-\frac{(\int_{A_n}a_k\,\mathrm dA)(\int_{A_n}a_l\,\mathrm dA)}{H_n}
\right].
\]

\(C\) is symmetric positive semidefinite because it is an integral of
projected-field squares. This gives the proposed form
\(P\simeq\rho I^\mathsf TGI+\dot I^\mathsf T C\dot I/\rho\), including
off-diagonal coupling between independent supplies. For one common normalized
linear trapezoid, rise \(T_r\), flat \(F\), fall \(T_f\), and zone current
vector \(I_*g(t)\), the full-pulse slow-limit heat is

\[
Q_{\rm slow}=\rho(I_*^\mathsf TGI_*)
\left(F+\frac{T_r+T_f}{3}\right)
+\frac{I_*^\mathsf TC I_*}{\rho}
\left(\frac1{T_r}+\frac1{T_f}\right).
\]

For symmetric ramps \(T_r=T_f=T\), the unique positive stationary point in
this restricted waveform family is

\[
T_*^{\rm slow}=\frac1\rho
\sqrt{\frac{3I_*^\mathsf TC I_*}{I_*^\mathsf TGI_*}}.
\]

The flat-top transport heat is constant in this derivative. At ramp times
comparable to a diffusion time, the \(1/T\) eddy term is **not** accurate;
the full modal response below must be integrated. The finite optimum above
applies to *linear* ramps with free duration and fixed peak current, not to
arbitrarily shaped waveforms.

The slow-limit stationary time is not guaranteed to be the global optimum of
the full diffusion problem. As \(T\to0\) at fixed peak current and fixed
flat-top duration, each modal state undergoes a finite jump rather than an
unbounded \(1/T\) dissipation. Eddy heat approaches a finite limit, while
ramp transport heat approaches zero. Depending on diffusion times and the
flat-top interval, the full heat minimum can therefore occur at the shortest
allowed ramp. A finite optimum must be checked with the modal solution and
any power-supply voltage, mechanical, and cooling constraints; the
\(\alpha\rho+\beta/\rho\) expression is an asymptotic screening model.

## Finite-diffusion Galerkin solver and reusable operator

Partition every turn into azimuthal current cells. Each cell uses the exact
\(1/r\) equilibrium profile within its rectangle. The finite-cross-section
loop Green function supplies the symmetric cell mutual matrix \(L\); the
diagonal \(R\) contains exact cell ring resistances. With \(C_t\) mapping
cell currents \(i\) to signed physical-turn currents \(I_t\), the weak
Maxwell/Ohm equations are

\[
L\dot i+Ri=C_t^\mathsf T v,\qquad C_t i=I_t.
\]

The inverse-resistance splitter \(G_t\) and within-turn zero-sum basis \(Q\)
give \(i=G_tI_t+Qq\), \(C_tQ=0\), and \(Q^\mathsf TRG_t=0\). Eliminating
the voltages yields

\[
L_q\dot q+R_q q=-F\dot I_t,\quad
L_q=Q^\mathsf TLQ,\quad R_q=Q^\mathsf TRQ,\quad
F=Q^\mathsf TLG_t.
\]

The generalized eigensystem \(R_q u_m=\lambda_m L_q u_m\), normalized with
\(u_m^\mathsf TL_q u_l=\delta_{ml}\), reduces this to
\(\dot y_m+\lambda_my_m=-f_m^\mathsf T\dot I_t\). The modal Joule increment
is \(P_e=\sum_m\lambda_m y_m^2\). Each piecewise-linear segment has an
analytic exponential response; after the final knot, the current-redistribution
tail is exactly \(Q_{e,\rm tail}=\tfrac12\sum_my_m(t_{\rm end})^2\).
For a completed pulse, set every turn current to zero at the final knot;
otherwise its transport Joule power continues beyond the reported interval.

For constant \(\rho\), \(R_q=\rho R_q^{(1)}\), so
\(\lambda_m=\rho\lambda_m^{(1)}\). The eigenvectors and forcing matrix depend
only on geometry. `compile_coupled_loss(...)` stores this operator once;
`operator.evaluate(waveform, resistivity=...)` reuses it for pulse and
material scans. Its `low_frequency_coefficients()` returns the discretized
\(G\) and \(C\) matrices: \(C\approx F_z^\mathsf T
\operatorname{diag}(1/\lambda_m^{(1)})F_z\). These are screening coefficients
for slow ramps, while `evaluate` retains all resolved diffusion modes.

`cell_distribution="edge_cosine"` concentrates cells at metal boundaries,
where skin/proximity current gradients can be strongest. It does not change
the physical model or fit to a COMSOL result. Mesh convergence must be checked
on **eddy heat and tail heat**, since a stable total can hide error in either.

## COMSOL comparison protocol

`validation/comsol/CoupledLossBenchmark.java` independently builds a 2D
axisymmetric current-controlled single-conductor model. The source cases and
runner are in `validation/comsol/run_validation.py`; `compare.py` evaluates
the corresponding Python geometry, current waveforms, and resistivity. A
second benchmark uses two isolated rings, two independently timed current
sources, and the same material law. Each COMSOL case is solved with two
metal-mesh/time-output resolutions; a larger air domain is checked separately.

**Axisymmetric integration audit:** The archived benchmark's `IntSurface`
node integrated only the \(r,z\) cross-section. Its output for `expr=1`
was \(0.0004\,\mathrm{m}^2\) for the 20×20 mm ring. The revolved volume
was \(\int_A 2\pi r\,\mathrm dA=0.00027646015\,\mathrm{m}^3\). Under
that **section** convention the physical power is
`IntSurface(2*pi*r*mf.Qh)`; integrating `mf.Qh` without the Jacobian in
that same convention caused a roughly 32% discrepancy in the first check.
The current benchmark explicitly sets `intvolume=on` and integrates
`mf.Qh`. It separately sets `intvolume=off` and integrates
`2*pi*r*mf.Qh`; the two agree at all 241 saved times in the new
`single_a_axisym_audit` run (reported maximum relative difference 0).
Its total heat is 0.08736684170182557 J, versus 0.08736684170181795 J
in the archived `single_a_coarse` run. The [new native-run summary](../validation/comsol/axisym_audit_result.json)
records the case and values; the ignored MPH and trace can be regenerated.
The historical benchmark record and
its source hashes are retained. Never multiply by `2*pi*r` a second time
when revolved-volume integration is enabled. An independent integration
of `rho*mf.Jphi^2` in the same volume convention checks `mf.Qh`; these
expressions describe the same loss and must not be added.

The benchmark record reports transport heat from the exact geometric ring
resistance, and inferred FEM eddy heat as total COMSOL heat minus that
transport value. This difference needs more time/mesh accuracy than the total
heat. A small discrepancy in this controlled continuous-metal benchmark does
not establish a universal bound for a 260-turn winding, variable \(\rho(T,B)\),
conductive support, or a composite CuNb strand topology.

## Numerical results and conditional mesh extrapolation

All energies below are joules for the deliberately small 1 kA or 1 kA/0.7 kA
verification cases, **not** the heat of the 20 T APEX magnet. The first three
are a single ring with inner radius 100 mm; the final case adds a second ring
starting at 150 mm with its pulse delayed 3 ms. Each flat top is 1 ms.
The high-resistivity row is a numerical sensitivity check at
\(1.2\times10^{-8}\ \Omega\,\mathrm m\); it does not change the specified
CuNb design value \(6\times10^{-9}\ \Omega\,\mathrm m\).

| Case; conductor \(w\times h\); ramp | COMSOL total | 24×24 edge-cell total | \(p=2\) extrapolated total | Total error | Eddy error | After-pulse tail error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| A; 20×20 mm; 10 ms | 0.08738472 | 0.08734048 | 0.08738147 | −0.004% | −0.040% | −0.232% |
| Narrow; 8×16 mm; 3 ms | 0.09930483 | 0.09927396 | 0.09931314 | +0.008% | +0.111% | −0.451% |
| Higher \(\rho\); 20×20 mm; 5 ms | 0.09796686 | 0.09793427 | 0.09797795 | +0.011% | +0.133% | −0.959% |
| Two independent rings; 20×20 mm; 10 ms | 0.15472642 | 0.15464108 | 0.15471535 | −0.007% | −0.059% | +0.332% |

The extrapolation is **not fitted to COMSOL**. Edge-graded 16×16, 20×20, and
24×24 cell results had observed convergence orders 1.988–1.990 for eddy
energy and 1.953–1.968 for tail energy. For each of those two quantities
separately, the stated estimate uses
\(Q_\infty^{(2)}=Q_{24}+(Q_{24}-Q_{20})/[(24/20)^2-1]\); exact transport
energy is added to the extrapolated eddy energy. The observed order is a
diagnostic, while exponent 2 was chosen from the smooth-cell discretization,
not adjusted to fit FEM. Do not apply this extrapolation if three successive
meshes do not show stable monotone convergence, especially under stronger
skin effects or for sharply different conductor topology.

The COMSOL reference itself is resolution dependent. Changing both metal
mesh and output interval from fine to ultra changed **total** heat by
0.008–0.053% across these four cases, and **tail** heat by 0.18–0.72%.
For case A, separate checks showed that halving the time-output interval
changed total heat by −0.025%, while halving only metal `hmax` changed it by
less than 0.00002%; increasing the air radius from 1 to 2 m changed it by
about 0.005%. These are measured refinement differences, not rigorous error
bounds. The after-pulse tail is especially sensitive because it is a small
and rapidly varying fraction of the total.

The 24×24 raw cell model is already below 0.55% eddy-energy error in these
four cases, but some raw tail errors exceed 1%. The conditional extrapolation
improves all four measured tail comparisons to within 1%; it gives no
universal 1% certificate. Every new conductor family and supply topology
needs its own mesh-order check, and selected full-magnet cases still need FEM
validation. The dense mutual matrix and generalized eigensystem scale poorly
with cell count, so applying 24×24 cells to every turn of a 260-turn magnet
would not be a resource-saving screening method.

For fixed geometry, the compiled operator avoids rebuilding the mutual
matrix and eigenmodes for every pulse. On this Windows host, the 16×16
single-ring edge-cell operator took 9.11 s to compile; 100 pulse evaluations
took 0.11 s for each tested resistivity. These are measurements of the
included [`compiled_loss_scan.py`](../examples/compiled_loss_scan.py), not
general timing guarantees. The scan's minimum was at its shortest sampled
ramp, illustrating why the slow-limit stationary formula alone cannot set
the pulse duration.

## Reproduce and inspect

From the repository root with SuperMagX installed and COMSOL 6.4 available
at the path configured in `run_validation.py`:

```powershell
python validation/comsol/run_validation.py single_a_coarse single_a_fine single_a_ultra single_a_timeonly single_a_meshonly single_a_farfield narrow_coarse narrow_fine narrow_ultra highrho_coarse highrho_fine highrho_ultra two_coarse two_fine two_ultra
python validation/comsol/compare.py single_a_ultra narrow_fine highrho_fine two_ultra --meshes 16 --distribution edge_cosine --output validation/comsol/output/edge_cosine_16.json
python validation/comsol/compare.py single_a_ultra narrow_fine highrho_fine two_fine --meshes 20 24 --distribution edge_cosine
python validation/comsol/summarize.py
```

The compact [benchmark record](../validation/comsol/benchmark_results.json)
includes per-case numerical values, the observed mesh orders, source SHA-256
values, and checksums of the local native MPH models. Raw time traces, logs,
and `.mph` files remain under `validation/comsol/output/` and are intentionally
excluded from Git because they are large. The local `.mph` files can be opened
directly in COMSOL for mesh, solver, physics, and material inspection.
