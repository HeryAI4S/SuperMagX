# First-principles equations and numerical reductions

All lengths are SI metres. Each metal turn occupies a rectangular area
`A_j` in the axisymmetric `r-z` plane and carries a prescribed transport
current `I_j(t)`. Its normalized cross-section current density is either

\[
p_j(r,z)=\frac{1}{r h_j\ln(r_{o,j}/r_{i,j})}
\quad\text{or}\quad
p_j(r,z)=\frac{1}{(r_{o,j}-r_{i,j})h_j}.
\]

The first expression follows from a constant azimuthal electric field around
each circular path, `E_phi(r) ∝ 1/r`, and Ohm's law at constant resistivity.
The second prescribes a uniform external current density, as some FEM models
do. Both satisfy `∫_A p_j dA = 1`.

## Field and magnetic energy

For free space, superpose the exact circular-loop Green kernel over the metal
area: `B(r,z,t) = Σ_j I_j(t) ∫_Aj p_j G_B(r,z;r',z') dA'`. The mutual
coefficient is `M_ij = ∫_Ai∫_Aj p_i p_j M_loop dA_i dA_j`, and
`W_m = Iᵀ M I/2`. This identity makes the field and inductance conventions
consistent for the same prescribed currents.

The loop kernel contains the complete elliptic integral `K(m)`. For almost
coincident loops, forming `1-m` after computing `m` rounds to zero even when
the separation is finite. Instead, the code forms its complementary
parameter directly:

\[
m'=\frac{(r-r')^2+(z-z')^2}{(r+r')^2+(z-z')^2},
\qquad K(m)=K(1-m').
\]

Near the symmetry axis, subtracting elliptic-integral terms before dividing
by `r` also loses digits. Maxwell's source-free expansion gives
`B_r(r,z) = -r B_z'(0,z)/2 + O(r³)` and
`B_z(r,z) = B_z(0,z) - r² B_z''(0,z)/4 + O(r⁴)`.

The axial part of a pair mutual integral originally had two quadrature
coordinates, one for each turn. Their kernel depends only on the difference.
For turn heights `h_a,h_b`, define `H=(h_a+h_b)/2`, `D=|h_a-h_b|/2` and
`u=(z_a-z_b)-(z_{c,a}-z_{c,b})`. The overlap weight is

\[
w(u)=\begin{cases}
\min(h_a,h_b), & |u|\le D,\\
H-|u|, & D<|u|<H,\\
0, & |u|\ge H.
\end{cases}
\]

Thus `∫∫f(z_a-z_b) dz_a dz_b = ∫_{-H}^{H} w(u)
f(z_{c,a}-z_{c,b}+u) du`. Three smooth intervals replace the two-dimensional
axial quadrature. The remaining radial quadratures retain the original turn
profiles. Repeated turn shapes and relative offsets reuse computed values;
the cache does not quantize geometry. Self terms retain a Duffy treatment of
the integrable logarithmic singularity.

## Transport and induced heat

Joule heat follows directly from `∫ρ J_phi² dV dt`. For a prescribed
piecewise-linear current, `∫ I² dt` is evaluated analytically on every time
segment. At fixed `ρ`, one turn has
`R = 2πρ/[h ln(r_o/r_i)]` for the `1/r` profile, or
`R = 2πρ r_c/[(r_o-r_i)h]` for uniform `J_phi`.

The optional **ideal-slab** screening model starts from the one-dimensional
magnetic diffusion equation for a uniform imposed field at both slab faces:

\[
\partial_t B=(\rho/\mu_0)\partial_x^2 B,\qquad
B(\pm d/2,t)=B_{\mathrm{app}}(t).
\]

For odd mode index `q=1,3,5,...`, let `f_q` solve
`τ_q df_q/dt + f_q = dB_app/dt`, with
`τ_q=μ₀d²/(ρπ²q²)`. The volume-averaged dissipated power of that slab axis is

\[
P_d(t)=\frac{8Vd^2}{\rho\pi^4}
\sum_{q\ \mathrm{odd}}\frac{f_q(t)^2}{q^4}.
\]

Because `Σ_{q odd}q⁻⁴=π⁴/96`, the slow-ramp limit is exactly
`P_d=Vd²(dB_app/dt)²/(12ρ)`. Each mode's `∫f_q²dt` is evaluated for the
piecewise-linear pulse; its remaining exponential decay after the last knot
is integrated to infinity. A stable short-step quadrature avoids subtracting
nearly equal terms in this integral.

The winding code adds radial and axial slab-axis estimates driven by the
turn-centre applied `B_z` and `B_r`. This **does not** solve the coupled 2D
current redistribution, terminal circuit, spatial field gradients, or
temperature-dependent resistivity. It is an unfitted alternative to the
existing calibrated one-pole screening model, not a verified FEM replacement.

## Mathematical references

- [NIST DLMF §19.12](https://dlmf.nist.gov/19.12) gives the near-singular
  elliptic-integral expansion in the complementary modulus.
- [SciPy `ellipkm1` documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipkm1.html)
  explains why supplying the complement preserves precision near `m=1`.
- [MIT OpenCourseWare, Electromagnetic Fields and Energy, Chapter 10](https://ocw.mit.edu/courses/res-6-001-electromagnetic-fields-and-energy-spring-2008/pages/chapter-10/)
  develops magnetic diffusion through thin conductors. The odd-mode heat
  formula above is derived for the stated 1D slab boundary-value problem.
