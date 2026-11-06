# APEX M30 — discrete-turn exact on-axis field and three-supply feasibility (2026-11-06)

**Status:** an independent closed-form electromagnetic design gate; not a structural model, not a transient circuit, not a Pareto optimum.

## Derivation

For each conductor occupying `a≤r≤b`, `zl≤z≤zu`, with prescribed normalized azimuthal profile `p(r,z)`, direct Biot–Savart integration yields:

`Bz(0,z0)/I = (μ0/2) ∫a^b ∫zl^zu p(r,z) r²/[r²+(z-z0)²]^(3/2) dz dr`.

For `uniform` source, `p=1/[(b-a)h]`, `h=zu-zl`. Set `G(c)=c[asinh(b/|c|)-asinh(a/|c|)]` (G(0)=0). Then

**`Bz/I=μ0[G(zu-z0)-G(zl-z0)]/[2(b-a)h]`.**

For `dc_equilibrium`, `p=1/[r h ln(b/a)]`. Set `H(c)=asinh(c/a)-asinh(c/b)`. Then

**`Bz/I=μ0[H(zu-z0)-H(zl-z0)]/[2h ln(b/a)]`.**

The turn direction signs enter linearly; group equal series-circuit identifiers to obtain `Bz=Σz K_z I_z`. These are exact for prescribed *rectangular* source distributions, not fully transient currents. Existing `AxisymmetricMagnetostatics` remains the full Green-kernel off-axis/inductance solver.

## Strict local current density

For a rectangular metal turn, `Icap=Jmax(b-a)h` for `uniform`. For `dc_equilibrium`, the peak J occurs at the **inner radial edge**, so `Icap=Jmax a h ln(b/a)`. All turns sharing a power supply carry the same current; the series circuit's cap is its most restrictive turn (and the supply limit).

With nonnegative currents in a *fixed* geometry, the largest positive on-axis field is the sum of positive per-zone field coefficients times their respective current caps. If `I<25 kA` is strict and that cap is active, the mathematical boundary is a **supremum**, not necessarily attainable exactly. A bound below 20 T proves that geometry infeasible under the assumed limits. A bound above 20 T does **not** prove structural, voltage or thermal feasibility.

For APEX require exactly three supplies; allow multiple mechanical modules per supply, but one cross-section per supply in the design workflow.

## M26-A independently reproduced

M26-A is rebuilt from its original 480 rectangular turn centers, 4×8 mm² section, 60 turns/layer, 8.3 mm axial pitch, and PS1=M1, PS2=M2, PS3=M3+M4. The validation is independent of an FEM field computation.

| Quantity | Closed-form result |
|---|---:|
| PS1 field coefficient | 0.266574835 T/kA |
| PS2 field coefficient | 0.258106714 T/kA |
| PS3 field coefficient | 0.489263177 T/kA |
| Historical peak currents | 14.600785 / 17.402973 / 23.740635 kA |
| Historical Bz(0,0) | 19.999444887 T |
| Historical Bz(0,+2.5 mm) | 19.998810354 T |
| Current cap per supply under J≤600 A/mm² | 19.200 kA |
| **Fixed-geometry B supremum under caps** | **19.467738741 T** |
| Target-to-supremum ratio | 1.027340682 |

The historical M26-A geometry therefore **cannot meet 20 T and the 600 A/mm² screening limit simultaneously**, even if all three supplies operate at their 19.2 kA rectangular-metal current ceiling. This is a *necessary magnetic* finding and does not certify any proposed replacement geometry.

**Rounded-corner caveat:** real 4×8 R0.5 mm² metal has area `32−(4−π)×0.5² = 31.7854 mm²`, giving `19.0712 kA` at 600 A/mm² under *uniform* current density. The closed-form model uses exact **rectangles**, not actual rounded boundaries: scaling the 19.4677 T field by this area ratio is only a screening approximation, not an exact rounded-metal magnetic field.

## Validation and scope

Run `python -m pytest -q tests/test_turn_axis_design.py` and `python examples/apex_m26_field_envelope.py`.

Tests use independent 2D numerical Biot–Savart quadrature for both source profiles and multiple axial positions, circuit sign mapping, series-current caps, and the full 480-turn M26-A regression. Field agreement does **not** validate mechanical stress, interface contact, pulse diffusion, terminal voltage or exact rounded conductor geometry.

**Next physics gate:** derive a finite-length, axially resolved force and composite reinforcement load-transfer model before structural Pareto ranking.
