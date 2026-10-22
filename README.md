# SuperMagX

SuperMagX is a modular Python project for pulsed-magnet calculations. Phase 1
contains a parameterized axisymmetric **magnetostatic** solver and separate
transport/induced loss models. Thermal and mechanical solvers are module
boundaries for subsequent phases; they are not implemented yet.

## Install and run

```powershell
python -m pip install -e ".[test]"
python -m pytest -q
python -m supermagx examples/two_zone_case.json result.json
```

The package accepts any number of named supply zones. Each physical turn has
its own rectangular metal cross-section, axial position, zone, and polarity.
`radial_stack()` generates regular groups with independently adjustable layer
counts, turns per layer, dimensions, pitches and radial gaps. For irregular
windings, construct `Magnet(tuple(RectTurn(...), ...))` or use JSON
`geometry.kind = "explicit_turns"`. Reusing a zone identifier represents
turns connected in series. Each zone's current may have a different
piecewise-linear waveform.

## Modules and data flow

```text
geometry + excitation
       |         |
       v         v
electromagnetics ---> losses
       |               |
       v               v
 mechanics          thermal
 (future)           (future)
```

| Module | Phase-1 function |
| --- | --- |
| `geometry` | Finite rectangular turns and arbitrary radial-group/zone mapping |
| `excitation` | Independent zone currents on a piecewise-linear time grid |
| `electromagnetics` | Free-space axisymmetric field at air points and finite-cross-section zone mutual matrix |
| `losses.joule` | Exact time integration of fixed-resistivity `I²R`; resistance follows each turn's selected current-density profile |
| `losses.eddy` | Turn-centre screening with either calibrated one-pole response or unfitted ideal-slab diffusion modes; both include the post-waveform tail |
| `losses.network` | Coupled conducting-cell magnetic diffusion and pulse heat, including the post-waveform eddy tail; compiled operator for repeated waveform and resistivity evaluations |
| `thermal`, `mechanics` | Reserved boundaries; no temperature, displacement or stress field is claimed |

`python -m supermagx` reads a JSON case and writes a JSON result. The output
states the field-probe currents, quadrature order/change, per-zone losses,
and the status of the unsolved modules. The sample input is
[`examples/two_zone_case.json`](examples/two_zone_case.json).
The coupled-loss example is
[`examples/coupled_loss_case.json`](examples/coupled_loss_case.json).

## Numerical model and precision

Each turn selects `source_profile = "dc_equilibrium"` (default,
`J_phi = I/(r H)`, `H = ∫_metal dr dz / r`) or `"uniform"`
(`J_phi = I/A`, as in a uniform external-current-density FEM feature).
The choice is available on `RectTurn` and `RadialGroup` JSON entries and
controls both magnetic sources and Joule resistance. Axisymmetric
circular-filament Green functions are
integrated over the finite metal rectangles. `field()` doubles Gauss order
until the field change is below the requested tolerance, and raises
`ConvergenceError` otherwise. The finite-turn mutual calculation uses
adaptive pair quadrature and a Duffy transformation for the logarithmic
self-term. The pair integral analytically collapses the two axial coordinates
to their trapezoidal difference distribution. Near-coincident loops use the
complementary elliptic parameter formed from physical separation; near-axis
radial fields use the Maxwell-series limit. There is no empirical
magnetic-field multiplier in this package.

The `uniform` option explicitly counts the prescribed metal current's
`I²R` heat in Python. COMSOL's External Current Density feature excludes
the external source's own Joule loss by default; a like-for-like comparison
must enable its loss contribution or provide that heat source separately.

For JSON losses, `"eddy": {"model": "slab_modes", "modes": 16}` selects an
unfitted ideal-slab diffusion expansion. The existing `"eddy": {"gamma":
...}` selects the calibrated one-pole screen. The slab model resolves odd
diffusion modes in each conductor dimension but still uses the applied field
at the turn centre and independent slab axes. Its tail heat is integrated
after the last waveform knot. See the [derivation](docs/first_principles.md).

`"eddy": {"model": "coupled_network", "radial_cells": 4,
"axial_cells": 4, "cell_distribution": "edge_cosine",
"refinement_meshes": [[2, 2], [3, 3], [4, 4]]}` selects
the coupled cell model. The same `resistivity_ohm_m` is used for transport
and induced currents. Output separates `joule_j`, `eddy_before_tail_j`,
`eddy_tail_j`, and `total_j`; `network_refinement` shows mesh-to-mesh changes.
Induced heat is not attributed to individual power zones because its modes
couple zones. Increase cell counts until **eddy heat itself** stabilizes for
the case of interest. The slab model remains useful for broad scans because
the network's dense mutual matrix and eigenproblem can be expensive. For
repeated scans at one geometry, use `compile_coupled_loss(magnet, ...)` once
and `operator.evaluate(waveform, resistivity=rho)` for each candidate; see
[`examples/compiled_loss_scan.py`](examples/compiled_loss_scan.py). The
`edge_cosine` setting places more cells near the metal boundary. Both mesh
families require a convergence check; the measured one/two-ring COMSOL
comparisons and full continuum derivation are in
[`docs/loss_theory_validation.md`](docs/loss_theory_validation.md).
The coupled model requires `source_profile = "dc_equilibrium"` for every turn;
it rejects a prescribed uniform-current source profile because that profile
is inconsistent with its resistive equilibrium and exact heat split.

The modal pulse integrals are evaluated in batches over diffusion modes.
The COMSOL benchmark pins the revolved-volume convention explicitly and
checks it against cross-section integration with the `2*pi*r` Jacobian.
The stand-alone [Maxwell-to-model theory reference](docs/theory_from_maxwell.md)
derives the field, transport and induced heat, thermal, and structural
equations and records the implementation and validation boundaries.

This is numerical convergence **within the stated free-space, prescribed-DC-
current model**. Consecutive quadrature agreement is not a mathematical error
certificate or a guarantee of 1% agreement with COMSOL. The precision-
controlled `field()` accepts observation points outside metal; the separate
`screening_coefficients()` used for eddy estimates samples turn centres and
does not have the same accuracy claim. Iron, leads, joints, actual supply
circuits, conductive supports, temperature-dependent resistivity and
mechanical preload are outside this phase. The one-pole and slab eddy models
remain screening approximations. The cell model has independent COMSOL 6.4
checks on representative one/two-ring cases, but the full winding and all
design variants are not certified by those comparisons.

See [accuracy notes](docs/accuracy.md) before using results as a design gate.

## APEX CuNb narrow-conductor study

The [op5 screening study](studies/apex_cunb_op5/op5_report.md) scans 70
rectangular CuNb conductors with radial metal width at least 5 mm, optimizes
the symmetric pulse ramps, and compares total pulse heat, force proxies,
zone voltage, and average-turn temperature rise. The
[reproducible script](studies/apex_cunb_op5/scan_narrow.py) uses the op3
calibrated one-pole loss approximation and writes the tabulated results in
the same directory. It is a screening study; the new geometry has no
structural or transient thermal FEM validation. A later op5 coupled COMSOL
run is recorded in the study report; its heat-source/thermal-energy balance
still needs resolution before the analytic screen is accepted as a validated
thermal prediction.

## Licensing

MIT; see [LICENSE](LICENSE).
