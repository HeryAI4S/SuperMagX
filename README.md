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
| `thermal`, `mechanics` | Reserved boundaries; no temperature, displacement or stress field is claimed |

`python -m supermagx` reads a JSON case and writes a JSON result. The output
states the field-probe currents, quadrature order/change, per-zone losses,
and the status of the unsolved modules. The sample input is
[`examples/two_zone_case.json`](examples/two_zone_case.json).

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

For JSON losses, `"eddy": {"model": "slab_modes", "modes": 16}` selects an
unfitted ideal-slab diffusion expansion. The existing `"eddy": {"gamma":
...}` selects the calibrated one-pole screen. The slab model resolves odd
diffusion modes in each conductor dimension but still uses the applied field
at the turn centre and independent slab axes. Its tail heat is integrated
after the last waveform knot. See the [derivation](docs/first_principles.md).

This is numerical convergence **within the stated free-space, prescribed-DC-
current model**. Consecutive quadrature agreement is not a mathematical error
certificate or a guarantee of 1% agreement with COMSOL. The precision-
controlled `field()` accepts observation points outside metal; the separate
`screening_coefficients()` used for eddy estimates samples turn centres and
does not have the same accuracy claim. Iron, leads, joints, actual supply
circuits, fully coupled magnetic diffusion, temperature-dependent resistivity
and mechanical preload are outside this phase. Both eddy models remain
screening approximations even if the magnetostatic field integral is highly
converged.

See [accuracy notes](docs/accuracy.md) before using results as a design gate.

## Licensing

MIT; see [LICENSE](LICENSE).
