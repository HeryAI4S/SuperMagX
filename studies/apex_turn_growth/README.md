# APEX M31: Deterministic turn-count growth, exact on-axis field

**2026-11-06.** This study reuses `supermagx.design.turn_axis.exact_turn_axis_coefficient` and **does not** implement a second electromagnetic solver.

## Configuration and assumptions

- Warm bore 200 mm, winding start radius 130 mm, 20 T target at origin.
- Four physical modules with three independent supply circuits: PS1=M1, PS2=M2, PS3=M3+M4 in series.
- Common peak current in all three circuits: the **minimum possible maximum current** for a fixed geometry if all per-zone coefficients are positive and conductor sections are identical. This is **not** the stress-optimal distribution, and does not optimize terminal voltage.
- All turns 4×8 mm² with uniform azimuthal current-density source; radial pitch 4.3 mm, axial pitch 8.3 mm.
- Intermodule radial Zylon gaps 10/12/8 mm and outer 16 mm, irrespective of layer count; these are geometry bookkeeping, **not validated sufficient reinforcement**.
- Nominal rounded wire R0.5 leaves 31.7854 mm² of metal. Current-density screening uses this area, but magnetic coefficients assume full ideal rectangles. **No exact rounded-wire field calculation is claimed.**
- Each module has the same axial turns-per-layer (60, 64, 70, 80, 100, 120, 160). The allowed winding dimensions and turns remain unrestricted; the grid is a small controlled study, not a global optimization.

## Results

| Radial layers | Turns/layer | Total turns | Equal current for 20 T | Rounded area J screen | Nominal coil OD |
|---|---:|---:|---:|---:|---:|
| 2+2+2+2 | 60 | 480 | 19.725 kA | 620.6 A/mm² | 420.8 mm |
| 2+2+2+2 | 80 | 640 | 18.398 kA | 578.8 A/mm² | 420.8 mm |
| 3+3+3+3 | 80 | 960 | 12.396 kA | 390.0 A/mm² | 455.2 mm |
| 4+4+4+4 | 80 | 1280 | 9.399 kA | 295.7 A/mm² | 489.6 mm |
| 5+5+5+5 | 80 | 1600 | 7.603 kA | 239.2 A/mm² | 524.0 mm |
| 6+6+6+6 | 80 | 1920 | 6.408 kA | 201.6 A/mm² | 558.4 mm |

The original historical M26-A **nonuniform** supply currents are different from this study's equal-current comparison. Therefore the 480-turn row is not a re-evaluation of historical structural loading.

## Physical interpretation

Adding axial turns at constant pitch elongates the coil, so its central-field contribution **saturates** toward the infinite-solenoid limit. It does *not* scale linearly with total turns. Adding radial layers increases the active radial current sheet and typically lowers the necessary current significantly more, but changes the structural lever arm and load sharing. Neither strategy guarantees low Zylon transverse stress, adequate insulation strength, or acceptable axial-force transfer.

**These results are not a structural Pareto front.** Before claiming any design is feasible, compare the exact discrete-field distribution, full Lorentz force, self/mutual inductances, source terminal voltages at feasible ramps, and an independently validated finite-length analytical structural solution.

### Reproduction

```sh
python -m pytest -q tests/test_apex_turn_growth.py
python studies/apex_turn_growth/scan.py
```

The script writes `studies/apex_turn_growth/results.csv`; results are deterministic and do not depend on COMSOL or an empirical stress proxy.
