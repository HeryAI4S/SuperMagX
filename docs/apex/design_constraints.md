# APEX Cu–Nb pulsed magnet — current constraint register (2026-10-28)

**Status: research baseline, not frozen engineering design.** Values below must not silently overwrite historical case data. Static M1 targets must not be confused with dynamically verified performance.

| Type | Quantity | Value | Scope |
|---|---|---|---|
| Target | On-axis center field | **20 T** | M1 static target; dynamic plateau checked in M2 |
| Geometry | Warm clear bore diameter | **200 mm** | Must include cryostat clearance |
| Baseline | Cu–Nb winding inner diameter | **260 mm** | Nominal, not fully frozen |
| Operating | Initial coolant temperature | **77 K**, LN2 immersion | Not a thermal calculation |
| Hard | Peak current per independent supply | **<25 kA** | Strict inequality per supply |
| Hard | Absolute terminal voltage per supply | **≤30 kV** | Entire pulse; M2 circuit model |
| Hard | Number of independent supplies | **≥3** | 3–6 preliminary search range, not global upper hard bound |
| Hard | Magnetic-field flat-top target | **1 ms** | Dynamic field plateau to be verified |
| Free | Rise/fall ramp durations | `t_rise`, `t_fall` variables | Previous fixed 3/1/3 ms discarded |
| Topologies | Internal support | S0 no ring / S1 slotted G10–GFRP ring | Separate discrete cases |
| S1 variable | Support thickness | **2–12 mm** preliminary | S0 uses 0 mm |
| S1 fixed | Coolant-channel orientation | Straight axial slots | No helical/interconnected slots in present S1 |
| Structural | Main support materials | Zylon, G10/GFRP | No metal load-bearing shell or reinforcement ring |
| Research goal | Cu–Nb peak von Mises stress | Approximately **700 MPa** | **Not** a verified material allowable or finalized safety factor |

**Undetermined:** detailed metal-section shape, total turns, axial and outer dimensions, layer/zone mapping, current waveform per supply, load paths, contact/prestress, 77 K Cu–Nb resistivity and strength as measured, thermal hot-spot and cycling limits. M1 should not silently supply numbers for them.

**Compatibility with existing code:** `RectTurn.source_profile` and radial group cases currently permit `dc_equilibrium` (default, `J_phi∝1/r`) and `uniform` (`J_phi=I/A`). M1 source selection, zone resistance and FEM comparisons must use the same profile. Existing `losses.network` requires `dc_equilibrium`; do not force `uniform` into its assumptions.

**Historical study separation:** [OP5](../../studies/apex_cunb_op5/op5_report.md) dates to 2026-10-21 and used a 29.785 kA, 8-supply candidate with former <40 kA / <40 kV selection logic. That candidate violates the current strict <25 kA limit. Its `r J_phi B_z` and local force-density estimates are **stress proxies, not computed structural von Mises**. Retain results for provenance; do not relabel as new constraints or verified design.

Module ownership: M1 magnetostatic field/inductance/force; M2 supply circuits, diffusion and Joule/eddy loss; M3 thermal/cryogenic; M4 preload/contact/structural failure; M5 materials and uncertainty. A static field or energy match does not establish system feasibility.
