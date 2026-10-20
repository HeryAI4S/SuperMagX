# Module architecture

The geometry is the shared source of turn locations, conductor dimensions,
polarity and electrical zone identity. It does not call an electromagnetic,
loss, thermal or mechanical solver. Excitation is independent of geometry.

Electromagnetics consumes geometry and a specified vector of zone currents.
Its zone coefficient and mutual matrix outputs keep the original zone order.
The transport-loss module consumes geometry, fixed resistivity and a current
waveform. The induced-loss module consumes geometry, waveform and an explicit
diffusion calibration, and computes its own turn-centre field coefficients.

The planned thermal module will consume a time- and space-resolved heat source,
not just total joules. The planned structural module will consume distributed
Lorentz loads plus mechanical material, boundary, contact and cooldown data.
The phase-one energy and field outputs alone do not determine a hotspot or a
von Mises peak.

The present coupling graph is one-way. A future temperature-dependent
resistivity model will require iteration between electromagnetic, loss and
thermal modules with a documented convergence criterion.

## Input conventions

- SI metres, seconds, amperes, ohm-metres, tesla, henries, joules.
- One `RectTurn` is one complete circular metal turn. A shared zone ID means
  the turns use the same series current; the user controls polarity per turn.
- `source_profile` chooses `dc_equilibrium` (`1/r`) or `uniform` current density
  per turn; this choice also determines its fixed-resistivity Joule resistance.
- `radial_stack()` generates regular layouts. `explicit_turns` allows
  nonuniform positions and cross-sections, with any number of supply zones.
- Conductors are normal metal with constant scalar resistivity. Passive
  insulation and reinforcement geometries are not yet represented as domains.
- A waveform is linear between knots; after its last knot, the current is
  held constant and the model's residual eddy state decays fully.
