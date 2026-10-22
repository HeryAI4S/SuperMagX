# Accuracy and validation status

## Verified numerically in the repository

The test suite checks the finite-rectangle on-axis field against an
independent two-dimensional adaptive integral, arbitrary-zone superposition,
mutual reciprocity and positive eigenvalues for a small winding, analytic
ring resistance for both `1/r` and uniform source profiles, exact
piecewise-linear `I²R` integration, near-coincident and near-axis kernel
limits, axial convolution against an independent double integral, the
short-ramp diffusion limit, and complete post-waveform eddy tail bookkeeping.
Run `python -m pytest -q` to repeat.

These tests check implementation consistency. They do not validate material
data or assert blanket agreement with a finite-element model.

## External FEM comparison

The coupled conducting-cell model has an independent COMSOL 6.4 benchmark on
one and two solid rectangular CuNb rings. The cases vary conductor size,
resistivity, ramp duration, and source timing. The benchmark keeps the
post-pulse diffusion tail and uses the physical axisymmetric power integral.
The archived benchmark used `IntSurface(2*pi*r*mf.Qh)` with cross-section
integration. The current benchmark sets `intvolume=on` and integrates
`mf.Qh`, checking the two conventions independently. See the
[derivation, cases, and reproducible commands](loss_theory_validation.md).
This is a bounded verification of the
stated free-space, constant-resistivity, continuous-ring equations; it is not
a validation of the complete APEX winding or composite conductor topology.

The slab-mode option remains a screening approximation for an **ideal 1D
slab**. Its centre-field and independent-axis approximations have no such
COMSOL comparison here. A total-heat match can conceal a larger error in the
smaller eddy or after-pulse contribution, so those quantities are reported
separately. The FEM heat itself depends on time-output and spatial refinement.

## Required engineering gates before claiming 1% vs FEM

1. Specify a bounded design family, probe locations, waveforms, material
   constitutive laws and the list of output quantities.
2. Refine the FEM mesh, time step and solver tolerances independently for
   each quantity, especially local field, voltage and peak stress.
3. Compare this package at identical metal geometry, current-density profile,
   current history, zone orientation, magnetic boundary conditions and
   resistivity. Keep source case identifiers and hashes.
4. Use held-out designs. If a requested quantity is not within its own
   tolerance, use the FEM result and improve the relevant module.

The field quadrature tolerance is a numerical stopping criterion; it is not a
bound on geometry/model error. Loss agreement does not establish temperature
or stress accuracy. The phase-one package does not compute those fields.
