# Accuracy and validation status

## Verified numerically in the repository

The test suite checks the finite-rectangle on-axis field against an
independent two-dimensional adaptive integral, arbitrary-zone superposition,
mutual reciprocity and positive eigenvalues for a small winding, analytic
ring resistance for both `1/r` and uniform source profiles, exact
piecewise-linear `I²R` integration, and complete
post-waveform eddy tail bookkeeping. Run `python -m pytest -q` to repeat.

These tests check implementation consistency. They do not validate material
data or assert blanket agreement with a finite-element model.

## External FEM comparison

The induced-loss model is explicitly labeled as screening. A match on one
geometry and waveform does not establish a transferable error bound across
conductor sizes, layouts or pulse durations. The package has not yet been
validated against converged FEM solutions over a specified design domain.

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
