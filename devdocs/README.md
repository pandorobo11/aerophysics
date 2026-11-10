# Developer documentation

These notes record implementation details that matter to contributors but are
too detailed for the user manual in `docs/`. The user manual explains what each
model calculates, its equations and theory, and when it applies. The files here
explain how the numerical solvers, range checks, and GUI internals behave.

- [Numerical implementation notes](numerical_implementation.md): root
  brackets, tolerances, convergence classification, and floating-point
  safeguards used by the calculation core.
- [GUI internals](gui_internals.md): settings JSON fields, settings
  validation, and how GUI pages classify and display solver failures.
- [Marching boundary-layer design](marching_boundary_layer_design.md): design
  and verification records of the variable-edge boundary-layer solver
  (written in Japanese).
- [Axisymmetric characteristics design](axisymmetric_characteristics_design.md):
  design, verification measurements, and review records of the axisymmetric
  method of characteristics.

The user manual keeps to what a user needs in order to choose and call a
model: what it calculates, its assumptions and limits, and the meaning of its
results. Solver algorithms, tolerances, grid construction, refinement and
timing measurements, and settings-file compatibility belong here.

The contributor workflow, validation gate, and release process are described
in [`DEVELOPMENT.md`](../DEVELOPMENT.md).

When behavior described here changes, update this directory in the same
change. When the change is visible to users (new inputs, new outputs, changed
ranges, or changed warnings), update the user manual as well.
