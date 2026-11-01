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

The contributor workflow, validation gate, and release process are described
in [`DEVELOPMENT.md`](../DEVELOPMENT.md).

When behavior described here changes, update this directory in the same
change. When the change is visible to users (new inputs, new outputs, changed
ranges, or changed warnings), update the user manual as well.
