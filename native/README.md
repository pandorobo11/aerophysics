# Experimental native laminar march

This optional prototype runs the laminar march in Rust: the leading-edge
solve, adaptive steps, BDF1/BDF2 stencils, Newton iterations, acceptance and
separation checks, and output collection. Scalar PPoly edge evaluation and
partial-panel quadrature also run in Rust. The march releases the GIL for
the entire computation. Sutherland and Keyes viscosity models are supported;
the Spalart–Allmaras solver and unsupported models remain in Python. From the
repository root, with the cached crates available, build and install the
separate extension:

```sh
export CARGO_HOME="$PWD/native/.cargo-home"
maturin build --release --offline -m native/Cargo.toml -o native/dist
uv pip install --offline --no-deps --python .venv/bin/python native/dist/aerophysics_native-*.whl
.venv/bin/python native/bench.py
```

`AEROPHYSICS_NATIVE=0` forces Python, `1` requires the extension, and `auto`
(or unset/empty) selects it when available. `uv run` and `uv sync` synchronize
against the main lockfile and remove the extension; use `.venv/bin/python` or
`uv run --no-sync` for native tests.

`AEROPHYSICS_NATIVE_JACOBIAN=analytic` (the default, also unset/empty) uses
an exact Jacobian and pivoted 3x3 block elimination, falling back to banded
LU if a block pivot is zero or a block calculation is non-finite.
`AEROPHYSICS_NATIVE_JACOBIAN=finite-difference` retains the original colored
finite-difference Jacobian and banded LU for comparison. Other values raise
`ValueError` when a native adapter is created. The extension constructor also
accepts `jacobian="analytic"` or `jacobian="finite-difference"`.

Both paths solve the same laminar equations with the same convergence and
damping rules. The analytic path evaluates viscosity at the current node
temperatures only; finite differences also evaluate perturbed temperatures,
which can cause a Keyes applicability warning at a fitted-range endpoint.
This remains an experimental, optional extension; the Python and SA station
solvers are unchanged.

`AEROPHYSICS_NATIVE_MARCH=0` keeps the march and edge evaluations in Python
while retaining the native Newton kernel. Unset, empty, and `1` enable the
native march when applicable; other values raise `ValueError` on the laminar
path. The Python interpolant still constructs the PCHIP coefficients and
integral panels. Rust receives those arrays, the derivative's own coefficients,
quadrature nodes/weights, and the Python march control constants. Output edge
values are reused by Python result construction. On separation, Python builds
the same partial result and extrapolates the separation position as before.
Keyes range excursions in a native call produce one applicability warning.

The physical applicability, units, and conventions are unchanged: planar
laminar perfect-gas flow with an isentropic prescribed edge, SI units, distance
from the sharp leading edge, and positive heat flux into the wall. This port
adds no physical model or validity range. Tests compare PPoly values bitwise
with SciPy (including knots and extrapolation), edge states at `rtol=1e-13`,
and all result fields and Howarth separation with the Python march at
`rtol=1e-8, atol=1e-10`. A Rust flat-plate test checks the Blasius wall slope;
three-path regression tests also compare accepted step positions exactly.

`native/bench.py` reports the best of five runs after one warm-up for Python,
native Newton only (analytic), and the full native march, plus full-march/Python
Cf and thermal differences. For sandbox testing without installation, unzip
the wheel under `tmp/` and prepend that directory and `src/` to `PYTHONPATH`.
Use `.venv/bin/python` to preserve the optional extension.
