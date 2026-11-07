# Experimental laminar Newton kernel

This optional prototype ports only the laminar station Newton solve. The
Spalart–Allmaras solver remains in Python. From the repository root, with the
cached crates available, build and install the separate extension:

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
solvers are unchanged. `native/bench.py` reports the best of five runs after
one warm-up for Python and both native modes, plus analytic/Python Cf and
thermal differences.
