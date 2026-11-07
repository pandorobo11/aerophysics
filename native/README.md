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
