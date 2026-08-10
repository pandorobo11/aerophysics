# Aerophysics agent instructions

## Completion gate

Before considering a change complete, run the repository's local validation
command from the repository root:

```console
scripts/check.sh
```

This command synchronizes the locked environment, formats Python files, and
then runs the lint, format, type, test, generated-asset, documentation, and
package-build checks. It must pass before handing work back to the user. Use
`scripts/check.sh --check-only` when a check-only run is needed; that mode does
not modify source files.

Generated and verification assets are deliberately excluded from the normal
cross-platform pytest matrix. Run `scripts/check-generated.sh` as part of the
local gate and keep the dedicated generated-assets CI job on Ubuntu with
Python 3.12. If an asset check fails, the failing file paths must remain
visible in the output.

CI must never run a formatter in write mode. CI uses `ruff format --check .`
and must not commit or otherwise apply automatic source changes.

Use the exact uv version required by `pyproject.toml` (`0.9.13`) and keep
dependency installation locked with:

```console
uv sync --all-groups --all-extras --locked
```

Do not weaken the existing Python/OS test matrix, coverage threshold, strict
lint/type settings, or documentation warnings-as-errors policy.
