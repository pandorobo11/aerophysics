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

### Generated-asset failure procedure

When a generated or verification asset check fails, follow this procedure:

1. Regenerate the assets with the existing generator that owns them. For the
   complete verification set, use:

   ```console
   uv run python docs/scripts/generate_verification.py
   ```

   Run the applicable individual `generate_*.py` script when only one asset
   family is involved.
2. Inspect the Git diff after regeneration. Confirm that every numerical
   value, table, and SVG change is intentional and supported by the source or
   model change being made.
3. If any change is unexpected, investigate the generator, dependencies,
   platform, or source data. Do not mechanically commit regenerated files just
   to make CI pass.
4. If the diff is understood and correct, rerun `scripts/check.sh` and confirm
   that every check succeeds before considering the work complete.

Never hand-edit generated files to make CI pass. Fix the owning generator or
its inputs, regenerate the assets, and review the resulting diff instead.

CI must never run a formatter in write mode. CI uses `ruff format --check .`
and must not commit or otherwise apply automatic source changes.

Use the exact uv version required by `pyproject.toml` (`0.9.13`) and keep
dependency installation locked with:

```console
uv sync --all-groups --all-extras --locked
```

Do not weaken the existing Python/OS test matrix, coverage threshold, strict
lint/type settings, or documentation warnings-as-errors policy.
