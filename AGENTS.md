# Aerophysics agent instructions

## Git workflow

Use a short-lived branch and pull request for Codex code changes. Start from
the latest `main`, use a descriptive prefix such as `feature/`, `fix/`,
`refactor/`, or `chore/`, and wait for CI to succeed. Do not merge a pull
request unless the user explicitly requests or approves that specific merge.
When a merge is authorized, use squash merge. Do not push code changes directly
to `main`; direct pushes are permitted only for clearly trivial edits such as
README typos or comment-only corrections.

Before opening a pull request, run the relevant tests, Ruff lint and format
checks, and mypy. For changes to a physical model or numerical correlation,
confirm and document its applicability or validity range, units, sign and
coordinate conventions, and analytical, literature, or regression
verification.

When adding, removing, renaming, or changing a required CI job in
`.github/workflows/ci.yml`, update the `ci-gate` job's `needs` dependencies and
result checks in the same change. `CI gate` must cover every required CI job.

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

The slow marching-solver tier runs automatically when solver-relevant files
changed relative to `origin/main`, including staged, unstaged, and untracked
files. `scripts/needs-verification.sh solver` applies the path rules listed in
`scripts/solver-verification-paths.txt` and runs the tier if the base cannot be
resolved. The list holds the transitive `aerophysics` import closure of the
`solver_verification` tests and the reference files they read;
`tests/test_verification_paths.py` fails when it drifts, so update the list
when a solver module gains an import or a solver test gains an input. Use `scripts/check.sh --full` to force
both solver and package verification, and use `--full` before changes whose
impact on the solver or distribution is unclear. The `--full` and
`--check-only` options can be combined.

By explicit user decision on 2026-11-07, package building and clean-install
verification are also conditional locally and in pull-request CI.
`scripts/needs-verification.sh package` owns these rules: packaging metadata,
the `hatch_build.py` hook, entry points/documentation lookup, packaging tests,
documentation configuration, README/license files, gate scripts/workflows,
and additions/deletions/renames under `src/` or `docs/`. The wheel hook's
documentation generators and reference data also trigger verification because
they are copied and executed in its isolated build tree. Ordinary source or
documentation content edits outside these rules do not trigger the package
tier. Main and release CI always run it; skipped PR jobs report the reason and
remain subject to `ci-gate`. Both modes use `origin/main` by default;
`VERIFICATION_BASE` overrides it (the former `SOLVER_VERIFICATION_BASE` remains
a fallback). Missing bases or failed comparisons run verification conservatively.

Generated and verification assets are deliberately excluded from the normal
cross-platform pytest matrix. Run `scripts/check-generated.sh` as part of the
local gate and keep the dedicated generated-assets CI job on Ubuntu with
Python 3.12. If an asset check fails, the failing file paths must remain
visible in the output.

By explicit user decision on 2026-11-07, slow `solver_verification` tests run
in CI only on Ubuntu / Python 3.12, separately from generated assets. They run
for solver-relevant pull-request changes and always on main and in releases.
Representative fast marching tests stay in the full Python/OS matrix. In the
local write-mode gate, verification assets are regenerated before the asset
tests; only their redundant currency test is skipped after successful
regeneration. Check-only, CI, and release runs retain that currency check.

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

Do not require an exact uv version in `pyproject.toml`. Local development may
use a current compatible uv version. CI and release workflows pin uv explicitly
to keep the canonical validation and release environments reproducible. Keep
dependency installation locked with:

```console
uv sync --all-groups --all-extras --locked
```

Do not change the CI or release uv pin as part of unrelated changes.

Do not weaken the existing Python/OS test matrix, coverage threshold, strict
lint/type settings, or documentation warnings-as-errors policy.

## Release readiness

After completing a user-facing feature, behavior change, or important fix and
after the completion gate passes, briefly evaluate whether the accumulated
changes since the latest version tag form a suitable release. Use the Git
history, public API changes, the `CHANGELOG.md` Unreleased section, and the
state of tests, verification, and documentation. Do not judge readiness from
commit counts or generated-file line counts alone.

Recommend a release only when the accumulated changes form a coherent,
user-visible release theme or contain an important calculation, public-API,
or distribution fix. When recommending one, state the suggested Semantic
Versioning bump, the reasons, and any remaining release work. If a release is
not warranted, do not mention release readiness in the final response.

This is an advisory check only. Never change the package version, create a
release commit or tag, push, or publish a GitHub Release unless the user
explicitly requests it.
