# Development

This document is the contributor workflow for `aerophysics`, including local
validation, generated documentation assets, and releases. User-facing model
documentation belongs in the English Sphinx manual under `docs/`. Solver
internals, numerical tolerances, and GUI settings-file details belong in
[`devdocs/`](devdocs/README.md) instead.

The user manual is opened from the GUI help page, so write it for users:
start each model page with a plain-language physical background, keep the
equations and assumptions, and avoid internal implementation detail.

## Contribution workflow

Use a lightweight GitHub Flow for normal feature, fix, and refactoring work:

1. Start a short-lived branch from the latest `main`, using a descriptive name
   such as `feature/...`, `fix/...`, `refactor/...`, or `chore/...`.
2. Before opening a pull request, run the relevant tests, Ruff lint and format
   checks, and mypy. Run the complete local gate described below before handing
   off the change.
3. Open a pull request to `main` and wait for CI to pass. An authorized
   maintainer then squash merges it.

Do not push normal code changes directly to `main`. Direct pushes are reserved
for clearly trivial edits such as README typos or comment-only corrections.
Reviewer approval is welcome but is not a required merge condition.
Codex must not merge a pull request unless the user explicitly requests or
approves that specific merge.

For changes to a physical model or numerical correlation, state and verify the
applicability or validity range, units, sign and coordinate conventions, and
the supporting analytical, literature, or regression verification in the pull
request.

Whenever a required CI job is added, removed, renamed, or changed in
`.github/workflows/ci.yml`, update the `ci-gate` job's `needs` dependencies and
result checks in the same pull request. `CI gate` must remain an exhaustive
summary of every required CI job.

## Environment

The project requires Python 3.12 or newer. Local development does not require
an exact `uv` version. CI and release workflows use a pinned uv version as the
canonical reproducible validation environment. Synchronize every locked
dependency group and extra from the repository root:

```console
uv sync --all-groups --all-extras --locked
```

Do not update the lock file as part of an ordinary environment sync.

## Dependency audit

Run the dependency audit explicitly after dependency updates and before a
release:

```console
scripts/check-dependencies.sh
```

The command checks every locked group and extra, including hashes, against
the current vulnerability database. It requires network access and is not
part of the local completion gate or required pull-request CI. Review any
findings for affected versions and reachable code paths, and update the lock
when a compatible fix is available.

## Local validation

Run the complete local gate before handing back a change:

```console
scripts/check.sh
```

The gate synchronizes the locked environment, formats Python in write mode,
and then runs Ruff lint and format checks, mypy, the normal test suite, the
generated-asset checks, and warning-as-error Sphinx HTML and doctest builds.
When package verification is selected, it builds the wheel and source
distribution, then installs the wheel with its GUI extra and the sdist into
separate clean virtual environments. Those
installed-package checks run from outside the checkout and cover metadata,
runtime dependencies, a public calculation, the console entry point, and
bundled documentation lookup. To run the same checks without rewriting Python
source, use:

```console
scripts/check.sh --check-only
```

Slow marching-solver verification runs automatically for solver-relevant
changes relative to `origin/main`, including staged, unstaged, and untracked
files. The path rules are in `scripts/needs-verification.sh solver`; a missing
base or failed comparison runs the tier conservatively. Force it with
`scripts/check.sh --full` (combinable with `--check-only`), and use `--full`
when a change's solver or distribution impact is unclear. `--full` forces
both solver and package verification.

By explicit user decision on 2026-11-07, package builds and clean-install
checks also run conditionally, using `scripts/needs-verification.sh package`.
They run for changes to packaging metadata, `hatch_build.py`, public package
initialization, GUI launcher/documentation lookup, package tests, `docs/conf.py`,
README/license files, and gate scripts/workflows. File additions, deletions,
and renames anywhere under `src/` or `docs/` also trigger them, including
untracked files. Ordinary content edits elsewhere in those trees do not.
`docs/scripts/` and `tests/reference_data/` also trigger the package tier:
the wheel hook copies these inputs into an isolated tree and runs the
generators before bundling HTML. Main and release CI always verify packages;
PR CI may skip the build/install steps with a reason in the job summary.

Both decision modes compare with `origin/main` and inspect the working tree.
Set `VERIFICATION_BASE` to use another base; `SOLVER_VERIFICATION_BASE` is
retained as a fallback for existing callers. Missing bases and Git failures
run verification conservatively.

The test tiers are mutually exclusive:

| Tier | Command / selection |
| --- | --- |
| Normal tests and coverage | `scripts/check-coverage.sh` (excludes the other three markers) |
| Generated assets | `scripts/check-generated.sh` |
| Slow solver verification | `uv run pytest -m solver_verification --no-cov` |
| Installed packages | `uv run pytest -m package_artifact --no-cov` (after building) |

By explicit user decision on 2026-11-07, slow solver verification runs in one
Ubuntu / Python 3.12 CI job for solver-relevant PR changes, and always on main
and in releases. Representative fast marching tests remain in the full
Python/OS matrix.

The normal test suite collects branch coverage once and applies independent
gates to the two product layers:

- the numerical core (top-level modules under ``src/aerophysics``) must remain
  at or above 95%;
- the GUI package (modules under ``src/aerophysics/gui``) must remain at or
  above 90%.

Run those same tests and gates on their own with:

```console
scripts/check-coverage.sh
```

Only ``gui/app.py`` is excluded. ``AppTest`` executes that declarative
bootstrap as a script through Streamlit's runner, and those executed lines are
not attributed to the importable module by the coverage collector. The page
renderers, shared components, launcher, and every other GUI module remain
measured by the GUI gate.

CI must only check formatting (`ruff format --check .`); it must never apply
formatting or commit generated changes.

For a focused documentation preview:

```console
uv run sphinx-build -W --keep-going -b html docs docs/_build/html
uv run sphinx-build -W --keep-going -b doctest docs docs/_build/doctest
```

When launched from a source checkout, `aerophysics-gui` automatically uses a
valid `docs/_build/html` tree. Wheel builds render the same manual into
`aerophysics/_docs/`, which lets the installed GUI serve it offline.

## Generated and verification assets

Generated RST, tables, and SVGs are owned by scripts under `docs/scripts/`.
Never hand-edit a generated file to make a check pass.

Regenerate the complete verification set with:

```console
uv run python docs/scripts/generate_verification.py
```

The viscosity comparison is generated separately:

```console
uv run python docs/scripts/generate_viscosity_comparison.py
```

When only one asset family changed, run its applicable `generate_*.py` script.
Each generator accepts `--check` to compare the committed files without
writing them.

| Verification record | Generator |
| --- | --- |
| Standard atmosphere | `docs/scripts/generate_standard_atmosphere_validation.py` |
| Compressible flow | `docs/scripts/generate_compressible_flow_validation.py` |
| Thermophysical properties | `docs/scripts/generate_thermophysical_validation.py` |
| Viscous flow | `docs/scripts/generate_viscous_flow_validation.py` |
| Viscosity model comparison | `docs/scripts/generate_viscosity_comparison.py` |

The fixtures behind the thermophysical record are rebuilt separately:

```console
uv run python docs/scripts/build_transport_reference.py
uv run python docs/scripts/generate_nist_transport_reference.py
```

Refresh external snapshots only in isolated environments:

```console
uv run --isolated --with cantera==3.2.0 python docs/scripts/capture_cantera_reference.py
```

The standard-atmosphere generator rewrites its RST record and SVGs; with
`--check` it checks only the deterministic SVGs.
After every regeneration, inspect the Git diff and confirm that each numerical
value, table, and SVG change follows from an intentional source or model
change. Unexpected output must be investigated in the generator, dependency,
platform, or source data. Once the diff is understood, rerun
`scripts/check.sh`.

The standard-atmosphere RST is a generated verification record of the measured
values from the current checkout, not a byte-for-byte numerical gate.
`scripts/check.sh`, documentation workflows, release workflows, and the
wheel build regenerate the verification set before their Sphinx build. The
check-only local gate remains non-destructive; it checks the checked-in
deterministic assets and builds the existing record without rewriting it.

The standalone generated-asset gate is:

```console
scripts/check-generated.sh
```

Generated verification tests remain outside the normal cross-platform pytest
matrix and run in their dedicated Ubuntu/Python 3.12 CI job. Failure output
must retain the paths of differing assets.

The local write-mode gate regenerates verification assets before this step,
then sets `AEROPHYSICS_SKIP_ASSET_CURRENCY=1` only after successful generation
to skip `test_deterministic_verification_assets_are_current`. Every other
asset test still runs. Check-only mode runs the currency test without writing;
CI and release workflows must never set the skip variable.

## Releases

Before proposing a release, review the changes since the latest version tag,
the public API, the `CHANGELOG.md` Unreleased section, and the state of tests,
verification, documentation, and packaging. Choose the Semantic Versioning
bump from the user-visible behavior rather than the number or size of changed
files.

For an approved release:

1. Update `CHANGELOG.md` and the version in `pyproject.toml`.
2. Run `scripts/check.sh` and commit the release changes.
3. Create and push a matching `vX.Y.Z` tag.

```console
git tag vX.Y.Z
git push origin vX.Y.Z
```

The tag-triggered release workflow verifies that the tag and project version
match, reruns tests, static checks, documentation builds, and package builds,
then publishes the wheel, source distribution, and
`aerophysics-docs-X.Y.Z.zip` to the GitHub Release. The documentation archive
opens at `aerophysics-docs-X.Y.Z/index.html` after extraction.

### Release permissions

The release workflow validates and builds with read-only repository access;
checkout credentials are not persisted. An isolated publish job downloads the
validated artifacts and is the only job with `contents: write`. It does not
check out or execute project code. Release Actions are pinned to immutable
commits; the existing Python and uv versions remain unchanged.

Publication uses `gh release create` only. An existing release causes a failure
instead of replacing its assets. Investigate partial publication failures
before retrying; the workflow does not overwrite an existing release to recover.
