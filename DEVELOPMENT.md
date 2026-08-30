# Development

This document is the contributor workflow for `aerophysics`, including local
validation, generated documentation assets, and releases. User-facing model
documentation belongs in the English Sphinx manual under `docs/`.

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

Build distributions with:

```console
scripts/build-distributions.sh
```

The script exports hash-locked constraints from every locked group and extra,
then passes them to both isolated PEP 517 builds with hash enforcement. Keep
the Hatchling backend in the `build` dependency group and every custom wheel
hook dependency in `uv.lock`; otherwise the constrained build must fail
rather than resolve an unreviewed build requirement.

The release build and CI package job initially synchronize with
`--no-install-project`. Their validation commands import the checkout through
`PYTHONPATH=src` and disable `uv run`'s implicit synchronization. This prevents
an earlier editable-project build from resolving backend requirements outside
the locked, hash-constrained distribution build. Those initial synchronizations
also use `--no-build`, so a missing third-party wheel fails closed instead of
running an unconstrained source-build backend.

## Local validation

Run the complete local gate before handing back a change:

```console
scripts/check.sh
```

The gate synchronizes the locked environment, formats Python in write mode,
and then runs Ruff lint and format checks, mypy, the normal test suite, the
generated-asset checks, warning-as-error Sphinx HTML and doctest builds, and
the wheel and source-distribution builds. To run the same checks without
rewriting Python source, use:

```console
scripts/check.sh --check-only
```

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

The tag-triggered release workflow first requires a successful `CI gate` for
the exact tagged commit and verifies that the tag and project version match.
Its repository-read-only build job does not persist checkout credentials and
reruns tests, static checks, documentation builds, and package builds. Release
Actions are pinned to immutable commit SHAs. The job checksums the exact
release bundle before an isolated provenance job attests it. Only then can the
narrowly write-enabled publish job run.

Publication is one-shot: the workflow fails if a GitHub Release already exists
for the tag and never replaces published assets. Immediately before creating
the release, it dereferences either a lightweight or annotated remote tag again
and requires it to still target the commit that produced the bundle. Repository
administrators must also protect `v*` tags from force updates or deletion with
a GitHub ruleset; that closes the unavoidable interval between the final API
check and release creation. A successful release contains the wheel, source
distribution, `aerophysics-docs-X.Y.Z.zip`, and `SHA256SUMS`; GitHub also
records build-provenance attestations for the bundle.
Verify downloaded bytes with `sha256sum --check SHA256SUMS`. The documentation
archive opens at `aerophysics-docs-X.Y.Z/index.html` after extraction.
