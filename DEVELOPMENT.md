# Development

This document is the contributor workflow for `aerophysics`, including local
validation, generated documentation assets, and releases. User-facing model
documentation belongs in the English Sphinx manual under `docs/`.

## Environment

The project requires Python 3.12 or newer and the exact `uv` version declared
in `pyproject.toml` (`0.9.13`). Synchronize every locked dependency group and
extra from the repository root:

```console
uv sync --all-groups --all-extras --locked
```

Do not update the lock file as part of an ordinary environment sync.

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

The tag-triggered release workflow verifies that the tag and project version
match, reruns tests, static checks, documentation builds, and package builds,
then publishes the wheel, source distribution, and
`aerophysics-docs-X.Y.Z.zip` to the GitHub Release. The documentation archive
opens at `aerophysics-docs-X.Y.Z/index.html` after extraction.
