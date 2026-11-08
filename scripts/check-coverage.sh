#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

if [[ "$#" -ne 0 ]]; then
    echo "usage: scripts/check-coverage.sh" >&2
    exit 2
fi

if [[ -z "${COVERAGE_CORE+x}" ]] && uv run python -c 'import sys; sys.exit(sys.version_info < (3, 14))'; then
    export COVERAGE_CORE=sysmon
fi
echo "Coverage core: ${COVERAGE_CORE-default}"

uv run coverage erase
uv run pytest -m "not generated_assets and not package_artifact and not solver_verification" -n auto

echo
echo "Enforcing numerical-core branch coverage >= 95%"
uv run coverage report \
    --include="src/aerophysics/*.py" \
    --fail-under=95

echo
echo "Enforcing GUI branch coverage >= 90%"
uv run coverage report \
    --include="src/aerophysics/gui/*.py" \
    --fail-under=90
