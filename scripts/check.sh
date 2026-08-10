#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

if [[ "${1:-}" == "--check-only" ]]; then
    FORMAT_COMMAND=(uv run ruff format --check .)
    shift
else
    FORMAT_COMMAND=(uv run ruff format .)
fi

if [[ "$#" -ne 0 ]]; then
    echo "usage: scripts/check.sh [--check-only]" >&2
    exit 2
fi

status=0

run_step() {
    local label="$1"
    shift
    printf '\n==> %s\n' "$label"
    if "$@"; then
        return 0
    fi
    status=1
    printf 'FAILED: %s\n' "$label" >&2
}

run_step "Sync locked dependencies" uv sync --all-groups --all-extras --locked
run_step "Format Python files" "${FORMAT_COMMAND[@]}"
run_step "Ruff lint" uv run ruff check .
run_step "Verify formatting" uv run ruff format --check .
run_step "Mypy" uv run mypy
run_step "Normal test suite" uv run pytest -m "not generated_assets"
run_step "Generated and verification assets" bash scripts/check-generated.sh
run_step "Build HTML documentation" \
    uv run sphinx-build -W --keep-going -b html docs docs/_build/html
run_step "Run documentation doctests" \
    uv run sphinx-build -W --keep-going -b doctest docs docs/_build/doctest
run_step "Build Python distributions" uv run python -m build

exit "$status"
