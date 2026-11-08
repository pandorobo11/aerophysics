#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

check_only=false
full=false
for argument in "$@"; do
    case "$argument" in
        --check-only) check_only=true ;;
        --full) full=true ;;
        *)
            echo "usage: scripts/check.sh [--check-only] [--full]" >&2
            exit 2
            ;;
    esac
done

if [[ "$check_only" == true ]]; then
    FORMAT_COMMAND=(uv run ruff format --check .)
else
    FORMAT_COMMAND=(uv run ruff format .)
fi

status=0

run_step() {
    local label="$1"
    shift
    printf '\n==> %s\n' "$label"
    if "$@"; then
        step_succeeded=true
        return 0
    fi
    step_succeeded=false
    status=1
    printf 'FAILED: %s\n' "$label" >&2
}

run_step "Sync locked dependencies" uv sync --all-groups --all-extras --locked
run_step "Format Python files" "${FORMAT_COMMAND[@]}"
run_step "Ruff lint" uv run ruff check .
run_step "Verify formatting" uv run ruff format --check .
run_step "Mypy" uv run mypy
run_step "Normal test suite and coverage gates" bash scripts/check-coverage.sh
# Only skip the currency test after successful regeneration in this invocation.
unset AEROPHYSICS_SKIP_ASSET_CURRENCY
if [[ "$check_only" == false ]]; then
    run_step "Generate verification records for documentation" \
        uv run python docs/scripts/generate_verification.py
    if [[ "$step_succeeded" == true ]]; then
        export AEROPHYSICS_SKIP_ASSET_CURRENCY=1
    fi
fi
run_step "Generated and verification assets" bash scripts/check-generated.sh
unset AEROPHYSICS_SKIP_ASSET_CURRENCY

if [[ "$full" == true ]]; then
    solver_decision=yes
    solver_reason="forced by --full"
else
    solver_check="$(bash scripts/needs-verification.sh solver 2>&1)"
    solver_decision="${solver_check%%$'\n'*}"
    solver_reason="${solver_check#*$'\n'}"
fi
if [[ "$solver_decision" == yes ]]; then
    printf '\nSolver verification: %s\n' "$solver_reason"
    run_step "Solver verification" uv run pytest -m solver_verification --no-cov -n auto
else
    printf '\nSKIPPED: Solver verification (%s)\n' "$solver_reason"
fi
run_step "Build HTML documentation" \
    uv run sphinx-build -E -a -W --keep-going -b html docs docs/_build/html
run_step "Run documentation doctests" \
    uv run sphinx-build -E -a -W --keep-going -b doctest docs docs/_build/doctest
if [[ "$full" == true ]]; then
    package_decision=yes
    package_reason="forced by --full"
else
    package_check="$(bash scripts/needs-verification.sh package 2>&1)"
    package_decision="${package_check%%$'\n'*}"
    package_reason="${package_check#*$'\n'}"
fi
if [[ "$package_decision" == yes ]]; then
    printf '\nPackage verification: %s\n' "$package_reason"
    run_step "Build Python distributions" uv run python -m build
    run_step "Verify clean wheel and sdist installations" \
        env AEROPHYSICS_TEST_ARTIFACT_DIRECTORY=dist \
        uv run pytest -m package_artifact --no-cov
else
    printf '\nSKIPPED: Build Python distributions (%s)\n' "$package_reason"
    printf 'SKIPPED: Verify clean wheel and sdist installations (%s)\n' "$package_reason"
fi

exit "$status"
