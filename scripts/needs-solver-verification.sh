#!/usr/bin/env bash
# Keep this compatible with macOS bash 3.2 as well as Linux bash.
set -Eeuo pipefail

PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"
base="${SOLVER_VERIFICATION_BASE-origin/main}"

answer() {
    printf '%s\n' "$1"
    printf '%s\n' "$2" >&2
    exit 0
}

git rev-parse --verify --quiet "${base}^{commit}" >/dev/null ||
    answer yes "base ref unavailable: $base (fail safe)"
git merge-base "$base" HEAD >/dev/null 2>&1 ||
    answer yes "no merge base with $base (fail safe)"

relevant_path() {
    case "$1" in
        __SOLVER_VERIFICATION_GIT_FAILED__)
            answer yes "cannot inspect changes relative to $base (fail safe)" ;;
        src/aerophysics/gui/*) ;;
        src/aerophysics/*|tests/test_marching_boundary_layer*.py|tests/test_turbulent_sa_reference.py|docs/scripts/*|tests/reference_data/*|pyproject.toml|uv.lock|scripts/check*.sh|scripts/needs-solver-verification.sh|.github/workflows/*)
            local display
            printf -v display '%q' "$1"
            answer yes "solver-relevant change: $display"
            ;;
    esac
}

# NUL delimiters preserve spaces, tabs, and newlines in filenames. Disable
# rename detection to check both old and new paths, and expand untracked dirs.
# A sentinel propagates Git errors out of process substitutions (fail safe).
while IFS= read -r -d '' path; do
    relevant_path "$path"
done < <(git diff --name-only --no-renames -z "$base...HEAD" 2>/dev/null ||
    printf '__SOLVER_VERIFICATION_GIT_FAILED__\0')
while IFS= read -r -d '' entry; do
    relevant_path "${entry:3}"
done < <(git status --porcelain=v1 --no-renames --untracked-files=all -z 2>/dev/null ||
    printf '?? __SOLVER_VERIFICATION_GIT_FAILED__\0')

answer no "no solver-relevant changes relative to $base or in the working tree"
