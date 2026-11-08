#!/usr/bin/env bash
# Keep this compatible with macOS bash 3.2 as well as Linux bash.
set -Eeuo pipefail

if [[ "$#" -ne 1 ]]; then
    echo "usage: scripts/needs-verification.sh <solver|package>" >&2
    exit 2
fi
case "$1" in
    solver|package) tier="$1" ;;
    *)
        echo "usage: scripts/needs-verification.sh <solver|package>" >&2
        exit 2
        ;;
esac

PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"
# Retain the former solver-specific override for existing local callers.
base="${VERIFICATION_BASE-${SOLVER_VERIFICATION_BASE-origin/main}}"

answer() {
    printf '%s\n' "$1"
    printf '%s\n' "$2" >&2
    exit 0
}

changed_path() {
    local display
    printf -v display '%q' "$1" "$1"
    answer yes "$tier-relevant $2: $display"
}

git rev-parse --verify --quiet "${base}^{commit}" >/dev/null 2>&1 ||
    answer yes "base ref unavailable: $base (fail safe)"
git merge-base "$base" HEAD >/dev/null 2>&1 ||
    answer yes "no merge base with $base (fail safe)"

relevant_path() {
    local path="$1" change="$2"
    if [[ "$path" == __VERIFICATION_GIT_FAILED__ ]]; then
        answer yes "cannot inspect changes relative to $base (fail safe)"
    fi
    if [[ "$tier" == solver ]]; then
        case "$path" in
            src/aerophysics/gui/*) ;;
            src/aerophysics/*|tests/test_marching_boundary_layer*.py|tests/test_turbulent_sa_reference.py|tests/test_characteristics.py|docs/scripts/*|tests/reference_data/*|pyproject.toml|uv.lock|scripts/check*.sh|scripts/needs-verification.sh|.github/workflows/*)
                changed_path "$path" change ;;
        esac
    else
        case "$path" in
            pyproject.toml|uv.lock|hatch_build.py|src/aerophysics/__init__.py|src/aerophysics/gui/launcher.py|src/aerophysics/gui/documentation.py|tests/test_package.py|docs/conf.py|README*|LICENSE*|scripts/check*.sh|scripts/needs-verification.sh|.github/workflows/*)
                changed_path "$path" change ;;
            # The wheel hook copies these inputs into an isolated tree and
            # executes the generators there before bundling HTML documentation.
            docs/scripts/*|tests/reference_data/*)
                changed_path "$path" "wheel documentation input change" ;;
        esac
        case "$change" in
            *A*|*D*|*R*|'??')
                case "$path" in
                    src/*|docs/*) changed_path "$path" "file addition/deletion/rename" ;;
                esac
                ;;
        esac
    fi
}

# NUL delimiters preserve spaces, tabs, and newlines in filenames. Renames
# become deletion + addition so both paths are inspected, including moves
# across src/docs boundaries. Expand untracked directories into file paths.
# Sentinels propagate Git errors out of process substitutions (fail safe).
while IFS= read -r -d '' change; do
    IFS= read -r -d '' path || answer yes "incomplete Git diff (fail safe)"
    relevant_path "$path" "$change"
done < <(git diff --name-status --no-renames -z "$base...HEAD" 2>/dev/null ||
    printf '!\0__VERIFICATION_GIT_FAILED__\0')
while IFS= read -r -d '' entry; do
    relevant_path "${entry:3}" "${entry:0:2}"
done < <(git status --porcelain=v1 --no-renames --untracked-files=all -z 2>/dev/null ||
    printf '?? __VERIFICATION_GIT_FAILED__\0')

answer no "no $tier-relevant changes relative to $base or in the working tree"
