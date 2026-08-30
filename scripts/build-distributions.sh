#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

if [[ "$#" -ne 0 ]]; then
    echo "usage: scripts/build-distributions.sh" >&2
    exit 2
fi

build_constraints="$(mktemp "${TMPDIR:-/tmp}/aerophysics-build.XXXXXX")"
trap 'rm -f "$build_constraints"' EXIT

echo "Exporting hash-locked build constraints"
uv export \
    --quiet \
    --all-groups \
    --all-extras \
    --no-emit-project \
    --locked \
    --output-file "$build_constraints"

echo "Building wheel from sdist with constrained isolated dependencies"
env -u PYTHONPATH uv build \
    --clear \
    --build-constraints "$build_constraints" \
    --require-hashes
