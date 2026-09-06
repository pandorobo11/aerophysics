#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

requirements_file="$(mktemp "${TMPDIR:-/tmp}/aerophysics-audit.XXXXXX")"
trap 'rm -f "$requirements_file"' EXIT

echo "Exporting all locked dependency groups and extras"
uv export \
    --quiet \
    --all-groups \
    --all-extras \
    --no-emit-project \
    --locked \
    --output-file "$requirements_file"

echo "Auditing locked dependencies"
uv run --quiet --all-groups --all-extras --locked pip-audit \
    --requirement "$requirements_file" \
    --require-hashes \
    --disable-pip \
    --progress-spinner=off
