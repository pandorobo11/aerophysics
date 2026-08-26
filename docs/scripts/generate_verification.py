"""Regenerate verification artifacts and check deterministic generated assets."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

DIRECTORY = Path(__file__).resolve().parent
GENERATORS = (
    "generate_standard_atmosphere_validation.py",
    "generate_compressible_flow_validation.py",
    "generate_nist_transport_reference.py",
    "generate_thermophysical_validation.py",
    "generate_viscous_flow_validation.py",
)
# The standard-atmosphere RST is a measured record and is regenerated before
# documentation builds; only the remaining deterministic assets are stale-checked.
DETERMINISTIC_ASSET_GENERATORS = GENERATORS[1:]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    generators = DETERMINISTIC_ASSET_GENERATORS if args.check else GENERATORS
    arguments = ["--check"] if args.check else []
    failures = []
    for generator in generators:
        result = subprocess.run(
            [sys.executable, str(DIRECTORY / generator), *arguments]
        )
        if result.returncode:
            failures.append(generator)
    if failures:
        action = "check" if args.check else "generation"
        joined = ", ".join(failures)
        raise SystemExit(f"verification asset {action} failed: {joined}")


if __name__ == "__main__":
    main()
