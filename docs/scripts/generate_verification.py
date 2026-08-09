"""Regenerate or check every committed verification artifact."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

DIRECTORY = Path(__file__).resolve().parent
GENERATORS = (
    "generate_standard_atmosphere_validation.py",
    "generate_compressible_flow_validation.py",
    "generate_thermophysical_validation.py",
    "generate_viscous_flow_validation.py",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    arguments = ["--check"] if args.check else []
    for generator in GENERATORS:
        subprocess.run(
            [sys.executable, str(DIRECTORY / generator), *arguments], check=True
        )


if __name__ == "__main__":
    main()
