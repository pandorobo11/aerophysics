"""Capture the pinned ``fluids`` standard-atmosphere snapshot.

Run this script in an isolated environment containing ``fluids==1.3.1``.
The resulting CSV and provenance JSON are committed so normal tests remain
offline and do not depend on the third-party package.
"""

from __future__ import annotations

import csv
import json
from importlib.metadata import version
from pathlib import Path

EXPECTED_VERSION = "1.3.1"
WHEEL_FILENAME = "fluids-1.3.1-py3-none-any.whl"
WHEEL_SHA256 = "d9097efe57c910ac14b89a1984d5e7062ee9df1afc84e089d944fdd85404e361"

PROJECT_ROOT = Path(__file__).resolve().parents[2]
REFERENCE_DIRECTORY = PROJECT_ROOT / "tests/reference_data/standard_atmosphere"
CSV_PATH = REFERENCE_DIRECTORY / "fluids-1.3.1.csv"
METADATA_PATH = REFERENCE_DIRECTORY / "fluids-1.3.1.json"


def main() -> None:
    """Write the comparison snapshot and its acquisition metadata."""
    installed_version = version("fluids")
    if installed_version != EXPECTED_VERSION:
        raise RuntimeError(
            f"fluids {EXPECTED_VERSION} is required; found {installed_version}"
        )

    from fluids.atmosphere import ATMOSPHERE_1976

    REFERENCE_DIRECTORY.mkdir(parents=True, exist_ok=True)
    fieldnames = (
        "geometric_altitude_m",
        "temperature_K",
        "pressure_Pa",
        "density_kg_m3",
        "gravity_m_s2",
        "speed_of_sound_m_s",
        "dynamic_viscosity_Pa_s",
        "kinematic_viscosity_m2_s",
        "thermal_conductivity_W_m_K",
    )
    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=fieldnames,
            lineterminator="\n",
        )
        writer.writeheader()
        for altitude in range(0, 86_001, 1_000):
            state = ATMOSPHERE_1976(float(altitude))
            writer.writerow(
                {
                    "geometric_altitude_m": altitude,
                    "temperature_K": format(state.T, ".17g"),
                    "pressure_Pa": format(state.P, ".17g"),
                    "density_kg_m3": format(state.rho, ".17g"),
                    "gravity_m_s2": format(state.g, ".17g"),
                    "speed_of_sound_m_s": format(state.v_sonic, ".17g"),
                    "dynamic_viscosity_Pa_s": format(state.mu, ".17g"),
                    "kinematic_viscosity_m2_s": format(state.mu / state.rho, ".17g"),
                    "thermal_conductivity_W_m_K": format(state.k, ".17g"),
                }
            )

    metadata = {
        "source": "fluids.atmosphere.ATMOSPHERE_1976",
        "version": EXPECTED_VERSION,
        "project_url": "https://pypi.org/project/fluids/1.3.1/",
        "documentation_url": ("https://fluids.readthedocs.io/fluids.atmosphere.html"),
        "capture_command": (
            "uv run --isolated --with fluids==1.3.1 python "
            "docs/scripts/capture_fluids_reference.py"
        ),
        "altitude_coordinate": "geometric altitude in metres",
        "altitude_grid": "0 through 86000 metres inclusive at 1000 metre steps",
        "wheel": {
            "filename": WHEEL_FILENAME,
            "sha256": WHEEL_SHA256,
            "hash_source": "https://pypi.org/pypi/fluids/1.3.1/json",
        },
        "variable_mapping": {
            "T": "temperature_K",
            "P": "pressure_Pa",
            "rho": "density_kg_m3",
            "g": "gravity_m_s2",
            "v_sonic": "speed_of_sound_m_s",
            "mu": "dynamic_viscosity_Pa_s",
            "mu / rho": "kinematic_viscosity_m2_s",
            "k": "thermal_conductivity_W_m_K",
        },
        "role": ("independent software cross-check, not experimental ground truth"),
    }
    METADATA_PATH.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
