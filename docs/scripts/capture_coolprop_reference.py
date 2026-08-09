"""Capture an observation-only CoolProp 8.0.0 air-transport snapshot."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import CoolProp
from CoolProp.CoolProp import PropsSI

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / "tests/reference_data/thermophysical"
CSV_PATH = DIRECTORY / "coolprop-8.0.0.csv"
JSON_PATH = DIRECTORY / "coolprop-8.0.0.json"
TEMPERATURES = (250.0, 300.0, 400.0, 500.0, 750.0, 1000.0, 1250.0, 1500.0)
PRESSURE = 101_325.0


def main() -> None:
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    rows = []
    for temperature in TEMPERATURES:
        viscosity = PropsSI("V", "T", temperature, "P", PRESSURE, "Air")
        conductivity = PropsSI("L", "T", temperature, "P", PRESSURE, "Air")
        rows.append(
            {
                "temperature_K": f"{temperature:.1f}",
                "pressure_Pa": f"{PRESSURE:.1f}",
                "dynamic_viscosity_Pa_s": f"{viscosity:.15g}",
                "thermal_conductivity_W_m_K": f"{conductivity:.15g}",
            }
        )
    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "source": "CoolProp Air",
        "version": CoolProp.__version__,
        "comparison_role": (
            "observation only; CoolProp and aerophysics use different "
            "transport correlations"
        ),
        "command": (
            "uv run --isolated --with CoolProp==8.0.0 python "
            "docs/scripts/capture_coolprop_reference.py"
        ),
        "wheel": {
            "filename": "coolprop-8.0.0-cp312-abi3-macosx_11_0_arm64.whl",
            "hash_source": "https://pypi.org/pypi/CoolProp/8.0.0/json",
            "sha256": (
                "46de36f51330ce8c7f8384a7360e1058739d20604736158aac6bd2abf6ae26bd"
            ),
        },
    }
    JSON_PATH.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
