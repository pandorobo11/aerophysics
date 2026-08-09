"""Capture an offline Cantera 3.2.0 frozen-air thermochemistry snapshot."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import cantera as ct

from aerophysics._nasa_data import (  # type: ignore[import-not-found]
    DRY_AIR_MOLE_FRACTIONS,
    NASA9_DATA,
)

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / "tests/reference_data/thermophysical"
CSV_PATH = DIRECTORY / "cantera-3.2.0.csv"
JSON_PATH = DIRECTORY / "cantera-3.2.0.json"
TEMPERATURES = (
    200.0,
    250.0,
    300.0,
    500.0,
    999.0,
    1000.0,
    1001.0,
    1500.0,
    2000.0,
    3000.0,
    4000.0,
    5000.0,
    6000.0,
)
PRESSURE = 101_325.0
MOLE_FRACTIONS = {
    name: value / sum(DRY_AIR_MOLE_FRACTIONS.values())
    for name, value in DRY_AIR_MOLE_FRACTIONS.items()
}


def _nasa9_species() -> list[ct.Species]:
    compositions = {
        "N2": {"N": 2},
        "O2": {"O": 2},
        "Ar": {"Ar": 1},
        "CO2": {"C": 1, "O": 2},
    }
    species: list[ct.Species] = []
    for name, composition in compositions.items():
        _, ranges, rows = NASA9_DATA[name]
        coefficients: list[float] = [float(len(rows))]
        for lower, upper, row in zip(ranges[:-1], ranges[1:], rows, strict=True):
            coefficients.extend((lower, upper, *row))
        item = ct.Species(name, composition)
        item.thermo = ct.Nasa9PolyMultiTempRegion(
            ranges[0], ranges[-1], 101_325.0, coefficients
        )
        species.append(item)
    return species


def _phase(model: str) -> ct.Solution:
    if model == "NASA7":
        wanted = set(MOLE_FRACTIONS)
        species = [
            item
            for item in ct.Species.list_from_file("nasa_gas.yaml")
            if item.name in wanted
        ]
    else:
        species = _nasa9_species()
    return ct.Solution(thermo="ideal-gas", species=species)


def main() -> None:
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    fixed_molar_mass = sum(
        MOLE_FRACTIONS[name] * NASA9_DATA[name][0] for name in MOLE_FRACTIONS
    )
    for model in ("NASA7", "NASA9"):
        phase = _phase(model)
        phase.X = MOLE_FRACTIONS
        for temperature in TEMPERATURES:
            phase.TP = temperature, PRESSURE
            cp = phase.cp_mole / 1000.0 / fixed_molar_mass
            cv = cp - ct.gas_constant / 1000.0 / fixed_molar_mass
            gamma = cp / cv
            gas_constant = ct.gas_constant / 1000.0 / fixed_molar_mass
            enthalpy = phase.enthalpy_mole / 1000.0 / fixed_molar_mass
            entropy = (
                phase.entropy_mole / 1000.0 / fixed_molar_mass
                - gas_constant * math.log(PRESSURE / 100_000.0)
            )
            sound_speed = math.sqrt(gamma * gas_constant * temperature)
            rows.append(
                {
                    "model": model,
                    "temperature_K": f"{temperature:.1f}",
                    "pressure_Pa": f"{PRESSURE:.1f}",
                    "cp_J_kg_K": f"{cp:.15g}",
                    "enthalpy_J_kg": f"{enthalpy:.15g}",
                    "entropy_J_kg_K": f"{entropy:.15g}",
                    "gamma": f"{gamma:.15g}",
                    "speed_of_sound_m_s": f"{sound_speed:.15g}",
                }
            )
    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "source": "Cantera ideal-gas thermodynamic evaluator",
        "version": ct.__version__,
        "composition": MOLE_FRACTIONS,
        "coefficient_mapping": {
            "NASA7": "Cantera nasa_gas.yaml; coefficients traced to NASA TM-4513",
            "NASA9": (
                "aerophysics pinned CEA coefficients evaluated by Cantera "
                "Nasa9PolyMultiTempRegion"
            ),
        },
        "entropy_reference": (
            "Cantera values normalized from 1 atm to the NASA 1 bar standard pressure"
        ),
        "command": (
            "uv run --isolated --with cantera==3.2.0 python "
            "docs/scripts/capture_cantera_reference.py"
        ),
        "wheel": {
            "filename": "cantera-3.2.0-cp312-cp312-macosx_11_0_arm64.whl",
            "hash_source": "https://pypi.org/pypi/cantera/3.2.0/json",
            "sha256": (
                "d58dd40112741423a4b9b95fbbd7789575250af1fa13c77d60fab47f472f694d"
            ),
        },
    }
    JSON_PATH.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
