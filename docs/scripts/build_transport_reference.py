"""Build independent source-equation values for dry-air transport models."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / "tests/reference_data/thermophysical"
CSV_PATH = DIRECTORY / "transport_source_equations.csv"
JSON_PATH = DIRECTORY / "transport_source_equations.json"

MOLE_FRACTIONS = (0.78084, 0.209476, 0.00934, 0.000314)
MOLAR_MASSES = (0.0280134, 0.0319988, 0.039948, 0.0440095)
BLOTTNER = (
    (0.0268142, 0.3177838, -11.3155513),
    (0.0449290, -0.0826158, -9.2019475),
    (-0.02201, 1.010, -13.42),
    (-0.041372, 1.3293, -15.016),
)


def _sutherland(temperature: float) -> float:
    ratio = temperature / 288.15
    return 1.7894e-5 * ratio**1.5 * (288.15 + 110.4) / (temperature + 110.4)


def _keyes(temperature: float) -> float:
    denominator = temperature + 122.1 * 10.0 ** (-5.0 / temperature)
    return 1.488e-6 * temperature**1.5 / denominator


def _blottner_wilke(temperature: float) -> float:
    fractions = tuple(value / sum(MOLE_FRACTIONS) for value in MOLE_FRACTIONS)
    logarithm = math.log(temperature)
    viscosities = tuple(
        0.1 * math.exp((a * logarithm + b) * logarithm + c) for a, b, c in BLOTTNER
    )
    mixture = 0.0
    for i, viscosity_i in enumerate(viscosities):
        denominator = 0.0
        for j, viscosity_j in enumerate(viscosities):
            phi = (
                1.0
                + math.sqrt(viscosity_i / viscosity_j)
                * (MOLAR_MASSES[j] / MOLAR_MASSES[i]) ** 0.25
            ) ** 2 / math.sqrt(8.0 * (1.0 + MOLAR_MASSES[i] / MOLAR_MASSES[j]))
            denominator += fractions[j] * phi
        mixture += fractions[i] * viscosity_i / denominator
    return mixture


def _ussa_conductivity(temperature: float) -> float:
    denominator = temperature + 245.4 * 10.0 ** (-12.0 / temperature)
    return 2.64638e-3 * temperature**1.5 / denominator


def main() -> None:
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    definitions = (
        (
            "Sutherland viscosity",
            "dynamic_viscosity_Pa_s",
            (200.0, 250.0, 288.15, 300.0, 500.0, 1000.0, 1500.0),
            _sutherland,
            "Sutherland (1893), equation and USSA sea-level normalization",
        ),
        (
            "Keyes viscosity",
            "dynamic_viscosity_Pa_s",
            (79.0, 100.0, 200.0, 300.0, 500.0, 1000.0, 1845.0),
            _keyes,
            "Bova et al. SAND2010-1168C, reproduced dry-air equation",
        ),
        (
            "Blottner/Wilke viscosity",
            "dynamic_viscosity_Pa_s",
            (1000.0, 1500.0, 2000.0, 5000.0, 10000.0, 30000.0),
            _blottner_wilke,
            "Blottner et al. SC-RR-70-754 coefficients; Wilke mixing equation",
        ),
        (
            "USSA conductivity",
            "thermal_conductivity_W_m_K",
            (200.0, 250.0, 288.15, 300.0, 500.0, 1000.0, 1500.0),
            _ussa_conductivity,
            "U.S. Standard Atmosphere 1976 corrected equation and constants",
        ),
    )
    for model, quantity, temperatures, equation, source in definitions:
        for temperature in temperatures:
            rows.append(
                {
                    "model": model,
                    "temperature_K": f"{temperature:.15g}",
                    "quantity": quantity,
                    "value_SI": f"{equation(temperature):.15g}",
                    "source": source,
                }
            )
    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "role": "independent direct evaluation of published equations and constants",
        "command": "python docs/scripts/build_transport_reference.py",
        "models": {
            "Sutherland viscosity": {
                "constants": {"mu0_Pa_s": 1.7894e-5, "T0_K": 288.15, "S_K": 110.4},
                "source": "Sutherland (1893); normalization from USSA 1976",
            },
            "Keyes viscosity": {
                "constants": {"coefficient": 1.488e-6, "A_K": 122.1, "B_K": 5.0},
                "source": "Bova et al., SAND2010-1168C",
            },
            "Blottner/Wilke viscosity": {
                "species_order": ["N2", "O2", "Ar", "CO2"],
                "blottner_coefficients": BLOTTNER,
                "mole_fractions_before_normalization": MOLE_FRACTIONS,
                "molar_masses_kg_mol": MOLAR_MASSES,
                "source": (
                    "SC-RR-70-754; Wilke (1950); Ar/CO2 audit source in references"
                ),
            },
            "USSA conductivity": {
                "constants": {"coefficient": 2.64638e-3, "A_K": 245.4, "B_K": 12.0},
                "source": "U.S. Standard Atmosphere 1976 corrected equation",
            },
        },
    }
    JSON_PATH.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
