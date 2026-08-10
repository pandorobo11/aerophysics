"""Generate dilute-air transport reference values from a NIST primary source."""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
from pathlib import Path

from _verification_common import write_or_check

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / "tests/reference_data/thermophysical"
CSV_PATH = DIRECTORY / "nist-lemmon-jacobsen-2004.csv"
JSON_PATH = DIRECTORY / "nist-lemmon-jacobsen-2004.json"

TEMPERATURES = (250.0, 300.0, 400.0, 500.0, 750.0, 1000.0, 1250.0, 1500.0)
COLLISION_INTEGRAL_COEFFICIENTS = (
    0.431,
    -0.4623,
    0.08406,
    0.005341,
    -0.00331,
)
MOLAR_MASS_G_MOL = 28.9586
LENNARD_JONES_ENERGY_K = 103.3
LENNARD_JONES_SIZE_NM = 0.360
CRITICAL_TEMPERATURE_K = 132.6312
CONDUCTIVITY_COEFFICIENTS = (1.308, 1.405, -1.1, -1.036, -0.3)


def dilute_air_transport(temperature: float) -> tuple[float, float]:
    """Return Lemmon--Jacobsen dilute-air viscosity and conductivity in SI."""
    logarithm = math.log(temperature / LENNARD_JONES_ENERGY_K)
    collision_integral = math.exp(
        sum(
            coefficient * logarithm**index
            for index, coefficient in enumerate(COLLISION_INTEGRAL_COEFFICIENTS)
        )
    )
    viscosity_micro_pa_s = (
        0.0266958
        * math.sqrt(MOLAR_MASS_G_MOL * temperature)
        / (LENNARD_JONES_SIZE_NM**2 * collision_integral)
    )
    reduced_inverse_temperature = CRITICAL_TEMPERATURE_K / temperature
    n1, n2, t2, n3, t3 = CONDUCTIVITY_COEFFICIENTS
    conductivity_mw_m_k = (
        n1 * viscosity_micro_pa_s
        + n2 * reduced_inverse_temperature**t2
        + n3 * reduced_inverse_temperature**t3
    )
    return viscosity_micro_pa_s * 1.0e-6, conductivity_mw_m_k * 1.0e-3


def _verify_paper_anchors() -> None:
    anchors = (
        (100.0, 7.09559e-6, 9.35902e-3, 5.0e-12, 5.0e-9),
        (300.0, 18.5230e-6, 26.3529e-3, 5.0e-11, 5.0e-8),
    )
    for (
        temperature,
        expected_viscosity,
        expected_conductivity,
        mu_atol,
        k_atol,
    ) in anchors:
        viscosity, conductivity = dilute_air_transport(temperature)
        if abs(viscosity - expected_viscosity) > mu_atol:
            raise ValueError(f"Table V viscosity anchor failed at {temperature:g} K")
        if abs(conductivity - expected_conductivity) > k_atol:
            raise ValueError(f"Table V conductivity anchor failed at {temperature:g} K")


def _csv_content() -> str:
    stream = io.StringIO()
    fieldnames = (
        "temperature_K",
        "molar_density_mol_dm3",
        "dynamic_viscosity_Pa_s",
        "thermal_conductivity_W_m_K",
    )
    writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    for temperature in TEMPERATURES:
        viscosity, conductivity = dilute_air_transport(temperature)
        writer.writerow(
            {
                "temperature_K": f"{temperature:.1f}",
                "molar_density_mol_dm3": "0.0",
                "dynamic_viscosity_Pa_s": f"{viscosity:.15g}",
                "thermal_conductivity_W_m_K": f"{conductivity:.15g}",
            }
        )
    return stream.getvalue()


def _metadata_content() -> str:
    metadata = {
        "source": (
            "Lemmon and Jacobsen (2004), Viscosity and Thermal Conductivity "
            "Equations for Nitrogen, Oxygen, Argon, and Air"
        ),
        "doi": "10.1023/B:IJOT.0000022327.04529.F3",
        "url": "https://trc.nist.gov/refprop/FAQ/NAO.PDF",
        "pdf_sha256": (
            "985428589472f6316af725b12d5ce1aaf3fed08cd29ff175ddf13e606da5d46e"
        ),
        "role": "non-gating evaluated-reference physical-accuracy assessment",
        "state": "dilute dry air at the zero-density limit",
        "source_locations": ["Equations (2) and (5)", "Tables I, II, IV, and V"],
        "units": {
            "dynamic_viscosity": "Pa s",
            "thermal_conductivity": "W/(m K)",
        },
        "estimated_relative_uncertainty": {
            "dynamic_viscosity_above_200_K": 0.01,
            "thermal_conductivity_dilute_gas": 0.02,
        },
        "paper_anchors": [
            {
                "temperature_K": 100.0,
                "molar_density_mol_dm3": 0.0,
                "dynamic_viscosity_microPa_s": 7.09559,
                "thermal_conductivity_mW_m_K": 9.35902,
            },
            {
                "temperature_K": 300.0,
                "molar_density_mol_dm3": 0.0,
                "dynamic_viscosity_microPa_s": 18.5230,
                "thermal_conductivity_mW_m_K": 26.3529,
            },
        ],
        "generation": {
            "command": "python docs/scripts/generate_nist_transport_reference.py",
            "equation_constants": {
                "molar_mass_g_mol": MOLAR_MASS_G_MOL,
                "lennard_jones_energy_K": LENNARD_JONES_ENERGY_K,
                "lennard_jones_size_nm": LENNARD_JONES_SIZE_NM,
                "critical_temperature_K": CRITICAL_TEMPERATURE_K,
                "collision_integral_coefficients": (COLLISION_INTEGRAL_COEFFICIENTS),
                "conductivity_coefficients": CONDUCTIVITY_COEFFICIENTS,
            },
        },
    }
    return json.dumps(metadata, indent=2, sort_keys=True)


def generate(*, check: bool) -> bool:
    _verify_paper_anchors()
    current = write_or_check(CSV_PATH, _csv_content(), check=check)
    current &= write_or_check(JSON_PATH, _metadata_content(), check=check)
    return current


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if not generate(check=args.check):
        raise SystemExit("NIST transport reference assets are stale")


if __name__ == "__main__":
    main()
