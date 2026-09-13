"""Offline third-party and invariant checks for thermophysical models."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

from aerophysics.real_gas import AIR_BEATTIE_BRIDGEMAN
from aerophysics.thermochemistry import AIR_NASA7, AIR_NASA9
from aerophysics.transport import (
    AIR_BLOTTNER_VISCOSITY,
    AIR_CONDUCTIVITY,
    AIR_KEYES_VISCOSITY,
    AIR_VISCOSITY,
)

REFERENCE = Path(__file__).parent / "reference_data/thermophysical"


def _rows(name: str | Path) -> list[dict[str, str]]:
    path = Path(name)
    if not path.is_absolute():
        path = REFERENCE / path
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def test_pinned_cantera_snapshot_matches_nasa_air_models() -> None:
    gases = {"NASA7": AIR_NASA7, "NASA9": AIR_NASA9}
    for row in _rows("cantera-3.2.0.csv"):
        gas = gases[row["model"]]
        temperature = float(row["temperature_K"])
        pressure = float(row["pressure_Pa"])
        actual = {
            "cp_J_kg_K": gas.cp(temperature),
            "enthalpy_J_kg": gas.standard_enthalpy(temperature),
            "entropy_J_kg_K": gas.entropy(temperature, pressure),
            "gamma": gas.heat_capacity_ratio(temperature),
            "speed_of_sound_m_s": gas.speed_of_sound(temperature),
        }
        for name, value in actual.items():
            expected = float(row[name])
            assert abs(float(value) - expected) <= 2.0e-6 * max(abs(expected), 1.0)


def test_nist_reference_covers_dilute_air_temperature_range() -> None:
    rows = _rows("nist-lemmon-jacobsen-2004.csv")
    assert len(rows) == 8
    assert all(float(row["molar_density_mol_dm3"]) == 0.0 for row in rows)
    assert float(rows[0]["temperature_K"]) == 250.0
    assert float(rows[-1]["temperature_K"]) == 1500.0


def test_transport_models_match_primary_ussa_references() -> None:
    metadata = json.loads(
        (REFERENCE / "ussa-1976-transport.json").read_text(encoding="utf-8")
    )
    viscosity_fixture = (REFERENCE / metadata["dynamic_viscosity"]["fixture"]).resolve()
    for row in _rows(viscosity_fixture):
        if not row["dynamic_viscosity_Pa_s"]:
            continue
        expected = float(row["dynamic_viscosity_Pa_s"])
        printed_tolerance = float(row["dynamic_viscosity_abs_tolerance_Pa_s"])
        tolerance = max(2.0 * printed_tolerance, 1.0e-4 * abs(expected))
        actual = AIR_VISCOSITY.dynamic_viscosity(float(row["temperature_K"]))
        assert abs(actual - expected) <= tolerance

    conductivity = metadata["thermal_conductivity"]
    actual_conductivity = AIR_CONDUCTIVITY.thermal_conductivity(
        float(conductivity["temperature_K"])
    )
    assert abs(actual_conductivity - float(conductivity["value_W_m_K"])) <= float(
        conductivity["printed_abs_tolerance_W_m_K"]
    )


def test_transport_models_match_direct_source_equation_reproductions() -> None:
    models = {
        "Sutherland viscosity": AIR_VISCOSITY.dynamic_viscosity,
        "Keyes viscosity": AIR_KEYES_VISCOSITY.dynamic_viscosity,
        "Blottner/Wilke viscosity": AIR_BLOTTNER_VISCOSITY.dynamic_viscosity,
        "USSA conductivity": AIR_CONDUCTIVITY.thermal_conductivity,
    }
    for row in _rows("transport_source_equations.csv"):
        actual = models[row["model"]](float(row["temperature_K"]))
        assert_allclose(float(actual), float(row["value_SI"]), rtol=1.0e-12)

    metadata = json.loads(
        (REFERENCE / "transport_source_equations.json").read_text(encoding="utf-8")
    )
    assert set(metadata["models"]) == set(models)


def test_nasa_polynomial_identities_and_derivative() -> None:
    # Endpoints, interiors and both sides of the polynomial region boundary.
    temperature = np.asarray([200.1, 500.0, 999.0, 1001.0, 3000.0, 5999.9])
    step = 1.0e-2
    for gas in (AIR_NASA7, AIR_NASA9):
        cp = np.asarray(gas.cp(temperature))
        cv = np.asarray(gas.cv(temperature))
        assert_allclose(cp - cv, gas.specific_gas_constant, rtol=1.0e-12)
        assert_allclose(gas.heat_capacity_ratio(temperature), cp / cv, rtol=1.0e-12)
        derivative = (
            np.asarray(gas.standard_enthalpy(temperature + step))
            - np.asarray(gas.standard_enthalpy(temperature - step))
        ) / (2.0 * step)
        assert_allclose(derivative, cp, rtol=1.0e-7)


def test_beattie_bridgeman_density_roots_are_stable_and_close_pressure() -> None:
    for pressure in (1.0e6, 10.0e6):
        for temperature in (400.0, 800.0, 1200.0):
            density = float(AIR_BEATTIE_BRIDGEMAN.density(temperature, pressure))
            delta = density * 1.0e-5
            lower = float(AIR_BEATTIE_BRIDGEMAN.pressure(temperature, density - delta))
            upper = float(AIR_BEATTIE_BRIDGEMAN.pressure(temperature, density + delta))
            assert upper > lower
            assert_allclose(
                AIR_BEATTIE_BRIDGEMAN.pressure(temperature, density),
                pressure,
                rtol=1.0e-12,
            )
