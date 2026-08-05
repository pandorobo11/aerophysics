"""Reproducible verification checks for the standard-atmosphere record."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics.atmosphere import (
    MAX_GEOMETRIC_ALTITUDE,
    MIN_GEOMETRIC_ALTITUDE,
    STANDARD_GRAVITY,
    geometric_to_geopotential,
    geopotential_to_geometric,
    standard_atmosphere,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_DIRECTORY = PROJECT_ROOT / "tests/reference_data/standard_atmosphere"
OFFICIAL_CSV = REFERENCE_DIRECTORY / "official_1976.csv"
FLUIDS_CSV = REFERENCE_DIRECTORY / "fluids-1.3.1.csv"
FLUIDS_METADATA = REFERENCE_DIRECTORY / "fluids-1.3.1.json"

OFFICIAL_PROPERTIES = (
    ("temperature_K", "temperature_abs_tolerance_K", "temperature"),
    ("pressure_Pa", "pressure_abs_tolerance_Pa", "pressure"),
    ("density_kg_m3", "density_abs_tolerance_kg_m3", "density"),
    ("gravity_m_s2", "gravity_abs_tolerance_m_s2", "gravity"),
    (
        "speed_of_sound_m_s",
        "speed_of_sound_abs_tolerance_m_s",
        "speed_of_sound",
    ),
    (
        "dynamic_viscosity_Pa_s",
        "dynamic_viscosity_abs_tolerance_Pa_s",
        "dynamic_viscosity",
    ),
    (
        "kinematic_viscosity_m2_s",
        "kinematic_viscosity_abs_tolerance_m2_s",
        "kinematic_viscosity",
    ),
)

KNOWN_OFFICIAL_DEVIATIONS = {
    ("geometric", -5_000.0, "dynamic_viscosity"),
    ("geometric", 5_000.0, "dynamic_viscosity"),
    ("geometric", 10_000.0, "pressure"),
    ("geometric", 15_000.0, "pressure"),
    ("geometric", 20_000.0, "kinematic_viscosity"),
    ("geometric", 30_000.0, "kinematic_viscosity"),
    ("geometric", 35_000.0, "kinematic_viscosity"),
    ("geometric", 40_000.0, "kinematic_viscosity"),
    ("geometric", 45_000.0, "kinematic_viscosity"),
    ("geometric", 55_000.0, "dynamic_viscosity"),
    ("geometric", 60_000.0, "pressure"),
    ("geometric", 65_000.0, "pressure"),
    ("geometric", 65_000.0, "kinematic_viscosity"),
    ("geometric", 70_000.0, "kinematic_viscosity"),
    ("geometric", 75_000.0, "kinematic_viscosity"),
    ("geometric", 80_000.0, "pressure"),
    ("geometric", 80_000.0, "kinematic_viscosity"),
    ("geometric", 85_000.0, "density"),
    ("geometric", 85_000.0, "kinematic_viscosity"),
}


def _csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _official_parameters() -> list[object]:
    parameters: list[object] = []
    for row in _csv_rows(OFFICIAL_CSV):
        for value_column, tolerance_column, state_field in OFFICIAL_PROPERTIES:
            if not row[value_column]:
                continue
            coordinate = row["coordinate_type"]
            altitude = float(row["altitude_m"])
            key = (coordinate, altitude, state_field)
            marks = (
                pytest.mark.xfail(
                    strict=True,
                    reason="documented strict printed-digit deviation",
                )
                if key in KNOWN_OFFICIAL_DEVIATIONS
                else ()
            )
            parameters.append(
                pytest.param(
                    coordinate,
                    altitude,
                    state_field,
                    float(row[value_column]),
                    float(row[tolerance_column]),
                    marks=marks,
                    id=f"{coordinate}-{altitude:g}-{state_field}",
                )
            )
    return parameters


@pytest.mark.parametrize(
    ("coordinate", "altitude", "state_field", "expected", "tolerance"),
    _official_parameters(),
)
def test_public_api_matches_official_printed_cells(
    coordinate: str,
    altitude: float,
    state_field: str,
    expected: float,
    tolerance: float,
) -> None:
    geometric_altitude = (
        altitude
        if coordinate == "geometric"
        else float(geopotential_to_geometric(altitude))
    )
    actual = float(getattr(standard_atmosphere(geometric_altitude), state_field))
    guard = 1.0e-14 * max(1.0, abs(expected))
    assert abs(actual - expected) <= tolerance + guard


def test_public_api_matches_official_acceptance_criterion() -> None:
    for row in _csv_rows(OFFICIAL_CSV):
        altitude = float(row["altitude_m"])
        geometric_altitude = (
            altitude
            if row["coordinate_type"] == "geometric"
            else float(geopotential_to_geometric(altitude))
        )
        state = standard_atmosphere(geometric_altitude)
        for value_column, tolerance_column, state_field in OFFICIAL_PROPERTIES:
            if not row[value_column]:
                continue
            expected = float(row[value_column])
            printed_half_digit = float(row[tolerance_column])
            tolerance = max(2.0 * printed_half_digit, 1.0e-4 * abs(expected))
            guard = 1.0e-14 * max(1.0, abs(expected))
            actual = float(getattr(state, state_field))
            assert abs(actual - expected) <= tolerance + guard, (
                row["coordinate_type"],
                altitude,
                state_field,
            )


EXTERNAL_PROPERTIES = (
    ("temperature_K", "temperature", "absolute", 1.0e-4),
    ("pressure_Pa", "pressure", "relative", 2.0e-5),
    ("density_kg_m3", "density", "relative", 2.0e-5),
    ("gravity_m_s2", "gravity", "relative", 2.0e-5),
    ("speed_of_sound_m_s", "speed_of_sound", "relative", 2.0e-5),
    (
        "dynamic_viscosity_Pa_s",
        "dynamic_viscosity",
        "relative",
        2.0e-5,
    ),
    (
        "kinematic_viscosity_m2_s",
        "kinematic_viscosity",
        "relative",
        2.0e-5,
    ),
    (
        "thermal_conductivity_W_m_K",
        "thermal_conductivity",
        "relative",
        2.0e-5,
    ),
)


@pytest.mark.parametrize(
    ("reference_column", "state_field", "criterion", "threshold"),
    EXTERNAL_PROPERTIES,
)
def test_public_api_matches_pinned_fluids_snapshot(
    reference_column: str,
    state_field: str,
    criterion: str,
    threshold: float,
) -> None:
    rows = _csv_rows(FLUIDS_CSV)
    altitude = np.asarray(
        [float(row["geometric_altitude_m"]) for row in rows], dtype=np.float64
    )
    reference = np.asarray(
        [float(row[reference_column]) for row in rows], dtype=np.float64
    )
    actual = np.asarray(
        getattr(standard_atmosphere(altitude), state_field), dtype=np.float64
    )
    difference = np.abs(actual - reference)
    if criterion == "relative":
        difference /= np.abs(reference)
    assert np.all(difference <= threshold)


def test_fluids_snapshot_provenance_is_pinned() -> None:
    metadata = json.loads(FLUIDS_METADATA.read_text(encoding="utf-8"))
    assert metadata["source"] == "fluids.atmosphere.ATMOSPHERE_1976"
    assert metadata["version"] == "1.3.1"
    assert metadata["wheel"] == {
        "filename": "fluids-1.3.1-py3-none-any.whl",
        "hash_source": "https://pypi.org/pypi/fluids/1.3.1/json",
        "sha256": ("d9097efe57c910ac14b89a1984d5e7062ee9df1afc84e089d944fdd85404e361"),
    }
    rows = _csv_rows(FLUIDS_CSV)
    assert len(rows) == 87
    assert float(rows[0]["geometric_altitude_m"]) == 0.0
    assert float(rows[-1]["geometric_altitude_m"]) == 86_000.0


def test_full_range_monotonic_physical_behaviour() -> None:
    altitude = np.arange(-5_000.0, 86_001.0, 1.0, dtype=np.float64)
    state = standard_atmosphere(altitude)
    for field in ("pressure", "density", "gravity"):
        values = np.asarray(getattr(state, field), dtype=np.float64)
        assert np.all(np.diff(values) < 0.0), field
    viscosity = np.asarray(state.kinematic_viscosity, dtype=np.float64)
    assert np.all(np.diff(viscosity) > 0.0)


def test_temperature_gradients_follow_the_seven_layers() -> None:
    boundaries = np.asarray(
        [0.0, 11_000.0, 20_000.0, 32_000.0, 47_000.0, 51_000.0, 71_000.0, 84_852.0]
    )
    expected = np.asarray([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0])
    midpoints = (boundaries[:-1] + boundaries[1:]) / 2.0
    offset = 100.0
    lower = np.asarray(geopotential_to_geometric(midpoints - offset))
    upper = np.asarray(geopotential_to_geometric(midpoints + offset))
    lower_temperature = np.asarray(standard_atmosphere(lower).temperature)
    upper_temperature = np.asarray(standard_atmosphere(upper).temperature)
    gradient = (upper_temperature - lower_temperature) / (2.0 * offset) * 1_000.0
    assert_allclose(gradient, expected, rtol=0.0, atol=2.0e-12)


def test_layer_boundary_properties_are_continuous() -> None:
    boundaries = np.asarray(
        [11_000.0, 20_000.0, 32_000.0, 47_000.0, 51_000.0, 71_000.0]
    )
    offset = 1.0e-6
    below = standard_atmosphere(geopotential_to_geometric(boundaries - offset))
    above = standard_atmosphere(geopotential_to_geometric(boundaries + offset))
    for field in (
        "temperature",
        "pressure",
        "density",
        "speed_of_sound",
        "dynamic_viscosity",
        "kinematic_viscosity",
        "thermal_conductivity",
        "prandtl_number",
    ):
        assert_allclose(
            getattr(below, field),
            getattr(above, field),
            rtol=1.0e-9,
            atol=1.0e-12,
            err_msg=field,
        )


def test_state_equations_close_over_the_full_range() -> None:
    altitude = np.arange(-5_000.0, 86_001.0, 1.0, dtype=np.float64)
    state = standard_atmosphere(altitude)
    pressure = np.asarray(state.pressure)
    density = np.asarray(state.density)
    temperature = np.asarray(state.temperature)
    gas_constant = np.asarray(state.specific_gas_constant)
    gamma = np.asarray(state.heat_capacity_ratio)
    dynamic_viscosity = np.asarray(state.dynamic_viscosity)
    conductivity = np.asarray(state.thermal_conductivity)
    assert_allclose(pressure, density * gas_constant * temperature, rtol=1.0e-12)
    assert_allclose(
        np.asarray(state.speed_of_sound) ** 2,
        gamma * gas_constant * temperature,
        rtol=1.0e-12,
    )
    assert_allclose(
        state.kinematic_viscosity,
        dynamic_viscosity / density,
        rtol=1.0e-12,
    )
    assert_allclose(
        state.prandtl_number,
        dynamic_viscosity * np.asarray(state.cp) / conductivity,
        rtol=1.0e-12,
    )
    assert_allclose(
        np.asarray(state.cp) - np.asarray(state.cv),
        gas_constant,
        rtol=1.0e-12,
    )
    assert_allclose(
        gamma,
        np.asarray(state.cp) / np.asarray(state.cv),
        rtol=1.0e-12,
    )


def test_numerical_hydrostatic_relation_over_the_full_range() -> None:
    altitude = np.arange(-5_000.0, 86_001.0, 1.0, dtype=np.float64)
    state = standard_atmosphere(altitude)
    geopotential = np.asarray(state.geopotential_altitude)
    temperature = np.asarray(state.temperature)
    pressure = np.asarray(state.pressure)
    gas_constant = np.asarray(state.specific_gas_constant)
    derivative = np.gradient(np.log(pressure), geopotential, edge_order=2)
    expected = -STANDARD_GRAVITY / (gas_constant * temperature)
    internal_boundaries = np.asarray(
        [11_000.0, 20_000.0, 32_000.0, 47_000.0, 51_000.0, 71_000.0]
    )
    distance = np.min(np.abs(geopotential[:, np.newaxis] - internal_boundaries), axis=1)
    mask = distance > 2.0
    mask[[0, -1]] = False
    assert_allclose(derivative[mask], expected[mask], rtol=1.0e-7, atol=0.0)


def test_altitude_round_trip_and_range_endpoints_on_one_metre_grid() -> None:
    geometric = np.arange(
        MIN_GEOMETRIC_ALTITUDE,
        MAX_GEOMETRIC_ALTITUDE + 1.0,
        1.0,
        dtype=np.float64,
    )
    geopotential = geometric_to_geopotential(geometric)
    assert_allclose(
        geopotential_to_geometric(geopotential), geometric, rtol=1.0e-14, atol=1.0e-10
    )
    state = standard_atmosphere([MIN_GEOMETRIC_ALTITUDE, MAX_GEOMETRIC_ALTITUDE])
    assert_allclose(
        state.geometric_altitude,
        [MIN_GEOMETRIC_ALTITUDE, MAX_GEOMETRIC_ALTITUDE],
        rtol=0.0,
        atol=0.0,
    )
