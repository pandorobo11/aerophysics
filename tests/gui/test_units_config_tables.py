"""Tests for GUI unit, configuration, and table helpers."""

import json

import numpy as np
import pytest

from aerophysics.gui.adapters import Row
from aerophysics.gui.config import (
    CONFIG_SCHEMA_VERSION,
    ConfigurationError,
    dump_configuration,
    load_configuration,
    make_configuration,
    validate_configuration,
)
from aerophysics.gui.tables import columns_for, display_rows, rows_to_csv
from aerophysics.gui.units import (
    UnitPreferences,
    from_si,
    selected_unit,
    to_si,
)


@pytest.mark.parametrize(
    ("kind", "unit", "values"),
    [
        ("length", "ft", [0.0, 3.0, 10.0]),
        ("speed", "kt", [0.0, 100.0]),
        ("pressure", "psi", [0.0, 14.7]),
        ("density", "slug/ft³", [0.0, 0.00237]),
        ("angle", "deg", [0.0, 10.0, 90.0]),
        ("temperature", "°F", [-40.0, 32.0, 100.0]),
    ],
)
def test_display_unit_round_trip(kind: str, unit: str, values: list[float]) -> None:
    si = to_si(values, kind, unit)  # type: ignore[arg-type]
    assert np.asarray(from_si(si, kind, unit)) == pytest.approx(values)  # type: ignore[arg-type]


def test_scalar_units_and_preferences_validation() -> None:
    assert to_si(1.0, "length", "m") == 1.0
    preferences = UnitPreferences(
        length="ft",
        speed="kt",
        pressure="psi",
        temperature="°F",
        density="slug/ft³",
        angle="rad",
    )
    assert selected_unit("speed", preferences) == "kt"
    assert UnitPreferences.from_dict(preferences.to_dict()) == preferences
    with pytest.raises(ValueError, match="unsupported length"):
        UnitPreferences(length="yard")
    with pytest.raises(ValueError, match="units must"):
        UnitPreferences.from_dict([])
    with pytest.raises(ValueError, match="missing"):
        UnitPreferences.from_dict({"length": "m"})
    with pytest.raises(ValueError, match="finite"):
        to_si(np.nan, "length", "m")
    with pytest.raises(ValueError, match="unsupported"):
        to_si(1.0, "pressure", "bar")


def test_configuration_round_trip() -> None:
    preferences = UnitPreferences(length="ft", angle="deg")
    configuration = make_configuration(
        calculator="flight",
        mode="sweep",
        inputs_si={
            "geometric_altitude": 1000.0,
            "motion": 0.8,
            "characteristic_length": None,
        },
        models={"motion_basis": "mach"},
        units=preferences,
        sweep_si={"field": "altitude", "start": 0.0, "stop": 2000.0, "points": 3},
    )
    serialized = dump_configuration(configuration)
    restored = load_configuration(serialized)
    assert restored == configuration
    assert restored["schema_version"] == CONFIG_SCHEMA_VERSION
    assert json.loads(serialized)["display_units"]["length"] == "ft"


def test_configuration_round_trip_with_embedded_profile_arrays() -> None:
    configuration = make_configuration(
        calculator="protrusion_drag",
        mode="single",
        inputs_si={
            "profile_height": [0.0, 0.01],
            "profile_velocity": [0.0, 100.0],
            "profile_density": [1.2, 1.0],
        },
        models={"profile_source": "csv"},
        units=UnitPreferences(),
    )
    assert load_configuration(dump_configuration(configuration)) == configuration


@pytest.mark.parametrize(
    "value",
    [
        [],
        {},
        {
            "schema_version": 99,
            "calculator": "flight",
            "mode": "single",
            "inputs_si": {},
            "models": {},
            "display_units": UnitPreferences().to_dict(),
        },
    ],
)
def test_configuration_rejects_invalid_objects(value: object) -> None:
    with pytest.raises(ConfigurationError):
        validate_configuration(value)


def test_configuration_rejects_invalid_json_and_fields() -> None:
    with pytest.raises(ConfigurationError, match="valid JSON"):
        load_configuration("{")
    base = {
        "schema_version": 1,
        "calculator": "unknown",
        "mode": "single",
        "inputs_si": {},
        "models": {},
        "display_units": UnitPreferences().to_dict(),
    }
    with pytest.raises(ConfigurationError, match="calculator"):
        validate_configuration(base)
    base["calculator"] = "flight"
    base["mode"] = "other"
    with pytest.raises(ConfigurationError, match="mode"):
        validate_configuration(base)
    base["mode"] = "sweep"
    with pytest.raises(ConfigurationError, match="sweep_si"):
        validate_configuration(base)


def test_display_table_and_csv() -> None:
    rows = (
        {
            "upstream_mach": 2.0,
            "deflection_angle": np.deg2rad(10.0),
            "maximum_deflection_angle": np.deg2rad(23.0),
            "shock_angle": np.deg2rad(39.0),
            "downstream_mach": 1.64,
            "static_pressure_ratio": 1.7,
            "static_density_ratio": 1.45,
            "static_temperature_ratio": 1.17,
            "total_pressure_ratio": 0.98,
            "status": "ok",
            "message": "",
        },
    )
    table = display_rows("oblique_shock", rows, UnitPreferences())
    assert table[0]["偏向角 θ [deg]"] == pytest.approx(10.0)
    csv_text = rows_to_csv(table)
    assert csv_text.startswith("\ufeff")
    assert "偏向角 θ [deg]" in csv_text
    with pytest.raises(ValueError, match="empty"):
        rows_to_csv([])
    with pytest.raises(ValueError, match="unsupported"):
        columns_for("missing")


def test_conical_shock_table_converts_angles() -> None:
    rows = (
        {
            "upstream_mach": 2.0,
            "cone_half_angle": np.deg2rad(10.0),
            "maximum_cone_half_angle": np.deg2rad(40.0),
            "shock_angle": np.deg2rad(31.0),
            "post_shock_mach": 1.95,
            "surface_mach": 1.83,
            "surface_pressure_ratio": 1.29,
            "surface_density_ratio": 1.20,
            "surface_temperature_ratio": 1.08,
            "total_pressure_ratio": 0.999,
            "status": "ok",
            "message": "",
        },
    )
    table = display_rows("conical_shock", rows, UnitPreferences())
    assert table[0]["円錐半頂角 θc [deg]"] == pytest.approx(10.0)
    assert table[0]["衝撃波角 β [deg]"] == pytest.approx(31.0)


@pytest.mark.parametrize(
    "calculator",
    [
        "isentropic",
        "normal_shock",
        "expansion",
        "boundary_layer_profile",
        "conical_shock",
        "protrusion_drag",
        "thermochemistry",
        "viscosity",
    ],
)
def test_additional_calculator_tables_and_config(calculator: str) -> None:
    assert columns_for(calculator)
    configuration = make_configuration(
        calculator=calculator,
        mode="single",
        inputs_si={},
        models={},
        units=UnitPreferences(),
    )
    assert configuration["calculator"] == calculator


def test_viscosity_table_converts_temperature_and_preserves_fixed_units() -> None:
    rows: tuple[Row, ...] = (
        {
            "model": "Keyes",
            "temperature": 300.0,
            "dynamic_viscosity": 1.8519327e-5,
            "relative_difference": 0.320,
            "status": "ok",
            "message": "",
        },
    )
    table = display_rows("viscosity", rows, UnitPreferences(temperature="°F"))
    assert table[0]["温度 T [°F]"] == pytest.approx(80.33)
    assert table[0]["粘性係数 μ [Pa·s]"] == pytest.approx(1.8519327e-5)
    assert table[0]["Sutherland基準相対差 [%]"] == pytest.approx(0.320)
