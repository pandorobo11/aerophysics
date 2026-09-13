"""Regression tests for calculator-discriminated GUI configurations."""

from __future__ import annotations

import pytest

from aerophysics.gui.config import (
    ConfigurationError,
    dump_configuration,
    load_configuration,
    validate_configuration,
)
from aerophysics.gui.units import UnitPreferences


def _flight_configuration(*, mode: str = "single") -> dict[str, object]:
    configuration: dict[str, object] = {
        "schema_version": 1,
        "calculator": "flight",
        "mode": mode,
        "inputs_si": {
            "geometric_altitude": 1_000.0,
            "motion": 0.8,
            "characteristic_length": None,
        },
        "models": {"motion_basis": "mach"},
        "display_units": UnitPreferences().to_dict(),
    }
    if mode == "sweep":
        configuration["sweep_si"] = {
            "field": "altitude",
            "start": 0.0,
            "stop": 2_000.0,
            "points": 101,
        }
    return configuration


@pytest.mark.parametrize(
    ("section", "field", "value", "message"),
    (
        ("inputs_si", "motion", "0.8", "must be a number"),
        ("inputs_si", "motion", True, "must be a number"),
        ("inputs_si", "motion", float("nan"), "must be finite"),
        ("inputs_si", "motion", float("inf"), "must be finite"),
        ("inputs_si", "characteristic_length", 0.0, "greater than"),
        ("models", "motion_basis", "airspeed", "must be one of"),
    ),
)
def test_field_types_finite_values_and_enums_are_validated(
    section: str, field: str, value: object, message: str
) -> None:
    configuration = _flight_configuration()
    payload = configuration[section]
    assert isinstance(payload, dict)
    payload[field] = value

    with pytest.raises(ConfigurationError, match=message):
        validate_configuration(configuration)


@pytest.mark.parametrize("operation", ("missing", "unknown"))
def test_required_and_unknown_payload_fields_are_rejected(operation: str) -> None:
    configuration = _flight_configuration()
    payload = configuration["inputs_si"]
    assert isinstance(payload, dict)
    if operation == "missing":
        payload.pop(next(iter(payload)))
        message = "missing required fields"
    else:
        payload["unexpected"] = 1.0
        message = "unsupported fields"

    with pytest.raises(ConfigurationError, match=message):
        validate_configuration(configuration)


@pytest.mark.parametrize(
    ("update", "message"),
    (
        ({"field": "velocity"}, "field must be one of"),
        ({"start": 2_000.0}, "sweep start must be less"),
        ({"points": True}, "points must be an integer"),
        ({"points": 1}, "points must be at least"),
        ({"points": 502}, "points must be at most"),
    ),
)
def test_sweep_contract_is_validated(update: dict[str, object], message: str) -> None:
    configuration = _flight_configuration(mode="sweep")
    sweep = configuration["sweep_si"]
    assert isinstance(sweep, dict)
    sweep.update(update)

    with pytest.raises(ConfigurationError, match=message):
        validate_configuration(configuration)


def test_mode_and_sweep_payload_must_be_consistent() -> None:
    single = _flight_configuration()
    single["sweep_si"] = {
        "field": "altitude",
        "start": 0.0,
        "stop": 1_000.0,
        "points": 11,
    }
    with pytest.raises(ConfigurationError, match="does not allow sweep_si"):
        validate_configuration(single)

    sweep = _flight_configuration(mode="sweep")
    del sweep["sweep_si"]
    with pytest.raises(ConfigurationError, match="requires sweep_si"):
        validate_configuration(sweep)


@pytest.mark.parametrize("invalid_sweep", (None, []))
def test_sweep_payload_must_be_an_object(invalid_sweep: object) -> None:
    configuration = _flight_configuration(mode="sweep")
    configuration["sweep_si"] = invalid_sweep

    with pytest.raises(ConfigurationError, match="sweep_si must be an object"):
        validate_configuration(configuration)


def test_calculator_specific_sweep_point_limit_is_validated() -> None:
    configuration = {
        "schema_version": 1,
        "calculator": "conical_shock",
        "mode": "sweep",
        "inputs_si": {"upstream_mach": 2.0, "cone_half_angle": 0.1},
        "models": {},
        "display_units": UnitPreferences().to_dict(),
        "sweep_si": {
            "field": "mach",
            "start": 1.1,
            "stop": 5.0,
            "points": 202,
        },
    }

    with pytest.raises(ConfigurationError, match="points must be at most 201"):
        validate_configuration(configuration)


def test_nonstandard_json_numbers_are_rejected_on_load() -> None:
    serialized = dump_configuration(_flight_configuration()).replace("0.8", "NaN")

    with pytest.raises(ConfigurationError, match="invalid JSON number"):
        load_configuration(serialized)


def test_duplicate_json_keys_are_rejected() -> None:
    serialized = dump_configuration(_flight_configuration()).replace(
        '"motion": 0.8\n', '"motion": 0.8,\n    "motion": 0.9\n'
    )

    with pytest.raises(ConfigurationError, match="duplicate key: motion"):
        load_configuration(serialized)


def test_dump_rejects_non_json_safe_values_as_configuration_errors() -> None:
    configuration = _flight_configuration()
    inputs = configuration["inputs_si"]
    assert isinstance(inputs, dict)
    inputs["motion"] = object()

    with pytest.raises(ConfigurationError, match=r"inputs_si\.motion"):
        dump_configuration(configuration)


def test_validation_returns_a_detached_normalized_payload() -> None:
    source = _flight_configuration()
    normalized = validate_configuration(source)
    source_inputs = source["inputs_si"]
    normalized_inputs = normalized["inputs_si"]
    assert isinstance(source_inputs, dict)
    assert isinstance(normalized_inputs, dict)
    source_inputs["motion"] = 2.0
    assert normalized_inputs["motion"] == 0.8
