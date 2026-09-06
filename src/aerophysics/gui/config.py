"""Versioned calculation-configuration interchange format."""

from __future__ import annotations

import json
from collections.abc import Mapping
from typing import Any

from aerophysics.gui._config_schema import (
    _SchemaError,
    calculator_names,
    validate_calculator_payload,
)
from aerophysics.gui.units import UnitPreferences

CONFIG_SCHEMA_VERSION = 1
CALCULATORS = calculator_names()
MODES = {"single", "sweep"}


class ConfigurationError(ValueError):
    """Raised when an imported GUI configuration is invalid."""


def make_configuration(
    *,
    calculator: str,
    mode: str,
    inputs_si: Mapping[str, object],
    models: Mapping[str, object],
    units: UnitPreferences,
    sweep_si: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """Build a canonical version-one configuration."""
    try:
        configuration: dict[str, object] = {
            "schema_version": CONFIG_SCHEMA_VERSION,
            "calculator": calculator,
            "mode": mode,
            "inputs_si": dict(inputs_si),
            "models": dict(models),
            "display_units": units.to_dict(),
        }
        if sweep_si is not None:
            configuration["sweep_si"] = dict(sweep_si)
    except (KeyError, OverflowError, TypeError, ValueError) as error:
        raise ConfigurationError(
            "configuration values could not be normalized"
        ) from error
    return validate_configuration(configuration)


def validate_configuration(value: object) -> dict[str, object]:
    """Validate and normalize an imported configuration."""
    if not isinstance(value, dict):
        raise ConfigurationError("configuration must be a JSON object")
    required = {
        "schema_version",
        "calculator",
        "mode",
        "inputs_si",
        "models",
        "display_units",
    }
    allowed = required | {"sweep_si"}
    if not required <= set(value) or set(value) - allowed:
        raise ConfigurationError("configuration fields are missing or unsupported")
    if (
        not isinstance(value["schema_version"], int)
        or isinstance(value["schema_version"], bool)
        or value["schema_version"] != CONFIG_SCHEMA_VERSION
    ):
        raise ConfigurationError(
            f"unsupported schema_version: {value['schema_version']!r}"
        )
    calculator = value["calculator"]
    if not isinstance(calculator, str) or calculator not in CALCULATORS:
        raise ConfigurationError("unsupported calculator")
    mode = value["mode"]
    if not isinstance(mode, str) or mode not in MODES:
        raise ConfigurationError("mode must be single or sweep")
    try:
        units = UnitPreferences.from_dict(value["display_units"])
        inputs, models, sweep = validate_calculator_payload(
            calculator=calculator,
            mode=mode,
            inputs_si=value["inputs_si"],
            models=value["models"],
            sweep_si=value.get("sweep_si"),
            has_sweep="sweep_si" in value,
        )
    except (_SchemaError, KeyError, OverflowError, TypeError, ValueError) as error:
        raise ConfigurationError(str(error)) from error
    normalized = dict(value)
    normalized["inputs_si"] = inputs
    normalized["models"] = models
    normalized["display_units"] = units.to_dict()
    if sweep is not None:
        normalized["sweep_si"] = sweep
    else:
        normalized.pop("sweep_si", None)
    return normalized


def dump_configuration(configuration: Mapping[str, Any]) -> str:
    """Serialize a validated configuration as stable, readable JSON."""
    try:
        normalized = validate_configuration(dict(configuration))
        return json.dumps(
            normalized,
            allow_nan=False,
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
    except ConfigurationError:
        raise
    except (KeyError, OverflowError, TypeError, ValueError) as error:
        raise ConfigurationError("configuration is not JSON serializable") from error


def _reject_json_constant(value: str) -> object:
    raise ConfigurationError(f"configuration contains invalid JSON number {value}")


def _object_without_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ConfigurationError(f"configuration contains duplicate key: {key}")
        result[key] = value
    return result


def load_configuration(serialized: str | bytes) -> dict[str, object]:
    """Deserialize and validate configuration JSON."""
    try:
        value = json.loads(
            serialized,
            parse_constant=_reject_json_constant,
            object_pairs_hook=_object_without_duplicate_keys,
        )
    except ConfigurationError:
        raise
    except (json.JSONDecodeError, TypeError, UnicodeDecodeError) as error:
        raise ConfigurationError("configuration is not valid JSON") from error
    return validate_configuration(value)
