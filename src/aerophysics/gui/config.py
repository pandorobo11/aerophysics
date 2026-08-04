"""Versioned calculation-configuration interchange format."""

from __future__ import annotations

import json
from collections.abc import Mapping
from typing import Any

from aerophysics.gui.units import UnitPreferences

CONFIG_SCHEMA_VERSION = 1
CALCULATORS = {
    "boundary_layer",
    "boundary_layer_profile",
    "conical_shock",
    "expansion",
    "flight",
    "isentropic",
    "normal_shock",
    "oblique_shock",
    "protrusion_drag",
    "thermochemistry",
}
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
    if value["schema_version"] != CONFIG_SCHEMA_VERSION:
        raise ConfigurationError(
            f"unsupported schema_version: {value['schema_version']!r}"
        )
    if value["calculator"] not in CALCULATORS:
        raise ConfigurationError("unsupported calculator")
    if value["mode"] not in MODES:
        raise ConfigurationError("mode must be single or sweep")
    if not isinstance(value["inputs_si"], dict) or not isinstance(
        value["models"], dict
    ):
        raise ConfigurationError("inputs_si and models must be objects")
    if value["mode"] == "sweep" and not isinstance(value.get("sweep_si"), dict):
        raise ConfigurationError("sweep mode requires sweep_si")
    try:
        units = UnitPreferences.from_dict(value["display_units"])
    except ValueError as error:
        raise ConfigurationError(str(error)) from error
    normalized = dict(value)
    normalized["inputs_si"] = dict(value["inputs_si"])
    normalized["models"] = dict(value["models"])
    normalized["display_units"] = units.to_dict()
    if "sweep_si" in value:
        normalized["sweep_si"] = dict(value["sweep_si"])
    return normalized


def dump_configuration(configuration: Mapping[str, Any]) -> str:
    """Serialize a validated configuration as stable, readable JSON."""
    normalized = validate_configuration(dict(configuration))
    return json.dumps(normalized, ensure_ascii=False, indent=2, sort_keys=True)


def load_configuration(serialized: str | bytes) -> dict[str, object]:
    """Deserialize and validate configuration JSON."""
    try:
        value = json.loads(serialized)
    except (json.JSONDecodeError, UnicodeDecodeError) as error:
        raise ConfigurationError("configuration is not valid JSON") from error
    return validate_configuration(value)
