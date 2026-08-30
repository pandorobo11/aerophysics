"""Calculator-specific validation for versioned GUI configurations."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

type _FieldKind = Literal["bool", "int", "number", "number_list", "string"]


class _SchemaError(ValueError):
    """Internal calculator-schema validation failure."""


@dataclass(frozen=True, slots=True)
class _FieldRule:
    kind: _FieldKind
    required: bool = True
    nullable: bool = False
    choices: frozenset[str] | None = None
    minimum: float | None = None
    maximum: float | None = None
    exclusive_minimum: bool = False


@dataclass(frozen=True, slots=True)
class _CalculatorSchema:
    inputs: dict[str, _FieldRule]
    models: dict[str, _FieldRule]
    sweep_variables: dict[str, _FieldRule] | None
    sweep_extras: dict[str, _FieldRule] | None = None
    maximum_sweep_points: int = 501


def _number(
    *,
    nullable: bool = False,
    minimum: float | None = None,
    maximum: float | None = None,
    exclusive_minimum: bool = False,
) -> _FieldRule:
    return _FieldRule(
        "number",
        nullable=nullable,
        minimum=minimum,
        maximum=maximum,
        exclusive_minimum=exclusive_minimum,
    )


def _integer(*, minimum: int, maximum: int) -> _FieldRule:
    return _FieldRule("int", minimum=float(minimum), maximum=float(maximum))


def _choice(*values: str, required: bool = True) -> _FieldRule:
    return _FieldRule("string", required=required, choices=frozenset(values))


def _nullable_number_list() -> _FieldRule:
    return _FieldRule("number_list", nullable=True)


_NUMBER = _number()
_NON_NEGATIVE = _number(minimum=0.0)
_POSITIVE = _number(minimum=0.0, exclusive_minimum=True)
_NULLABLE_POSITIVE = _number(nullable=True, minimum=0.0, exclusive_minimum=True)
_BOOL = _FieldRule("bool")
_SUPERSONIC = _number(minimum=1.0, exclusive_minimum=True)


_SCHEMAS: dict[str, _CalculatorSchema] = {
    "flight": _CalculatorSchema(
        inputs={
            "geometric_altitude": _number(minimum=-5000.0, maximum=86_000.0),
            "motion": _NON_NEGATIVE,
            "characteristic_length": _NULLABLE_POSITIVE,
        },
        models={"motion_basis": _choice("mach", "velocity")},
        sweep_variables={
            "altitude": _number(minimum=-5000.0, maximum=86_000.0),
            "motion": _NON_NEGATIVE,
        },
    ),
    "oblique_shock": _CalculatorSchema(
        inputs={
            "upstream_mach": _SUPERSONIC,
            "deflection_angle": _NON_NEGATIVE,
        },
        models={"branch": _choice("weak", "strong")},
        sweep_variables={
            "deflection": _NON_NEGATIVE,
            "mach": _SUPERSONIC,
        },
    ),
    "conical_shock": _CalculatorSchema(
        inputs={
            "upstream_mach": _SUPERSONIC,
            "cone_half_angle": _NON_NEGATIVE,
        },
        models={},
        sweep_variables={
            "cone_half_angle": _NON_NEGATIVE,
            "mach": _SUPERSONIC,
        },
        maximum_sweep_points=201,
    ),
    "boundary_layer": _CalculatorSchema(
        inputs={
            "distance": _POSITIVE,
            "edge_velocity": _POSITIVE,
            "edge_density": _POSITIVE,
            "edge_dynamic_viscosity": _POSITIVE,
            "transition_reynolds": _NULLABLE_POSITIVE,
            "mach": _number(nullable=True, minimum=0.0),
            "edge_temperature": _POSITIVE,
            "wall_temperature": _NULLABLE_POSITIVE,
        },
        models={
            "source": _choice("manual", "flight"),
            "regime": _choice("laminar", "turbulent", "transitional"),
            "turbulent_correlation": _choice("power_law", "schlichting"),
            "compressibility_correction": _choice("none", "eckert", "van_driest_ii"),
        },
        sweep_variables={"distance": _POSITIVE},
        sweep_extras={"logarithmic": _BOOL},
    ),
    "isentropic": _CalculatorSchema(
        inputs={
            "input_value": _NON_NEGATIVE,
            "total_pressure": _NULLABLE_POSITIVE,
            "total_temperature": _NULLABLE_POSITIVE,
        },
        models={
            "input_basis": _choice(
                "mach",
                "temperature_ratio",
                "pressure_ratio",
                "density_ratio",
                "area_ratio",
            ),
            "branch": _choice("subsonic", "supersonic"),
            "gas_model": _choice(
                "AIR",
                "NASA7",
                "NASA9",
                "HARMONIC_OSCILLATOR",
                "BEATTIE_BRIDGEMAN",
                required=False,
            ),
            "with_mass_flux": _BOOL,
            "allow_extrapolation": _FieldRule("bool", required=False),
        },
        sweep_variables={
            "input_value": _NUMBER,
        },
    ),
    "normal_shock": _CalculatorSchema(
        inputs={"upstream_mach": _number(minimum=1.0)},
        models={},
        sweep_variables={"upstream_mach": _number(minimum=1.0)},
    ),
    "expansion": _CalculatorSchema(
        inputs={
            "upstream_mach": _number(minimum=1.0),
            "turn_angle": _NON_NEGATIVE,
        },
        models={},
        sweep_variables={
            "turn_angle": _NON_NEGATIVE,
            "mach": _number(minimum=1.0),
        },
    ),
    "detached_shock": _CalculatorSchema(
        inputs={"upstream_mach": _SUPERSONIC, "nose_radius": _POSITIVE},
        models={
            "geometry": _choice("axisymmetric_sphere", "cylindrical_nose_2d"),
            "model": _choice("ambrosio_wortman", "seiff", "comparison"),
        },
        sweep_variables={"upstream_mach": _SUPERSONIC},
    ),
    "boundary_layer_profile": _CalculatorSchema(
        inputs={
            "edge_velocity": _POSITIVE,
            "edge_density": _POSITIVE,
            "edge_temperature": _POSITIVE,
            "boundary_layer_thickness": _POSITIVE,
            "wall_shear_stress": _POSITIVE,
            "wall_temperature": _NULLABLE_POSITIVE,
            "wake_parameter": _number(nullable=True, minimum=0.0),
            "points": _integer(minimum=51, maximum=501),
        },
        models={
            "source": _choice("manual", "boundary"),
            "transformation": _choice("van_driest", "volpiani", "compare"),
            "temperature_velocity_relation": _choice(
                "generalized_reynolds_analogy", "walz"
            ),
        },
        sweep_variables=None,
    ),
    "protrusion_drag": _CalculatorSchema(
        inputs={
            "drag_coefficient": _NON_NEGATIVE,
            "height": _POSITIVE,
            "base_width": _POSITIVE,
            "edge_velocity": _POSITIVE,
            "edge_density": _POSITIVE,
            "boundary_layer_thickness": _POSITIVE,
            "mach": _number(nullable=True, minimum=0.0),
            "edge_temperature": _NULLABLE_POSITIVE,
            "wall_temperature": _NULLABLE_POSITIVE,
            "profile_height": _nullable_number_list(),
            "profile_velocity": _nullable_number_list(),
            "profile_density": _nullable_number_list(),
            "shape_height": _nullable_number_list(),
            "shape_width": _nullable_number_list(),
        },
        models={
            "profile_source": _choice("power_law", "saved", "csv"),
            "shape": _choice("rectangle", "triangle", "ellipse", "csv"),
            "compressible": _BOOL,
        },
        sweep_variables={
            "height": _POSITIVE,
            "drag_coefficient": _NON_NEGATIVE,
            "base_width": _POSITIVE,
            "boundary_layer_thickness": _POSITIVE,
            "mach": _NON_NEGATIVE,
        },
    ),
    "thermochemistry": _CalculatorSchema(
        inputs={
            "temperature": _POSITIVE,
            "pressure": _POSITIVE,
            "reference_temperature": _POSITIVE,
        },
        models={
            "selection": _choice("NASA7", "NASA9", "compare"),
            "allow_extrapolation": _BOOL,
        },
        sweep_variables={"temperature": _POSITIVE},
    ),
    "viscosity": _CalculatorSchema(
        inputs={"temperature": _POSITIVE},
        models={
            "selection": _choice("Sutherland", "Keyes", "Blottner/Wilke", "compare"),
            "allow_extrapolation": _BOOL,
        },
        sweep_variables={"temperature": _POSITIVE},
        sweep_extras={"scale": _choice("linear", "log")},
    ),
}


def calculator_names() -> frozenset[str]:
    """Return the calculator discriminators supported by schema version one."""
    return frozenset(_SCHEMAS)


def _validate_number(value: object, rule: _FieldRule, path: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise _SchemaError(f"{path} must be a number")
    try:
        normalized = float(value)
    except (OverflowError, TypeError, ValueError) as error:
        raise _SchemaError(f"{path} must be a finite number") from error
    if not math.isfinite(normalized):
        raise _SchemaError(f"{path} must be finite")
    if rule.minimum is not None:
        below = (
            normalized <= rule.minimum
            if rule.exclusive_minimum
            else normalized < rule.minimum
        )
        if below:
            comparison = "greater than" if rule.exclusive_minimum else "at least"
            raise _SchemaError(f"{path} must be {comparison} {rule.minimum:g}")
    if rule.maximum is not None and normalized > rule.maximum:
        raise _SchemaError(f"{path} must be at most {rule.maximum:g}")
    return normalized


def _validate_field(value: object, rule: _FieldRule, path: str) -> object:
    if value is None:
        if rule.nullable:
            return None
        raise _SchemaError(f"{path} must not be null")
    if rule.kind == "number":
        return _validate_number(value, rule, path)
    if rule.kind == "int":
        if isinstance(value, bool) or not isinstance(value, int):
            raise _SchemaError(f"{path} must be an integer")
        normalized = int(value)
        if rule.minimum is not None and normalized < rule.minimum:
            raise _SchemaError(f"{path} must be at least {int(rule.minimum)}")
        if rule.maximum is not None and normalized > rule.maximum:
            raise _SchemaError(f"{path} must be at most {int(rule.maximum)}")
        return normalized
    if rule.kind == "bool":
        if not isinstance(value, bool):
            raise _SchemaError(f"{path} must be a boolean")
        return value
    if rule.kind == "string":
        if not isinstance(value, str):
            raise _SchemaError(f"{path} must be a string")
        if rule.choices is not None and value not in rule.choices:
            choices = ", ".join(sorted(rule.choices))
            raise _SchemaError(f"{path} must be one of: {choices}")
        return value
    if not isinstance(value, list):
        raise _SchemaError(f"{path} must be an array of numbers")
    return [
        _validate_number(item, _NUMBER, f"{path}[{index}]")
        for index, item in enumerate(value)
    ]


def _validate_object(
    value: object, rules: dict[str, _FieldRule], path: str
) -> dict[str, object]:
    if not isinstance(value, dict):
        raise _SchemaError(f"{path} must be an object")
    unknown = set(value) - set(rules)
    if unknown:
        names = ", ".join(sorted(str(name) for name in unknown))
        raise _SchemaError(f"{path} contains unsupported fields: {names}")
    missing = {name for name, rule in rules.items() if rule.required} - set(value)
    if missing:
        names = ", ".join(sorted(missing))
        raise _SchemaError(f"{path} is missing required fields: {names}")
    return {
        name: _validate_field(item, rules[name], f"{path}.{name}")
        for name, item in value.items()
    }


def _validate_sweep(value: object, schema: _CalculatorSchema) -> dict[str, object]:
    assert schema.sweep_variables is not None
    if not isinstance(value, dict):
        raise _SchemaError("sweep_si must be an object")
    field = value.get("field")
    if not isinstance(field, str) or field not in schema.sweep_variables:
        choices = ", ".join(sorted(schema.sweep_variables))
        raise _SchemaError(f"sweep_si.field must be one of: {choices}")
    rules = {
        "field": _choice(*schema.sweep_variables),
        "start": schema.sweep_variables[field],
        "stop": schema.sweep_variables[field],
        "points": _integer(minimum=2, maximum=schema.maximum_sweep_points),
        **(schema.sweep_extras or {}),
    }
    normalized = _validate_object(value, rules, "sweep_si")
    start = normalized["start"]
    stop = normalized["stop"]
    assert isinstance(start, float)
    assert isinstance(stop, float)
    if start >= stop:
        raise _SchemaError("sweep_si.start must be less than sweep_si.stop")
    if normalized.get("logarithmic") is True and start <= 0.0:
        raise _SchemaError("a logarithmic sweep requires a positive start")
    if normalized.get("scale") == "log" and start <= 0.0:
        raise _SchemaError("a logarithmic sweep requires a positive start")
    return normalized


def _require_matching_lists(
    inputs: dict[str, object], names: tuple[str, ...], *, path: str
) -> None:
    values = [inputs[name] for name in names]
    present = [value is not None for value in values]
    if any(present) and not all(present):
        raise _SchemaError(f"{path} fields must be supplied together")
    if not all(present):
        return
    arrays = [value for value in values if isinstance(value, list)]
    if len(arrays) != len(values):
        raise _SchemaError(f"{path} fields must be numeric arrays")
    if len(arrays[0]) < 2 or len({len(value) for value in arrays}) != 1:
        raise _SchemaError(f"{path} arrays must have equal lengths of at least two")


def _validate_cross_fields(
    calculator: str,
    inputs: dict[str, object],
    models: dict[str, object],
    sweep: dict[str, object] | None,
) -> None:
    if calculator == "isentropic":
        basis = models["input_basis"]
        input_value = inputs["input_value"]
        assert isinstance(input_value, float)
        if basis != "mach" and input_value < 1.0:
            raise _SchemaError(
                "inputs_si.input_value must be at least one for ratio inputs"
            )
        if sweep is not None and basis == "mach":
            start = sweep["start"]
            assert isinstance(start, float)
            if start < 0.0:
                raise _SchemaError("sweep_si.start must be non-negative for Mach")
        if sweep is not None and basis != "mach":
            start = sweep["start"]
            assert isinstance(start, float)
            if start < 1.0:
                raise _SchemaError(
                    "sweep_si.start must be at least one for ratio inputs"
                )
        gas_model = models.get("gas_model", "AIR")
        temperature = inputs["total_temperature"]
        pressure = inputs["total_pressure"]
        with_mass_flux = models["with_mass_flux"]
        if gas_model != "AIR" and temperature is None:
            raise _SchemaError(
                "inputs_si.total_temperature is required for this gas model"
            )
        if with_mass_flux and (temperature is None or pressure is None):
            raise _SchemaError(
                "total_temperature and total_pressure are required for mass flux"
            )
        if gas_model == "BEATTIE_BRIDGEMAN" and not with_mass_flux:
            raise _SchemaError(
                "Beattie-Bridgeman configurations require mass-flux conditions"
            )
    elif calculator == "boundary_layer":
        regime = models["regime"]
        transition = inputs["transition_reynolds"]
        if regime == "transitional" and transition is None:
            raise _SchemaError(
                "inputs_si.transition_reynolds is required for transitional flow"
            )
        correction = models["compressibility_correction"]
        if correction != "none" and inputs["mach"] is None:
            raise _SchemaError(
                "inputs_si.mach is required for a compressibility correction"
            )
    elif calculator == "detached_shock":
        if (
            models["geometry"] == "cylindrical_nose_2d"
            and models["model"] != "ambrosio_wortman"
        ):
            raise _SchemaError(
                "cylindrical geometry supports only the ambrosio_wortman model"
            )
    elif calculator == "protrusion_drag":
        profile_names = ("profile_height", "profile_velocity", "profile_density")
        shape_names = ("shape_height", "shape_width")
        _require_matching_lists(inputs, profile_names, path="embedded profile")
        _require_matching_lists(inputs, shape_names, path="embedded shape")
        profile_present = inputs["profile_height"] is not None
        shape_present = inputs["shape_height"] is not None
        if models["profile_source"] in {"saved", "csv"} and not profile_present:
            raise _SchemaError("the selected profile source requires embedded arrays")
        if models["profile_source"] == "power_law" and profile_present:
            raise _SchemaError("power_law profiles must not include embedded arrays")
        if models["shape"] == "csv" and not shape_present:
            raise _SchemaError("the csv shape requires embedded shape arrays")
        if models["shape"] != "csv" and shape_present:
            raise _SchemaError(
                "a built-in shape must not include embedded shape arrays"
            )
        if models["compressible"]:
            if models["profile_source"] != "power_law":
                raise _SchemaError("compressibility is available only for power_law")
            if inputs["mach"] is None or inputs["edge_temperature"] is None:
                raise _SchemaError(
                    "compressible protrusion profiles require mach and edge_temperature"
                )
        if sweep is not None:
            field = sweep["field"]
            if field == "mach" and not models["compressible"]:
                raise _SchemaError("a Mach sweep requires compressibility")
            if field == "base_width" and models["shape"] == "csv":
                raise _SchemaError("base_width cannot be swept for a csv shape")


def validate_calculator_payload(
    *,
    calculator: str,
    mode: str,
    inputs_si: object,
    models: object,
    sweep_si: object | None,
    has_sweep: bool,
) -> tuple[dict[str, object], dict[str, object], dict[str, object] | None]:
    """Validate and normalize one calculator-discriminated payload."""
    schema = _SCHEMAS[calculator]
    normalized_inputs = _validate_object(inputs_si, schema.inputs, "inputs_si")
    normalized_models = _validate_object(models, schema.models, "models")
    if mode == "single":
        if has_sweep:
            raise _SchemaError("single mode does not allow sweep_si")
        normalized_sweep = None
    else:
        if schema.sweep_variables is None:
            raise _SchemaError(f"{calculator} does not support sweep mode")
        if not has_sweep:
            raise _SchemaError("sweep mode requires sweep_si")
        normalized_sweep = _validate_sweep(sweep_si, schema)
    _validate_cross_fields(
        calculator, normalized_inputs, normalized_models, normalized_sweep
    )
    return normalized_inputs, normalized_models, normalized_sweep


__all__ = ["_SchemaError", "calculator_names", "validate_calculator_payload"]
