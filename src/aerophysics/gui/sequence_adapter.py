"""Validated sequence configuration and tabular SI output for the GUI."""

from dataclasses import asdict
from math import isfinite
from typing import Any

from aerophysics.flow_sequence import (
    ConicalShockStep,
    ExpansionStep,
    FlowProperties,
    FlowSequenceResult,
    FlowState,
    FlowStep,
    IsentropicStep,
    NormalShockStep,
    ObliqueShockStep,
    StepResult,
    solve_flow_sequence,
)
from aerophysics.gui.adapters import _SHOCK_GASES
from aerophysics.isentropic import MachBranch
from aerophysics.shocks import ShockBranch

KINDS = ("normal", "oblique", "conical", "expansion", "isentropic")
BASES = (
    "static_mach",
    "static_velocity",
    "total",
    "atmosphere_mach",
    "atmosphere_velocity",
)
INITIAL_FIELDS = {
    "static_mach": ("pressure", "temperature", "mach"),
    "static_velocity": ("pressure", "temperature", "velocity"),
    "total": ("pressure", "temperature", "mach"),
    "atmosphere_mach": ("altitude", "mach"),
    "atmosphere_velocity": ("altitude", "velocity"),
}


def _number(value: object) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (float, int))
        or not isfinite(value)
    ):
        raise ValueError("sequence inputs must be finite numbers")
    return float(value)


def decode_step(data: dict[str, Any]) -> FlowStep:
    """Validate one stage, rejecting ignored or unsupported parameters."""
    kind = data.get("kind")
    name = data.get("name")
    if kind not in KINDS or not isinstance(name, str):
        raise ValueError("invalid stage kind or name")
    expected = {"kind", "name"}
    if kind in {"oblique", "conical", "expansion"}:
        expected.add("angle")
    if kind == "oblique":
        expected.add("branch")
    if kind == "isentropic":
        expected |= {"basis", "value", "branch"}
    if set(data) != expected:
        raise ValueError("stage contains missing or unsupported fields")
    if kind == "normal":
        return NormalShockStep(name=name)
    if kind == "isentropic":
        if data["basis"] not in {"mach", "pressure_ratio", "area_ratio"}:
            raise ValueError("unknown isentropic basis")
        branch = None if data["branch"] is None else MachBranch(data["branch"])
        if data["basis"] == "area_ratio" and branch is None:
            raise ValueError("area ratio requires an exit branch")
        return IsentropicStep(_number(data["value"]), data["basis"], branch, name)
    angle = _number(data["angle"])
    if kind == "oblique":
        return ObliqueShockStep(angle, ShockBranch(data["branch"]), name)
    if kind == "conical":
        return ConicalShockStep(angle, name)
    return ExpansionStep(angle, name)


def validate_sequence_payload(
    inputs: object, models: object
) -> tuple[dict[str, object], dict[str, object]]:
    """Validate JSON structure without invoking any physical solver."""
    if (
        not isinstance(models, dict)
        or set(models) != {"gas_model"}
        or models["gas_model"] not in _SHOCK_GASES
    ):
        raise ValueError("unknown sequence gas model")
    if not isinstance(inputs, dict) or set(inputs) != {"initial", "steps"}:
        raise ValueError("sequence requires initial and steps")
    initial, steps = inputs["initial"], inputs["steps"]
    if not isinstance(initial, dict) or initial.get("basis") not in BASES:
        raise ValueError("unknown initial state basis")
    fields = INITIAL_FIELDS[initial["basis"]]
    if set(initial) != {"basis", *fields}:
        raise ValueError("initial state contains missing or unsupported fields")
    normalized = {
        "basis": initial["basis"],
        **{key: _number(initial[key]) for key in fields},
    }
    if not isinstance(steps, list) or not all(isinstance(step, dict) for step in steps):
        raise ValueError("steps must be a list of stage objects")
    for step in steps:
        decode_step(step)
    return {"initial": normalized, "steps": steps}, dict(models)


def calculate_sequence(
    inputs: dict[str, Any], models: dict[str, Any]
) -> FlowSequenceResult:
    """Resolve the initial state and call the common numerical API."""
    validate_sequence_payload(inputs, models)
    data = inputs["initial"]
    gas = _SHOCK_GASES[models["gas_model"]]
    basis = data["basis"]
    if basis.startswith("atmosphere"):
        initial = FlowState.from_atmosphere(data["altitude"], data.get("mach", 0.0))
        if basis == "atmosphere_velocity":
            initial = FlowState.from_velocity(
                initial.pressure, initial.temperature, data["velocity"], gas
            )
    elif basis == "total":
        initial = FlowState.from_total(
            data["pressure"], data["temperature"], data["mach"], gas
        )
    elif basis == "static_velocity":
        initial = FlowState.from_velocity(
            data["pressure"], data["temperature"], data["velocity"], gas
        )
    else:
        initial = FlowState(data["pressure"], data["temperature"], data["mach"])
    return solve_flow_sequence(
        initial, [decode_step(step) for step in inputs["steps"]], gas
    )


def sequence_rows(result: FlowSequenceResult) -> list[dict[str, Any]]:
    """Return initial and stage rows, with failed stages retained."""
    rows: list[dict[str, Any]] = []
    records: list[
        tuple[
            str,
            str,
            FlowProperties | None,
            float | None,
            float | None,
            StepResult | None,
        ]
    ] = [("Initial", "ok", result.initial, 1.0, 0.0, None)]
    records.extend(
        (
            stage.step.name,
            stage.status,
            stage.properties,
            stage.total_pressure_recovery,
            stage.entropy_change,
            stage,
        )
        for stage in result.steps
    )
    for index, (name, status, properties, recovery, entropy, stage) in enumerate(
        records
    ):
        row: dict[str, Any] = {"stage": index, "name": name, "status": status}
        if properties is not None:
            values = asdict(properties)
            row.update(values.pop("state"))
            row.update(values)
        row.update(total_pressure_recovery=recovery, entropy_change=entropy)
        if stage is not None:
            row.update(
                pressure_ratio=stage.pressure_ratio,
                temperature_ratio=stage.temperature_ratio,
                shock_angle=stage.shock_angle,
                post_shock_mach=stage.post_shock_mach,
                message=stage.message,
            )
        rows.append(row)
    return rows
