"""Derived flow-state columns, independent of the shock solvers."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from aerophysics._temperature import (
    restore_static_temperature as restore_static_temperature,
)
from aerophysics.gas import PerfectGas
from aerophysics.real_gas import BeattieBridgemanGas
from aerophysics.shocks import ShockGasModel
from aerophysics.transport import (
    AIR_BLOTTNER_VISCOSITY,
    AIR_KEYES_VISCOSITY,
    AIR_VISCOSITY,
    DynamicViscosityModel,
)

if TYPE_CHECKING:
    from aerophysics.gui.adapters import Row

VISCOSITY_MODELS: dict[str, DynamicViscosityModel] = {
    "Sutherland": AIR_VISCOSITY,
    "Keyes": AIR_KEYES_VISCOSITY,
    "Blottner/Wilke": AIR_BLOTTNER_VISCOSITY,
}
VISCOSITY_RANGES: dict[str, tuple[float, float] | None] = {
    "Sutherland": None,
    "Keyes": (79.0, 1845.0),
    "Blottner/Wilke": (1000.0, 30_000.0),
}


def validate_output_inputs(
    pressure: float | None, viscosity_model: str, characteristic_length: float | None
) -> None:
    """Reject invalid optional inputs before computing any flow state."""
    if viscosity_model not in VISCOSITY_MODELS:
        raise ValueError("viscosity_model must be Sutherland, Keyes, or Blottner/Wilke")
    for name, value in (
        ("pressure", pressure),
        ("characteristic_length", characteristic_length),
    ):
        if value is not None and (not np.isfinite(value) or value <= 0.0):
            raise ValueError(f"{name} must be finite and greater than zero")
    if characteristic_length is not None and pressure is None:
        raise ValueError("pressure is required when characteristic_length is specified")


def heat_capacities(
    gas: ShockGasModel | BeattieBridgemanGas,
    temperature: float | None,
    pressure: float | None = None,
    *,
    allow_extrapolation: bool = False,
) -> tuple[float, float, float]:
    """Return local gamma, cp and cv with the selected thermodynamic model."""
    if isinstance(gas, PerfectGas):
        return gas.heat_capacity_ratio, gas.cp, gas.cv
    assert temperature is not None
    if isinstance(gas, BeattieBridgemanGas):
        assert pressure is not None
        cp = float(
            gas.cp(temperature, pressure, allow_extrapolation=allow_extrapolation)
        )
        cv = float(
            gas.cv(temperature, pressure, allow_extrapolation=allow_extrapolation)
        )
    else:
        cp = float(gas.cp(temperature, allow_extrapolation=allow_extrapolation))
        cv = float(gas.cv(temperature, allow_extrapolation=allow_extrapolation))
    return cp / cv, cp, cv


def add_transport(
    row: Row,
    *,
    prefix: str = "",
    temperature: float,
    density: float,
    velocity: float,
    viscosity_model: str,
    characteristic_length: float | None,
) -> str:
    """Add local mu/Re; omit only transport outputs outside a fitted range."""
    row["viscosity_model"] = viscosity_model
    row[f"{prefix}dynamic_viscosity"] = None
    row[f"{prefix}reynolds_number_per_length"] = None
    if characteristic_length is not None:
        row[f"{prefix}reynolds_number"] = None
    nominal_range = VISCOSITY_RANGES[viscosity_model]
    if nominal_range is not None and not (
        nominal_range[0] <= temperature <= nominal_range[1]
    ):
        return (
            f"{prefix or 'static_'}T={temperature:g} K: {viscosity_model} "
            f"({nominal_range[0]:g}–{nominal_range[1]:g} K) の範囲外のため μ・Re は空欄"
        )
    viscosity = float(VISCOSITY_MODELS[viscosity_model].dynamic_viscosity(temperature))
    unit_reynolds = density * velocity / viscosity
    row[f"{prefix}dynamic_viscosity"] = viscosity
    row[f"{prefix}reynolds_number_per_length"] = unit_reynolds
    if characteristic_length is not None:
        row[f"{prefix}reynolds_number"] = unit_reynolds * characteristic_length
    return ""


def add_shock_outputs(
    row: Row,
    gas: ShockGasModel,
    upstream_temperature: float | None,
    upstream_pressure: float | None,
    viscosity_model: str,
    characteristic_length: float | None,
    with_heat_capacities: bool,
    *,
    surface: bool = False,
) -> tuple[str, ...]:
    """Add full velocity magnitudes and upstream/downstream or surface states."""
    if upstream_temperature is not None and (
        not np.isfinite(upstream_temperature) or upstream_temperature <= 0.0
    ):
        raise ValueError("upstream_temperature must be finite and greater than zero")
    end = "surface" if surface else "downstream"
    ratio_prefix = "surface" if surface else "static"
    temperature_ratio = float(row[f"{ratio_prefix}_temperature_ratio"])  # type: ignore[arg-type]
    end_temperature = (
        restore_static_temperature(upstream_temperature * temperature_ratio, gas)
        if upstream_temperature is not None
        else None
    )
    gamma1, cp1, cv1 = heat_capacities(gas, upstream_temperature)
    gamma2, cp2, cv2 = heat_capacities(gas, end_temperature)
    row["upstream_heat_capacity_ratio"] = gamma1
    row[f"{end}_heat_capacity_ratio"] = gamma2
    if with_heat_capacities:
        row.update({"upstream_cp": cp1, "upstream_cv": cv1})
        row.update({f"{end}_cp": cp2, f"{end}_cv": cv2})
    mach1 = float(row["upstream_mach"])  # type: ignore[arg-type]
    mach2 = float(row[f"{end}_mach"])  # type: ignore[arg-type]
    velocity_ratio = mach2 / mach1 * np.sqrt(gamma2 / gamma1 * temperature_ratio)
    row["velocity_ratio"] = float(velocity_ratio)
    row["dynamic_pressure_ratio"] = float(
        float(row[f"{ratio_prefix}_density_ratio"]) * velocity_ratio**2  # type: ignore[arg-type]
    )
    row["entropy_change_over_r"] = float(
        -np.log(float(row["total_pressure_ratio"]))  # type: ignore[arg-type]
    )
    if upstream_temperature is None:
        if upstream_pressure is not None:
            raise ValueError("upstream_temperature is required with upstream_pressure")
        return ()
    assert end_temperature is not None
    messages: list[str] = []
    for prefix, temperature, gamma, mach, pressure_ratio in (
        ("upstream", upstream_temperature, gamma1, mach1, 1.0),
        (
            end,
            end_temperature,
            gamma2,
            mach2,
            float(row[f"{ratio_prefix}_pressure_ratio"]),  # type: ignore[arg-type]
        ),
    ):
        sound_speed = float(np.sqrt(gamma * gas.specific_gas_constant * temperature))
        velocity = mach * sound_speed
        row[f"{prefix}_temperature"] = temperature
        row[f"{prefix}_speed_of_sound"] = sound_speed
        row[f"{prefix}_velocity"] = velocity
        if upstream_pressure is not None:
            pressure = upstream_pressure * pressure_ratio
            density = pressure / (gas.specific_gas_constant * temperature)
            row[f"{prefix}_pressure"] = pressure
            row[f"{prefix}_density"] = density
            row[f"{prefix}_dynamic_pressure"] = 0.5 * density * velocity**2
            message = add_transport(
                row,
                prefix=f"{prefix}_",
                temperature=temperature,
                density=density,
                velocity=velocity,
                viscosity_model=viscosity_model,
                characteristic_length=characteristic_length,
            )
            if message:
                messages.append(message)
    if messages:
        row["message"] = "; ".join(filter(None, (str(row["message"]), *messages)))
    return tuple(messages)
