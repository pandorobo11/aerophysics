"""Pure GUI adapters for viscous-profile, protrusion, and property models."""

from __future__ import annotations

import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np

from aerophysics import AIR_NASA7, AIR_NASA9
from aerophysics.boundary_layer_profile import (
    CompressibleVelocityTransformation,
    TemperatureVelocityRelation,
    compressible_turbulent_boundary_layer_profile,
)
from aerophysics.gui.adapters import CalculationResult, Row, sweep_values
from aerophysics.protrusion import protrusion_drag
from aerophysics.transport import (
    AIR_BLOTTNER_VISCOSITY,
    AIR_KEYES_VISCOSITY,
    AIR_VISCOSITY,
    DynamicViscosityModel,
)


@dataclass(frozen=True, slots=True)
class BoundaryLayerCase:
    """Scalar turbulent flat-plate state passed to the profile calculator."""

    edge_velocity: float
    edge_density: float
    edge_temperature: float
    boundary_layer_thickness: float
    wall_shear_stress: float


@dataclass(frozen=True, slots=True)
class BoundaryProfileCase:
    """SI profile passed from profile prediction to protrusion drag."""

    transformation: CompressibleVelocityTransformation
    wall_distance: np.ndarray
    velocity: np.ndarray
    density: np.ndarray
    edge_velocity: float
    edge_density: float
    edge_temperature: float
    boundary_layer_thickness: float


@dataclass(frozen=True, slots=True)
class ProfileCalculation:
    """Tabular result plus one transferable profile per transformation."""

    result: CalculationResult
    profiles: tuple[BoundaryProfileCase, ...]


def wall_normal_grid(thickness: float, points: int) -> np.ndarray:
    """Create a wall-resolving grid from zero through ``delta_99``."""
    if not np.isfinite(thickness) or thickness <= 0.0:
        raise ValueError("boundary_layer_thickness must be finite and positive")
    if not 51 <= points <= 501:
        raise ValueError("profile points must be between 51 and 501")
    positive = np.geomspace(thickness * 1.0e-7, thickness, points - 1)
    return np.concatenate(([0.0], positive)).astype(np.float64)


def boundary_layer_profiles(
    *,
    edge_velocity: float,
    edge_density: float,
    edge_temperature: float,
    boundary_layer_thickness: float,
    wall_shear_stress: float,
    transformations: Sequence[CompressibleVelocityTransformation],
    temperature_velocity_relation: TemperatureVelocityRelation,
    wall_temperature: float | None,
    wake_parameter: float | None,
    points: int,
) -> ProfileCalculation:
    """Predict one or two compressible turbulent boundary-layer profiles."""
    selected = tuple(transformations)
    if not selected or len(set(selected)) != len(selected):
        raise ValueError("transformations must contain unique model choices")
    grid = wall_normal_grid(boundary_layer_thickness, points)
    rows: list[Row] = []
    profiles: list[BoundaryProfileCase] = []
    messages: list[str] = []
    for transformation in selected:
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            output = compressible_turbulent_boundary_layer_profile(
                grid,
                edge_velocity,
                edge_density,
                edge_temperature,
                boundary_layer_thickness,
                wall_shear_stress,
                transformation=transformation,
                wall_temperature=wall_temperature,
                temperature_velocity_relation=temperature_velocity_relation,
                wake_parameter=wake_parameter,
            )
        messages.extend(str(item.message) for item in captured)
        profiles.append(
            BoundaryProfileCase(
                transformation=transformation,
                wall_distance=output.wall_distance.copy(),
                velocity=output.velocity.copy(),
                density=output.density.copy(),
                edge_velocity=edge_velocity,
                edge_density=edge_density,
                edge_temperature=edge_temperature,
                boundary_layer_thickness=boundary_layer_thickness,
            )
        )
        for index, height in enumerate(output.wall_distance):
            rows.append(
                {
                    "model": transformation.value,
                    "wall_distance": float(height),
                    "outer_coordinate": float(height / boundary_layer_thickness),
                    "wall_distance_plus": float(output.wall_distance_plus[index]),
                    "transformed_wall_coordinate": float(
                        output.transformed_wall_coordinate[index]
                    ),
                    "velocity": float(output.velocity[index]),
                    "velocity_ratio": float(output.velocity[index] / edge_velocity),
                    "velocity_plus": float(output.velocity_plus[index]),
                    "transformed_velocity_plus": float(
                        output.transformed_velocity_plus[index]
                    ),
                    "temperature": float(output.temperature[index]),
                    "density": float(output.density[index]),
                    "dynamic_viscosity": float(output.dynamic_viscosity[index]),
                    "local_mach_number": float(output.local_mach_number[index]),
                    "dynamic_pressure": float(output.dynamic_pressure[index]),
                    "friction_velocity": output.friction_velocity,
                    "friction_reynolds_number": output.friction_reynolds_number,
                    "recovery_temperature": output.recovery_temperature,
                    "wall_temperature": output.wall_temperature,
                    "wake_parameter": output.wake_parameter,
                    "displacement_thickness": output.displacement_thickness,
                    "momentum_thickness": output.momentum_thickness,
                    "shape_factor": output.shape_factor,
                    "local_skin_friction_coefficient": (
                        output.local_skin_friction_coefficient
                    ),
                    "status": "ok",
                    "message": "",
                }
            )
    return ProfileCalculation(
        CalculationResult(tuple(rows), tuple(dict.fromkeys(messages))),
        tuple(profiles),
    )


def _representative_width(
    shape: str, height: float, base_width: float
) -> float | Callable[[np.ndarray], np.ndarray]:
    if shape == "rectangle":
        return base_width
    if shape == "triangle":
        return lambda y: base_width * np.maximum(1.0 - y / height, 0.0)
    if shape == "ellipse":
        return lambda y: base_width * np.sqrt(np.maximum(1.0 - (y / height) ** 2, 0.0))
    raise ValueError("shape must be rectangle, triangle, ellipse, or csv")


def _width_input(
    *,
    shape: str,
    height: float,
    base_width: float,
    shape_height: np.ndarray | None,
    shape_width: np.ndarray | None,
) -> float | Callable[[np.ndarray], np.ndarray]:
    if shape != "csv":
        return _representative_width(shape, height, base_width)
    if shape_height is None or shape_width is None:
        raise ValueError("CSV shape requires height and width arrays")
    if shape_height[-1] < height:
        raise ValueError("CSV shape must cover the protrusion height")
    return lambda y: np.interp(y, shape_height, shape_width)


def protrusion_condition(
    *,
    drag_coefficient: float,
    height: float,
    base_width: float,
    shape: str,
    edge_velocity: float,
    edge_density: float,
    boundary_layer_thickness: float,
    profile_height: np.ndarray | None = None,
    profile_velocity: np.ndarray | None = None,
    profile_density: np.ndarray | None = None,
    mach: float | None = None,
    edge_temperature: float | None = None,
    wall_temperature: float | None = None,
    shape_height: np.ndarray | None = None,
    shape_width: np.ndarray | None = None,
) -> CalculationResult:
    """Calculate direct drag for a representative or uploaded frontal shape."""
    width = _width_input(
        shape=shape,
        height=height,
        base_width=base_width,
        shape_height=shape_height,
        shape_width=shape_width,
    )
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        output = protrusion_drag(
            drag_coefficient,
            height,
            width,
            edge_velocity,
            edge_density,
            boundary_layer_thickness,
            profile_height=profile_height,
            profile_velocity=profile_velocity,
            profile_density=profile_density,
            mach=mach,
            edge_temperature=edge_temperature,
            wall_temperature=wall_temperature,
        )
    row: Row = {
        "drag_coefficient": drag_coefficient,
        "height": height,
        "base_width": base_width,
        "boundary_layer_thickness": boundary_layer_thickness,
        "shape": shape,
        "mach": mach,
        "direct_drag": output.direct_drag,
        "effective_dynamic_pressure": output.effective_dynamic_pressure,
        "shielding_factor": output.shielding_factor,
        "frontal_area": output.frontal_area,
        "edge_dynamic_pressure": output.edge_dynamic_pressure,
        "height_to_boundary_layer_thickness": (
            output.height_to_boundary_layer_thickness
        ),
        "profile": output.profile.value,
        "compressibility_applied": str(output.compressibility_applied),
        "status": "ok",
        "message": "",
    }
    return CalculationResult((row,), tuple(str(item.message) for item in captured))


def protrusion_sweep(
    *,
    sweep_field: str,
    start: float,
    stop: float,
    points: int,
    drag_coefficient: float,
    height: float,
    base_width: float,
    shape: str,
    edge_velocity: float,
    edge_density: float,
    boundary_layer_thickness: float,
    profile_height: np.ndarray | None = None,
    profile_velocity: np.ndarray | None = None,
    profile_density: np.ndarray | None = None,
    mach: float | None = None,
    edge_temperature: float | None = None,
    wall_temperature: float | None = None,
    shape_height: np.ndarray | None = None,
    shape_width: np.ndarray | None = None,
) -> CalculationResult:
    """Sweep one protrusion input and retain invalid points as gaps."""
    allowed = {
        "height",
        "drag_coefficient",
        "base_width",
        "boundary_layer_thickness",
        "mach",
    }
    if sweep_field not in allowed:
        raise ValueError("unsupported protrusion sweep field")
    if sweep_field == "mach" and mach is None:
        raise ValueError("Mach sweep requires compressible mode")
    rows: list[Row] = []
    messages: list[str] = []
    for value in sweep_values(start, stop, points):
        current_drag = (
            float(value) if sweep_field == "drag_coefficient" else drag_coefficient
        )
        current_height = float(value) if sweep_field == "height" else height
        current_width = float(value) if sweep_field == "base_width" else base_width
        current_thickness = (
            float(value)
            if sweep_field == "boundary_layer_thickness"
            else boundary_layer_thickness
        )
        current_mach = float(value) if sweep_field == "mach" else mach
        try:
            result = protrusion_condition(
                drag_coefficient=current_drag,
                height=current_height,
                base_width=current_width,
                shape=shape,
                edge_velocity=edge_velocity,
                edge_density=edge_density,
                boundary_layer_thickness=current_thickness,
                profile_height=profile_height,
                profile_velocity=profile_velocity,
                profile_density=profile_density,
                mach=current_mach,
                edge_temperature=edge_temperature,
                wall_temperature=wall_temperature,
                shape_height=shape_height,
                shape_width=shape_width,
            )
        except ValueError as error:
            rows.append(
                {
                    "drag_coefficient": current_drag,
                    "height": current_height,
                    "base_width": current_width,
                    "boundary_layer_thickness": current_thickness,
                    "shape": shape,
                    "mach": current_mach,
                    "direct_drag": None,
                    "effective_dynamic_pressure": None,
                    "shielding_factor": None,
                    "frontal_area": None,
                    "edge_dynamic_pressure": None,
                    "height_to_boundary_layer_thickness": None,
                    "profile": "",
                    "compressibility_applied": "",
                    "status": "invalid",
                    "message": str(error),
                }
            )
        else:
            rows.extend(result.rows)
            messages.extend(result.warnings)
    return CalculationResult(tuple(rows), tuple(dict.fromkeys(messages)))


_THERMO_MODELS = {"NASA7": AIR_NASA7, "NASA9": AIR_NASA9}


def thermochemistry_condition(
    *,
    temperature: float | np.ndarray,
    pressure: float,
    reference_temperature: float,
    models: Sequence[str],
    allow_extrapolation: bool,
) -> CalculationResult:
    """Calculate frozen dry-air thermodynamic properties."""
    selected = tuple(models)
    if not selected or len(set(selected)) != len(selected):
        raise ValueError("models must contain unique NASA model choices")
    temperatures = np.atleast_1d(np.asarray(temperature, dtype=np.float64))
    rows: list[Row] = []
    messages: list[str] = []
    for name in selected:
        try:
            gas = _THERMO_MODELS[name]
        except KeyError as error:
            raise ValueError("model must be NASA7 or NASA9") from error
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            cp = np.atleast_1d(
                gas.cp(temperatures, allow_extrapolation=allow_extrapolation)
            )
            cv = np.atleast_1d(
                gas.cv(temperatures, allow_extrapolation=allow_extrapolation)
            )
            gamma = np.atleast_1d(
                gas.heat_capacity_ratio(
                    temperatures, allow_extrapolation=allow_extrapolation
                )
            )
            sound = np.atleast_1d(
                gas.speed_of_sound(
                    temperatures, allow_extrapolation=allow_extrapolation
                )
            )
            standard_h = np.atleast_1d(
                gas.standard_enthalpy(
                    temperatures, allow_extrapolation=allow_extrapolation
                )
            )
            standard_u = np.atleast_1d(
                gas.standard_internal_energy(
                    temperatures, allow_extrapolation=allow_extrapolation
                )
            )
            sensible_h = np.atleast_1d(
                gas.sensible_enthalpy(
                    temperatures,
                    reference_temperature,
                    allow_extrapolation=allow_extrapolation,
                )
            )
            sensible_u = np.atleast_1d(
                gas.sensible_internal_energy(
                    temperatures,
                    reference_temperature,
                    allow_extrapolation=allow_extrapolation,
                )
            )
            entropy = np.atleast_1d(
                gas.entropy(
                    temperatures,
                    pressure,
                    allow_extrapolation=allow_extrapolation,
                )
            )
        messages.extend(str(item.message) for item in captured)
        for index, value in enumerate(temperatures):
            rows.append(
                {
                    "model": name,
                    "temperature": float(value),
                    "pressure": pressure,
                    "reference_temperature": reference_temperature,
                    "molar_mass": gas.molar_mass,
                    "specific_gas_constant": gas.specific_gas_constant,
                    "cp": float(cp[index]),
                    "cv": float(cv[index]),
                    "heat_capacity_ratio": float(gamma[index]),
                    "speed_of_sound": float(sound[index]),
                    "standard_enthalpy": float(standard_h[index]),
                    "standard_internal_energy": float(standard_u[index]),
                    "sensible_enthalpy": float(sensible_h[index]),
                    "sensible_internal_energy": float(sensible_u[index]),
                    "entropy": float(entropy[index]),
                    "status": "ok",
                    "message": "",
                }
            )
    return CalculationResult(tuple(rows), tuple(dict.fromkeys(messages)))


def thermochemistry_sweep(
    *,
    start: float,
    stop: float,
    points: int,
    pressure: float,
    reference_temperature: float,
    models: Sequence[str],
    allow_extrapolation: bool,
) -> CalculationResult:
    """Calculate a bounded temperature sweep for frozen dry air."""
    return thermochemistry_condition(
        temperature=sweep_values(start, stop, points),
        pressure=pressure,
        reference_temperature=reference_temperature,
        models=models,
        allow_extrapolation=allow_extrapolation,
    )


_VISCOSITY_MODELS: dict[str, DynamicViscosityModel] = {
    "Sutherland": AIR_VISCOSITY,
    "Keyes": AIR_KEYES_VISCOSITY,
    "Blottner/Wilke": AIR_BLOTTNER_VISCOSITY,
}
_VISCOSITY_RANGES: dict[str, tuple[float, float] | None] = {
    "Sutherland": None,
    "Keyes": (79.0, 1845.0),
    "Blottner/Wilke": (1000.0, 30_000.0),
}


def viscosity_condition(
    *,
    temperature: float | np.ndarray,
    models: Sequence[str],
    allow_extrapolation: bool,
) -> CalculationResult:
    """Calculate dry-air dynamic viscosity with one or more models."""
    selected = tuple(models)
    if not selected or len(set(selected)) != len(selected):
        raise ValueError("models must contain unique viscosity model choices")
    unknown = tuple(name for name in selected if name not in _VISCOSITY_MODELS)
    if unknown:
        raise ValueError("model must be Sutherland, Keyes, or Blottner/Wilke")

    temperatures = np.atleast_1d(np.asarray(temperature, dtype=np.float64)).reshape(-1)
    if not np.all(np.isfinite(temperatures)) or np.any(temperatures <= 0.0):
        raise ValueError("temperature must be finite and greater than zero")
    baseline = np.atleast_1d(
        np.asarray(AIR_VISCOSITY.dynamic_viscosity(temperatures), dtype=np.float64)
    )
    rows: list[Row] = []
    messages: list[str] = []

    for name in selected:
        model = _VISCOSITY_MODELS[name]
        nominal_range = _VISCOSITY_RANGES[name]
        valid = np.ones(temperatures.shape, dtype=np.bool_)
        if nominal_range is not None:
            minimum, maximum = nominal_range
            valid = (temperatures >= minimum) & (temperatures <= maximum)

        values = np.full(temperatures.shape, np.nan, dtype=np.float64)
        if allow_extrapolation:
            with warnings.catch_warnings(record=True) as captured:
                warnings.simplefilter("always")
                values = np.atleast_1d(
                    np.asarray(model.dynamic_viscosity(temperatures), dtype=np.float64)
                )
            messages.extend(str(item.message) for item in captured)
        elif np.any(valid):
            values[valid] = np.atleast_1d(
                np.asarray(
                    model.dynamic_viscosity(temperatures[valid]), dtype=np.float64
                )
            )

        omitted = int(np.count_nonzero(~valid))
        if omitted and not allow_extrapolation and nominal_range is not None:
            minimum, maximum = nominal_range
            messages.append(
                f"{name}: {omitted} temperature value(s) outside the nominal "
                f"range {minimum:g}--{maximum:g} K were omitted"
            )

        for index, value in enumerate(temperatures):
            if valid[index]:
                viscosity = float(values[index])
                relative = float((viscosity / baseline[index] - 1.0) * 100.0)
                status = "ok"
                message = ""
            elif allow_extrapolation and nominal_range is not None:
                minimum, maximum = nominal_range
                viscosity = float(values[index])
                relative = float((viscosity / baseline[index] - 1.0) * 100.0)
                status = "extrapolated"
                message = (
                    f"temperature is outside the nominal range "
                    f"{minimum:g}--{maximum:g} K"
                )
            else:
                minimum, maximum = nominal_range or (0.0, 0.0)
                viscosity = None
                relative = None
                status = "out_of_range"
                message = (
                    f"temperature is outside the nominal range "
                    f"{minimum:g}--{maximum:g} K"
                )
            rows.append(
                {
                    "model": name,
                    "temperature": float(value),
                    "dynamic_viscosity": viscosity,
                    "relative_difference": relative,
                    "status": status,
                    "message": message,
                }
            )

    return CalculationResult(tuple(rows), tuple(dict.fromkeys(messages)))


def viscosity_sweep(
    *,
    start: float,
    stop: float,
    points: int,
    models: Sequence[str],
    allow_extrapolation: bool,
    log_temperature: bool,
) -> CalculationResult:
    """Calculate a linear or logarithmic dry-air viscosity sweep."""
    return viscosity_condition(
        temperature=sweep_values(start, stop, points, log=log_temperature),
        models=models,
        allow_extrapolation=allow_extrapolation,
    )
