"""Pure adapters between GUI requests and the public calculation API."""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

from aerophysics import FlightCondition
from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    CompressibilityCorrection,
    TurbulentCorrelation,
    flat_plate_boundary_layer,
)
from aerophysics.shocks import (
    ShockBranch,
    maximum_attached_deflection,
    oblique_shock,
)

type CellValue = float | str | None
type Row = dict[str, CellValue]


@dataclass(frozen=True, slots=True)
class CalculationResult:
    """Tabular SI results and captured model warnings."""

    rows: tuple[Row, ...]
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class FlightCase:
    """Scalar edge conditions passed from the flight page to boundary layers."""

    geometric_altitude: float
    mach: float
    velocity: float
    density: float
    dynamic_viscosity: float
    temperature: float

    @classmethod
    def from_row(cls, row: Row) -> FlightCase:
        """Create a case from a successful scalar flight-result row."""

        def number(name: str) -> float:
            value = row[name]
            if not isinstance(value, float):
                raise ValueError(f"{name} is not numeric")
            return value

        return cls(
            geometric_altitude=number("geometric_altitude"),
            mach=number("mach"),
            velocity=number("velocity"),
            density=number("density"),
            dynamic_viscosity=number("dynamic_viscosity"),
            temperature=number("temperature"),
        )


def sweep_values(
    start: float, stop: float, points: int, *, log: bool = False
) -> np.ndarray:
    """Return a validated bounded sweep grid."""
    if not np.isfinite(start) or not np.isfinite(stop):
        raise ValueError("sweep bounds must be finite")
    if start >= stop:
        raise ValueError("sweep start must be less than sweep stop")
    if not 2 <= points <= 501:
        raise ValueError("sweep points must be between 2 and 501")
    if log:
        if start <= 0.0:
            raise ValueError("logarithmic sweep start must be positive")
        return np.geomspace(start, stop, points, dtype=np.float64)
    return np.linspace(start, stop, points, dtype=np.float64)


def _array(value: object) -> np.ndarray:
    return np.atleast_1d(np.asarray(value, dtype=np.float64))


def _optional_at(value: object, index: int) -> float | None:
    if value is None:
        return None
    return float(_array(value)[index])


def flight_condition(
    *,
    geometric_altitude: float | np.ndarray,
    motion: float | np.ndarray,
    motion_basis: str,
    characteristic_length: float | None,
) -> CalculationResult:
    """Calculate scalar or broadcast flight conditions."""
    if motion_basis == "mach":
        condition = FlightCondition.from_mach(
            geometric_altitude, motion, characteristic_length
        )
    elif motion_basis == "velocity":
        condition = FlightCondition.from_velocity(
            geometric_altitude, motion, characteristic_length
        )
    else:
        raise ValueError("motion_basis must be mach or velocity")
    atmosphere = condition.atmosphere
    altitude = _array(atmosphere.geometric_altitude)
    rows: list[Row] = []
    for index in range(altitude.size):
        rows.append(
            {
                "geometric_altitude": float(altitude[index]),
                "geopotential_altitude": float(
                    _array(atmosphere.geopotential_altitude)[index]
                ),
                "temperature": float(_array(atmosphere.temperature)[index]),
                "pressure": float(_array(atmosphere.pressure)[index]),
                "density": float(_array(atmosphere.density)[index]),
                "speed_of_sound": float(_array(atmosphere.speed_of_sound)[index]),
                "dynamic_viscosity": float(_array(atmosphere.dynamic_viscosity)[index]),
                "mach": float(_array(condition.mach)[index]),
                "velocity": float(_array(condition.velocity)[index]),
                "dynamic_pressure": float(_array(condition.dynamic_pressure)[index]),
                "reynolds_number_per_length": float(
                    _array(condition.reynolds_number_per_length)[index]
                ),
                "reynolds_number": _optional_at(condition.reynolds_number, index),
                "total_temperature": float(_array(condition.total_temperature)[index]),
                "total_pressure": float(_array(condition.total_pressure)[index]),
                "total_density": float(_array(condition.total_density)[index]),
                "status": "ok",
                "message": "",
            }
        )
    return CalculationResult(tuple(rows))


def flight_sweep(
    *,
    fixed_altitude: float,
    fixed_motion: float,
    motion_basis: str,
    sweep_field: str,
    start: float,
    stop: float,
    points: int,
    characteristic_length: float | None,
) -> CalculationResult:
    """Sweep altitude or the selected motion quantity."""
    values = sweep_values(start, stop, points)
    if sweep_field == "altitude":
        altitude: float | np.ndarray = values
        motion: float | np.ndarray = fixed_motion
    elif sweep_field == "motion":
        altitude = fixed_altitude
        motion = values
    else:
        raise ValueError("sweep_field must be altitude or motion")
    return flight_condition(
        geometric_altitude=altitude,
        motion=motion,
        motion_basis=motion_basis,
        characteristic_length=characteristic_length,
    )


def oblique_shock_condition(
    *,
    upstream_mach: float,
    deflection_angle: float,
    branch: ShockBranch,
) -> CalculationResult:
    """Calculate one oblique-shock state."""
    result = oblique_shock(upstream_mach, deflection_angle, branch)
    limit = maximum_attached_deflection(upstream_mach)
    row: Row = {
        "upstream_mach": float(result.upstream_mach),
        "downstream_mach": float(result.downstream_mach),
        "deflection_angle": float(result.deflection_angle),
        "shock_angle": float(result.shock_angle),
        "maximum_deflection_angle": float(limit.deflection_angle),
        "upstream_normal_mach": float(result.upstream_normal_mach),
        "downstream_normal_mach": float(result.downstream_normal_mach),
        "static_pressure_ratio": float(result.static_pressure_ratio),
        "static_density_ratio": float(result.static_density_ratio),
        "static_temperature_ratio": float(result.static_temperature_ratio),
        "total_pressure_ratio": float(result.total_pressure_ratio),
        "status": "ok",
        "message": "",
    }
    return CalculationResult((row,))


def oblique_shock_sweep(
    *,
    fixed_mach: float,
    fixed_deflection: float,
    branch: ShockBranch,
    sweep_field: str,
    start: float,
    stop: float,
    points: int,
) -> CalculationResult:
    """Sweep Mach or deflection while retaining non-attached rows."""
    values = sweep_values(start, stop, points)
    rows: list[Row] = []
    for value in values:
        mach = float(value) if sweep_field == "mach" else fixed_mach
        theta = float(value) if sweep_field == "deflection" else fixed_deflection
        try:
            result = oblique_shock_condition(
                upstream_mach=mach,
                deflection_angle=theta,
                branch=branch,
            )
        except ValueError as error:
            maximum: float | None = None
            if mach > 1.0:
                maximum = float(maximum_attached_deflection(mach).deflection_angle)
            rows.append(
                {
                    "upstream_mach": mach,
                    "downstream_mach": None,
                    "deflection_angle": theta,
                    "shock_angle": None,
                    "maximum_deflection_angle": maximum,
                    "upstream_normal_mach": None,
                    "downstream_normal_mach": None,
                    "static_pressure_ratio": None,
                    "static_density_ratio": None,
                    "static_temperature_ratio": None,
                    "total_pressure_ratio": None,
                    "status": "no_attached_shock",
                    "message": str(error),
                }
            )
        else:
            rows.append(result.rows[0])
    if sweep_field not in {"mach", "deflection"}:
        raise ValueError("sweep_field must be mach or deflection")
    return CalculationResult(tuple(rows))


def flat_plate(
    *,
    distance: float | np.ndarray,
    edge_velocity: float,
    edge_density: float,
    edge_dynamic_viscosity: float,
    regime: BoundaryLayerRegime,
    turbulent_correlation: TurbulentCorrelation,
    transition_reynolds: float | None,
    compressibility_correction: CompressibilityCorrection,
    mach: float | None,
    edge_temperature: float | None,
    wall_temperature: float | None,
) -> CalculationResult:
    """Calculate a scalar or distance-array flat-plate boundary layer."""
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        result = flat_plate_boundary_layer(
            distance,
            edge_velocity,
            edge_density,
            edge_dynamic_viscosity,
            regime=regime,
            turbulent_correlation=turbulent_correlation,
            transition_reynolds=transition_reynolds,
            compressibility_correction=compressibility_correction,
            mach=mach,
            edge_temperature=edge_temperature,
            wall_temperature=wall_temperature,
        )
    distances = _array(result.distance)
    rows: list[Row] = []
    for index in range(distances.size):
        rows.append(
            {
                "distance": float(distances[index]),
                "reynolds_number": float(_array(result.reynolds_number)[index]),
                "effective_reynolds_number": float(
                    _array(result.effective_reynolds_number)[index]
                ),
                "boundary_layer_thickness": float(
                    _array(result.boundary_layer_thickness)[index]
                ),
                "displacement_thickness": float(
                    _array(result.displacement_thickness)[index]
                ),
                "momentum_thickness": float(_array(result.momentum_thickness)[index]),
                "local_skin_friction_coefficient": float(
                    _array(result.local_skin_friction_coefficient)[index]
                ),
                "average_skin_friction_coefficient": float(
                    _array(result.average_skin_friction_coefficient)[index]
                ),
                "wall_shear_stress": float(_array(result.wall_shear_stress)[index]),
                "drag_per_unit_width": float(_array(result.drag_per_unit_width)[index]),
                "recovery_temperature": _optional_at(
                    result.recovery_temperature, index
                ),
                "wall_temperature": _optional_at(result.wall_temperature, index),
                "status": "ok",
                "message": "",
            }
        )
    return CalculationResult(tuple(rows), tuple(str(item.message) for item in captured))


def flat_plate_sweep(
    *,
    start: float,
    stop: float,
    points: int,
    logarithmic: bool,
    edge_velocity: float,
    edge_density: float,
    edge_dynamic_viscosity: float,
    regime: BoundaryLayerRegime,
    turbulent_correlation: TurbulentCorrelation,
    transition_reynolds: float | None,
    compressibility_correction: CompressibilityCorrection,
    mach: float | None,
    edge_temperature: float | None,
    wall_temperature: float | None,
) -> CalculationResult:
    """Calculate a bounded distance sweep."""
    return flat_plate(
        distance=sweep_values(start, stop, points, log=logarithmic),
        edge_velocity=edge_velocity,
        edge_density=edge_density,
        edge_dynamic_viscosity=edge_dynamic_viscosity,
        regime=regime,
        turbulent_correlation=turbulent_correlation,
        transition_reynolds=transition_reynolds,
        compressibility_correction=compressibility_correction,
        mach=mach,
        edge_temperature=edge_temperature,
        wall_temperature=wall_temperature,
    )
