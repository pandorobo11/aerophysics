"""Pure adapters between GUI requests and the public calculation API."""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

from aerophysics import (
    AIR_BEATTIE_BRIDGEMAN,
    AIR_HARMONIC_OSCILLATOR,
    AIR_NASA7,
    AIR_NASA9,
    FlightCondition,
)
from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    CompressibilityCorrection,
    TurbulentCorrelation,
    flat_plate_boundary_layer,
)
from aerophysics.detached_shock import (
    BilligShockShapeResult,
    DetachedShockGeometry,
    billig_shock_shape,
    compare_standoff_distances,
    seiff_standoff_distance_from_mach,
    shock_standoff_distance,
)
from aerophysics.expansion import (
    maximum_prandtl_meyer_angle,
    prandtl_meyer_angle,
    prandtl_meyer_expansion,
)
from aerophysics.gas import AIR, PerfectGas
from aerophysics.isentropic import (
    MachBranch,
    area_ratio,
    choked_mass_flux,
    critical_ratios,
    isentropic_ratios,
    isentropic_state,
    mach_from_area_ratio,
    mach_from_total_density_ratio,
    mach_from_total_pressure_ratio,
    mach_from_total_temperature_ratio,
    mass_flow_parameter,
    mass_flux,
)
from aerophysics.real_gas import BeattieBridgemanGas, HarmonicOscillatorGas
from aerophysics.shocks import (
    ShockBranch,
    conical_shock,
    maximum_attached_cone_angle,
    maximum_attached_deflection,
    normal_shock,
    oblique_shock,
    supersonic_pitot_pressure_ratio,
)
from aerophysics.thermochemistry import ThermallyPerfectGas

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


type _GuiIsentropicGas = (
    PerfectGas | ThermallyPerfectGas | HarmonicOscillatorGas | BeattieBridgemanGas
)


_ISENTROPIC_GASES: dict[str, _GuiIsentropicGas] = {
    "AIR": AIR,
    "NASA7": AIR_NASA7,
    "NASA9": AIR_NASA9,
    "HARMONIC_OSCILLATOR": AIR_HARMONIC_OSCILLATOR,
    "BEATTIE_BRIDGEMAN": AIR_BEATTIE_BRIDGEMAN,
}


def _mach_from_isentropic_input(
    value: float,
    basis: str,
    branch: MachBranch,
    gas: _GuiIsentropicGas,
    total_temperature: float | None,
    total_pressure: float | None,
    *,
    allow_extrapolation: bool,
) -> float:
    if basis == "mach":
        if value < 0.0:
            raise ValueError("mach must be non-negative")
        return value
    inverse_functions = {
        "temperature_ratio": mach_from_total_temperature_ratio,
        "pressure_ratio": mach_from_total_pressure_ratio,
        "density_ratio": mach_from_total_density_ratio,
    }
    if basis == "area_ratio":
        return float(
            mach_from_area_ratio(
                value,
                branch,
                gas,
                total_temperature=total_temperature,
                total_pressure=total_pressure,
                allow_extrapolation=allow_extrapolation,
            )
        )
    try:
        inverse = inverse_functions[basis]
    except KeyError as error:
        raise ValueError("unsupported isentropic input basis") from error
    return float(
        inverse(
            value,
            gas,
            total_temperature=total_temperature,
            total_pressure=total_pressure,
            allow_extrapolation=allow_extrapolation,
        )
    )


def isentropic_condition(
    *,
    input_value: float | np.ndarray,
    input_basis: str,
    branch: MachBranch = MachBranch.SUBSONIC,
    gas_model: str = "AIR",
    total_pressure: float | None = None,
    total_temperature: float | None = None,
    allow_extrapolation: bool = True,
) -> CalculationResult:
    """Calculate isentropic state, inverse, area, and mass-flow relations."""
    try:
        gas = _ISENTROPIC_GASES[gas_model]
    except KeyError as error:
        raise ValueError(
            "gas_model must be AIR, NASA7, NASA9, HARMONIC_OSCILLATOR, "
            "or BEATTIE_BRIDGEMAN"
        ) from error
    if isinstance(gas, (ThermallyPerfectGas, HarmonicOscillatorGas)) and (
        total_temperature is None
    ):
        raise ValueError("total_temperature is required for a thermally perfect gas")
    if isinstance(gas, BeattieBridgemanGas) and (
        total_temperature is None or total_pressure is None
    ):
        raise ValueError(
            "total_temperature and total_pressure are required for a "
            "Beattie--Bridgeman gas"
        )
    if total_pressure is not None and total_temperature is None:
        raise ValueError(
            "total_temperature is required when total_pressure is specified"
        )

    rows: list[Row] = []
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        critical = critical_ratios(
            gas,
            total_temperature=total_temperature,
            total_pressure=total_pressure,
            allow_extrapolation=allow_extrapolation,
        )
        for raw_value in _array(input_value):
            value = float(raw_value)
            mach = _mach_from_isentropic_input(
                value,
                input_basis,
                branch,
                gas,
                total_temperature,
                total_pressure,
                allow_extrapolation=allow_extrapolation,
            )
            ratios = isentropic_ratios(
                mach,
                gas,
                total_temperature=total_temperature,
                total_pressure=total_pressure,
                allow_extrapolation=allow_extrapolation,
            )
            flux = (
                float(
                    mass_flux(
                        total_pressure,
                        total_temperature,
                        mach,
                        gas,
                        allow_extrapolation=allow_extrapolation,
                    )
                )
                if total_pressure is not None and total_temperature is not None
                else None
            )
            choked = (
                float(
                    choked_mass_flux(
                        total_pressure,
                        total_temperature,
                        gas,
                        allow_extrapolation=allow_extrapolation,
                    )
                )
                if total_pressure is not None and total_temperature is not None
                else None
            )
            absolute_state = (
                isentropic_state(
                    mach,
                    gas,
                    total_temperature=total_temperature,
                    total_pressure=total_pressure,
                    allow_extrapolation=allow_extrapolation,
                )
                if total_pressure is not None and total_temperature is not None
                else None
            )
            rows.append(
                {
                    "gas_model": gas_model,
                    "input_value": value,
                    "input_basis": input_basis,
                    "mach": mach,
                    "total_temperature": total_temperature,
                    "static_temperature": (
                        total_temperature / float(ratios.total_temperature_ratio)
                        if total_temperature is not None
                        else None
                    ),
                    "static_pressure": (
                        float(absolute_state.static_pressure)
                        if absolute_state is not None
                        else None
                    ),
                    "static_density": (
                        float(absolute_state.static_density)
                        if absolute_state is not None
                        else None
                    ),
                    "velocity": (
                        float(absolute_state.velocity)
                        if absolute_state is not None
                        else None
                    ),
                    "speed_of_sound": (
                        float(absolute_state.speed_of_sound)
                        if absolute_state is not None
                        else None
                    ),
                    "dynamic_pressure": (
                        float(absolute_state.dynamic_pressure)
                        if absolute_state is not None
                        else None
                    ),
                    "total_temperature_ratio": float(ratios.total_temperature_ratio),
                    "total_pressure_ratio": float(ratios.total_pressure_ratio),
                    "total_density_ratio": float(ratios.total_density_ratio),
                    "area_ratio": (
                        float(
                            area_ratio(
                                mach,
                                gas,
                                total_temperature=total_temperature,
                                total_pressure=total_pressure,
                                allow_extrapolation=allow_extrapolation,
                            )
                        )
                        if mach > 0.0
                        else None
                    ),
                    "mass_flow_parameter": float(
                        mass_flow_parameter(
                            mach,
                            gas,
                            total_temperature=total_temperature,
                            total_pressure=total_pressure,
                            allow_extrapolation=allow_extrapolation,
                        )
                    ),
                    "mass_flux": flux,
                    "choked_mass_flux": choked,
                    "critical_temperature_ratio": float(
                        critical.total_temperature_ratio
                    ),
                    "critical_pressure_ratio": float(critical.total_pressure_ratio),
                    "critical_density_ratio": float(critical.total_density_ratio),
                    "status": "ok",
                    "message": "",
                }
            )
    messages = tuple(dict.fromkeys(str(item.message) for item in captured))
    return CalculationResult(tuple(rows), messages)


def isentropic_sweep(
    *,
    input_basis: str,
    branch: MachBranch,
    start: float,
    stop: float,
    points: int,
    gas_model: str = "AIR",
    total_pressure: float | None = None,
    total_temperature: float | None = None,
    allow_extrapolation: bool = True,
) -> CalculationResult:
    """Sweep the selected isentropic input quantity."""
    return isentropic_condition(
        input_value=sweep_values(start, stop, points),
        input_basis=input_basis,
        branch=branch,
        gas_model=gas_model,
        total_pressure=total_pressure,
        total_temperature=total_temperature,
        allow_extrapolation=allow_extrapolation,
    )


def normal_shock_condition(*, upstream_mach: float | np.ndarray) -> CalculationResult:
    """Calculate normal-shock ratios and the Rayleigh pitot relation."""
    result = normal_shock(upstream_mach)
    mach_values = _array(result.upstream_mach)
    pitot_values = _array(supersonic_pitot_pressure_ratio(mach_values))
    rows: list[Row] = []
    for index, mach in enumerate(mach_values):
        rows.append(
            {
                "upstream_mach": float(mach),
                "downstream_mach": float(_array(result.downstream_mach)[index]),
                "static_pressure_ratio": float(
                    _array(result.static_pressure_ratio)[index]
                ),
                "static_density_ratio": float(
                    _array(result.static_density_ratio)[index]
                ),
                "static_temperature_ratio": float(
                    _array(result.static_temperature_ratio)[index]
                ),
                "total_pressure_ratio": float(
                    _array(result.total_pressure_ratio)[index]
                ),
                "pitot_pressure_ratio": float(pitot_values[index]),
                "status": "ok",
                "message": "",
            }
        )
    return CalculationResult(tuple(rows))


def normal_shock_sweep(*, start: float, stop: float, points: int) -> CalculationResult:
    """Sweep upstream Mach number through a normal shock."""
    return normal_shock_condition(upstream_mach=sweep_values(start, stop, points))


def detached_shock_condition(
    *,
    upstream_mach: float | np.ndarray,
    nose_radius: float,
    geometry: DetachedShockGeometry,
    selection: str,
) -> CalculationResult:
    """Calculate detached-shock standoff and Billig curvature data."""
    if selection not in {"ambrosio_wortman", "seiff", "comparison"}:
        raise ValueError(
            "selection must be ambrosio_wortman, seiff, or comparison"
        )
    if (
        geometry is DetachedShockGeometry.CYLINDRICAL_NOSE_2D
        and selection != "ambrosio_wortman"
    ):
        raise ValueError("Seiff and comparison are available only for sphere geometry")

    aw = shock_standoff_distance(
        upstream_mach,
        nose_radius,
        geometry=geometry,
    )
    curvature = billig_shock_shape(
        upstream_mach,
        nose_radius,
        [0.0],
        geometry=geometry,
    )
    mach_values = _array(aw.upstream_mach)
    aw_normalized = _array(aw.normalized_standoff_distance)
    aw_distance = _array(aw.standoff_distance)
    curvature_radius = _array(curvature.vertex_curvature_radius)

    seiff_normalized: np.ndarray | None = None
    seiff_distance: np.ndarray | None = None
    density_ratio: np.ndarray | None = None
    normalized_difference: np.ndarray | None = None
    distance_difference: np.ndarray | None = None
    relative_difference: np.ndarray | None = None
    if geometry is DetachedShockGeometry.AXISYMMETRIC_SPHERE:
        seiff = seiff_standoff_distance_from_mach(upstream_mach, nose_radius)
        seiff_normalized = _array(seiff.normalized_standoff_distance)
        seiff_distance = _array(seiff.standoff_distance)
        density_ratio = _array(seiff.density_ratio)
        comparison = compare_standoff_distances(upstream_mach, nose_radius)
        normalized_difference = _array(
            comparison.normalized_standoff_difference
        )
        distance_difference = _array(comparison.standoff_distance_difference)
        relative_difference = _array(comparison.relative_difference)

    rows: list[Row] = []
    for index, mach in enumerate(mach_values):
        selected_normalized = (
            float(seiff_normalized[index])
            if selection == "seiff" and seiff_normalized is not None
            else float(aw_normalized[index])
            if selection == "ambrosio_wortman"
            else None
        )
        selected_distance = (
            float(seiff_distance[index])
            if selection == "seiff" and seiff_distance is not None
            else float(aw_distance[index])
            if selection == "ambrosio_wortman"
            else None
        )
        rows.append(
            {
                "upstream_mach": float(mach),
                "nose_radius": nose_radius,
                "geometry": geometry.value,
                "selection": selection,
                "normalized_standoff_distance": selected_normalized,
                "standoff_distance": selected_distance,
                "aw_normalized_standoff_distance": float(aw_normalized[index]),
                "aw_standoff_distance": float(aw_distance[index]),
                "seiff_density_ratio": _optional_at(density_ratio, index),
                "seiff_normalized_standoff_distance": _optional_at(
                    seiff_normalized, index
                ),
                "seiff_standoff_distance": _optional_at(seiff_distance, index),
                "normalized_standoff_difference": _optional_at(
                    normalized_difference, index
                ),
                "standoff_distance_difference": _optional_at(
                    distance_difference, index
                ),
                "relative_difference": _optional_at(relative_difference, index),
                "billig_vertex_curvature_radius": float(curvature_radius[index]),
                "status": "ok",
                "message": "",
            }
        )
    return CalculationResult(tuple(rows))


def detached_shock_sweep(
    *,
    start: float,
    stop: float,
    points: int,
    nose_radius: float,
    geometry: DetachedShockGeometry,
    selection: str,
) -> CalculationResult:
    """Sweep Mach number for detached-shock correlations."""
    return detached_shock_condition(
        upstream_mach=sweep_values(start, stop, points),
        nose_radius=nose_radius,
        geometry=geometry,
        selection=selection,
    )


def detached_shock_shape(
    *,
    upstream_mach: float,
    nose_radius: float,
    geometry: DetachedShockGeometry,
    points: int = 401,
    span_radii: float = 2.0,
) -> BilligShockShapeResult:
    """Return the GUI's symmetric Billig shock-shape sampling."""
    if not 3 <= points <= 2001 or points % 2 == 0:
        raise ValueError("points must be an odd integer between 3 and 2001")
    if not np.isfinite(span_radii) or span_radii <= 0.0:
        raise ValueError("span_radii must be greater than zero")
    transverse = np.linspace(
        -span_radii * nose_radius,
        span_radii * nose_radius,
        points,
        dtype=np.float64,
    )
    return billig_shock_shape(
        upstream_mach,
        nose_radius,
        transverse,
        geometry=geometry,
    )


def expansion_condition(
    *, upstream_mach: float, turn_angle: float
) -> CalculationResult:
    """Calculate one centered Prandtl-Meyer expansion."""
    result = prandtl_meyer_expansion(upstream_mach, turn_angle)
    maximum_turn = maximum_prandtl_meyer_angle() - float(
        prandtl_meyer_angle(upstream_mach)
    )
    return CalculationResult(
        (
            {
                "upstream_mach": float(result.upstream_mach),
                "downstream_mach": float(result.downstream_mach),
                "turn_angle": float(result.turn_angle),
                "maximum_turn_angle": maximum_turn,
                "upstream_prandtl_meyer_angle": float(
                    result.upstream_prandtl_meyer_angle
                ),
                "downstream_prandtl_meyer_angle": float(
                    result.downstream_prandtl_meyer_angle
                ),
                "static_temperature_ratio": float(result.static_temperature_ratio),
                "static_pressure_ratio": float(result.static_pressure_ratio),
                "static_density_ratio": float(result.static_density_ratio),
                "status": "ok",
                "message": "",
            },
        )
    )


def expansion_sweep(
    *,
    fixed_mach: float,
    fixed_turn_angle: float,
    sweep_field: str,
    start: float,
    stop: float,
    points: int,
) -> CalculationResult:
    """Sweep expansion Mach or turn angle while retaining limit failures."""
    if sweep_field not in {"mach", "turn_angle"}:
        raise ValueError("sweep_field must be mach or turn_angle")
    rows: list[Row] = []
    for value in sweep_values(start, stop, points):
        mach = float(value) if sweep_field == "mach" else fixed_mach
        turn = float(value) if sweep_field == "turn_angle" else fixed_turn_angle
        try:
            rows.append(
                expansion_condition(upstream_mach=mach, turn_angle=turn).rows[0]
            )
        except ValueError as error:
            maximum_turn = (
                maximum_prandtl_meyer_angle() - float(prandtl_meyer_angle(mach))
                if mach >= 1.0
                else None
            )
            rows.append(
                {
                    "upstream_mach": mach,
                    "downstream_mach": None,
                    "turn_angle": turn,
                    "maximum_turn_angle": maximum_turn,
                    "upstream_prandtl_meyer_angle": None,
                    "downstream_prandtl_meyer_angle": None,
                    "static_temperature_ratio": None,
                    "static_pressure_ratio": None,
                    "static_density_ratio": None,
                    "status": "limit_exceeded",
                    "message": str(error),
                }
            )
    return CalculationResult(tuple(rows))


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


def conical_shock_condition(
    *, upstream_mach: float, cone_half_angle: float
) -> CalculationResult:
    """Calculate one axisymmetric Taylor-Maccoll conical-shock state."""
    result = conical_shock(upstream_mach, cone_half_angle)
    limit = maximum_attached_cone_angle(upstream_mach)
    row: Row = {
        "upstream_mach": float(result.upstream_mach),
        "cone_half_angle": float(result.cone_half_angle),
        "maximum_cone_half_angle": float(limit.cone_half_angle),
        "shock_angle": float(result.shock_angle),
        "post_shock_mach": float(result.post_shock_mach),
        "surface_mach": float(result.surface_mach),
        "surface_pressure_ratio": float(result.surface_pressure_ratio),
        "surface_density_ratio": float(result.surface_density_ratio),
        "surface_temperature_ratio": float(result.surface_temperature_ratio),
        "total_pressure_ratio": float(result.total_pressure_ratio),
        "status": "ok",
        "message": "",
    }
    return CalculationResult((row,))


def conical_shock_sweep(
    *,
    fixed_mach: float,
    fixed_cone_half_angle: float,
    sweep_field: str,
    start: float,
    stop: float,
    points: int,
) -> CalculationResult:
    """Sweep Mach or cone half-angle while retaining non-attached rows."""
    if sweep_field not in {"mach", "cone_half_angle"}:
        raise ValueError("sweep_field must be mach or cone_half_angle")
    values = sweep_values(start, stop, points)
    rows: list[Row] = []
    for value in values:
        mach = float(value) if sweep_field == "mach" else fixed_mach
        angle = (
            float(value) if sweep_field == "cone_half_angle" else fixed_cone_half_angle
        )
        try:
            result = conical_shock_condition(upstream_mach=mach, cone_half_angle=angle)
        except ValueError as error:
            maximum: float | None = None
            if mach > 1.0:
                maximum = float(maximum_attached_cone_angle(mach).cone_half_angle)
            rows.append(
                {
                    "upstream_mach": mach,
                    "cone_half_angle": angle,
                    "maximum_cone_half_angle": maximum,
                    "shock_angle": None,
                    "post_shock_mach": None,
                    "surface_mach": None,
                    "surface_pressure_ratio": None,
                    "surface_density_ratio": None,
                    "surface_temperature_ratio": None,
                    "total_pressure_ratio": None,
                    "status": "no_attached_shock",
                    "message": str(error),
                }
            )
        else:
            rows.append(result.rows[0])
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
