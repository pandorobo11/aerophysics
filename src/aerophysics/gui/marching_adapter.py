"""Pure SI boundaries for the marching boundary-layer calculator."""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

import numpy as np

from aerophysics.boundary_layer import BoundaryLayerRegime
from aerophysics.characteristics import AxisymmetricBodyFlowResult
from aerophysics.exceptions import BoundaryLayerSeparationError
from aerophysics.gas import AIR
from aerophysics.gui.adapters import CalculationResult, Row
from aerophysics.marching_boundary_layer import (
    BoundaryLayerEdge,
    MarchingBoundaryLayerResult,
    marching_boundary_layer,
)


@dataclass(frozen=True, slots=True)
class EdgeDefinition:
    """Edge pressure along the wall, optional wall radius and total conditions."""

    distance: np.ndarray
    pressure: np.ndarray
    radius: np.ndarray | None
    total_pressure: float
    total_temperature: float


def total_conditions(
    mach: float, pressure: float, temperature: float
) -> tuple[float, float]:
    """Return AIR total pressure [Pa] and total temperature [K]."""
    gamma = AIR.heat_capacity_ratio
    ratio = 1.0 + 0.5 * (gamma - 1.0) * mach**2
    return pressure * ratio ** (gamma / (gamma - 1.0)), temperature * ratio


def uniform_edge(
    *,
    mach: float = 2.0,
    pressure: float = 10000.0,
    temperature: float = 220.0,
    length: float = 1.0,
    station_count: int = 50,
    geometry: str = "planar",
    cone_half_angle: float | None = None,
) -> EdgeDefinition:
    """Build a constant-pressure edge on a flat plate or a sharp cone.

    ``mach``, ``pressure`` and ``temperature`` are the local edge state, which
    on a cone is the cone-surface state behind the conical shock. ``length``
    is the wall distance in metres, divided into ``station_count`` intervals.
    """
    if not all(
        math.isfinite(value) and value > 0
        for value in (mach, pressure, temperature, length)
    ):
        raise ValueError(
            "mach, pressure, temperature and length must be finite and positive"
        )
    if (
        isinstance(station_count, bool)
        or not isinstance(station_count, int)
        or not 2 <= station_count <= 1000
    ):
        raise ValueError("station_count must be an integer from 2 to 1000")
    distance = np.linspace(0.0, length, station_count + 1)
    if geometry == "planar":
        radius = None
    elif geometry == "cone":
        if (
            cone_half_angle is None
            or not math.isfinite(cone_half_angle)
            or not 0 < cone_half_angle < math.pi / 2
        ):
            raise ValueError("cone_half_angle must satisfy 0 < angle < pi/2")
        radius = distance * math.sin(cone_half_angle)
    else:
        raise ValueError("unsupported geometry")
    total_pressure, total_temperature = total_conditions(mach, pressure, temperature)
    return EdgeDefinition(
        distance,
        np.full(distance.size, pressure),
        radius,
        total_pressure,
        total_temperature,
    )


def body_case_edge(case: AxisymmetricBodyFlowResult) -> EdgeDefinition:
    """Take the body-streamline edge of a saved case with default corner ramps."""
    edge, radius = case.boundary_layer_inputs()
    return EdgeDefinition(
        np.array(edge.distance),
        np.array(edge.pressure),
        np.array(radius),
        case.total_pressure,
        case.total_temperature,
    )


@dataclass(frozen=True, slots=True)
class MarchingCalculation:
    """Station rows, the public solver result and the edge it was marched on.

    ``separation_distance`` is set when the march stopped at separation; the
    rows and ``raw`` then hold only the stations solved upstream of it.
    """

    stations: CalculationResult
    raw: MarchingBoundaryLayerResult
    edge: BoundaryLayerEdge
    separation_distance: float | None = None


def calculate_marching_boundary_layer(
    definition: EdgeDefinition,
    *,
    regime: str = "laminar",
    wall_temperature: float | None = None,
    transition_distance: float | None = None,
    prandtl_number: float = 0.72,
    turbulent_prandtl_number: float = 0.9,
    edge_sa_ratio: float = 3.0,
) -> MarchingCalculation:
    """March the AIR/Sutherland boundary layer once and tabulate each station.

    A separation that leaves at least one solved station is returned with
    ``separation_distance`` set; other solver errors propagate.
    """
    edge = BoundaryLayerEdge.from_pressure(
        definition.distance,
        definition.pressure,
        total_pressure=definition.total_pressure,
        total_temperature=definition.total_temperature,
    )
    separation = None
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        try:
            raw = marching_boundary_layer(
                edge,
                regime=BoundaryLayerRegime(regime),
                wall_temperature=wall_temperature,
                prandtl_number=prandtl_number,
                turbulent_prandtl_number=turbulent_prandtl_number,
                edge_sa_ratio=edge_sa_ratio,
                wall_radius=definition.radius,
                transition_distance=transition_distance,
            )
        except BoundaryLayerSeparationError as error:
            if error.result.distance.size == 0:
                raise
            raw, separation = error.result, error.separation_distance
    stations = slice(1, raw.distance.size + 1)
    temperature = edge.temperature[stations]
    velocity = edge.velocity[stations]
    values = {
        "surface_distance": raw.distance,
        "edge_pressure": edge.pressure[stations],
        "edge_temperature": temperature,
        "edge_mach": velocity
        / np.sqrt(AIR.heat_capacity_ratio * AIR.specific_gas_constant * temperature),
        "edge_velocity": velocity,
        "reynolds_number": raw.reynolds_number,
        "boundary_layer_thickness": raw.boundary_layer_thickness,
        "displacement_thickness": raw.displacement_thickness,
        "momentum_thickness": raw.momentum_thickness,
        "shape_factor": raw.displacement_thickness / raw.momentum_thickness,
        "local_skin_friction_coefficient": raw.local_skin_friction_coefficient,
        "wall_shear_stress": raw.wall_shear_stress,
        "wall_temperature": raw.wall_temperature,
        "wall_heat_flux": raw.wall_heat_flux,
        "pressure_gradient_parameter": raw.pressure_gradient_parameter,
    }
    if definition.radius is not None:
        values["radius"] = np.asarray(definition.radius, dtype=np.float64)[stations]
    rows: tuple[Row, ...] = tuple(
        {key: float(value[index]) for key, value in values.items()}
        for index in range(raw.distance.size)
    )
    return MarchingCalculation(
        CalculationResult(rows, tuple(str(item.message) for item in captured)),
        raw,
        edge,
        separation,
    )


def profile_rows(calculation: MarchingCalculation, station: int) -> tuple[Row, ...]:
    """Return the wall-normal profile at one solved station in SI units."""
    raw = calculation.raw
    if not 0 <= station < raw.distance.size:
        raise ValueError("station is outside the solved range")
    temperature = float(calculation.edge.temperature[station + 1])
    velocity = float(calculation.edge.velocity[station + 1])
    rows: list[Row] = []
    for node in range(raw.eta.size):
        row: Row = {
            "wall_distance": float(raw.wall_distance[station, node]),
            "eta": float(raw.eta[node]),
            "velocity_ratio": float(raw.velocity_ratio[station, node]),
            "temperature_ratio": float(raw.temperature_ratio[station, node]),
            "velocity": float(raw.velocity_ratio[station, node]) * velocity,
            "temperature": float(raw.temperature_ratio[station, node]) * temperature,
        }
        if raw.nu_tilde_ratio is not None:
            row["nu_tilde_ratio"] = float(raw.nu_tilde_ratio[station, node])
        rows.append(row)
    return tuple(rows)
