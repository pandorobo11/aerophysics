"""Pure SI boundaries for the axisymmetric body-flow calculator."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from aerophysics.atmosphere import geopotential_to_geometric, standard_atmosphere
from aerophysics.characteristics import (
    AxisymmetricBodyFlowResult,
    axisymmetric_body_flow,
)
from aerophysics.gas import AIR
from aerophysics.gui.adapters import CalculationResult, Row
from aerophysics.gui.csv_inputs import MeridianCSV


def body_flow_inputs_match(
    current: Mapping[str, object], saved: Mapping[str, object]
) -> bool:
    """Ignore SI roundoff from display conversions when checking saved inputs."""

    def matches(left: object, right: object) -> bool:
        if isinstance(left, float) and isinstance(right, float):
            return math.isclose(left, right, rel_tol=1e-12, abs_tol=0.0)
        if isinstance(left, list) and isinstance(right, list):
            return len(left) == len(right) and all(
                matches(a, b) for a, b in zip(left, right, strict=True)
            )
        return type(left) is type(right) and left == right

    return current.keys() == saved.keys() and all(
        matches(value, saved[key]) for key, value in current.items()
    )


def build_meridian(
    *,
    nose: str = "cone",
    nose_length: float = 0.3,
    base_radius: float = 0.05,
    ogive_segments: int = 32,
    cylinder_length: float = 1.0,
    boattail: bool = False,
    boattail_length: float | None = None,
    boattail_radius: float | None = None,
) -> MeridianCSV:
    """Build sharp convex meridian vertices in metres, including exact endpoints."""
    if not all(math.isfinite(v) and v > 0 for v in (nose_length, base_radius)):
        raise ValueError("nose_length and base_radius must be finite and positive")
    if not math.isfinite(cylinder_length) or cylinder_length < 0:
        raise ValueError("cylinder_length must be finite and non-negative")
    if nose == "cone":
        x = [0.0, nose_length]
        r = [0.0, base_radius]
    elif nose == "tangent_ogive":
        if base_radius > nose_length:
            raise ValueError("tangent ogive requires base_radius <= nose_length")
        if (
            isinstance(ogive_segments, bool)
            or not isinstance(ogive_segments, int)
            or not 2 <= ogive_segments <= 200
        ):
            raise ValueError("ogive_segments must be an integer from 2 to 200")
        rho = (base_radius**2 + nose_length**2) / (2 * base_radius)
        phi = np.linspace(math.asin(nose_length / rho), 0.0, ogive_segments + 1)
        x = (nose_length - rho * np.sin(phi)).tolist()
        r = (rho * np.cos(phi) + base_radius - rho).tolist()
        x[0], r[0] = 0.0, 0.0
        x[-1], r[-1] = nose_length, base_radius
    else:
        raise ValueError("unsupported nose type")
    if cylinder_length > 0:
        x.append(nose_length + cylinder_length)
        r.append(base_radius)
    if boattail:
        if (
            boattail_length is None
            or not math.isfinite(boattail_length)
            or boattail_length <= 0
        ):
            raise ValueError("boattail_length must be finite and positive")
        if (
            boattail_radius is None
            or not math.isfinite(boattail_radius)
            or not 0 < boattail_radius < base_radius
        ):
            raise ValueError(
                "boattail_radius must satisfy 0 < boattail_radius < base_radius"
            )
        x.append(x[-1] + boattail_length)
        r.append(boattail_radius)
    return MeridianCSV(np.asarray(x), np.asarray(r))


def resolve_freestream(
    *,
    freestream_basis: str = "static",
    pressure: float = 101325.0,
    temperature: float = 288.15,
    altitude: float | None = None,
    altitude_basis: str = "geometric",
) -> tuple[float, float]:
    """Resolve static pressure [Pa] and temperature [K] from the selected basis."""
    if freestream_basis == "static":
        return pressure, temperature
    if freestream_basis != "atmosphere":
        raise ValueError("unsupported freestream basis")
    if altitude is None:
        raise ValueError("altitude is required for atmosphere mode")
    if altitude_basis == "geopotential":
        altitude = float(geopotential_to_geometric(altitude))
    elif altitude_basis != "geometric":
        raise ValueError("unsupported altitude basis")
    state = standard_atmosphere(altitude)
    return float(state.pressure), float(state.temperature)


@dataclass(frozen=True, slots=True)
class BodyFlowCalculation:
    """Surface rows, shock rows and the unmodified public solver result."""

    surface: CalculationResult
    shock_rows: tuple[Row, ...]
    raw: AxisymmetricBodyFlowResult


def calculate_body_flow(
    meridian: MeridianCSV,
    *,
    mach: float = 3.0,
    pressure: float = 101325.0,
    temperature: float = 288.15,
    characteristic_count: int = 88,
    fan_angle_step: float = math.pi / 180,
) -> BodyFlowCalculation:
    """Call the AIR solver once and derive dimensionless surface pressures."""
    raw = axisymmetric_body_flow(
        meridian.axial_distance,
        meridian.radius,
        mach=mach,
        pressure=pressure,
        temperature=temperature,
        gas=AIR,
        characteristic_count=characteristic_count,
        fan_angle_step=fan_angle_step,
    )
    values = {
        "axial_distance": raw.axial_distance,
        "radius": raw.radius,
        "surface_distance": raw.surface_distance,
        "pressure": raw.pressure,
        "pressure_ratio": raw.pressure / pressure,
        "pressure_rise": (raw.pressure - pressure) / pressure,
        "cp": (raw.pressure - pressure)
        / (0.5 * AIR.heat_capacity_ratio * pressure * mach**2),
        "temperature": raw.temperature,
        "mach": raw.mach,
        "flow_angle": raw.flow_angle,
    }
    rows: tuple[Row, ...] = tuple(
        {key: float(value[i]) for key, value in values.items()}
        for i in range(len(raw.axial_distance))
    )
    shock: tuple[Row, ...] = tuple(
        {
            "shock_axial_distance": float(x),
            "shock_radius": float(r),
            "shock_angle": float(beta),
        }
        for x, r, beta in zip(
            raw.shock_axial_distance, raw.shock_radius, raw.shock_angle, strict=True
        )
    )
    return BodyFlowCalculation(CalculationResult(rows), shock, raw)
