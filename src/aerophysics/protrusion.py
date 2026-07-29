"""Direct drag of an isolated protrusion immersed in a boundary layer.

The engineering model in this module scales a separately supplied free-stream
drag coefficient by the dynamic pressure averaged over the protrusion frontal
area.  It estimates direct protrusion drag only; wall interference, transition,
and downstream drag increments are outside the model.

References
----------
Young, A. D. and Paterson, J. H., *Aircraft Excrescence Drag*,
AGARD-AG-264, 1981.
Hoerner, S. F., *Fluid-Dynamic Drag*, 1965, Chapter 5.
Johnson, D. F. and Mitchell, G. A., *Experimental Investigation of Two
Methods for Generating an Artificially Thickened Boundary Layer*,
NASA TM X-2238, 1971.
"""

import warnings
from dataclasses import dataclass
from enum import StrEnum
from typing import Protocol

import numpy as np
from numpy.typing import ArrayLike

from aerophysics._array import FloatArray, as_float_array
from aerophysics.exceptions import ApplicabilityWarning
from aerophysics.gas import AIR, PerfectGas


class FrontalWidthFunction(Protocol):
    """Callable returning projected protrusion width at supplied heights."""

    def __call__(self, height: FloatArray, /) -> ArrayLike: ...


class ProtrusionProfile(StrEnum):
    """Boundary-layer profile used for the effective-dynamic-pressure integral."""

    TURBULENT_ONE_SEVENTH_POWER = "turbulent_one_seventh_power"
    PROVIDED = "provided"


@dataclass(frozen=True, slots=True)
class ProtrusionDragResult:
    """Direct protrusion drag and effective-dynamic-pressure metadata."""

    direct_drag: float
    effective_dynamic_pressure: float
    shielding_factor: float
    frontal_area: float
    edge_dynamic_pressure: float
    height_to_boundary_layer_thickness: float
    profile: ProtrusionProfile
    compressibility_applied: bool


def _positive_scalar(value: float, *, name: str, allow_zero: bool = False) -> float:
    array, scalar = as_float_array(value, name=name)
    if not scalar:
        raise ValueError(f"{name} must be a scalar")
    result = float(array)
    if result < 0.0 if allow_zero else result <= 0.0:
        qualifier = "non-negative" if allow_zero else "greater than zero"
        raise ValueError(f"{name} must be {qualifier}")
    return result


def _integration_grid(
    height: float, boundary_layer_thickness: float, points: int
) -> FloatArray:
    immersed_height = min(height, boundary_layer_thickness)
    immersed = np.linspace(0.0, immersed_height, points, dtype=np.float64)
    if height <= boundary_layer_thickness:
        return immersed
    outer = np.linspace(
        boundary_layer_thickness, height, points, dtype=np.float64
    )[1:]
    return np.concatenate((immersed, outer))


def _frontal_width_values(
    width: float | FrontalWidthFunction, height: FloatArray
) -> FloatArray:
    if callable(width):
        try:
            raw_width = width(height)
        except Exception as error:
            raise ValueError("frontal_width could not be evaluated") from error
    else:
        raw_width = width
    values, _ = as_float_array(raw_width, name="frontal_width")
    try:
        broadcast = np.broadcast_to(values, height.shape)
    except ValueError as error:
        raise ValueError(
            "frontal_width must return a scalar or values matching profile height"
        ) from error
    result = np.asarray(broadcast, dtype=np.float64)
    if np.any(result < 0.0):
        raise ValueError("frontal_width must be non-negative")
    return result


def _provided_profile(
    evaluation_height: FloatArray,
    *,
    protrusion_height: float,
    boundary_layer_thickness: float,
    edge_velocity: float,
    edge_density: float,
    profile_height: ArrayLike,
    profile_velocity: ArrayLike,
    profile_density: ArrayLike,
) -> tuple[FloatArray, FloatArray]:
    sample_height, _ = as_float_array(profile_height, name="profile_height")
    sample_velocity, _ = as_float_array(profile_velocity, name="profile_velocity")
    sample_density, _ = as_float_array(profile_density, name="profile_density")
    if (
        sample_height.ndim != 1
        or sample_velocity.ndim != 1
        or sample_density.ndim != 1
    ):
        raise ValueError("provided profile arrays must be one-dimensional")
    if not (
        sample_height.size == sample_velocity.size == sample_density.size
        and sample_height.size >= 2
    ):
        raise ValueError(
            "provided profile arrays must have the same length of at least two"
        )
    if sample_height[0] != 0.0 or np.any(np.diff(sample_height) <= 0.0):
        raise ValueError(
            "profile_height must start at zero and be strictly increasing"
        )
    required_height = min(protrusion_height, boundary_layer_thickness)
    if sample_height[-1] < required_height:
        raise ValueError(
            "provided profile must cover the immersed protrusion height"
        )
    if np.any(sample_velocity < 0.0):
        raise ValueError("profile_velocity must be non-negative")
    if np.any(sample_density <= 0.0):
        raise ValueError("profile_density must be greater than zero")

    immersed = evaluation_height <= boundary_layer_thickness
    velocity = np.full_like(evaluation_height, edge_velocity)
    density = np.full_like(evaluation_height, edge_density)
    velocity[immersed] = np.interp(
        evaluation_height[immersed], sample_height, sample_velocity
    )
    density[immersed] = np.interp(
        evaluation_height[immersed], sample_height, sample_density
    )
    return velocity, density


def _one_seventh_power_profile(
    height: FloatArray,
    *,
    edge_velocity: float,
    edge_density: float,
    boundary_layer_thickness: float,
    mach: float | None,
    edge_temperature: float | None,
    wall_temperature: float | None,
    prandtl_number: float,
    gas: PerfectGas,
) -> tuple[FloatArray, FloatArray, bool]:
    velocity_ratio = np.minimum(height / boundary_layer_thickness, 1.0) ** (
        1.0 / 7.0
    )
    velocity = edge_velocity * velocity_ratio
    if mach is None:
        return velocity, np.full_like(height, edge_density), False

    assert edge_temperature is not None
    recovery_factor = np.cbrt(prandtl_number)
    recovery_temperature = edge_temperature * (
        1.0
        + recovery_factor
        * 0.5
        * (gas.heat_capacity_ratio - 1.0)
        * mach**2
    )
    wall = recovery_temperature if wall_temperature is None else wall_temperature
    temperature = (
        wall
        + (recovery_temperature - wall) * velocity_ratio
        + (edge_temperature - recovery_temperature) * velocity_ratio**2
    )
    density = edge_density * edge_temperature / temperature
    density[height > boundary_layer_thickness] = edge_density
    return velocity, np.asarray(density, dtype=np.float64), True


def protrusion_drag(
    drag_coefficient: float,
    height: float,
    frontal_width: float | FrontalWidthFunction,
    edge_velocity: float,
    edge_density: float,
    boundary_layer_thickness: float,
    *,
    profile_height: ArrayLike | None = None,
    profile_velocity: ArrayLike | None = None,
    profile_density: ArrayLike | None = None,
    mach: float | None = None,
    edge_temperature: float | None = None,
    wall_temperature: float | None = None,
    prandtl_number: float = 0.72,
    gas: PerfectGas = AIR,
    integration_points: int = 4097,
) -> ProtrusionDragResult:
    """Estimate the direct drag of one boundary-layer-immersed protrusion.

    Parameters
    ----------
    drag_coefficient:
        Free-stream drag coefficient based on protrusion frontal area.
    height:
        Protrusion height normal to the wall in m.
    frontal_width:
        Constant projected width in m, or a callable ``b(y)`` returning the
        projected width at heights in m.
    edge_velocity, edge_density:
        Boundary-layer edge velocity in m/s and density in kg/m³.
    boundary_layer_thickness:
        Boundary-layer thickness in m.
    profile_height, profile_velocity, profile_density:
        Optional measured or computed one-dimensional profile. All three must
        be supplied together. Height is in m, velocity in m/s, and density in
        kg/m³. The height array must start at zero and cover
        ``min(height, boundary_layer_thickness)``.
    mach, edge_temperature:
        Optional edge Mach number and temperature in K for the compressible
        Walz-profile approximation. They must be supplied together and cannot
        be combined with a provided profile.
    wall_temperature:
        Wall temperature in K. If omitted in compressible mode, the wall is
        adiabatic at the turbulent recovery temperature.
    prandtl_number:
        Constant Prandtl number used for turbulent recovery temperature.
    gas:
        Perfect gas supplying the heat-capacity ratio.
    integration_points:
        Points used in each of the immersed and outer integration intervals.

    Notes
    -----
    The returned force is ``C_D`` times the local dynamic pressure integrated
    over frontal area. It excludes wall interference, horseshoe vortices,
    transition, downstream skin-friction changes, and interactions between
    multiple protrusions. A scalar free-stream ``C_D`` is especially uncertain
    in transonic or shock/boundary-layer-interaction regimes.
    """
    coefficient = _positive_scalar(
        drag_coefficient, name="drag_coefficient", allow_zero=True
    )
    protrusion_height = _positive_scalar(height, name="height")
    velocity_edge = _positive_scalar(edge_velocity, name="edge_velocity")
    density_edge = _positive_scalar(edge_density, name="edge_density")
    thickness = _positive_scalar(
        boundary_layer_thickness, name="boundary_layer_thickness"
    )
    prandtl = _positive_scalar(prandtl_number, name="prandtl_number")
    if (
        isinstance(integration_points, bool)
        or not isinstance(integration_points, int)
        or integration_points < 32
    ):
        raise ValueError("integration_points must be an integer of at least 32")

    profile_inputs = (profile_height, profile_velocity, profile_density)
    has_profile = all(value is not None for value in profile_inputs)
    if any(value is not None for value in profile_inputs) and not has_profile:
        raise ValueError(
            "profile_height, profile_velocity, and profile_density "
            "must be supplied together"
        )

    thermal_inputs = mach is not None or edge_temperature is not None
    if thermal_inputs and (mach is None or edge_temperature is None):
        raise ValueError("mach and edge_temperature must be supplied together")
    if has_profile and (
        thermal_inputs or wall_temperature is not None
    ):
        raise ValueError(
            "thermal inputs cannot be combined with a provided density profile"
        )
    if not thermal_inputs and wall_temperature is not None:
        raise ValueError(
            "wall_temperature requires mach and edge_temperature"
        )

    mach_value: float | None = None
    edge_temperature_value: float | None = None
    wall_temperature_value: float | None = None
    if mach is not None:
        mach_value = _positive_scalar(mach, name="mach", allow_zero=True)
        assert edge_temperature is not None
        edge_temperature_value = _positive_scalar(
            edge_temperature, name="edge_temperature"
        )
        if wall_temperature is not None:
            wall_temperature_value = _positive_scalar(
                wall_temperature, name="wall_temperature"
            )
        if 0.8 <= mach_value <= 1.4:
            warnings.warn(
                "a scalar free-stream drag coefficient may not capture "
                "transonic shock and boundary-layer interaction",
                ApplicabilityWarning,
                stacklevel=2,
            )

    grid = _integration_grid(protrusion_height, thickness, integration_points)
    width = _frontal_width_values(frontal_width, grid)
    frontal_area = float(np.trapezoid(width, grid))
    if frontal_area <= 0.0:
        raise ValueError("frontal_width must define a positive frontal area")

    if has_profile:
        assert profile_height is not None
        assert profile_velocity is not None
        assert profile_density is not None
        velocity, density = _provided_profile(
            grid,
            protrusion_height=protrusion_height,
            boundary_layer_thickness=thickness,
            edge_velocity=velocity_edge,
            edge_density=density_edge,
            profile_height=profile_height,
            profile_velocity=profile_velocity,
            profile_density=profile_density,
        )
        profile = ProtrusionProfile.PROVIDED
        compressibility_applied = False
    else:
        velocity, density, compressibility_applied = _one_seventh_power_profile(
            grid,
            edge_velocity=velocity_edge,
            edge_density=density_edge,
            boundary_layer_thickness=thickness,
            mach=mach_value,
            edge_temperature=edge_temperature_value,
            wall_temperature=wall_temperature_value,
            prandtl_number=prandtl,
            gas=gas,
        )
        profile = ProtrusionProfile.TURBULENT_ONE_SEVENTH_POWER

    dynamic_pressure = 0.5 * density * velocity**2
    dynamic_pressure_area_integral = float(
        np.trapezoid(dynamic_pressure * width, grid)
    )
    effective_dynamic_pressure = dynamic_pressure_area_integral / frontal_area
    edge_dynamic_pressure = 0.5 * density_edge * velocity_edge**2
    shielding_factor = effective_dynamic_pressure / edge_dynamic_pressure

    return ProtrusionDragResult(
        direct_drag=coefficient * dynamic_pressure_area_integral,
        effective_dynamic_pressure=effective_dynamic_pressure,
        shielding_factor=shielding_factor,
        frontal_area=frontal_area,
        edge_dynamic_pressure=edge_dynamic_pressure,
        height_to_boundary_layer_thickness=protrusion_height / thickness,
        profile=profile,
        compressibility_applied=compressibility_applied,
    )


__all__ = [
    "FrontalWidthFunction",
    "ProtrusionDragResult",
    "ProtrusionProfile",
    "protrusion_drag",
]
