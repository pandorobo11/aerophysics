"""Compressible turbulent boundary-layer mean-property profiles.

This module combines an incompressible Spalding--Coles composite velocity
profile with either the Van Driest or Volpiani compressibility
transformation.  The inverse model is intended for smooth, fully turbulent,
approximately zero-pressure-gradient boundary layers.

All dimensional inputs and outputs use SI units.  ``boundary_layer_thickness``
is the conventional 99-percent thickness.

References
----------
Van Driest, E. R., *Turbulent Boundary Layer in Compressible Fluids*,
Journal of the Aeronautical Sciences, 18(3), 1951.
Volpiani, P. S., Iyer, P. S., Pirozzoli, S., and Larsson, J.,
*Data-Driven Compressibility Transformation for Turbulent Wall Layers*,
Physical Review Fluids, 5, 052602(R), 2020.
Coles, D., *The Law of the Wake in the Turbulent Boundary Layer*,
Journal of Fluid Mechanics, 1(2), 191--226, 1956.
Spalding, D. B., *A Single Formula for the Law of the Wall*,
Journal of Applied Mechanics, 28(3), 455--458, 1961.
Zhang, Y.-S., Bi, W.-T., Hussain, F., and She, Z.-S.,
*A Generalized Reynolds Analogy for Compressible Wall-Bounded Turbulent
Flows*, Journal of Fluid Mechanics, 739, 392--420, 2014.
"""

import warnings
from dataclasses import dataclass
from enum import StrEnum

import numpy as np
from numpy.typing import ArrayLike
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.optimize import brentq

from aerophysics._array import FloatArray, as_float_array
from aerophysics.exceptions import ApplicabilityWarning, ModelRangeError
from aerophysics.gas import AIR, AIR_VISCOSITY, PerfectGas, SutherlandModel


class CompressibleVelocityTransformation(StrEnum):
    """Compressible mean-velocity transformation used by a profile model."""

    VAN_DRIEST = "van_driest"
    VOLPIANI = "volpiani"


class TemperatureVelocityRelation(StrEnum):
    """Algebraic temperature--velocity relation used by the inverse model."""

    GENERALIZED_REYNOLDS_ANALOGY = "generalized_reynolds_analogy"
    WALZ = "walz"


@dataclass(frozen=True, slots=True)
class TransformedVelocityProfileResult:
    """A supplied compressible profile mapped to incompressible variables."""

    transformation: CompressibleVelocityTransformation
    wall_distance: FloatArray
    wall_distance_plus: FloatArray
    velocity_plus: FloatArray
    transformed_wall_coordinate: FloatArray
    transformed_velocity_plus: FloatArray
    friction_velocity: float


@dataclass(frozen=True, slots=True)
class CompressibleBoundaryLayerProfileResult:
    """Predicted mean properties from the wall through ``delta_99``.

    The displacement and momentum thicknesses are integrated only through
    ``boundary_layer_thickness``.  They therefore neglect the small remaining
    defect outside the 99-percent thickness.
    """

    transformation: CompressibleVelocityTransformation
    temperature_velocity_relation: TemperatureVelocityRelation
    wall_distance: FloatArray
    wall_distance_plus: FloatArray
    transformed_wall_coordinate: FloatArray
    velocity: FloatArray
    velocity_plus: FloatArray
    transformed_velocity_plus: FloatArray
    temperature: FloatArray
    density: FloatArray
    dynamic_viscosity: FloatArray
    local_mach_number: FloatArray
    dynamic_pressure: FloatArray
    friction_velocity: float
    friction_reynolds_number: float
    recovery_temperature: float
    wall_temperature: float
    wake_parameter: float
    edge_velocity_ratio: float
    displacement_thickness: float
    momentum_thickness: float
    shape_factor: float
    local_skin_friction_coefficient: float


@dataclass(frozen=True, slots=True)
class _InverseParameters:
    transformation: CompressibleVelocityTransformation
    temperature_velocity_relation: TemperatureVelocityRelation
    edge_velocity: float
    edge_density: float
    edge_temperature: float
    wall_temperature: float
    recovery_temperature: float
    wall_density: float
    wall_viscosity: float
    friction_velocity: float
    friction_reynolds_number: float
    prandtl_number: float
    reynolds_analogy_factor: float
    von_karman_constant: float
    log_law_intercept: float
    gas: PerfectGas
    viscosity_model: SutherlandModel


def _positive_scalar(value: float, *, name: str) -> float:
    array, scalar = as_float_array(value, name=name)
    if not scalar:
        raise ValueError(f"{name} must be a scalar")
    result = float(array)
    if result <= 0.0:
        raise ValueError(f"{name} must be greater than zero")
    return result


def _validate_transformation(
    transformation: CompressibleVelocityTransformation,
) -> None:
    if not isinstance(transformation, CompressibleVelocityTransformation):
        raise ValueError(
            "transformation must be a CompressibleVelocityTransformation value"
        )


def _validate_temperature_relation(
    relation: TemperatureVelocityRelation,
) -> None:
    if not isinstance(relation, TemperatureVelocityRelation):
        raise ValueError(
            "temperature_velocity_relation must be a TemperatureVelocityRelation value"
        )


def _profile_array(value: ArrayLike, *, name: str) -> FloatArray:
    array, _ = as_float_array(value, name=name)
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    return np.asarray(array, dtype=np.float64)


def _validate_wall_grid(wall_distance: FloatArray, *, minimum_size: int = 2) -> None:
    if wall_distance.size < minimum_size:
        raise ValueError(f"wall_distance must contain at least {minimum_size} points")
    if wall_distance[0] != 0.0:
        raise ValueError("wall_distance must start at zero")
    if np.any(np.diff(wall_distance) <= 0.0):
        raise ValueError("wall_distance must be strictly increasing")


def _spalding_wall_coordinate(
    velocity_plus: FloatArray, *, von_karman_constant: float, intercept: float
) -> FloatArray:
    scaled = von_karman_constant * velocity_plus
    remainder = np.expm1(scaled) - scaled - 0.5 * scaled**2 - scaled**3 / 6.0
    return np.asarray(
        velocity_plus + np.exp(-von_karman_constant * intercept) * remainder,
        dtype=np.float64,
    )


def _spalding_coordinate_derivative(
    velocity_plus: float, *, von_karman_constant: float, intercept: float
) -> float:
    scaled = von_karman_constant * velocity_plus
    remainder_derivative = np.expm1(scaled) - scaled - 0.5 * scaled**2
    return float(
        1.0
        + np.exp(-von_karman_constant * intercept)
        * von_karman_constant
        * remainder_derivative
    )


def _wake_function(outer_coordinate: FloatArray) -> FloatArray:
    return np.asarray(
        2.0 * np.sin(0.5 * np.pi * outer_coordinate) ** 2,
        dtype=np.float64,
    )


def _temperature_profile(
    velocity_ratio: FloatArray,
    *,
    relation: TemperatureVelocityRelation,
    edge_temperature: float,
    wall_temperature: float,
    recovery_temperature: float,
    prandtl_number: float,
    reynolds_analogy_factor: float,
) -> FloatArray:
    ratio = np.asarray(velocity_ratio, dtype=np.float64)
    if relation is TemperatureVelocityRelation.WALZ:
        temperature = (
            wall_temperature
            + (recovery_temperature - wall_temperature) * ratio
            + (edge_temperature - recovery_temperature) * ratio**2
        )
    else:
        temperature = (
            wall_temperature
            + reynolds_analogy_factor
            * prandtl_number
            * (recovery_temperature - wall_temperature)
            * ratio
            * (1.0 - ratio)
            + (edge_temperature - wall_temperature) * ratio**2
        )
    result = np.asarray(temperature, dtype=np.float64)
    if np.any(~np.isfinite(result) | (result <= 0.0)):
        raise ModelRangeError(
            "temperature--velocity relation produced a non-physical temperature"
        )
    return result


def _property_ratios(
    physical_velocity_plus: float,
    parameters: _InverseParameters,
) -> tuple[float, float]:
    velocity = parameters.friction_velocity * physical_velocity_plus
    ratio = np.asarray([velocity / parameters.edge_velocity], dtype=np.float64)
    temperature = float(
        _temperature_profile(
            ratio,
            relation=parameters.temperature_velocity_relation,
            edge_temperature=parameters.edge_temperature,
            wall_temperature=parameters.wall_temperature,
            recovery_temperature=parameters.recovery_temperature,
            prandtl_number=parameters.prandtl_number,
            reynolds_analogy_factor=parameters.reynolds_analogy_factor,
        )[0]
    )
    density = parameters.edge_density * parameters.edge_temperature / temperature
    viscosity = float(parameters.viscosity_model.dynamic_viscosity(temperature))
    return (
        density / parameters.wall_density,
        viscosity / parameters.wall_viscosity,
    )


def _mapping_factors(
    physical_velocity_plus: float,
    parameters: _InverseParameters,
) -> tuple[float, float]:
    density_ratio, viscosity_ratio = _property_ratios(
        physical_velocity_plus, parameters
    )
    density_root = np.sqrt(density_ratio)
    if parameters.transformation is CompressibleVelocityTransformation.VAN_DRIEST:
        return 1.0, float(density_root)
    coordinate_factor = density_root * viscosity_ratio**-1.5
    velocity_factor = density_root * viscosity_ratio**-0.5
    return float(coordinate_factor), float(velocity_factor)


def _integrate_inverse_profile(
    wake_parameter: float,
    evaluation_wall_distance_plus: FloatArray,
    parameters: _InverseParameters,
) -> tuple[FloatArray, FloatArray]:
    reynolds = parameters.friction_reynolds_number
    kappa = parameters.von_karman_constant
    intercept = parameters.log_law_intercept

    def right_hand_side(wall_distance_plus: float, state: FloatArray) -> FloatArray:
        physical_velocity_plus = float(state[0])
        wall_velocity_plus = float(state[1])
        coordinate_factor, velocity_factor = _mapping_factors(
            physical_velocity_plus, parameters
        )
        spalding_derivative = _spalding_coordinate_derivative(
            wall_velocity_plus,
            von_karman_constant=kappa,
            intercept=intercept,
        )
        wall_velocity_gradient = coordinate_factor / spalding_derivative
        outer_coordinate = wall_distance_plus / reynolds
        wake_gradient = (
            wake_parameter
            * np.pi
            * np.sin(np.pi * outer_coordinate)
            / (kappa * reynolds)
        )
        physical_velocity_gradient = (
            wall_velocity_gradient + wake_gradient
        ) / velocity_factor
        return np.asarray(
            [physical_velocity_gradient, wall_velocity_gradient],
            dtype=np.float64,
        )

    solution = solve_ivp(
        right_hand_side,
        (0.0, reynolds),
        np.zeros(2, dtype=np.float64),
        t_eval=evaluation_wall_distance_plus,
        rtol=1e-9,
        atol=1e-11,
    )
    if not solution.success:
        raise ModelRangeError(
            f"inverse velocity integration failed: {solution.message}"
        )
    values = np.asarray(solution.y, dtype=np.float64)
    if values.shape != (2, evaluation_wall_distance_plus.size):
        raise ModelRangeError("inverse velocity integration returned incomplete data")
    if np.any(~np.isfinite(values)):
        raise ModelRangeError("inverse velocity integration returned non-finite data")
    return values[0], values[1]


def _edge_velocity_ratio(
    wake_parameter: float, parameters: _InverseParameters
) -> float:
    endpoint = np.asarray([parameters.friction_reynolds_number], dtype=np.float64)
    velocity_plus, _ = _integrate_inverse_profile(wake_parameter, endpoint, parameters)
    return float(
        velocity_plus[-1] * parameters.friction_velocity / parameters.edge_velocity
    )


def _choose_wake_parameter(
    specified_wake_parameter: float | None,
    parameters: _InverseParameters,
) -> tuple[float, float]:
    target = 0.99
    tolerance = 1e-3
    if specified_wake_parameter is not None:
        array, scalar = as_float_array(specified_wake_parameter, name="wake_parameter")
        if not scalar:
            raise ValueError("wake_parameter must be a scalar")
        wake = float(array)
        if wake < 0.0:
            raise ValueError("wake_parameter must be non-negative")
        ratio = _edge_velocity_ratio(wake, parameters)
        if abs(ratio - target) > tolerance:
            raise ModelRangeError(
                "specified wake_parameter does not satisfy "
                "U(delta_99) / U_e = 0.99 within 1e-3"
            )
        return wake, ratio

    residual_at_zero = _edge_velocity_ratio(0.0, parameters) - target
    residual_at_one = _edge_velocity_ratio(1.0, parameters) - target
    if abs(residual_at_zero) <= 1e-10:
        return 0.0, residual_at_zero + target
    if abs(residual_at_one) <= 1e-10:
        return 1.0, residual_at_one + target
    if residual_at_zero * residual_at_one > 0.0:
        raise ModelRangeError(
            "no Coles wake parameter in 0 <= Pi <= 1 satisfies "
            "U(delta_99) / U_e = 0.99; edge velocity, wall shear stress, "
            "and boundary-layer thickness are mutually inconsistent"
        )
    wake = float(
        brentq(
            lambda value: _edge_velocity_ratio(value, parameters) - target,
            0.0,
            1.0,
            xtol=1e-10,
            rtol=1e-10,
        )
    )
    return wake, _edge_velocity_ratio(wake, parameters)


def _integration_grid(
    requested_wall_distance_plus: FloatArray, friction_reynolds_number: float
) -> FloatArray:
    outer = np.linspace(0.0, friction_reynolds_number, 2049, dtype=np.float64)
    first_positive = max(
        np.finfo(np.float64).eps,
        friction_reynolds_number * 1e-12,
    )
    inner = np.geomspace(
        first_positive,
        friction_reynolds_number,
        2048,
        dtype=np.float64,
    )
    return np.unique(np.concatenate((outer, inner, requested_wall_distance_plus)))


def transform_compressible_velocity_profile(
    wall_distance: ArrayLike,
    velocity: ArrayLike,
    density: ArrayLike,
    dynamic_viscosity: ArrayLike,
    wall_shear_stress: float,
    *,
    transformation: CompressibleVelocityTransformation,
) -> TransformedVelocityProfileResult:
    """Map a known compressible mean-velocity profile to wall variables.

    Parameters
    ----------
    wall_distance:
        One-dimensional wall-normal coordinates in m, starting at zero.
    velocity:
        Favre mean streamwise velocity in m/s. Values must start at zero and
        be non-negative and non-decreasing.
    density:
        Mean density in kg/m³.
    dynamic_viscosity:
        Mean dynamic viscosity in Pa s.
    wall_shear_stress:
        Positive wall shear stress in Pa.
    transformation:
        Van Driest or Volpiani transformation.
    """
    _validate_transformation(transformation)
    height = _profile_array(wall_distance, name="wall_distance")
    speed = _profile_array(velocity, name="velocity")
    density_values = _profile_array(density, name="density")
    viscosity_values = _profile_array(dynamic_viscosity, name="dynamic_viscosity")
    _validate_wall_grid(height)
    if not (height.size == speed.size == density_values.size == viscosity_values.size):
        raise ValueError("profile arrays must have the same length")
    if speed[0] != 0.0:
        raise ValueError("velocity must start at zero at the wall")
    if np.any(speed < 0.0) or np.any(np.diff(speed) < 0.0):
        raise ValueError("velocity must be non-negative and non-decreasing")
    if np.any(density_values <= 0.0):
        raise ValueError("density must be greater than zero")
    if np.any(viscosity_values <= 0.0):
        raise ValueError("dynamic_viscosity must be greater than zero")
    shear = _positive_scalar(wall_shear_stress, name="wall_shear_stress")

    wall_density = float(density_values[0])
    wall_viscosity = float(viscosity_values[0])
    friction_velocity = np.sqrt(shear / wall_density)
    wall_distance_plus = wall_density * friction_velocity * height / wall_viscosity
    velocity_plus = speed / friction_velocity
    density_ratio = density_values / wall_density
    viscosity_ratio = viscosity_values / wall_viscosity

    if transformation is CompressibleVelocityTransformation.VAN_DRIEST:
        transformed_coordinate = wall_distance_plus.copy()
        velocity_integrand = np.sqrt(density_ratio)
    else:
        transformed_coordinate = cumulative_trapezoid(
            np.sqrt(density_ratio) * viscosity_ratio**-1.5,
            x=wall_distance_plus,
            initial=0,
        )
        velocity_integrand = np.sqrt(density_ratio) * viscosity_ratio**-0.5
    transformed_velocity = cumulative_trapezoid(
        velocity_integrand,
        x=velocity_plus,
        initial=0,
    )
    return TransformedVelocityProfileResult(
        transformation=transformation,
        wall_distance=height.copy(),
        wall_distance_plus=np.asarray(wall_distance_plus, dtype=np.float64),
        velocity_plus=np.asarray(velocity_plus, dtype=np.float64),
        transformed_wall_coordinate=np.asarray(
            transformed_coordinate, dtype=np.float64
        ),
        transformed_velocity_plus=np.asarray(transformed_velocity, dtype=np.float64),
        friction_velocity=float(friction_velocity),
    )


def compressible_turbulent_boundary_layer_profile(
    wall_distance: ArrayLike,
    edge_velocity: float,
    edge_density: float,
    edge_temperature: float,
    boundary_layer_thickness: float,
    wall_shear_stress: float,
    *,
    transformation: CompressibleVelocityTransformation,
    wall_temperature: float | None = None,
    temperature_velocity_relation: TemperatureVelocityRelation = (
        TemperatureVelocityRelation.GENERALIZED_REYNOLDS_ANALOGY
    ),
    wake_parameter: float | None = None,
    prandtl_number: float = 0.72,
    reynolds_analogy_factor: float = 1.14,
    von_karman_constant: float = 0.41,
    log_law_intercept: float = 5.2,
    gas: PerfectGas = AIR,
    viscosity_model: SutherlandModel = AIR_VISCOSITY,
) -> CompressibleBoundaryLayerProfileResult:
    """Predict a compressible Spalding--Coles mean-property profile.

    Parameters
    ----------
    wall_distance:
        Strictly increasing one-dimensional output coordinates in m, starting
        at zero and not exceeding ``boundary_layer_thickness``.
    edge_velocity, edge_density, edge_temperature:
        Boundary-layer-edge velocity in m/s, density in kg/m³, and
        temperature in K.
    boundary_layer_thickness:
        Conventional 99-percent boundary-layer thickness in m.
    wall_shear_stress:
        Positive wall shear stress in Pa.
    transformation:
        Van Driest or Volpiani inverse velocity transformation.
    wall_temperature:
        Wall temperature in K. If omitted, the turbulent recovery temperature
        is used.
    temperature_velocity_relation:
        Generalized Reynolds analogy (default) or Walz relation.
    wake_parameter:
        Optional Coles wake strength. If omitted, it is solved in the range
        zero to one so that ``U(delta_99) / U_e = 0.99``.
    prandtl_number:
        Constant Prandtl number used by the recovery and temperature models.
    reynolds_analogy_factor:
        Coefficient ``s`` in the generalized Reynolds analogy. It is ignored
        by the Walz relation.
    von_karman_constant, log_law_intercept:
        Constants in the Spalding law of the wall and Coles wake amplitude.
    gas, viscosity_model:
        Perfect-gas and Sutherland models used for thermodynamic properties.

    Notes
    -----
    The Coles wake is added in transformed velocity and uses physical
    ``y / delta_99``. This is an engineering coupling of incompressible and
    compressible models, not a new consequence of either transformation.
    Pressure-gradient, rough-wall, transitional, and real-gas flows are
    outside the implemented scope.
    """
    _validate_transformation(transformation)
    _validate_temperature_relation(temperature_velocity_relation)
    height = _profile_array(wall_distance, name="wall_distance")
    _validate_wall_grid(height)
    velocity_edge = _positive_scalar(edge_velocity, name="edge_velocity")
    density_edge = _positive_scalar(edge_density, name="edge_density")
    temperature_edge = _positive_scalar(edge_temperature, name="edge_temperature")
    thickness = _positive_scalar(
        boundary_layer_thickness, name="boundary_layer_thickness"
    )
    shear = _positive_scalar(wall_shear_stress, name="wall_shear_stress")
    prandtl = _positive_scalar(prandtl_number, name="prandtl_number")
    analogy_factor = _positive_scalar(
        reynolds_analogy_factor, name="reynolds_analogy_factor"
    )
    kappa = _positive_scalar(von_karman_constant, name="von_karman_constant")
    intercept_array, intercept_scalar = as_float_array(
        log_law_intercept, name="log_law_intercept"
    )
    if not intercept_scalar:
        raise ValueError("log_law_intercept must be a scalar")
    intercept = float(intercept_array)
    if height[-1] > thickness:
        raise ValueError("wall_distance must not exceed boundary_layer_thickness")

    recovery_factor = np.cbrt(prandtl)
    recovery_temperature = temperature_edge + recovery_factor * velocity_edge**2 / (
        2.0 * gas.cp
    )
    wall = (
        recovery_temperature
        if wall_temperature is None
        else _positive_scalar(wall_temperature, name="wall_temperature")
    )
    wall_density = density_edge * temperature_edge / wall
    wall_viscosity = float(viscosity_model.dynamic_viscosity(wall))
    friction_velocity = np.sqrt(shear / wall_density)
    friction_reynolds_number = (
        wall_density * friction_velocity * thickness / wall_viscosity
    )
    if friction_reynolds_number < 500.0:
        warnings.warn(
            "Spalding--Coles scale separation is weak for Re_tau < 500",
            ApplicabilityWarning,
            stacklevel=2,
        )

    parameters = _InverseParameters(
        transformation=transformation,
        temperature_velocity_relation=temperature_velocity_relation,
        edge_velocity=velocity_edge,
        edge_density=density_edge,
        edge_temperature=temperature_edge,
        wall_temperature=wall,
        recovery_temperature=float(recovery_temperature),
        wall_density=float(wall_density),
        wall_viscosity=wall_viscosity,
        friction_velocity=float(friction_velocity),
        friction_reynolds_number=float(friction_reynolds_number),
        prandtl_number=prandtl,
        reynolds_analogy_factor=analogy_factor,
        von_karman_constant=kappa,
        log_law_intercept=intercept,
        gas=gas,
        viscosity_model=viscosity_model,
    )
    wake, edge_ratio = _choose_wake_parameter(wake_parameter, parameters)

    requested_plus = wall_density * friction_velocity * height / wall_viscosity
    integration_plus = _integration_grid(requested_plus, friction_reynolds_number)
    velocity_plus_all, wall_velocity_plus_all = _integrate_inverse_profile(
        wake, integration_plus, parameters
    )
    velocity_all = velocity_plus_all * friction_velocity
    if np.any(velocity_all > velocity_edge * (1.0 + 1e-10)):
        raise ModelRangeError(
            "inverse model produced velocity greater than edge_velocity"
        )
    ratio_all = velocity_all / velocity_edge
    temperature_all = _temperature_profile(
        ratio_all,
        relation=temperature_velocity_relation,
        edge_temperature=temperature_edge,
        wall_temperature=wall,
        recovery_temperature=float(recovery_temperature),
        prandtl_number=prandtl,
        reynolds_analogy_factor=analogy_factor,
    )
    density_all = density_edge * temperature_edge / temperature_all
    physical_height_all = (
        integration_plus * wall_viscosity / (wall_density * friction_velocity)
    )
    mass_velocity_ratio = density_all * velocity_all / (density_edge * velocity_edge)
    displacement_thickness = float(
        np.trapezoid(1.0 - mass_velocity_ratio, physical_height_all)
    )
    momentum_thickness = float(
        np.trapezoid(
            mass_velocity_ratio * (1.0 - ratio_all),
            physical_height_all,
        )
    )
    if momentum_thickness <= 0.0:
        raise ModelRangeError(
            "inverse model produced a non-positive momentum thickness"
        )

    velocity_plus = np.interp(requested_plus, integration_plus, velocity_plus_all)
    wall_velocity_plus = np.interp(
        requested_plus, integration_plus, wall_velocity_plus_all
    )
    velocity = velocity_plus * friction_velocity
    velocity_ratio = velocity / velocity_edge
    temperature = _temperature_profile(
        velocity_ratio,
        relation=temperature_velocity_relation,
        edge_temperature=temperature_edge,
        wall_temperature=wall,
        recovery_temperature=float(recovery_temperature),
        prandtl_number=prandtl,
        reynolds_analogy_factor=analogy_factor,
    )
    density = density_edge * temperature_edge / temperature
    viscosity = np.asarray(
        viscosity_model.dynamic_viscosity(temperature), dtype=np.float64
    )
    outer_coordinate = height / thickness
    transformed_coordinate = _spalding_wall_coordinate(
        np.asarray(wall_velocity_plus, dtype=np.float64),
        von_karman_constant=kappa,
        intercept=intercept,
    )
    transformed_velocity = wall_velocity_plus + wake / kappa * _wake_function(
        outer_coordinate
    )
    speed_of_sound = np.asarray(gas.speed_of_sound(temperature), dtype=np.float64)
    mach = velocity / speed_of_sound
    dynamic_pressure = 0.5 * density * velocity**2
    skin_friction = 2.0 * shear / (density_edge * velocity_edge**2)

    return CompressibleBoundaryLayerProfileResult(
        transformation=transformation,
        temperature_velocity_relation=temperature_velocity_relation,
        wall_distance=height.copy(),
        wall_distance_plus=np.asarray(requested_plus, dtype=np.float64),
        transformed_wall_coordinate=transformed_coordinate,
        velocity=np.asarray(velocity, dtype=np.float64),
        velocity_plus=np.asarray(velocity_plus, dtype=np.float64),
        transformed_velocity_plus=np.asarray(transformed_velocity, dtype=np.float64),
        temperature=temperature,
        density=np.asarray(density, dtype=np.float64),
        dynamic_viscosity=viscosity,
        local_mach_number=np.asarray(mach, dtype=np.float64),
        dynamic_pressure=np.asarray(dynamic_pressure, dtype=np.float64),
        friction_velocity=float(friction_velocity),
        friction_reynolds_number=float(friction_reynolds_number),
        recovery_temperature=float(recovery_temperature),
        wall_temperature=float(wall),
        wake_parameter=wake,
        edge_velocity_ratio=edge_ratio,
        displacement_thickness=displacement_thickness,
        momentum_thickness=momentum_thickness,
        shape_factor=displacement_thickness / momentum_thickness,
        local_skin_friction_coefficient=float(skin_friction),
    )


__all__ = [
    "CompressibleBoundaryLayerProfileResult",
    "CompressibleVelocityTransformation",
    "TemperatureVelocityRelation",
    "TransformedVelocityProfileResult",
    "compressible_turbulent_boundary_layer_profile",
    "transform_compressible_velocity_profile",
]
