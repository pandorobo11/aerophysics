"""Smooth zero-pressure-gradient flat-plate boundary layers.

The correlations in this module describe one side of a smooth flat plate with
constant edge conditions. Distance is measured from the sharp leading edge.
All inputs and outputs use SI units.

References
----------
Schlichting, H. and Gersten, K., *Boundary-Layer Theory*, 8th ed.,
Springer, 2000.
White, F. M., *Viscous Fluid Flow*, 3rd ed., McGraw-Hill, 2006.
Glass, C. E. and Hunt, L. R., *Aerothermal Tests of Quilted Dome Models on
a Flat Plate at a Mach Number of 6.5*, NASA TP-2804, 1988.
Gnoffo, P. A., Berry, S. A., and Van Norman, J. W., *Uncertainty Assessments
of 2D and Axisymmetric Hypersonic Shock Wave--Turbulent Boundary Layer
Interaction Simulations at Compression Corners*, 2011, Appendix A.
"""

import warnings
from dataclasses import dataclass
from enum import StrEnum

import numpy as np
from numpy.typing import ArrayLike

from aerophysics._array import FloatArray, FloatResult, as_float_array, return_float
from aerophysics.exceptions import ApplicabilityWarning, ModelRangeError
from aerophysics.gas import AIR, AIR_VISCOSITY, PerfectGas, SutherlandModel


class BoundaryLayerRegime(StrEnum):
    """Boundary-layer state prescribed by the caller."""

    LAMINAR = "laminar"
    TURBULENT = "turbulent"
    TRANSITIONAL = "transitional"


class TurbulentCorrelation(StrEnum):
    """Smooth-plate turbulent skin-friction correlation."""

    POWER_LAW = "power_law"
    SCHLICHTING = "schlichting"


class CompressibilityCorrection(StrEnum):
    """Property and skin-friction treatment for compressible flow."""

    NONE = "none"
    ECKERT = "eckert"
    VAN_DRIEST_II = "van_driest_ii"


@dataclass(frozen=True, slots=True)
class FlatPlateBoundaryLayerResult:
    """Local boundary-layer properties and accumulated one-sided drag."""

    distance: FloatResult
    reynolds_number: FloatResult
    effective_reynolds_number: FloatResult
    boundary_layer_thickness: FloatResult
    displacement_thickness: FloatResult
    momentum_thickness: FloatResult
    local_skin_friction_coefficient: FloatResult
    average_skin_friction_coefficient: FloatResult
    wall_shear_stress: FloatResult
    drag_per_unit_width: FloatResult
    recovery_temperature: FloatResult | None
    wall_temperature: FloatResult | None


def _laminar_correlations(
    distance: FloatArray, reynolds: FloatArray
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]:
    root_reynolds = np.sqrt(reynolds)
    thickness = 5.0 * distance / root_reynolds
    displacement = 1.7208 * distance / root_reynolds
    momentum = 0.664 * distance / root_reynolds
    local_friction = 0.664 / root_reynolds
    average_friction = 1.328 / root_reynolds
    return (
        thickness,
        displacement,
        momentum,
        local_friction,
        average_friction,
    )


def _turbulent_friction(
    reynolds: FloatArray,
    correlation: TurbulentCorrelation,
) -> tuple[FloatArray, FloatArray]:
    if correlation is TurbulentCorrelation.POWER_LAW:
        local = 0.0592 * reynolds**-0.2
        average = 0.074 * reynolds**-0.2
        return local, average

    if np.any(reynolds <= 1.0):
        raise ModelRangeError(
            "Schlichting correlation requires Reynolds number greater than one"
        )
    logarithm = np.log10(reynolds)
    average = 0.455 / logarithm**2.58
    local = average * (1.0 - 2.58 / np.log(reynolds))
    return local, average


def _turbulent_correlations(
    distance: FloatArray,
    reynolds: FloatArray,
    correlation: TurbulentCorrelation,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]:
    thickness = 0.37 * distance * reynolds**-0.2
    displacement = 0.125 * thickness
    momentum = (7.0 / 72.0) * thickness
    local_friction, average_friction = _turbulent_friction(reynolds, correlation)
    return (
        thickness,
        displacement,
        momentum,
        local_friction,
        average_friction,
    )


def _validate_model_choices(
    regime: BoundaryLayerRegime,
    turbulent_correlation: TurbulentCorrelation,
    transition_reynolds: float | None,
) -> float | None:
    if not isinstance(regime, BoundaryLayerRegime):
        raise ValueError("regime must be a BoundaryLayerRegime value")
    if not isinstance(turbulent_correlation, TurbulentCorrelation):
        raise ValueError("turbulent_correlation must be a TurbulentCorrelation value")
    if regime is BoundaryLayerRegime.TRANSITIONAL:
        if transition_reynolds is None:
            raise ValueError(
                "transition_reynolds is required for a transitional boundary layer"
            )
        if not np.isfinite(transition_reynolds) or transition_reynolds <= 0.0:
            raise ValueError("transition_reynolds must be finite and greater than zero")
        return float(transition_reynolds)
    if transition_reynolds is not None:
        raise ValueError(
            "transition_reynolds is only valid for a transitional boundary layer"
        )
    return None


def _warn_outside_correlation_range(
    reynolds: FloatArray,
    turbulent_mask: np.ndarray,
) -> None:
    if np.any(turbulent_mask & ((reynolds < 5e5) | (reynolds > 1e9))):
        warnings.warn(
            "turbulent flat-plate correlations are nominally valid for "
            "5e5 <= Reynolds number <= 1e9",
            ApplicabilityWarning,
            stacklevel=3,
        )


def _recovery_temperature(
    edge_temperature: FloatArray,
    mach: FloatArray,
    *,
    turbulent: bool,
    prandtl_number: float,
    gas: PerfectGas,
) -> FloatArray:
    recovery_factor = np.cbrt(prandtl_number) if turbulent else np.sqrt(prandtl_number)
    recovery = edge_temperature * (
        1.0 + recovery_factor * 0.5 * (gas.heat_capacity_ratio - 1.0) * mach**2
    )
    return np.asarray(recovery, dtype=np.float64)


def _eckert_effective_reynolds(
    reynolds: FloatArray,
    edge_temperature: FloatArray,
    mach: FloatArray,
    specified_wall_temperature: FloatArray | None,
    *,
    turbulent: bool,
    prandtl_number: float,
    gas: PerfectGas,
    viscosity_model: SutherlandModel,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    recovery = _recovery_temperature(
        edge_temperature,
        mach,
        turbulent=turbulent,
        prandtl_number=prandtl_number,
        gas=gas,
    )
    wall = (
        recovery if specified_wall_temperature is None else specified_wall_temperature
    )
    reference = 0.22 * recovery + 0.28 * edge_temperature + 0.50 * wall
    edge_viscosity = np.asarray(
        viscosity_model.dynamic_viscosity(edge_temperature), dtype=np.float64
    )
    reference_viscosity = np.asarray(
        viscosity_model.dynamic_viscosity(reference), dtype=np.float64
    )
    effective_reynolds = (
        reynolds
        * (edge_temperature / reference)
        * (edge_viscosity / reference_viscosity)
    )
    return (
        np.asarray(effective_reynolds, dtype=np.float64),
        recovery,
        np.asarray(wall, dtype=np.float64),
    )


def _van_driest_ii_state(
    reynolds: FloatArray,
    edge_temperature: FloatArray,
    mach: FloatArray,
    specified_wall_temperature: FloatArray | None,
    *,
    prandtl_number: float,
    gas: PerfectGas,
    viscosity_model: SutherlandModel,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    recovery = _recovery_temperature(
        edge_temperature,
        mach,
        turbulent=True,
        prandtl_number=prandtl_number,
        gas=gas,
    )
    wall = (
        recovery if specified_wall_temperature is None else specified_wall_temperature
    )
    edge_viscosity = np.asarray(
        viscosity_model.dynamic_viscosity(edge_temperature), dtype=np.float64
    )
    with np.errstate(over="ignore", invalid="ignore"):
        wall_viscosity = np.asarray(
            viscosity_model.dynamic_viscosity(wall), dtype=np.float64
        )
    theta_factor = edge_viscosity / wall_viscosity
    effective_reynolds = reynolds * theta_factor

    recovery_rise = recovery / edge_temperature - 1.0
    a = np.sqrt(
        0.5
        * prandtl_number ** (1.0 / 3.0)
        * (gas.heat_capacity_ratio - 1.0)
        * mach**2
        * edge_temperature
        / wall
    )
    b = recovery / wall - 1.0
    denominator = np.sqrt(b**2 + 4.0 * a**2)
    alpha = np.clip((2.0 * a**2 - b) / denominator, -1.0, 1.0)
    beta = np.clip(b / denominator, -1.0, 1.0)
    angle = np.arcsin(alpha) + np.arcsin(beta)
    with np.errstate(divide="ignore", invalid="ignore"):
        friction_factor = recovery_rise / angle**2
    if np.any(~np.isfinite(friction_factor) | (friction_factor <= 0.0)):
        raise ModelRangeError(
            "Van Driest II transformation is undefined for this thermal state"
        )
    return (
        np.asarray(effective_reynolds, dtype=np.float64),
        np.asarray(friction_factor, dtype=np.float64),
        recovery,
        np.asarray(wall, dtype=np.float64),
    )


def flat_plate_boundary_layer(
    distance: ArrayLike,
    edge_velocity: ArrayLike,
    edge_density: ArrayLike,
    edge_dynamic_viscosity: ArrayLike,
    *,
    regime: BoundaryLayerRegime,
    turbulent_correlation: TurbulentCorrelation = (TurbulentCorrelation.SCHLICHTING),
    transition_reynolds: float | None = None,
    compressibility_correction: CompressibilityCorrection = (
        CompressibilityCorrection.NONE
    ),
    mach: ArrayLike | None = None,
    edge_temperature: ArrayLike | None = None,
    wall_temperature: ArrayLike | None = None,
    prandtl_number: float = 0.72,
    gas: PerfectGas = AIR,
    viscosity_model: SutherlandModel = AIR_VISCOSITY,
) -> FlatPlateBoundaryLayerResult:
    """Return a smooth flat-plate boundary-layer estimate.

    Parameters
    ----------
    distance:
        Distance from the leading edge in m. Values must be positive.
    edge_velocity:
        Boundary-layer edge velocity in m/s.
    edge_density:
        Boundary-layer edge density in kg/m³.
    edge_dynamic_viscosity:
        Boundary-layer edge dynamic viscosity in Pa s.
    regime:
        Explicit laminar, turbulent, or transitional model selection.
    turbulent_correlation:
        Skin-friction correlation used wherever the boundary layer is
        turbulent. The default is the Schlichting logarithmic correlation.
    transition_reynolds:
        Reynolds number at transition. Required only for ``TRANSITIONAL``;
        there is deliberately no default transition location.
    compressibility_correction:
        ``NONE``, Eckert's reference-temperature method, or Van Driest II.
        With Van Driest II, laminar portions use Eckert and turbulent portions
        use Van Driest II.
    mach, edge_temperature:
        Required for a compressibility correction. Temperature is in K.
    wall_temperature:
        Wall temperature in K. If omitted for a compressibility correction,
        an adiabatic wall at the regime-appropriate recovery temperature is
        used.
    prandtl_number:
        Constant Prandtl number used for the recovery factors.
    gas, viscosity_model:
        Perfect-gas and Sutherland models used by compressibility corrections.

    Notes
    -----
    ``boundary_layer_thickness`` is the conventional 99-percent thickness.
    Turbulent thicknesses use the one-fifth-power smooth-plate thickness
    correlation and a one-seventh-power velocity-profile approximation.
    ``drag_per_unit_width`` is accumulated from the leading edge to
    ``distance`` on one side of the plate.
    """
    transition = _validate_model_choices(
        regime, turbulent_correlation, transition_reynolds
    )
    if not isinstance(compressibility_correction, CompressibilityCorrection):
        raise ValueError(
            "compressibility_correction must be a CompressibilityCorrection value"
        )
    if not np.isfinite(prandtl_number) or prandtl_number <= 0.0:
        raise ValueError("prandtl_number must be finite and greater than zero")
    raw_values = (
        ("distance", distance),
        ("edge_velocity", edge_velocity),
        ("edge_density", edge_density),
        ("edge_dynamic_viscosity", edge_dynamic_viscosity),
    )
    arrays: list[FloatArray] = []
    scalar = True
    for name, value in raw_values:
        array, is_scalar = as_float_array(value, name=name)
        arrays.append(array)
        scalar = scalar and is_scalar
    thermal_values: list[tuple[str, ArrayLike]] = []
    if compressibility_correction is CompressibilityCorrection.NONE:
        if (
            mach is not None
            or edge_temperature is not None
            or wall_temperature is not None
        ):
            raise ValueError(
                "thermal inputs require an explicit compressibility correction"
            )
    else:
        if mach is None or edge_temperature is None:
            raise ValueError(
                "mach and edge_temperature are required for compressibility correction"
            )
        thermal_values.extend((("mach", mach), ("edge_temperature", edge_temperature)))
        if wall_temperature is not None:
            thermal_values.append(("wall_temperature", wall_temperature))
    for name, value in thermal_values:
        array, is_scalar = as_float_array(value, name=name)
        arrays.append(array)
        scalar = scalar and is_scalar
    try:
        broadcast = np.broadcast_arrays(*arrays)
    except ValueError as error:
        raise ValueError("flat-plate inputs must be broadcastable") from error
    x, velocity, density, viscosity = broadcast[:4]
    if np.any(x <= 0.0):
        raise ValueError("distance must be greater than zero")
    if np.any(velocity <= 0.0):
        raise ValueError("edge_velocity must be greater than zero")
    if np.any(density <= 0.0):
        raise ValueError("edge_density must be greater than zero")
    if np.any(viscosity <= 0.0):
        raise ValueError("edge_dynamic_viscosity must be greater than zero")

    reynolds = density * velocity * x / viscosity
    if regime is BoundaryLayerRegime.LAMINAR:
        turbulent_mask = np.zeros_like(reynolds, dtype=np.bool_)
    elif regime is BoundaryLayerRegime.TURBULENT:
        turbulent_mask = np.ones_like(reynolds, dtype=np.bool_)
    else:
        assert transition is not None
        turbulent_mask = reynolds > transition

    laminar_reynolds = np.asarray(reynolds, dtype=np.float64)
    turbulent_reynolds = np.asarray(reynolds, dtype=np.float64)
    turbulent_friction_factor = np.ones_like(reynolds)
    recovery_temperature: FloatArray | None = None
    used_wall_temperature: FloatArray | None = None
    if compressibility_correction is not CompressibilityCorrection.NONE:
        mach_values = np.asarray(broadcast[4], dtype=np.float64)
        edge_temperature_values = np.asarray(broadcast[5], dtype=np.float64)
        specified_wall = (
            np.asarray(broadcast[6], dtype=np.float64)
            if wall_temperature is not None
            else None
        )
        if np.any(mach_values < 0.0):
            raise ValueError("mach must be non-negative")
        if np.any(edge_temperature_values <= 0.0):
            raise ValueError("edge_temperature must be greater than zero")
        if specified_wall is not None and np.any(specified_wall <= 0.0):
            raise ValueError("wall_temperature must be greater than zero")
        if (
            compressibility_correction is CompressibilityCorrection.VAN_DRIEST_II
            and np.any(turbulent_mask & (mach_values <= 0.0))
        ):
            raise ModelRangeError(
                "Van Driest II requires Mach number greater than zero"
            )

        (
            laminar_reynolds,
            laminar_recovery,
            laminar_wall,
        ) = _eckert_effective_reynolds(
            reynolds,
            edge_temperature_values,
            mach_values,
            specified_wall,
            turbulent=False,
            prandtl_number=prandtl_number,
            gas=gas,
            viscosity_model=viscosity_model,
        )
        if regime is BoundaryLayerRegime.LAMINAR:
            recovery_temperature = laminar_recovery
            used_wall_temperature = laminar_wall
        elif compressibility_correction is CompressibilityCorrection.ECKERT:
            (
                turbulent_reynolds,
                turbulent_recovery,
                turbulent_wall,
            ) = _eckert_effective_reynolds(
                reynolds,
                edge_temperature_values,
                mach_values,
                specified_wall,
                turbulent=True,
                prandtl_number=prandtl_number,
                gas=gas,
                viscosity_model=viscosity_model,
            )
            recovery_temperature = np.where(
                turbulent_mask, turbulent_recovery, laminar_recovery
            )
            used_wall_temperature = np.where(
                turbulent_mask, turbulent_wall, laminar_wall
            )
        else:
            (
                turbulent_reynolds,
                turbulent_friction_factor,
                turbulent_recovery,
                turbulent_wall,
            ) = _van_driest_ii_state(
                reynolds,
                edge_temperature_values,
                mach_values,
                specified_wall,
                prandtl_number=prandtl_number,
                gas=gas,
                viscosity_model=viscosity_model,
            )
            recovery_temperature = np.where(
                turbulent_mask, turbulent_recovery, laminar_recovery
            )
            used_wall_temperature = np.where(
                turbulent_mask, turbulent_wall, laminar_wall
            )

    laminar = _laminar_correlations(x, laminar_reynolds)
    if regime is BoundaryLayerRegime.LAMINAR:
        selected = laminar
    elif regime is BoundaryLayerRegime.TURBULENT:
        turbulent_base = _turbulent_correlations(
            x, turbulent_reynolds, turbulent_correlation
        )
        turbulent = (
            turbulent_base[0],
            turbulent_base[1],
            turbulent_base[2],
            turbulent_base[3] / turbulent_friction_factor,
            turbulent_base[4] / turbulent_friction_factor,
        )
        selected = turbulent
    else:
        assert transition is not None
        turbulent_scale = turbulent_reynolds / reynolds
        turbulent_base = _turbulent_correlations(
            x,
            turbulent_scale * np.maximum(reynolds, transition),
            turbulent_correlation,
        )
        turbulent = (
            turbulent_base[0],
            turbulent_base[1],
            turbulent_base[2],
            turbulent_base[3] / turbulent_friction_factor,
            turbulent_base[4] / turbulent_friction_factor,
        )
        laminar_transition_reynolds = laminar_reynolds / reynolds * transition
        turbulent_transition_reynolds = turbulent_scale * transition
        laminar_average = 1.328 / np.sqrt(laminar_transition_reynolds)
        _, turbulent_average = _turbulent_friction(
            turbulent_transition_reynolds, turbulent_correlation
        )
        turbulent_average = turbulent_average / turbulent_friction_factor
        mixed_average = turbulent[4] + (transition / reynolds) * (
            laminar_average - turbulent_average
        )
        selected = (
            np.where(turbulent_mask, turbulent[0], laminar[0]),
            np.where(turbulent_mask, turbulent[1], laminar[1]),
            np.where(turbulent_mask, turbulent[2], laminar[2]),
            np.where(turbulent_mask, turbulent[3], laminar[3]),
            np.where(turbulent_mask, mixed_average, laminar[4]),
        )

    _warn_outside_correlation_range(reynolds, turbulent_mask)
    thickness, displacement, momentum, local_friction, average_friction = selected
    effective_reynolds = np.where(turbulent_mask, turbulent_reynolds, laminar_reynolds)
    dynamic_pressure = 0.5 * density * velocity**2
    shear = dynamic_pressure * local_friction
    drag = dynamic_pressure * x * average_friction

    def output(values: FloatArray) -> FloatResult:
        return return_float(np.asarray(values, dtype=np.float64), scalar=scalar)

    return FlatPlateBoundaryLayerResult(
        distance=output(x),
        reynolds_number=output(reynolds),
        effective_reynolds_number=output(effective_reynolds),
        boundary_layer_thickness=output(thickness),
        displacement_thickness=output(displacement),
        momentum_thickness=output(momentum),
        local_skin_friction_coefficient=output(local_friction),
        average_skin_friction_coefficient=output(average_friction),
        wall_shear_stress=output(shear),
        drag_per_unit_width=output(drag),
        recovery_temperature=(
            output(recovery_temperature) if recovery_temperature is not None else None
        ),
        wall_temperature=(
            output(used_wall_temperature) if used_wall_temperature is not None else None
        ),
    )


__all__ = [
    "BoundaryLayerRegime",
    "CompressibilityCorrection",
    "FlatPlateBoundaryLayerResult",
    "TurbulentCorrelation",
    "flat_plate_boundary_layer",
]
