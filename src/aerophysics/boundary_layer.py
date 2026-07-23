"""Smooth zero-pressure-gradient flat-plate boundary layers.

The correlations in this module describe one side of a smooth flat plate with
constant edge conditions. Distance is measured from the sharp leading edge.
All inputs and outputs use SI units.

References
----------
Schlichting, H. and Gersten, K., *Boundary-Layer Theory*, 8th ed.,
Springer, 2000.
White, F. M., *Viscous Fluid Flow*, 3rd ed., McGraw-Hill, 2006.
"""

import warnings
from dataclasses import dataclass
from enum import StrEnum

import numpy as np
from numpy.typing import ArrayLike

from aerophysics._array import FloatArray, FloatResult, as_float_array, return_float
from aerophysics.exceptions import ApplicabilityWarning, ModelRangeError


class BoundaryLayerRegime(StrEnum):
    """Boundary-layer state prescribed by the caller."""

    LAMINAR = "laminar"
    TURBULENT = "turbulent"
    TRANSITIONAL = "transitional"


class TurbulentCorrelation(StrEnum):
    """Smooth-plate turbulent skin-friction correlation."""

    POWER_LAW = "power_law"
    SCHLICHTING = "schlichting"


@dataclass(frozen=True, slots=True)
class FlatPlateBoundaryLayerResult:
    """Local boundary-layer properties and accumulated one-sided drag."""

    distance: FloatResult
    reynolds_number: FloatResult
    boundary_layer_thickness: FloatResult
    displacement_thickness: FloatResult
    momentum_thickness: FloatResult
    local_skin_friction_coefficient: FloatResult
    average_skin_friction_coefficient: FloatResult
    wall_shear_stress: FloatResult
    drag_per_unit_width: FloatResult


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


def flat_plate_boundary_layer(
    distance: ArrayLike,
    edge_velocity: ArrayLike,
    edge_density: ArrayLike,
    edge_dynamic_viscosity: ArrayLike,
    *,
    regime: BoundaryLayerRegime,
    turbulent_correlation: TurbulentCorrelation = (TurbulentCorrelation.SCHLICHTING),
    transition_reynolds: float | None = None,
) -> FlatPlateBoundaryLayerResult:
    """Return an incompressible smooth flat-plate boundary-layer estimate.

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
    try:
        x, velocity, density, viscosity = np.broadcast_arrays(*arrays)
    except ValueError as error:
        raise ValueError("flat-plate inputs must be broadcastable") from error
    if np.any(x <= 0.0):
        raise ValueError("distance must be greater than zero")
    if np.any(velocity <= 0.0):
        raise ValueError("edge_velocity must be greater than zero")
    if np.any(density <= 0.0):
        raise ValueError("edge_density must be greater than zero")
    if np.any(viscosity <= 0.0):
        raise ValueError("edge_dynamic_viscosity must be greater than zero")

    reynolds = density * velocity * x / viscosity
    laminar = _laminar_correlations(x, reynolds)

    if regime is BoundaryLayerRegime.LAMINAR:
        selected = laminar
        turbulent_mask = np.zeros_like(reynolds, dtype=np.bool_)
    elif regime is BoundaryLayerRegime.TURBULENT:
        turbulent = _turbulent_correlations(x, reynolds, turbulent_correlation)
        selected = turbulent
        turbulent_mask = np.ones_like(reynolds, dtype=np.bool_)
    else:
        assert transition is not None
        turbulent_mask = reynolds > transition
        turbulent = _turbulent_correlations(
            x,
            np.maximum(reynolds, transition),
            turbulent_correlation,
        )
        transition_array = np.asarray(transition, dtype=np.float64)
        laminar_average = np.asarray(
            1.328 / np.sqrt(transition_array), dtype=np.float64
        )
        _, turbulent_average = _turbulent_friction(
            transition_array, turbulent_correlation
        )
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
    dynamic_pressure = 0.5 * density * velocity**2
    shear = dynamic_pressure * local_friction
    drag = dynamic_pressure * x * average_friction

    def output(values: FloatArray) -> FloatResult:
        return return_float(np.asarray(values, dtype=np.float64), scalar=scalar)

    return FlatPlateBoundaryLayerResult(
        distance=output(x),
        reynolds_number=output(reynolds),
        boundary_layer_thickness=output(thickness),
        displacement_thickness=output(displacement),
        momentum_thickness=output(momentum),
        local_skin_friction_coefficient=output(local_friction),
        average_skin_friction_coefficient=output(average_friction),
        wall_shear_stress=output(shear),
        drag_per_unit_width=output(drag),
    )


__all__ = [
    "BoundaryLayerRegime",
    "FlatPlateBoundaryLayerResult",
    "TurbulentCorrelation",
    "flat_plate_boundary_layer",
]
