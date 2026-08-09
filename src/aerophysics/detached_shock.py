"""Engineering correlations for detached bow shocks over blunt noses.

The Ambrosio--Wortman and Billig correlations implemented here describe
continuum, low-temperature perfect-gas data for spherical/hemispherical and
two-dimensional cylindrical noses. For a hemispherical or cylindrical nose
followed by a body parallel to the free stream, Billig's far-field angle is
the Mach angle.

Coordinates use the nose-curvature centre as the origin, with positive
``x`` pointing upstream.  The body vertex is therefore at ``x = Rn`` and the
shock vertex at ``x = Rn + Delta``.

References
----------
Ambrosio, A. and Wortman, A., *Stagnation-Point Shock-Detachment Distance
for Flow around Spheres and Cylinders*, ARS Journal, 32(2), 281, 1962.
DOI: 10.2514/8.5988.
Billig, F. S., *Shock-Wave Shapes Around Spherical- and Cylindrical-Nosed
Bodies*, Journal of Spacecraft and Rockets, 4(6), 822--823, 1967.
DOI: 10.2514/3.28969.
Seiff, A., *Recent Information on Hypersonic Flow Fields*, in *The High
Temperature Aspects of Hypersonic Flow*, NASA SP-24, 1964.
Inouye, M., *Blunt Body Solutions for Spheres and Ellipsoids in Equilibrium
Gas Mixtures*, NASA TN D-2780, 1965.
"""

from dataclasses import dataclass
from enum import StrEnum

import numpy as np
from numpy.typing import ArrayLike

from aerophysics._array import FloatArray, FloatResult, as_float_array, return_float
from aerophysics.gas import AIR, PerfectGas
from aerophysics.shocks import normal_shock


class DetachedShockGeometry(StrEnum):
    """Blunt-nose geometry used by a detached-shock correlation."""

    AXISYMMETRIC_SPHERE = "axisymmetric_sphere"
    CYLINDRICAL_NOSE_2D = "cylindrical_nose_2d"


class DetachedShockModel(StrEnum):
    """Correlation used to predict detached-shock position or shape."""

    AMBROSIO_WORTMAN = "ambrosio_wortman"
    SEIFF = "seiff"
    BILLIG = "billig"


@dataclass(frozen=True, slots=True)
class DetachedShockStandoffResult:
    """Detached-shock standoff distance and correlation metadata."""

    upstream_mach: FloatResult | None
    nose_radius: FloatResult
    normalized_standoff_distance: FloatResult
    standoff_distance: FloatResult
    model: DetachedShockModel
    geometry: DetachedShockGeometry
    density_ratio: FloatResult | None = None


@dataclass(frozen=True, slots=True)
class BilligShockShapeResult:
    """Billig hyperbolic shock shape in nose-centred coordinates."""

    upstream_mach: FloatResult
    nose_radius: FloatResult
    normalized_standoff_distance: FloatResult
    standoff_distance: FloatResult
    vertex_curvature_radius: FloatResult
    transverse_coordinates: FloatArray
    shock_x: FloatArray
    shock_y: FloatArray
    model: DetachedShockModel
    standoff_model: DetachedShockModel
    geometry: DetachedShockGeometry


@dataclass(frozen=True, slots=True)
class DetachedShockComparisonResult:
    """Ambrosio--Wortman and Seiff sphere standoff comparison."""

    upstream_mach: FloatResult
    nose_radius: FloatResult
    ambrosio_wortman: DetachedShockStandoffResult
    seiff: DetachedShockStandoffResult
    normalized_standoff_difference: FloatResult
    standoff_distance_difference: FloatResult
    relative_difference: FloatResult
    geometry: DetachedShockGeometry


def _require_geometry(geometry: DetachedShockGeometry) -> None:
    if not isinstance(geometry, DetachedShockGeometry):
        raise ValueError("geometry must be a DetachedShockGeometry")


def _broadcast_positive_cases(
    first: ArrayLike,
    nose_radius: ArrayLike,
    *,
    first_name: str,
    first_minimum: float,
) -> tuple[FloatArray, FloatArray, bool]:
    values, values_scalar = as_float_array(first, name=first_name)
    radius, radius_scalar = as_float_array(nose_radius, name="nose_radius")
    if np.any(values <= first_minimum):
        raise ValueError(f"{first_name} must be greater than {first_minimum:g}")
    if np.any(radius <= 0.0):
        raise ValueError("nose_radius must be greater than zero")
    try:
        values, radius = np.broadcast_arrays(values, radius)
    except ValueError as error:
        raise ValueError(
            f"{first_name} and nose_radius must be broadcastable"
        ) from error
    return (
        np.asarray(values, dtype=np.float64),
        np.asarray(radius, dtype=np.float64),
        values_scalar and radius_scalar,
    )


def _standoff_coefficients(
    geometry: DetachedShockGeometry,
) -> tuple[float, float]:
    if geometry is DetachedShockGeometry.AXISYMMETRIC_SPHERE:
        return 0.143, 3.24
    return 0.386, 4.67


def shock_standoff_distance(
    upstream_mach: ArrayLike,
    nose_radius: ArrayLike,
    *,
    geometry: DetachedShockGeometry,
) -> DetachedShockStandoffResult:
    """Return the Ambrosio--Wortman detached-shock standoff distance."""
    _require_geometry(geometry)
    mach, radius, scalar = _broadcast_positive_cases(
        upstream_mach,
        nose_radius,
        first_name="upstream_mach",
        first_minimum=1.0,
    )
    coefficient, exponent = _standoff_coefficients(geometry)
    normalized = coefficient * np.exp(exponent / mach**2)
    distance = normalized * radius
    return DetachedShockStandoffResult(
        upstream_mach=return_float(mach, scalar=scalar),
        nose_radius=return_float(radius, scalar=scalar),
        normalized_standoff_distance=return_float(normalized, scalar=scalar),
        standoff_distance=return_float(distance, scalar=scalar),
        model=DetachedShockModel.AMBROSIO_WORTMAN,
        geometry=geometry,
    )


def seiff_standoff_distance(
    density_ratio: ArrayLike,
    nose_radius: ArrayLike,
) -> DetachedShockStandoffResult:
    """Return Seiff sphere standoff distance from ``rho2/rho1`` directly."""
    ratio, radius, scalar = _broadcast_positive_cases(
        density_ratio,
        nose_radius,
        first_name="density_ratio",
        first_minimum=1.0,
    )
    normalized = 0.78 / ratio
    distance = normalized * radius
    return DetachedShockStandoffResult(
        upstream_mach=None,
        nose_radius=return_float(radius, scalar=scalar),
        normalized_standoff_distance=return_float(normalized, scalar=scalar),
        standoff_distance=return_float(distance, scalar=scalar),
        model=DetachedShockModel.SEIFF,
        geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
        density_ratio=return_float(ratio, scalar=scalar),
    )


def seiff_standoff_distance_from_mach(
    upstream_mach: ArrayLike,
    nose_radius: ArrayLike,
    gas: PerfectGas = AIR,
) -> DetachedShockStandoffResult:
    """Return Seiff sphere standoff using a perfect-gas normal-shock ratio."""
    mach, radius, scalar = _broadcast_positive_cases(
        upstream_mach,
        nose_radius,
        first_name="upstream_mach",
        first_minimum=1.0,
    )
    normal = normal_shock(mach, gas)
    ratio = np.asarray(normal.static_density_ratio, dtype=np.float64)
    normalized = 0.78 / ratio
    distance = normalized * radius
    return DetachedShockStandoffResult(
        upstream_mach=return_float(mach, scalar=scalar),
        nose_radius=return_float(radius, scalar=scalar),
        normalized_standoff_distance=return_float(normalized, scalar=scalar),
        standoff_distance=return_float(distance, scalar=scalar),
        model=DetachedShockModel.SEIFF,
        geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
        density_ratio=return_float(ratio, scalar=scalar),
    )


def billig_shock_shape(
    upstream_mach: ArrayLike,
    nose_radius: ArrayLike,
    transverse_coordinates: ArrayLike,
    *,
    geometry: DetachedShockGeometry,
) -> BilligShockShapeResult:
    """Return Billig's hyperbolic shock coordinates at supplied ``y`` values."""
    _require_geometry(geometry)
    mach, radius, scalar = _broadcast_positive_cases(
        upstream_mach,
        nose_radius,
        first_name="upstream_mach",
        first_minimum=1.0,
    )
    transverse, _ = as_float_array(
        transverse_coordinates, name="transverse_coordinates"
    )
    if transverse.ndim != 1:
        raise ValueError("transverse_coordinates must be one-dimensional")
    if transverse.size == 0:
        raise ValueError("transverse_coordinates must not be empty")

    standoff = shock_standoff_distance(mach, radius, geometry=geometry)
    normalized = np.asarray(
        standoff.normalized_standoff_distance, dtype=np.float64
    )
    distance = np.asarray(standoff.standoff_distance, dtype=np.float64)
    if geometry is DetachedShockGeometry.AXISYMMETRIC_SPHERE:
        curvature = 1.143 * radius * np.exp(0.54 / (mach - 1.0) ** 1.2)
    else:
        curvature = 1.386 * radius * np.exp(1.8 / (mach - 1.0) ** 0.75)

    case_dimensions = (1,) * mach.ndim
    y = transverse.reshape((*case_dimensions, transverse.size))
    expanded_shape = (*mach.shape, transverse.size)
    shock_y = np.broadcast_to(y, expanded_shape).astype(np.float64, copy=False)
    mach_expanded = mach[..., np.newaxis]
    radius_expanded = radius[..., np.newaxis]
    distance_expanded = distance[..., np.newaxis]
    curvature_expanded = curvature[..., np.newaxis]
    beta = np.arcsin(1.0 / mach_expanded)
    tangent = np.tan(beta)
    shock_x = (
        radius_expanded
        + distance_expanded
        - curvature_expanded
        / tangent**2
        * (
            np.sqrt(1.0 + shock_y**2 * tangent**2 / curvature_expanded**2)
            - 1.0
        )
    )

    return BilligShockShapeResult(
        upstream_mach=return_float(mach, scalar=scalar),
        nose_radius=return_float(radius, scalar=scalar),
        normalized_standoff_distance=return_float(normalized, scalar=scalar),
        standoff_distance=return_float(distance, scalar=scalar),
        vertex_curvature_radius=return_float(curvature, scalar=scalar),
        transverse_coordinates=np.asarray(shock_y, dtype=np.float64),
        shock_x=np.asarray(shock_x, dtype=np.float64),
        shock_y=np.asarray(shock_y, dtype=np.float64),
        model=DetachedShockModel.BILLIG,
        standoff_model=DetachedShockModel.AMBROSIO_WORTMAN,
        geometry=geometry,
    )


def compare_standoff_distances(
    upstream_mach: ArrayLike,
    nose_radius: ArrayLike,
    gas: PerfectGas = AIR,
) -> DetachedShockComparisonResult:
    """Compare sphere standoff distances from Ambrosio--Wortman and Seiff."""
    ambrosio_wortman = shock_standoff_distance(
        upstream_mach,
        nose_radius,
        geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
    )
    seiff = seiff_standoff_distance_from_mach(upstream_mach, nose_radius, gas)
    aw_normalized = np.asarray(
        ambrosio_wortman.normalized_standoff_distance, dtype=np.float64
    )
    seiff_normalized = np.asarray(
        seiff.normalized_standoff_distance, dtype=np.float64
    )
    aw_distance = np.asarray(ambrosio_wortman.standoff_distance, dtype=np.float64)
    seiff_distance = np.asarray(seiff.standoff_distance, dtype=np.float64)
    difference = seiff_normalized - aw_normalized
    distance_difference = seiff_distance - aw_distance
    relative = difference / aw_normalized
    assert ambrosio_wortman.upstream_mach is not None
    scalar = np.ndim(ambrosio_wortman.upstream_mach) == 0
    return DetachedShockComparisonResult(
        upstream_mach=ambrosio_wortman.upstream_mach,
        nose_radius=ambrosio_wortman.nose_radius,
        ambrosio_wortman=ambrosio_wortman,
        seiff=seiff,
        normalized_standoff_difference=return_float(difference, scalar=scalar),
        standoff_distance_difference=return_float(
            distance_difference, scalar=scalar
        ),
        relative_difference=return_float(relative, scalar=scalar),
        geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
    )


__all__ = [
    "BilligShockShapeResult",
    "DetachedShockComparisonResult",
    "DetachedShockGeometry",
    "DetachedShockModel",
    "DetachedShockStandoffResult",
    "billig_shock_shape",
    "compare_standoff_distances",
    "seiff_standoff_distance",
    "seiff_standoff_distance_from_mach",
    "shock_standoff_distance",
]
