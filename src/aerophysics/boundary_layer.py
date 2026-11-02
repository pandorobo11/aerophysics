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
Van Driest, E. R., *Turbulent Boundary Layer in Compressible Fluids*, 1951,
doi:10.2514/8.1895.
Hopkins, E. J. and Inouye, M., *An Evaluation of Theories for Predicting
Turbulent Skin Friction and Heat Transfer on Flat Plates at Supersonic and
Hypersonic Mach Numbers*, AIAA Journal, 1971.
Hopkins, E. J., *Charts for Predicting Turbulent Skin Friction from the Van
Driest Method (II)*, NASA TN D-6945, 1972.
Willems, S. and Gülhan, A., *Experiments on Shock Induced Laminar-Turbulent
Transition on a Flat Plate at Mach 6*, EUCASS 2013, equation (7).
"""

import warnings
from dataclasses import dataclass
from enum import StrEnum

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import brentq

from aerophysics._array import FloatArray, FloatResult, as_float_array, return_float
from aerophysics._laminar_similarity import similarity_coefficients
from aerophysics.boundary_layer_profile import (
    _van_driest_walz_profile,
)
from aerophysics.exceptions import ApplicabilityWarning, ModelRangeError
from aerophysics.gas import AIR, PerfectGas
from aerophysics.transport import AIR_VISCOSITY, DynamicViscosityModel

_VD2_WALL_TEMPERATURE_EXPONENT = 0.064
_DELTA99_MOMENTUM_FRACTION = 0.978
# Blasius velocity reaches 99 percent at this similarity coordinate.
_BLASIUS_DELTA99 = 4.9100
_VISCOUS_INTERACTION_THRESHOLD = 1.0


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
    VAN_DRIEST_II_WALL_TEMPERATURE = "van_driest_ii_wall_temperature"


@dataclass(frozen=True, slots=True)
class FlatPlateBoundaryLayerResult:
    """Local boundary-layer properties and accumulated one-sided drag.

    ``local_skin_friction_coefficient`` and
    ``average_skin_friction_coefficient`` are always referred to the edge
    dynamic pressure ``rho_e U_e**2/2`` for every compressibility correction,
    so ``wall_shear_stress`` is the local coefficient multiplied by the edge
    dynamic pressure, while ``drag_per_unit_width`` is the average coefficient
    multiplied by the edge dynamic pressure and ``distance``. With Eckert
    turbulent portions, ``effective_reynolds_number`` is the reference Reynolds number
    ``Re*`` at which the incompressible friction correlations are evaluated.
    Laminar similarity portions return the edge Reynolds number ``Re_x``.
    """

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
    thickness = _BLASIUS_DELTA99 * distance / root_reynolds
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


def _warn_viscous_interaction(
    reynolds: FloatArray,
    mach: FloatArray,
    edge_temperature: FloatArray,
    edge_viscosity: FloatArray,
    laminar_wall: FloatArray,
    turbulent_mask: np.ndarray,
    viscosity_model: DynamicViscosityModel,
) -> None:
    laminar_mask = ~turbulent_mask
    if not np.any(laminar_mask):
        return
    wall = laminar_wall[laminar_mask]
    wall_viscosity = np.asarray(viscosity_model.dynamic_viscosity(wall))
    wall_chapman_rubesin = (
        edge_temperature[laminar_mask]
        / wall
        * wall_viscosity
        / edge_viscosity[laminar_mask]
    )
    interaction = mach[laminar_mask] ** 3 * np.sqrt(
        wall_chapman_rubesin / reynolds[laminar_mask]
    )
    affected = interaction > _VISCOUS_INTERACTION_THRESHOLD
    if np.any(affected):
        warnings.warn(
            "laminar viscous interaction: leading-edge induced pressure invalidates "
            "the zero-pressure-gradient flat-plate estimate "
            f"(chi_bar > 1 at {np.count_nonzero(affected)} point(s); "
            f"maximum chi_bar = {np.max(interaction):.6g})",
            ApplicabilityWarning,
            stacklevel=3,
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
    thickness, displacement, momentum = _turbulent_thicknesses(distance, reynolds)
    local_friction, average_friction = _turbulent_friction(reynolds, correlation)
    return (
        thickness,
        displacement,
        momentum,
        local_friction,
        average_friction,
    )


def _turbulent_thicknesses(
    distance: FloatArray,
    reynolds: FloatArray,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    thickness = 0.37 * distance * reynolds**-0.2
    displacement = 0.125 * thickness
    momentum = (7.0 / 72.0) * thickness
    return thickness, displacement, momentum


def _edge_basis_friction(
    correlations: tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray],
    density_ratio: FloatArray,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]:
    """Rescale friction coefficients by ``rho*/rho_e``; thicknesses are lengths."""
    thickness, displacement, momentum, local_friction, average_friction = correlations
    return (
        thickness,
        displacement,
        momentum,
        density_ratio * local_friction,
        density_ratio * average_friction,
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
    prandtl_number: float,
    gas: PerfectGas,
) -> FloatArray:
    recovery_factor = np.cbrt(prandtl_number)
    recovery = edge_temperature * (
        1.0 + recovery_factor * 0.5 * (gas.heat_capacity_ratio - 1.0) * mach**2
    )
    return np.asarray(recovery, dtype=np.float64)


def _eckert_effective_reynolds(
    reynolds: FloatArray,
    edge_temperature: FloatArray,
    edge_viscosity: FloatArray,
    mach: FloatArray,
    specified_wall_temperature: FloatArray | None,
    *,
    prandtl_number: float,
    gas: PerfectGas,
    viscosity_model: DynamicViscosityModel,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    """Return ``Re*``, recovery and wall temperatures, and ``rho*/rho_e``.

    The edge pressure is constant through the layer, so the perfect-gas
    reference-to-edge density ratio is ``T_e/T*``.  Correlations evaluated at
    ``Re*`` give coefficients on the reference dynamic pressure
    ``rho* U_e**2/2``; multiplying by this ratio converts them to the edge
    dynamic-pressure basis.
    """
    recovery = _recovery_temperature(
        edge_temperature,
        mach,
        prandtl_number=prandtl_number,
        gas=gas,
    )
    wall = (
        recovery if specified_wall_temperature is None else specified_wall_temperature
    )
    reference = 0.22 * recovery + 0.28 * edge_temperature + 0.50 * wall
    density_ratio = edge_temperature / reference
    reference_viscosity = np.asarray(
        viscosity_model.dynamic_viscosity(reference), dtype=np.float64
    )
    effective_reynolds = (
        reynolds * density_ratio * (edge_viscosity / reference_viscosity)
    )
    return (
        np.asarray(effective_reynolds, dtype=np.float64),
        recovery,
        np.asarray(wall, dtype=np.float64),
        np.asarray(density_ratio, dtype=np.float64),
    )


@dataclass(frozen=True, slots=True)
class _VanDriestIIState:
    friction_reynolds: FloatArray
    thickness_reynolds: FloatArray
    friction_factor: FloatArray
    momentum_factor: FloatArray
    reynolds_factor: FloatArray
    recovery_temperature: FloatArray
    wall_temperature: FloatArray


def _van_driest_ii_state(
    reynolds: FloatArray,
    edge_temperature: FloatArray,
    edge_viscosity: FloatArray,
    mach: FloatArray,
    specified_wall_temperature: FloatArray | None,
    *,
    prandtl_number: float,
    gas: PerfectGas,
    viscosity_model: DynamicViscosityModel,
) -> _VanDriestIIState:
    """Return the Hopkins--Inouye form of the Van Driest II factors."""
    recovery = _recovery_temperature(
        edge_temperature,
        mach,
        prandtl_number=prandtl_number,
        gas=gas,
    )
    wall = (
        recovery if specified_wall_temperature is None else specified_wall_temperature
    )
    with np.errstate(over="ignore", invalid="ignore"):
        wall_viscosity = np.asarray(
            viscosity_model.dynamic_viscosity(wall), dtype=np.float64
        )
    momentum_factor = edge_viscosity / wall_viscosity

    # Hopkins--Inouye / Willems notation.  This is the Van Driest II
    # skin-friction transformation, not the 1951 velocity transformation.
    recovery_rise = (
        np.cbrt(prandtl_number) * 0.5 * (gas.heat_capacity_ratio - 1.0) * mach**2
    )
    temperature_factor = wall / edge_temperature
    discriminant = (
        recovery_rise + 1.0 + temperature_factor
    ) ** 2 - 4.0 * temperature_factor
    incompressible_limit = recovery_rise <= 1e-12
    safe_discriminant = np.where(incompressible_limit, 1.0, discriminant)
    denominator = np.sqrt(safe_discriminant)
    alpha = np.clip(
        (recovery_rise - 1.0 + temperature_factor) / denominator,
        -1.0,
        1.0,
    )
    beta = np.clip(
        (recovery_rise + 1.0 - temperature_factor) / denominator,
        -1.0,
        1.0,
    )
    angle = np.arcsin(alpha) + np.arcsin(beta)
    with np.errstate(divide="ignore", invalid="ignore"):
        friction_factor = recovery_rise / angle**2
    friction_factor = np.where(incompressible_limit, 1.0, friction_factor)
    momentum_factor = np.where(incompressible_limit, 1.0, momentum_factor)
    reynolds_factor = momentum_factor / friction_factor
    friction_reynolds = reynolds * reynolds_factor
    thickness_reynolds = reynolds * momentum_factor
    invalid = (
        ~np.isfinite(friction_factor)
        | (friction_factor <= 0.0)
        | ~np.isfinite(momentum_factor)
        | (momentum_factor <= 0.0)
        | ~np.isfinite(reynolds_factor)
        | (reynolds_factor <= 0.0)
        | ~np.isfinite(friction_reynolds)
        | (friction_reynolds <= 0.0)
    )
    if np.any(invalid):
        raise ModelRangeError(
            "Van Driest II transformation is undefined for this thermal state"
        )
    return _VanDriestIIState(
        friction_reynolds=np.asarray(friction_reynolds, dtype=np.float64),
        thickness_reynolds=np.asarray(thickness_reynolds, dtype=np.float64),
        friction_factor=np.asarray(friction_factor, dtype=np.float64),
        momentum_factor=np.asarray(momentum_factor, dtype=np.float64),
        reynolds_factor=np.asarray(reynolds_factor, dtype=np.float64),
        recovery_temperature=recovery,
        wall_temperature=np.asarray(wall, dtype=np.float64),
    )


def _positive_decreasing_root(
    reynolds: FloatArray,
    *,
    right_hand_offset: float,
    description: str,
) -> FloatArray:
    """Solve independent positive monotonic equations with Brent's method."""
    lower = np.finfo(np.float64).eps
    upper = 1.0
    roots = np.empty_like(reynolds)
    for index in np.ndindex(reynolds.shape):
        reynolds_value = float(reynolds[index])

        def residual(
            friction: float,
            reynolds_value: float = reynolds_value,
        ) -> float:
            return float(
                0.242 / np.sqrt(friction)
                - right_hand_offset
                - np.log10(reynolds_value * friction)
            )

        lower_residual = residual(lower)
        upper_residual = residual(upper)
        if (
            not np.isfinite(lower_residual)
            or lower_residual <= 0.0
            or not np.isfinite(upper_residual)
            or upper_residual >= 0.0
        ):
            raise ModelRangeError(
                f"could not bracket a positive {description} solution"
            )
        try:
            roots[index] = brentq(
                residual,
                lower,
                upper,
                xtol=np.finfo(np.float64).tiny,
                rtol=4.0 * np.finfo(np.float64).eps,
            )
        except (ValueError, RuntimeError) as error:
            raise ModelRangeError(
                f"could not solve the {description} equation"
            ) from error
    return np.asarray(roots, dtype=np.float64)


def _willems_friction(
    reynolds: FloatArray,
) -> tuple[FloatArray, FloatArray]:
    """Solve the Willems et al. (2013) local and mean implicit equations."""
    if np.any(~np.isfinite(reynolds) | (reynolds <= 0.0)):
        raise ModelRangeError(
            "Willems skin-friction equations require a positive finite Reynolds number"
        )

    local = _positive_decreasing_root(
        reynolds,
        right_hand_offset=0.41,
        description="local Willems skin-friction",
    )
    average = _positive_decreasing_root(
        reynolds,
        right_hand_offset=0.0,
        description="average Willems skin-friction",
    )
    return local, average


def _wall_temperature_friction(
    reynolds: FloatArray,
) -> tuple[FloatArray, FloatArray]:
    return _turbulent_friction(reynolds, TurbulentCorrelation.SCHLICHTING)


def _wall_temperature_thickness(
    theta: float,
    velocity: float,
    density: float,
    temperature: float,
    wall: float,
    shear: float,
    *,
    prandtl_number: float,
    gas: PerfectGas,
    viscosity_model: DynamicViscosityModel,
) -> tuple[float, float]:
    profile = _van_driest_walz_profile(
        velocity,
        density,
        temperature,
        wall,
        shear,
        prandtl_number=prandtl_number,
        gas=gas,
        viscosity_model=viscosity_model,
    )
    # Pi decreases with thickness. Its two admissible endpoints give the
    # entire valid bracket, including roots arbitrarily close to Pi=0 or 1.
    lower = profile.thickness_at_wake(1.0)
    upper = profile.thickness_at_wake(0.0)
    target = _DELTA99_MOMENTUM_FRACTION * theta

    def residual(d: float) -> float:
        return profile.integrals(d)[0] - target

    if residual(lower) > 0.0:
        raise ModelRangeError("required thickness is below the Pi=1 bound (Pi > 1)")
    if residual(upper) < 0.0:
        raise ModelRangeError("required thickness is above the Pi=0 bound (Pi < 0)")
    root = brentq(residual, lower, upper, xtol=np.finfo(float).tiny, rtol=1e-6)
    return root, profile.integrals(root)[1]


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
    viscosity_model: DynamicViscosityModel = AIR_VISCOSITY,
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
        Boundary-layer edge dynamic viscosity in Pa s. Compressibility
        corrections use this same edge-viscosity value.
    regime:
        Explicit laminar, turbulent, or transitional model selection.
    turbulent_correlation:
        Skin-friction correlation used wherever the boundary layer is
        turbulent. The default is the Schlichting logarithmic correlation.
        This option is ignored for turbulent portions when
        ``compressibility_correction=VAN_DRIEST_II`` (Willems equations), or
        ``VAN_DRIEST_II_WALL_TEMPERATURE`` (Schlichting pair).
    transition_reynolds:
        Reynolds number at transition. Required only for ``TRANSITIONAL``;
        there is deliberately no default transition location.
    compressibility_correction:
        ``NONE``, ``ECKERT``, ``VAN_DRIEST_II``, or
        ``VAN_DRIEST_II_WALL_TEMPERATURE``. All three compressible modes solve
        laminar similarity equations. Turbulent Eckert friction at ``Re*`` gives
        coefficients on the reference dynamic pressure ``rho* U_e**2/2``;
        they are converted to the edge basis with ``rho*/rho_e = T_e/T*``.
    mach, edge_temperature:
        Required for a compressibility correction. Temperature is in K.
    wall_temperature:
        Wall temperature in K. If omitted for a compressibility correction,
        an adiabatic wall at the regime-appropriate recovery temperature is
        used.
    prandtl_number:
        Constant Prandtl number used for the recovery factors.
    gas, viscosity_model:
        Perfect-gas and dynamic-viscosity models used by compressibility
        corrections. ``viscosity_model`` evaluates reference or wall
        viscosity; it does not replace ``edge_dynamic_viscosity``.

    Notes
    -----
    ``boundary_layer_thickness`` is the conventional 99-percent thickness.
    With ``NONE``, turbulent thicknesses use the one-fifth-power correlation
    and one-seventh-power profile approximation. ``ECKERT`` and
    ``VAN_DRIEST_II`` return NaN for all three thicknesses at turbulent
    points, including post-transition points, because these approximations
    are inaccurate for compressible flow.

    ``VAN_DRIEST_II_WALL_TEMPERATURE`` multiplies the Schlichting friction
    pair divided by the VD2 friction factor by ``(Tw/Tr)**0.064``. Its
    momentum thickness is ``distance * average_skin_friction_coefficient/2``
    on the edge basis (mixed average for transitional flow). A Van Driest
    velocity transformation with Walz temperature relation, kappa=0.41,
    B=5.2, and 0 <= Pi <= 1 determines delta99 and displacement thickness
    through delta99, matching profile momentum thickness to 0.978 times the
    returned, full-edge momentum thickness. An inadmissible closure raises
    ModelRangeError with the broadcast array index and wake-bound direction.
    External validation covers smooth perfect-gas ZPG flow, Me=0.3--13.3,
    Tw/Taw=0.07--1, freestream T=47--271 K, Re_x,e about 1e6--1e9 (Me<2 only
    to 1.6e8). Against SU2 SA, the delta99 difference falls by 5--6 points per
    decade of Re_x,e (about -9 to -12% at 1e9 on cold supersonic walls).
    Supply local edge states; extrapolation beyond that envelope is unverified.

    Skin-friction coefficients use edge dynamic pressure for every mode.
    ``drag_per_unit_width`` is accumulated from the leading edge on one side.
    The integral momentum relation assumes constant edge and wall conditions
    and zero pressure gradient.

    Compressible laminar cf, theta, delta*, and delta99 come from the
    Illingworth/Levy--Lees similarity equations ``(C f'')' + f f'' = 0`` and
    ``(C g'/Pr)' + f g' + (gamma-1)*Me**2*C*f''**2 = 0``, with
    g=T/Te and C=mu(g*Te)/(mu_model(Te)*g). The solver uses conservative
    finite differences, Newton iteration, grid refinement and Richardson
    extrapolation, with residual-controlled solve_bvp fallback. Agreement
    with the independent reference is within 0.1 percent; displacement near
    zero uses an absolute coefficient tolerance of 5e-5. Numerical accuracy
    does not establish physical validity.

    Omitted wall temperature imposes zero wall heat flux, g'(0)=0.
    Recovery temperature is the solved adiabatic wall temperature, even for
    a prescribed isothermal wall. Effective Reynolds number is the edge Re_x.
    theta=x*cf_local and average_cf=2*cf_local follow from the ZPG integral
    balance. Negative displacement thickness on strongly cooled walls is
    retained. If supplied mu_e differs from mu_model(Te), all four normalized
    coefficients scale by sqrt(mu_model(Te)/mu_e).

    These estimates assume a perfect gas, constant cp and Pr, zero pressure gradient,
    and no induced pressure; at Taw roughly above 2000 K, dissociation makes
    them indicative only. In compressible modes, ApplicabilityWarning is
    emitted once if any laminar point has chi_bar=Me**3*sqrt(Cw/Re_x)>1,
    where Cw=(Te/Tw)*mu(Tw)/mu_e uses the laminar wall temperature and the
    selected viscosity model. The warning reports the maximum and point count.
    Turbulent points and NONE are excluded from this warning.
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
    turbulent_friction_reynolds = np.asarray(reynolds, dtype=np.float64)
    turbulent_thickness_reynolds = np.asarray(reynolds, dtype=np.float64)
    wall_corrected = (
        compressibility_correction
        is CompressibilityCorrection.VAN_DRIEST_II_WALL_TEMPERATURE
    )
    vd2 = wall_corrected or (
        compressibility_correction is CompressibilityCorrection.VAN_DRIEST_II
    )
    vd2_friction = _wall_temperature_friction if wall_corrected else _willems_friction
    turbulent_friction_factor = np.ones_like(reynolds)
    turbulent_reynolds_factor = np.ones_like(reynolds)
    # Reference-to-edge density ratios that convert Eckert coefficients from
    # the reference dynamic-pressure basis to the edge basis.
    laminar_coefficient = np.full_like(reynolds, 0.664)
    laminar = _laminar_correlations(x, reynolds)
    turbulent_density_ratio = np.ones_like(reynolds)
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
        laminar_recovery = np.zeros_like(reynolds)
        laminar_wall = np.zeros_like(reynolds)
        if regime is not BoundaryLayerRegime.TURBULENT:
            # Transitional turbulent points also need the laminar coefficient
            # at transition to accumulate drag from the leading edge.
            states = np.column_stack(
                (
                    mach_values.ravel(),
                    edge_temperature_values.ravel(),
                    np.full(reynolds.size, np.nan)
                    if specified_wall is None
                    else (specified_wall / edge_temperature_values).ravel(),
                    np.full(reynolds.size, prandtl_number),
                    np.full(reynolds.size, gas.heat_capacity_ratio),
                )
            )
            if specified_wall is None:
                solved = similarity_coefficients(
                    states, viscosity_model, _input_shape=reynolds.shape
                )
                adiabatic = solved
            else:
                adiabatic_states = states.copy()
                adiabatic_states[:, 2] = np.nan
                combined = similarity_coefficients(
                    np.concatenate((states, adiabatic_states)),
                    viscosity_model,
                    _input_shape=reynolds.shape,
                )
                solved, adiabatic = np.split(combined, 2)
            model_mu = np.asarray(
                viscosity_model.dynamic_viscosity(edge_temperature_values)
            )
            coefficients = (
                solved[:, :4].reshape((*reynolds.shape, 4))
                * np.sqrt(model_mu / viscosity)[..., None]
            )
            laminar_coefficient = coefficients[..., 0]
            length_scale = x / np.sqrt(reynolds)
            local = laminar_coefficient / np.sqrt(reynolds)
            laminar = (
                coefficients[..., 1] * length_scale,
                coefficients[..., 2] * length_scale,
                coefficients[..., 3] * length_scale,
                local,
                2 * local,
            )
            laminar_recovery = (
                adiabatic[:, 4].reshape(reynolds.shape) * edge_temperature_values
            )
            laminar_wall = (
                laminar_recovery if specified_wall is None else specified_wall
            )
        if regime is BoundaryLayerRegime.LAMINAR:
            recovery_temperature = laminar_recovery
            used_wall_temperature = laminar_wall
        elif compressibility_correction is CompressibilityCorrection.ECKERT:
            (
                turbulent_friction_reynolds,
                turbulent_recovery,
                turbulent_wall,
                turbulent_density_ratio,
            ) = _eckert_effective_reynolds(
                reynolds,
                edge_temperature_values,
                np.asarray(viscosity, dtype=np.float64),
                mach_values,
                specified_wall,
                prandtl_number=prandtl_number,
                gas=gas,
                viscosity_model=viscosity_model,
            )
            turbulent_thickness_reynolds = turbulent_friction_reynolds
            turbulent_reynolds_factor = turbulent_friction_reynolds / reynolds
            recovery_temperature = np.where(
                turbulent_mask, turbulent_recovery, laminar_recovery
            )
            used_wall_temperature = np.where(
                turbulent_mask, turbulent_wall, laminar_wall
            )
        else:
            van_driest = _van_driest_ii_state(
                reynolds,
                edge_temperature_values,
                np.asarray(viscosity, dtype=np.float64),
                mach_values,
                specified_wall,
                prandtl_number=prandtl_number,
                gas=gas,
                viscosity_model=viscosity_model,
            )
            turbulent_friction_reynolds = van_driest.friction_reynolds
            turbulent_thickness_reynolds = van_driest.thickness_reynolds
            turbulent_friction_factor = van_driest.friction_factor
            if wall_corrected:
                weight = (
                    van_driest.wall_temperature / van_driest.recovery_temperature
                ) ** _VD2_WALL_TEMPERATURE_EXPONENT
                turbulent_friction_factor = turbulent_friction_factor / weight
            turbulent_reynolds_factor = van_driest.reynolds_factor
            recovery_temperature = np.where(
                turbulent_mask,
                van_driest.recovery_temperature,
                laminar_recovery,
            )
            used_wall_temperature = np.where(
                turbulent_mask,
                van_driest.wall_temperature,
                laminar_wall,
            )

    if compressibility_correction is not CompressibilityCorrection.NONE:
        _warn_viscous_interaction(
            reynolds,
            mach_values,
            edge_temperature_values,
            viscosity,
            laminar_wall,
            turbulent_mask,
            viscosity_model,
        )
    if regime is BoundaryLayerRegime.LAMINAR:
        selected = laminar
    elif regime is BoundaryLayerRegime.TURBULENT:
        if vd2:
            turbulent_thickness = _turbulent_thicknesses(
                x, turbulent_thickness_reynolds
            )
            local_i, average_i = vd2_friction(turbulent_friction_reynolds)
            turbulent = (
                *turbulent_thickness,
                local_i / turbulent_friction_factor,
                average_i / turbulent_friction_factor,
            )
        else:
            turbulent = _edge_basis_friction(
                _turbulent_correlations(
                    x, turbulent_friction_reynolds, turbulent_correlation
                ),
                turbulent_density_ratio,
            )
        selected = turbulent
    else:
        assert transition is not None
        turbulent_thickness_scale = turbulent_thickness_reynolds / reynolds
        turbulent_thickness = _turbulent_thicknesses(
            x,
            turbulent_thickness_scale * np.maximum(reynolds, transition),
        )
        if vd2:
            friction_reynolds = turbulent_reynolds_factor * np.maximum(
                reynolds, transition
            )
            local_i, average_i = vd2_friction(friction_reynolds)
            turbulent = (
                *turbulent_thickness,
                local_i / turbulent_friction_factor,
                average_i / turbulent_friction_factor,
            )
        else:
            local, average = _turbulent_friction(
                turbulent_reynolds_factor * np.maximum(reynolds, transition),
                turbulent_correlation,
            )
            turbulent = (
                *turbulent_thickness,
                turbulent_density_ratio * local,
                turbulent_density_ratio * average,
            )
        laminar_average = 2 * laminar_coefficient / np.sqrt(transition)
        if vd2:
            _, turbulent_average_i = vd2_friction(
                turbulent_reynolds_factor * transition
            )
            turbulent_average = turbulent_average_i / turbulent_friction_factor
        else:
            _, turbulent_average_reference = _turbulent_friction(
                turbulent_reynolds_factor * transition,
                turbulent_correlation,
            )
            turbulent_average = turbulent_density_ratio * turbulent_average_reference
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
    effective_reynolds = np.where(
        turbulent_mask, turbulent_friction_reynolds, laminar_reynolds
    )
    dynamic_pressure = 0.5 * density * velocity**2
    shear = dynamic_pressure * local_friction
    drag = dynamic_pressure * x * average_friction
    if compressibility_correction in (
        CompressibilityCorrection.ECKERT,
        CompressibilityCorrection.VAN_DRIEST_II,
    ):
        thickness = np.where(turbulent_mask, np.nan, thickness)
        displacement = np.where(turbulent_mask, np.nan, displacement)
        momentum = np.where(turbulent_mask, np.nan, momentum)
    elif wall_corrected:
        thickness = np.array(thickness, copy=True)
        displacement = np.array(displacement, copy=True)
        momentum = np.where(turbulent_mask, x * average_friction / 2.0, momentum)
        assert used_wall_temperature is not None
        for index in np.ndindex(reynolds.shape):
            if not turbulent_mask[index]:
                continue
            try:
                thickness[index], displacement[index] = _wall_temperature_thickness(
                    float(momentum[index]),
                    float(velocity[index]),
                    float(density[index]),
                    float(broadcast[5][index]),
                    float(used_wall_temperature[index]),
                    float(shear[index]),
                    prandtl_number=prandtl_number,
                    gas=gas,
                    viscosity_model=viscosity_model,
                )
            except ModelRangeError as error:
                raise ModelRangeError(
                    f"thickness at array index {index}: {error}"
                ) from error

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
