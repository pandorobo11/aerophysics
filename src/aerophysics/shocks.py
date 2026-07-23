"""Normal and oblique shocks for a calorically perfect gas.

Angles are expressed in radians. State ratios use downstream over upstream
static quantities, while ``total_pressure_ratio`` is ``p02/p01``.

References
----------
Ames Research Staff, *Equations, Tables, and Charts for Compressible Flow*,
NACA Report 1135, 1953.
"""

from dataclasses import dataclass
from enum import StrEnum
from typing import Final

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import brentq, minimize_scalar

from aerophysics._array import FloatArray, FloatResult, as_float_array, return_float
from aerophysics.exceptions import NoAttachedShockError
from aerophysics.gas import AIR, PerfectGas

_ROOT_XTOL: Final = 1e-12
_ROOT_RTOL: Final = 4.0 * np.finfo(np.float64).eps
_ROOT_MAXITER: Final = 100


class ShockBranch(StrEnum):
    """Branch of the theta-beta-Mach relation."""

    WEAK = "weak"
    STRONG = "strong"


@dataclass(frozen=True, slots=True)
class NormalShockResult:
    """State ratios across a normal shock."""

    upstream_mach: FloatResult
    downstream_mach: FloatResult
    static_pressure_ratio: FloatResult
    static_density_ratio: FloatResult
    static_temperature_ratio: FloatResult
    total_pressure_ratio: FloatResult


@dataclass(frozen=True, slots=True)
class AttachedShockLimit:
    """Maximum attached-shock deflection and corresponding shock angle."""

    upstream_mach: FloatResult
    deflection_angle: FloatResult
    shock_angle: FloatResult


@dataclass(frozen=True, slots=True)
class ObliqueShockResult:
    """Flow state across an attached oblique shock."""

    upstream_mach: FloatResult
    downstream_mach: FloatResult
    deflection_angle: FloatResult
    shock_angle: FloatResult
    upstream_normal_mach: FloatResult
    downstream_normal_mach: FloatResult
    static_pressure_ratio: FloatResult
    static_density_ratio: FloatResult
    static_temperature_ratio: FloatResult
    total_pressure_ratio: FloatResult


def _validate_supersonic_mach(mach: ArrayLike) -> tuple[FloatArray, bool]:
    values, scalar = as_float_array(mach, name="upstream_mach")
    if np.any(values < 1.0):
        raise ValueError("upstream_mach must be greater than or equal to one")
    return values, scalar


def _normal_shock_arrays(
    upstream_mach: FloatArray, gas: PerfectGas
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]:
    gamma = gas.heat_capacity_ratio
    mach_squared = upstream_mach**2
    downstream_mach = np.sqrt(
        (1.0 + 0.5 * (gamma - 1.0) * mach_squared)
        / (gamma * mach_squared - 0.5 * (gamma - 1.0))
    )
    pressure_ratio = 1.0 + 2.0 * gamma * (mach_squared - 1.0) / (gamma + 1.0)
    density_ratio = (gamma + 1.0) * mach_squared / ((gamma - 1.0) * mach_squared + 2.0)
    temperature_ratio = pressure_ratio / density_ratio
    total_pressure_ratio = (
        (gamma + 1.0) * mach_squared / ((gamma - 1.0) * mach_squared + 2.0)
    ) ** (gamma / (gamma - 1.0)) * (
        (gamma + 1.0) / (2.0 * gamma * mach_squared - (gamma - 1.0))
    ) ** (1.0 / (gamma - 1.0))
    return (
        downstream_mach,
        pressure_ratio,
        density_ratio,
        temperature_ratio,
        total_pressure_ratio,
    )


def normal_shock(upstream_mach: ArrayLike, gas: PerfectGas = AIR) -> NormalShockResult:
    """Return perfect-gas normal-shock state ratios."""
    mach, scalar = _validate_supersonic_mach(upstream_mach)
    downstream, pressure, density, temperature, total_pressure = _normal_shock_arrays(
        mach, gas
    )

    def output(values: FloatArray) -> FloatResult:
        return return_float(values, scalar=scalar)

    return NormalShockResult(
        upstream_mach=output(mach),
        downstream_mach=output(downstream),
        static_pressure_ratio=output(pressure),
        static_density_ratio=output(density),
        static_temperature_ratio=output(temperature),
        total_pressure_ratio=output(total_pressure),
    )


def _broadcast_mach_angle(
    upstream_mach: ArrayLike,
    angle: ArrayLike,
    *,
    angle_name: str,
) -> tuple[FloatArray, FloatArray, bool]:
    mach, mach_scalar = _validate_supersonic_mach(upstream_mach)
    angle_values, angle_scalar = as_float_array(angle, name=angle_name)
    try:
        mach, angle_values = np.broadcast_arrays(mach, angle_values)
    except ValueError as error:
        raise ValueError(
            f"upstream_mach and {angle_name} must be broadcastable"
        ) from error
    return mach, angle_values, mach_scalar and angle_scalar


def _theta_from_beta_scalar(mach: float, beta: float, gas: PerfectGas) -> float:
    gamma = gas.heat_capacity_ratio
    numerator = 2.0 * (mach**2 * np.sin(beta) ** 2 - 1.0)
    denominator = np.tan(beta) * (mach**2 * (gamma + np.cos(2.0 * beta)) + 2.0)
    return float(np.arctan2(numerator, denominator))


def theta_from_shock_angle(
    upstream_mach: ArrayLike,
    shock_angle: ArrayLike,
    gas: PerfectGas = AIR,
) -> FloatResult:
    """Return flow-deflection angle from shock angle and upstream Mach."""
    mach, beta, scalar = _broadcast_mach_angle(
        upstream_mach, shock_angle, angle_name="shock_angle"
    )
    mach_angle = np.arcsin(1.0 / mach)
    if np.any((beta < mach_angle) | (beta > 0.5 * np.pi)):
        raise ValueError("shock_angle must lie between the Mach angle and pi/2")
    gamma = gas.heat_capacity_ratio
    numerator = 2.0 * (mach**2 * np.sin(beta) ** 2 - 1.0)
    denominator = np.tan(beta) * (mach**2 * (gamma + np.cos(2.0 * beta)) + 2.0)
    theta = np.maximum(np.arctan2(numerator, denominator), 0.0)
    return return_float(theta, scalar=scalar)


def _attached_limit_scalar(mach: float, gas: PerfectGas) -> tuple[float, float]:
    mach_angle = float(np.arcsin(1.0 / mach))
    if mach == 1.0:
        return 0.0, 0.5 * np.pi
    optimum = minimize_scalar(
        lambda beta: -_theta_from_beta_scalar(mach, float(beta), gas),
        bounds=(mach_angle, 0.5 * np.pi),
        method="bounded",
        options={"xatol": _ROOT_XTOL, "maxiter": _ROOT_MAXITER},
    )
    beta = float(optimum.x)
    return _theta_from_beta_scalar(mach, beta, gas), beta


def maximum_attached_deflection(
    upstream_mach: ArrayLike, gas: PerfectGas = AIR
) -> AttachedShockLimit:
    """Return the maximum deflection permitting an attached oblique shock."""
    mach, scalar = _validate_supersonic_mach(upstream_mach)
    theta = np.empty_like(mach)
    beta = np.empty_like(mach)
    for index, value in np.ndenumerate(mach):
        theta[index], beta[index] = _attached_limit_scalar(float(value), gas)

    def output(values: FloatArray) -> FloatResult:
        return return_float(values, scalar=scalar)

    return AttachedShockLimit(
        upstream_mach=output(mach),
        deflection_angle=output(theta),
        shock_angle=output(beta),
    )


def _shock_angle_scalar(
    mach: float,
    theta: float,
    branch: ShockBranch,
    gas: PerfectGas,
) -> float:
    mach_angle = float(np.arcsin(1.0 / mach))
    theta_max, beta_peak = _attached_limit_scalar(mach, gas)
    tolerance = 1e-11
    if theta > theta_max + tolerance:
        raise NoAttachedShockError(
            f"no attached shock for Mach {mach:g} and deflection {theta:g} rad"
        )
    if abs(theta - theta_max) <= tolerance:
        return beta_peak
    if theta == 0.0:
        return mach_angle if branch is ShockBranch.WEAK else 0.5 * np.pi

    def residual(beta: float) -> float:
        return _theta_from_beta_scalar(mach, beta, gas) - theta

    lower, upper = (
        (mach_angle, beta_peak)
        if branch is ShockBranch.WEAK
        else (beta_peak, 0.5 * np.pi)
    )
    return float(
        brentq(
            residual,
            lower,
            upper,
            xtol=_ROOT_XTOL,
            rtol=_ROOT_RTOL,
            maxiter=_ROOT_MAXITER,
        )
    )


def shock_angle(
    upstream_mach: ArrayLike,
    deflection_angle: ArrayLike,
    branch: ShockBranch,
    gas: PerfectGas = AIR,
) -> FloatResult:
    """Solve the theta-beta-Mach relation on a selected branch."""
    if not isinstance(branch, ShockBranch):
        raise ValueError("branch must be a ShockBranch value")
    mach, theta, scalar = _broadcast_mach_angle(
        upstream_mach, deflection_angle, angle_name="deflection_angle"
    )
    if np.any(theta < 0.0):
        raise ValueError("deflection_angle must be non-negative")
    beta = np.empty_like(mach)
    for index, value in np.ndenumerate(mach):
        beta[index] = _shock_angle_scalar(
            float(value), float(theta[index]), branch, gas
        )
    return return_float(beta, scalar=scalar)


def oblique_shock(
    upstream_mach: ArrayLike,
    deflection_angle: ArrayLike,
    branch: ShockBranch = ShockBranch.WEAK,
    gas: PerfectGas = AIR,
) -> ObliqueShockResult:
    """Return the state across an attached oblique shock."""
    mach, theta, scalar = _broadcast_mach_angle(
        upstream_mach, deflection_angle, angle_name="deflection_angle"
    )
    beta = np.asarray(shock_angle(mach, theta, branch, gas), dtype=np.float64)
    upstream_normal = mach * np.sin(beta)
    (
        downstream_normal,
        pressure,
        density,
        temperature,
        total_pressure,
    ) = _normal_shock_arrays(upstream_normal, gas)
    downstream = downstream_normal / np.sin(beta - theta)

    def output(values: FloatArray) -> FloatResult:
        return return_float(np.asarray(values, dtype=np.float64), scalar=scalar)

    return ObliqueShockResult(
        upstream_mach=output(mach),
        downstream_mach=output(downstream),
        deflection_angle=output(theta),
        shock_angle=output(beta),
        upstream_normal_mach=output(upstream_normal),
        downstream_normal_mach=output(downstream_normal),
        static_pressure_ratio=output(pressure),
        static_density_ratio=output(density),
        static_temperature_ratio=output(temperature),
        total_pressure_ratio=output(total_pressure),
    )


def supersonic_pitot_pressure_ratio(
    upstream_mach: ArrayLike, gas: PerfectGas = AIR
) -> FloatResult:
    """Return post-shock total pressure over upstream static pressure."""
    mach, scalar = _validate_supersonic_mach(upstream_mach)
    gamma = gas.heat_capacity_ratio
    mach_squared = mach**2
    ratio = ((gamma + 1.0) * mach_squared / 2.0) ** (gamma / (gamma - 1.0)) * (
        (gamma + 1.0) / (2.0 * gamma * mach_squared - (gamma - 1.0))
    ) ** (1.0 / (gamma - 1.0))
    return return_float(ratio, scalar=scalar)


__all__ = [
    "AttachedShockLimit",
    "NormalShockResult",
    "ObliqueShockResult",
    "ShockBranch",
    "maximum_attached_deflection",
    "normal_shock",
    "oblique_shock",
    "shock_angle",
    "supersonic_pitot_pressure_ratio",
    "theta_from_shock_angle",
]
