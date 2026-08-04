"""Normal, oblique, and conical shocks for a calorically perfect gas.

Angles are expressed in radians. Normal- and oblique-shock state ratios use
downstream over upstream static quantities. Conical-shock static ratios use
cone surface over free stream. ``total_pressure_ratio`` is always post-shock
over upstream total pressure.

References
----------
Ames Research Staff, *Equations, Tables, and Charts for Compressible Flow*,
NACA Report 1135, 1953.
Sims, *Tables for Supersonic Flow Around Right Circular Cones at Zero Angle
of Attack*, NASA SP-3004, 1964.
"""

from dataclasses import dataclass
from enum import StrEnum
from functools import lru_cache
from typing import Final

import numpy as np
from numpy.typing import ArrayLike
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

from aerophysics._array import FloatArray, FloatResult, as_float_array, return_float
from aerophysics.exceptions import NoAttachedShockError
from aerophysics.gas import AIR, PerfectGas

_ROOT_XTOL: Final = 1e-12
_ROOT_RTOL: Final = 4.0 * np.finfo(np.float64).eps
_ROOT_MAXITER: Final = 100
_ODE_RTOL: Final = 1e-10
_ODE_ATOL: Final = 1e-12
_CONE_ANGLE_TOLERANCE: Final = 1e-10
_CONE_AXIS_FLOOR: Final = 1e-8


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


@dataclass(frozen=True, slots=True)
class AttachedConicalShockLimit:
    """Maximum cone half-angle permitting an attached conical shock."""

    upstream_mach: FloatResult
    cone_half_angle: FloatResult
    shock_angle: FloatResult


@dataclass(frozen=True, slots=True)
class ConicalShockResult:
    """Surface state behind an attached axisymmetric conical shock."""

    upstream_mach: FloatResult
    cone_half_angle: FloatResult
    shock_angle: FloatResult
    post_shock_mach: FloatResult
    surface_mach: FloatResult
    surface_pressure_ratio: FloatResult
    surface_density_ratio: FloatResult
    surface_temperature_ratio: FloatResult
    total_pressure_ratio: FloatResult


@dataclass(frozen=True, slots=True)
class _ConicalSurfaceState:
    cone_half_angle: float
    surface_mach: float
    post_shock_mach: float
    total_pressure_ratio: float


class _ConeSurfaceEvent:
    """Locate the ray where the polar velocity becomes tangent to the cone."""

    terminal: bool = True
    direction: float = 0.0

    def __call__(self, _angle: float, velocity: FloatArray) -> float:
        return float(velocity[1])


_CONE_SURFACE_EVENT = _ConeSurfaceEvent()


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


def _mach_from_limiting_velocity(speed_squared: float, gas: PerfectGas) -> float:
    gamma = gas.heat_capacity_ratio
    sound_speed_squared = 0.5 * (gamma - 1.0) * (1.0 - speed_squared)
    if sound_speed_squared <= 0.0:
        raise RuntimeError("Taylor-Maccoll integration left the physical state space")
    return float(np.sqrt(speed_squared / sound_speed_squared))


def _taylor_maccoll_rhs(
    angle: float, velocity: FloatArray, gas: PerfectGas
) -> FloatArray:
    radial, polar = velocity
    gamma = gas.heat_capacity_ratio
    sound_speed_squared = 0.5 * (gamma - 1.0) * (1.0 - radial**2 - polar**2)
    denominator = sound_speed_squared - polar**2
    polar_derivative = (
        radial * polar**2 - sound_speed_squared * (2.0 * radial + polar / np.tan(angle))
    ) / denominator
    return np.asarray([polar, polar_derivative], dtype=np.float64)


def _conical_surface_state_scalar(
    mach: float, shock_angle_value: float, gas: PerfectGas
) -> _ConicalSurfaceState | None:
    mach_angle = float(np.arcsin(1.0 / mach))
    if shock_angle_value <= mach_angle + np.finfo(np.float64).eps:
        return _ConicalSurfaceState(0.0, mach, mach, 1.0)

    gamma = gas.heat_capacity_ratio
    limiting_velocity = mach / np.sqrt(mach**2 + 2.0 / (gamma - 1.0))
    upstream_normal_mach = mach * np.sin(shock_angle_value)
    _, _, density_ratio, _, total_pressure_ratio = _normal_shock_arrays(
        np.asarray(upstream_normal_mach, dtype=np.float64), gas
    )
    radial = limiting_velocity * np.cos(shock_angle_value)
    polar = -limiting_velocity * np.sin(shock_angle_value) / float(density_ratio)
    initial = np.asarray([radial, polar], dtype=np.float64)
    post_shock_mach = _mach_from_limiting_velocity(float(initial @ initial), gas)

    solution = solve_ivp(
        lambda angle, velocity: _taylor_maccoll_rhs(angle, velocity, gas),
        (shock_angle_value, _CONE_AXIS_FLOOR),
        initial,
        method="DOP853",
        events=_CONE_SURFACE_EVENT,
        rtol=_ODE_RTOL,
        atol=_ODE_ATOL,
    )
    event_times = solution.t_events
    event_states = solution.y_events
    if (
        not solution.success
        or event_times is None
        or event_states is None
        or not event_times[0].size
    ):
        return None
    cone_half_angle = float(event_times[0][0])
    surface_velocity = event_states[0][0]
    surface_mach = _mach_from_limiting_velocity(
        float(surface_velocity @ surface_velocity), gas
    )
    return _ConicalSurfaceState(
        cone_half_angle=cone_half_angle,
        surface_mach=surface_mach,
        post_shock_mach=post_shock_mach,
        total_pressure_ratio=float(total_pressure_ratio),
    )


@lru_cache(maxsize=512)
def _attached_conical_limit_scalar(mach: float, gas: PerfectGas) -> tuple[float, float]:
    mach_angle = float(np.arcsin(1.0 / mach))
    beta_values = np.linspace(
        mach_angle + 1e-7,
        0.5 * np.pi - 1e-7,
        33,
        dtype=np.float64,
    )
    cone_angles = np.asarray(
        [
            -np.inf if state is None else state.cone_half_angle
            for state in (
                _conical_surface_state_scalar(mach, float(beta), gas)
                for beta in beta_values
            )
        ],
        dtype=np.float64,
    )
    peak_index = int(np.argmax(cone_angles))
    if not np.isfinite(cone_angles[peak_index]):
        raise RuntimeError(
            "Taylor-Maccoll integration could not find an attached shock"
        )
    lower = float(beta_values[max(peak_index - 1, 0)])
    upper = float(beta_values[min(peak_index + 1, beta_values.size - 1)])

    def objective(beta: float) -> float:
        state = _conical_surface_state_scalar(mach, beta, gas)
        return 1.0 if state is None else -state.cone_half_angle

    optimum = minimize_scalar(
        objective,
        bounds=(lower, upper),
        method="bounded",
        options={"xatol": _ROOT_XTOL, "maxiter": _ROOT_MAXITER},
    )
    beta_peak = float(optimum.x)
    peak_state = _conical_surface_state_scalar(mach, beta_peak, gas)
    if peak_state is None:
        raise RuntimeError("Taylor-Maccoll integration failed at the attached limit")
    return peak_state.cone_half_angle, beta_peak


def _validate_conical_mach(upstream_mach: ArrayLike) -> tuple[FloatArray, bool]:
    mach, scalar = _validate_supersonic_mach(upstream_mach)
    if np.any(mach <= 1.0):
        raise ValueError("upstream_mach must be greater than one")
    return mach, scalar


def maximum_attached_cone_angle(
    upstream_mach: ArrayLike, gas: PerfectGas = AIR
) -> AttachedConicalShockLimit:
    """Return the largest half-angle permitting an attached conical shock."""
    mach, scalar = _validate_conical_mach(upstream_mach)
    cone_angle = np.empty_like(mach)
    beta = np.empty_like(mach)
    for index, value in np.ndenumerate(mach):
        cone_angle[index], beta[index] = _attached_conical_limit_scalar(
            float(value), gas
        )

    def output(values: FloatArray) -> FloatResult:
        return return_float(values, scalar=scalar)

    return AttachedConicalShockLimit(
        upstream_mach=output(mach),
        cone_half_angle=output(cone_angle),
        shock_angle=output(beta),
    )


def _conical_shock_scalar(
    mach: float, cone_half_angle: float, gas: PerfectGas
) -> tuple[float, _ConicalSurfaceState]:
    mach_angle = float(np.arcsin(1.0 / mach))
    if cone_half_angle == 0.0:
        return mach_angle, _ConicalSurfaceState(0.0, mach, mach, 1.0)

    maximum_angle, beta_peak = _attached_conical_limit_scalar(mach, gas)
    if cone_half_angle > maximum_angle + _CONE_ANGLE_TOLERANCE:
        raise NoAttachedShockError(
            f"no attached conical shock for Mach {mach:g} and cone half-angle "
            f"{cone_half_angle:g} rad"
        )
    if abs(cone_half_angle - maximum_angle) <= _CONE_ANGLE_TOLERANCE:
        state = _conical_surface_state_scalar(mach, beta_peak, gas)
        if state is None:
            raise RuntimeError(
                "Taylor-Maccoll integration failed at the attached limit"
            )
        return beta_peak, state

    def residual(beta: float) -> float:
        state = _conical_surface_state_scalar(mach, beta, gas)
        if state is None:
            raise RuntimeError("Taylor-Maccoll integration failed during root solving")
        return state.cone_half_angle - cone_half_angle

    beta = float(
        brentq(
            residual,
            mach_angle,
            beta_peak,
            xtol=_ROOT_XTOL,
            rtol=_ROOT_RTOL,
            maxiter=_ROOT_MAXITER,
        )
    )
    state = _conical_surface_state_scalar(mach, beta, gas)
    if state is None:
        raise RuntimeError("Taylor-Maccoll integration failed for the conical shock")
    return beta, state


def conical_shock(
    upstream_mach: ArrayLike,
    cone_half_angle: ArrayLike,
    gas: PerfectGas = AIR,
) -> ConicalShockResult:
    """Return the cone-surface state behind an attached Taylor-Maccoll shock."""
    mach, mach_scalar = _validate_conical_mach(upstream_mach)
    angle, angle_scalar = as_float_array(cone_half_angle, name="cone_half_angle")
    try:
        mach, angle = np.broadcast_arrays(mach, angle)
    except ValueError as error:
        raise ValueError(
            "upstream_mach and cone_half_angle must be broadcastable"
        ) from error
    if np.any((angle < 0.0) | (angle >= 0.5 * np.pi)):
        raise ValueError("cone_half_angle must be between zero and pi/2")
    scalar = mach_scalar and angle_scalar

    beta = np.empty_like(mach)
    post_shock_mach = np.empty_like(mach)
    surface_mach = np.empty_like(mach)
    surface_pressure = np.empty_like(mach)
    surface_density = np.empty_like(mach)
    surface_temperature = np.empty_like(mach)
    total_pressure = np.empty_like(mach)
    gamma = gas.heat_capacity_ratio

    for index, value in np.ndenumerate(mach):
        beta_value, state = _conical_shock_scalar(
            float(value), float(angle[index]), gas
        )
        beta[index] = beta_value
        post_shock_mach[index] = state.post_shock_mach
        surface_mach[index] = state.surface_mach
        total_pressure[index] = state.total_pressure_ratio
        upstream_factor = 1.0 + 0.5 * (gamma - 1.0) * float(value) ** 2
        surface_factor = 1.0 + 0.5 * (gamma - 1.0) * state.surface_mach**2
        temperature_ratio = upstream_factor / surface_factor
        pressure_ratio = state.total_pressure_ratio * temperature_ratio ** (
            gamma / (gamma - 1.0)
        )
        surface_temperature[index] = temperature_ratio
        surface_pressure[index] = pressure_ratio
        surface_density[index] = pressure_ratio / temperature_ratio

    def output(values: FloatArray) -> FloatResult:
        return return_float(np.asarray(values, dtype=np.float64), scalar=scalar)

    return ConicalShockResult(
        upstream_mach=output(mach),
        cone_half_angle=output(angle),
        shock_angle=output(beta),
        post_shock_mach=output(post_shock_mach),
        surface_mach=output(surface_mach),
        surface_pressure_ratio=output(surface_pressure),
        surface_density_ratio=output(surface_density),
        surface_temperature_ratio=output(surface_temperature),
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
    "AttachedConicalShockLimit",
    "AttachedShockLimit",
    "ConicalShockResult",
    "NormalShockResult",
    "ObliqueShockResult",
    "ShockBranch",
    "conical_shock",
    "maximum_attached_cone_angle",
    "maximum_attached_deflection",
    "normal_shock",
    "oblique_shock",
    "shock_angle",
    "supersonic_pitot_pressure_ratio",
    "theta_from_shock_angle",
]
