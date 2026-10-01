"""Frozen thermally perfect Taylor--Maccoll flow at zero incidence.

Velocities are scaled by free-stream speed, angles are radians, and static
states remain in the gas temperature range. The shock fixes entropy; total
enthalpy closes the smooth conical flow without a constant-gamma substitution.
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

from aerophysics._array import FloatArray
from aerophysics._thermal_shocks import (
    ThermalShockGas,
    _hugoniot,
    normal_state,
    properties,
)
from aerophysics.exceptions import (
    ModelRangeError,
    NoAttachedShockError,
    ShockConvergenceError,
)


@dataclass(frozen=True, slots=True)
class ThermalConicalState:
    cone_half_angle: float
    post_shock_mach: float
    surface_mach: float
    pressure_ratio: float
    density_ratio: float
    temperature_ratio: float
    total_pressure_ratio: float


@lru_cache(maxsize=4096)
def surface_state(
    mach: float, beta: float, temperature: float, gas: ThermalShockGas
) -> ThermalConicalState | None:
    """Integrate a candidate shock to the tangent-flow surface, or range limit."""
    upstream = properties(temperature, gas)
    mu = float(np.arcsin(1.0 / mach))
    if beta <= mu + np.finfo(float).eps:
        return ThermalConicalState(0.0, mach, mach, 1.0, 1.0, 1.0, 1.0)
    shock = normal_state(float(mach * np.sin(beta)), temperature, gas)
    t2 = temperature * shock.temperature_ratio
    downstream = properties(t2, gas)
    speed_squared = mach**2 * upstream.sound_speed_squared
    initial = np.array([np.cos(beta), -np.sin(beta) / shock.density_ratio])
    initial_squared = float(initial @ initial)
    post_shock_mach = float(
        np.sqrt(speed_squared * initial_squared / downstream.sound_speed_squared)
    )
    maximum = gas.temperature_range[1]
    finite_range = maximum < 1e100
    maximum_enthalpy = properties(maximum, gas).enthalpy if finite_range else np.inf

    def static_temperature(velocity: FloatArray) -> float:
        target = downstream.enthalpy + 0.5 * speed_squared * (
            initial_squared - float(velocity @ velocity)
        )
        # RK trial stages can cross a terminal boundary before event location.
        # Evaluate only at in-range temperatures for those stages; the enthalpy
        # event below stops the physical trajectory at Tmax. No extrapolated
        # thermodynamic state is evaluated or returned.
        if target <= upstream.enthalpy:
            return temperature
        if target >= maximum_enthalpy:
            return maximum
        upper = min(maximum, max(t2, 2.0 * temperature))
        while properties(upper, gas).enthalpy < target:
            upper = min(maximum, 2.0 * upper)
        return float(
            brentq(
                lambda t: properties(t, gas).enthalpy - target,
                temperature,
                upper,
                xtol=1e-10,
                rtol=1e-14,
            )
        )

    def rhs(angle: float, velocity: FloatArray) -> FloatArray:
        radial, polar = velocity
        state = properties(static_temperature(velocity), gas)
        a_squared = state.sound_speed_squared / speed_squared
        derivative = (
            radial * polar**2 - a_squared * (2.0 * radial + polar / np.tan(angle))
        ) / (a_squared - polar**2)
        return np.array([polar, derivative])

    def cone_surface(_angle: float, velocity: FloatArray) -> float:
        return float(velocity[1])

    def range_boundary(_angle: float, velocity: FloatArray) -> float:
        return (maximum_enthalpy - downstream.enthalpy) / speed_squared - 0.5 * (
            initial_squared - float(velocity @ velocity)
        )

    cone_surface.terminal = True  # type: ignore[attr-defined]
    range_boundary.terminal = True  # type: ignore[attr-defined]
    solution = solve_ivp(
        rhs,
        (beta, 1e-8),
        initial,
        method="DOP853",
        events=(cone_surface, range_boundary) if finite_range else cone_surface,
        rtol=1e-10,
        atol=1e-12,
    )
    times, states = solution.t_events, solution.y_events
    if times is not None and finite_range and times[1].size:
        raise ModelRangeError(
            "conical flow temperature exceeds the gas temperature range"
        )
    if not solution.success or times is None or states is None or not times[0].size:
        return None
    surface_velocity = states[0][0]
    # Tmax can be crossed twice in one RK step that passes the surface: the
    # temperature reaches its maximum at V_theta=0 and falls beyond it. An
    # endpoint sign-change event alone can therefore miss the range crossing.
    # Check the physical surface enthalpy before adopting a boundary T.
    if finite_range:
        roundoff = (
            4.0
            * np.finfo(float).eps
            * max(
                1.0,
                (abs(maximum_enthalpy) + abs(downstream.enthalpy)) / speed_squared,
            )
        )
        if range_boundary(float(times[0][0]), surface_velocity) < -roundoff:
            raise ModelRangeError(
                "conical flow temperature exceeds the gas temperature range"
            )
    ts = static_temperature(surface_velocity)
    surface = properties(ts, gas)
    pressure_ratio = shock.pressure_ratio * np.exp(
        (surface.entropy - downstream.entropy) / gas.specific_gas_constant
    )
    return ThermalConicalState(
        float(times[0][0]),
        post_shock_mach,
        float(
            np.sqrt(
                speed_squared
                * float(surface_velocity @ surface_velocity)
                / surface.sound_speed_squared
            )
        ),
        float(pressure_ratio),
        float(pressure_ratio * temperature / ts),
        ts / temperature,
        shock.total_pressure_ratio,
    )


@lru_cache(maxsize=512)
def cone_limit(
    mach: float, temperature: float, gas: ThermalShockGas
) -> tuple[float, float, bool]:
    """Return the available cone-angle peak, beta and physical-peak flag."""
    upstream = properties(temperature, gas)
    mu = float(np.arcsin(1.0 / mach))
    upper = 0.5 * np.pi - 1e-7
    maximum = gas.temperature_range[1]
    if temperature == maximum:
        raise ModelRangeError(
            "gas temperature range contains no compressive conical flow"
        )
    if maximum < 1e100:
        normal_squared, _ = _hugoniot(maximum, temperature, upstream, gas)
        upper = min(
            upper,
            float(
                np.arcsin(
                    min(
                        1.0,
                        np.sqrt(normal_squared / upstream.sound_speed_squared) / mach,
                    )
                )
            ),
        )
    if upper <= mu:
        raise ModelRangeError(
            "gas temperature range contains no compressive conical flow"
        )
    beta_values = np.linspace(mu, upper, 33)
    angles = [0.0]
    available = [mu]
    truncated = False
    for beta in beta_values[1:]:
        try:
            state = surface_state(mach, float(beta), temperature, gas)
        except ModelRangeError:
            # Find the surface-temperature endpoint, not just the shock's Tmax.
            lower, boundary = available[-1], float(beta)
            for _ in range(40):
                middle = 0.5 * (lower + boundary)
                try:
                    middle_state = surface_state(mach, middle, temperature, gas)
                except ModelRangeError:
                    boundary = middle
                else:
                    if middle_state is None:
                        raise ShockConvergenceError(
                            "Taylor-Maccoll integration failed at range boundary"
                        )
                    lower = middle
            state = surface_state(mach, lower, temperature, gas)
            if state is None:
                raise ShockConvergenceError(
                    "Taylor-Maccoll integration failed at range boundary"
                ) from None
            available.append(lower)
            angles.append(state.cone_half_angle)
            truncated = True
            break
        available.append(float(beta))
        angles.append(-np.inf if state is None else state.cone_half_angle)
    peak = int(np.argmax(angles))
    if peak == 0:
        raise ShockConvergenceError(
            "Taylor-Maccoll integration could not find an attached shock"
        )

    def objective(beta: float) -> float:
        state = surface_state(mach, beta, temperature, gas)
        return 1.0 if state is None else -state.cone_half_angle

    optimum = minimize_scalar(
        objective,
        bounds=(available[peak - 1], available[min(peak + 1, len(available) - 1)]),
        method="bounded",
        options={"xatol": 1e-12, "maxiter": 100},
    )
    if not optimum.success:
        raise ShockConvergenceError("thermally perfect cone-angle maximization failed")
    beta_peak = float(optimum.x)
    state = surface_state(mach, beta_peak, temperature, gas)
    if state is None:
        raise ShockConvergenceError(
            "Taylor-Maccoll integration failed at the attached limit"
        )
    # The last sample can straddle a physical peak just inside the range.
    # Near a peak, the angle drop is quadratic and can be smaller than the
    # angle-root tolerance even with an in-range maximum. Use beta separation
    # instead, allowing for bounded minimization's sqrt(eps) relative stopping
    # term and its xatol. A resolved interior maximum is a physical peak.
    if peak == len(angles) - 1 and truncated:
        beta_uncertainty = 4.0 * (np.sqrt(np.finfo(float).eps) * abs(beta_peak) + 1e-12)
        if available[-1] - beta_peak <= beta_uncertainty:
            return angles[-1], available[-1], False
    return state.cone_half_angle, beta_peak, True


def conical_state(
    mach: float, angle: float, temperature: float, gas: ThermalShockGas
) -> tuple[float, ThermalConicalState]:
    """Solve the weak conical branch inside the available static-state range."""
    mu = float(np.arcsin(1.0 / mach))
    if angle == 0.0:
        state = surface_state(mach, mu, temperature, gas)
        assert state is not None
        return mu, state
    peak, beta_peak, full_peak = cone_limit(mach, temperature, gas)
    if angle > peak + 1e-10:
        if not full_peak:
            raise ModelRangeError(
                "requested conical shock exceeds the gas temperature range"
            )
        raise NoAttachedShockError(
            f"no attached conical shock for Mach {mach:g} "
            f"and cone half-angle {angle:g} rad"
        )
    if abs(angle - peak) <= 1e-10:
        beta = beta_peak
    else:

        def residual(beta: float) -> float:
            state = surface_state(mach, beta, temperature, gas)
            if state is None:
                raise ShockConvergenceError(
                    "Taylor-Maccoll integration failed during root solving"
                )
            return state.cone_half_angle - angle

        beta, root_result = brentq(
            residual,
            mu,
            beta_peak,
            xtol=1e-12,
            rtol=1e-14,
            full_output=True,
            disp=False,
        )
        if not root_result.converged:
            raise ShockConvergenceError(
                "thermally perfect conical angle root did not converge"
            )
        if abs(residual(beta)) > 1e-10:
            # Slender cones have a steep angle residual near the Mach wave.
            beta, root_result = brentq(
                residual,
                mu,
                beta_peak,
                xtol=np.finfo(float).eps,
                rtol=4.0 * np.finfo(float).eps,
                full_output=True,
                disp=False,
            )
            if not root_result.converged:
                raise ShockConvergenceError(
                    "thermally perfect conical angle root did not converge"
                )
    state = surface_state(mach, beta, temperature, gas)
    if state is None:
        raise ShockConvergenceError(
            "Taylor-Maccoll integration failed for the conical shock"
        )
    if abs(state.cone_half_angle - angle) > 1e-10:
        raise ShockConvergenceError("conical shock angle is below numerical resolution")
    return beta, state
