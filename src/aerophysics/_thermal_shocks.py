"""Frozen ideal-gas Rankine--Hugoniot solver and shock polar.

The temperature Hugoniot removes the zero-strength root before solving.
No constant-gamma substitution or thermodynamic extrapolation is used.
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from aerophysics.exceptions import ModelRangeError, NoAttachedShockError
from aerophysics.isentropic import _thermal_properties, _ThermalProperties
from aerophysics.real_gas import HarmonicOscillatorGas
from aerophysics.thermochemistry import ThermallyPerfectGas

type ThermalShockGas = ThermallyPerfectGas | HarmonicOscillatorGas


def properties(temperature: float, gas: ThermalShockGas) -> _ThermalProperties:
    """Validate all thermal models before evaluating their shared properties."""
    minimum, maximum = gas.temperature_range
    if not minimum <= temperature <= maximum:
        raise ModelRangeError(
            f"shock temperature must be within {minimum:g} K and {maximum:g} K"
        )
    return _thermal_properties(temperature, gas, allow_extrapolation=False)


@dataclass(frozen=True, slots=True)
class ThermalNormalState:
    downstream_mach: float
    pressure_ratio: float
    density_ratio: float
    temperature_ratio: float
    total_pressure_ratio: float


def _enthalpy_jump(
    temperature: float,
    upstream_temperature: float,
    upstream: _ThermalProperties,
    downstream: _ThermalProperties,
    gas: ThermalShockGas,
) -> float:
    delta = temperature - upstream_temperature
    # Avoid subtracting nearly equal enthalpies for vanishing shocks. Use
    # Simpson integration only within a single polynomial region; across a
    # region boundary retain the database's actual enthalpy offsets.
    crosses_boundary = isinstance(gas, ThermallyPerfectGas) and any(
        upstream_temperature <= boundary < temperature
        for species in gas.species
        for boundary in species.thermo.temperature_ranges[1:-1]
    )
    if delta < 1e-4 * upstream_temperature and not crosses_boundary:
        middle = properties(upstream_temperature + 0.5 * delta, gas)
        r = gas.specific_gas_constant

        def cp(state: _ThermalProperties) -> float:
            gamma = state.heat_capacity_ratio
            return r * gamma / (gamma - 1.0)

        return delta * (cp(upstream) + 4.0 * cp(middle) + cp(downstream)) / 6.0
    return downstream.enthalpy - upstream.enthalpy


def _hugoniot(
    temperature: float,
    upstream_temperature: float,
    upstream: _ThermalProperties,
    gas: ThermalShockGas,
) -> tuple[float, float]:
    """Return required normal speed squared and velocity ratio at a given T2."""
    if temperature == upstream_temperature:
        return upstream.sound_speed_squared, 1.0
    downstream = properties(temperature, gas)
    rt = gas.specific_gas_constant * upstream_temperature
    tau_minus_one = (temperature - upstream_temperature) / upstream_temperature
    q = (
        _enthalpy_jump(temperature, upstream_temperature, upstream, downstream, gas)
        / rt
    )
    c = 2.0 * q + 2.0 - tau_minus_one
    # d = 1 - u_n2/u_n1, from energy + momentum + p=rho*R*T.
    d = 4.0 * (q - tau_minus_one) / (c + np.sqrt(c * c - 8.0 * (q - tau_minus_one)))
    if not 0.0 < d < 1.0:
        raise ModelRangeError("thermally perfect shock Hugoniot is non-physical")
    speed_squared = rt * (tau_minus_one + d) / (d * (1.0 - d))
    return float(speed_squared), float(1.0 - d)


@lru_cache(maxsize=4096)
def normal_state(
    normal_mach: float, upstream_temperature: float, gas: ThermalShockGas
) -> ThermalNormalState:
    upstream = properties(upstream_temperature, gas)
    if normal_mach <= 1.0:
        return ThermalNormalState(1.0, 1.0, 1.0, 1.0, 1.0)
    speed_squared = normal_mach**2 * upstream.sound_speed_squared

    def residual(temperature: float) -> float:
        required, _ = _hugoniot(temperature, upstream_temperature, upstream, gas)
        return required / speed_squared - 1.0

    maximum = gas.temperature_range[1]
    upper = min(maximum, 2.0 * upstream_temperature)
    while residual(upper) < 0.0:
        if upper == maximum:
            # The polar endpoint is computed from Tmax; its reconstructed
            # normal speed can round outward by a few ulps. Adopt Tmax itself
            # only for a machine-precision residual, never a same-sign bracket.
            if residual(upper) >= -4.0 * np.finfo(float).eps:
                temperature = upper
                break
            raise ModelRangeError(
                "downstream shock temperature exceeds the gas temperature range"
            )
        upper = min(maximum, 2.0 * upper)
    else:
        temperature = float(
            brentq(residual, upstream_temperature, upper, xtol=1e-10, rtol=1e-14)
        )
    _, velocity_ratio = _hugoniot(temperature, upstream_temperature, upstream, gas)
    downstream = properties(temperature, gas)
    temperature_ratio = temperature / upstream_temperature
    density_ratio = 1.0 / velocity_ratio
    pressure_ratio = temperature_ratio * density_ratio
    entropy_jump_over_r = (
        downstream.entropy - upstream.entropy
    ) / gas.specific_gas_constant - np.log(pressure_ratio)
    entropy_roundoff = (
        32.0
        * np.finfo(float).eps
        * max(
            1.0,
            (abs(downstream.entropy) + abs(upstream.entropy))
            / gas.specific_gas_constant,
        )
    )
    if entropy_jump_over_r < -entropy_roundoff:
        raise ModelRangeError("thermally perfect shock has negative entropy production")
    return ThermalNormalState(
        downstream_mach=float(
            np.sqrt(speed_squared / downstream.sound_speed_squared) * velocity_ratio
        ),
        pressure_ratio=pressure_ratio,
        density_ratio=density_ratio,
        temperature_ratio=temperature_ratio,
        total_pressure_ratio=float(np.exp(-max(0.0, entropy_jump_over_r))),
    )


def theta_from_beta(
    mach: float, beta: float, temperature: float, gas: ThermalShockGas
) -> float:
    mu = float(np.arcsin(1.0 / mach))
    if beta == mu or beta == 0.5 * np.pi:
        return 0.0
    state = normal_state(float(mach * np.sin(beta)), temperature, gas)
    return float(beta - np.arctan2(np.sin(beta) / state.density_ratio, np.cos(beta)))


@lru_cache(maxsize=512)
def polar_limit(
    mach: float, temperature: float, gas: ThermalShockGas
) -> tuple[float, float, float, bool]:
    """Return theta peak, beta peak, available beta endpoint, full-peak flag."""
    upstream = properties(temperature, gas)
    mu = float(np.arcsin(1.0 / mach))
    if mach == 1.0:
        return 0.0, 0.5 * np.pi, 0.5 * np.pi, True
    maximum = gas.temperature_range[1]
    # Determine the available polar without demanding that a normal shock
    # fit the database: a weak oblique shock may still lie well inside it.
    if maximum < 1e100:
        maximum_speed_squared, _ = _hugoniot(maximum, temperature, upstream, gas)
        maximum_normal_mach = np.sqrt(
            maximum_speed_squared / upstream.sound_speed_squared
        )
        upper = float(np.arcsin(min(1.0, maximum_normal_mach / mach)))
    else:
        upper = 0.5 * np.pi
    if upper <= mu:
        raise ModelRangeError(
            "gas temperature range contains no compressive shock polar"
        )
    optimum = minimize_scalar(
        lambda beta: -theta_from_beta(mach, float(beta), temperature, gas),
        bounds=(mu, upper),
        method="bounded",
        options={"xatol": 1e-12},
    )
    if not optimum.success:
        raise RuntimeError("thermally perfect shock-polar maximization failed")
    beta = float(optimum.x)
    full_peak = upper == 0.5 * np.pi or upper - beta > 1e-6
    if not full_peak:
        beta = upper
    return theta_from_beta(mach, beta, temperature, gas), beta, upper, full_peak


def beta_from_theta(
    mach: float, theta: float, strong: bool, temperature: float, gas: ThermalShockGas
) -> float:
    properties(temperature, gas)
    mu = float(np.arcsin(1.0 / mach))
    if theta == 0.0:
        beta = 0.5 * np.pi if strong else mu
        normal_state(mach if strong else 1.0, temperature, gas)
        return beta
    peak, beta_peak, upper, full_peak = polar_limit(mach, temperature, gas)
    if theta > peak + 1e-11:
        if not full_peak:
            raise ModelRangeError(
                "requested shock lies beyond the available temperature range"
            )
        raise NoAttachedShockError(
            f"no attached shock for Mach {mach:g} and deflection {theta:g} rad"
        )
    if strong and (
        not full_peak or theta < theta_from_beta(mach, upper, temperature, gas) - 1e-11
    ):
        raise ModelRangeError(
            "strong shock lies beyond the available temperature range"
        )
    if abs(theta - peak) <= 1e-11:
        return beta_peak
    lower, upper = (beta_peak, upper) if strong else (mu, beta_peak)
    return float(
        brentq(
            lambda beta: theta_from_beta(mach, beta, temperature, gas) - theta,
            lower,
            upper,
            xtol=1e-12,
            rtol=1e-14,
        )
    )
