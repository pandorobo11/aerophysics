"""Frozen, variable-specific-heat centered expansion at constant enthalpy."""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from aerophysics.exceptions import ExpansionConvergenceError, ModelRangeError
from aerophysics.isentropic import _thermal_properties, _ThermalProperties
from aerophysics.real_gas import HarmonicOscillatorGas
from aerophysics.thermochemistry import ThermallyPerfectGas

type ThermalExpansionGas = ThermallyPerfectGas | HarmonicOscillatorGas

_ROOT_XTOL = 1e-12
_ROOT_RTOL = 4.0 * np.finfo(np.float64).eps


def _properties(temperature: float, gas: ThermalExpansionGas) -> _ThermalProperties:
    minimum, maximum = gas.temperature_range
    if not minimum <= temperature <= maximum:
        raise ModelRangeError(
            f"expansion temperature must be within {minimum:g} K and {maximum:g} K"
        )
    return _thermal_properties(temperature, gas, allow_extrapolation=False)


@dataclass(frozen=True, slots=True)
class ThermalExpansionState:
    downstream_mach: float
    temperature_ratio: float
    pressure_ratio: float
    density_ratio: float
    upstream_angle: float
    available_turn_angle: float


def expansion_state(
    mach: float, turn: float, temperature: float, gas: ThermalExpansionGas
) -> ThermalExpansionState:
    """Integrate turning from the upstream state, without a sonic reference."""
    upstream = _properties(temperature, gas)
    speed_squared = mach**2 * upstream.sound_speed_squared
    minimum, maximum = gas.temperature_range
    unbounded = isinstance(gas, HarmonicOscillatorGas) and (
        gas.applicable_temperature_range is None
    )
    boundaries = (
        sorted(
            {
                boundary
                for species in gas.species
                for boundary in species.thermo.temperature_ranges[1:-1]
            }
        )
        if isinstance(gas, ThermallyPerfectGas)
        else []
    )

    def velocity_squared(state: _ThermalProperties) -> float:
        return speed_squared + 2.0 * (upstream.enthalpy - state.enthalpy)

    def angle(lower: float, upper: float) -> float:
        # x = sqrt(T/T1) removes the integrable 1/sqrt(T) vacuum singularity.
        # QUADPACK never samples the endpoints, including x=0 for unbounded gas.
        def integrand(x: float) -> float:
            state = _properties(temperature * x * x, gas)
            velocity = velocity_squared(state)
            m_squared = velocity / state.sound_speed_squared
            if m_squared < 1.0 - 32.0 * np.finfo(np.float64).eps:
                raise ExpansionConvergenceError("expansion path became subsonic")
            gamma = state.heat_capacity_ratio
            cp = gas.specific_gas_constant * gamma / (gamma - 1.0)
            return float(
                2.0
                * temperature
                * x
                * cp
                / velocity
                * np.sqrt(max(m_squared - 1.0, 0.0))
            )

        points = [np.sqrt(t / temperature) for t in boundaries if lower < t < upper]
        value, error = quad(
            integrand,
            np.sqrt(lower / temperature),
            np.sqrt(upper / temperature),
            points=points,
            epsabs=_ROOT_XTOL,
            epsrel=_ROOT_XTOL,
            limit=100,
        )
        if not np.isfinite(value) or error > 10.0 * _ROOT_XTOL * max(1.0, abs(value)):
            raise ExpansionConvergenceError(
                "expansion angle integration did not converge"
            )
        return float(value)

    available = angle(0.0 if unbounded else minimum, temperature)
    endpoint_tolerance = 32.0 * np.finfo(np.float64).eps * max(1.0, available)
    if unbounded and turn >= available:
        raise ValueError("expansion reaches or exceeds the limiting angle")
    if turn > available + endpoint_tolerance or (minimum == temperature and turn > 0.0):
        raise ModelRangeError(
            "downstream expansion temperature is below the gas temperature range"
        )
    if turn == 0.0:
        downstream_temperature = temperature
    elif not unbounded and abs(turn - available) <= endpoint_tolerance:
        # Adopt the endpoint itself; never pass a same-sign bracket to brentq.
        downstream_temperature = minimum
    else:
        try:
            root = brentq(
                lambda x: angle(temperature * x * x, temperature) - turn,
                0.0 if unbounded else np.sqrt(minimum / temperature),
                1.0,
                xtol=_ROOT_XTOL,
                rtol=_ROOT_RTOL,
            )
        except RuntimeError as error:
            raise ExpansionConvergenceError(
                "expansion temperature root did not converge"
            ) from error
        downstream_temperature = temperature * root * root
    downstream = _properties(downstream_temperature, gas)
    temperature_ratio = downstream_temperature / temperature
    pressure_ratio = float(
        np.exp((downstream.entropy - upstream.entropy) / gas.specific_gas_constant)
    )

    def sonic_residual(t: float) -> float:
        state = _properties(t, gas)
        return (
            state.enthalpy
            - upstream.enthalpy
            + 0.5 * (state.sound_speed_squared - speed_squared)
        )

    upstream_angle = np.nan
    if mach == 1.0:
        upstream_angle = 0.0
    else:
        upper = min(2.0 * temperature, maximum)
        while sonic_residual(upper) < 0.0 and upper < maximum:
            upper = min(2.0 * upper, maximum)
        if sonic_residual(upper) >= 0.0:
            try:
                sonic_temperature = brentq(
                    sonic_residual, temperature, upper, xtol=_ROOT_XTOL, rtol=_ROOT_RTOL
                )
            except RuntimeError as error:
                raise ExpansionConvergenceError(
                    "expansion sonic reference did not converge"
                ) from error
            upstream_angle = angle(temperature, sonic_temperature)

    return ThermalExpansionState(
        downstream_mach=float(
            np.sqrt(velocity_squared(downstream) / downstream.sound_speed_squared)
        ),
        temperature_ratio=temperature_ratio,
        pressure_ratio=pressure_ratio,
        density_ratio=pressure_ratio / temperature_ratio,
        upstream_angle=float(upstream_angle),
        available_turn_angle=available,
    )
