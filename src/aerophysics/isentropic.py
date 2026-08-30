"""Isentropic relations for perfect and Beattie--Bridgeman gases.

State ratios use the total-to-static convention: ``T0/T``, ``p0/p``, and
``rho0/rho``. All inputs and outputs are dimensionless except mass flux,
pressure, and temperature arguments explicitly documented in SI units.

Calorically perfect gases use the closed-form constant-``gamma`` relations.
Thermally perfect gases use temperature-dependent enthalpy and entropy and
therefore require total temperature. Beattie--Bridgeman gases additionally
require total pressure because their properties depend on density.

References
----------
Ames Research Staff, *Equations, Tables, and Charts for Compressible Flow*,
NACA Report 1135, 1953.
"""

import warnings
from collections.abc import Callable
from dataclasses import dataclass
from enum import StrEnum
from typing import Final, Literal

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import brentq

from aerophysics._array import FloatArray, FloatResult, as_float_array, return_float
from aerophysics.exceptions import ApplicabilityWarning, ModelRangeError
from aerophysics.gas import AIR, PerfectGas
from aerophysics.real_gas import (
    BeattieBridgemanGas,
    HarmonicOscillatorGas,
    _ScalarThermodynamicState,
)
from aerophysics.thermochemistry import ThermallyPerfectGas

_ROOT_XTOL: Final = 1e-12
_ROOT_RTOL: Final = 4.0 * np.finfo(np.float64).eps
_ROOT_MAXITER: Final = 100
_MAX_SUPERSONIC_BRACKET: Final = 1e6
_MAX_TEMPERATURE_BRACKET_STEPS: Final = 80
_MAX_BRANCH_LOG_TEMPERATURE_STEP: Final = 1.0e-2
_BRANCH_GAP_SAFETY_FACTOR: Final = 0.25
_BRANCH_GAP_DIFFERENCE_STEP: Final = 1.0e-4
_MAX_BRANCH_CONTINUATION_STEPS: Final = 3000

type _IdealThermalGas = ThermallyPerfectGas | HarmonicOscillatorGas
type IsentropicGasModel = PerfectGas | _IdealThermalGas | BeattieBridgemanGas


class MachBranch(StrEnum):
    """Branch of the area-Mach relation."""

    SUBSONIC = "subsonic"
    SUPERSONIC = "supersonic"


@dataclass(frozen=True, slots=True)
class IsentropicRatios:
    """Total-to-static state ratios at one or more Mach numbers."""

    mach: FloatResult
    total_temperature_ratio: FloatResult
    total_pressure_ratio: FloatResult
    total_density_ratio: FloatResult


@dataclass(frozen=True, slots=True)
class CriticalRatios:
    """Total-to-critical ratios at Mach one."""

    total_temperature_ratio: FloatResult
    total_pressure_ratio: FloatResult
    total_density_ratio: FloatResult


@dataclass(frozen=True, slots=True)
class IsentropicFlowState:
    """Absolute total and static properties for isentropic flow."""

    mach: FloatResult
    total_temperature: FloatResult
    total_pressure: FloatResult
    total_density: FloatResult
    static_temperature: FloatResult
    static_pressure: FloatResult
    static_density: FloatResult
    velocity: FloatResult
    speed_of_sound: FloatResult
    dynamic_pressure: FloatResult
    mass_flux: FloatResult


@dataclass(frozen=True, slots=True)
class IsentropicAnalysis:
    """Fused forward-flow results evaluated from one set of solved states.

    ``state``, ``mass_flux``, and ``choked_mass_flux`` are available when both
    total temperature and total pressure are supplied. ``area_ratio`` uses the
    limiting value ``inf`` at Mach zero; the standalone :func:`area_ratio`
    function continues to require strictly positive Mach numbers.
    """

    ratios: IsentropicRatios
    critical_ratios: CriticalRatios
    area_ratio: FloatResult
    mass_flow_parameter: FloatResult
    critical_mass_flow_parameter: FloatResult
    mass_flux: FloatResult | None
    choked_mass_flux: FloatResult | None
    state: IsentropicFlowState | None


@dataclass(frozen=True, slots=True)
class _ThermalProperties:
    enthalpy: float
    entropy: float
    heat_capacity_ratio: float
    sound_speed_squared: float
    extrapolated: bool


@dataclass(frozen=True, slots=True)
class _ThermalFlowState:
    static_temperature: float
    total_temperature_ratio: float
    total_pressure_ratio: float
    total_density_ratio: float
    mass_flow_parameter: float
    extrapolated: bool


@dataclass(frozen=True, slots=True)
class _RealFlowState:
    total: _ScalarThermodynamicState
    static: _ScalarThermodynamicState
    mach: float
    velocity: float
    total_temperature_ratio: float
    total_pressure_ratio: float
    total_density_ratio: float
    mass_flow_parameter: float


@dataclass(frozen=True, slots=True)
class _AbsoluteFlowArrays:
    total_density: FloatArray
    static_temperature: FloatArray
    static_pressure: FloatArray
    static_density: FloatArray
    velocity: FloatArray
    speed_of_sound: FloatArray


@dataclass(frozen=True, slots=True)
class _AnalysisArrays:
    temperature_ratio: FloatArray
    pressure_ratio: FloatArray
    density_ratio: FloatArray
    mass_flow_parameter: FloatArray
    critical_temperature_ratio: FloatArray
    critical_pressure_ratio: FloatArray
    critical_density_ratio: FloatArray
    critical_mass_flow_parameter: FloatArray
    absolute: _AbsoluteFlowArrays | None


def _validate_mach(
    mach: ArrayLike, *, positive: bool = False
) -> tuple[FloatArray, bool]:
    values, scalar = as_float_array(mach, name="mach")
    lower_invalid = values <= 0.0 if positive else values < 0.0
    if np.any(lower_invalid):
        condition = "greater than zero" if positive else "non-negative"
        raise ValueError(f"mach must be {condition}")
    return values, scalar


def _validate_ratio(ratio: ArrayLike, *, name: str) -> tuple[FloatArray, bool]:
    values, scalar = as_float_array(ratio, name=name)
    if np.any(values < 1.0):
        raise ValueError(f"{name} must be greater than or equal to one")
    return values, scalar


def _temperature_factor(mach: FloatArray, gas: PerfectGas) -> FloatArray:
    return 1.0 + 0.5 * (gas.heat_capacity_ratio - 1.0) * mach**2


def _thermal_properties(
    temperature: float,
    gas: _IdealThermalGas,
    *,
    allow_extrapolation: bool,
) -> _ThermalProperties:
    value = np.asarray(temperature, dtype=np.float64)
    if isinstance(gas, ThermallyPerfectGas):
        validated, _ = gas._validated_temperature(
            value,
            allow_extrapolation=allow_extrapolation,
            warn=False,
        )
        molar_cp, molar_enthalpy, molar_entropy = (
            gas._standard_molar_properties_from_values(validated)
        )
        cp = float(molar_cp) / gas.molar_mass
        enthalpy = float(molar_enthalpy) / gas.molar_mass
        entropy = float(molar_entropy) / gas.molar_mass
        minimum, maximum = gas.temperature_range
        extrapolated = temperature < minimum or temperature > maximum
    else:
        cp_values, cv_values, _, enthalpy_values, entropy_values = (
            gas._properties_from_values(value)
        )
        cp = float(cp_values)
        enthalpy = float(enthalpy_values)
        entropy = float(entropy_values)
        cv = float(cv_values)
        extrapolated = False
    gas_constant = gas.specific_gas_constant
    if isinstance(gas, ThermallyPerfectGas):
        cv = cp - gas_constant
    if not np.isfinite(cp) or not np.isfinite(cv) or cv <= 0.0:
        raise ModelRangeError(
            "thermally perfect gas has non-physical heat capacity at the "
            f"required temperature {temperature:g} K"
        )
    gamma = cp / cv
    sound_speed_squared = gamma * gas_constant * temperature
    if (
        not np.isfinite(gamma)
        or gamma <= 1.0
        or not np.isfinite(sound_speed_squared)
        or sound_speed_squared <= 0.0
    ):
        raise ModelRangeError(
            "thermally perfect gas has non-physical sound speed at the "
            f"required temperature {temperature:g} K"
        )
    return _ThermalProperties(
        enthalpy=enthalpy,
        entropy=entropy,
        heat_capacity_ratio=gamma,
        sound_speed_squared=sound_speed_squared,
        extrapolated=extrapolated,
    )


def _warn_if_extrapolated(gas: _IdealThermalGas, extrapolated: bool) -> None:
    if not extrapolated:
        return
    minimum, maximum = gas.temperature_range
    detail = (
        "; the nearest polynomial region was extrapolated"
        if isinstance(gas, ThermallyPerfectGas)
        else ""
    )
    warnings.warn(
        "isentropic solution uses temperature outside the documented range "
        f"{minimum:g}--{maximum:g} K{detail}",
        ApplicabilityWarning,
        stacklevel=3,
    )


def _check_harmonic_total_temperature(
    gas: HarmonicOscillatorGas,
    total_temperature: FloatArray,
    *,
    allow_extrapolation: bool,
) -> None:
    applicable_range = gas.applicable_temperature_range
    if applicable_range is None:
        return
    minimum, maximum = applicable_range
    outside = bool(
        np.any((total_temperature < minimum) | (total_temperature > maximum))
    )
    if outside and not allow_extrapolation:
        raise ModelRangeError(
            f"total_temperature must be within {minimum:g}--{maximum:g} K"
        )
    if outside:
        warnings.warn(
            "total_temperature is outside the documented harmonic-oscillator "
            f"air range {minimum:g}--{maximum:g} K",
            ApplicabilityWarning,
            stacklevel=3,
        )


def _broadcast_thermal_inputs(
    values: FloatArray,
    values_scalar: bool,
    total_temperature: ArrayLike | None,
    *,
    value_name: str,
) -> tuple[FloatArray, FloatArray, bool]:
    if total_temperature is None:
        raise ValueError("total_temperature is required for a thermally perfect gas")
    temperatures, temperature_scalar = as_float_array(
        total_temperature, name="total_temperature"
    )
    if np.any(temperatures <= 0.0):
        raise ValueError("total_temperature must be greater than zero")
    try:
        broadcast_values, broadcast_temperatures = np.broadcast_arrays(
            values, temperatures
        )
    except ValueError as error:
        raise ValueError(
            f"{value_name} and total_temperature must be broadcastable"
        ) from error
    return (
        np.asarray(broadcast_values, dtype=np.float64),
        np.asarray(broadcast_temperatures, dtype=np.float64),
        values_scalar and temperature_scalar,
    )


def _broadcast_real_inputs(
    values: FloatArray,
    values_scalar: bool,
    total_temperature: ArrayLike | None,
    total_pressure: ArrayLike | None,
    *,
    value_name: str,
) -> tuple[FloatArray, FloatArray, FloatArray, bool]:
    if total_temperature is None:
        raise ValueError("total_temperature is required for a Beattie--Bridgeman gas")
    if total_pressure is None:
        raise ValueError("total_pressure is required for a Beattie--Bridgeman gas")
    temperatures, temperature_scalar = as_float_array(
        total_temperature, name="total_temperature"
    )
    pressures, pressure_scalar = as_float_array(total_pressure, name="total_pressure")
    if np.any(temperatures <= 0.0):
        raise ValueError("total_temperature must be greater than zero")
    if np.any(pressures <= 0.0):
        raise ValueError("total_pressure must be greater than zero")
    try:
        broadcast_values, broadcast_temperatures, broadcast_pressures = (
            np.broadcast_arrays(values, temperatures, pressures)
        )
    except ValueError as error:
        raise ValueError(
            f"{value_name}, total_temperature, and total_pressure must be broadcastable"
        ) from error
    return (
        np.asarray(broadcast_values, dtype=np.float64),
        np.asarray(broadcast_temperatures, dtype=np.float64),
        np.asarray(broadcast_pressures, dtype=np.float64),
        values_scalar and temperature_scalar and pressure_scalar,
    )


def _check_real_total_conditions(
    gas: BeattieBridgemanGas,
    total_temperature: FloatArray,
    total_pressure: FloatArray,
    *,
    allow_extrapolation: bool,
) -> bool:
    outside = gas._check_applicability(
        total_temperature,
        total_pressure,
        allow_extrapolation=allow_extrapolation,
        warn=False,
    )
    return outside


def _check_real_static_conditions(
    gas: BeattieBridgemanGas,
    total_outside: bool,
    static_temperature: FloatArray,
    static_pressure: FloatArray,
    *,
    allow_extrapolation: bool,
) -> None:
    """Validate all states and emit one extrapolation warning per public call."""
    static_outside = gas._check_applicability(
        static_temperature,
        static_pressure,
        allow_extrapolation=allow_extrapolation,
        warn=False,
    )
    if total_outside or static_outside:
        warnings.warn(
            f"isentropic solution is outside the {gas._applicability_description()}",
            ApplicabilityWarning,
            stacklevel=3,
        )


def _find_lower_temperature(
    residual: Callable[[float], float],
    total_temperature: float,
    gas: _IdealThermalGas,
    *,
    allow_extrapolation: bool,
) -> float:
    minimum = gas.temperature_range[0]
    if not allow_extrapolation and isinstance(gas, ThermallyPerfectGas):
        lower = minimum
        lower_residual = float(residual(lower))
        if lower >= total_temperature or (
            lower_residual < 0.0
            and not np.isclose(lower_residual, 0.0, rtol=0.0, atol=1e-8)
        ):
            raise ModelRangeError(
                "isentropic solution requires temperature below the fitted "
                f"range minimum {minimum:g} K"
            )
        return lower

    lower = min(minimum, 0.5 * total_temperature)
    floor = max(np.finfo(np.float64).tiny, total_temperature * 1e-12)
    for _ in range(_MAX_TEMPERATURE_BRACKET_STEPS):
        if float(residual(lower)) >= 0.0:
            return lower
        if lower <= floor:
            break
        lower = max(0.5 * lower, floor)
    raise ModelRangeError(
        "could not bracket a positive-temperature thermally perfect isentropic solution"
    )


def _thermal_flow_state(
    mach: float,
    total_temperature: float,
    gas: _IdealThermalGas,
    *,
    allow_extrapolation: bool,
) -> _ThermalFlowState:
    total = _thermal_properties(
        total_temperature, gas, allow_extrapolation=allow_extrapolation
    )
    if mach == 0.0:
        return _ThermalFlowState(
            static_temperature=total_temperature,
            total_temperature_ratio=1.0,
            total_pressure_ratio=1.0,
            total_density_ratio=1.0,
            mass_flow_parameter=0.0,
            extrapolated=total.extrapolated,
        )

    def residual(temperature: float) -> float:
        static = _thermal_properties(
            temperature, gas, allow_extrapolation=allow_extrapolation
        )
        return (
            total.enthalpy
            - static.enthalpy
            - 0.5 * mach**2 * static.sound_speed_squared
        )

    lower = _find_lower_temperature(
        residual,
        total_temperature,
        gas,
        allow_extrapolation=allow_extrapolation,
    )
    if np.isclose(residual(lower), 0.0, rtol=0.0, atol=1e-8):
        static_temperature = lower
    else:
        static_temperature = float(
            brentq(
                residual,
                lower,
                total_temperature,
                xtol=_ROOT_XTOL,
                rtol=_ROOT_RTOL,
                maxiter=_ROOT_MAXITER,
            )
        )
    static = _thermal_properties(
        static_temperature, gas, allow_extrapolation=allow_extrapolation
    )
    log_pressure_ratio = (total.entropy - static.entropy) / gas.specific_gas_constant
    pressure_ratio = float(np.exp(log_pressure_ratio))
    density_ratio = pressure_ratio * static_temperature / total_temperature
    temperature_ratio = total_temperature / static_temperature
    parameter = (
        mach * np.sqrt(static.heat_capacity_ratio * temperature_ratio) / pressure_ratio
    )
    values = (pressure_ratio, density_ratio, temperature_ratio, parameter)
    if not np.all(np.isfinite(values)) or any(value <= 0.0 for value in values):
        raise ModelRangeError("thermally perfect isentropic state is non-physical")
    return _ThermalFlowState(
        static_temperature=static_temperature,
        total_temperature_ratio=temperature_ratio,
        total_pressure_ratio=pressure_ratio,
        total_density_ratio=density_ratio,
        mass_flow_parameter=float(parameter),
        extrapolated=total.extrapolated or static.extrapolated,
    )


def _thermal_flow_states(
    mach: FloatArray,
    total_temperature: FloatArray,
    gas: _IdealThermalGas,
    *,
    allow_extrapolation: bool,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, bool]:
    temperature_ratio = np.empty_like(mach)
    pressure_ratio = np.empty_like(mach)
    density_ratio = np.empty_like(mach)
    parameter = np.empty_like(mach)
    extrapolated = False
    for index in np.ndindex(mach.shape):
        state = _thermal_flow_state(
            float(mach[index]),
            float(total_temperature[index]),
            gas,
            allow_extrapolation=allow_extrapolation,
        )
        temperature_ratio[index] = state.total_temperature_ratio
        pressure_ratio[index] = state.total_pressure_ratio
        density_ratio[index] = state.total_density_ratio
        parameter[index] = state.mass_flow_parameter
        extrapolated = extrapolated or state.extrapolated
    return temperature_ratio, pressure_ratio, density_ratio, parameter, extrapolated


def _real_isentropic_state_at_temperature(
    temperature: float,
    total: _ScalarThermodynamicState,
    gas: BeattieBridgemanGas,
) -> _ScalarThermodynamicState:
    if temperature == total.temperature:
        return total
    density_floor = max(
        np.finfo(np.float64).tiny,
        total.density * 1e-14,
    )

    spinodal = gas._first_spinodal_density(temperature)
    density_upper = total.density
    if spinodal is not None:
        # ``np.roots`` is accurate enough to locate the boundary, but a small
        # relative retreat keeps endpoint and Brent iterates away from roundoff
        # that can otherwise change the sign of dp/drho at the spinodal.
        stable_upper = np.nextafter(spinodal * (1.0 - 1.0e-10), 0.0)
        for _ in range(64):
            if gas._dp_drho_scalar(temperature, stable_upper) > 0.0:
                break
            stable_upper = np.nextafter(stable_upper, 0.0)
        else:
            raise ModelRangeError(
                "could not establish the mechanically stable "
                "Beattie--Bridgeman gas branch"
            )
        density_upper = min(density_upper, stable_upper)
    if density_upper <= density_floor:
        raise ModelRangeError(
            "could not bracket the Beattie--Bridgeman isentropic gas branch"
        )

    def entropy_residual(density: float) -> float:
        state = gas._scalar_state_from_density(temperature, density)
        return state.entropy - total.entropy

    lower_residual = entropy_residual(density_floor)
    upper_residual = entropy_residual(density_upper)
    if lower_residual < 0.0 or upper_residual > 0.0:
        raise ModelRangeError(
            "could not bracket the Beattie--Bridgeman isentropic gas branch"
        )
    if upper_residual == 0.0:
        density = density_upper
    else:
        log_density = float(
            brentq(
                lambda value: entropy_residual(float(np.exp(value))),
                float(np.log(density_floor)),
                float(np.log(density_upper)),
                xtol=_ROOT_XTOL,
                rtol=_ROOT_RTOL,
                maxiter=_ROOT_MAXITER,
            )
        )
        density = float(np.exp(log_density))
    return gas._scalar_state_from_density(temperature, density)


def _real_flow_state(
    mach: float,
    total_temperature: float,
    total_pressure: float,
    gas: BeattieBridgemanGas,
) -> _RealFlowState:
    total = gas._scalar_state(total_temperature, total_pressure)
    if mach == 0.0:
        return _RealFlowState(
            total=total,
            static=total,
            mach=0.0,
            velocity=0.0,
            total_temperature_ratio=1.0,
            total_pressure_ratio=1.0,
            total_density_ratio=1.0,
            mass_flow_parameter=0.0,
        )

    def state_and_residual(
        temperature: float,
    ) -> tuple[_ScalarThermodynamicState, float]:
        static = _real_isentropic_state_at_temperature(temperature, total, gas)
        residual = (
            total.enthalpy - static.enthalpy - 0.5 * mach**2 * static.speed_of_sound**2
        )
        return static, residual

    floor = max(np.finfo(np.float64).tiny, total_temperature * 1e-12)
    upper = total_temperature
    upper_static, upper_residual = state_and_residual(upper)

    def continuation_temperature(
        current_temperature: float, current_static: _ScalarThermodynamicState
    ) -> float:
        """Advance the entropy branch without stepping over its spinodal event."""
        spinodal = gas._first_spinodal_density(current_temperature)
        if spinodal is None:
            log_step = _MAX_BRANCH_LOG_TEMPERATURE_STEP
        else:
            gap = spinodal - current_static.density
            if gap <= 0.0:
                raise ModelRangeError(
                    "no connected Beattie--Bridgeman isentropic solution "
                    "above the non-physical temperature boundary"
                )
            probe_temperature = float(
                np.exp(np.log(current_temperature) - _BRANCH_GAP_DIFFERENCE_STEP)
            )
            try:
                probe_static = _real_isentropic_state_at_temperature(
                    probe_temperature, total, gas
                )
            except ModelRangeError:
                # The infinitesimal continuation probe is itself the first
                # observed branch event; let the outer boundary refinement
                # locate it from the last connected state.
                return probe_temperature
            probe_spinodal = gas._first_spinodal_density(probe_temperature)
            if probe_spinodal is None:
                log_step = _MAX_BRANCH_LOG_TEMPERATURE_STEP
            else:
                probe_gap = probe_spinodal - probe_static.density
                gap_rate = abs((probe_gap - gap) / _BRANCH_GAP_DIFFERENCE_STEP)
                log_step = min(
                    _MAX_BRANCH_LOG_TEMPERATURE_STEP,
                    _BRANCH_GAP_SAFETY_FACTOR
                    * gap
                    / max(gap_rate, np.finfo(np.float64).tiny),
                )
        return max(float(np.exp(np.log(current_temperature) - log_step)), floor)

    lower = continuation_temperature(upper, upper_static)
    for _ in range(_MAX_BRANCH_CONTINUATION_STEPS):
        try:
            lower_static, lower_residual = state_and_residual(lower)
        except ModelRangeError:
            invalid_lower = lower
            connected_lower = upper
            connected_residual = upper_residual
            for _ in range(_MAX_TEMPERATURE_BRACKET_STEPS):
                midpoint = float(
                    np.exp(0.5 * (np.log(connected_lower) + np.log(invalid_lower)))
                )
                try:
                    _, midpoint_residual = state_and_residual(midpoint)
                except ModelRangeError:
                    invalid_lower = midpoint
                else:
                    connected_lower = midpoint
                    connected_residual = midpoint_residual
            if connected_residual < 0.0:
                raise ModelRangeError(
                    "no connected Beattie--Bridgeman isentropic solution "
                    "above the non-physical temperature boundary"
                ) from None
            lower = connected_lower
            lower_residual = connected_residual
            break
        if lower_residual >= 0.0:
            break
        if lower <= floor:
            raise ModelRangeError(
                "could not bracket a positive-temperature "
                "Beattie--Bridgeman isentropic solution"
            )
        upper = lower
        upper_residual = lower_residual
        upper_static = lower_static
        if lower <= floor:
            raise ModelRangeError(
                "could not bracket a positive-temperature "
                "Beattie--Bridgeman isentropic solution"
            )
        next_lower = continuation_temperature(upper, upper_static)
        if not np.isfinite(next_lower) or next_lower <= 0.0 or next_lower >= upper:
            raise ModelRangeError(
                "could not make positive-temperature progress on the "
                "Beattie--Bridgeman isentropic branch"
            )
        lower = next_lower
    else:
        raise ModelRangeError(
            "could not reach a Beattie--Bridgeman isentropic branch boundary"
        )

    log_lower = float(np.log(lower))
    log_upper = float(np.log(upper))

    def temperature_from_log(log_temperature: float) -> float:
        if log_temperature == log_upper:
            return upper
        if log_temperature == log_lower:
            return lower
        return float(np.exp(log_temperature))

    def log_temperature_residual(log_temperature: float) -> float:
        return state_and_residual(temperature_from_log(log_temperature))[1]

    static_temperature = float(
        brentq(
            log_temperature_residual,
            log_lower,
            log_upper,
            xtol=_ROOT_XTOL,
            rtol=_ROOT_RTOL,
            maxiter=_ROOT_MAXITER,
        )
    )
    static_temperature = temperature_from_log(static_temperature)
    static = _real_isentropic_state_at_temperature(static_temperature, total, gas)
    velocity = mach * static.speed_of_sound
    parameter = (
        static.density
        * velocity
        * np.sqrt(gas.specific_gas_constant * total_temperature)
        / total_pressure
    )
    values = (
        total_temperature / static.temperature,
        total_pressure / static.pressure,
        total.density / static.density,
        parameter,
    )
    if not np.all(np.isfinite(values)) or any(value <= 0.0 for value in values):
        raise ModelRangeError("Beattie--Bridgeman isentropic state is non-physical")
    return _RealFlowState(
        total=total,
        static=static,
        mach=mach,
        velocity=velocity,
        total_temperature_ratio=float(values[0]),
        total_pressure_ratio=float(values[1]),
        total_density_ratio=float(values[2]),
        mass_flow_parameter=float(parameter),
    )


def _real_flow_states(
    mach: FloatArray,
    total_temperature: FloatArray,
    total_pressure: FloatArray,
    gas: BeattieBridgemanGas,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]:
    temperature_ratio = np.empty_like(mach)
    pressure_ratio = np.empty_like(mach)
    density_ratio = np.empty_like(mach)
    parameter = np.empty_like(mach)
    static_temperature = np.empty_like(mach)
    static_pressure = np.empty_like(mach)
    for index in np.ndindex(mach.shape):
        state = _real_flow_state(
            float(mach[index]),
            float(total_temperature[index]),
            float(total_pressure[index]),
            gas,
        )
        temperature_ratio[index] = state.total_temperature_ratio
        pressure_ratio[index] = state.total_pressure_ratio
        density_ratio[index] = state.total_density_ratio
        parameter[index] = state.mass_flow_parameter
        static_temperature[index] = state.static.temperature
        static_pressure[index] = state.static.pressure
    return (
        temperature_ratio,
        pressure_ratio,
        density_ratio,
        parameter,
        static_temperature,
        static_pressure,
    )


def _check_real_mach_states(
    gas: BeattieBridgemanGas,
    total_outside: bool,
    mach: FloatArray,
    total_temperature: FloatArray,
    total_pressure: FloatArray,
    *,
    allow_extrapolation: bool,
) -> None:
    """Validate final static states without warning during root iterations."""
    _, _, _, _, static_temperature, static_pressure = _real_flow_states(
        mach, total_temperature, total_pressure, gas
    )
    _check_real_static_conditions(
        gas,
        total_outside,
        static_temperature,
        static_pressure,
        allow_extrapolation=allow_extrapolation,
    )


def _check_real_area_states(
    gas: BeattieBridgemanGas,
    total_outside: bool,
    mach: FloatArray,
    total_temperature: FloatArray,
    total_pressure: FloatArray,
    *,
    allow_extrapolation: bool,
) -> None:
    """Validate requested and critical states with one aggregate warning."""
    _, _, _, _, static_temperature, static_pressure = _real_flow_states(
        mach, total_temperature, total_pressure, gas
    )
    _, _, _, _, critical_temperature, critical_pressure = _real_flow_states(
        np.ones_like(mach), total_temperature, total_pressure, gas
    )
    _check_real_static_conditions(
        gas,
        total_outside,
        np.concatenate((static_temperature.ravel(), critical_temperature.ravel())),
        np.concatenate((static_pressure.ravel(), critical_pressure.ravel())),
        allow_extrapolation=allow_extrapolation,
    )


def _analysis_inputs(
    mach: ArrayLike,
    gas: IsentropicGasModel,
    total_temperature: ArrayLike | None,
    total_pressure: ArrayLike | None,
) -> tuple[FloatArray, FloatArray | None, FloatArray | None, bool]:
    mach_values, mach_scalar = _validate_mach(mach)
    if isinstance(gas, BeattieBridgemanGas) and (
        total_temperature is None or total_pressure is None
    ):
        raise ValueError(
            "total_temperature and total_pressure are required for a "
            "Beattie--Bridgeman gas"
        )
    if isinstance(gas, (ThermallyPerfectGas, HarmonicOscillatorGas)) and (
        total_temperature is None
    ):
        raise ValueError("total_temperature is required for a thermally perfect gas")
    if total_pressure is not None and total_temperature is None:
        raise ValueError(
            "total_temperature is required when total_pressure is specified"
        )

    temperatures: FloatArray | None = None
    pressures: FloatArray | None = None
    temperature_scalar = True
    pressure_scalar = True
    if total_temperature is not None:
        temperatures, temperature_scalar = as_float_array(
            total_temperature, name="total_temperature"
        )
        if np.any(temperatures <= 0.0):
            raise ValueError("total_temperature must be greater than zero")
    if total_pressure is not None:
        pressures, pressure_scalar = as_float_array(
            total_pressure, name="total_pressure"
        )
        if np.any(pressures <= 0.0):
            raise ValueError("total_pressure must be greater than zero")

    try:
        if pressures is not None:
            assert temperatures is not None
            mach_values, temperatures, pressures = np.broadcast_arrays(
                mach_values, temperatures, pressures
            )
        elif temperatures is not None:
            mach_values, temperatures = np.broadcast_arrays(mach_values, temperatures)
    except ValueError as error:
        raise ValueError(
            "mach, total_temperature, and total_pressure must be broadcastable"
        ) from error

    return (
        np.asarray(mach_values, dtype=np.float64),
        (
            np.asarray(temperatures, dtype=np.float64)
            if temperatures is not None
            else None
        ),
        np.asarray(pressures, dtype=np.float64) if pressures is not None else None,
        mach_scalar and temperature_scalar and pressure_scalar,
    )


def _perfect_analysis_arrays(
    mach: FloatArray,
    total_temperature: FloatArray | None,
    total_pressure: FloatArray | None,
    gas: PerfectGas,
) -> _AnalysisArrays:
    gamma = gas.heat_capacity_ratio
    temperature_ratio = _temperature_factor(mach, gas)
    pressure_ratio = temperature_ratio ** (gamma / (gamma - 1.0))
    density_ratio = temperature_ratio ** (1.0 / (gamma - 1.0))
    exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    parameter = np.sqrt(gamma) * mach / temperature_ratio**exponent

    critical_temperature = 0.5 * (gamma + 1.0)
    critical_pressure = critical_temperature ** (gamma / (gamma - 1.0))
    critical_density = critical_temperature ** (1.0 / (gamma - 1.0))
    critical_parameter = np.sqrt(gamma) / critical_temperature**exponent
    critical_temperature_ratio = np.full_like(mach, critical_temperature)
    critical_pressure_ratio = np.full_like(mach, critical_pressure)
    critical_density_ratio = np.full_like(mach, critical_density)
    critical_mass_flow_parameter = np.full_like(mach, critical_parameter)

    absolute = None
    if total_pressure is not None:
        assert total_temperature is not None
        static_temperature = total_temperature / temperature_ratio
        static_pressure = total_pressure / pressure_ratio
        total_density = total_pressure / (gas.specific_gas_constant * total_temperature)
        static_density = static_pressure / (
            gas.specific_gas_constant * static_temperature
        )
        speed_of_sound = np.sqrt(gamma * gas.specific_gas_constant * static_temperature)
        absolute = _AbsoluteFlowArrays(
            total_density=total_density,
            static_temperature=static_temperature,
            static_pressure=static_pressure,
            static_density=static_density,
            velocity=mach * speed_of_sound,
            speed_of_sound=speed_of_sound,
        )

    return _AnalysisArrays(
        temperature_ratio=temperature_ratio,
        pressure_ratio=pressure_ratio,
        density_ratio=density_ratio,
        mass_flow_parameter=parameter,
        critical_temperature_ratio=critical_temperature_ratio,
        critical_pressure_ratio=critical_pressure_ratio,
        critical_density_ratio=critical_density_ratio,
        critical_mass_flow_parameter=critical_mass_flow_parameter,
        absolute=absolute,
    )


def _thermal_analysis_arrays(
    mach: FloatArray,
    total_temperature: FloatArray,
    total_pressure: FloatArray | None,
    gas: _IdealThermalGas,
    *,
    allow_extrapolation: bool,
) -> _AnalysisArrays:
    if isinstance(gas, HarmonicOscillatorGas):
        _check_harmonic_total_temperature(
            gas,
            total_temperature,
            allow_extrapolation=allow_extrapolation,
        )

    temperature_ratio = np.empty_like(mach)
    pressure_ratio = np.empty_like(mach)
    density_ratio = np.empty_like(mach)
    parameter = np.empty_like(mach)
    critical_temperature_ratio = np.empty_like(mach)
    critical_pressure_ratio = np.empty_like(mach)
    critical_density_ratio = np.empty_like(mach)
    critical_parameter = np.empty_like(mach)

    total_density: FloatArray | None = None
    static_temperature: FloatArray | None = None
    static_pressure: FloatArray | None = None
    static_density: FloatArray | None = None
    velocity: FloatArray | None = None
    speed_of_sound: FloatArray | None = None
    if total_pressure is not None:
        total_density = np.empty_like(mach)
        static_temperature = np.empty_like(mach)
        static_pressure = np.empty_like(mach)
        static_density = np.empty_like(mach)
        velocity = np.empty_like(mach)
        speed_of_sound = np.empty_like(mach)

    solved: dict[tuple[float, float], _ThermalFlowState] = {}

    def solve(mach_value: float, temperature: float) -> _ThermalFlowState:
        key = (mach_value, temperature)
        state = solved.get(key)
        if state is None:
            state = _thermal_flow_state(
                mach_value,
                temperature,
                gas,
                allow_extrapolation=allow_extrapolation,
            )
            solved[key] = state
        return state

    extrapolated = False
    for index in np.ndindex(mach.shape):
        total_temperature_value = float(total_temperature[index])
        critical = solve(1.0, total_temperature_value)
        state = solve(float(mach[index]), total_temperature_value)
        temperature_ratio[index] = state.total_temperature_ratio
        pressure_ratio[index] = state.total_pressure_ratio
        density_ratio[index] = state.total_density_ratio
        parameter[index] = state.mass_flow_parameter
        critical_temperature_ratio[index] = critical.total_temperature_ratio
        critical_pressure_ratio[index] = critical.total_pressure_ratio
        critical_density_ratio[index] = critical.total_density_ratio
        critical_parameter[index] = critical.mass_flow_parameter
        extrapolated = extrapolated or state.extrapolated or critical.extrapolated

        if total_pressure is not None:
            assert total_density is not None
            assert static_temperature is not None
            assert static_pressure is not None
            assert static_density is not None
            assert velocity is not None
            assert speed_of_sound is not None
            properties = _thermal_properties(
                state.static_temperature,
                gas,
                allow_extrapolation=allow_extrapolation,
            )
            static_temperature[index] = state.static_temperature
            static_pressure[index] = total_pressure[index] / state.total_pressure_ratio
            total_density[index] = total_pressure[index] / (
                gas.specific_gas_constant * total_temperature[index]
            )
            static_density[index] = static_pressure[index] / (
                gas.specific_gas_constant * static_temperature[index]
            )
            speed_of_sound[index] = np.sqrt(properties.sound_speed_squared)
            velocity[index] = mach[index] * speed_of_sound[index]

    _warn_if_extrapolated(gas, extrapolated)
    absolute = None
    if total_pressure is not None:
        assert total_density is not None
        assert static_temperature is not None
        assert static_pressure is not None
        assert static_density is not None
        assert velocity is not None
        assert speed_of_sound is not None
        absolute = _AbsoluteFlowArrays(
            total_density=total_density,
            static_temperature=static_temperature,
            static_pressure=static_pressure,
            static_density=static_density,
            velocity=velocity,
            speed_of_sound=speed_of_sound,
        )
    return _AnalysisArrays(
        temperature_ratio=temperature_ratio,
        pressure_ratio=pressure_ratio,
        density_ratio=density_ratio,
        mass_flow_parameter=parameter,
        critical_temperature_ratio=critical_temperature_ratio,
        critical_pressure_ratio=critical_pressure_ratio,
        critical_density_ratio=critical_density_ratio,
        critical_mass_flow_parameter=critical_parameter,
        absolute=absolute,
    )


def _real_analysis_arrays(
    mach: FloatArray,
    total_temperature: FloatArray,
    total_pressure: FloatArray,
    gas: BeattieBridgemanGas,
    *,
    allow_extrapolation: bool,
) -> _AnalysisArrays:
    total_outside = _check_real_total_conditions(
        gas,
        total_temperature,
        total_pressure,
        allow_extrapolation=allow_extrapolation,
    )
    temperature_ratio = np.empty_like(mach)
    pressure_ratio = np.empty_like(mach)
    density_ratio = np.empty_like(mach)
    parameter = np.empty_like(mach)
    critical_temperature_ratio = np.empty_like(mach)
    critical_pressure_ratio = np.empty_like(mach)
    critical_density_ratio = np.empty_like(mach)
    critical_parameter = np.empty_like(mach)
    total_density = np.empty_like(mach)
    static_temperature = np.empty_like(mach)
    static_pressure = np.empty_like(mach)
    static_density = np.empty_like(mach)
    velocity = np.empty_like(mach)
    speed_of_sound = np.empty_like(mach)
    critical_static_temperature = np.empty_like(mach)
    critical_static_pressure = np.empty_like(mach)
    solved: dict[tuple[float, float, float], _RealFlowState] = {}

    def solve(mach_value: float, temperature: float, pressure: float) -> _RealFlowState:
        key = (mach_value, temperature, pressure)
        state = solved.get(key)
        if state is None:
            state = _real_flow_state(mach_value, temperature, pressure, gas)
            solved[key] = state
        return state

    for index in np.ndindex(mach.shape):
        total_temperature_value = float(total_temperature[index])
        total_pressure_value = float(total_pressure[index])
        critical = solve(1.0, total_temperature_value, total_pressure_value)
        state = solve(float(mach[index]), total_temperature_value, total_pressure_value)
        temperature_ratio[index] = state.total_temperature_ratio
        pressure_ratio[index] = state.total_pressure_ratio
        density_ratio[index] = state.total_density_ratio
        parameter[index] = state.mass_flow_parameter
        critical_temperature_ratio[index] = critical.total_temperature_ratio
        critical_pressure_ratio[index] = critical.total_pressure_ratio
        critical_density_ratio[index] = critical.total_density_ratio
        critical_parameter[index] = critical.mass_flow_parameter
        total_density[index] = state.total.density
        static_temperature[index] = state.static.temperature
        static_pressure[index] = state.static.pressure
        static_density[index] = state.static.density
        velocity[index] = state.velocity
        speed_of_sound[index] = state.static.speed_of_sound
        critical_static_temperature[index] = critical.static.temperature
        critical_static_pressure[index] = critical.static.pressure

    _check_real_static_conditions(
        gas,
        total_outside,
        np.concatenate(
            (static_temperature.ravel(), critical_static_temperature.ravel())
        ),
        np.concatenate((static_pressure.ravel(), critical_static_pressure.ravel())),
        allow_extrapolation=allow_extrapolation,
    )
    return _AnalysisArrays(
        temperature_ratio=temperature_ratio,
        pressure_ratio=pressure_ratio,
        density_ratio=density_ratio,
        mass_flow_parameter=parameter,
        critical_temperature_ratio=critical_temperature_ratio,
        critical_pressure_ratio=critical_pressure_ratio,
        critical_density_ratio=critical_density_ratio,
        critical_mass_flow_parameter=critical_parameter,
        absolute=_AbsoluteFlowArrays(
            total_density=total_density,
            static_temperature=static_temperature,
            static_pressure=static_pressure,
            static_density=static_density,
            velocity=velocity,
            speed_of_sound=speed_of_sound,
        ),
    )


def _make_isentropic_analysis(
    mach: FloatArray,
    total_temperature: FloatArray | None,
    total_pressure: FloatArray | None,
    gas: IsentropicGasModel,
    values: _AnalysisArrays,
    *,
    scalar: bool,
) -> IsentropicAnalysis:
    def output(result: FloatArray) -> FloatResult:
        return return_float(result, scalar=scalar)

    area = np.full_like(mach, np.inf)
    np.divide(
        values.critical_mass_flow_parameter,
        values.mass_flow_parameter,
        out=area,
        where=values.mass_flow_parameter > 0.0,
    )

    mass_flux_value: FloatResult | None = None
    choked_mass_flux_value: FloatResult | None = None
    state = None
    if total_pressure is not None:
        assert total_temperature is not None
        assert values.absolute is not None
        normalization = np.sqrt(gas.specific_gas_constant * total_temperature)
        mass_flux_value = output(
            total_pressure * values.mass_flow_parameter / normalization
        )
        choked_mass_flux_value = output(
            total_pressure * values.critical_mass_flow_parameter / normalization
        )
        dynamic_pressure = (
            0.5 * values.absolute.static_density * values.absolute.velocity**2
        )
        state = IsentropicFlowState(
            mach=output(mach),
            total_temperature=output(total_temperature),
            total_pressure=output(total_pressure),
            total_density=output(values.absolute.total_density),
            static_temperature=output(values.absolute.static_temperature),
            static_pressure=output(values.absolute.static_pressure),
            static_density=output(values.absolute.static_density),
            velocity=output(values.absolute.velocity),
            speed_of_sound=output(values.absolute.speed_of_sound),
            dynamic_pressure=output(dynamic_pressure),
            mass_flux=output(values.absolute.static_density * values.absolute.velocity),
        )

    return IsentropicAnalysis(
        ratios=IsentropicRatios(
            mach=output(mach),
            total_temperature_ratio=output(values.temperature_ratio),
            total_pressure_ratio=output(values.pressure_ratio),
            total_density_ratio=output(values.density_ratio),
        ),
        critical_ratios=CriticalRatios(
            total_temperature_ratio=output(values.critical_temperature_ratio),
            total_pressure_ratio=output(values.critical_pressure_ratio),
            total_density_ratio=output(values.critical_density_ratio),
        ),
        area_ratio=output(area),
        mass_flow_parameter=output(values.mass_flow_parameter),
        critical_mass_flow_parameter=output(values.critical_mass_flow_parameter),
        mass_flux=mass_flux_value,
        choked_mass_flux=choked_mass_flux_value,
        state=state,
    )


def isentropic_analysis(
    mach: ArrayLike,
    gas: IsentropicGasModel = AIR,
    *,
    total_temperature: ArrayLike | None = None,
    total_pressure: ArrayLike | None = None,
    allow_extrapolation: bool = True,
) -> IsentropicAnalysis:
    """Evaluate coupled forward isentropic quantities with shared flow states.

    Use this fused API when ratios, area, mass-flow, and absolute-state values
    are needed together. Thermally perfect and Beattie--Bridgeman numerical
    states are solved once per distinct input state, while each distinct
    reservoir condition shares one Mach-one critical state.

    ``total_temperature`` is required for thermally perfect gases. A
    Beattie--Bridgeman gas also requires ``total_pressure``. For other gases,
    absolute state and mass-flux results are included only when both reservoir
    quantities are supplied.
    """
    mach_values, temperatures, pressures, scalar = _analysis_inputs(
        mach,
        gas,
        total_temperature,
        total_pressure,
    )
    if isinstance(gas, BeattieBridgemanGas):
        assert temperatures is not None
        assert pressures is not None
        values = _real_analysis_arrays(
            mach_values,
            temperatures,
            pressures,
            gas,
            allow_extrapolation=allow_extrapolation,
        )
    elif isinstance(gas, (ThermallyPerfectGas, HarmonicOscillatorGas)):
        assert temperatures is not None
        values = _thermal_analysis_arrays(
            mach_values,
            temperatures,
            pressures,
            gas,
            allow_extrapolation=allow_extrapolation,
        )
    else:
        values = _perfect_analysis_arrays(
            mach_values,
            temperatures,
            pressures,
            gas,
        )
    return _make_isentropic_analysis(
        mach_values,
        temperatures,
        pressures,
        gas,
        values,
        scalar=scalar,
    )


def isentropic_ratios(
    mach: ArrayLike,
    gas: IsentropicGasModel = AIR,
    *,
    total_temperature: ArrayLike | None = None,
    total_pressure: ArrayLike | None = None,
    allow_extrapolation: bool = True,
) -> IsentropicRatios:
    """Return total-to-static temperature, pressure, and density ratios.

    ``total_temperature`` is required for thermally perfect gases. A
    Beattie--Bridgeman gas also requires ``total_pressure``.
    """
    values, scalar = _validate_mach(mach)
    if isinstance(gas, BeattieBridgemanGas):
        values, temperatures, pressures, scalar = _broadcast_real_inputs(
            values,
            scalar,
            total_temperature,
            total_pressure,
            value_name="mach",
        )
        total_outside = _check_real_total_conditions(
            gas,
            temperatures,
            pressures,
            allow_extrapolation=allow_extrapolation,
        )
        (
            temperature_ratio,
            pressure_ratio,
            density_ratio,
            _,
            static_temperature,
            static_pressure,
        ) = _real_flow_states(values, temperatures, pressures, gas)
        _check_real_static_conditions(
            gas,
            total_outside,
            static_temperature,
            static_pressure,
            allow_extrapolation=allow_extrapolation,
        )
    elif isinstance(gas, (ThermallyPerfectGas, HarmonicOscillatorGas)):
        assert isinstance(gas, (ThermallyPerfectGas, HarmonicOscillatorGas))
        values, temperatures, scalar = _broadcast_thermal_inputs(
            values, scalar, total_temperature, value_name="mach"
        )
        if isinstance(gas, HarmonicOscillatorGas):
            _check_harmonic_total_temperature(
                gas,
                temperatures,
                allow_extrapolation=allow_extrapolation,
            )
        temperature_ratio, pressure_ratio, density_ratio, _, extrapolated = (
            _thermal_flow_states(
                values,
                temperatures,
                gas,
                allow_extrapolation=allow_extrapolation,
            )
        )
        _warn_if_extrapolated(gas, extrapolated)
    else:
        gamma = gas.heat_capacity_ratio
        temperature_ratio = _temperature_factor(values, gas)
        pressure_ratio = temperature_ratio ** (gamma / (gamma - 1.0))
        density_ratio = temperature_ratio ** (1.0 / (gamma - 1.0))

    def output(result: FloatArray) -> FloatResult:
        return return_float(result, scalar=scalar)

    return IsentropicRatios(
        mach=output(values),
        total_temperature_ratio=output(temperature_ratio),
        total_pressure_ratio=output(pressure_ratio),
        total_density_ratio=output(density_ratio),
    )


def _thermal_mach_from_static_temperature(
    static_temperature: float,
    total_temperature: float,
    gas: _IdealThermalGas,
    *,
    allow_extrapolation: bool,
) -> tuple[float, bool]:
    total = _thermal_properties(
        total_temperature, gas, allow_extrapolation=allow_extrapolation
    )
    static = _thermal_properties(
        static_temperature, gas, allow_extrapolation=allow_extrapolation
    )
    enthalpy_difference = total.enthalpy - static.enthalpy
    tolerance = 1e-12 * max(1.0, abs(total.enthalpy), abs(static.enthalpy))
    if enthalpy_difference < -tolerance:
        raise ModelRangeError(
            "thermally perfect gas enthalpy is not monotonic over the required "
            "temperature interval"
        )
    mach = np.sqrt(2.0 * max(0.0, enthalpy_difference) / static.sound_speed_squared)
    if not np.isfinite(mach):
        raise ModelRangeError("thermally perfect isentropic Mach number is invalid")
    return float(mach), total.extrapolated or static.extrapolated


def mach_from_total_temperature_ratio(
    ratio: ArrayLike,
    gas: IsentropicGasModel = AIR,
    *,
    total_temperature: ArrayLike | None = None,
    total_pressure: ArrayLike | None = None,
    allow_extrapolation: bool = True,
) -> FloatResult:
    """Return Mach number from ``T0/T``."""
    values, scalar = _validate_ratio(ratio, name="total_temperature_ratio")
    if isinstance(gas, PerfectGas):
        result = np.sqrt(2.0 * (values - 1.0) / (gas.heat_capacity_ratio - 1.0))
        return return_float(result, scalar=scalar)

    if isinstance(gas, BeattieBridgemanGas):
        values, temperatures, pressures, scalar = _broadcast_real_inputs(
            values,
            scalar,
            total_temperature,
            total_pressure,
            value_name="total_temperature_ratio",
        )
        total_outside = _check_real_total_conditions(
            gas,
            temperatures,
            pressures,
            allow_extrapolation=allow_extrapolation,
        )
        result = np.empty_like(values)
        for index in np.ndindex(values.shape):
            result[index] = _real_mach_from_ratio_scalar(
                float(values[index]),
                float(temperatures[index]),
                float(pressures[index]),
                gas,
                basis="temperature",
            )
        _check_real_mach_states(
            gas,
            total_outside,
            result,
            temperatures,
            pressures,
            allow_extrapolation=allow_extrapolation,
        )
        return return_float(result, scalar=scalar)

    assert isinstance(gas, (ThermallyPerfectGas, HarmonicOscillatorGas))
    values, temperatures, scalar = _broadcast_thermal_inputs(
        values, scalar, total_temperature, value_name="total_temperature_ratio"
    )
    if isinstance(gas, HarmonicOscillatorGas):
        _check_harmonic_total_temperature(
            gas, temperatures, allow_extrapolation=allow_extrapolation
        )
    result = np.empty_like(values)
    extrapolated = False
    for index in np.ndindex(values.shape):
        total = float(temperatures[index])
        result[index], outside = _thermal_mach_from_static_temperature(
            total / float(values[index]),
            total,
            gas,
            allow_extrapolation=allow_extrapolation,
        )
        extrapolated = extrapolated or outside
    _warn_if_extrapolated(gas, extrapolated)
    return return_float(result, scalar=scalar)


def _thermal_static_temperature_from_ratio(
    ratio: float,
    total_temperature: float,
    gas: _IdealThermalGas,
    *,
    basis: Literal["pressure", "density"],
    allow_extrapolation: bool,
) -> tuple[float, bool]:
    total = _thermal_properties(
        total_temperature, gas, allow_extrapolation=allow_extrapolation
    )
    if ratio == 1.0:
        return total_temperature, total.extrapolated
    target = np.log(ratio)

    def residual(temperature: float) -> float:
        static = _thermal_properties(
            temperature, gas, allow_extrapolation=allow_extrapolation
        )
        value = (total.entropy - static.entropy) / gas.specific_gas_constant
        if basis == "density":
            value += np.log(temperature / total_temperature)
        return float(value - target)

    lower = _find_lower_temperature(
        residual,
        total_temperature,
        gas,
        allow_extrapolation=allow_extrapolation,
    )
    if np.isclose(residual(lower), 0.0, rtol=0.0, atol=1e-8):
        temperature = lower
    else:
        temperature = float(
            brentq(
                residual,
                lower,
                total_temperature,
                xtol=_ROOT_XTOL,
                rtol=_ROOT_RTOL,
                maxiter=_ROOT_MAXITER,
            )
        )
    static = _thermal_properties(
        temperature, gas, allow_extrapolation=allow_extrapolation
    )
    return temperature, total.extrapolated or static.extrapolated


def _thermal_mach_from_state_ratio(
    ratio: FloatArray,
    total_temperature: FloatArray,
    gas: _IdealThermalGas,
    *,
    basis: Literal["pressure", "density"],
    allow_extrapolation: bool,
) -> tuple[FloatArray, bool]:
    result = np.empty_like(ratio)
    extrapolated = False
    for index in np.ndindex(ratio.shape):
        total = float(total_temperature[index])
        static, outside = _thermal_static_temperature_from_ratio(
            float(ratio[index]),
            total,
            gas,
            basis=basis,
            allow_extrapolation=allow_extrapolation,
        )
        result[index], mach_outside = _thermal_mach_from_static_temperature(
            static,
            total,
            gas,
            allow_extrapolation=allow_extrapolation,
        )
        extrapolated = extrapolated or outside or mach_outside
    return result, extrapolated


def _real_mach_from_ratio_scalar(
    target: float,
    total_temperature: float,
    total_pressure: float,
    gas: BeattieBridgemanGas,
    *,
    basis: Literal["temperature", "pressure", "density"],
) -> float:
    if target == 1.0:
        return 0.0

    attribute = {
        "temperature": "total_temperature_ratio",
        "pressure": "total_pressure_ratio",
        "density": "total_density_ratio",
    }[basis]

    def residual(mach: float) -> float:
        state = _real_flow_state(mach, total_temperature, total_pressure, gas)
        return float(getattr(state, attribute)) - target

    upper = 1.0
    for _ in range(_MAX_TEMPERATURE_BRACKET_STEPS):
        if residual(upper) >= 0.0:
            return float(
                brentq(
                    residual,
                    0.0,
                    upper,
                    xtol=_ROOT_XTOL,
                    rtol=_ROOT_RTOL,
                    maxiter=_ROOT_MAXITER,
                )
            )
        upper *= 2.0
    raise ModelRangeError(
        f"could not bracket Mach from Beattie--Bridgeman {basis} ratio"
    )


def mach_from_total_pressure_ratio(
    ratio: ArrayLike,
    gas: IsentropicGasModel = AIR,
    *,
    total_temperature: ArrayLike | None = None,
    total_pressure: ArrayLike | None = None,
    allow_extrapolation: bool = True,
) -> FloatResult:
    """Return Mach number from ``p0/p``."""
    values, scalar = _validate_ratio(ratio, name="total_pressure_ratio")
    if isinstance(gas, PerfectGas):
        gamma = gas.heat_capacity_ratio
        temperature_ratio = values ** ((gamma - 1.0) / gamma)
        result = np.sqrt(2.0 * (temperature_ratio - 1.0) / (gamma - 1.0))
        return return_float(result, scalar=scalar)

    if isinstance(gas, BeattieBridgemanGas):
        values, temperatures, pressures, scalar = _broadcast_real_inputs(
            values,
            scalar,
            total_temperature,
            total_pressure,
            value_name="total_pressure_ratio",
        )
        total_outside = _check_real_total_conditions(
            gas,
            temperatures,
            pressures,
            allow_extrapolation=allow_extrapolation,
        )
        result = np.empty_like(values)
        for index in np.ndindex(values.shape):
            result[index] = _real_mach_from_ratio_scalar(
                float(values[index]),
                float(temperatures[index]),
                float(pressures[index]),
                gas,
                basis="pressure",
            )
        _check_real_mach_states(
            gas,
            total_outside,
            result,
            temperatures,
            pressures,
            allow_extrapolation=allow_extrapolation,
        )
        return return_float(result, scalar=scalar)

    assert isinstance(gas, (ThermallyPerfectGas, HarmonicOscillatorGas))
    values, temperatures, scalar = _broadcast_thermal_inputs(
        values, scalar, total_temperature, value_name="total_pressure_ratio"
    )
    if isinstance(gas, HarmonicOscillatorGas):
        _check_harmonic_total_temperature(
            gas, temperatures, allow_extrapolation=allow_extrapolation
        )
    result, extrapolated = _thermal_mach_from_state_ratio(
        values,
        temperatures,
        gas,
        basis="pressure",
        allow_extrapolation=allow_extrapolation,
    )
    _warn_if_extrapolated(gas, extrapolated)
    return return_float(result, scalar=scalar)


def mach_from_total_density_ratio(
    ratio: ArrayLike,
    gas: IsentropicGasModel = AIR,
    *,
    total_temperature: ArrayLike | None = None,
    total_pressure: ArrayLike | None = None,
    allow_extrapolation: bool = True,
) -> FloatResult:
    """Return Mach number from ``rho0/rho``."""
    values, scalar = _validate_ratio(ratio, name="total_density_ratio")
    if isinstance(gas, PerfectGas):
        gamma = gas.heat_capacity_ratio
        temperature_ratio = values ** (gamma - 1.0)
        result = np.sqrt(2.0 * (temperature_ratio - 1.0) / (gamma - 1.0))
        return return_float(result, scalar=scalar)

    if isinstance(gas, BeattieBridgemanGas):
        values, temperatures, pressures, scalar = _broadcast_real_inputs(
            values,
            scalar,
            total_temperature,
            total_pressure,
            value_name="total_density_ratio",
        )
        total_outside = _check_real_total_conditions(
            gas,
            temperatures,
            pressures,
            allow_extrapolation=allow_extrapolation,
        )
        result = np.empty_like(values)
        for index in np.ndindex(values.shape):
            result[index] = _real_mach_from_ratio_scalar(
                float(values[index]),
                float(temperatures[index]),
                float(pressures[index]),
                gas,
                basis="density",
            )
        _check_real_mach_states(
            gas,
            total_outside,
            result,
            temperatures,
            pressures,
            allow_extrapolation=allow_extrapolation,
        )
        return return_float(result, scalar=scalar)

    assert isinstance(gas, (ThermallyPerfectGas, HarmonicOscillatorGas))
    values, temperatures, scalar = _broadcast_thermal_inputs(
        values, scalar, total_temperature, value_name="total_density_ratio"
    )
    if isinstance(gas, HarmonicOscillatorGas):
        _check_harmonic_total_temperature(
            gas, temperatures, allow_extrapolation=allow_extrapolation
        )
    result, extrapolated = _thermal_mach_from_state_ratio(
        values,
        temperatures,
        gas,
        basis="density",
        allow_extrapolation=allow_extrapolation,
    )
    _warn_if_extrapolated(gas, extrapolated)
    return return_float(result, scalar=scalar)


def _area_ratio_formula(mach: FloatArray, gas: PerfectGas) -> FloatArray:
    gamma = gas.heat_capacity_ratio
    factor = (2.0 / (gamma + 1.0)) * _temperature_factor(mach, gas)
    exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    return factor**exponent / mach


def _thermal_area_ratio_scalar(
    mach: float,
    total_temperature: float,
    gas: _IdealThermalGas,
    *,
    allow_extrapolation: bool,
) -> tuple[float, bool]:
    state = _thermal_flow_state(
        mach,
        total_temperature,
        gas,
        allow_extrapolation=allow_extrapolation,
    )
    critical = _thermal_flow_state(
        1.0,
        total_temperature,
        gas,
        allow_extrapolation=allow_extrapolation,
    )
    return (
        critical.mass_flow_parameter / state.mass_flow_parameter,
        state.extrapolated or critical.extrapolated,
    )


def _real_area_ratio_scalar(
    mach: float,
    total_temperature: float,
    total_pressure: float,
    gas: BeattieBridgemanGas,
) -> float:
    state = _real_flow_state(mach, total_temperature, total_pressure, gas)
    critical = _real_flow_state(1.0, total_temperature, total_pressure, gas)
    return critical.mass_flow_parameter / state.mass_flow_parameter


def area_ratio(
    mach: ArrayLike,
    gas: IsentropicGasModel = AIR,
    *,
    total_temperature: ArrayLike | None = None,
    total_pressure: ArrayLike | None = None,
    allow_extrapolation: bool = True,
) -> FloatResult:
    """Return the isentropic area ratio ``A/A*``."""
    values, scalar = _validate_mach(mach, positive=True)
    if isinstance(gas, PerfectGas):
        return return_float(_area_ratio_formula(values, gas), scalar=scalar)

    if isinstance(gas, BeattieBridgemanGas):
        values, temperatures, pressures, scalar = _broadcast_real_inputs(
            values,
            scalar,
            total_temperature,
            total_pressure,
            value_name="mach",
        )
        total_outside = _check_real_total_conditions(
            gas,
            temperatures,
            pressures,
            allow_extrapolation=allow_extrapolation,
        )
        result = np.empty_like(values)
        for index in np.ndindex(values.shape):
            result[index] = _real_area_ratio_scalar(
                float(values[index]),
                float(temperatures[index]),
                float(pressures[index]),
                gas,
            )
        _check_real_area_states(
            gas,
            total_outside,
            values,
            temperatures,
            pressures,
            allow_extrapolation=allow_extrapolation,
        )
        return return_float(result, scalar=scalar)

    assert isinstance(gas, (ThermallyPerfectGas, HarmonicOscillatorGas))
    values, temperatures, scalar = _broadcast_thermal_inputs(
        values, scalar, total_temperature, value_name="mach"
    )
    if isinstance(gas, HarmonicOscillatorGas):
        _check_harmonic_total_temperature(
            gas, temperatures, allow_extrapolation=allow_extrapolation
        )
    result = np.empty_like(values)
    extrapolated = False
    for index in np.ndindex(values.shape):
        result[index], outside = _thermal_area_ratio_scalar(
            float(values[index]),
            float(temperatures[index]),
            gas,
            allow_extrapolation=allow_extrapolation,
        )
        extrapolated = extrapolated or outside
    _warn_if_extrapolated(gas, extrapolated)
    return return_float(result, scalar=scalar)


def _mach_from_area_scalar(target: float, branch: MachBranch, gas: PerfectGas) -> float:
    if target == 1.0:
        return 1.0

    def residual(mach: float) -> float:
        value = _area_ratio_formula(np.asarray(mach, dtype=np.float64), gas)
        return float(value) - target

    if branch is MachBranch.SUBSONIC:
        lower = float(np.finfo(np.float64).tiny)
        upper = 1.0
    else:
        lower = 1.0
        upper = 2.0
        while residual(upper) < 0.0 and upper < _MAX_SUPERSONIC_BRACKET:
            upper *= 2.0
        if residual(upper) < 0.0:
            raise ValueError("area_ratio is too large for the supersonic solver")

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


def _thermal_mach_from_area_scalar(
    target: float,
    branch: MachBranch,
    total_temperature: float,
    gas: _IdealThermalGas,
    *,
    allow_extrapolation: bool,
) -> tuple[float, bool]:
    if target == 1.0:
        state = _thermal_flow_state(
            1.0,
            total_temperature,
            gas,
            allow_extrapolation=allow_extrapolation,
        )
        return 1.0, state.extrapolated

    extrapolated = False

    def residual(mach: float) -> float:
        nonlocal extrapolated
        value, outside = _thermal_area_ratio_scalar(
            mach,
            total_temperature,
            gas,
            allow_extrapolation=allow_extrapolation,
        )
        extrapolated = extrapolated or outside
        return value - target

    if branch is MachBranch.SUBSONIC:
        lower = float(np.finfo(np.float64).tiny)
        upper = 1.0
    elif not allow_extrapolation and isinstance(gas, ThermallyPerfectGas):
        lower = 1.0
        upper, outside = _thermal_mach_from_static_temperature(
            gas.temperature_range[0],
            total_temperature,
            gas,
            allow_extrapolation=False,
        )
        extrapolated = extrapolated or outside
        if upper <= 1.0 or residual(upper) < 0.0:
            raise ModelRangeError(
                "supersonic area-ratio solution requires temperature below "
                f"the fitted range minimum {gas.temperature_range[0]:g} K"
            )
    else:
        lower = 1.0
        upper = 2.0
        while residual(upper) < 0.0 and upper < _MAX_SUPERSONIC_BRACKET:
            upper *= 2.0
        if residual(upper) < 0.0:
            raise ValueError("area_ratio is too large for the supersonic solver")
    result = float(
        brentq(
            residual,
            lower,
            upper,
            xtol=_ROOT_XTOL,
            rtol=_ROOT_RTOL,
            maxiter=_ROOT_MAXITER,
        )
    )
    return result, extrapolated


def _real_mach_from_area_scalar(
    target: float,
    branch: MachBranch,
    total_temperature: float,
    total_pressure: float,
    gas: BeattieBridgemanGas,
) -> float:
    if target == 1.0:
        return 1.0

    def residual(mach: float) -> float:
        return (
            _real_area_ratio_scalar(mach, total_temperature, total_pressure, gas)
            - target
        )

    if branch is MachBranch.SUBSONIC:
        lower = float(np.finfo(np.float64).tiny)
        upper = 1.0
    else:
        lower = 1.0
        upper = 2.0
        while residual(upper) < 0.0 and upper < _MAX_SUPERSONIC_BRACKET:
            upper *= 2.0
        if residual(upper) < 0.0:
            raise ValueError("area_ratio is too large for the supersonic solver")
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


def mach_from_area_ratio(
    ratio: ArrayLike,
    branch: MachBranch,
    gas: IsentropicGasModel = AIR,
    *,
    total_temperature: ArrayLike | None = None,
    total_pressure: ArrayLike | None = None,
    allow_extrapolation: bool = True,
) -> FloatResult:
    """Invert ``A/A*`` on an explicitly selected Mach branch."""
    if not isinstance(branch, MachBranch):
        raise ValueError("branch must be a MachBranch value")
    values, scalar = _validate_ratio(ratio, name="area_ratio")
    if isinstance(gas, PerfectGas):
        result = np.empty_like(values)
        for index, target in np.ndenumerate(values):
            result[index] = _mach_from_area_scalar(float(target), branch, gas)
        return return_float(result, scalar=scalar)

    if isinstance(gas, BeattieBridgemanGas):
        values, temperatures, pressures, scalar = _broadcast_real_inputs(
            values,
            scalar,
            total_temperature,
            total_pressure,
            value_name="area_ratio",
        )
        total_outside = _check_real_total_conditions(
            gas,
            temperatures,
            pressures,
            allow_extrapolation=allow_extrapolation,
        )
        result = np.empty_like(values)
        for index in np.ndindex(values.shape):
            result[index] = _real_mach_from_area_scalar(
                float(values[index]),
                branch,
                float(temperatures[index]),
                float(pressures[index]),
                gas,
            )
        _check_real_area_states(
            gas,
            total_outside,
            result,
            temperatures,
            pressures,
            allow_extrapolation=allow_extrapolation,
        )
        return return_float(result, scalar=scalar)

    assert isinstance(gas, (ThermallyPerfectGas, HarmonicOscillatorGas))
    values, temperatures, scalar = _broadcast_thermal_inputs(
        values, scalar, total_temperature, value_name="area_ratio"
    )
    if isinstance(gas, HarmonicOscillatorGas):
        _check_harmonic_total_temperature(
            gas, temperatures, allow_extrapolation=allow_extrapolation
        )
    result = np.empty_like(values)
    extrapolated = False
    for index in np.ndindex(values.shape):
        result[index], outside = _thermal_mach_from_area_scalar(
            float(values[index]),
            branch,
            float(temperatures[index]),
            gas,
            allow_extrapolation=allow_extrapolation,
        )
        extrapolated = extrapolated or outside
    _warn_if_extrapolated(gas, extrapolated)
    return return_float(result, scalar=scalar)


def critical_ratios(
    gas: IsentropicGasModel = AIR,
    *,
    total_temperature: ArrayLike | None = None,
    total_pressure: ArrayLike | None = None,
    allow_extrapolation: bool = True,
) -> CriticalRatios:
    """Return total-to-critical state ratios at Mach one."""
    if isinstance(gas, PerfectGas):
        gamma = gas.heat_capacity_ratio
        temperature_ratio = 0.5 * (gamma + 1.0)
        return CriticalRatios(
            total_temperature_ratio=temperature_ratio,
            total_pressure_ratio=temperature_ratio ** (gamma / (gamma - 1.0)),
            total_density_ratio=temperature_ratio ** (1.0 / (gamma - 1.0)),
        )

    placeholder = np.asarray(1.0, dtype=np.float64)
    if isinstance(gas, BeattieBridgemanGas):
        mach, temperatures, pressures, scalar = _broadcast_real_inputs(
            placeholder,
            True,
            total_temperature,
            total_pressure,
            value_name="critical state",
        )
        total_outside = _check_real_total_conditions(
            gas,
            temperatures,
            pressures,
            allow_extrapolation=allow_extrapolation,
        )
        (
            real_temperature_ratio,
            real_pressure_ratio,
            real_density_ratio,
            _,
            static_temperature,
            static_pressure,
        ) = _real_flow_states(np.ones_like(mach), temperatures, pressures, gas)
        _check_real_static_conditions(
            gas,
            total_outside,
            static_temperature,
            static_pressure,
            allow_extrapolation=allow_extrapolation,
        )
        return CriticalRatios(
            total_temperature_ratio=return_float(real_temperature_ratio, scalar=scalar),
            total_pressure_ratio=return_float(real_pressure_ratio, scalar=scalar),
            total_density_ratio=return_float(real_density_ratio, scalar=scalar),
        )

    assert isinstance(gas, (ThermallyPerfectGas, HarmonicOscillatorGas))
    _, temperatures, scalar = _broadcast_thermal_inputs(
        placeholder, True, total_temperature, value_name="critical state"
    )
    if isinstance(gas, HarmonicOscillatorGas):
        _check_harmonic_total_temperature(
            gas, temperatures, allow_extrapolation=allow_extrapolation
        )
    mach = np.ones_like(temperatures)
    thermal_temperature_ratio, pressure_ratio, density_ratio, _, extrapolated = (
        _thermal_flow_states(
            mach,
            temperatures,
            gas,
            allow_extrapolation=allow_extrapolation,
        )
    )
    _warn_if_extrapolated(gas, extrapolated)
    return CriticalRatios(
        total_temperature_ratio=return_float(thermal_temperature_ratio, scalar=scalar),
        total_pressure_ratio=return_float(pressure_ratio, scalar=scalar),
        total_density_ratio=return_float(density_ratio, scalar=scalar),
    )


def mass_flow_parameter(
    mach: ArrayLike,
    gas: IsentropicGasModel = AIR,
    *,
    total_temperature: ArrayLike | None = None,
    total_pressure: ArrayLike | None = None,
    allow_extrapolation: bool = True,
) -> FloatResult:
    """Return the dimensionless mass-flow parameter.

    The normalization is chosen so that mass flux equals
    ``p0 / sqrt(R T0)`` multiplied by the returned value.
    """
    values, scalar = _validate_mach(mach)
    if isinstance(gas, PerfectGas):
        gamma = gas.heat_capacity_ratio
        exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
        result = np.sqrt(gamma) * values / _temperature_factor(values, gas) ** exponent
        return return_float(result, scalar=scalar)

    if isinstance(gas, BeattieBridgemanGas):
        values, temperatures, pressures, scalar = _broadcast_real_inputs(
            values,
            scalar,
            total_temperature,
            total_pressure,
            value_name="mach",
        )
        total_outside = _check_real_total_conditions(
            gas,
            temperatures,
            pressures,
            allow_extrapolation=allow_extrapolation,
        )
        _, _, _, result, static_temperature, static_pressure = _real_flow_states(
            values, temperatures, pressures, gas
        )
        _check_real_static_conditions(
            gas,
            total_outside,
            static_temperature,
            static_pressure,
            allow_extrapolation=allow_extrapolation,
        )
        return return_float(result, scalar=scalar)

    assert isinstance(gas, (ThermallyPerfectGas, HarmonicOscillatorGas))
    values, temperatures, scalar = _broadcast_thermal_inputs(
        values, scalar, total_temperature, value_name="mach"
    )
    if isinstance(gas, HarmonicOscillatorGas):
        _check_harmonic_total_temperature(
            gas, temperatures, allow_extrapolation=allow_extrapolation
        )
    _, _, _, result, extrapolated = _thermal_flow_states(
        values,
        temperatures,
        gas,
        allow_extrapolation=allow_extrapolation,
    )
    _warn_if_extrapolated(gas, extrapolated)
    return return_float(result, scalar=scalar)


def mass_flux(
    total_pressure: ArrayLike,
    total_temperature: ArrayLike,
    mach: ArrayLike,
    gas: IsentropicGasModel = AIR,
    *,
    allow_extrapolation: bool = True,
) -> FloatResult:
    """Return isentropic mass flux in kg/(m² s)."""
    pressure, pressure_scalar = as_float_array(total_pressure, name="total_pressure")
    temperature, temperature_scalar = as_float_array(
        total_temperature, name="total_temperature"
    )
    mach_values, mach_scalar = _validate_mach(mach)
    try:
        pressure, temperature, mach_values = np.broadcast_arrays(
            pressure, temperature, mach_values
        )
    except ValueError as error:
        raise ValueError(
            "total_pressure, total_temperature, and mach must be broadcastable"
        ) from error
    if np.any(pressure <= 0.0):
        raise ValueError("total_pressure must be greater than zero")
    if np.any(temperature <= 0.0):
        raise ValueError("total_temperature must be greater than zero")

    if isinstance(gas, BeattieBridgemanGas):
        total_outside = _check_real_total_conditions(
            gas,
            np.asarray(temperature, dtype=np.float64),
            np.asarray(pressure, dtype=np.float64),
            allow_extrapolation=allow_extrapolation,
        )
        _, _, _, parameter, static_temperature, static_pressure = _real_flow_states(
            np.asarray(mach_values, dtype=np.float64),
            np.asarray(temperature, dtype=np.float64),
            np.asarray(pressure, dtype=np.float64),
            gas,
        )
        _check_real_static_conditions(
            gas,
            total_outside,
            np.asarray(static_temperature, dtype=np.float64),
            np.asarray(static_pressure, dtype=np.float64),
            allow_extrapolation=allow_extrapolation,
        )
    elif isinstance(gas, (ThermallyPerfectGas, HarmonicOscillatorGas)):
        assert isinstance(gas, (ThermallyPerfectGas, HarmonicOscillatorGas))
        if isinstance(gas, HarmonicOscillatorGas):
            _check_harmonic_total_temperature(
                gas,
                np.asarray(temperature, dtype=np.float64),
                allow_extrapolation=allow_extrapolation,
            )
        _, _, _, parameter, extrapolated = _thermal_flow_states(
            np.asarray(mach_values, dtype=np.float64),
            np.asarray(temperature, dtype=np.float64),
            gas,
            allow_extrapolation=allow_extrapolation,
        )
        _warn_if_extrapolated(gas, extrapolated)
    else:
        gamma = gas.heat_capacity_ratio
        exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
        parameter = (
            np.sqrt(gamma)
            * mach_values
            / _temperature_factor(mach_values, gas) ** exponent
        )
    result = pressure * parameter / np.sqrt(gas.specific_gas_constant * temperature)
    scalar = pressure_scalar and temperature_scalar and mach_scalar
    return return_float(np.asarray(result, dtype=np.float64), scalar=scalar)


def choked_mass_flux(
    total_pressure: ArrayLike,
    total_temperature: ArrayLike,
    gas: IsentropicGasModel = AIR,
    *,
    allow_extrapolation: bool = True,
) -> FloatResult:
    """Return the maximum isentropic mass flux at Mach one."""
    return mass_flux(
        total_pressure,
        total_temperature,
        1.0,
        gas,
        allow_extrapolation=allow_extrapolation,
    )


def isentropic_state(
    mach: ArrayLike,
    gas: IsentropicGasModel = AIR,
    *,
    total_temperature: ArrayLike,
    total_pressure: ArrayLike,
    allow_extrapolation: bool = True,
) -> IsentropicFlowState:
    """Return absolute total and static properties for isentropic flow."""
    mach_values, mach_scalar = _validate_mach(mach)
    temperatures, temperature_scalar = as_float_array(
        total_temperature, name="total_temperature"
    )
    pressures, pressure_scalar = as_float_array(total_pressure, name="total_pressure")
    if np.any(temperatures <= 0.0):
        raise ValueError("total_temperature must be greater than zero")
    if np.any(pressures <= 0.0):
        raise ValueError("total_pressure must be greater than zero")
    try:
        mach_values, temperatures, pressures = np.broadcast_arrays(
            mach_values, temperatures, pressures
        )
    except ValueError as error:
        raise ValueError(
            "mach, total_temperature, and total_pressure must be broadcastable"
        ) from error
    mach_values = np.asarray(mach_values, dtype=np.float64)
    temperatures = np.asarray(temperatures, dtype=np.float64)
    pressures = np.asarray(pressures, dtype=np.float64)
    scalar = mach_scalar and temperature_scalar and pressure_scalar

    total_density = np.empty_like(mach_values)
    static_temperature = np.empty_like(mach_values)
    static_pressure = np.empty_like(mach_values)
    static_density = np.empty_like(mach_values)
    velocity = np.empty_like(mach_values)
    sound_speed = np.empty_like(mach_values)

    if isinstance(gas, BeattieBridgemanGas):
        total_outside = _check_real_total_conditions(
            gas,
            temperatures,
            pressures,
            allow_extrapolation=allow_extrapolation,
        )
        for index in np.ndindex(mach_values.shape):
            real_state = _real_flow_state(
                float(mach_values[index]),
                float(temperatures[index]),
                float(pressures[index]),
                gas,
            )
            total_density[index] = real_state.total.density
            static_temperature[index] = real_state.static.temperature
            static_pressure[index] = real_state.static.pressure
            static_density[index] = real_state.static.density
            velocity[index] = real_state.velocity
            sound_speed[index] = real_state.static.speed_of_sound
        _check_real_static_conditions(
            gas,
            total_outside,
            static_temperature,
            static_pressure,
            allow_extrapolation=allow_extrapolation,
        )
    elif isinstance(gas, (ThermallyPerfectGas, HarmonicOscillatorGas)):
        assert isinstance(gas, (ThermallyPerfectGas, HarmonicOscillatorGas))
        if isinstance(gas, HarmonicOscillatorGas):
            _check_harmonic_total_temperature(
                gas, temperatures, allow_extrapolation=allow_extrapolation
            )
        extrapolated = False
        for index in np.ndindex(mach_values.shape):
            thermal_state = _thermal_flow_state(
                float(mach_values[index]),
                float(temperatures[index]),
                gas,
                allow_extrapolation=allow_extrapolation,
            )
            properties = _thermal_properties(
                thermal_state.static_temperature,
                gas,
                allow_extrapolation=allow_extrapolation,
            )
            static_temperature[index] = thermal_state.static_temperature
            static_pressure[index] = (
                pressures[index] / thermal_state.total_pressure_ratio
            )
            total_density[index] = pressures[index] / (
                gas.specific_gas_constant * temperatures[index]
            )
            static_density[index] = static_pressure[index] / (
                gas.specific_gas_constant * static_temperature[index]
            )
            sound_speed[index] = np.sqrt(properties.sound_speed_squared)
            velocity[index] = mach_values[index] * sound_speed[index]
            extrapolated = extrapolated or thermal_state.extrapolated
        _warn_if_extrapolated(gas, extrapolated)
    else:
        ratios = _temperature_factor(mach_values, gas)
        static_temperature = temperatures / ratios
        pressure_ratios = ratios ** (
            gas.heat_capacity_ratio / (gas.heat_capacity_ratio - 1.0)
        )
        static_pressure = pressures / pressure_ratios
        total_density = pressures / (gas.specific_gas_constant * temperatures)
        static_density = static_pressure / (
            gas.specific_gas_constant * static_temperature
        )
        sound_speed = np.sqrt(
            gas.heat_capacity_ratio * gas.specific_gas_constant * static_temperature
        )
        velocity = mach_values * sound_speed

    dynamic_pressure = 0.5 * static_density * velocity**2
    flux = static_density * velocity
    return IsentropicFlowState(
        mach=return_float(mach_values, scalar=scalar),
        total_temperature=return_float(temperatures, scalar=scalar),
        total_pressure=return_float(pressures, scalar=scalar),
        total_density=return_float(total_density, scalar=scalar),
        static_temperature=return_float(static_temperature, scalar=scalar),
        static_pressure=return_float(static_pressure, scalar=scalar),
        static_density=return_float(static_density, scalar=scalar),
        velocity=return_float(velocity, scalar=scalar),
        speed_of_sound=return_float(sound_speed, scalar=scalar),
        dynamic_pressure=return_float(dynamic_pressure, scalar=scalar),
        mass_flux=return_float(flux, scalar=scalar),
    )


__all__ = [
    "CriticalRatios",
    "IsentropicAnalysis",
    "IsentropicFlowState",
    "IsentropicGasModel",
    "IsentropicRatios",
    "MachBranch",
    "area_ratio",
    "choked_mass_flux",
    "critical_ratios",
    "isentropic_analysis",
    "isentropic_ratios",
    "isentropic_state",
    "mach_from_area_ratio",
    "mach_from_total_density_ratio",
    "mach_from_total_pressure_ratio",
    "mach_from_total_temperature_ratio",
    "mass_flow_parameter",
    "mass_flux",
]
