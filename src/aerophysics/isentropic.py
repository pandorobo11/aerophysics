"""Isentropic relations for calorically and thermally perfect gases.

State ratios use the total-to-static convention: ``T0/T``, ``p0/p``, and
``rho0/rho``. All inputs and outputs are dimensionless except mass flux,
pressure, and temperature arguments explicitly documented in SI units.

Calorically perfect gases use the closed-form constant-``gamma`` relations.
Thermally perfect gases use temperature-dependent NASA-polynomial enthalpy and
entropy and therefore require total temperature.

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
from aerophysics.thermochemistry import ThermallyPerfectGas

_ROOT_XTOL: Final = 1e-12
_ROOT_RTOL: Final = 4.0 * np.finfo(np.float64).eps
_ROOT_MAXITER: Final = 100
_MAX_SUPERSONIC_BRACKET: Final = 1e6
_MAX_TEMPERATURE_BRACKET_STEPS: Final = 80

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
    gas: ThermallyPerfectGas,
    *,
    allow_extrapolation: bool,
) -> _ThermalProperties:
    value = np.asarray(temperature, dtype=np.float64)
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
    gas_constant = gas.specific_gas_constant
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
    minimum, maximum = gas.temperature_range
    return _ThermalProperties(
        enthalpy=enthalpy,
        entropy=entropy,
        heat_capacity_ratio=gamma,
        sound_speed_squared=sound_speed_squared,
        extrapolated=temperature < minimum or temperature > maximum,
    )


def _warn_if_extrapolated(gas: ThermallyPerfectGas, extrapolated: bool) -> None:
    if not extrapolated:
        return
    minimum, maximum = gas.temperature_range
    warnings.warn(
        (
            "isentropic solution uses temperature outside the fitted range "
            f"{minimum:g}--{maximum:g} K; the nearest polynomial region was "
            "extrapolated"
        ),
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
        raise ValueError(
            "total_temperature is required for a thermally perfect gas"
        )
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


def _find_lower_temperature(
    residual: Callable[[float], float],
    total_temperature: float,
    gas: ThermallyPerfectGas,
    *,
    allow_extrapolation: bool,
) -> float:
    minimum = gas.temperature_range[0]
    if not allow_extrapolation:
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
        "could not bracket a positive-temperature thermally perfect "
        "isentropic solution"
    )


def _thermal_flow_state(
    mach: float,
    total_temperature: float,
    gas: ThermallyPerfectGas,
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
        mach
        * np.sqrt(static.heat_capacity_ratio * temperature_ratio)
        / pressure_ratio
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
    gas: ThermallyPerfectGas,
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


def isentropic_ratios(
    mach: ArrayLike,
    gas: PerfectGas | ThermallyPerfectGas = AIR,
    *,
    total_temperature: ArrayLike | None = None,
    allow_extrapolation: bool = True,
) -> IsentropicRatios:
    """Return total-to-static temperature, pressure, and density ratios.

    ``total_temperature`` is required when ``gas`` is thermally perfect. NASA
    polynomial temperatures are extrapolated by default with an
    :class:`~aerophysics.exceptions.ApplicabilityWarning`; pass
    ``allow_extrapolation=False`` to enforce the fitted range.
    """
    values, scalar = _validate_mach(mach)
    if isinstance(gas, ThermallyPerfectGas):
        values, temperatures, scalar = _broadcast_thermal_inputs(
            values, scalar, total_temperature, value_name="mach"
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
    gas: ThermallyPerfectGas,
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
    mach = np.sqrt(
        2.0 * max(0.0, enthalpy_difference) / static.sound_speed_squared
    )
    if not np.isfinite(mach):
        raise ModelRangeError("thermally perfect isentropic Mach number is invalid")
    return float(mach), total.extrapolated or static.extrapolated


def mach_from_total_temperature_ratio(
    ratio: ArrayLike,
    gas: PerfectGas | ThermallyPerfectGas = AIR,
    *,
    total_temperature: ArrayLike | None = None,
    allow_extrapolation: bool = True,
) -> FloatResult:
    """Return Mach number from ``T0/T``."""
    values, scalar = _validate_ratio(ratio, name="total_temperature_ratio")
    if not isinstance(gas, ThermallyPerfectGas):
        result = np.sqrt(2.0 * (values - 1.0) / (gas.heat_capacity_ratio - 1.0))
        return return_float(result, scalar=scalar)

    values, temperatures, scalar = _broadcast_thermal_inputs(
        values, scalar, total_temperature, value_name="total_temperature_ratio"
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
    gas: ThermallyPerfectGas,
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
    gas: ThermallyPerfectGas,
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


def mach_from_total_pressure_ratio(
    ratio: ArrayLike,
    gas: PerfectGas | ThermallyPerfectGas = AIR,
    *,
    total_temperature: ArrayLike | None = None,
    allow_extrapolation: bool = True,
) -> FloatResult:
    """Return Mach number from ``p0/p``."""
    values, scalar = _validate_ratio(ratio, name="total_pressure_ratio")
    if not isinstance(gas, ThermallyPerfectGas):
        gamma = gas.heat_capacity_ratio
        temperature_ratio = values ** ((gamma - 1.0) / gamma)
        result = np.sqrt(2.0 * (temperature_ratio - 1.0) / (gamma - 1.0))
        return return_float(result, scalar=scalar)

    values, temperatures, scalar = _broadcast_thermal_inputs(
        values, scalar, total_temperature, value_name="total_pressure_ratio"
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
    gas: PerfectGas | ThermallyPerfectGas = AIR,
    *,
    total_temperature: ArrayLike | None = None,
    allow_extrapolation: bool = True,
) -> FloatResult:
    """Return Mach number from ``rho0/rho``."""
    values, scalar = _validate_ratio(ratio, name="total_density_ratio")
    if not isinstance(gas, ThermallyPerfectGas):
        gamma = gas.heat_capacity_ratio
        temperature_ratio = values ** (gamma - 1.0)
        result = np.sqrt(2.0 * (temperature_ratio - 1.0) / (gamma - 1.0))
        return return_float(result, scalar=scalar)

    values, temperatures, scalar = _broadcast_thermal_inputs(
        values, scalar, total_temperature, value_name="total_density_ratio"
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
    gas: ThermallyPerfectGas,
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


def area_ratio(
    mach: ArrayLike,
    gas: PerfectGas | ThermallyPerfectGas = AIR,
    *,
    total_temperature: ArrayLike | None = None,
    allow_extrapolation: bool = True,
) -> FloatResult:
    """Return the isentropic area ratio ``A/A*``."""
    values, scalar = _validate_mach(mach, positive=True)
    if not isinstance(gas, ThermallyPerfectGas):
        return return_float(_area_ratio_formula(values, gas), scalar=scalar)

    values, temperatures, scalar = _broadcast_thermal_inputs(
        values, scalar, total_temperature, value_name="mach"
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
    gas: ThermallyPerfectGas,
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
    elif not allow_extrapolation:
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


def mach_from_area_ratio(
    ratio: ArrayLike,
    branch: MachBranch,
    gas: PerfectGas | ThermallyPerfectGas = AIR,
    *,
    total_temperature: ArrayLike | None = None,
    allow_extrapolation: bool = True,
) -> FloatResult:
    """Invert ``A/A*`` on an explicitly selected Mach branch."""
    if not isinstance(branch, MachBranch):
        raise ValueError("branch must be a MachBranch value")
    values, scalar = _validate_ratio(ratio, name="area_ratio")
    if not isinstance(gas, ThermallyPerfectGas):
        result = np.empty_like(values)
        for index, target in np.ndenumerate(values):
            result[index] = _mach_from_area_scalar(float(target), branch, gas)
        return return_float(result, scalar=scalar)

    values, temperatures, scalar = _broadcast_thermal_inputs(
        values, scalar, total_temperature, value_name="area_ratio"
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
    gas: PerfectGas | ThermallyPerfectGas = AIR,
    *,
    total_temperature: ArrayLike | None = None,
    allow_extrapolation: bool = True,
) -> CriticalRatios:
    """Return total-to-critical state ratios at Mach one."""
    if not isinstance(gas, ThermallyPerfectGas):
        gamma = gas.heat_capacity_ratio
        temperature_ratio = 0.5 * (gamma + 1.0)
        return CriticalRatios(
            total_temperature_ratio=temperature_ratio,
            total_pressure_ratio=temperature_ratio ** (gamma / (gamma - 1.0)),
            total_density_ratio=temperature_ratio ** (1.0 / (gamma - 1.0)),
        )

    placeholder = np.asarray(1.0, dtype=np.float64)
    _, temperatures, scalar = _broadcast_thermal_inputs(
        placeholder, True, total_temperature, value_name="critical state"
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
        total_temperature_ratio=return_float(
            thermal_temperature_ratio, scalar=scalar
        ),
        total_pressure_ratio=return_float(pressure_ratio, scalar=scalar),
        total_density_ratio=return_float(density_ratio, scalar=scalar),
    )


def mass_flow_parameter(
    mach: ArrayLike,
    gas: PerfectGas | ThermallyPerfectGas = AIR,
    *,
    total_temperature: ArrayLike | None = None,
    allow_extrapolation: bool = True,
) -> FloatResult:
    """Return the dimensionless ideal-gas mass-flow parameter.

    The normalization is chosen so that mass flux equals
    ``p0 / sqrt(R T0)`` multiplied by the returned value.
    """
    values, scalar = _validate_mach(mach)
    if not isinstance(gas, ThermallyPerfectGas):
        gamma = gas.heat_capacity_ratio
        exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
        result = (
            np.sqrt(gamma)
            * values
            / _temperature_factor(values, gas) ** exponent
        )
        return return_float(result, scalar=scalar)

    values, temperatures, scalar = _broadcast_thermal_inputs(
        values, scalar, total_temperature, value_name="mach"
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
    gas: PerfectGas | ThermallyPerfectGas = AIR,
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

    if isinstance(gas, ThermallyPerfectGas):
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
    result = (
        pressure
        * parameter
        / np.sqrt(gas.specific_gas_constant * temperature)
    )
    scalar = pressure_scalar and temperature_scalar and mach_scalar
    return return_float(np.asarray(result, dtype=np.float64), scalar=scalar)


def choked_mass_flux(
    total_pressure: ArrayLike,
    total_temperature: ArrayLike,
    gas: PerfectGas | ThermallyPerfectGas = AIR,
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


__all__ = [
    "CriticalRatios",
    "IsentropicRatios",
    "MachBranch",
    "area_ratio",
    "choked_mass_flux",
    "critical_ratios",
    "isentropic_ratios",
    "mach_from_area_ratio",
    "mach_from_total_density_ratio",
    "mach_from_total_pressure_ratio",
    "mach_from_total_temperature_ratio",
    "mass_flow_parameter",
    "mass_flux",
]
