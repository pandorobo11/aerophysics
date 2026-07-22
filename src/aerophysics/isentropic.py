"""Isentropic relations for a calorically perfect gas.

State ratios use the total-to-static convention: ``T0/T``, ``p0/p``, and
``rho0/rho``. All inputs and outputs are dimensionless except mass flux,
pressure, and temperature arguments explicitly documented in SI units.

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
from scipy.optimize import brentq

from aerophysics._array import FloatArray, FloatResult, as_float_array, return_float
from aerophysics.gas import AIR, PerfectGas

_ROOT_XTOL: Final = 1e-12
_ROOT_RTOL: Final = 4.0 * np.finfo(np.float64).eps
_ROOT_MAXITER: Final = 100
_MAX_SUPERSONIC_BRACKET: Final = 1e6


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

    total_temperature_ratio: float
    total_pressure_ratio: float
    total_density_ratio: float


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


def isentropic_ratios(mach: ArrayLike, gas: PerfectGas = AIR) -> IsentropicRatios:
    """Return total-to-static temperature, pressure, and density ratios."""
    values, scalar = _validate_mach(mach)
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


def mach_from_total_temperature_ratio(
    ratio: ArrayLike, gas: PerfectGas = AIR
) -> FloatResult:
    """Return Mach number from ``T0/T``."""
    values, scalar = _validate_ratio(ratio, name="total_temperature_ratio")
    result = np.sqrt(2.0 * (values - 1.0) / (gas.heat_capacity_ratio - 1.0))
    return return_float(result, scalar=scalar)


def mach_from_total_pressure_ratio(
    ratio: ArrayLike, gas: PerfectGas = AIR
) -> FloatResult:
    """Return Mach number from ``p0/p``."""
    values, scalar = _validate_ratio(ratio, name="total_pressure_ratio")
    gamma = gas.heat_capacity_ratio
    temperature_ratio = values ** ((gamma - 1.0) / gamma)
    result = np.sqrt(2.0 * (temperature_ratio - 1.0) / (gamma - 1.0))
    return return_float(result, scalar=scalar)


def mach_from_total_density_ratio(
    ratio: ArrayLike, gas: PerfectGas = AIR
) -> FloatResult:
    """Return Mach number from ``rho0/rho``."""
    values, scalar = _validate_ratio(ratio, name="total_density_ratio")
    gamma = gas.heat_capacity_ratio
    temperature_ratio = values ** (gamma - 1.0)
    result = np.sqrt(2.0 * (temperature_ratio - 1.0) / (gamma - 1.0))
    return return_float(result, scalar=scalar)


def _area_ratio_formula(mach: FloatArray, gas: PerfectGas) -> FloatArray:
    gamma = gas.heat_capacity_ratio
    factor = (2.0 / (gamma + 1.0)) * _temperature_factor(mach, gas)
    exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    return factor**exponent / mach


def area_ratio(mach: ArrayLike, gas: PerfectGas = AIR) -> FloatResult:
    """Return the isentropic area ratio ``A/A*``."""
    values, scalar = _validate_mach(mach, positive=True)
    return return_float(_area_ratio_formula(values, gas), scalar=scalar)


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


def mach_from_area_ratio(
    ratio: ArrayLike,
    branch: MachBranch,
    gas: PerfectGas = AIR,
) -> FloatResult:
    """Invert ``A/A*`` on an explicitly selected Mach branch."""
    if not isinstance(branch, MachBranch):
        raise ValueError("branch must be a MachBranch value")
    values, scalar = _validate_ratio(ratio, name="area_ratio")
    result = np.empty_like(values)
    for index, target in np.ndenumerate(values):
        result[index] = _mach_from_area_scalar(float(target), branch, gas)
    return return_float(result, scalar=scalar)


def critical_ratios(gas: PerfectGas = AIR) -> CriticalRatios:
    """Return total-to-critical state ratios at Mach one."""
    gamma = gas.heat_capacity_ratio
    temperature_ratio = 0.5 * (gamma + 1.0)
    return CriticalRatios(
        total_temperature_ratio=temperature_ratio,
        total_pressure_ratio=temperature_ratio ** (gamma / (gamma - 1.0)),
        total_density_ratio=temperature_ratio ** (1.0 / (gamma - 1.0)),
    )


def mass_flow_parameter(mach: ArrayLike, gas: PerfectGas = AIR) -> FloatResult:
    """Return the dimensionless perfect-gas mass-flow parameter.

    The normalization is chosen so that mass flux equals
    ``p0 / sqrt(R T0)`` multiplied by the returned value.
    """
    values, scalar = _validate_mach(mach)
    gamma = gas.heat_capacity_ratio
    exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    result = np.sqrt(gamma) * values / _temperature_factor(values, gas) ** exponent
    return return_float(result, scalar=scalar)


def mass_flux(
    total_pressure: ArrayLike,
    total_temperature: ArrayLike,
    mach: ArrayLike,
    gas: PerfectGas = AIR,
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

    gamma = gas.heat_capacity_ratio
    exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    parameter = (
        np.sqrt(gamma) * mach_values / _temperature_factor(mach_values, gas) ** exponent
    )
    result = pressure * parameter / np.sqrt(gas.specific_gas_constant * temperature)
    scalar = pressure_scalar and temperature_scalar and mach_scalar
    return return_float(np.asarray(result, dtype=np.float64), scalar=scalar)


def choked_mass_flux(
    total_pressure: ArrayLike,
    total_temperature: ArrayLike,
    gas: PerfectGas = AIR,
) -> FloatResult:
    """Return the maximum isentropic mass flux at Mach one."""
    return mass_flux(total_pressure, total_temperature, 1.0, gas)


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
