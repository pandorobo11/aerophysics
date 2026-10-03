"""Prandtl-Meyer expansion of calorically or thermally perfect ideal gases.

Angles are expressed in radians. Expansion state ratios are downstream over
upstream static quantities; total temperature and total pressure remain
constant.

References
----------
Ames Research Staff, *Equations, Tables, and Charts for Compressible Flow*,
NACA Report 1135, 1953.
"""

from dataclasses import dataclass
from typing import Final

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import brentq

from aerophysics._array import FloatArray, FloatResult, as_float_array, return_float
from aerophysics._thermal_expansion import expansion_state
from aerophysics.gas import AIR, PerfectGas
from aerophysics.real_gas import HarmonicOscillatorGas
from aerophysics.thermochemistry import ThermallyPerfectGas

type ExpansionGasModel = PerfectGas | ThermallyPerfectGas | HarmonicOscillatorGas

_ROOT_XTOL: Final = 1e-12
_ROOT_RTOL: Final = 4.0 * np.finfo(np.float64).eps
_ROOT_MAXITER: Final = 100
_MAX_MACH_BRACKET: Final = 1e8


@dataclass(frozen=True, slots=True)
class PrandtlMeyerExpansionResult:
    """State change through a centered Prandtl-Meyer expansion."""

    upstream_mach: FloatResult
    downstream_mach: FloatResult
    turn_angle: FloatResult
    upstream_prandtl_meyer_angle: FloatResult | None
    downstream_prandtl_meyer_angle: FloatResult | None
    static_temperature_ratio: FloatResult
    static_pressure_ratio: FloatResult
    static_density_ratio: FloatResult
    available_turn_angle: FloatResult | None = None


def _validate_supersonic_mach(mach: ArrayLike) -> tuple[FloatArray, bool]:
    values, scalar = as_float_array(mach, name="mach")
    if np.any(values < 1.0):
        raise ValueError("mach must be greater than or equal to one")
    return values, scalar


def _prandtl_meyer_array(mach: FloatArray, gas: PerfectGas) -> FloatArray:
    gamma = gas.heat_capacity_ratio
    root = np.sqrt(mach**2 - 1.0)
    coefficient = np.sqrt((gamma + 1.0) / (gamma - 1.0))
    angle = coefficient * np.arctan(
        root * np.sqrt((gamma - 1.0) / (gamma + 1.0))
    ) - np.arctan(root)
    return np.asarray(angle, dtype=np.float64)


def maximum_prandtl_meyer_angle(gas: PerfectGas = AIR) -> float:
    """Return the limiting Prandtl-Meyer angle as Mach tends to infinity."""
    gamma = gas.heat_capacity_ratio
    return float(0.5 * np.pi * (np.sqrt((gamma + 1.0) / (gamma - 1.0)) - 1.0))


def prandtl_meyer_angle(mach: ArrayLike, gas: PerfectGas = AIR) -> FloatResult:
    """Return the Prandtl-Meyer function ``nu(M)`` in radians."""
    values, scalar = _validate_supersonic_mach(mach)
    return return_float(_prandtl_meyer_array(values, gas), scalar=scalar)


def _mach_from_angle_scalar(angle: float, gas: PerfectGas) -> float:
    if angle == 0.0:
        return 1.0

    def residual(mach: float) -> float:
        value = _prandtl_meyer_array(np.asarray(mach, dtype=np.float64), gas)
        return float(value) - angle

    upper = 2.0
    while residual(upper) < 0.0 and upper < _MAX_MACH_BRACKET:
        upper *= 2.0
    if residual(upper) < 0.0:
        raise ValueError("Prandtl-Meyer angle is too close to its limiting value")
    return float(
        brentq(
            residual,
            1.0,
            upper,
            xtol=_ROOT_XTOL,
            rtol=_ROOT_RTOL,
            maxiter=_ROOT_MAXITER,
        )
    )


def mach_from_prandtl_meyer(angle: ArrayLike, gas: PerfectGas = AIR) -> FloatResult:
    """Invert the Prandtl-Meyer function for a non-negative angle."""
    values, scalar = as_float_array(angle, name="prandtl_meyer_angle")
    maximum = maximum_prandtl_meyer_angle(gas)
    if np.any(values < 0.0):
        raise ValueError("prandtl_meyer_angle must be non-negative")
    if np.any(values >= maximum):
        raise ValueError("prandtl_meyer_angle must be less than the limiting angle")
    result = np.empty_like(values)
    for index, target in np.ndenumerate(values):
        result[index] = _mach_from_angle_scalar(float(target), gas)
    return return_float(result, scalar=scalar)


def prandtl_meyer_expansion(
    upstream_mach: ArrayLike,
    turn_angle: ArrayLike,
    gas: ExpansionGasModel = AIR,
    *,
    upstream_temperature: ArrayLike | None = None,
) -> PrandtlMeyerExpansionResult:
    """Return the state after an isentropic centered expansion.

    Thermal models require upstream static temperature in K, broadcast with
    Mach and turn angle (rad). Composition is frozen and no extrapolation is
    performed. ``available_turn_angle`` is the turn to the temperature lower
    bound for thermal gases, or the exclusive vacuum limit for perfect gases.
    Absolute nu uses the sonic state at the same total enthalpy. If that state
    is outside the model range, nu is None for scalars or NaN in array entries;
    this does not discard valid downstream states.
    """
    mach, mach_scalar = _validate_supersonic_mach(upstream_mach)
    turn, turn_scalar = as_float_array(turn_angle, name="turn_angle")
    try:
        mach, turn = np.broadcast_arrays(mach, turn)
    except ValueError as error:
        raise ValueError(
            "upstream_mach and turn_angle must be broadcastable"
        ) from error
    if np.any(turn < 0.0):
        raise ValueError("turn_angle must be non-negative")

    if not isinstance(gas, PerfectGas):
        if not isinstance(gas, (ThermallyPerfectGas, HarmonicOscillatorGas)):
            raise TypeError(
                "gas must be a PerfectGas, ThermallyPerfectGas, "
                "or HarmonicOscillatorGas"
            )
        if upstream_temperature is None:
            raise ValueError(
                "upstream_temperature is required for a thermally perfect gas"
            )
        temperature, temperature_scalar = as_float_array(
            upstream_temperature, name="upstream_temperature"
        )
        if np.any(temperature <= 0.0):
            raise ValueError("upstream_temperature must be greater than zero")
        try:
            mach, turn, temperature = np.broadcast_arrays(mach, turn, temperature)
        except ValueError as error:
            raise ValueError(
                "Mach, angle, and upstream_temperature must be broadcastable"
            ) from error
        scalar = mach_scalar and turn_scalar and temperature_scalar
        values = np.empty((6, *mach.shape), dtype=np.float64)
        for index in np.ndindex(mach.shape):
            state = expansion_state(
                float(mach[index]), float(turn[index]), float(temperature[index]), gas
            )
            values[(slice(None), *index)] = (
                state.downstream_mach,
                state.temperature_ratio,
                state.pressure_ratio,
                state.density_ratio,
                state.upstream_angle,
                state.available_turn_angle,
            )

        def thermal_output(value: FloatArray) -> FloatResult:
            return return_float(value, scalar=scalar)

        def angle_output(value: FloatArray) -> FloatResult | None:
            return None if scalar and np.isnan(value) else thermal_output(value)

        return PrandtlMeyerExpansionResult(
            upstream_mach=thermal_output(mach),
            downstream_mach=thermal_output(values[0]),
            turn_angle=thermal_output(turn),
            upstream_prandtl_meyer_angle=angle_output(values[4]),
            downstream_prandtl_meyer_angle=angle_output(values[4] + turn),
            static_temperature_ratio=thermal_output(values[1]),
            static_pressure_ratio=thermal_output(values[2]),
            static_density_ratio=thermal_output(values[3]),
            available_turn_angle=thermal_output(values[5]),
        )

    upstream_angle = _prandtl_meyer_array(mach, gas)
    downstream_angle = upstream_angle + turn
    maximum = maximum_prandtl_meyer_angle(gas)
    if np.any(downstream_angle >= maximum):
        raise ValueError("expansion reaches or exceeds the limiting angle")
    downstream_mach = np.asarray(
        mach_from_prandtl_meyer(downstream_angle, gas), dtype=np.float64
    )

    gamma = gas.heat_capacity_ratio
    upstream_factor = 1.0 + 0.5 * (gamma - 1.0) * mach**2
    downstream_factor = 1.0 + 0.5 * (gamma - 1.0) * downstream_mach**2
    temperature_ratio = upstream_factor / downstream_factor
    pressure_ratio = temperature_ratio ** (gamma / (gamma - 1.0))
    density_ratio = temperature_ratio ** (1.0 / (gamma - 1.0))
    scalar = mach_scalar and turn_scalar

    def output(values: FloatArray) -> FloatResult:
        return return_float(np.asarray(values, dtype=np.float64), scalar=scalar)

    return PrandtlMeyerExpansionResult(
        upstream_mach=output(mach),
        downstream_mach=output(downstream_mach),
        turn_angle=output(turn),
        upstream_prandtl_meyer_angle=output(upstream_angle),
        downstream_prandtl_meyer_angle=output(downstream_angle),
        static_temperature_ratio=output(temperature_ratio),
        static_pressure_ratio=output(pressure_ratio),
        static_density_ratio=output(density_ratio),
        available_turn_angle=output(maximum - upstream_angle),
    )


__all__ = [
    "ExpansionGasModel",
    "PrandtlMeyerExpansionResult",
    "mach_from_prandtl_meyer",
    "maximum_prandtl_meyer_angle",
    "prandtl_meyer_angle",
    "prandtl_meyer_expansion",
]
