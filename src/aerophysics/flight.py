"""Integrated standard-atmosphere flight conditions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Self

import numpy as np
from numpy.typing import ArrayLike

from aerophysics._array import FloatArray, FloatResult, as_float_array, return_float
from aerophysics.atmosphere import AtmosphereState, standard_atmosphere
from aerophysics.isentropic import isentropic_ratios


def _broadcast_inputs(
    geometric_altitude: ArrayLike,
    motion: ArrayLike,
    characteristic_length: ArrayLike | None,
    *,
    motion_name: str,
) -> tuple[FloatArray, FloatArray, FloatArray | None, bool]:
    altitude, altitude_scalar = as_float_array(
        geometric_altitude, name="geometric_altitude"
    )
    motion_values, motion_scalar = as_float_array(motion, name=motion_name)
    if characteristic_length is None:
        try:
            altitude, motion_values = np.broadcast_arrays(altitude, motion_values)
        except ValueError as error:
            raise ValueError(
                f"geometric_altitude and {motion_name} must be broadcastable"
            ) from error
        return altitude, motion_values, None, altitude_scalar and motion_scalar

    length, length_scalar = as_float_array(
        characteristic_length, name="characteristic_length"
    )
    if np.any(length <= 0.0):
        raise ValueError("characteristic_length must be greater than zero")
    try:
        altitude, motion_values, length = np.broadcast_arrays(
            altitude, motion_values, length
        )
    except ValueError as error:
        raise ValueError(
            "geometric_altitude, motion input, and characteristic_length "
            "must be broadcastable"
        ) from error
    return (
        altitude,
        motion_values,
        length,
        altitude_scalar and motion_scalar and length_scalar,
    )


@dataclass(frozen=True, slots=True)
class FlightCondition:
    """Atmospheric and aerodynamic state at one or more flight conditions."""

    atmosphere: AtmosphereState
    mach: FloatResult
    velocity: FloatResult
    dynamic_pressure: FloatResult
    reynolds_number_per_length: FloatResult
    characteristic_length: FloatResult | None
    reynolds_number: FloatResult | None
    total_temperature: FloatResult
    total_pressure: FloatResult
    total_density: FloatResult

    @classmethod
    def from_mach(
        cls,
        geometric_altitude: ArrayLike,
        mach: ArrayLike,
        characteristic_length: ArrayLike | None = None,
    ) -> Self:
        """Build a condition from geometric altitude and Mach number.

        Parameters
        ----------
        geometric_altitude:
            Geometric altitude in metres.
        mach:
            Non-negative Mach number.
        characteristic_length:
            Optional characteristic length in metres. When omitted,
            ``reynolds_number`` is ``None`` while
            ``reynolds_number_per_length`` remains available.
        """
        altitude, mach_values, length, scalar = _broadcast_inputs(
            geometric_altitude,
            mach,
            characteristic_length,
            motion_name="mach",
        )
        if np.any(mach_values < 0.0):
            raise ValueError("mach must be non-negative")
        atmosphere = standard_atmosphere(altitude)
        return cls._from_state(atmosphere, mach_values, length, scalar=scalar)

    @classmethod
    def from_velocity(
        cls,
        geometric_altitude: ArrayLike,
        velocity: ArrayLike,
        characteristic_length: ArrayLike | None = None,
    ) -> Self:
        """Build a condition from geometric altitude and velocity in m/s."""
        altitude, velocity_values, length, scalar = _broadcast_inputs(
            geometric_altitude,
            velocity,
            characteristic_length,
            motion_name="velocity",
        )
        if np.any(velocity_values < 0.0):
            raise ValueError("velocity must be non-negative")
        atmosphere = standard_atmosphere(altitude)
        speed_of_sound = np.asarray(atmosphere.speed_of_sound, dtype=np.float64)
        mach_values = velocity_values / speed_of_sound
        return cls._from_state(atmosphere, mach_values, length, scalar=scalar)

    @classmethod
    def _from_state(
        cls,
        atmosphere: AtmosphereState,
        mach: FloatArray,
        characteristic_length: FloatArray | None,
        *,
        scalar: bool,
    ) -> Self:
        density = np.asarray(atmosphere.density, dtype=np.float64)
        temperature = np.asarray(atmosphere.temperature, dtype=np.float64)
        pressure = np.asarray(atmosphere.pressure, dtype=np.float64)
        speed_of_sound = np.asarray(atmosphere.speed_of_sound, dtype=np.float64)
        dynamic_viscosity = np.asarray(atmosphere.dynamic_viscosity, dtype=np.float64)
        velocity = mach * speed_of_sound
        dynamic_pressure = 0.5 * density * velocity**2
        reynolds_per_length = density * velocity / dynamic_viscosity

        ratios = isentropic_ratios(mach)
        temperature_ratio = np.asarray(ratios.total_temperature_ratio, dtype=np.float64)
        pressure_ratio = np.asarray(ratios.total_pressure_ratio, dtype=np.float64)
        density_ratio = np.asarray(ratios.total_density_ratio, dtype=np.float64)

        def output(values: FloatArray) -> FloatResult:
            return return_float(np.asarray(values, dtype=np.float64), scalar=scalar)

        length_output: FloatResult | None = None
        reynolds_number: FloatResult | None = None
        if characteristic_length is not None:
            length_output = output(characteristic_length)
            reynolds_number = output(reynolds_per_length * characteristic_length)

        return cls(
            atmosphere=atmosphere,
            mach=output(mach),
            velocity=output(velocity),
            dynamic_pressure=output(dynamic_pressure),
            reynolds_number_per_length=output(reynolds_per_length),
            characteristic_length=length_output,
            reynolds_number=reynolds_number,
            total_temperature=output(temperature * temperature_ratio),
            total_pressure=output(pressure * pressure_ratio),
            total_density=output(density * density_ratio),
        )


__all__ = ["FlightCondition"]
