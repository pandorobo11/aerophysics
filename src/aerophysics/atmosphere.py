"""U.S. Standard Atmosphere 1976 below 86 km geometric altitude.

The public altitude coordinate is geometric altitude in metres. Layer
hydrostatics use geopotential altitude as specified by the standard. Returned
properties are SI values and use the standard's source-specific constants.

References
----------
U.S. Standard Atmosphere, 1976, NOAA-S/T 76-1562, NASA-TM-X-74335.
"""

from dataclasses import dataclass
from typing import Final, overload

import numpy as np
from numpy.typing import ArrayLike

from aerophysics._array import FloatArray, FloatResult, as_float_array, return_float
from aerophysics.exceptions import ModelRangeError
from aerophysics.gas import AIR
from aerophysics.transport import AIR_CONDUCTIVITY, AIR_VISCOSITY

EARTH_RADIUS: Final = 6_356_766.0
"""Effective Earth radius in metres used by U.S. Standard Atmosphere 1976."""

STANDARD_GRAVITY: Final = 9.80665
"""Standard sea-level acceleration due to gravity in m/s²."""

MIN_GEOMETRIC_ALTITUDE: Final = -5_000.0
MAX_GEOMETRIC_ALTITUDE: Final = 86_000.0

_BASE_GEOPOTENTIAL_ALTITUDES = np.array(
    [0.0, 11_000.0, 20_000.0, 32_000.0, 47_000.0, 51_000.0, 71_000.0],
    dtype=np.float64,
)
_BASE_TEMPERATURES = np.array(
    [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65],
    dtype=np.float64,
)
_LAPSE_RATES = np.array(
    [-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002],
    dtype=np.float64,
)


def _base_pressures() -> FloatArray:
    pressures = np.empty_like(_BASE_GEOPOTENTIAL_ALTITUDES)
    pressures[0] = 101_325.0
    gas_constant = AIR.specific_gas_constant
    for layer in range(1, pressures.size):
        previous = layer - 1
        delta_height = (
            _BASE_GEOPOTENTIAL_ALTITUDES[layer] - _BASE_GEOPOTENTIAL_ALTITUDES[previous]
        )
        base_temperature = _BASE_TEMPERATURES[previous]
        lapse_rate = _LAPSE_RATES[previous]
        if lapse_rate == 0.0:
            factor = np.exp(
                -STANDARD_GRAVITY * delta_height / (gas_constant * base_temperature)
            )
        else:
            top_temperature = base_temperature + lapse_rate * delta_height
            factor = (base_temperature / top_temperature) ** (
                STANDARD_GRAVITY / (gas_constant * lapse_rate)
            )
        pressures[layer] = pressures[previous] * factor
    return pressures


_BASE_PRESSURES = _base_pressures()


@dataclass(frozen=True, slots=True)
class AtmosphereState:
    """Properties of the U.S. Standard Atmosphere at one or more altitudes."""

    geometric_altitude: FloatResult
    geopotential_altitude: FloatResult
    temperature: FloatResult
    pressure: FloatResult
    density: FloatResult
    speed_of_sound: FloatResult
    gravity: FloatResult
    dynamic_viscosity: FloatResult
    kinematic_viscosity: FloatResult
    cp: FloatResult
    cv: FloatResult
    heat_capacity_ratio: FloatResult
    specific_gas_constant: FloatResult
    thermal_conductivity: FloatResult
    prandtl_number: FloatResult


@overload
def geometric_to_geopotential(geometric_altitude: float) -> float: ...


@overload
def geometric_to_geopotential(geometric_altitude: ArrayLike) -> FloatArray: ...


def geometric_to_geopotential(geometric_altitude: ArrayLike) -> FloatResult:
    """Convert geometric altitude in metres to geopotential altitude."""
    altitude, scalar = as_float_array(geometric_altitude, name="geometric_altitude")
    if np.any(altitude <= -EARTH_RADIUS):
        raise ValueError("geometric_altitude must be greater than -EARTH_RADIUS")
    result = EARTH_RADIUS * altitude / (EARTH_RADIUS + altitude)
    return return_float(result, scalar=scalar)


@overload
def geopotential_to_geometric(geopotential_altitude: float) -> float: ...


@overload
def geopotential_to_geometric(geopotential_altitude: ArrayLike) -> FloatArray: ...


def geopotential_to_geometric(geopotential_altitude: ArrayLike) -> FloatResult:
    """Convert geopotential altitude in metres to geometric altitude."""
    altitude, scalar = as_float_array(
        geopotential_altitude, name="geopotential_altitude"
    )
    if np.any(altitude >= EARTH_RADIUS):
        raise ValueError("geopotential_altitude must be less than EARTH_RADIUS")
    result = EARTH_RADIUS * altitude / (EARTH_RADIUS - altitude)
    return return_float(result, scalar=scalar)


def _layer_state(geopotential_altitude: FloatArray) -> tuple[FloatArray, FloatArray]:
    layer = np.searchsorted(
        _BASE_GEOPOTENTIAL_ALTITUDES[1:], geopotential_altitude, side="right"
    )
    base_height = _BASE_GEOPOTENTIAL_ALTITUDES[layer]
    base_temperature = _BASE_TEMPERATURES[layer]
    base_pressure = _BASE_PRESSURES[layer]
    lapse_rate = _LAPSE_RATES[layer]
    delta_height = geopotential_altitude - base_height
    temperature = base_temperature + lapse_rate * delta_height

    gas_constant = AIR.specific_gas_constant
    pressure = np.empty_like(temperature)
    isothermal = lapse_rate == 0.0
    pressure[isothermal] = base_pressure[isothermal] * np.exp(
        -STANDARD_GRAVITY
        * delta_height[isothermal]
        / (gas_constant * base_temperature[isothermal])
    )
    gradient = ~isothermal
    pressure[gradient] = base_pressure[gradient] * (
        base_temperature[gradient] / temperature[gradient]
    ) ** (STANDARD_GRAVITY / (gas_constant * lapse_rate[gradient]))
    return temperature, pressure


def standard_atmosphere(geometric_altitude: ArrayLike) -> AtmosphereState:
    """Evaluate U.S. Standard Atmosphere 1976 from -5 to 86 km.

    Parameters
    ----------
    geometric_altitude:
        Geometric altitude in metres.

    Returns
    -------
    AtmosphereState
        Standard-atmosphere properties. Scalar input produces scalar fields;
        array-like input produces float64 arrays of matching shape.

    Raises
    ------
    ModelRangeError
        If any altitude lies outside the implemented geometric range.
    ValueError
        If an altitude is non-numeric or non-finite.
    """
    altitude, scalar = as_float_array(geometric_altitude, name="geometric_altitude")
    if np.any(
        (altitude < MIN_GEOMETRIC_ALTITUDE) | (altitude > MAX_GEOMETRIC_ALTITUDE)
    ):
        raise ModelRangeError(
            "geometric_altitude must be between -5000 and 86000 metres"
        )

    geopotential = np.asarray(geometric_to_geopotential(altitude), dtype=np.float64)
    temperature, pressure = _layer_state(geopotential)
    density = pressure / (AIR.specific_gas_constant * temperature)
    speed_of_sound = np.sqrt(
        AIR.heat_capacity_ratio * AIR.specific_gas_constant * temperature
    )
    gravity = STANDARD_GRAVITY * (EARTH_RADIUS / (EARTH_RADIUS + altitude)) ** 2
    dynamic_viscosity = np.asarray(
        AIR_VISCOSITY.dynamic_viscosity(temperature), dtype=np.float64
    )
    thermal_conductivity = np.asarray(
        AIR_CONDUCTIVITY.thermal_conductivity(temperature), dtype=np.float64
    )
    kinematic_viscosity = dynamic_viscosity / density
    prandtl_number = dynamic_viscosity * AIR.cp / thermal_conductivity

    def output(values: FloatArray) -> FloatResult:
        return return_float(values, scalar=scalar)

    shape = temperature.shape
    cp = np.full(shape, AIR.cp, dtype=np.float64)
    cv = np.full(shape, AIR.cv, dtype=np.float64)
    gamma = np.full(shape, AIR.heat_capacity_ratio, dtype=np.float64)
    gas_constant = np.full(shape, AIR.specific_gas_constant, dtype=np.float64)
    return AtmosphereState(
        geometric_altitude=output(altitude),
        geopotential_altitude=output(geopotential),
        temperature=output(temperature),
        pressure=output(pressure),
        density=output(density),
        speed_of_sound=output(speed_of_sound),
        gravity=output(gravity),
        dynamic_viscosity=output(dynamic_viscosity),
        kinematic_viscosity=output(kinematic_viscosity),
        cp=output(cp),
        cv=output(cv),
        heat_capacity_ratio=output(gamma),
        specific_gas_constant=output(gas_constant),
        thermal_conductivity=output(thermal_conductivity),
        prandtl_number=output(prandtl_number),
    )


__all__ = [
    "EARTH_RADIUS",
    "MAX_GEOMETRIC_ALTITUDE",
    "MIN_GEOMETRIC_ALTITUDE",
    "STANDARD_GRAVITY",
    "AtmosphereState",
    "geometric_to_geopotential",
    "geopotential_to_geometric",
    "standard_atmosphere",
]
