"""Explicit aviation and SI unit conversions.

Functions accept finite scalar or array-like values. Scalar input returns a
Python ``float``; array-like input returns a float64 NumPy array.

References
----------
NIST Special Publication 811, Appendix B, Conversion Factors.
"""

from typing import Final

import numpy as np
from numpy.typing import ArrayLike

from aerophysics._array import FloatResult, as_float_array, return_float

FOOT_TO_METER: Final = 0.3048
KNOT_TO_METER_PER_SECOND: Final = 1_852.0 / 3_600.0
POUND_MASS_TO_KILOGRAM: Final = 0.45359237
POUND_FORCE_TO_NEWTON: Final = 4.4482216152605
PSI_TO_PASCAL: Final = POUND_FORCE_TO_NEWTON / 0.0254**2
PSF_TO_PASCAL: Final = PSI_TO_PASCAL / 144.0
SLUG_TO_KILOGRAM: Final = POUND_FORCE_TO_NEWTON / FOOT_TO_METER


def _scale(value: ArrayLike, *, factor: float, name: str) -> FloatResult:
    array, scalar = as_float_array(value, name=name)
    return return_float(array * factor, scalar=scalar)


def feet_to_meters(value: ArrayLike) -> FloatResult:
    """Convert international feet to metres."""
    return _scale(value, factor=FOOT_TO_METER, name="feet")


def meters_to_feet(value: ArrayLike) -> FloatResult:
    """Convert metres to international feet."""
    return _scale(value, factor=1.0 / FOOT_TO_METER, name="meters")


def knots_to_meters_per_second(value: ArrayLike) -> FloatResult:
    """Convert knots to metres per second."""
    return _scale(value, factor=KNOT_TO_METER_PER_SECOND, name="knots")


def meters_per_second_to_knots(value: ArrayLike) -> FloatResult:
    """Convert metres per second to knots."""
    return _scale(
        value,
        factor=1.0 / KNOT_TO_METER_PER_SECOND,
        name="meters_per_second",
    )


def fahrenheit_to_kelvin(value: ArrayLike) -> FloatResult:
    """Convert absolute degrees Fahrenheit to kelvin."""
    array, scalar = as_float_array(value, name="fahrenheit")
    if np.any(array < -459.67):
        raise ValueError("fahrenheit must not be below absolute zero")
    result = np.maximum((array - 32.0) * (5.0 / 9.0) + 273.15, 0.0)
    return return_float(result, scalar=scalar)


def kelvin_to_fahrenheit(value: ArrayLike) -> FloatResult:
    """Convert kelvin to absolute degrees Fahrenheit."""
    array, scalar = as_float_array(value, name="kelvin")
    if np.any(array < 0.0):
        raise ValueError("kelvin must not be below zero")
    result = (array - 273.15) * (9.0 / 5.0) + 32.0
    return return_float(result, scalar=scalar)


def psi_to_pascals(value: ArrayLike) -> FloatResult:
    """Convert pounds-force per square inch to pascals."""
    return _scale(value, factor=PSI_TO_PASCAL, name="psi")


def pascals_to_psi(value: ArrayLike) -> FloatResult:
    """Convert pascals to pounds-force per square inch."""
    return _scale(value, factor=1.0 / PSI_TO_PASCAL, name="pascals")


def psf_to_pascals(value: ArrayLike) -> FloatResult:
    """Convert pounds-force per square foot to pascals."""
    return _scale(value, factor=PSF_TO_PASCAL, name="psf")


def pascals_to_psf(value: ArrayLike) -> FloatResult:
    """Convert pascals to pounds-force per square foot."""
    return _scale(value, factor=1.0 / PSF_TO_PASCAL, name="pascals")


def pounds_mass_to_kilograms(value: ArrayLike) -> FloatResult:
    """Convert pounds mass to kilograms."""
    return _scale(value, factor=POUND_MASS_TO_KILOGRAM, name="pounds_mass")


def kilograms_to_pounds_mass(value: ArrayLike) -> FloatResult:
    """Convert kilograms to pounds mass."""
    return _scale(
        value,
        factor=1.0 / POUND_MASS_TO_KILOGRAM,
        name="kilograms",
    )


def slugs_to_kilograms(value: ArrayLike) -> FloatResult:
    """Convert slugs to kilograms."""
    return _scale(value, factor=SLUG_TO_KILOGRAM, name="slugs")


def kilograms_to_slugs(value: ArrayLike) -> FloatResult:
    """Convert kilograms to slugs."""
    return _scale(value, factor=1.0 / SLUG_TO_KILOGRAM, name="kilograms")


def degrees_to_radians(value: ArrayLike) -> FloatResult:
    """Convert degrees to radians."""
    return _scale(value, factor=np.pi / 180.0, name="degrees")


def radians_to_degrees(value: ArrayLike) -> FloatResult:
    """Convert radians to degrees."""
    return _scale(value, factor=180.0 / np.pi, name="radians")


__all__ = [
    "FOOT_TO_METER",
    "KNOT_TO_METER_PER_SECOND",
    "POUND_FORCE_TO_NEWTON",
    "POUND_MASS_TO_KILOGRAM",
    "PSF_TO_PASCAL",
    "PSI_TO_PASCAL",
    "SLUG_TO_KILOGRAM",
    "degrees_to_radians",
    "fahrenheit_to_kelvin",
    "feet_to_meters",
    "kelvin_to_fahrenheit",
    "kilograms_to_pounds_mass",
    "kilograms_to_slugs",
    "knots_to_meters_per_second",
    "meters_per_second_to_knots",
    "meters_to_feet",
    "pascals_to_psf",
    "pascals_to_psi",
    "pounds_mass_to_kilograms",
    "psf_to_pascals",
    "psi_to_pascals",
    "radians_to_degrees",
    "slugs_to_kilograms",
]
