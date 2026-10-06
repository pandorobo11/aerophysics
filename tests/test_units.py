"""Tests for explicit SI and aviation unit conversions."""

from collections.abc import Callable

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics import units
from aerophysics._array import FloatResult


@pytest.mark.parametrize(
    ("forward", "inverse", "expected_factor"),
    [
        (units.feet_to_meters, units.meters_to_feet, 0.3048),
        (units.inches_to_meters, units.meters_to_inches, 0.0254),
        (
            units.knots_to_meters_per_second,
            units.meters_per_second_to_knots,
            pytest.approx(1852 / 3600),
        ),
        (units.psi_to_pascals, units.pascals_to_psi, pytest.approx(6894.757293168)),
        (units.psf_to_pascals, units.pascals_to_psf, pytest.approx(47.8802589803)),
        (
            units.pounds_mass_to_kilograms,
            units.kilograms_to_pounds_mass,
            0.45359237,
        ),
        (
            units.slugs_to_kilograms,
            units.kilograms_to_slugs,
            pytest.approx(14.5939029372),
        ),
        (
            units.pounds_force_to_newtons,
            units.newtons_to_pounds_force,
            pytest.approx(4.4482216152605),
        ),
        (
            units.degrees_to_radians,
            units.radians_to_degrees,
            pytest.approx(np.pi / 180.0),
        ),
    ],
)
def test_multiplicative_conversion_factors_and_round_trip(
    forward: Callable[[float], FloatResult],
    inverse: Callable[[float], FloatResult],
    expected_factor: object,
) -> None:
    assert forward(1.0) == expected_factor
    value = 123.456
    converted = forward(value)
    assert isinstance(converted, float)
    assert inverse(converted) == pytest.approx(value, rel=1e-14)


def test_temperature_reference_points_and_round_trip() -> None:
    assert units.fahrenheit_to_kelvin(-459.67) == pytest.approx(0.0, abs=1e-12)
    assert units.fahrenheit_to_kelvin(32.0) == pytest.approx(273.15)
    assert units.kelvin_to_fahrenheit(373.15) == pytest.approx(212.0)
    values = np.array([0.0, 273.15, 373.15])
    assert_allclose(
        units.fahrenheit_to_kelvin(units.kelvin_to_fahrenheit(values)),
        values,
        atol=1e-12,
    )
    assert units.celsius_to_kelvin(-273.15) == pytest.approx(0.0, abs=1e-12)
    assert units.celsius_to_kelvin(100.0) == pytest.approx(373.15)
    assert units.kelvin_to_celsius(273.15) == pytest.approx(0.0)
    assert units.rankine_to_kelvin(491.67) == pytest.approx(273.15)
    assert units.kelvin_to_rankine(373.15) == pytest.approx(671.67)
    assert_allclose(
        units.rankine_to_kelvin(units.kelvin_to_rankine(values)),
        values,
        atol=1e-12,
    )


def test_array_conversion_preserves_shape_and_float64() -> None:
    result = units.feet_to_meters([[1, 2], [3, 4]])
    assert isinstance(result, np.ndarray)
    assert result.shape == (2, 2)
    assert result.dtype == np.float64
    assert_allclose(result, [[0.3048, 0.6096], [0.9144, 1.2192]])


def test_temperature_rejects_values_below_absolute_zero() -> None:
    with pytest.raises(ValueError):
        units.fahrenheit_to_kelvin([-459.67, -500.0])
    with pytest.raises(ValueError):
        units.kelvin_to_fahrenheit(-1.0)
    with pytest.raises(ValueError):
        units.celsius_to_kelvin(-273.16)
    with pytest.raises(ValueError):
        units.rankine_to_kelvin(-1.0)
    with pytest.raises(ValueError):
        units.kelvin_to_celsius(-1.0)
    with pytest.raises(ValueError):
        units.kelvin_to_rankine(-1.0)


def test_conversions_reject_non_finite_or_non_numeric_input() -> None:
    with pytest.raises(ValueError):
        units.feet_to_meters([1.0, np.inf])
    with pytest.raises(ValueError):
        units.degrees_to_radians(["right angle"])
