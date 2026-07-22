"""Tests for the integrated FlightCondition API."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics import FlightCondition
from aerophysics.exceptions import ModelRangeError


def test_scalar_condition_from_mach() -> None:
    condition = FlightCondition.from_mach(0.0, 0.8, characteristic_length=1.5)
    assert condition.mach == 0.8
    assert condition.velocity == pytest.approx(272.235, rel=1e-5)
    assert condition.dynamic_pressure == pytest.approx(45_393.6, rel=1e-5)
    assert condition.reynolds_number_per_length == pytest.approx(18_637_500.0, rel=1e-3)
    assert condition.characteristic_length == 1.5
    assert condition.reynolds_number == pytest.approx(
        1.5 * float(condition.reynolds_number_per_length)
    )
    assert condition.total_temperature == pytest.approx(325.0332)
    assert condition.total_pressure == pytest.approx(154_453.75, rel=1e-5)
    assert condition.total_density == pytest.approx(1.6554213, rel=1e-7)


def test_velocity_and_mach_constructors_are_equivalent() -> None:
    from_mach = FlightCondition.from_mach(12_000.0, 2.2, 0.75)
    from_velocity = FlightCondition.from_velocity(
        12_000.0, from_mach.velocity, characteristic_length=0.75
    )
    for field in (
        "mach",
        "velocity",
        "dynamic_pressure",
        "reynolds_number_per_length",
        "reynolds_number",
        "total_temperature",
        "total_pressure",
        "total_density",
    ):
        assert getattr(from_velocity, field) == pytest.approx(
            getattr(from_mach, field), rel=1e-13
        )


def test_characteristic_length_is_optional() -> None:
    condition = FlightCondition.from_mach(1_000.0, 0.3)
    assert condition.characteristic_length is None
    assert condition.reynolds_number is None
    assert float(condition.reynolds_number_per_length) > 0.0


def test_zero_velocity_condition() -> None:
    condition = FlightCondition.from_velocity(0.0, 0.0)
    assert condition.mach == 0.0
    assert condition.dynamic_pressure == 0.0
    assert condition.reynolds_number_per_length == 0.0
    assert condition.total_temperature == condition.atmosphere.temperature
    assert condition.total_pressure == condition.atmosphere.pressure
    assert condition.total_density == condition.atmosphere.density


def test_array_inputs_broadcast_and_return_float64() -> None:
    condition = FlightCondition.from_mach(
        [[0.0], [10_000.0]], [0.5, 1.0, 2.0], characteristic_length=2.0
    )
    for value in (
        condition.mach,
        condition.velocity,
        condition.dynamic_pressure,
        condition.reynolds_number_per_length,
        condition.characteristic_length,
        condition.reynolds_number,
        condition.total_temperature,
        condition.total_pressure,
        condition.total_density,
        condition.atmosphere.temperature,
    ):
        assert isinstance(value, np.ndarray)
        assert value.shape == (2, 3)
        assert value.dtype == np.float64
    assert isinstance(condition.reynolds_number, np.ndarray)
    assert isinstance(condition.reynolds_number_per_length, np.ndarray)
    assert_allclose(
        condition.reynolds_number,
        2.0 * condition.reynolds_number_per_length,
    )


def test_flight_condition_is_frozen() -> None:
    condition = FlightCondition.from_mach(0.0, 1.0)
    with pytest.raises(FrozenInstanceError):
        condition.mach = 2.0  # type: ignore[misc]


@pytest.mark.parametrize("mach", [-0.1, np.nan, np.inf])
def test_from_mach_rejects_invalid_mach(mach: float) -> None:
    with pytest.raises(ValueError):
        FlightCondition.from_mach(0.0, mach)


@pytest.mark.parametrize("velocity", [-1.0, np.nan, np.inf])
def test_from_velocity_rejects_invalid_velocity(velocity: float) -> None:
    with pytest.raises(ValueError):
        FlightCondition.from_velocity(0.0, velocity)


@pytest.mark.parametrize("length", [0.0, -1.0, np.nan, np.inf])
def test_rejects_invalid_characteristic_length(length: float) -> None:
    with pytest.raises(ValueError):
        FlightCondition.from_mach(0.0, 1.0, length)


def test_rejects_incompatible_input_shapes() -> None:
    with pytest.raises(ValueError, match="geometric_altitude and mach"):
        FlightCondition.from_mach([0.0, 1.0], [0.5, 1.0, 2.0])
    with pytest.raises(ValueError, match="characteristic_length"):
        FlightCondition.from_velocity([0.0, 1.0], [100.0, 200.0], [1.0, 2.0, 3.0])


def test_atmosphere_range_errors_are_preserved() -> None:
    with pytest.raises(ModelRangeError):
        FlightCondition.from_mach(90_000.0, 1.0)
