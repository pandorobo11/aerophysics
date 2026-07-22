"""Tests for the U.S. Standard Atmosphere 1976 implementation."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics.atmosphere import (
    EARTH_RADIUS,
    AtmosphereState,
    geometric_to_geopotential,
    geopotential_to_geometric,
    standard_atmosphere,
)
from aerophysics.exceptions import ModelRangeError


def test_sea_level_reference_state() -> None:
    state = standard_atmosphere(0.0)
    assert isinstance(state, AtmosphereState)
    assert state.geometric_altitude == 0.0
    assert state.geopotential_altitude == 0.0
    assert state.temperature == pytest.approx(288.15)
    assert state.pressure == pytest.approx(101_325.0)
    assert state.density == pytest.approx(1.225, rel=1e-6)
    assert state.speed_of_sound == pytest.approx(340.294, rel=1e-6)
    assert state.gravity == pytest.approx(9.80665)
    assert state.dynamic_viscosity == pytest.approx(1.7894e-5)
    assert state.kinematic_viscosity == pytest.approx(1.46074e-5, rel=1e-5)
    assert state.thermal_conductivity == pytest.approx(0.025325884264264)
    assert state.prandtl_number == pytest.approx(0.70986, rel=1e-5)


@pytest.mark.parametrize(
    ("geopotential_altitude", "temperature", "pressure"),
    [
        (11_000.0, 216.65, 22_632.1),
        (20_000.0, 216.65, 5_474.89),
        (32_000.0, 228.65, 868.019),
        (47_000.0, 270.65, 110.906),
        (51_000.0, 270.65, 66.9389),
        (71_000.0, 214.65, 3.95642),
    ],
)
def test_geopotential_layer_boundary_reference_values(
    geopotential_altitude: float, temperature: float, pressure: float
) -> None:
    geometric_altitude = geopotential_to_geometric(geopotential_altitude)
    state = standard_atmosphere(geometric_altitude)
    assert state.temperature == pytest.approx(temperature, abs=1e-10)
    assert state.pressure == pytest.approx(pressure, rel=1e-5)


def test_geometric_altitude_endpoints() -> None:
    state = standard_atmosphere([-5_000.0, 86_000.0])
    assert isinstance(state.temperature, np.ndarray)
    assert_allclose(state.temperature, [320.676, 186.946], rtol=5e-6)
    assert_allclose(state.pressure, [177_762.0, 0.37338], rtol=5e-4)


def test_layer_boundaries_are_continuous() -> None:
    boundaries = np.array([11_000.0, 20_000.0, 32_000.0, 47_000.0, 51_000.0, 71_000.0])
    geometric = geopotential_to_geometric(boundaries)
    below = standard_atmosphere(geometric - 1e-5)
    above = standard_atmosphere(geometric + 1e-5)
    assert_allclose(below.temperature, above.temperature, rtol=1e-9, atol=1e-6)
    assert_allclose(below.pressure, above.pressure, rtol=1e-8)


def test_array_shape_and_derived_properties() -> None:
    state = standard_atmosphere([[0.0, 1_000.0], [10_000.0, 20_000.0]])
    for value in (
        state.geometric_altitude,
        state.geopotential_altitude,
        state.temperature,
        state.pressure,
        state.density,
        state.speed_of_sound,
        state.gravity,
        state.dynamic_viscosity,
        state.kinematic_viscosity,
        state.cp,
        state.cv,
        state.heat_capacity_ratio,
        state.specific_gas_constant,
        state.thermal_conductivity,
        state.prandtl_number,
    ):
        assert isinstance(value, np.ndarray)
        assert value.shape == (2, 2)
        assert value.dtype == np.float64
    assert_allclose(
        state.prandtl_number,
        np.asarray(state.dynamic_viscosity)
        * np.asarray(state.cp)
        / np.asarray(state.thermal_conductivity),
    )


def test_altitude_conversions_round_trip() -> None:
    geometric = np.array([-5_000.0, 0.0, 11_000.0, 86_000.0])
    geopotential = geometric_to_geopotential(geometric)
    assert isinstance(geopotential, np.ndarray)
    assert_allclose(geopotential_to_geometric(geopotential), geometric, atol=1e-10)
    assert geometric_to_geopotential(0.0) == 0.0
    assert geopotential_to_geometric(0.0) == 0.0


@pytest.mark.parametrize("altitude", [-5_000.01, 86_000.01])
def test_standard_atmosphere_rejects_out_of_range_altitude(altitude: float) -> None:
    with pytest.raises(ModelRangeError):
        standard_atmosphere(altitude)


def test_standard_atmosphere_rejects_invalid_array_as_a_whole() -> None:
    with pytest.raises(ModelRangeError):
        standard_atmosphere([0.0, 90_000.0])
    with pytest.raises(ValueError):
        standard_atmosphere([0.0, np.nan])


def test_altitude_conversion_rejects_singular_coordinates() -> None:
    with pytest.raises(ValueError):
        geometric_to_geopotential(-EARTH_RADIUS)
    with pytest.raises(ValueError):
        geopotential_to_geometric(EARTH_RADIUS)
