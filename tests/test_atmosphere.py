"""Tests for the U.S. Standard Atmosphere 1976 implementation."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics.atmosphere import (
    EARTH_RADIUS,
    geometric_to_geopotential,
    geopotential_to_geometric,
    standard_atmosphere,
)
from aerophysics.exceptions import ModelRangeError


def test_sea_level_scalar_contract_and_conductivity() -> None:
    # Tabulated sea-level values are covered by the official-table verification.
    state = standard_atmosphere(0.0)
    assert isinstance(state.pressure, float)
    assert state.geometric_altitude == 0.0
    assert state.geopotential_altitude == 0.0
    assert state.thermal_conductivity == pytest.approx(0.025325884264264)


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


def test_altitude_conversions_round_trip() -> None:
    geometric = np.array([-5_000.0, 0.0, 11_000.0, 86_000.0])
    geopotential = geometric_to_geopotential(geometric)
    assert isinstance(geopotential, np.ndarray)
    assert_allclose(
        geopotential_to_geometric(geopotential), geometric, rtol=1e-14, atol=1e-10
    )
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
