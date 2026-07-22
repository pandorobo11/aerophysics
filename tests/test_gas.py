"""Tests for perfect-gas and transport-property models."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics.gas import (
    AIR,
    AIR_CONDUCTIVITY,
    AIR_VISCOSITY,
    PerfectGas,
    SutherlandModel,
    USSAConductivityModel,
)


def test_air_heat_capacities_are_consistent() -> None:
    assert AIR.cp / AIR.cv == pytest.approx(AIR.heat_capacity_ratio)
    assert AIR.cp - AIR.cv == pytest.approx(AIR.specific_gas_constant)


def test_speed_of_sound_scalar_and_array() -> None:
    assert AIR.speed_of_sound(288.15) == pytest.approx(340.293988, rel=1e-8)
    result = AIR.speed_of_sound([[288.15, 216.65]])
    assert isinstance(result, np.ndarray)
    assert result.dtype == np.float64
    assert result.shape == (1, 2)
    assert_allclose(result, [[340.293988, 295.069494]], rtol=1e-8)


@pytest.mark.parametrize(
    ("gas_constant", "gamma"),
    [(0.0, 1.4), (-1.0, 1.4), (np.inf, 1.4), (287.0, 1.0), (287.0, np.nan)],
)
def test_perfect_gas_rejects_invalid_constants(
    gas_constant: float, gamma: float
) -> None:
    with pytest.raises(ValueError):
        PerfectGas(gas_constant, gamma)


@pytest.mark.parametrize("temperature", [0.0, -1.0, np.nan, np.inf])
def test_speed_of_sound_rejects_invalid_temperature(temperature: float) -> None:
    with pytest.raises(ValueError):
        AIR.speed_of_sound(temperature)


def test_array_conversion_rejects_non_numeric_values() -> None:
    with pytest.raises(ValueError, match="real numeric"):
        AIR.speed_of_sound(["cold"])


def test_sutherland_viscosity_reference_and_array() -> None:
    assert AIR_VISCOSITY.dynamic_viscosity(288.15) == pytest.approx(1.7894e-5)
    result = AIR_VISCOSITY.dynamic_viscosity([200.0, 300.0])
    assert isinstance(result, np.ndarray)
    assert_allclose(result, [1.328573517871709e-5, 1.846021864575035e-5])


def test_ussa_conductivity_reference_and_array() -> None:
    assert AIR_CONDUCTIVITY.thermal_conductivity(288.15) == pytest.approx(
        0.025325884264264
    )
    result = AIR_CONDUCTIVITY.thermal_conductivity([200.0, 300.0])
    assert isinstance(result, np.ndarray)
    assert_allclose(result, [0.018091538433556, 0.026252000780939])


@pytest.mark.parametrize(
    ("viscosity", "temperature", "sutherland"),
    [(-1.0, 288.15, 110.4), (1.0, -1.0, 110.4), (1.0, 288.15, np.inf)],
)
def test_sutherland_model_rejects_invalid_constants(
    viscosity: float, temperature: float, sutherland: float
) -> None:
    with pytest.raises(ValueError):
        SutherlandModel(viscosity, temperature, sutherland)


@pytest.mark.parametrize(
    ("coefficient", "additive", "exponential"),
    [(0.0, 245.4, 12.0), (1.0, -1.0, 12.0), (1.0, 245.4, np.nan)],
)
def test_conductivity_model_rejects_invalid_constants(
    coefficient: float, additive: float, exponential: float
) -> None:
    with pytest.raises(ValueError):
        USSAConductivityModel(coefficient, additive, exponential)


@pytest.mark.parametrize("temperature", [0.0, -10.0, np.nan, np.inf])
def test_transport_models_reject_invalid_temperature(temperature: float) -> None:
    with pytest.raises(ValueError):
        AIR_VISCOSITY.dynamic_viscosity(temperature)
    with pytest.raises(ValueError):
        AIR_CONDUCTIVITY.thermal_conductivity(temperature)
