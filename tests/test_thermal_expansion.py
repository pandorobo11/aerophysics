"""Variable-cp expansion conservation, independent integration, and bounds."""

from dataclasses import fields

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.integrate import solve_ivp

from aerophysics import AIR_HARMONIC_OSCILLATOR, AIR_NASA7, AIR_NASA9
from aerophysics.exceptions import ExpansionConvergenceError, ModelRangeError
from aerophysics.expansion import prandtl_meyer_angle, prandtl_meyer_expansion
from aerophysics.gas import PerfectGas
from aerophysics.real_gas import HarmonicOscillatorGas
from aerophysics.thermochemistry import (
    IdealGasSpecies,
    NASA7Polynomial,
    ThermallyPerfectGas,
)


@pytest.mark.parametrize("gamma", [1.2, 1.4, 5.0 / 3.0])
def test_constant_cp_recovers_perfect_gas(gamma: float) -> None:
    thermal = HarmonicOscillatorGas(287.0, gamma)
    actual = prandtl_meyer_expansion(
        [[1.0], [2.0], [5.0]],
        np.deg2rad([0.0, 5.0, 15.0]),
        thermal,
        upstream_temperature=500.0,
    )
    expected = prandtl_meyer_expansion(
        [[1.0], [2.0], [5.0]],
        np.deg2rad([0.0, 5.0, 15.0]),
        PerfectGas(287.0, gamma),
    )
    for field in fields(actual):
        assert_allclose(
            getattr(actual, field.name),
            getattr(expected, field.name),
            rtol=2e-10,
            atol=2e-11,
        )
    with pytest.raises(ValueError, match="limiting"):
        prandtl_meyer_expansion(
            2.0,
            float(np.asarray(actual.available_turn_angle)[1, 0]) + 1e-6,
            thermal,
            upstream_temperature=500.0,
        )


def test_variable_cp_against_independent_turn_coordinate_ode() -> None:
    # Manufactured cp/R=3.5+0.001*T has analytic h/R and s/R. Integrate
    # dT/dtheta directly, independently of the production temperature integral.
    gas = ThermallyPerfectGas(
        (
            IdealGasSpecies(
                "test",
                0.028,
                NASA7Polynomial(
                    (100.0, 30000.0), ((3.5, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0),)
                ),
            ),
        ),
        (1.0,),
    )
    t1, m1, turn = 1000.0, 3.0, float(np.deg2rad(15.0))
    cp1 = 3.5 + 0.001 * t1
    u1_squared = m1**2 * cp1 / (cp1 - 1.0) * t1
    h1 = 3.5 * t1 + 0.0005 * t1**2

    def derivative(_angle: float, temperature: np.ndarray) -> list[float]:
        t = float(temperature[0])
        cp = 3.5 + 0.001 * t
        velocity = u1_squared + 2.0 * (h1 - 3.5 * t - 0.0005 * t**2)
        mach_squared = velocity / (cp / (cp - 1.0) * t)
        return [-velocity / (cp * np.sqrt(mach_squared - 1.0))]

    reference = solve_ivp(
        derivative, (0.0, turn), [t1], method="DOP853", rtol=2e-12, atol=1e-9
    )
    assert reference.success
    t2 = reference.y[0, -1]
    expected_pressure = np.exp(3.5 * np.log(t2 / t1) + 0.001 * (t2 - t1))
    actual = prandtl_meyer_expansion(m1, turn, gas, upstream_temperature=t1)
    assert actual.static_temperature_ratio == pytest.approx(t2 / t1, rel=2e-11)
    assert actual.static_pressure_ratio == pytest.approx(expected_pressure, rel=2e-11)
    assert actual.static_density_ratio == pytest.approx(
        expected_pressure * t1 / t2, rel=2e-11
    )


@pytest.mark.parametrize("gas", [AIR_NASA7, AIR_NASA9, AIR_HARMONIC_OSCILLATOR])
def test_presets_conserve_energy_entropy_across_polynomial_boundaries(
    gas: ThermallyPerfectGas | HarmonicOscillatorGas,
) -> None:
    t1, m1 = 1200.0, 3.0
    result = prandtl_meyer_expansion(m1, np.deg2rad(20.0), gas, upstream_temperature=t1)
    t2 = t1 * float(result.static_temperature_ratio)
    assert gas.temperature_range[0] < t2 < 1000.0  # NASA7 and NASA9 region boundary.
    u1 = m1 * float(gas.speed_of_sound(t1))
    u2 = float(result.downstream_mach) * float(gas.speed_of_sound(t2))
    assert float(gas.standard_enthalpy(t1)) + 0.5 * u1**2 == pytest.approx(
        float(gas.standard_enthalpy(t2)) + 0.5 * u2**2,
        rel=2e-11,
    )
    entropy_change = float(gas.entropy(t2, 100000.0)) - float(gas.entropy(t1, 100000.0))
    assert entropy_change == pytest.approx(
        gas.specific_gas_constant * np.log(result.static_pressure_ratio), rel=2e-11
    )
    assert result.downstream_mach > m1
    assert 0.0 < result.static_pressure_ratio < result.static_temperature_ratio < 1.0


def test_valid_expansion_without_in_range_sonic_reference() -> None:
    result = prandtl_meyer_expansion(
        3.0, np.deg2rad(10.0), AIR_HARMONIC_OSCILLATOR, upstream_temperature=1000.0
    )
    # Independent DOP853 integration in theta with rtol=2e-13.
    assert result.downstream_mach == pytest.approx(3.495015031609171, rel=2e-11)
    assert 1000.0 * result.static_temperature_ratio == pytest.approx(
        815.5794702524529, rel=2e-11
    )
    assert result.upstream_prandtl_meyer_angle is None
    assert result.downstream_prandtl_meyer_angle is None
    assert result.available_turn_angle is not None


def test_finite_constant_cp_temperature_endpoint_is_inclusive() -> None:
    gas = HarmonicOscillatorGas(
        287.0, 1.4, applicable_temperature_range=(400.0, 2000.0)
    )
    zero = prandtl_meyer_expansion(2.0, 0.0, gas, upstream_temperature=500.0)
    assert zero.available_turn_angle is not None
    minimum_mach = np.sqrt(((1.0 + 0.2 * 2.0**2) / 0.8 - 1.0) / 0.2)
    analytic_limit = float(prandtl_meyer_angle(minimum_mach)) - float(
        prandtl_meyer_angle(2.0)
    )
    assert zero.available_turn_angle == pytest.approx(analytic_limit, abs=2e-14)
    result = prandtl_meyer_expansion(
        2.0, analytic_limit, gas, upstream_temperature=500.0
    )
    assert result.static_temperature_ratio == 0.8
    with pytest.raises(ModelRangeError, match="below"):
        prandtl_meyer_expansion(
            2.0,
            float(zero.available_turn_angle) + 1e-8,
            gas,
            upstream_temperature=500.0,
        )
    lower = prandtl_meyer_expansion(1.0, 0.0, gas, upstream_temperature=400.0)
    assert lower.static_pressure_ratio == lower.static_temperature_ratio == 1.0
    assert lower.available_turn_angle == 0.0
    with pytest.raises(ModelRangeError, match="below"):
        prandtl_meyer_expansion(1.0, 1e-10, gas, upstream_temperature=400.0)
    assert prandtl_meyer_expansion(
        3.0, 0.0, gas, upstream_temperature=2000.0
    ).downstream_mach == pytest.approx(3.0)


def test_thermal_broadcasting_and_owned_read_only_results() -> None:
    temperatures = np.array([500.0, 1000.0])
    result = prandtl_meyer_expansion(
        [[1.0], [3.0]],
        np.deg2rad([0.0, 5.0]),
        AIR_HARMONIC_OSCILLATOR,
        upstream_temperature=temperatures,
    )
    for field in fields(result):
        value = getattr(result, field.name)
        assert isinstance(value, np.ndarray) and value.shape == (2, 2)
        assert value.dtype == np.float64 and value.flags.owndata
        assert not value.flags.writeable
    angles = np.asarray(result.upstream_prandtl_meyer_angle)
    assert angles[0, 0] == 0.0 and np.isnan(angles[1, 1])
    expected_mach = np.asarray(result.downstream_mach).copy()
    temperatures[:] = 400.0
    assert_allclose(result.downstream_mach, expected_mach)


@pytest.mark.parametrize("temperature", [399.0, 2001.0])
def test_upstream_out_of_range_is_not_extrapolated(temperature: float) -> None:
    with pytest.raises(ModelRangeError):
        prandtl_meyer_expansion(
            2.0, 0.0, AIR_HARMONIC_OSCILLATOR, upstream_temperature=temperature
        )


@pytest.mark.parametrize("temperature", [0.0, -1.0, np.inf, np.nan, None])
def test_upstream_temperature_validation(temperature: float | None) -> None:
    with pytest.raises(ValueError, match="upstream_temperature"):
        prandtl_meyer_expansion(
            2.0, 0.1, AIR_HARMONIC_OSCILLATOR, upstream_temperature=temperature
        )
    with pytest.raises(ValueError, match="broadcastable"):
        prandtl_meyer_expansion(
            [2.0, 3.0],
            0.0,
            AIR_HARMONIC_OSCILLATOR,
            upstream_temperature=[500.0, 600.0, 700.0],
        )


def test_numerical_integration_failure_has_dedicated_exception(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import aerophysics._thermal_expansion as core

    monkeypatch.setattr(core, "quad", lambda *args, **kwargs: (np.nan, np.inf))
    with pytest.raises(ExpansionConvergenceError, match="integration"):
        prandtl_meyer_expansion(
            2.0, 0.1, AIR_HARMONIC_OSCILLATOR, upstream_temperature=500.0
        )
