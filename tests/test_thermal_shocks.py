"""Conservation, independent variable-cp solutions, and perfect-gas limits."""

from dataclasses import fields

import numpy as np
import pytest
from scipy.optimize import root

from aerophysics import AIR_HARMONIC_OSCILLATOR, AIR_NASA7, AIR_NASA9
from aerophysics.exceptions import ModelRangeError, NoAttachedShockError
from aerophysics.gas import PerfectGas
from aerophysics.real_gas import HarmonicOscillatorGas
from aerophysics.shocks import (
    ShockBranch,
    maximum_attached_deflection,
    normal_shock,
    oblique_shock,
    shock_angle,
    theta_from_shock_angle,
)
from aerophysics.thermochemistry import (
    IdealGasSpecies,
    NASA7Polynomial,
    ThermallyPerfectGas,
)


def polynomial_gas(slope: float = 0.0) -> ThermallyPerfectGas:
    """Manufactured cp/R = 3.5 + slope*T over 100--30000 K."""
    return ThermallyPerfectGas(
        (
            IdealGasSpecies(
                "test",
                0.028,
                NASA7Polynomial(
                    (100.0, 30000.0), ((3.5, slope, 0.0, 0.0, 0.0, 0.0, 0.0),)
                ),
            ),
        ),
        (1.0,),
    )


@pytest.mark.parametrize("gamma", [1.2, 1.4, 5.0 / 3.0])
@pytest.mark.parametrize("branch", list(ShockBranch))
def test_constant_cp_recovers_closed_form(gamma: float, branch: ShockBranch) -> None:
    thermal = HarmonicOscillatorGas(287.0, gamma)
    perfect = PerfectGas(heat_capacity_ratio=gamma, specific_gas_constant=287.0)
    for mach in (1.2, 2.0, 5.0, 10.0):
        expected_normal = normal_shock(mach, perfect)
        actual_normal = normal_shock(mach, thermal, upstream_temperature=300.0)
        for field in fields(expected_normal):
            assert getattr(actual_normal, field.name) == pytest.approx(
                getattr(expected_normal, field.name), rel=2e-10, abs=1e-12
            )
        limit = maximum_attached_deflection(mach, perfect)
        thermal_limit = maximum_attached_deflection(
            mach, thermal, upstream_temperature=300.0
        )
        assert thermal_limit.deflection_angle == pytest.approx(
            limit.deflection_angle, abs=1e-11
        )
        theta = 0.5 * float(limit.deflection_angle)
        expected = oblique_shock(mach, theta, branch, perfect)
        actual = oblique_shock(mach, theta, branch, thermal, upstream_temperature=300.0)
        for field in fields(expected):
            assert getattr(actual, field.name) == pytest.approx(
                getattr(expected, field.name), rel=2e-9, abs=1e-11
            )


def test_constant_nasa_polynomial_recovers_closed_form() -> None:
    thermal = polynomial_gas()
    expected = oblique_shock(3.0, np.deg2rad(20.0))
    actual = oblique_shock(
        3.0, np.deg2rad(20.0), gas=thermal, upstream_temperature=300.0
    )
    for field in fields(expected):
        assert getattr(actual, field.name) == pytest.approx(
            getattr(expected, field.name), rel=1e-10
        )


@pytest.mark.parametrize("gas", [AIR_NASA7, AIR_NASA9, AIR_HARMONIC_OSCILLATOR])
@pytest.mark.parametrize("branch", list(ShockBranch))
@pytest.mark.parametrize(
    "mach,temperature,theta", [(2.5, 500.0, 15.0), (3.0, 500.0, 20.0)]
)
def test_shock_conserves_mass_momentum_energy_and_tangent_velocity(
    gas: ThermallyPerfectGas | HarmonicOscillatorGas,
    branch: ShockBranch,
    mach: float,
    temperature: float,
    theta: float,
) -> None:
    result = oblique_shock(
        mach, np.deg2rad(theta), branch, gas, upstream_temperature=temperature
    )
    beta = float(result.shock_angle)
    deflection = float(result.deflection_angle)
    t2 = temperature * float(result.static_temperature_ratio)
    r = gas.specific_gas_constant
    p1 = 43000.0
    p2 = p1 * float(result.static_pressure_ratio)
    rho1, rho2 = p1 / (r * temperature), p2 / (r * t2)
    v1 = mach * float(gas.speed_of_sound(temperature))
    v2 = float(result.downstream_mach) * float(gas.speed_of_sound(t2))
    u1, u2 = v1 * np.sin(beta), v2 * np.sin(beta - deflection)
    assert rho1 * u1 == pytest.approx(rho2 * u2, rel=2e-11)
    assert p1 + rho1 * u1**2 == pytest.approx(p2 + rho2 * u2**2, rel=2e-11)
    assert gas.standard_enthalpy(temperature) + 0.5 * v1**2 == pytest.approx(
        gas.standard_enthalpy(t2) + 0.5 * v2**2, rel=2e-11
    )
    assert v1 * np.cos(beta) == pytest.approx(v2 * np.cos(beta - deflection), rel=2e-11)
    ds_over_r = (float(gas.entropy(t2, p2)) - float(gas.entropy(temperature, p1))) / r
    assert ds_over_r > 0.0
    assert result.total_pressure_ratio == pytest.approx(np.exp(-ds_over_r), rel=2e-12)
    assert theta_from_shock_angle(
        mach, beta, gas, upstream_temperature=temperature
    ) == pytest.approx(deflection, abs=2e-12)


@pytest.mark.parametrize("branch", list(ShockBranch))
def test_variable_cp_against_independent_coupled_equations(branch: ShockBranch) -> None:
    gas = polynomial_gas(0.001)
    mach, t1, theta = 4.0, 500.0, np.deg2rad(20.0)
    gamma1 = 4.0 / 3.0
    speed_squared_over_r = mach**2 * gamma1 * t1

    def equations(values: np.ndarray) -> list[float]:
        beta, tau, density = values
        normal_speed_squared_over_r = speed_squared_over_r * np.sin(beta) ** 2
        return [
            (
                density * tau
                - 1.0
                - normal_speed_squared_over_r / t1 * (1.0 - 1.0 / density)
            ),
            (
                3.5 * t1 * (tau - 1.0)
                + 0.0005 * t1**2 * (tau**2 - 1.0)
                - 0.5 * normal_speed_squared_over_r * (1.0 - 1.0 / density**2)
            )
            / t1,
            np.tan(beta - theta) - np.tan(beta) / density,
        ]

    guess = [0.6, 1.5, 3.0] if branch is ShockBranch.WEAK else [1.45, 3.0, 6.0]
    reference = root(equations, guess, tol=1e-12)
    assert reference.success
    np.testing.assert_allclose(equations(reference.x), 0.0, atol=2e-12)
    beta, tau, density = reference.x
    actual = oblique_shock(mach, theta, branch, gas, upstream_temperature=t1)
    assert actual.shock_angle == pytest.approx(beta, abs=2e-11)
    assert actual.static_temperature_ratio == pytest.approx(tau, rel=2e-11)
    assert actual.static_density_ratio == pytest.approx(density, rel=2e-11)
    # Independently integrate this manufactured cp/R, then include p2/p1.
    ds_over_r = 3.5 * np.log(tau) + 0.001 * t1 * (tau - 1.0) - np.log(density * tau)
    assert actual.total_pressure_ratio == pytest.approx(np.exp(-ds_over_r), rel=2e-11)


def test_thermal_broadcast_and_scalar_contract() -> None:
    result = oblique_shock(
        [[2.0], [3.0]],
        np.deg2rad(10.0),
        gas=AIR_NASA9,
        upstream_temperature=[300.0, 600.0, 1000.0],
    )
    assert np.shape(result.shock_angle) == (2, 3)
    for i, mach in enumerate((2.0, 3.0)):
        for j, temperature in enumerate((300.0, 600.0, 1000.0)):
            single = oblique_shock(
                mach, np.deg2rad(10.0), gas=AIR_NASA9, upstream_temperature=temperature
            )
            assert isinstance(single.shock_angle, float)
            for field in fields(single):
                assert np.asarray(getattr(result, field.name))[i, j] == pytest.approx(
                    getattr(single, field.name)
                )
    assert np.shape(
        normal_shock(
            2.0, AIR_NASA9, upstream_temperature=[300.0, 500.0]
        ).downstream_mach
    ) == (2,)
    assert np.shape(
        theta_from_shock_angle(
            3.0, np.deg2rad(40.0), AIR_NASA9, upstream_temperature=[300.0, 500.0]
        )
    ) == (2,)
    assert np.shape(
        maximum_attached_deflection(
            3.0, AIR_NASA9, upstream_temperature=[300.0, 500.0]
        ).shock_angle
    ) == (2,)


@pytest.mark.parametrize("gas", [AIR_NASA7, AIR_NASA9, AIR_HARMONIC_OSCILLATOR])
def test_endpoints_sonic_and_detachment(
    gas: ThermallyPerfectGas | HarmonicOscillatorGas,
) -> None:
    for mach in (1.0, 3.0):
        weak = oblique_shock(mach, 0.0, gas=gas, upstream_temperature=500.0)
        assert weak.shock_angle == pytest.approx(np.arcsin(1.0 / mach))
        assert weak.downstream_mach == pytest.approx(mach)
        assert weak.static_temperature_ratio == 1.0
        assert weak.total_pressure_ratio == 1.0
        strong = oblique_shock(
            mach, 0.0, ShockBranch.STRONG, gas, upstream_temperature=500.0
        )
        normal = normal_shock(mach, gas, upstream_temperature=500.0)
        assert strong.shock_angle == 0.5 * np.pi
        assert strong.downstream_mach == normal.downstream_mach
        assert strong.total_pressure_ratio == normal.total_pressure_ratio
        assert (
            theta_from_shock_angle(mach, 0.5 * np.pi, gas, upstream_temperature=500.0)
            == 0.0
        )
    limit = maximum_attached_deflection(3.0, gas, upstream_temperature=500.0)
    for branch in ShockBranch:
        assert (
            shock_angle(
                3.0, limit.deflection_angle, branch, gas, upstream_temperature=500.0
            )
            == limit.shock_angle
        )
        with pytest.raises(NoAttachedShockError):
            oblique_shock(
                3.0,
                float(limit.deflection_angle) + 0.001,
                branch,
                gas,
                upstream_temperature=500.0,
            )
    assert (
        maximum_attached_deflection(
            1.0, gas, upstream_temperature=500.0
        ).deflection_angle
        == 0.0
    )


def test_near_sonic_and_vanishing_deflection() -> None:
    result = normal_shock(1.000001, AIR_NASA9, upstream_temperature=300.0)
    assert 0.999 < float(result.downstream_mach) < 1.0
    assert 1.0 < float(result.static_temperature_ratio) < 1.00001
    weak = oblique_shock(3.0, 1e-7, gas=AIR_NASA9, upstream_temperature=300.0)
    assert float(weak.static_pressure_ratio) > 1.0
    assert theta_from_shock_angle(
        3.0, weak.shock_angle, AIR_NASA9, upstream_temperature=300.0
    ) == pytest.approx(1e-7, abs=1e-12)


def test_temperature_range_does_not_discard_valid_weak_shock() -> None:
    result = oblique_shock(
        20.0, np.deg2rad(10.0), gas=AIR_NASA7, upstream_temperature=300.0
    )
    assert 300.0 < 300.0 * float(result.static_temperature_ratio) < 6000.0
    with pytest.raises(ModelRangeError):
        normal_shock(20.0, AIR_NASA7, upstream_temperature=300.0)
    with pytest.raises(ModelRangeError):
        maximum_attached_deflection(20.0, AIR_NASA7, upstream_temperature=300.0)
    with pytest.raises(ModelRangeError):
        oblique_shock(
            20.0,
            np.deg2rad(10.0),
            ShockBranch.STRONG,
            AIR_NASA7,
            upstream_temperature=300.0,
        )
    with pytest.raises(ModelRangeError):
        oblique_shock(20.0, np.deg2rad(40.0), gas=AIR_NASA7, upstream_temperature=300.0)
    with pytest.raises(ModelRangeError):
        oblique_shock(3.0, 0.1, gas=AIR_NASA9, upstream_temperature=6000.0)
    with pytest.raises(ModelRangeError):
        theta_from_shock_angle(20.0, 0.5 * np.pi, AIR_NASA9, upstream_temperature=300.0)


@pytest.mark.parametrize("temperature", [0.0, -1.0, np.nan, np.inf, 199.0, 6001.0])
def test_invalid_temperature(temperature: float) -> None:
    with pytest.raises(ValueError):
        oblique_shock(3.0, 0.1, gas=AIR_NASA9, upstream_temperature=temperature)


def test_missing_temperature_and_invalid_inputs() -> None:
    with pytest.raises(ValueError, match="upstream_temperature is required"):
        oblique_shock(3.0, 0.1, gas=AIR_NASA9)
    with pytest.raises(ValueError, match="broadcastable"):
        oblique_shock(
            [2.0, 3.0], 0.1, gas=AIR_NASA9, upstream_temperature=[300.0, 400.0, 500.0]
        )
    with pytest.raises(ValueError, match="non-negative"):
        oblique_shock(3.0, -0.1, gas=AIR_NASA9, upstream_temperature=300.0)
    with pytest.raises(ValueError, match="ShockBranch"):
        shock_angle(3.0, 0.1, "weak", AIR_NASA9, upstream_temperature=300.0)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="shock_angle"):
        theta_from_shock_angle(3.0, 0.1, AIR_NASA9, upstream_temperature=300.0)
    with pytest.raises(TypeError, match="gas must"):
        normal_shock(3.0, object(), upstream_temperature=300.0)  # type: ignore[arg-type]
    nonphysical = ThermallyPerfectGas(
        (
            IdealGasSpecies(
                "bad",
                0.028,
                NASA7Polynomial(
                    (100.0, 30000.0), ((0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),)
                ),
            ),
        ),
        (1.0,),
    )
    with pytest.raises(ModelRangeError, match="heat capacity"):
        normal_shock(3.0, nonphysical, upstream_temperature=300.0)
    limited = HarmonicOscillatorGas(
        287.0, 1.4, applicable_temperature_range=(250.0, 500.0)
    )
    with pytest.raises(ModelRangeError):
        normal_shock(3.0, limited, upstream_temperature=300.0)


def test_static_temperature_changes_solution_and_defaults_stay_perfect() -> None:
    cold = oblique_shock(
        3.0, np.deg2rad(20.0), gas=AIR_NASA9, upstream_temperature=300.0
    )
    hot = oblique_shock(
        3.0, np.deg2rad(20.0), gas=AIR_NASA9, upstream_temperature=1500.0
    )
    assert float(hot.shock_angle) < float(cold.shock_angle)
    assert float(hot.static_temperature_ratio) < float(cold.static_temperature_ratio)
    assert oblique_shock(3.0, 0.1) == oblique_shock(
        3.0, 0.1, upstream_temperature=500.0
    )
