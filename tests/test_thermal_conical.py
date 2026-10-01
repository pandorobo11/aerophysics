"""Thermal cone verification: perfect limits and independent variable-cp flow."""

from dataclasses import fields

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root

from aerophysics import AIR_HARMONIC_OSCILLATOR, AIR_NASA7, AIR_NASA9
from aerophysics.exceptions import (
    ModelRangeError,
    NoAttachedShockError,
    ShockConvergenceError,
)
from aerophysics.gas import PerfectGas
from aerophysics.real_gas import HarmonicOscillatorGas
from aerophysics.shocks import (
    conical_shock,
    maximum_attached_cone_angle,
    normal_shock,
)
from aerophysics.thermochemistry import (
    IdealGasSpecies,
    NASA7Polynomial,
    ThermallyPerfectGas,
)


def polynomial_gas() -> ThermallyPerfectGas:
    """cp/R = 3.5 + 0.001*T, with a deliberately nonzero enthalpy offset."""
    return ThermallyPerfectGas(
        (
            IdealGasSpecies(
                "test",
                0.028,
                NASA7Polynomial(
                    (100.0, 30000.0), ((3.5, 0.001, 0.0, 0.0, 0.0, -12000.0, 0.0),)
                ),
            ),
        ),
        (1.0,),
    )


@pytest.mark.parametrize("gamma", [1.2, 1.4, 5.0 / 3.0])
@pytest.mark.parametrize("mach", [1.2, 3.0, 10.0])
def test_constant_cp_recovers_perfect_cone(gamma: float, mach: float) -> None:
    thermal = HarmonicOscillatorGas(287.0, gamma)
    perfect = PerfectGas(287.0, gamma)
    limit = maximum_attached_cone_angle(mach, perfect)
    angle = 0.4 * float(limit.cone_half_angle)
    actual = conical_shock(mach, angle, thermal, upstream_temperature=500.0)
    expected = conical_shock(mach, angle, perfect)
    for field in fields(expected):
        assert getattr(actual, field.name) == pytest.approx(
            getattr(expected, field.name), rel=3e-9, abs=2e-10
        )
    thermal_limit = maximum_attached_cone_angle(
        mach, thermal, upstream_temperature=500.0
    )
    assert thermal_limit.cone_half_angle == pytest.approx(
        limit.cone_half_angle, abs=2e-10
    )
    assert thermal_limit.shock_angle == pytest.approx(limit.shock_angle, abs=5e-8)


def test_variable_cp_against_independent_temperature_ode() -> None:
    # The reference uses a coupled shock solve and a third ODE for temperature,
    # with analytic caloric properties; no production Hugoniot or h(T) inversion.
    mach, t1, angle = 4.0, 500.0, np.deg2rad(15.0)
    gamma1 = 4.0 / 3.0
    speed_squared_over_r = mach**2 * gamma1 * t1

    def reference(beta: float) -> tuple[float, np.ndarray, float, float]:
        un_squared_over_r = speed_squared_over_r * np.sin(beta) ** 2

        def shock_equations(values: np.ndarray) -> list[float]:
            tau, density = values
            return [
                density * tau - 1.0 - un_squared_over_r / t1 * (1.0 - 1.0 / density),
                (
                    3.5 * t1 * (tau - 1.0)
                    + 0.0005 * t1**2 * (tau**2 - 1.0)
                    - 0.5 * un_squared_over_r * (1.0 - 1.0 / density**2)
                )
                / t1,
            ]

        shock = root(shock_equations, [1.5, 3.0], tol=1e-11)
        assert shock.success
        tau, density = shock.x
        assert tau > 1.0 and density > 1.0
        initial = [np.cos(beta), -np.sin(beta) / density, t1 * tau]

        def rhs(theta: float, state: np.ndarray) -> list[float]:
            radial, polar, temperature = state
            cp_over_r = 3.5 + 0.001 * temperature
            a_squared = (
                cp_over_r / (cp_over_r - 1.0) * temperature / speed_squared_over_r
            )
            derivative = (
                radial * polar**2 - a_squared * (2 * radial + polar / np.tan(theta))
            ) / (a_squared - polar**2)
            return [
                polar,
                derivative,
                -speed_squared_over_r * polar * (radial + derivative) / cp_over_r,
            ]

        def surface(_theta: float, state: np.ndarray) -> float:
            return float(state[1])

        surface.terminal = True  # type: ignore[attr-defined]
        flow = solve_ivp(
            rhs,
            (beta, 1e-8),
            initial,
            events=surface,
            method="RK45",
            rtol=2e-11,
            atol=2e-12,
            max_step=0.002,
        )
        assert flow.success and flow.t_events is not None and flow.y_events is not None
        return float(flow.t_events[0][0]), flow.y_events[0][0], tau, density

    beta = brentq(lambda b: reference(b)[0] - angle, 0.3, 0.6, xtol=1e-12)
    _, state, tau, density = reference(beta)
    radial, polar, ts = state
    cp_over_r = 3.5 + 0.001 * ts
    surface_mach = np.sqrt(
        speed_squared_over_r
        * (radial**2 + polar**2)
        / (cp_over_r / (cp_over_r - 1.0) * ts)
    )
    pressure_ratio = (
        density * tau * np.exp(3.5 * np.log(ts / (t1 * tau)) + 0.001 * (ts - t1 * tau))
    )
    ds_over_r = 3.5 * np.log(tau) + 0.001 * t1 * (tau - 1.0) - np.log(density * tau)
    actual = conical_shock(mach, angle, polynomial_gas(), upstream_temperature=t1)
    assert actual.shock_angle == pytest.approx(beta, abs=2e-10)
    assert actual.surface_temperature_ratio == pytest.approx(ts / t1, rel=3e-9)
    assert actual.surface_pressure_ratio == pytest.approx(pressure_ratio, rel=3e-9)
    assert actual.surface_mach == pytest.approx(surface_mach, rel=3e-9)
    assert actual.total_pressure_ratio == pytest.approx(np.exp(-ds_over_r), rel=3e-9)


@pytest.mark.parametrize("gas", [AIR_NASA7, AIR_NASA9, AIR_HARMONIC_OSCILLATOR])
def test_presets_conserve_energy_entropy_and_equation_of_state(
    gas: ThermallyPerfectGas | HarmonicOscillatorGas,
) -> None:
    t1, mach = 1000.0, 3.0
    actual = conical_shock(mach, np.deg2rad(10.0), gas, upstream_temperature=t1)
    ts = t1 * float(actual.surface_temperature_ratio)
    v1 = mach * float(gas.speed_of_sound(t1))
    vs = float(actual.surface_mach) * float(gas.speed_of_sound(ts))
    assert gas.standard_enthalpy(t1) + 0.5 * v1**2 == pytest.approx(
        gas.standard_enthalpy(ts) + 0.5 * vs**2, rel=3e-10
    )
    pressure = 43000.0
    ds_over_r = (
        float(gas.entropy(ts, pressure * float(actual.surface_pressure_ratio)))
        - float(gas.entropy(t1, pressure))
    ) / gas.specific_gas_constant
    assert ds_over_r > 0.0
    assert actual.total_pressure_ratio == pytest.approx(np.exp(-ds_over_r), rel=3e-10)
    assert actual.surface_pressure_ratio == pytest.approx(
        float(actual.surface_density_ratio) * float(actual.surface_temperature_ratio),
        rel=2e-12,
    )
    assert t1 < ts <= gas.temperature_range[1]
    assert actual.surface_mach < actual.post_shock_mach < mach


def test_range_limited_polar_keeps_valid_weak_cone() -> None:
    gas = HarmonicOscillatorGas(287.0, 1.4, applicable_temperature_range=(200.0, 600.0))
    expected = conical_shock(3.0, np.deg2rad(10.0), PerfectGas(287.0, 1.4))
    actual = conical_shock(3.0, np.deg2rad(10.0), gas, upstream_temperature=500.0)
    assert actual.shock_angle == pytest.approx(expected.shock_angle, abs=2e-10)
    assert 500.0 * float(actual.surface_temperature_ratio) < 600.0
    with pytest.raises(ModelRangeError):
        normal_shock(3.0, gas, upstream_temperature=500.0)
    with pytest.raises(ModelRangeError):
        maximum_attached_cone_angle(3.0, gas, upstream_temperature=500.0)
    with pytest.raises(ModelRangeError):
        conical_shock(3.0, np.deg2rad(20.0), gas, upstream_temperature=500.0)


def test_surface_can_exceed_range_when_post_shock_temperature_is_valid() -> None:
    # A 10-degree Mach-3 cone has T_shock~545 K but T_surface~567 K.
    gas = HarmonicOscillatorGas(287.0, 1.4, applicable_temperature_range=(200.0, 550.0))
    perfect = conical_shock(3.0, np.deg2rad(10.0), PerfectGas(287.0, 1.4))
    shock = normal_shock(
        3.0 * np.sin(float(perfect.shock_angle)), PerfectGas(287.0, 1.4)
    )
    assert 500.0 * float(shock.static_temperature_ratio) < 550.0
    assert 500.0 * float(perfect.surface_temperature_ratio) > 550.0
    with pytest.raises(ModelRangeError):
        conical_shock(3.0, np.deg2rad(10.0), gas, upstream_temperature=500.0)


def test_thermal_attached_limit_and_detachment() -> None:
    gas = polynomial_gas()
    limit = maximum_attached_cone_angle(3.0, gas, upstream_temperature=500.0)
    actual = conical_shock(3.0, limit.cone_half_angle, gas, upstream_temperature=500.0)
    assert actual.shock_angle == limit.shock_angle
    with pytest.raises(NoAttachedShockError):
        conical_shock(
            3.0, float(limit.cone_half_angle) + 1e-6, gas, upstream_temperature=500.0
        )


def test_broadcast_scalar_and_zero_cone_contract() -> None:
    gas = AIR_HARMONIC_OSCILLATOR
    result = conical_shock(
        [[2.0], [3.0]],
        np.deg2rad([0.0, 10.0]),
        gas,
        upstream_temperature=[500.0, 1000.0],
    )
    assert np.shape(result.shock_angle) == (2, 2)
    for i, mach in enumerate((2.0, 3.0)):
        for j, (angle, temperature) in enumerate(((0.0, 500.0), (10.0, 1000.0))):
            scalar = conical_shock(
                mach, np.deg2rad(angle), gas, upstream_temperature=temperature
            )
            for field in fields(scalar):
                assert isinstance(getattr(scalar, field.name), float)
                assert getattr(result, field.name)[i, j] == getattr(scalar, field.name)
    assert np.all(np.asarray(result.surface_temperature_ratio)[:, 0] == 1.0)
    assert np.all(np.asarray(result.surface_pressure_ratio)[:, 0] == 1.0)
    # An upstream state at Tmax allows the zero-strength endpoint only.
    zero = conical_shock(3.0, 0.0, gas, upstream_temperature=2000.0)
    assert zero.surface_temperature_ratio == 1.0
    with pytest.raises(ModelRangeError):
        conical_shock(3.0, 0.01, gas, upstream_temperature=2000.0)
    limit = maximum_attached_cone_angle(
        [[2.0], [3.0]],
        HarmonicOscillatorGas(287.0, 1.4),
        upstream_temperature=[300.0, 600.0],
    )
    assert np.shape(limit.cone_half_angle) == (2, 2)


@pytest.mark.parametrize("solver", [conical_shock, maximum_attached_cone_angle])
def test_thermal_input_validation(solver: object) -> None:
    # Keep calls explicit so both public signatures are type checked.
    if solver is conical_shock:
        with pytest.raises(ValueError, match="upstream_temperature"):
            conical_shock(3.0, 0.1, AIR_NASA9)
        with pytest.raises(ModelRangeError):
            conical_shock(3.0, 0.0, AIR_NASA9, upstream_temperature=100.0)
        with pytest.raises(ValueError, match="broadcast"):
            conical_shock(
                [2.0, 3.0], 0.1, AIR_NASA9, upstream_temperature=[300.0, 400.0, 500.0]
            )
    else:
        with pytest.raises(ValueError, match="upstream_temperature"):
            maximum_attached_cone_angle(3.0, AIR_NASA9)
        with pytest.raises(ModelRangeError):
            maximum_attached_cone_angle(3.0, AIR_NASA9, upstream_temperature=100.0)


@pytest.mark.parametrize("temperature_margin", [0.001, 0.1])
def test_physical_peak_just_inside_surface_temperature_range(
    temperature_margin: float,
) -> None:
    perfect = PerfectGas(287.0, 1.4)
    limit = maximum_attached_cone_angle(3.0, perfect)
    peak = conical_shock(3.0, limit.cone_half_angle, perfect)
    tmax = 500.0 * float(peak.surface_temperature_ratio) + temperature_margin
    gas = HarmonicOscillatorGas(287.0, 1.4, applicable_temperature_range=(200.0, tmax))
    actual = maximum_attached_cone_angle(3.0, gas, upstream_temperature=500.0)
    assert actual.cone_half_angle == pytest.approx(limit.cone_half_angle, abs=2e-10)
    with pytest.raises(NoAttachedShockError):
        conical_shock(
            3.0, float(actual.cone_half_angle) + 1e-6, gas, upstream_temperature=500.0
        )


def test_surface_temperature_endpoint_and_true_exceed() -> None:
    unbounded = HarmonicOscillatorGas(287.0, 1.4)
    angle = np.deg2rad(10.0)
    expected = conical_shock(3.0, angle, unbounded, upstream_temperature=500.0)
    tmax = 500.0 * float(expected.surface_temperature_ratio)
    bounded = HarmonicOscillatorGas(
        287.0, 1.4, applicable_temperature_range=(200.0, tmax)
    )
    actual = conical_shock(3.0, angle, bounded, upstream_temperature=500.0)
    assert 500.0 * float(actual.surface_temperature_ratio) <= tmax
    assert actual.shock_angle == pytest.approx(expected.shock_angle, abs=2e-10)
    with pytest.raises(ModelRangeError):
        conical_shock(3.0, angle + 1e-6, bounded, upstream_temperature=500.0)


@pytest.mark.parametrize("temperature_margin", [0.001, 0.1])
def test_temperature_boundary_before_physical_peak_is_still_range_limited(
    temperature_margin: float,
) -> None:
    unbounded = HarmonicOscillatorGas(287.0, 1.4)
    peak = maximum_attached_cone_angle(3.0, unbounded, upstream_temperature=500.0)
    state = conical_shock(
        3.0, peak.cone_half_angle, unbounded, upstream_temperature=500.0
    )
    gas = HarmonicOscillatorGas(
        287.0,
        1.4,
        applicable_temperature_range=(
            200.0,
            500.0 * float(state.surface_temperature_ratio) - temperature_margin,
        ),
    )
    with pytest.raises(ModelRangeError):
        maximum_attached_cone_angle(3.0, gas, upstream_temperature=500.0)
    with pytest.raises(ModelRangeError):
        conical_shock(
            3.0, float(peak.cone_half_angle) + 1e-6, gas, upstream_temperature=500.0
        )


@pytest.mark.parametrize("angle", [0.01, 0.1])
def test_slender_cone_returns_solution_or_typed_numerical_failure(angle: float) -> None:
    # Near-Mach-wave roots depend on floating-point/backend details. A failure
    # must be identifiable without treating it as range failure or detachment.
    try:
        result = conical_shock(
            3.0, np.deg2rad(angle), AIR_HARMONIC_OSCILLATOR, upstream_temperature=500.0
        )
    except ShockConvergenceError as error:
        assert str(error)
    else:
        assert (
            np.arcsin(1.0 / 3.0)
            < result.shock_angle
            < np.deg2rad(10.0) + np.arcsin(1.0 / 3.0)
        )
        assert 1.0 < result.surface_temperature_ratio < 1.01
    regular = conical_shock(
        3.0, np.deg2rad(10.0), AIR_HARMONIC_OSCILLATOR, upstream_temperature=500.0
    )
    assert regular.surface_temperature_ratio > 1.1
