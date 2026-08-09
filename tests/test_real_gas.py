"""Tests for harmonic-oscillator and Beattie--Bridgeman gas models."""

import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics import (
    AIR_BEATTIE_BRIDGEMAN,
    AIR_HARMONIC_OSCILLATOR,
    BeattieBridgemanGas,
    HarmonicOscillatorGas,
    ThermodynamicState,
    VibrationalMode,
)
from aerophysics.exceptions import ApplicabilityWarning, ModelRangeError
from aerophysics.isentropic import (
    MachBranch,
    _real_flow_state,
    area_ratio,
    choked_mass_flux,
    critical_ratios,
    isentropic_ratios,
    isentropic_state,
    mach_from_area_ratio,
    mach_from_total_density_ratio,
    mach_from_total_pressure_ratio,
    mach_from_total_temperature_ratio,
    mass_flow_parameter,
    mass_flux,
)


def _zero_correction_pair() -> tuple[HarmonicOscillatorGas, BeattieBridgemanGas]:
    modes = (VibrationalMode(0.8, 3000.0), VibrationalMode(0.2, 2200.0))
    harmonic = HarmonicOscillatorGas(287.0, 1.4, modes)
    real = BeattieBridgemanGas(287.0, 1.4, 0.0, 0.0, 0.0, 0.0, 0.0, modes)
    return harmonic, real


def test_vibrational_mode_and_model_validation() -> None:
    assert VibrationalMode(1, 3000) == VibrationalMode(1.0, 3000.0)
    with pytest.raises(ValueError, match="weight"):
        VibrationalMode(-1.0, 3000.0)
    with pytest.raises(ValueError, match="characteristic_temperature"):
        VibrationalMode(1.0, 0.0)
    with pytest.raises(ValueError, match="specific_gas_constant"):
        HarmonicOscillatorGas(0.0, 1.4)
    with pytest.raises(ValueError, match="base_heat_capacity_ratio"):
        HarmonicOscillatorGas(287.0, 1.0)
    with pytest.raises(TypeError, match="VibrationalMode"):
        HarmonicOscillatorGas(287.0, 1.4, (object(),))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="strictly increasing"):
        HarmonicOscillatorGas(287.0, 1.4, applicable_temperature_range=(2, 1))
    with pytest.raises(ValueError, match="a0"):
        BeattieBridgemanGas(287.0, 1.4, np.inf, 0, 0, 0, 0)


def test_harmonic_oscillator_matches_independent_kennard_equation() -> None:
    gas = AIR_HARMONIC_OSCILLATOR
    temperature = np.array([400.0, 1000.0, 2000.0])
    theta = 3055.56
    x = theta / temperature
    expected_cv = 287.05287 / 0.4 + 287.05287 * (x**2 * np.exp(x) / np.expm1(x) ** 2)
    assert_allclose(gas.cv(temperature), expected_cv, rtol=2e-15)
    assert_allclose(gas.cp(temperature), expected_cv + 287.05287, rtol=2e-15)
    assert_allclose(
        gas.heat_capacity_ratio(temperature),
        (expected_cv + 287.05287) / expected_cv,
    )


def test_harmonic_oscillator_energy_entropy_and_state_identities() -> None:
    gas = AIR_HARMONIC_OSCILLATOR
    temperature = 1200.0
    pressure = 6.0e6
    state = gas.state(temperature, pressure)
    assert isinstance(state, ThermodynamicState)
    assert state.density == pytest.approx(
        pressure / (gas.specific_gas_constant * temperature)
    )
    assert state.enthalpy - state.internal_energy == pytest.approx(
        gas.specific_gas_constant * temperature
    )
    assert state.cp - state.cv == pytest.approx(gas.specific_gas_constant)
    assert state.heat_capacity_ratio == pytest.approx(state.cp / state.cv)
    assert state.speed_of_sound == pytest.approx(
        np.sqrt(state.heat_capacity_ratio * gas.specific_gas_constant * temperature)
    )
    delta = 1.0e-3
    derivative = (
        gas.standard_enthalpy(temperature + delta)
        - gas.standard_enthalpy(temperature - delta)
    ) / (2.0 * delta)
    assert derivative == pytest.approx(state.cp, rel=2e-9)
    assert gas.standard_internal_energy(temperature) == state.internal_energy
    assert gas.entropy(temperature, pressure) == state.entropy
    assert gas.speed_of_sound(temperature) == state.speed_of_sound


def test_harmonic_oscillator_arrays_warnings_and_validation() -> None:
    gas = AIR_HARMONIC_OSCILLATOR
    state = gas.state(np.array([[800.0], [1200.0]]), np.array([2.0e6, 4.0e6]))
    assert isinstance(state.temperature, np.ndarray)
    assert state.temperature.shape == (2, 2)
    assert state.temperature.dtype == np.float64
    with pytest.warns(ApplicabilityWarning, match="documented range"):
        gas.cp(300.0, allow_extrapolation=True)
    with pytest.raises(ModelRangeError, match="400--2000"):
        gas.cv(300.0)
    with pytest.raises(ValueError, match="temperature"):
        gas.state(0.0, 1.0e6)
    with pytest.raises(ValueError, match="pressure"):
        gas.entropy(1000.0, 0.0)
    with pytest.raises(ValueError, match="broadcast"):
        gas.state(np.ones(2) * 1000.0, np.ones(3) * 1.0e6)


def test_beattie_bridgeman_eos_and_thermodynamic_identities() -> None:
    gas = AIR_BEATTIE_BRIDGEMAN
    temperature = 1200.0
    pressure = 6.0e6
    density = gas.density(temperature, pressure)
    assert gas.pressure(temperature, density) == pytest.approx(pressure, rel=2e-15)
    state = gas.state(temperature, pressure)
    assert state.enthalpy == pytest.approx(state.internal_energy + pressure / density)
    assert state.cp > state.cv > 0.0
    assert state.heat_capacity_ratio == pytest.approx(state.cp / state.cv)
    assert state.speed_of_sound > 0.0
    assert gas.cp(temperature, pressure) == state.cp
    assert gas.cv(temperature, pressure) == state.cv
    assert gas.heat_capacity_ratio(temperature, pressure) == state.heat_capacity_ratio
    assert gas.internal_energy(temperature, pressure) == state.internal_energy
    assert gas.enthalpy(temperature, pressure) == state.enthalpy
    assert gas.entropy(temperature, pressure) == state.entropy
    assert gas.speed_of_sound(temperature, pressure) == state.speed_of_sound


def test_beattie_bridgeman_arrays_applicability_and_validation() -> None:
    gas = AIR_BEATTIE_BRIDGEMAN
    temperatures = np.array([[800.0], [1200.0]])
    pressures = np.array([2.0e6, 8.0e6])
    densities = gas.density(temperatures, pressures)
    assert isinstance(densities, np.ndarray)
    assert densities.shape == (2, 2)
    assert_allclose(
        gas.pressure(temperatures, densities), np.broadcast_to(pressures, (2, 2))
    )
    with pytest.warns(ApplicabilityWarning, match="outside"):
        gas.state(1300.0, 20.0e6, allow_extrapolation=True)
    with pytest.raises(ModelRangeError, match="temperature"):
        gas.state(30.0, 6.0e6)
    with pytest.raises(ModelRangeError, match="pressure"):
        gas.density(1200.0, 100.0)
    with pytest.raises(ValueError, match="density"):
        gas.pressure(1200.0, 0.0)
    with pytest.raises(ValueError, match="broadcast"):
        gas.pressure(np.ones(2) * 1200.0, np.ones(3))


def test_beattie_bridgeman_randall_tabulated_range_endpoints() -> None:
    gas = AIR_BEATTIE_BRIDGEMAN
    temperature_range = gas.applicable_temperature_range
    pressure_range = gas.applicable_pressure_range
    assert temperature_range is not None
    assert pressure_range is not None
    temperature_minimum, temperature_maximum = temperature_range
    pressure_minimum, pressure_maximum = pressure_range

    # The low-temperature endpoint has a stable gas root at low pressure.
    for temperature in (temperature_minimum, temperature_maximum):
        gas.state(temperature, pressure_minimum)
    for pressure in (pressure_minimum, pressure_maximum):
        gas.state(900.0, pressure)

    for temperature in (
        np.nextafter(temperature_minimum, 0.0),
        np.nextafter(temperature_maximum, np.inf),
    ):
        with pytest.raises(ModelRangeError, match="temperature"):
            gas.state(temperature, pressure_minimum, allow_extrapolation=False)
        with pytest.warns(
            ApplicabilityWarning, match=r"Randall, AEDC-TR-57-8.*tabulated range"
        ):
            gas.state(temperature, pressure_minimum, allow_extrapolation=True)

    for pressure in (
        np.nextafter(pressure_minimum, 0.0),
        np.nextafter(pressure_maximum, np.inf),
    ):
        with pytest.raises(ModelRangeError, match="pressure"):
            gas.state(900.0, pressure, allow_extrapolation=False)
        with pytest.warns(
            ApplicabilityWarning, match=r"Randall, AEDC-TR-57-8.*tabulated range"
        ):
            gas.state(900.0, pressure, allow_extrapolation=True)


def test_beattie_bridgeman_pressure_checks_computed_randall_range() -> None:
    gas = AIR_BEATTIE_BRIDGEMAN
    pressure_range = gas.applicable_pressure_range
    assert pressure_range is not None
    pressure_minimum, pressure_maximum = pressure_range

    def density_for_output_pressure(
        pressure: float, *, comparison: str, direction: float
    ) -> float:
        """Invert pressure, then choose the representable density on one side."""
        density = gas._density_scalar(900.0, pressure)
        for _ in range(32):
            output_pressure = gas._pressure_scalar(900.0, density)
            if (
                comparison == "at_or_above_minimum"
                and output_pressure >= pressure_minimum
            ):
                return density
            if (
                comparison == "at_or_below_maximum"
                and output_pressure <= pressure_maximum
            ):
                return density
            if comparison == "below_minimum" and output_pressure < pressure_minimum:
                return density
            if comparison == "above_maximum" and output_pressure > pressure_maximum:
                return density
            density = np.nextafter(density, direction)
        raise AssertionError("could not construct a pressure output at the boundary")

    for pressure, comparison, direction in (
        (pressure_minimum, "at_or_above_minimum", np.inf),
        (pressure_maximum, "at_or_below_maximum", 0.0),
    ):
        density = density_for_output_pressure(
            pressure, comparison=comparison, direction=direction
        )
        assert gas.pressure(900.0, density) == pytest.approx(pressure)

    for pressure, comparison, direction in (
        (np.nextafter(pressure_minimum, 0.0), "below_minimum", 0.0),
        (np.nextafter(pressure_maximum, np.inf), "above_maximum", np.inf),
    ):
        density = density_for_output_pressure(
            pressure, comparison=comparison, direction=direction
        )
        with pytest.raises(ModelRangeError, match="pressure"):
            gas.pressure(900.0, density, allow_extrapolation=False)
        with pytest.warns(
            ApplicabilityWarning, match=r"Randall, AEDC-TR-57-8.*tabulated range"
        ) as captured:
            assert gas.pressure(
                900.0, density, allow_extrapolation=True
            ) == pytest.approx(pressure)
        assert len(captured) == 1


@pytest.mark.parametrize("density", [1.0, 100.0])
def test_beattie_bridgeman_pressure_combines_temperature_and_pressure_range(
    density: float,
) -> None:
    gas = AIR_BEATTIE_BRIDGEMAN
    with pytest.warns(
        ApplicabilityWarning, match=r"Randall, AEDC-TR-57-8.*tabulated range"
    ) as captured:
        pressure = gas.pressure(1300.0, density, allow_extrapolation=True)
    assert pressure > 0.0
    assert len(captured) == 1
    with pytest.raises(ModelRangeError):
        gas.pressure(1300.0, density, allow_extrapolation=False)


def test_beattie_bridgeman_exact_total_temperature_endpoint() -> None:
    state = isentropic_state(
        1.0e-10,
        AIR_BEATTIE_BRIDGEMAN,
        total_temperature=45.0,
        total_pressure=1.0e5,
        allow_extrapolation=False,
    )
    total = AIR_BEATTIE_BRIDGEMAN.state(45.0, 1.0e5)
    static = AIR_BEATTIE_BRIDGEMAN._scalar_state_from_density(
        float(state.static_temperature), float(state.static_density)
    )
    assert state.static_temperature == pytest.approx(45.0, abs=2e-12)
    assert state.static_pressure > 0.0
    assert state.static_density > 0.0
    assert static.entropy == pytest.approx(total.entropy, abs=2e-8)
    assert static.enthalpy + 0.5 * state.velocity**2 == pytest.approx(
        total.enthalpy, rel=2e-13
    )


@pytest.mark.parametrize("mach", [0.0, 1.0e-10, 1.0e-6])
def test_beattie_bridgeman_low_mach_log_endpoint_regression(mach: float) -> None:
    total_temperature = 565.0
    total_pressure = 1.0e5
    state = isentropic_state(
        mach,
        AIR_BEATTIE_BRIDGEMAN,
        total_temperature=total_temperature,
        total_pressure=total_pressure,
        allow_extrapolation=False,
    )
    assert 0.0 < state.static_temperature <= total_temperature
    assert state.static_pressure > 0.0
    assert state.static_density > 0.0
    if mach > 0.0:
        for value in (
            isentropic_ratios(
                mach,
                AIR_BEATTIE_BRIDGEMAN,
                total_temperature=total_temperature,
                total_pressure=total_pressure,
            ).total_temperature_ratio,
            mass_flow_parameter(
                mach,
                AIR_BEATTIE_BRIDGEMAN,
                total_temperature=total_temperature,
                total_pressure=total_pressure,
            ),
            mass_flux(total_pressure, total_temperature, mach, AIR_BEATTIE_BRIDGEMAN),
            area_ratio(
                mach,
                AIR_BEATTIE_BRIDGEMAN,
                total_temperature=total_temperature,
                total_pressure=total_pressure,
            ),
        ):
            assert np.isfinite(value)


def test_beattie_bridgeman_area_checks_critical_applicability() -> None:
    with pytest.raises(ModelRangeError, match="temperature"):
        area_ratio(
            0.1,
            AIR_BEATTIE_BRIDGEMAN,
            total_temperature=45.0,
            total_pressure=1.0e5,
            allow_extrapolation=False,
        )
    with pytest.warns(
        ApplicabilityWarning, match=r"Randall, AEDC-TR-57-8.*tabulated range"
    ) as captured:
        target = area_ratio(
            0.1,
            AIR_BEATTIE_BRIDGEMAN,
            total_temperature=45.0,
            total_pressure=1.0e5,
        )
    assert len(captured) == 1
    with pytest.raises(ModelRangeError, match="temperature"):
        mach_from_area_ratio(
            target,
            MachBranch.SUBSONIC,
            AIR_BEATTIE_BRIDGEMAN,
            total_temperature=45.0,
            total_pressure=1.0e5,
            allow_extrapolation=False,
        )
    with pytest.warns(
        ApplicabilityWarning, match=r"Randall, AEDC-TR-57-8.*tabulated range"
    ) as captured:
        mach = mach_from_area_ratio(
            target,
            MachBranch.SUBSONIC,
            AIR_BEATTIE_BRIDGEMAN,
            total_temperature=45.0,
            total_pressure=1.0e5,
        )
    assert len(captured) == 1
    assert mach == pytest.approx(0.1, abs=2e-11)


def test_custom_beattie_bridgeman_warning_does_not_claim_randall_air() -> None:
    gas = BeattieBridgemanGas(
        287.0,
        1.4,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        applicable_temperature_range=(100.0, 200.0),
        applicable_pressure_range=(1.0e3, 2.0e3),
    )
    with pytest.warns(
        ApplicabilityWarning, match=r"configured.*tabulated range"
    ) as captured:
        gas.state(300.0, 3.0e3, allow_extrapolation=True)
    assert len(captured) == 1
    assert "Randall" not in str(captured[0].message)

    with pytest.warns(
        ApplicabilityWarning, match=r"configured.*tabulated range"
    ) as captured:
        isentropic_ratios(
            1.0,
            gas,
            total_temperature=300.0,
            total_pressure=3.0e3,
            allow_extrapolation=True,
        )
    assert len(captured) == 1
    assert "Randall" not in str(captured[0].message)
    with pytest.raises(ModelRangeError):
        isentropic_ratios(
            1.0,
            gas,
            total_temperature=300.0,
            total_pressure=3.0e3,
            allow_extrapolation=False,
        )


@pytest.mark.parametrize(
    ("coefficients", "mach", "total_temperature"),
    [
        (
            (
                12825.936657414117,
                2.0739052512358574e-06,
                0.21038314102297515,
                -0.000816567523215019,
                28906.11424191098,
            ),
            32.0,
            221.48207105566667,
        ),
        (
            (
                4918.005459004671,
                6.453841220102751e-06,
                0.12006443547655751,
                -0.006863118307642916,
                25215.558965841552,
            ),
            20.344423140819163,
            221.48207105566667,
        ),
        (
            (
                4918.005459004671,
                6.453841220102751e-06,
                0.12006443547655751,
                -0.006863118307642916,
                25215.558965841552 * 0.9205,
            ),
            20.0,
            221.49,
        ),
    ],
)
def test_real_isentropic_solver_rejects_disconnected_validity_gap(
    coefficients: tuple[float, float, float, float, float],
    mach: float,
    total_temperature: float,
) -> None:
    disconnected = BeattieBridgemanGas(287.05287, 1.4, *coefficients)
    with pytest.raises(ModelRangeError, match=r"no connected.*solution"):
        _real_flow_state(mach, total_temperature, 55.716069666e6, disconnected)


def test_zero_beattie_bridgeman_correction_is_harmonic_oscillator_limit() -> None:
    harmonic, real = _zero_correction_pair()
    temperature = np.array([600.0, 1200.0, 1800.0])
    pressure = np.array([1.0e6, 5.0e6, 10.0e6])
    ideal_state = harmonic.state(temperature, pressure)
    real_state = real.state(temperature, pressure)
    for name in (
        "density",
        "internal_energy",
        "enthalpy",
        "cp",
        "cv",
        "heat_capacity_ratio",
        "speed_of_sound",
    ):
        assert_allclose(
            getattr(real_state, name), getattr(ideal_state, name), rtol=2e-14
        )


@pytest.mark.parametrize(
    "gas,total_pressure",
    [
        (AIR_HARMONIC_OSCILLATOR, None),
        (AIR_BEATTIE_BRIDGEMAN, 6.0e6),
    ],
)
def test_new_models_isentropic_forward_inverse_and_area_branches(
    gas: HarmonicOscillatorGas | BeattieBridgemanGas,
    total_pressure: float | None,
) -> None:
    mach = np.array([0.0, 0.4, 1.0, 2.0])
    ratios = isentropic_ratios(
        mach,
        gas,
        total_temperature=1200.0,
        total_pressure=total_pressure,
    )
    assert_allclose(
        mach_from_total_temperature_ratio(
            ratios.total_temperature_ratio,
            gas,
            total_temperature=1200.0,
            total_pressure=total_pressure,
        ),
        mach,
        atol=3e-11,
    )
    assert_allclose(
        mach_from_total_pressure_ratio(
            ratios.total_pressure_ratio,
            gas,
            total_temperature=1200.0,
            total_pressure=total_pressure,
        ),
        mach,
        atol=3e-11,
    )
    assert_allclose(
        mach_from_total_density_ratio(
            ratios.total_density_ratio,
            gas,
            total_temperature=1200.0,
            total_pressure=total_pressure,
        ),
        mach,
        atol=3e-11,
    )
    critical = critical_ratios(
        gas, total_temperature=1200.0, total_pressure=total_pressure
    )
    at_one = isentropic_ratios(
        1.0,
        gas,
        total_temperature=1200.0,
        total_pressure=total_pressure,
    )
    assert critical.total_pressure_ratio == pytest.approx(at_one.total_pressure_ratio)
    assert area_ratio(
        1.0,
        gas,
        total_temperature=1200.0,
        total_pressure=total_pressure,
    ) == pytest.approx(1.0)
    target = area_ratio(
        2.0,
        gas,
        total_temperature=1200.0,
        total_pressure=total_pressure,
    )
    assert mach_from_area_ratio(
        target,
        MachBranch.SUPERSONIC,
        gas,
        total_temperature=1200.0,
        total_pressure=total_pressure,
    ) == pytest.approx(2.0, abs=2e-11)
    assert (
        0.0
        < mach_from_area_ratio(
            target,
            MachBranch.SUBSONIC,
            gas,
            total_temperature=1200.0,
            total_pressure=total_pressure,
        )
        < 1.0
    )


@pytest.mark.parametrize(
    "gas",
    [
        AIR_HARMONIC_OSCILLATOR,
        AIR_BEATTIE_BRIDGEMAN,
    ],
)
def test_isentropic_state_conserves_energy_and_entropy(
    gas: HarmonicOscillatorGas | BeattieBridgemanGas,
) -> None:
    total_temperature = 1200.0
    total_pressure = 6.0e6
    state = isentropic_state(
        np.array([0.0, 1.0, 2.0]),
        gas,
        total_temperature=total_temperature,
        total_pressure=total_pressure,
    )
    assert isinstance(state.static_temperature, np.ndarray)
    assert_allclose(state.velocity, np.asarray(state.mach) * state.speed_of_sound)
    assert_allclose(state.mass_flux, state.static_density * state.velocity)
    assert_allclose(
        state.dynamic_pressure,
        0.5 * state.static_density * np.asarray(state.velocity) ** 2,
    )
    total = gas.state(total_temperature, total_pressure)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ApplicabilityWarning)
        static = gas.state(
            state.static_temperature,
            state.static_pressure,
            allow_extrapolation=True,
        )
    assert_allclose(static.entropy, total.entropy, atol=2e-8)
    assert_allclose(
        static.enthalpy + 0.5 * np.asarray(state.velocity) ** 2,
        total.enthalpy,
        rtol=2e-13,
    )


def test_beattie_bridgeman_high_mach_gas_branch_conserves_state() -> None:
    gas = AIR_BEATTIE_BRIDGEMAN
    total_temperature = 900.0
    total_pressure = 1.0e7
    mach = np.arange(11.0)
    state = isentropic_state(
        mach,
        gas,
        total_temperature=total_temperature,
        total_pressure=total_pressure,
        allow_extrapolation=False,
    )
    total = gas.state(total_temperature, total_pressure)
    static_temperatures = np.asarray(state.static_temperature, dtype=np.float64)
    static_densities = np.asarray(state.static_density, dtype=np.float64)
    velocities = np.asarray(state.velocity, dtype=np.float64)
    scalar_states = [
        gas._scalar_state_from_density(float(temperature), float(density))
        for temperature, density in zip(
            static_temperatures, static_densities, strict=True
        )
    ]
    assert np.all(np.asarray(state.static_pressure) > 0.0)
    assert np.all(np.asarray(state.static_density) > 0.0)
    assert all(
        gas._dp_drho_scalar(item.temperature, item.density) > 0.0
        for item in scalar_states
    )
    assert all(item.cp > 0.0 and item.cv > 0.0 for item in scalar_states)
    assert all(item.speed_of_sound**2 > 0.0 for item in scalar_states)
    assert_allclose([item.entropy for item in scalar_states], total.entropy, atol=2e-8)
    assert_allclose(
        [
            item.enthalpy + 0.5 * velocity**2
            for item, velocity in zip(scalar_states, velocities, strict=True)
        ],
        total.enthalpy,
        rtol=2e-13,
    )


def test_beattie_bridgeman_continuation_reaches_very_low_temperature() -> None:
    _, gas = _zero_correction_pair()
    total_temperature = 1000.0
    total_pressure = 1.0e5
    state = isentropic_state(
        50.0,
        gas,
        total_temperature=total_temperature,
        total_pressure=total_pressure,
    )
    static = gas._scalar_state_from_density(
        float(state.static_temperature), float(state.static_density)
    )
    total = gas.state(total_temperature, total_pressure)
    assert state.static_temperature < total_temperature * np.exp(-0.8)
    assert static.pressure > 0.0
    assert static.density > 0.0
    assert gas._dp_drho_scalar(static.temperature, static.density) > 0.0
    assert static.cp > 0.0
    assert static.cv > 0.0
    assert static.entropy == pytest.approx(total.entropy, abs=2e-8)
    assert static.enthalpy + 0.5 * state.velocity**2 == pytest.approx(
        total.enthalpy, rel=2e-13
    )


def test_beattie_bridgeman_mach_ten_regression_and_static_range_handling() -> None:
    state = isentropic_state(
        10.0,
        AIR_BEATTIE_BRIDGEMAN,
        total_temperature=900.0,
        total_pressure=1.0e7,
        allow_extrapolation=False,
    )
    static = AIR_BEATTIE_BRIDGEMAN._scalar_state_from_density(
        float(state.static_temperature), float(state.static_density)
    )
    assert np.isfinite(state.static_temperature)
    assert state.static_pressure > 0.0
    assert state.static_density > 0.0
    assert static.speed_of_sound > 0.0
    assert (
        AIR_BEATTIE_BRIDGEMAN._dp_drho_scalar(static.temperature, static.density) > 0.0
    )

    with pytest.raises(ModelRangeError, match="temperature"):
        isentropic_state(
            11.0,
            AIR_BEATTIE_BRIDGEMAN,
            total_temperature=900.0,
            total_pressure=1.0e7,
            allow_extrapolation=False,
        )
    with pytest.warns(
        ApplicabilityWarning, match=r"Randall, AEDC-TR-57-8.*tabulated range"
    ):
        isentropic_state(
            11.0,
            AIR_BEATTIE_BRIDGEMAN,
            total_temperature=900.0,
            total_pressure=1.0e7,
            allow_extrapolation=True,
        )


@pytest.mark.parametrize(
    "gas,total_pressure",
    [
        (AIR_HARMONIC_OSCILLATOR, None),
        (AIR_BEATTIE_BRIDGEMAN, 6.0e6),
    ],
)
def test_mass_flux_chokes_at_mach_one(
    gas: HarmonicOscillatorGas | BeattieBridgemanGas,
    total_pressure: float | None,
) -> None:
    mach = np.linspace(0.05, 2.0, 80)
    parameter = mass_flow_parameter(
        mach,
        gas,
        total_temperature=1200.0,
        total_pressure=total_pressure,
    )
    assert int(np.argmax(parameter)) in {38, 39}
    flux = mass_flux(6.0e6, 1200.0, 1.0, gas)
    assert flux == pytest.approx(choked_mass_flux(6.0e6, 1200.0, gas))


def test_isentropic_requirements_broadcast_and_range_handling() -> None:
    with pytest.raises(ValueError, match="total_temperature"):
        isentropic_ratios(1.0, AIR_HARMONIC_OSCILLATOR)
    with pytest.raises(ValueError, match="total_pressure"):
        isentropic_ratios(1.0, AIR_BEATTIE_BRIDGEMAN, total_temperature=1200.0)
    with pytest.warns(ApplicabilityWarning) as captured:
        isentropic_ratios(
            [1.0, 2.0],
            AIR_BEATTIE_BRIDGEMAN,
            total_temperature=[300.0, 2200.0],
            total_pressure=20.0e6,
        )
    assert len(captured) == 1
    with pytest.raises(ModelRangeError, match="temperature"):
        isentropic_ratios(
            1.0,
            AIR_BEATTIE_BRIDGEMAN,
            total_temperature=30.0,
            total_pressure=6.0e6,
            allow_extrapolation=False,
        )
    with warnings.catch_warnings(record=True) as captured_harmonic:
        warnings.simplefilter("always")
        isentropic_ratios(
            [1.0, 2.0],
            AIR_HARMONIC_OSCILLATOR,
            total_temperature=300.0,
        )
    assert len(captured_harmonic) == 1
    with pytest.raises(ModelRangeError, match="total_temperature"):
        critical_ratios(
            AIR_HARMONIC_OSCILLATOR,
            total_temperature=300.0,
            allow_extrapolation=False,
        )
    ratios = isentropic_ratios(
        np.array([[0.5], [1.0]]),
        AIR_BEATTIE_BRIDGEMAN,
        total_temperature=np.array([1000.0, 1200.0]),
        total_pressure=6.0e6,
    )
    assert np.shape(ratios.total_pressure_ratio) == (2, 2)


def test_unstable_beattie_bridgeman_state_fails_explicitly() -> None:
    pathological = BeattieBridgemanGas(
        287.0,
        1.4,
        a0=1.0e8,
        b0=0.0,
        a=0.0,
        b=0.0,
        c=0.0,
    )
    with pytest.raises(ModelRangeError, match="gas-phase density root"):
        pathological.density(300.0, 1.0e6)
