"""Tests for perfect-gas isentropic relations."""

import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics import AIR_NASA7, AIR_NASA9
from aerophysics.exceptions import ApplicabilityWarning, ModelRangeError
from aerophysics.gas import PerfectGas
from aerophysics.isentropic import (
    MachBranch,
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
from aerophysics.thermochemistry import (
    UNIVERSAL_GAS_CONSTANT,
    IdealGasSpecies,
    NASA7Polynomial,
    ThermallyPerfectGas,
)


def _constant_cp_gas() -> tuple[PerfectGas, ThermallyPerfectGas]:
    gas_constant = 287.0
    caloric = PerfectGas(gas_constant, 1.4)
    polynomial = NASA7Polynomial(
        (100.0, 10_000.0),
        ((3.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),),
    )
    species = IdealGasSpecies(
        "constant-cp",
        UNIVERSAL_GAS_CONSTANT / gas_constant,
        polynomial,
    )
    return caloric, ThermallyPerfectGas((species,), (1.0,))


def test_isentropic_reference_values() -> None:
    ratios = isentropic_ratios(2.0)
    assert ratios.mach == 2.0
    assert ratios.total_temperature_ratio == pytest.approx(1.8)
    assert ratios.total_pressure_ratio == pytest.approx(7.8244490669)
    assert ratios.total_density_ratio == pytest.approx(4.3469161483)


def test_zero_mach_has_unity_state_ratios() -> None:
    ratios = isentropic_ratios(0.0)
    assert ratios.total_temperature_ratio == 1.0
    assert ratios.total_pressure_ratio == 1.0
    assert ratios.total_density_ratio == 1.0
    assert mass_flow_parameter(0.0) == 0.0


def test_ratio_functions_vectorize_and_invert() -> None:
    mach = np.array([[0.0, 0.5], [1.0, 5.0]])
    ratios = isentropic_ratios(mach)
    assert isinstance(ratios.total_pressure_ratio, np.ndarray)
    assert ratios.total_pressure_ratio.shape == (2, 2)
    assert ratios.total_pressure_ratio.dtype == np.float64
    assert_allclose(
        mach_from_total_temperature_ratio(ratios.total_temperature_ratio), mach
    )
    assert_allclose(mach_from_total_pressure_ratio(ratios.total_pressure_ratio), mach)
    assert_allclose(mach_from_total_density_ratio(ratios.total_density_ratio), mach)


def test_area_ratio_reference_values() -> None:
    assert area_ratio(0.5) == pytest.approx(1.33984375)
    assert area_ratio(1.0) == 1.0
    assert area_ratio(2.0) == pytest.approx(1.6875)


def test_area_ratio_inverse_branches() -> None:
    assert mach_from_area_ratio(1.0, MachBranch.SUBSONIC) == 1.0
    assert mach_from_area_ratio(1.0, MachBranch.SUPERSONIC) == 1.0
    assert mach_from_area_ratio(1.6875, MachBranch.SUBSONIC) == pytest.approx(
        0.3722444862, rel=1e-10
    )
    assert mach_from_area_ratio(1.6875, MachBranch.SUPERSONIC) == pytest.approx(
        2.0, rel=1e-12
    )


def test_area_ratio_inverse_vectorizes_and_expands_bracket() -> None:
    targets = area_ratio([0.2, 1.0, 2.0, 5.0])
    assert isinstance(targets, np.ndarray)
    subsonic = mach_from_area_ratio(targets[:2], MachBranch.SUBSONIC)
    supersonic = mach_from_area_ratio(targets[1:], MachBranch.SUPERSONIC)
    assert_allclose(subsonic, [0.2, 1.0], atol=1e-11)
    assert_allclose(supersonic, [1.0, 2.0, 5.0], atol=1e-11)


def test_critical_ratios_match_mach_one() -> None:
    critical = critical_ratios()
    at_one = isentropic_ratios(1.0)
    assert critical.total_temperature_ratio == at_one.total_temperature_ratio
    assert critical.total_pressure_ratio == at_one.total_pressure_ratio
    assert critical.total_density_ratio == at_one.total_density_ratio


def test_custom_gas_is_supported() -> None:
    helium = PerfectGas(2_077.1, 5.0 / 3.0)
    ratios = isentropic_ratios(1.0, helium)
    assert ratios.total_temperature_ratio == pytest.approx(4.0 / 3.0)
    assert mach_from_total_pressure_ratio(ratios.total_pressure_ratio, helium) == (
        pytest.approx(1.0)
    )


def test_mass_flow_parameter_and_flux() -> None:
    parameter = mass_flow_parameter(1.0)
    assert parameter == pytest.approx(0.6847314564)
    flux = mass_flux(101_325.0, 288.15, 1.0)
    assert flux == pytest.approx(241.2384229, rel=1e-10)
    assert choked_mass_flux(101_325.0, 288.15) == flux
    assert mass_flow_parameter(1.0) > mass_flow_parameter(0.9)
    assert mass_flow_parameter(1.0) > mass_flow_parameter(1.1)


def test_mass_flux_broadcasts() -> None:
    result = mass_flux([[100_000.0], [200_000.0]], [250.0, 300.0], 1.0)
    assert isinstance(result, np.ndarray)
    assert result.shape == (2, 2)
    assert result.dtype == np.float64
    assert_allclose(result[1], 2.0 * result[0])


def test_absolute_isentropic_state_supports_caloric_and_nasa_air() -> None:
    caloric = isentropic_state(
        [0.0, 2.0],
        total_temperature=300.0,
        total_pressure=101_325.0,
    )
    assert_allclose(caloric.static_temperature, [300.0, 300.0 / 1.8])
    assert_allclose(caloric.velocity, np.asarray(caloric.mach) * caloric.speed_of_sound)
    assert_allclose(caloric.mass_flux, caloric.static_density * caloric.velocity)
    assert_allclose(
        caloric.dynamic_pressure,
        0.5 * caloric.static_density * np.asarray(caloric.velocity) ** 2,
    )

    thermal = isentropic_state(
        2.0,
        AIR_NASA9,
        total_temperature=1000.0,
        total_pressure=100_000.0,
        allow_extrapolation=False,
    )
    ratios = isentropic_ratios(
        2.0,
        AIR_NASA9,
        total_temperature=1000.0,
        allow_extrapolation=False,
    )
    assert thermal.static_temperature == pytest.approx(
        1000.0 / ratios.total_temperature_ratio
    )
    assert thermal.static_pressure == pytest.approx(
        100_000.0 / ratios.total_pressure_ratio
    )
    assert thermal.velocity == pytest.approx(2.0 * thermal.speed_of_sound)


def test_absolute_isentropic_state_validates_total_state() -> None:
    with pytest.raises(ValueError, match="total_temperature"):
        isentropic_state(1.0, total_temperature=0.0, total_pressure=1.0)
    with pytest.raises(ValueError, match="total_pressure"):
        isentropic_state(1.0, total_temperature=300.0, total_pressure=0.0)
    with pytest.raises(ValueError, match="broadcastable"):
        isentropic_state(
            [1.0, 2.0],
            total_temperature=[300.0, 400.0, 500.0],
            total_pressure=1.0,
        )


@pytest.mark.parametrize("mach", [-1.0, np.nan, np.inf])
def test_relations_reject_invalid_mach(mach: float) -> None:
    with pytest.raises(ValueError):
        isentropic_ratios(mach)
    with pytest.raises(ValueError):
        mass_flow_parameter(mach)


def test_area_ratio_rejects_zero_mach_and_invalid_ratio() -> None:
    with pytest.raises(ValueError):
        area_ratio(0.0)
    with pytest.raises(ValueError):
        mach_from_area_ratio(0.9, MachBranch.SUBSONIC)
    with pytest.raises(ValueError, match="MachBranch"):
        mach_from_area_ratio(2.0, "subsonic")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="too large"):
        mach_from_area_ratio(1e100, MachBranch.SUPERSONIC)


@pytest.mark.parametrize(
    "inverse",
    [
        mach_from_total_temperature_ratio,
        mach_from_total_pressure_ratio,
        mach_from_total_density_ratio,
    ],
)
def test_ratio_inverses_reject_values_below_one(inverse: object) -> None:
    assert callable(inverse)
    with pytest.raises(ValueError):
        inverse(0.99)


def test_mass_flux_rejects_invalid_or_incompatible_inputs() -> None:
    with pytest.raises(ValueError, match="total_pressure"):
        mass_flux(0.0, 300.0, 1.0)
    with pytest.raises(ValueError, match="total_temperature"):
        mass_flux(100_000.0, 0.0, 1.0)
    with pytest.raises(ValueError, match="broadcastable"):
        mass_flux([1.0, 2.0], [1.0, 2.0, 3.0], 1.0)


def test_constant_cp_thermal_model_matches_closed_form_relations() -> None:
    caloric, thermal = _constant_cp_gas()
    total_temperature = 1000.0
    mach = np.asarray([0.0, 0.5, 1.0, 2.0])

    expected = isentropic_ratios(mach, caloric)
    actual = isentropic_ratios(
        mach,
        thermal,
        total_temperature=total_temperature,
        allow_extrapolation=False,
    )
    assert_allclose(actual.total_temperature_ratio, expected.total_temperature_ratio)
    assert_allclose(actual.total_pressure_ratio, expected.total_pressure_ratio)
    assert_allclose(actual.total_density_ratio, expected.total_density_ratio)
    assert_allclose(
        mass_flow_parameter(
            mach,
            thermal,
            total_temperature=total_temperature,
            allow_extrapolation=False,
        ),
        mass_flow_parameter(mach, caloric),
    )
    assert_allclose(
        area_ratio(
            mach[1:],
            thermal,
            total_temperature=total_temperature,
            allow_extrapolation=False,
        ),
        area_ratio(mach[1:], caloric),
    )

    expected_critical = critical_ratios(caloric)
    actual_critical = critical_ratios(
        thermal,
        total_temperature=total_temperature,
        allow_extrapolation=False,
    )
    assert actual_critical.total_temperature_ratio == pytest.approx(
        expected_critical.total_temperature_ratio
    )
    assert actual_critical.total_pressure_ratio == pytest.approx(
        expected_critical.total_pressure_ratio
    )
    assert actual_critical.total_density_ratio == pytest.approx(
        expected_critical.total_density_ratio
    )
    assert mass_flux(
        200_000.0,
        total_temperature,
        2.0,
        thermal,
        allow_extrapolation=False,
    ) == pytest.approx(mass_flux(200_000.0, total_temperature, 2.0, caloric))


def test_nasa9_forward_relations_conserve_energy_and_entropy() -> None:
    total_temperature = 1000.0
    total_pressure = 100_000.0
    mach = 2.0
    result = isentropic_ratios(
        mach,
        AIR_NASA9,
        total_temperature=total_temperature,
        allow_extrapolation=False,
    )
    static_temperature = total_temperature / result.total_temperature_ratio
    static_pressure = total_pressure / result.total_pressure_ratio
    velocity = mach * AIR_NASA9.speed_of_sound(static_temperature)

    assert result.total_temperature_ratio == pytest.approx(1.7221397147)
    assert result.total_pressure_ratio == pytest.approx(7.8946725067)
    assert result.total_density_ratio == pytest.approx(4.5842230101)
    assert AIR_NASA9.standard_enthalpy(total_temperature) == pytest.approx(
        AIR_NASA9.standard_enthalpy(static_temperature) + 0.5 * velocity**2
    )
    assert AIR_NASA9.entropy(total_temperature, total_pressure) == pytest.approx(
        AIR_NASA9.entropy(static_temperature, static_pressure)
    )
    assert result.total_density_ratio == pytest.approx(
        result.total_pressure_ratio * static_temperature / total_temperature
    )


@pytest.mark.parametrize("gas", [AIR_NASA7, AIR_NASA9])
def test_thermal_relations_vectorize_and_all_inverses_round_trip(
    gas: ThermallyPerfectGas,
) -> None:
    mach = np.asarray([[0.0], [0.5], [1.0], [2.0]])
    total_temperature = np.asarray([[1000.0, 1500.0]])
    ratios = isentropic_ratios(
        mach,
        gas,
        total_temperature=total_temperature,
        allow_extrapolation=False,
    )
    expected_mach = np.broadcast_to(mach, (4, 2))
    assert np.asarray(ratios.total_pressure_ratio).shape == (4, 2)
    assert_allclose(
        mach_from_total_temperature_ratio(
            ratios.total_temperature_ratio,
            gas,
            total_temperature=total_temperature,
            allow_extrapolation=False,
        ),
        expected_mach,
        atol=2e-12,
    )
    assert_allclose(
        mach_from_total_pressure_ratio(
            ratios.total_pressure_ratio,
            gas,
            total_temperature=total_temperature,
            allow_extrapolation=False,
        ),
        expected_mach,
        atol=2e-12,
    )
    assert_allclose(
        mach_from_total_density_ratio(
            ratios.total_density_ratio,
            gas,
            total_temperature=total_temperature,
            allow_extrapolation=False,
        ),
        expected_mach,
        atol=2e-12,
    )


@pytest.mark.parametrize("gas", [AIR_NASA7, AIR_NASA9])
def test_thermal_area_branches_critical_state_and_choking(
    gas: ThermallyPerfectGas,
) -> None:
    total_temperature = 1000.0
    subsonic_ratio = area_ratio(
        0.4,
        gas,
        total_temperature=total_temperature,
        allow_extrapolation=False,
    )
    supersonic_ratio = area_ratio(
        2.0,
        gas,
        total_temperature=total_temperature,
        allow_extrapolation=False,
    )
    assert mach_from_area_ratio(
        subsonic_ratio,
        MachBranch.SUBSONIC,
        gas,
        total_temperature=total_temperature,
        allow_extrapolation=False,
    ) == pytest.approx(0.4)
    assert mach_from_area_ratio(
        supersonic_ratio,
        MachBranch.SUPERSONIC,
        gas,
        total_temperature=total_temperature,
        allow_extrapolation=False,
    ) == pytest.approx(2.0)

    critical = critical_ratios(
        gas,
        total_temperature=total_temperature,
        allow_extrapolation=False,
    )
    at_one = isentropic_ratios(
        1.0,
        gas,
        total_temperature=total_temperature,
        allow_extrapolation=False,
    )
    assert critical.total_temperature_ratio == at_one.total_temperature_ratio
    assert critical.total_pressure_ratio == at_one.total_pressure_ratio
    assert critical.total_density_ratio == at_one.total_density_ratio
    choked = choked_mass_flux(
        100_000.0,
        total_temperature,
        gas,
        allow_extrapolation=False,
    )
    assert choked > mass_flux(
        100_000.0,
        total_temperature,
        0.9,
        gas,
        allow_extrapolation=False,
    )
    assert choked > mass_flux(
        100_000.0,
        total_temperature,
        1.1,
        gas,
        allow_extrapolation=False,
    )


def test_thermal_critical_ratios_and_mass_flux_broadcast() -> None:
    temperatures = np.asarray([800.0, 1200.0])
    critical = critical_ratios(
        AIR_NASA9,
        total_temperature=temperatures,
        allow_extrapolation=False,
    )
    assert isinstance(critical.total_pressure_ratio, np.ndarray)
    assert critical.total_pressure_ratio.shape == (2,)
    flux = mass_flux(
        [[100_000.0], [200_000.0]],
        temperatures,
        1.0,
        AIR_NASA9,
        allow_extrapolation=False,
    )
    assert isinstance(flux, np.ndarray)
    assert flux.shape == (2, 2)
    assert_allclose(flux[1], 2.0 * flux[0])


def test_strict_supersonic_area_inverse_uses_available_temperature_range() -> None:
    target = area_ratio(
        1.2,
        AIR_NASA9,
        total_temperature=300.0,
        allow_extrapolation=False,
    )
    assert mach_from_area_ratio(
        target,
        MachBranch.SUPERSONIC,
        AIR_NASA9,
        total_temperature=300.0,
        allow_extrapolation=False,
    ) == pytest.approx(1.2)


@pytest.mark.parametrize("gas", [AIR_NASA7, AIR_NASA9])
@pytest.mark.parametrize("mach", [1.01, 1.1, 1.2, 1.5])
def test_thermal_area_inverse_ignores_extrapolated_solver_probes(
    gas: ThermallyPerfectGas,
    mach: float,
) -> None:
    target = area_ratio(
        mach,
        gas,
        total_temperature=300.0,
        allow_extrapolation=False,
    )

    with warnings.catch_warnings():
        warnings.simplefilter("error", ApplicabilityWarning)
        result = mach_from_area_ratio(
            target,
            MachBranch.SUPERSONIC,
            gas,
            total_temperature=300.0,
        )

    assert result == pytest.approx(mach)


@pytest.mark.parametrize("gas", [AIR_NASA7, AIR_NASA9])
def test_thermal_area_inverse_reports_converged_out_of_range_state(
    gas: ThermallyPerfectGas,
) -> None:
    with pytest.warns(ApplicabilityWarning):
        target = area_ratio(2.0, gas, total_temperature=300.0)

    with pytest.warns(ApplicabilityWarning) as caught:
        result = mach_from_area_ratio(
            target,
            MachBranch.SUPERSONIC,
            gas,
            total_temperature=300.0,
        )

    assert len(caught) == 1
    assert result == pytest.approx(2.0)

    with pytest.raises(ModelRangeError, match="below the fitted range"):
        mach_from_area_ratio(
            target,
            MachBranch.SUPERSONIC,
            gas,
            total_temperature=300.0,
            allow_extrapolation=False,
        )


def test_thermal_default_extrapolation_warns_once_and_strict_mode_rejects() -> None:
    with pytest.warns(ApplicabilityWarning) as caught:
        ratios = isentropic_ratios(
            [2.0, 3.0],
            AIR_NASA9,
            total_temperature=300.0,
        )
    assert len(caught) == 1
    assert np.all(np.isfinite(ratios.total_pressure_ratio))

    with pytest.raises(ModelRangeError, match="below the fitted range"):
        isentropic_ratios(
            2.0,
            AIR_NASA9,
            total_temperature=300.0,
            allow_extrapolation=False,
        )


@pytest.mark.parametrize(
    "function,arguments",
    [
        (isentropic_ratios, (1.0, AIR_NASA9)),
        (mach_from_total_temperature_ratio, (1.2, AIR_NASA9)),
        (mach_from_total_pressure_ratio, (2.0, AIR_NASA9)),
        (mach_from_total_density_ratio, (2.0, AIR_NASA9)),
        (area_ratio, (1.0, AIR_NASA9)),
        (mass_flow_parameter, (1.0, AIR_NASA9)),
        (critical_ratios, (AIR_NASA9,)),
    ],
)
def test_thermal_relations_require_total_temperature(
    function: object, arguments: tuple[object, ...]
) -> None:
    assert callable(function)
    with pytest.raises(ValueError, match="total_temperature is required"):
        function(*arguments)


def test_thermal_relations_reject_incompatible_and_nonphysical_models() -> None:
    with pytest.raises(ValueError, match="broadcastable"):
        isentropic_ratios(
            [1.0, 2.0],
            AIR_NASA9,
            total_temperature=[500.0, 600.0, 700.0],
        )

    polynomial = NASA7Polynomial(
        (100.0, 1000.0),
        ((0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),),
    )
    species = IdealGasSpecies("invalid-cv", 0.03, polynomial)
    nonphysical = ThermallyPerfectGas((species,), (1.0,))
    with pytest.raises(ModelRangeError, match="heat capacity"):
        isentropic_ratios(
            1.0,
            nonphysical,
            total_temperature=500.0,
            allow_extrapolation=False,
        )
