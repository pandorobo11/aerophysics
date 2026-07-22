"""Tests for perfect-gas isentropic relations."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics.gas import PerfectGas
from aerophysics.isentropic import (
    MachBranch,
    area_ratio,
    choked_mass_flux,
    critical_ratios,
    isentropic_ratios,
    mach_from_area_ratio,
    mach_from_total_density_ratio,
    mach_from_total_pressure_ratio,
    mach_from_total_temperature_ratio,
    mass_flow_parameter,
    mass_flux,
)


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
