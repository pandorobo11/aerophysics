"""Tests for Prandtl-Meyer expansion relations."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics.expansion import (
    mach_from_prandtl_meyer,
    maximum_prandtl_meyer_angle,
    prandtl_meyer_angle,
    prandtl_meyer_expansion,
)
from aerophysics.gas import PerfectGas
from aerophysics.units import degrees_to_radians, radians_to_degrees


def test_prandtl_meyer_reference_values() -> None:
    assert float(radians_to_degrees(prandtl_meyer_angle(1.0))) == 0.0
    assert float(radians_to_degrees(prandtl_meyer_angle(2.0))) == pytest.approx(
        26.3797608, rel=1e-9
    )
    assert float(radians_to_degrees(prandtl_meyer_angle(3.0))) == pytest.approx(
        49.7573467, rel=1e-9
    )


def test_limiting_angle() -> None:
    assert float(radians_to_degrees(maximum_prandtl_meyer_angle())) == pytest.approx(
        130.4540769, rel=1e-9
    )


def test_inverse_round_trip_vectorizes_and_expands_bracket() -> None:
    mach = np.array([[1.0, 2.0], [5.0, 20.0]])
    angle = prandtl_meyer_angle(mach)
    result = mach_from_prandtl_meyer(angle)
    assert isinstance(result, np.ndarray)
    assert result.shape == (2, 2)
    assert result.dtype == np.float64
    assert_allclose(result, mach, rtol=1e-11, atol=1e-11)


def test_ten_degree_expansion_from_mach_two() -> None:
    result = prandtl_meyer_expansion(2.0, float(degrees_to_radians(10.0)))
    assert result.downstream_mach == pytest.approx(2.38488715, rel=1e-8)
    assert result.static_temperature_ratio == pytest.approx(0.84209055, rel=1e-8)
    assert result.static_pressure_ratio == pytest.approx(0.54796873, rel=1e-8)
    assert result.static_density_ratio == pytest.approx(0.65072424, rel=1e-8)
    assert result.downstream_prandtl_meyer_angle == pytest.approx(
        float(result.upstream_prandtl_meyer_angle) + float(result.turn_angle)
    )


def test_zero_turn_is_identity() -> None:
    result = prandtl_meyer_expansion(3.0, 0.0)
    assert result.downstream_mach == pytest.approx(3.0)
    assert result.static_temperature_ratio == pytest.approx(1.0)
    assert result.static_pressure_ratio == pytest.approx(1.0)
    assert result.static_density_ratio == pytest.approx(1.0)


def test_expansion_broadcasts() -> None:
    result = prandtl_meyer_expansion(
        [[2.0], [3.0]], degrees_to_radians([0.0, 5.0, 10.0])
    )
    for value in (
        result.upstream_mach,
        result.downstream_mach,
        result.turn_angle,
        result.upstream_prandtl_meyer_angle,
        result.downstream_prandtl_meyer_angle,
        result.static_temperature_ratio,
        result.static_pressure_ratio,
        result.static_density_ratio,
    ):
        assert isinstance(value, np.ndarray)
        assert value.shape == (2, 3)
        assert value.dtype == np.float64
    assert np.all(
        np.asarray(result.downstream_mach) >= np.asarray(result.upstream_mach)
    )


def test_custom_gas() -> None:
    helium = PerfectGas(2_077.1, 5.0 / 3.0)
    angle = prandtl_meyer_angle(2.0, helium)
    assert mach_from_prandtl_meyer(angle, helium) == pytest.approx(2.0)
    assert maximum_prandtl_meyer_angle(helium) < maximum_prandtl_meyer_angle()


@pytest.mark.parametrize("mach", [0.9, np.nan, np.inf])
def test_relations_reject_invalid_mach(mach: float) -> None:
    with pytest.raises(ValueError):
        prandtl_meyer_angle(mach)
    with pytest.raises(ValueError):
        prandtl_meyer_expansion(mach, 0.1)


def test_inverse_rejects_invalid_angles() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        mach_from_prandtl_meyer(-0.1)
    with pytest.raises(ValueError, match="limiting"):
        mach_from_prandtl_meyer(maximum_prandtl_meyer_angle())
    with pytest.raises(ValueError, match="too close"):
        mach_from_prandtl_meyer(np.nextafter(maximum_prandtl_meyer_angle(), 0.0))
    with pytest.raises(ValueError):
        mach_from_prandtl_meyer(np.nan)


def test_expansion_rejects_invalid_turn_or_shape() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        prandtl_meyer_expansion(2.0, -0.1)
    with pytest.raises(ValueError, match="limiting"):
        prandtl_meyer_expansion(2.0, maximum_prandtl_meyer_angle())
    with pytest.raises(ValueError, match="broadcastable"):
        prandtl_meyer_expansion([2.0, 3.0], [0.1, 0.2, 0.3])
