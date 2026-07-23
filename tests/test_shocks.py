"""Tests for normal and oblique perfect-gas shocks."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics.exceptions import NoAttachedShockError
from aerophysics.shocks import (
    ShockBranch,
    maximum_attached_deflection,
    normal_shock,
    oblique_shock,
    shock_angle,
    supersonic_pitot_pressure_ratio,
    theta_from_shock_angle,
)
from aerophysics.units import degrees_to_radians, radians_to_degrees


def test_normal_shock_mach_two_reference_values() -> None:
    result = normal_shock(2.0)
    assert result.downstream_mach == pytest.approx(0.5773502692)
    assert result.static_pressure_ratio == pytest.approx(4.5)
    assert result.static_density_ratio == pytest.approx(8.0 / 3.0)
    assert result.static_temperature_ratio == pytest.approx(1.6875)
    assert result.total_pressure_ratio == pytest.approx(0.7208738615)


def test_normal_shock_at_mach_one_is_degenerate() -> None:
    result = normal_shock(1.0)
    assert result.downstream_mach == 1.0
    assert result.static_pressure_ratio == 1.0
    assert result.static_density_ratio == 1.0
    assert result.static_temperature_ratio == 1.0
    assert result.total_pressure_ratio == 1.0


def test_normal_shock_vectorizes() -> None:
    result = normal_shock([[1.0, 2.0], [3.0, 5.0]])
    for value in (
        result.upstream_mach,
        result.downstream_mach,
        result.static_pressure_ratio,
        result.static_density_ratio,
        result.static_temperature_ratio,
        result.total_pressure_ratio,
    ):
        assert isinstance(value, np.ndarray)
        assert value.shape == (2, 2)
        assert value.dtype == np.float64
    assert np.all(np.asarray(result.total_pressure_ratio) <= 1.0)


def test_theta_beta_mach_weak_reference_value() -> None:
    mach = 2.5
    theta = float(degrees_to_radians(15.0))
    beta = shock_angle(mach, theta, ShockBranch.WEAK)
    assert float(radians_to_degrees(beta)) == pytest.approx(36.9449, rel=1e-5)
    assert theta_from_shock_angle(mach, beta) == pytest.approx(theta, rel=1e-12)


def test_weak_and_strong_oblique_shocks() -> None:
    theta = float(degrees_to_radians(15.0))
    weak = oblique_shock(2.5, theta, ShockBranch.WEAK)
    strong = oblique_shock(2.5, theta, ShockBranch.STRONG)
    assert float(weak.shock_angle) < float(strong.shock_angle)
    assert float(weak.downstream_mach) > 1.0
    assert float(strong.downstream_mach) < 1.0
    assert float(weak.static_pressure_ratio) < float(strong.static_pressure_ratio)
    assert float(weak.total_pressure_ratio) > float(strong.total_pressure_ratio)


def test_zero_deflection_branch_limits() -> None:
    weak = oblique_shock(2.0, 0.0, ShockBranch.WEAK)
    strong = oblique_shock(2.0, 0.0, ShockBranch.STRONG)
    assert weak.shock_angle == pytest.approx(np.arcsin(0.5))
    assert weak.downstream_mach == pytest.approx(2.0)
    assert weak.static_pressure_ratio == pytest.approx(1.0)
    assert strong.shock_angle == pytest.approx(0.5 * np.pi)
    assert strong.downstream_mach == pytest.approx(normal_shock(2.0).downstream_mach)


def test_maximum_attached_deflection_is_turning_point() -> None:
    limit = maximum_attached_deflection(2.0)
    assert float(radians_to_degrees(limit.deflection_angle)) == pytest.approx(
        22.9735, rel=1e-5
    )
    weak = shock_angle(2.0, limit.deflection_angle, ShockBranch.WEAK)
    strong = shock_angle(2.0, limit.deflection_angle, ShockBranch.STRONG)
    assert weak == pytest.approx(strong, rel=1e-10)
    assert theta_from_shock_angle(2.0, limit.shock_angle) == pytest.approx(
        limit.deflection_angle
    )


def test_mach_one_attached_limit() -> None:
    limit = maximum_attached_deflection(1.0)
    assert limit.deflection_angle == 0.0
    assert limit.shock_angle == pytest.approx(0.5 * np.pi)


def test_oblique_shock_broadcasts() -> None:
    theta = degrees_to_radians([5.0, 10.0, 15.0])
    result = oblique_shock([[2.0], [3.0]], theta)
    assert isinstance(result.shock_angle, np.ndarray)
    assert result.shock_angle.shape == (2, 3)
    assert_allclose(
        theta_from_shock_angle(result.upstream_mach, result.shock_angle),
        np.broadcast_to(theta, (2, 3)),
    )


def test_supersonic_pitot_pressure_ratio() -> None:
    assert supersonic_pitot_pressure_ratio(1.0) == pytest.approx(1.8929291587)
    assert supersonic_pitot_pressure_ratio(2.0) == pytest.approx(5.6404408128)
    result = supersonic_pitot_pressure_ratio([1.0, 2.0, 3.0])
    assert isinstance(result, np.ndarray)
    assert np.all(np.diff(result) > 0.0)


@pytest.mark.parametrize("mach", [0.9, np.nan, np.inf])
def test_shock_relations_reject_invalid_mach(mach: float) -> None:
    with pytest.raises(ValueError):
        normal_shock(mach)
    with pytest.raises(ValueError):
        maximum_attached_deflection(mach)


def test_theta_beta_mach_rejects_invalid_inputs() -> None:
    with pytest.raises(ValueError, match="Mach angle"):
        theta_from_shock_angle(2.0, 0.1)
    with pytest.raises(ValueError, match="Mach angle"):
        theta_from_shock_angle(2.0, 2.0)
    with pytest.raises(ValueError, match="non-negative"):
        shock_angle(2.0, -0.1, ShockBranch.WEAK)
    with pytest.raises(ValueError, match="ShockBranch"):
        shock_angle(2.0, 0.1, "weak")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="broadcastable"):
        shock_angle([2.0, 3.0], [0.1, 0.2, 0.3], ShockBranch.WEAK)


def test_detached_shock_raises_dedicated_error() -> None:
    with pytest.raises(NoAttachedShockError):
        oblique_shock(2.0, float(degrees_to_radians(30.0)))
