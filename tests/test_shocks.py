"""Tests for normal and oblique perfect-gas shocks."""

import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics.exceptions import NoAttachedShockError
from aerophysics.gas import PerfectGas
from aerophysics.shocks import (
    ShockBranch,
    conical_shock,
    maximum_attached_cone_angle,
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


def test_conical_shock_matches_nasa_sp_3004() -> None:
    result = conical_shock(2.0, float(degrees_to_radians(10.0)))
    assert result.shock_angle == pytest.approx(0.54464827, rel=2e-6)
    # NASA SP-3004 tabulates critical Mach numbers M*=V/a*.  Converting its
    # surface and post-shock values gives the ordinary local Mach numbers.
    assert result.surface_mach == pytest.approx(1.83403, rel=2e-5)
    assert result.post_shock_mach == pytest.approx(1.94679, rel=2e-5)
    assert result.surface_pressure_ratio == pytest.approx(1.2924832, rel=3e-5)
    assert result.surface_density_ratio == pytest.approx(1.2011081, rel=3e-5)
    assert result.surface_temperature_ratio == pytest.approx(1.0760757, rel=3e-5)


def test_zero_angle_conical_shock_is_a_mach_wave() -> None:
    result = conical_shock(2.0, 0.0)
    assert result.shock_angle == pytest.approx(np.arcsin(0.5))
    assert result.post_shock_mach == 2.0
    assert result.surface_mach == 2.0
    assert result.surface_pressure_ratio == 1.0
    assert result.surface_density_ratio == 1.0
    assert result.surface_temperature_ratio == 1.0
    assert result.total_pressure_ratio == 1.0


def test_near_sonic_zero_angle_conical_shock_remains_a_mach_wave() -> None:
    mach = np.nextafter(1.0, np.inf)
    result = conical_shock(mach, 0.0)
    assert result.shock_angle == pytest.approx(np.arcsin(1.0 / mach))
    assert 0.0 < result.shock_angle < 0.5 * np.pi
    assert result.surface_mach == mach
    assert result.total_pressure_ratio == 1.0


def test_near_sonic_conical_limit_uses_public_degeneracy_error() -> None:
    mach = np.nextafter(1.0, np.inf)
    with pytest.raises(
        NoAttachedShockError,
        match="conical-shock interval is physically or numerically degenerate",
    ):
        maximum_attached_cone_angle(mach)


@pytest.mark.parametrize(
    ("mach", "cone_angle"),
    [
        (1.01, 1e-4),
        (1.0 + 1e-10, 1e-12),
        (np.nextafter(1.0, np.inf), 1e-14),
    ],
)
def test_tiny_cone_solver_failure_uses_public_degeneracy_error(
    mach: float, cone_angle: float
) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        with pytest.raises(
            NoAttachedShockError,
            match=(
                r"conical-shock solution is physically or numerically degenerate.*"
                r"cone half-angle"
            ),
        ):
            conical_shock(mach, cone_angle)


def test_near_sonic_conical_array_normalizes_failing_element() -> None:
    with pytest.raises(NoAttachedShockError, match=r"Mach 1\.01"):
        conical_shock(
            [2.0, 1.01],
            [float(degrees_to_radians(10.0)), 1e-4],
        )


def test_resolvable_near_sonic_conical_root_remains_physical() -> None:
    mach = 1.0001
    cone_angle = 1e-4
    result = conical_shock(mach, cone_angle)
    assert cone_angle < result.shock_angle < 0.5 * np.pi
    assert result.post_shock_mach > result.surface_mach > 1.0
    assert result.surface_pressure_ratio > 1.0


def test_conical_shock_vectorizes_and_broadcasts() -> None:
    result = conical_shock([[2.0], [3.0]], degrees_to_radians([0.0, 5.0]))
    for value in (
        result.upstream_mach,
        result.cone_half_angle,
        result.shock_angle,
        result.post_shock_mach,
        result.surface_mach,
        result.surface_pressure_ratio,
        result.surface_density_ratio,
        result.surface_temperature_ratio,
        result.total_pressure_ratio,
    ):
        assert isinstance(value, np.ndarray)
        assert value.shape == (2, 2)
        assert value.dtype == np.float64
    assert_allclose(
        np.asarray(result.surface_pressure_ratio)
        / np.asarray(result.surface_density_ratio),
        result.surface_temperature_ratio,
    )
    limit = maximum_attached_cone_angle([2.0, 3.0])
    assert isinstance(limit.cone_half_angle, np.ndarray)
    assert np.asarray(limit.cone_half_angle).shape == (2,)


def test_conical_shock_accepts_custom_perfect_gas() -> None:
    helium = PerfectGas(specific_gas_constant=2077.1, heat_capacity_ratio=5.0 / 3.0)
    air = conical_shock(3.0, float(degrees_to_radians(10.0)))
    result = conical_shock(3.0, float(degrees_to_radians(10.0)), helium)
    assert result.shock_angle != pytest.approx(air.shock_angle)
    assert result.surface_pressure_ratio == pytest.approx(
        result.surface_density_ratio * result.surface_temperature_ratio
    )


def test_conical_attached_limit_and_detached_error() -> None:
    limit = maximum_attached_cone_angle(2.0)
    attached = conical_shock(2.0, limit.cone_half_angle)
    assert attached.shock_angle == pytest.approx(limit.shock_angle)
    with pytest.raises(NoAttachedShockError):
        conical_shock(2.0, float(limit.cone_half_angle) + 1e-6)


@pytest.mark.parametrize("mach", [1.0, 0.9, np.nan, np.inf])
def test_conical_shock_rejects_invalid_mach(mach: float) -> None:
    with pytest.raises(ValueError):
        conical_shock(mach, 0.1)
    with pytest.raises(ValueError):
        maximum_attached_cone_angle(mach)


@pytest.mark.parametrize("angle", [-0.1, 0.5 * np.pi, np.nan, np.inf])
def test_conical_shock_rejects_invalid_angle(angle: float) -> None:
    with pytest.raises(ValueError):
        conical_shock(2.0, angle)


def test_conical_shock_rejects_nonbroadcastable_inputs() -> None:
    with pytest.raises(ValueError, match="broadcastable"):
        conical_shock([2.0, 3.0], [0.1, 0.2, 0.3])
