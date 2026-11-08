"""Analytical, conservation and grid verification for rotational external MOC."""

import json
import time
import warnings
from dataclasses import fields, replace
from itertools import pairwise
from math import exp, log, radians, tan
from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

import aerophysics.characteristics as moc
from aerophysics import (
    AIR,
    AxisymmetricBodyFlowResult,
    BoundaryLayerEdge,
    BoundaryLayerRegime,
    PerfectGas,
    axisymmetric_body_flow,
    conical_shock,
    marching_boundary_layer,
    oblique_shock,
    prandtl_meyer_expansion,
)
from aerophysics._array import FloatArray
from aerophysics.exceptions import ModelRangeError
from aerophysics.transport import AIR_VISCOSITY


def solve(
    mach: float = 3, angle: float = 9.46, count: int = 88, step: float = 1.0
) -> AxisymmetricBodyFlowResult:
    radius = tan(radians(angle))
    return axisymmetric_body_flow(
        [0, 1, 1 + 20 * radius],
        [0, radius, radius],
        mach=mach,
        pressure=10000,
        temperature=220,
        characteristic_count=count,
        fan_angle_step=radians(step),
    )


@pytest.mark.parametrize(("mach", "angle"), [(2, 10), (6, 5)])
def test_m1_cone(mach: float, angle: float) -> None:
    rad = tan(radians(angle))
    result = axisymmetric_body_flow(
        [0, 1], [0, rad], mach=mach, pressure=10000, temperature=220
    )
    exact = conical_shock(mach, radians(angle))
    p_error = float(
        np.max(abs(result.pressure / (10000 * exact.surface_pressure_ratio) - 1))
    )
    m_error = float(np.max(abs(result.mach / exact.surface_mach - 1)))
    shock_error = float(
        np.max(
            abs(
                result.shock_radius
                - result.shock_axial_distance * tan(float(exact.shock_angle))
            )
        )
    )
    print("M1", mach, p_error, m_error, shock_error)
    assert max(p_error, m_error) < 1e-4
    assert shock_error < 1e-10


def test_collinear_cone_vertices() -> None:
    result = axisymmetric_body_flow(
        [0, 1, 2],
        [0, 0.1, 0.2],
        mach=3,
        pressure=10000,
        temperature=220,
        characteristic_count=48,
    )
    assert np.max(abs(result.pressure / result.pressure[0] - 1)) < 1e-4
    assert result.axial_distance[-1] == 2


def test_m2_planar_wedge_and_expansion() -> None:
    angle = radians(10)
    result, _ = moc._solve(
        [0, 1, 2],
        [0, tan(angle), tan(angle)],
        mach=3,
        pressure=10000,
        temperature=220,
        axisymmetric=False,
        fan_angle_step=radians(0.5),
    )
    wedge = oblique_shock(3, angle)
    expansion = prandtl_meyer_expansion(wedge.downstream_mach, angle)
    np.testing.assert_allclose(
        result.mach[result.axial_distance < 1], wedge.downstream_mach, rtol=1e-12
    )
    np.testing.assert_allclose(
        result.pressure[result.axial_distance < 1],
        10000 * wedge.static_pressure_ratio,
        rtol=1e-12,
    )
    error = float(
        np.max(
            abs(result.mach[result.axial_distance >= 1] / expansion.downstream_mach - 1)
        )
    )
    print("M2 planar post-corner Mach relative error", error)
    assert error < 2e-6


@pytest.mark.parametrize(("mach", "angle"), [(3, 9.46), (6, 4.76)])
def test_m5_corner_and_edge(mach: float, angle: float) -> None:
    result = solve(mach, angle)
    cone = conical_shock(mach, radians(angle))
    expansion = prandtl_meyer_expansion(cone.surface_mach, radians(angle))
    corner = int(np.flatnonzero(result.axial_distance == 1)[0])
    assert result.mach[corner] == pytest.approx(expansion.downstream_mach, rel=2e-12)
    assert result.pressure[corner] == pytest.approx(
        10000 * cone.surface_pressure_ratio * expansion.static_pressure_ratio, rel=2e-12
    )
    assert np.all(np.diff(result.surface_distance) > 0)
    assert result.surface_distance[0] == 0
    edge = result.boundary_layer_inputs().edge
    original = np.searchsorted(edge.distance, result.surface_distance)
    np.testing.assert_allclose(
        edge.temperature[original], result.temperature, rtol=2e-14
    )
    np.testing.assert_allclose(
        edge.velocity[original] / AIR.speed_of_sound(edge.temperature[original]),
        result.mach,
        rtol=2e-14,
    )
    assert (
        np.ptp(result.pressure / (result.temperature / result.total_temperature) ** 3.5)
        / result.total_pressure
        < 2e-14
    )
    print("M5", mach, "Mach", result.mach[corner], "pressure", result.pressure[corner])


def test_owned_readonly_and_nonair_edge() -> None:
    x, r = np.array([0.0, 1.0]), np.array([0.0, 0.1])
    gas = PerfectGas(300, 1.3)
    result = axisymmetric_body_flow(
        x, r, mach=3, pressure=10000, temperature=220, gas=gas
    )
    x[1], r[1] = 4, 9
    assert result.axial_distance[-1] == 1
    assert result.radius[-1] == pytest.approx(0.1)
    for field in fields(result):
        value = getattr(result, field.name)
        if isinstance(value, np.ndarray):
            assert value.flags.owndata
            assert not value.flags.writeable
            with pytest.raises(ValueError):
                value.flat[0] = 0
    copied = replace(result)
    assert not np.shares_memory(copied.pressure, result.pressure)
    edge = result.boundary_layer_inputs(viscosity_model=AIR_VISCOSITY).edge
    assert edge.gas is gas
    assert edge.viscosity_model is AIR_VISCOSITY
    np.testing.assert_allclose(edge.temperature, result.temperature, rtol=1e-14)


@pytest.mark.parametrize(
    ("x", "r"),
    [
        ([0, 1, 2], [0, 0.1, 0.3]),
        ([0, 1], [0, 4]),
        ([0, 1, 1], [0, 0.1, 0.2]),
        ([0, 1], [0, -0.1]),
        ([0, 1], [0, 0]),
        ([1, 2], [0, 0.1]),
        ([0, 1], [0.1, 0.2]),
        ([0], [0]),
        ([[0, 1]], [[0, 0.1]]),
        ([0, 1], [0, 0.1, 0.2]),
        ([0, np.nan], [0, 0.1]),
    ],
)
def test_m6_geometry(x: list[float], r: list[float]) -> None:
    with pytest.raises(ValueError):
        axisymmetric_body_flow(x, r, mach=3, pressure=10000, temperature=220)


@pytest.mark.parametrize("mach", [0.5, 1, 1.00000001])
def test_m6_mach(mach: float) -> None:
    with pytest.raises(ModelRangeError, match=r"x=0"):
        axisymmetric_body_flow(
            [0, 1], [0, 0.1], mach=mach, pressure=10000, temperature=220
        )


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("pressure", 0),
        ("temperature", -1),
        ("fan_angle_step", 0),
        ("mach", np.inf),
        ("mach", [3, 4]),
        ("characteristic_count", 3),
        ("characteristic_count", True),
        ("characteristic_count", 4.5),
    ],
)
def test_invalid_state_and_resolution(key: str, value: object) -> None:
    kw: dict[str, Any] = dict(mach=3, pressure=10000, temperature=220)
    kw[key] = value
    with pytest.raises(ValueError):
        axisymmetric_body_flow([0, 1], [0, 0.1], **kw)


def test_invalid_gas() -> None:
    with pytest.raises(TypeError, match=r"PerfectGas"):
        axisymmetric_body_flow(
            [0, 1],
            [0, 0.1],
            mach=3,
            pressure=10000,
            temperature=220,
            gas=cast(PerfectGas, None),
        )


def test_boattail_and_multiple_corners() -> None:
    result = axisymmetric_body_flow(
        [0, 1, 1.2, 1.5], [0, 0.2, 0.2, 0.17], mach=3, pressure=10000, temperature=220
    )
    assert result.flow_angle[-1] < 0
    for vertex in (1, 1.2, 1.5):
        assert np.count_nonzero(result.axial_distance == vertex) == 1
    np.testing.assert_allclose(
        result.radius,
        np.interp(result.axial_distance, [0, 1, 1.2, 1.5], [0, 0.2, 0.2, 0.17]),
        atol=1e-14,
    )


def mass_error(row: list[moc._Point], mach: float) -> float:
    """Independent integral of rho (u dr - v dx) 2 pi r along the polygon."""
    flow = moc._Flow(mach, AIR)
    data = np.array([[p.x, p.r, p.theta, p.p, flow.mach_at(p)] for p in row])
    x, r, theta, logp, mm = data.T
    tr = (1 + 0.2 * mach**2) / (1 + 0.2 * mm**2)
    density = np.exp(logp) / tr  # rho / rho_inf
    velocity = mm / mach * np.sqrt(tr)  # q / U_inf
    fx, fy = (
        density * velocity * r * np.cos(theta),
        density * velocity * r * np.sin(theta),
    )
    flux = np.sum((fx[1:] + fx[:-1]) * np.diff(r) - (fy[1:] + fy[:-1]) * np.diff(x))
    return float(flux / r[-1] ** 2 - 1)


@pytest.mark.parametrize(("mach", "angle"), [(3, 9.46), (6, 4.76)])
def test_m4_mass(mach: float, angle: float) -> None:
    r = tan(radians(angle))
    _, row = moc._solve(
        [0, 1, 1 + 20 * r], [0, r, r], mach=mach, pressure=10000, temperature=220
    )
    error = mass_error(row, mach)
    print("M4", mach, "relative mass error", error)
    assert abs(error) < 1e-4


def test_rotational_shock_entropy_and_mass() -> None:
    result, row = moc._solve(
        [0, 1, 5],
        [0, 0.1, 0.1],
        mach=6,
        pressure=10000,
        temperature=220,
        characteristic_count=48,
        fan_angle_step=radians(0.25),
    )
    totals = np.array([p.p0 for p in row])
    assert np.ptp(totals) > 1e-4
    assert totals.min() == pytest.approx(log(result.total_pressure / 10000))
    assert np.all(totals <= moc._Flow(6, AIR).upstream_p0)
    assert np.ptp(result.shock_angle) > 1e-3
    assert abs(mass_error(row, 6)) < 1e-4


@pytest.mark.solver_verification
@pytest.mark.parametrize(("mach", "angle"), [(3, 9.46), (6, 4.76)])
def test_m3_refinement(mach: float, angle: float) -> None:
    results = []
    times = []
    for n in (24, 48, 96, 192):
        start = time.perf_counter()
        results.append(solve(mach, angle, n, angle / (n / 2) * (1 + 1e-7)))
        times.append(time.perf_counter() - start)
    stations = np.linspace(1, results[-1].axial_distance[-1], 501)
    pressures: list[FloatArray] = []
    for result in results:
        mask = result.axial_distance >= 1
        pressures.append(
            CubicSpline(result.axial_distance[mask], result.pressure[mask])(stations)
        )
    differences = [float(np.max(abs(a - b))) for a, b in pairwise(pressures)]
    relative_differences = [
        float(np.max(abs(a / b - 1))) for a, b in pairwise(pressures)
    ]
    order = np.log2(np.array(differences[:-1]) / differences[1:])
    estimate = differences[-1] / 3
    print(
        "M3 "
        + json.dumps(
            dict(
                mach=mach,
                counts=[24, 48, 96, 192],
                differences_Pa=differences,
                relative_differences=relative_differences,
                orders=order.tolist(),
                fine_error_estimate_Pa=estimate,
                seconds=times,
            )
        )
    )
    assert 1.8 < order[-1] < 2.2
    assert differences[0] / 10000 < 1e-4
    assert estimate / 10000 < 2e-6
    default = solve(mach, angle)
    for region, positions in (
        ("shoulder_1_to_1.1", np.linspace(1, 1.1, 1001)),
        ("downstream_1.1_to_end", np.linspace(1.1, stations[-1], 1501)),
    ):
        sampled = []
        for result in (default, results[-2], results[-1]):
            mask = result.axial_distance >= 1
            sampled.append(
                CubicSpline(result.axial_distance[mask], result.pressure[mask])(
                    positions
                )
            )
        ordinary, medium, fine = sampled
        reference = fine + (fine - medium) / 3
        relative = float(np.max(abs(ordinary / reference - 1)))
        print(
            "M3_DEFAULT "
            + json.dumps(
                dict(
                    mach=mach,
                    region=region,
                    max_absolute_Pa=float(np.max(abs(ordinary - reference))),
                    max_relative=relative,
                    fine_correction_Pa=float(np.max(abs((fine - medium) / 3))),
                )
            )
        )
        assert relative < 1.2e-4


@pytest.mark.solver_verification
@pytest.mark.parametrize(
    "regime", [BoundaryLayerRegime.LAMINAR, BoundaryLayerRegime.TURBULENT]
)
def test_m7_boundary_layers(regime: BoundaryLayerRegime) -> None:
    flow = axisymmetric_body_flow(
        [0, 0.3, 1.3], [0, 0.05, 0.05], mach=3, pressure=10000, temperature=220
    )
    edge, radius = flow.boundary_layer_inputs()
    start = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        result = marching_boundary_layer(
            edge, regime=regime, wall_radius=radius, wall_temperature=300
        )
    assert any("pressure" in str(w.message) for w in caught)
    assert result.distance[-1] == pytest.approx(flow.surface_distance[-1])
    assert np.all(result.local_skin_friction_coefficient > 0)
    print(
        "M7",
        regime,
        "seconds",
        time.perf_counter() - start,
        "end_s",
        result.distance[-1],
        "warnings",
        [str(w.message) for w in caught],
    )
    if regime is BoundaryLayerRegime.LAMINAR:
        cone = flow.surface_distance < 0.25
        plate_edge = BoundaryLayerEdge.from_pressure(
            flow.surface_distance[cone],
            flow.pressure[cone],
            total_pressure=flow.total_pressure,
            total_temperature=flow.total_temperature,
        )
        plate = marching_boundary_layer(plate_edge, regime=regime, wall_temperature=300)
        ratio = (
            result.local_skin_friction_coefficient[: len(plate.distance)]
            / plate.local_skin_friction_coefficient
        )
        error = float(np.max(abs(ratio / np.sqrt(3) - 1)))
        print("M7 Mangler relative error", error)
        assert error < 5e-4


def test_rotational_parallel_shear() -> None:
    # Exact Euler solution: constant pressure and direction, entropy varying
    # normally to straight streamlines, with a uniform stagnation enthalpy.
    a = moc._Point(0, 1, 0, 0, 2)
    b = moc._Point(0, 2, 0, 0, 4)
    q = moc._interior(a, b, moc._Flow(3, AIR))
    assert q.theta == pytest.approx(0, abs=1e-12)
    assert q.p == pytest.approx(0, abs=1e-12)
    assert q.p0 == pytest.approx(2 + 2 * (q.r - 1), abs=1e-12)


@pytest.mark.parametrize("angle", [10, 20])
def test_local_fitted_shock_state(angle: float) -> None:
    expected = oblique_shock(3, radians(angle))
    flow = moc._Flow(3, AIR)
    actual = flow.shock(1, 0.5, float(expected.shock_angle))
    assert actual.theta == pytest.approx(radians(angle), abs=1e-12)
    assert exp(actual.p) == pytest.approx(expected.static_pressure_ratio, rel=1e-12)
    assert flow.mach_at(actual) == pytest.approx(expected.downstream_mach, rel=1e-12)
    assert exp(actual.p0 - flow.upstream_p0) == pytest.approx(
        expected.total_pressure_ratio, rel=1e-12
    )


def test_characteristic_domain_errors() -> None:
    flow = moc._Flow(3, AIR)
    q = moc._Point(1, 1, 0, 0, 2)
    with pytest.raises(ModelRangeError, match=r"subsonic.*x=1"):
        flow.mach_at(replace(q, p0=0.01))
    with pytest.raises(ModelRangeError, match=r"cannot march.*x=1"):
        flow.coefficients(replace(q, theta=1.5), 1)
    with pytest.raises(ModelRangeError, match=r"invalid shock angle.*x=1"):
        flow.shock(1, 1, 0.1)
    with pytest.raises(ModelRangeError, match=r"parallel.*x=1"):
        moc._intersection(q, replace(q, r=2), 1, 1)
    with pytest.raises(ModelRangeError, match=r"crossed.*x="):
        moc._intersection(q, replace(q, r=2), -1, 1)
    with pytest.raises(ModelRangeError, match=r"degenerate.*x=1"):
        moc._entropy(q, replace(q, x=0), q)
    with pytest.raises(ModelRangeError, match=r"left its.*x=1"):
        moc._entropy(q, replace(q, r=2), replace(q, r=3))


@pytest.mark.parametrize("operation", ["interior", "wall", "shock"])
def test_unit_process_iteration_exhaustion(
    monkeypatch: pytest.MonkeyPatch, operation: str
) -> None:
    flow = moc._Flow(3, AIR, axisymmetric=False)
    beta = float(oblique_shock(3, radians(10)).shock_angle)
    a = flow.shock(0, 1, beta)
    b = flow.shock(0, 2, beta)
    monkeypatch.setattr(moc, "_MAX_ITERATIONS", 1)
    with pytest.raises(ModelRangeError, match=r"corrector did not converge.*x="):
        if operation == "interior":
            moc._interior(a, b, flow)
        elif operation == "wall":
            moc._wall(b, a, a.theta, flow)
        else:
            moc._shock_point(a, b, beta, flow)


def test_shock_root_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    from types import SimpleNamespace

    flow = moc._Flow(3, AIR)
    row, beta, _ = moc._initial_line(1, radians(10), 24, flow)
    monkeypatch.setattr(
        moc, "root_scalar", lambda *a, **k: SimpleNamespace(converged=False)
    )
    with pytest.raises(ModelRangeError, match=r"shock corrector.*x="):
        moc._shock_point(row[-2], row[-1], beta, flow)


@pytest.mark.parametrize("failed_call", [1, 2])
def test_starting_line_ode_failure(
    monkeypatch: pytest.MonkeyPatch, failed_call: int
) -> None:
    from types import SimpleNamespace

    original = solve_ivp
    count = 0

    def fail(*args: Any, **kwargs: Any) -> Any:
        nonlocal count
        count += 1
        if count == failed_call:
            return SimpleNamespace(success=False, sol=None)
        return original(*args, **kwargs)

    monkeypatch.setattr(moc, "solve_ivp", fail)
    with pytest.raises(
        ModelRangeError, match=r"(starting line failed|integration failed).*x="
    ):
        solve()


@pytest.mark.parametrize("budget", ["_MAX_STATIONS", "_MAX_FAN_STEPS"])
def test_work_limit_failure(monkeypatch: pytest.MonkeyPatch, budget: str) -> None:
    monkeypatch.setattr(moc, budget, 1)
    with pytest.raises(ModelRangeError, match=r"budget exceeded.*x="):
        solve()


def test_expansion_limit_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    def fail(*args: Any, **kwargs: Any) -> Any:
        raise ValueError("expansion reaches vacuum")

    monkeypatch.setattr(moc, "prandtl_meyer_expansion", fail)
    with pytest.raises(ModelRangeError, match=r"unresolvable expansion.*x=1"):
        solve()


def test_shock_updates_after_fan_and_collinear_segment() -> None:
    result = axisymmetric_body_flow(
        [0, 1, 5, 10],
        [0, 0.2, 0.2 + 4 * (1.6 / 9), 1.8],
        mach=3,
        pressure=10000,
        temperature=220,
    )
    assert len(result.shock_angle) > 10
    assert np.all(np.diff(result.shock_axial_distance) > 0)
    assert result.axial_distance[-1] == 10


def test_short_remaining_characteristic_front() -> None:
    # Exercises the two-point inverse-wall stencil after the outer front has
    # left the downstream dependency region; no fictitious outflow state.
    result = axisymmetric_body_flow(
        [0, 1, 5],
        [0, tan(radians(4.76)), tan(radians(4.76))],
        mach=6,
        pressure=10000,
        temperature=220,
        fan_angle_step=radians(0.25),
    )
    assert result.axial_distance[-1] == 5


def test_unresolved_long_body_reports_location() -> None:
    with pytest.raises(ModelRangeError, match=r"crossed characteristics at x=.*r="):
        axisymmetric_body_flow(
            [0, 1, 10], [0, 0.2, 1.1], mach=3, pressure=10000, temperature=220
        )


def test_streamline_iteration_exhaustion(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(moc, "_MAX_ITERATIONS", 1)
    with pytest.raises(ModelRangeError, match=r"streamline corrector.*x="):
        moc._entropy(
            moc._Point(0, 1, 0.1, 0, 2),
            moc._Point(0, 2, 0.2, 0, 3),
            moc._Point(1, 1.3, 0.15, 0, 2.5),
        )


@pytest.mark.parametrize("failure", ["unconverged", "bracket", "physical"])
def test_inverse_wall_failure_reporting(
    monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    from types import SimpleNamespace

    flow = moc._Flow(3, AIR)
    row, _, wall = moc._initial_line(1, radians(10), 24, flow)

    def fail(*args: Any, **kwargs: Any) -> Any:
        if failure == "bracket":
            raise ValueError("not bracketed")
        if failure == "physical":
            raise ModelRangeError("subsonic at x=1 m, r=0.2 m")
        return SimpleNamespace(converged=False)

    monkeypatch.setattr(moc, "root_scalar", fail)
    with pytest.raises(ModelRangeError, match=r"(inverse wall root|subsonic).*x="):
        moc._wall_at(row, wall, wall.theta, 1.01, flow)


def test_forced_cone_uses_a_marching_net() -> None:
    result, _ = moc._solve(
        [0, 1],
        [0, tan(radians(10))],
        mach=2,
        pressure=10000,
        temperature=220,
        starting_distance=0.05,
    )
    assert result.axial_distance[-1] == 1
    assert len(result.pressure[result.axial_distance > 0.05]) > 20
    # Numerical truncation is nonzero: returning the analytic initializer
    # everywhere cannot satisfy this independent check of the forced path.
    assert np.ptp(result.pressure) / result.pressure[0] > 1e-4
    assert len(result.shock_angle) > 10


@pytest.mark.parametrize("start", [0, 1, np.nan])
def test_invalid_internal_start(start: float) -> None:
    with pytest.raises(ValueError, match="starting_distance"):
        moc._solve(
            [0, 1],
            [0, 0.1],
            mach=3,
            pressure=10000,
            temperature=220,
            starting_distance=start,
        )


@pytest.mark.solver_verification
@pytest.mark.parametrize(
    ("mach", "angle", "limits"),
    [
        (2, 10, [0.015, 0.005, 0.0025]),
        (6, 5, [0.001, 0.00016, 0.001]),
    ],
)
def test_m1_marched_cone(mach: float, angle: float, limits: list[float]) -> None:
    exact = conical_shock(mach, radians(angle))
    errors = []
    for count in (24, 48, 96, 192):
        result, _ = moc._solve(
            [0, 1],
            [0, tan(radians(angle))],
            mach=mach,
            pressure=10000,
            temperature=220,
            characteristic_count=count,
            starting_distance=0.05,
        )
        wall = result.axial_distance > 0.05
        shock = result.shock_axial_distance > 0
        errors.append(
            [
                float(
                    np.max(
                        abs(
                            result.pressure[wall]
                            / (10000 * exact.surface_pressure_ratio)
                            - 1
                        )
                    )
                ),
                float(np.max(abs(result.mach[wall] / exact.surface_mach - 1))),
                float(
                    np.max(
                        abs(
                            result.shock_radius[shock]
                            / (
                                result.shock_axial_distance[shock]
                                * tan(float(exact.shock_angle))
                            )
                            - 1
                        )
                    )
                ),
            ]
        )
        # All compared surface stations beyond x0 came from wall unit processes.
        assert len(result.pressure[wall]) > 20
    error = np.array(errors)
    orders = np.log2(error[:-1] / error[1:])
    print(
        "M1_MOC "
        + json.dumps(
            dict(
                mach=mach,
                angle_deg=angle,
                start=0.05,
                end=1.0,
                counts=[24, 48, 96, 192],
                errors=errors,
                orders=orders.tolist(),
            )
        )
    )
    assert np.all(error[0] < np.asarray(limits))
    assert np.all((orders[-1] > 1.8) & (orders[-1] < 2.2))


@pytest.mark.parametrize("case_index", [0, 1])
def test_naca_tn2250(case_index: int) -> None:
    path = (
        Path(__file__).parent / "reference_data/inviscid/naca_tn2250_cone_cylinder.json"
    )
    case = json.loads(path.read_text())["cases"][case_index]
    mach, radius = case["free_stream_mach"], 0.5 / case["nose_fineness_ratio"]
    result = axisymmetric_body_flow(
        [0, 1, 3.4], [0, radius, radius], mach=mach, pressure=10000, temperature=220
    )
    data = np.asarray(case["points"])
    mask = result.axial_distance >= 1
    predicted = CubicSpline(
        result.axial_distance[mask], result.pressure[mask] / 10000 - 1
    )(data[:, 0] / 100)
    error = float(np.max(abs(predicted - data[:, 1])))
    print(
        "TN2250 "
        + json.dumps(
            dict(
                mach=mach,
                predicted=predicted.tolist(),
                reference=data[:, 1].tolist(),
                max_absolute_difference=error,
            )
        )
    )
    assert error < 0.02


@pytest.mark.solver_verification
@pytest.mark.parametrize(("mach", "angle"), [(3, 9.46), (6, 4.76)])
def test_terminal_total_pressure_range(
    monkeypatch: pytest.MonkeyPatch, mach: float, angle: float
) -> None:
    radius = tan(radians(angle))
    end = 1 + 20 * radius
    original = moc._row
    samples: list[float] = []

    def record(*args: Any, **kwargs: Any) -> tuple[list[moc._Point], float, bool]:
        value = original(*args, **kwargs)
        points = value[0]
        x = [p.x for p in points]
        if x[0] <= end <= x[-1]:
            samples.append(float(np.interp(end, x, [p.p0 for p in points])))
        return value

    monkeypatch.setattr(moc, "_row", record)
    for count, step in ((24, 1.0), (192, 0.125)):
        samples.clear()
        result, _ = moc._solve(
            [0, 1, end],
            [0, radius, radius],
            mach=mach,
            pressure=10000,
            temperature=220,
            characteristic_count=count,
            fan_angle_step=radians(step),
        )
        # The returned last complete characteristic is NOT the terminal x cut.
        # Sample every active row's intersection with the terminal plane and
        # close the cut with the body and the local interpolated fitted shock.
        beta = float(np.interp(end, result.shock_axial_distance, result.shock_angle))
        shock_r = float(
            np.interp(end, result.shock_axial_distance, result.shock_radius)
        )
        shock = moc._Flow(mach, AIR).shock(end, shock_r, beta)
        totals = 10000 * np.exp(
            [log(result.total_pressure / 10000), shock.p0, *samples]
        )
        span = float(np.ptp(totals))
        print(
            "TERMINAL_P0 "
            + json.dumps(
                dict(
                    mach=mach,
                    count=count,
                    x=end,
                    samples=len(samples),
                    minimum_Pa=float(totals.min()),
                    maximum_Pa=float(totals.max()),
                    range_Pa=span,
                    range_over_body=span / result.total_pressure,
                )
            )
        )
        assert samples
        assert span / result.total_pressure < 1.5e-4
        assert span > 0 if mach == 3 else span == 0


@pytest.mark.solver_verification
@pytest.mark.parametrize(
    "regime", [BoundaryLayerRegime.LAMINAR, BoundaryLayerRegime.TURBULENT]
)
def test_shoulder_pressure_width_sensitivity(
    monkeypatch: pytest.MonkeyPatch, regime: BoundaryLayerRegime
) -> None:
    # The fixed-output V20 harness preserves each edge interpolant. Resampling
    # pressure onto the output grid would instead change the experiment.
    from test_marching_boundary_layer_axisymmetric import axisymmetric_march

    import aerophysics._boundary_layer_march as march_module
    from aerophysics.marching_boundary_layer import (
        _EdgeInterpolant,
        _warn_applicability,
    )

    flow = axisymmetric_body_flow(
        [0, 0.3, 1.3], [0, 0.05, 0.05], mach=3, pressure=10000, temperature=220
    )
    s, p, r = flow.surface_distance, flow.pressure, flow.radius
    corner = int(np.flatnonzero(flow.axial_distance == 0.3)[0])
    shoulder = float(s[corner])
    widths = [0.001 * shoulder, 0.001, 0.0002]
    cone_edge = BoundaryLayerEdge.from_pressure(
        s[: corner + 1],
        np.full(corner + 1, p[0]),
        total_pressure=flow.total_pressure,
        total_temperature=flow.total_temperature,
    )
    with warnings.catch_warnings(record=True):
        cone = marching_boundary_layer(
            cone_edge, regime=regime, wall_temperature=300, wall_radius=r[: corner + 1]
        )
    delta = float(cone.boundary_layer_thickness[-1])
    targets = shoulder + np.array([5, 20, 100]) * delta
    # One uniform ramp grid avoids mixing three differently spaced output
    # grids. Edge interpolation is unchanged, and the half-step check below
    # separately measures streamwise-discretization error.
    common = np.unique(
        np.concatenate(
            (
                s[1:],
                targets,
                np.linspace(shoulder - max(widths), shoulder, 65),
            )
        )
    )
    indices = np.searchsorted(common, targets)
    original = march_module._step_is_resolved
    accepted = [0]

    def count_steps(*args: Any, **kwargs: Any) -> bool:
        resolved = original(*args, **kwargs)
        accepted[0] += int(resolved)
        return resolved

    monkeypatch.setattr(march_module, "_step_is_resolved", count_steps)
    results: list[dict[str, Any]] = []
    for index, width in enumerate(widths):
        edge, radius = flow.boundary_layer_inputs(
            corner_length=None if index == 0 else width
        )
        accepted[0] = 0
        start = time.perf_counter()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            out = axisymmetric_march(edge, radius, regime, 300, common)
            ip = _EdgeInterpolant(
                edge,
                np.full(edge.distance.size, 300.0),
                3.0 if regime is BoundaryLayerRegime.TURBULENT else None,
                radius,
            )
            _warn_applicability(edge, ip, out)
        record = dict(
            regime=str(regime),
            width_m=width,
            delta_upstream_m=delta,
            target_s=targets.tolist(),
            Cf=out.local_skin_friction_coefficient[indices].tolist(),
            heat_flux=out.wall_heat_flux[indices].tolist(),
            delta99=out.boundary_layer_thickness[indices].tolist(),
            input_stations=edge.distance.size,
            output_stations=common.size,
            accepted_steps=accepted[0],
            seconds=time.perf_counter() - start,
            pressure_length_max=float(
                np.max(
                    out.boundary_layer_thickness * abs(ip.log_pressure_slope(common))
                )
            ),
            warnings=[str(w.message) for w in caught],
        )
        results.append(record)
        print("SHOULDER_WIDTH " + json.dumps(record))
        assert out.distance[-1] == flow.surface_distance[-1]
        assert np.all(out.local_skin_friction_coefficient > 0)
    # Check marching-step contamination of the narrowest reference separately
    # from the deliberately varied edge interpolation width.
    dense = np.unique(np.r_[common, (np.r_[0.0, common[:-1]] + common) / 2])
    refined_indices = np.searchsorted(dense, targets)
    accepted[0] = 0
    refined = axisymmetric_march(edge, radius, regime, 300, dense, relative_step=0.015)
    print("SHOULDER_REFINED_STEPS", str(regime), accepted[0])
    for quantity, values in (
        ("Cf", refined.local_skin_friction_coefficient[refined_indices]),
        ("heat_flux", refined.wall_heat_flux[refined_indices]),
        ("delta99", refined.boundary_layer_thickness[refined_indices]),
    ):
        reference = np.asarray(results[-1][quantity])
        relative = [
            100 * (np.asarray(record[quantity]) / reference - 1)
            for record in results[:-1]
        ]
        stream_error = float(np.max(abs(values / reference - 1)))
        print(
            "SHOULDER_DIFFERENCE "
            + json.dumps(
                dict(
                    regime=str(regime),
                    quantity=quantity,
                    default_percent=relative[0].tolist(),
                    one_mm_percent=relative[1].tolist(),
                    default_vs_one_mm_percent=(
                        100
                        * (
                            np.asarray(results[0][quantity])
                            / np.asarray(results[1][quantity])
                            - 1
                        )
                    ).tolist(),
                    half_step_max_relative=stream_error,
                )
            )
        )
        assert stream_error < 0.002


@pytest.mark.parametrize("width", [0, -0.001, np.nan, np.inf, True, [0.001]])
def test_invalid_corner_length(width: Any) -> None:
    flow = axisymmetric_body_flow(
        [0, 1], [0, 0.1], mach=3, pressure=10000, temperature=220
    )
    with pytest.raises(ValueError, match="corner_length"):
        flow.boundary_layer_inputs(corner_length=width)


def test_corner_inputs_preserve_limits_and_geometry() -> None:
    x, r = [0, 1, 1.2, 1.5], [0, 0.2, 0.2, 0.17]
    flow = axisymmetric_body_flow(x, r, mach=3, pressure=10000, temperature=220)
    inputs = flow.boundary_layer_inputs()
    edge, radius = inputs
    assert not hasattr(flow, "boundary_layer_edge")
    assert edge.distance.size == flow.surface_distance.size + 4
    assert not radius.flags.writeable
    assert radius.flags.owndata
    assert not np.shares_memory(radius, flow.radius)
    np.testing.assert_allclose(
        radius,
        np.interp(
            edge.distance,
            np.r_[0, np.cumsum(np.hypot(np.diff(x), np.diff(r)))],
            r,
        ),
        atol=2e-16,
    )
    indices = np.searchsorted(edge.distance, flow.surface_distance)
    np.testing.assert_array_equal(edge.pressure[indices], flow.pressure)
    for vertex in (1, 2):
        truncated = axisymmetric_body_flow(
            x[: vertex + 1], r[: vertex + 1], mach=3, pressure=10000, temperature=220
        )
        corner = int(np.flatnonzero(flow.axial_distance == x[vertex])[0])
        sc = flow.surface_distance[corner]
        j = int(np.flatnonzero(edge.distance == sc)[0])
        assert edge.distance[j - 2] == pytest.approx(sc * 0.999, abs=1e-15)
        assert edge.distance[j - 1] == pytest.approx(sc * 0.9995, abs=1e-15)
        assert edge.pressure[j - 2] == pytest.approx(truncated.pressure[-1], rel=2e-12)
        assert edge.pressure[j - 1] == pytest.approx(
            np.sqrt(edge.pressure[j - 2] * edge.pressure[j]), rel=2e-15
        )
    gap = flow.surface_distance[corner] - flow.surface_distance[corner - 1]
    for width in (gap, gap * 2, 1e-30):
        with pytest.raises(ValueError, match="preceding output interval"):
            flow.boundary_layer_inputs(corner_length=width)


@pytest.mark.solver_verification
def test_default_marched_cone_accuracy_and_minimum_count() -> None:
    exact = conical_shock(2, radians(10))
    errors = []
    options_list: list[dict[str, Any]] = [{"characteristic_count": 87}, {}]
    for options in options_list:
        result, _ = moc._solve(
            [0, 1],
            [0, tan(radians(10))],
            mach=2,
            pressure=10000,
            temperature=220,
            starting_distance=0.05,
            **options,
        )
        errors.append(
            float(
                np.max(
                    abs(result.pressure / (10000 * exact.surface_pressure_ratio) - 1)
                )
            )
        )
    print("DEFAULT_CONE", errors)
    assert errors[0] >= 0.001
    assert errors[1] < 0.001
    # The public entry point uses the same default as the private verification.
    public = axisymmetric_body_flow(
        [0, 1], [0, 0.1], mach=3, pressure=10000, temperature=220
    )
    assert public.axial_distance.size == 89
