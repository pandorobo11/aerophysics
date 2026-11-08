"""Specified SA production onset: numerical verification, not transition prediction.

Run with -s to emit the T1--T8 measurements documented in the design note.
"""

import warnings
from dataclasses import replace

import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator
from test_marching_boundary_layer import smooth_edge, tighten_streamwise_control
from test_marching_boundary_layer_sa import constant_edge, discretization

from aerophysics import (
    AIR,
    BoundaryLayerEdge,
    BoundaryLayerRegime,
    MarchingBoundaryLayerResult,
    TurbulentModel,
    flat_plate_boundary_layer,
    marching_boundary_layer,
)
from aerophysics import _boundary_layer_march as march_module
from aerophysics import _boundary_layer_sa as sa_module
from aerophysics._array import FloatArray
from aerophysics._boundary_layer_march import StationEdge, eta_grid, march
from aerophysics.exceptions import ApplicabilityWarning, BoundaryLayerSeparationError
from aerophysics.marching_boundary_layer import (
    _EdgeInterpolant,
    _result,
    _turbulent_eta_grid,
    _warn_applicability,
)

LAMINAR = BoundaryLayerRegime.LAMINAR
TURBULENT = BoundaryLayerRegime.TURBULENT
TRANSITIONAL = BoundaryLayerRegime.TRANSITIONAL
CASES = [
    pytest.param(0.3, 300.0, 250.0, id="M0.3-cold"),
    pytest.param(2.0, 220.0, 300.0, id="M2-cold"),
    pytest.param(2.0, 220.0, None, id="M2-adiabatic"),
    pytest.param(6.0, 65.0, 125.0, id="M6-cold"),
]


def relative_difference(a: FloatArray, b: FloatArray) -> float:
    return float(np.max(abs(a / b - 1.0)))


@pytest.mark.parametrize(
    "transition", [None, 0.0, -1.0, 1.0, 2.0, np.nan, np.inf, -np.inf, [0.1]]
)
def test_transition_distance_is_required_and_bounded(transition: float | None) -> None:
    edge = constant_edge(2.0, 220.0, [1.0])
    with pytest.raises(ValueError, match="transition_distance"):
        marching_boundary_layer(
            edge, regime=TRANSITIONAL, transition_distance=transition
        )


@pytest.mark.parametrize("regime", [LAMINAR, TURBULENT])
def test_transition_distance_is_rejected_in_other_regimes(
    regime: BoundaryLayerRegime,
) -> None:
    with pytest.raises(ValueError, match="only valid for TRANSITIONAL"):
        marching_boundary_layer(
            constant_edge(2.0, 220.0, [1.0]), regime=regime, transition_distance=0.1
        )


def test_zero_seed_explains_why_transition_cannot_develop() -> None:
    with pytest.raises(ValueError, match="proportional to nu_tilde"):
        marching_boundary_layer(
            constant_edge(2.0, 220.0, [1.0]),
            regime=TRANSITIONAL,
            transition_distance=0.1,
            edge_sa_ratio=0.0,
        )


def test_switch_changes_only_production_and_includes_the_transition_station() -> None:
    edge = constant_edge(2.0, 220.0, [0.05, 0.1, 0.2])
    ordinary = _EdgeInterpolant(edge, None, 3.0)
    transition = _EdgeInterpolant(edge, None, 3.0, transition_distance=0.1)
    for s in [0.0, 0.05, 0.1, np.nextafter(0.1, np.inf), 0.2]:
        state = transition.state(s)
        assert state.turbulence_production == float(s > 0.1)
        assert replace(state, turbulence_production=1.0) == ordinary.state(s)


def test_residual_suppresses_only_the_sa_production_row() -> None:
    solver = discretization()
    edge = StationEdge(0.01, 1.0, 0.1, 1.6, 220.0, 1.3, 2000.0, 0.05, 3.0)
    state = solver.starting_guess(edge)
    lagged = -100.0 * state
    on = solver.residual(state, edge, 100.0, lagged)
    off = solver.residual(
        state, replace(edge, turbulence_production=0.0), 100.0, lagged
    )
    np.testing.assert_array_equal(on[[0, 1, 2, 4]], off[[0, 1, 2, 4]])
    np.testing.assert_array_equal(on[3, [0, -1]], off[3, [0, -1]])
    g, n = state[2, 1:-1], state[3, 1:-1]
    c = solver.chapman_rubesin(state[2], edge.temperature)[1:-1]
    modified, _, _ = sa_module._source_terms(
        edge.reynolds_root * abs(solver.derivative(state[1])) / g,
        n,
        c,
        g,
        state[4, 1:-1],
    )
    # Subtraction of full residuals loses the tiny outer-source digits.
    # Bound cancellation by machine precision of the terms being subtracted.
    error = abs((on[3, 1:-1] - off[3, 1:-1]) - 0.1355 * modified * n)
    roundoff = (
        8.0
        * np.finfo(float).eps
        * np.maximum(abs(on[3, 1:-1]) + abs(off[3, 1:-1]), 1.0)
    )
    assert np.all(error <= roundoff)


@pytest.mark.parametrize("production", [0.0, 1.0])
def test_picard_and_newton_solve_the_same_suppressed_equations(
    production: float,
) -> None:
    solver = discretization()
    edge = StationEdge(0.01, 1.0, 0.15, 1.6, 220.0, 1.3, 2000.0, 0.05, 3.0, production)
    initial = solver.solve(
        solver.starting_guess(edge), edge, 0.0, np.zeros((5, len(solver.eta)))
    )
    assert initial is not None
    edge = replace(edge, xi=0.0101, beta=-0.05, viscosity_growth=-0.02)
    coefficient = 10000.0
    lagged = -coefficient * initial.state
    newton = solver.newton(initial.state, edge, coefficient, lagged)
    picard = solver.picard(initial.state, edge, coefficient, lagged)
    assert newton is not None and picard is not None
    np.testing.assert_allclose(picard.state, newton.state, atol=1e-6, rtol=1e-8)
    residual = solver.residual(picard.state, edge, coefficient, lagged)
    residual[1:4, 1:-1] *= solver.width**2
    residual[3] /= max(1.0, float(np.max(picard.state[3])))
    assert np.max(abs(residual)) < 1e-5


@pytest.fixture(scope="module")
def short_transition() -> tuple[BoundaryLayerEdge, MarchingBoundaryLayerResult]:
    edge = constant_edge(2.0, 220.0, [0.05, 0.1, 0.11])
    result = marching_boundary_layer(
        edge, regime=TRANSITIONAL, transition_distance=0.1, wall_temperature=300.0
    )
    return edge, result


def test_public_result_uses_sa_grid_and_owned_read_only_profiles(
    short_transition: tuple[BoundaryLayerEdge, MarchingBoundaryLayerResult],
) -> None:
    edge, result = short_transition
    np.testing.assert_array_equal(result.distance, edge.distance[1:])
    np.testing.assert_array_equal(result.eta, _turbulent_eta_grid(1.1e6))
    assert result.nu_tilde_ratio is not None
    assert result.nu_tilde_ratio.shape == result.velocity_ratio.shape
    assert not result.nu_tilde_ratio.flags.writeable
    assert not np.shares_memory(result.nu_tilde_ratio, result.velocity_ratio)
    assert np.all(result.nu_tilde_ratio[:, 0] == 0.0)
    np.testing.assert_allclose(result.nu_tilde_ratio[:, -1], 3.0)


@pytest.mark.parametrize("downstream", [False, True])
def test_y_plus_warning_ignores_upstream_and_transition_node(
    short_transition: tuple[BoundaryLayerEdge, MarchingBoundaryLayerResult],
    downstream: bool,
) -> None:
    edge, result = short_transition
    interp = _EdgeInterpolant(edge, np.full(4, 300.0), 3.0, transition_distance=0.1)
    y = result.wall_distance.copy()
    y[:, 1] = [1e6, 1e6, 1e6 if downstream else 0.0]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ApplicabilityWarning)
        _warn_applicability(edge, interp, replace(result, wall_distance=y))
    assert any("y+ > 1" in str(w.message) for w in caught) is downstream
    # An upstream-only partial result must not reduce an empty downstream array.
    interp.transition_distance = 0.2
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ApplicabilityWarning)
        _warn_applicability(edge, interp, replace(result, wall_distance=y))
    assert not any("y+ > 1" in str(w.message) for w in caught)


@pytest.mark.solver_verification
@pytest.mark.parametrize(("mach", "temperature", "wall"), CASES)
def test_t1_small_seed_upstream_laminar_limit(
    mach: float,
    temperature: float,
    wall: float | None,
) -> None:
    edge = constant_edge(mach, temperature, [0.01, 0.05, 0.1, 0.10001])
    laminar = marching_boundary_layer(edge, regime=LAMINAR, wall_temperature=wall)
    transition = marching_boundary_layer(
        edge,
        regime=TRANSITIONAL,
        wall_temperature=wall,
        edge_sa_ratio=1e-3,
        transition_distance=0.1,
    )
    zero_sa = marching_boundary_layer(
        edge, regime=TURBULENT, wall_temperature=wall, edge_sa_ratio=0.0
    )
    np.testing.assert_array_equal(transition.eta, zero_sa.eta)
    sa_cf = relative_difference(
        transition.local_skin_friction_coefficient[:3],
        zero_sa.local_skin_friction_coefficient[:3],
    )
    if wall is None:
        np.testing.assert_array_equal(transition.wall_heat_flux[:3], np.zeros(3))
        np.testing.assert_array_equal(zero_sa.wall_heat_flux[:3], np.zeros(3))
        sa_heat = 0.0
    else:
        sa_heat = relative_difference(
            transition.wall_heat_flux[:3], zero_sa.wall_heat_flux[:3]
        )
    print(f"T1(a) M={mach} Tw={wall}: Cf={100 * sa_cf:.9g}%, q={100 * sa_heat:.9g}%")
    # Same SA grid/operators isolate the production switch from discretization error.
    # 1e-10 relative = 1e-8%; measured differences are at roundoff (design note 13).
    assert sa_cf < 1e-10
    assert sa_heat < 1e-10

    cf = relative_difference(
        transition.local_skin_friction_coefficient[:3],
        laminar.local_skin_friction_coefficient[:3],
    )
    heat = (
        0.0
        if wall is None
        else relative_difference(
            transition.wall_heat_flux[:3], laminar.wall_heat_flux[:3]
        )
    )
    print(f"T1(b) M={mach} Tw={wall}: Cf={100 * cf:.9g}%, q={100 * heat:.9g}%")
    # V6 judges Cf only; heat-flux differences between discretizations are diagnostic.
    assert cf < 5e-4


@pytest.mark.solver_verification
@pytest.mark.parametrize(("mach", "temperature", "wall"), CASES)
def test_t2_t3_t4_uniform_edge(
    mach: float,
    temperature: float,
    wall: float | None,
) -> None:
    x = np.unique(np.r_[np.geomspace(0.01, 10.0, 350), 0.1, 0.05, 0.2, 0.5, 1, 2, 5])
    edge = constant_edge(mach, temperature, x.tolist())
    laminar = marching_boundary_layer(edge, regime=LAMINAR, wall_temperature=wall)
    transition = marching_boundary_layer(
        edge, regime=TRANSITIONAL, wall_temperature=wall, transition_distance=0.1
    )
    turbulent = marching_boundary_layer(edge, regime=TURBULENT, wall_temperature=wall)
    upstream = x <= 0.1
    cf = relative_difference(
        transition.local_skin_friction_coefficient[upstream],
        laminar.local_skin_friction_coefficient[upstream],
    )
    heat = (
        0.0
        if wall is None
        else relative_difference(
            transition.wall_heat_flux[upstream], laminar.wall_heat_flux[upstream]
        )
    )
    print(f"T2 M={mach} Tw={wall}: Cf={100 * cf:.9g}%, q={100 * heat:.9g}%")
    # Regression bounds only; residual eddy viscosity is not physical validation.
    assert cf < 0.0045
    assert heat < 0.0045
    rt = turbulent.momentum_thickness * 1e7
    rr = transition.momentum_thickness * 1e7
    developed = x >= 0.1
    valid = (x > 0.1) & (rr >= rt[developed].min()) & (rr <= rt.max())
    assert np.all(np.diff(rt[developed]) > 0.0)
    difference = np.full(x.size, np.inf)
    difference[valid] = abs(
        transition.local_skin_friction_coefficient[valid]
        / PchipInterpolator(
            rt[developed],
            turbulent.local_skin_friction_coefficient[developed],
            extrapolate=False,
        )(rr[valid])
        - 1.0
    )
    settled = np.flatnonzero(np.maximum.accumulate(difference[::-1])[::-1] < 0.01)
    assert settled.size
    first = int(settled[0])
    delta = transition.boundary_layer_thickness[x == 0.1][0]
    far_error = float(np.max(difference[x >= 5.0]))
    print(
        f"T3 M={mach} Tw={wall}: s={x[first]:.9g} m, "
        f"ds/delta99(tr)={(x[first] - 0.1) / delta:.9g}, "
        f"Re_theta={rr[first]:.9g}..{rr[-1]:.9g}, far={100 * far_error:.9g}%"
    )
    assert far_error < 2e-5  # 0.002%, rounded above the measured 0.001172%.
    flat = flat_plate_boundary_layer(
        x,
        edge.velocity[0],
        edge.pressure[0] / (AIR.specific_gas_constant * temperature),
        edge.viscosity_model.dynamic_viscosity(temperature),
        regime=TRANSITIONAL,
        transition_reynolds=1e6,
        compressible=True,
        turbulent_model=TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE,
        mach=mach,
        edge_temperature=temperature,
        wall_temperature=wall,
    )
    for field in ["local_skin_friction_coefficient", "wall_heat_flux"]:
        if wall is None and field == "wall_heat_flux":
            continue
        comparison = 100.0 * (getattr(transition, field) / getattr(flat, field) - 1.0)
        print(
            f"T4 M={mach} Tw={wall} {field}: "
            f"at s=.2,1,10 m = {comparison[np.isin(x, [0.2, 1, 10])].tolist()} %"
        )


def refined_march(
    edge: BoundaryLayerEdge,
    targets: FloatArray,
    *,
    normal: bool = False,
    relative_step: float = 0.03,
) -> MarchingBoundaryLayerResult:
    interp = _EdgeInterpolant(
        edge, np.full(edge.distance.size, 300.0), 3.0, transition_distance=0.1
    )
    eta = eta_grid(0.00035, 1.02, 0.5, 150.0) if normal else _turbulent_eta_grid(2e7)
    solver = sa_module.TurbulentDiscretization(eta, 0.72, 0.9, edge.viscosity_model)
    return _result(
        edge,
        interp,
        march(
            interp.state,
            targets,
            solver,
            relative_step=relative_step,
            internal_stations=(0.1,),
        ),
        0.72,
    )


@pytest.mark.solver_verification
def test_t5_streamwise_and_wall_normal_refinement(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    x = np.unique(
        np.r_[
            np.geomspace(0.01, 2, 180),
            0.1,
            0.1001,
            0.101,
            0.102,
            0.105,
            0.11,
            0.12,
            0.15,
            0.2,
            0.3,
            0.5,
            1,
        ]
    )
    edge = constant_edge(2.0, 220.0, x.tolist())
    base = refined_march(edge, x)
    normal = refined_march(edge, x, normal=True)
    fine = np.unique(
        np.concatenate(
            [
                np.linspace(left, right, 9)[1:]
                for left, right in zip(
                    edge.distance[:-1], edge.distance[1:], strict=True
                )
            ]
        )
    )
    tighten_streamwise_control(monkeypatch)
    stream = refined_march(edge, fine, relative_step=0.0075)
    for label, result, ix, cf_limit in [
        ("normal", normal, np.arange(x.size), 0.003),
        ("stream", stream, np.searchsorted(fine, x), 0.001),
    ]:
        differences = []
        for field in ["local_skin_friction_coefficient", "wall_heat_flux"]:
            difference = abs(getattr(base, field) / getattr(result, field)[ix] - 1.0)
            near = float(np.max(difference[(x > 0.1) & (x < 0.2)]))
            far = float(np.max(difference[x >= 0.2]))
            print(f"T5 {label} {field}: near={100 * near:.9g}%, far={100 * far:.9g}%")
            differences.append(far)
        assert differences[0] < cf_limit
        assert differences[1] < 0.005


@pytest.mark.solver_verification
@pytest.mark.parametrize(("start", "end"), [(2.0, 3.0), (3.0, 2.5)])
@pytest.mark.parametrize("cone", [False, True])
def test_t6_pressure_gradient_and_cone_complete(
    start: float, end: float, cone: bool
) -> None:
    edge = smooth_edge(start, end)
    radius = edge.distance * np.sin(np.deg2rad(10.0)) if cone else None
    try:
        laminar = marching_boundary_layer(
            edge, regime=LAMINAR, wall_temperature=300.0, wall_radius=radius
        )
    except BoundaryLayerSeparationError as separation:
        # Compression can separate the laminar reference downstream of the trip.
        # Its upstream result still supplies the same input interpolation/history.
        laminar = separation.result
    transition = marching_boundary_layer(
        edge,
        regime=TRANSITIONAL,
        wall_temperature=300.0,
        wall_radius=radius,
        transition_distance=0.5,
    )
    assert transition.distance[-1] == edge.distance[-1]
    upstream = transition.distance <= 0.5
    reference = np.searchsorted(laminar.distance, transition.distance[upstream])
    np.testing.assert_array_equal(
        laminar.distance[reference], transition.distance[upstream]
    )
    for field in ["local_skin_friction_coefficient", "wall_heat_flux"]:
        error = relative_difference(
            getattr(transition, field)[upstream], getattr(laminar, field)[reference]
        )
        print(f"T6 M{start}->{end} cone={cone} {field}={100 * error:.9g}%")
        assert error < 0.002  # Numerical proximity, not an accuracy claim.


@pytest.mark.solver_verification
@pytest.mark.parametrize("before_first_output", [False, True])
def test_t8_transition_on_and_between_output_stations(
    monkeypatch: pytest.MonkeyPatch,
    before_first_output: bool,
) -> None:
    x = np.unique(np.r_[np.geomspace(0.01, 2, 140), 0.1001, 0.101, 0.11, 0.2, 0.5, 1])
    if before_first_output:
        x = x[x > 0.1]
    accepted_xi: list[float] = []
    original = march_module._step_is_resolved

    def observe(
        previous: march_module.Station,
        station: march_module.Station,
        state: StationEdge,
    ) -> bool:
        resolved = original(previous, station, state)
        if resolved:
            accepted_xi.append(state.xi)
        return resolved

    monkeypatch.setattr(march_module, "_step_is_resolved", observe)
    y = np.union1d(x, [0.1])
    a = marching_boundary_layer(
        constant_edge(2.0, 220.0, x.tolist()),
        regime=TRANSITIONAL,
        wall_temperature=300.0,
        transition_distance=0.1,
    )
    interp = _EdgeInterpolant(constant_edge(2.0, 220.0, x.tolist()), None)
    assert interp.xi(0.1) in accepted_xi
    b = marching_boundary_layer(
        constant_edge(2.0, 220.0, y.tolist()),
        regime=TRANSITIONAL,
        wall_temperature=300.0,
        transition_distance=0.1,
    )
    np.testing.assert_array_equal(a.distance, x)
    assert 0.1 not in a.distance and 0.1 in b.distance
    ix = np.searchsorted(y, x)
    for field in ["local_skin_friction_coefficient", "wall_heat_flux"]:
        error = relative_difference(getattr(a, field), getattr(b, field)[ix])
        print(f"T8 {field}={100 * error:.9g}%")
        assert error < 1e-8
