"""Exact derivatives of the SA discrete residual, including its active branches."""

import warnings
from importlib import import_module

import numpy as np
import pytest
from test_marching_boundary_layer import plate_edge
from test_marching_boundary_layer_jacobian import (
    MODELS,
    assert_results_close,
)
from test_marching_boundary_layer_jacobian import (
    SA_ETA as ETA,
)

from aerophysics import (
    BoundaryLayerRegime,
    MarchingBoundaryLayerResult,
    marching_boundary_layer,
)
from aerophysics import _boundary_layer_sa as sa
from aerophysics._array import FloatArray
from aerophysics._boundary_layer_march import StationEdge, eta_grid
from aerophysics._boundary_layer_sa import TurbulentDiscretization
from aerophysics.exceptions import ApplicabilityWarning
from aerophysics.transport import (
    AIR_VISCOSITY,
    DynamicViscosityModel,
)


def dense(bands: FloatArray, band: int = 9) -> FloatArray:
    size = bands.shape[1]
    rows = np.arange(size)[None, :] + np.arange(3 * band + 1)[:, None] - 2 * band
    columns = np.broadcast_to(np.arange(size), rows.shape)
    valid = (rows >= 0) & (rows < size)
    result = np.zeros((size, size))
    result[rows[valid], columns[valid]] = bands[valid]
    return result


@pytest.mark.parametrize(
    ("branch_index", "branch"),
    list(enumerate(["regular", "limited", "capped", "negative"])),
    ids=["regular", "limited", "capped", "negative"],
)
@pytest.mark.parametrize(
    ("state_index", "xi", "beta", "heating", "coefficient", "production", "growth"),
    [
        (0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (1, 0.37, 0.23, 1.8, 1.3, 1.0, 0.2),
        (2, 2.1, -0.15, 4.0, 23.0, 0.4, -0.17),
    ],
    ids=["starting", "favorable", "adverse"],
)
def test_sa_analytic_matches_both_difference_jacobians(
    branch_index: int,
    branch: str,
    state_index: int,
    xi: float,
    beta: float,
    heating: float,
    coefficient: float,
    production: float,
    growth: float,
) -> None:
    model = MODELS[(branch_index + state_index) % 2]
    wall = (None, 1.4)[state_index % 2]
    solver = TurbulentDiscretization(ETA, 0.72, 0.9, model)
    edge = StationEdge(
        xi, 1.0, beta, heating, 290.0, wall, 2.8, growth, 3.0, production
    )
    rng = np.random.default_rng(173 + int(xi * 10))
    state = rng.uniform(0.6, 1.4, (5, len(ETA)))
    state[1] = (0.4 * ETA + rng.uniform(-0.02, 0.02, len(ETA))) * (
        -1.0 if xi > 1.0 else 1.0
    )
    state[3] = rng.uniform(0.3, 0.7, len(ETA))
    state[4] = rng.uniform(2.0, 3.0, len(ETA))
    if branch in ("limited", "capped"):
        state[2] = rng.uniform(0.95, 1.05, len(ETA))
        state[3] = rng.uniform(3.0, 3.3, len(ETA))
        c = solver.chapman_rubesin(state[2], edge.temperature)
        chi = state[3] / (c * state[2] ** 2)
        fv2 = 1.0 - chi / (1.0 + chi * sa._fv1(chi))
        vort = edge.reynolds_root * np.abs(solver.derivative(state[1])) / state[2, 1:-1]
        factor = 0.72 if branch == "limited" else 5.0
        state[4, 1:-1] = np.sqrt(
            -state[3, 1:-1] * fv2[1:-1] / (sa._KAPPA**2 * vort * factor)
        )
    elif branch == "negative":
        state[3, 1::3] = -0.2
    inner_g = state[2, 1:-1]
    c = solver.chapman_rubesin(state[2], edge.temperature)[1:-1]
    vort = edge.reynolds_root * np.abs(solver.derivative(state[1])) / inner_g
    modified, extra, _ = sa._source_terms(
        vort, np.maximum(state[3, 1:-1], 0.0), c, inner_g, state[4, 1:-1]
    )
    r = np.maximum(state[3, 1:-1], 0.0) / (
        modified * sa._KAPPA**2 * state[4, 1:-1] ** 2
    )
    assert np.all((extra >= -sa._CV2 * vort) == (branch not in ("limited", "capped")))
    assert np.all(r > sa._MAX_R) if branch == "capped" else np.all(r < sa._MAX_R)
    lagged = rng.uniform(-0.3, 0.3, state.shape)
    residual = solver.residual(state, edge, coefficient, lagged)
    analytic = solver.jacobian(state, residual, edge, coefficient, lagged)
    finite = solver._finite_difference_jacobian(
        state, residual, edge, coefficient, lagged
    )
    assert analytic.shape == (28, state.size)
    assert analytic.flags.f_contiguous
    np.testing.assert_array_equal(analytic != 0.0, finite != 0.0)
    matrix, forward = dense(analytic), dense(finite)
    scale = np.max(np.abs(forward), axis=1, keepdims=True)
    # sqrt(eps) forward steps lose digits to residual cancellation; use the
    # same equation-row scaling as the laminar test, including off-stencil zeros.
    assert np.all(np.abs(matrix - forward) <= 2e-6 * scale)
    central = np.empty_like(matrix)
    for node in range(len(ETA)):
        for variable in range(5):
            step = 2e-6 * max(abs(state[variable, node]), 1.0)
            plus, minus = state.copy(), state.copy()
            plus[variable, node] += step
            minus[variable, node] -= step
            central[:, 5 * node + variable] = (
                solver.residual(plus, edge, coefficient, lagged)
                - solver.residual(minus, edge, coefficient, lagged)
            ).T.ravel() / (2.0 * step)
    # Central steps trade O(h^2) truncation against O(eps/h) cancellation.
    assert np.all(np.abs(matrix - central) <= 2e-8 * scale)


@pytest.mark.parametrize("case", ["pressure", "transitional"])
def test_sa_complete_marches_match_finite_differences(
    monkeypatch: pytest.MonkeyPatch,
    case: str,
) -> None:
    # Short, low-Re marches on a compact test grid keep both comparisons in
    # the normal tier. This tests equality of the discrete solves, not grid
    # accuracy; solver_verification retains the production grids and high Re.
    edge = plate_edge([0.0, 0.0001, 0.00011, 0.0002])
    if case == "pressure":
        edge = type(edge).from_pressure(
            edge.distance,
            [20000.0, 20000.0, 19900.0, 19800.0],
            total_pressure=1.6e5,
            total_temperature=500.0,
        )
    solvers: list[TurbulentDiscretization] = []

    class CountingSolver(TurbulentDiscretization):
        def __init__(
            self,
            eta: FloatArray,
            prandtl_number: float,
            turbulent_prandtl_number: float,
            viscosity_model: DynamicViscosityModel,
        ) -> None:
            super().__init__(
                eta, prandtl_number, turbulent_prandtl_number, viscosity_model
            )
            solvers.append(self)

    module = import_module("aerophysics.marching_boundary_layer")
    monkeypatch.setattr(module, "TurbulentDiscretization", CountingSolver)
    monkeypatch.setattr(
        module, "_turbulent_eta_grid", lambda _: eta_grid(0.01, 1.25, 10.0, 150.0)
    )

    def run() -> MarchingBoundaryLayerResult:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ApplicabilityWarning)
            return marching_boundary_layer(
                edge,
                regime=(
                    BoundaryLayerRegime.TRANSITIONAL
                    if case == "transitional"
                    else BoundaryLayerRegime.TURBULENT
                ),
                wall_temperature=300.0,
                transition_distance=0.0001 if case == "transitional" else None,
            )

    analytic = run()
    analytic_picard = sum(s.picard_fallbacks for s in solvers)
    solvers.clear()
    with monkeypatch.context() as patch:
        patch.setattr(TurbulentDiscretization, "_analytic_jacobian_enabled", False)
        finite = run()
    finite_picard = sum(s.picard_fallbacks for s in solvers)
    assert_results_close(analytic, finite)
    assert analytic_picard == finite_picard


def test_sa_abs_and_positive_n_ties_match_forward_columns() -> None:
    solver = TurbulentDiscretization(ETA, 0.72, 0.9, AIR_VISCOSITY)
    edge = StationEdge(0.2, 1.0, 0.1, 0.5, 290.0, None, 2.8, 0.2, 0.3)
    state = solver.starting_guess(edge)
    state[1] = 0.0  # Exactly zero vorticity, with both signs of stencil weights.
    state[3] = 0.3
    state[3, 3] = 0.0
    state[3, 5] = -0.2
    state[4] = 2.0
    lagged = -0.8 * state
    residual = solver.residual(state, edge, 1.3, lagged)
    analytic = dense(solver._analytic_jacobian(state, edge, 1.3, lagged))
    scale = np.max(np.abs(analytic), axis=1)
    for node in range(len(ETA)):
        for variable in range(5):
            step = 1e-7
            plus = state.copy()
            plus[variable, node] += step
            forward = (
                solver.residual(plus, edge, 1.3, lagged) - residual
            ).T.ravel() / step
            assert np.all(
                np.abs(analytic[:, 5 * node + variable] - forward) < 2e-5 * scale
            )
