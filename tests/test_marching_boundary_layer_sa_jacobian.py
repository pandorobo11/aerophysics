"""Exact derivatives of the SA discrete residual, including its active branches."""

import warnings
from importlib import import_module

import numpy as np
import pytest
from test_marching_boundary_layer import LinearViscosity, plate_edge
from test_marching_boundary_layer_jacobian import (
    MODELS,
    KeyesSubclass,
    SutherlandSubclass,
    assert_results_close,
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
    BlottnerModel,
    DynamicViscosityModel,
    KeyesModel,
    SutherlandModel,
    WilkeMixtureViscosityModel,
)

# Small, nonuniform grid: every column can be checked independently, including
# both boundaries, without the cancellation of tiny near-wall verification cells.
ETA = np.array([0.0, 0.17, 0.38, 0.64, 0.95, 1.32, 1.76, 2.27, 2.85])


def dense(bands: FloatArray, band: int = 9) -> FloatArray:
    size = bands.shape[1]
    rows = np.arange(size)[None, :] + np.arange(3 * band + 1)[:, None] - 2 * band
    columns = np.broadcast_to(np.arange(size), rows.shape)
    valid = (rows >= 0) & (rows < size)
    result = np.zeros((size, size))
    result[rows[valid], columns[valid]] = bands[valid]
    return result


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.parametrize("wall", [None, 1.4])
@pytest.mark.parametrize("branch", ["regular", "limited", "capped", "negative"])
@pytest.mark.parametrize(
    ("xi", "beta", "heating", "coefficient", "production", "growth"),
    [
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.37, 0.23, 1.8, 1.3, 1.0, 0.2),
        (2.1, -0.15, 4.0, 23.0, 0.4, -0.17),
    ],
)
def test_sa_analytic_matches_both_difference_jacobians(
    model: SutherlandModel | KeyesModel,
    wall: float | None,
    branch: str,
    xi: float,
    beta: float,
    heating: float,
    coefficient: float,
    production: float,
    growth: float,
) -> None:
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


@pytest.mark.parametrize(
    "model",
    [
        BlottnerModel(0.026, 0.3, -11.3, temperature_range=(79.0, 30000.0)),
        LinearViscosity(),
        SutherlandSubclass(1.7894e-5, 288.15, 110.4),
        KeyesSubclass(1.488e-6, 122.1, 5.0),
        WilkeMixtureViscosityModel((AIR_VISCOSITY,), (0.028,), (1.0,)),
    ],
)
def test_sa_unsupported_models_keep_finite_differences(
    model: DynamicViscosityModel,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    solver = TurbulentDiscretization(ETA, 0.72, 0.9, model)
    edge = StationEdge(0.37, 1.0, 0.23, 1.8, 290.0, None, 1000.0, 0.1, 3.0)
    state = solver.starting_guess(edge)
    lagged = -0.8 * state
    residual = solver.residual(state, edge, 1.3, lagged)

    def forbidden(*args: object) -> FloatArray:
        pytest.fail("unsupported model reached the analytic assembly")

    monkeypatch.setattr(solver, "_analytic_jacobian", forbidden)
    np.testing.assert_array_equal(
        solver.jacobian(state, residual, edge, 1.3, lagged),
        solver._finite_difference_jacobian(state, residual, edge, 1.3, lagged),
    )


@pytest.mark.parametrize("model", MODELS)
def test_sa_analytic_preserves_viscosity_validation_and_warnings(
    model: SutherlandModel | KeyesModel,
) -> None:
    solver = TurbulentDiscretization(ETA, 0.72, 0.9, model)
    edge = StationEdge(0.0, 1.0, 0.0, 0.0, 290.0, 1.0, 2.0, 0.0, 3.0)
    state = solver.starting_guess(edge)
    for bad in (0.0, -1.0, np.nan, np.inf):
        state[2, 2] = bad
        with pytest.raises(ValueError):
            solver._analytic_jacobian(state, edge, 0.0, np.zeros_like(state))
    if isinstance(model, KeyesModel):
        state[2, 2] = 1900.0 / edge.temperature
        with pytest.warns(ApplicabilityWarning, match="fitted range"):
            solver._analytic_jacobian(state, edge, 0.0, np.zeros_like(state))


@pytest.mark.parametrize(
    "case",
    ["plate", "pressure", "wall-change", "adiabatic", "axisymmetric", "transitional"],
)
def test_sa_complete_marches_match_finite_differences(
    monkeypatch: pytest.MonkeyPatch,
    case: str,
) -> None:
    # Short, low-Re marches on a compact test grid keep all six comparisons in
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
    wall: float | FloatArray | None = 300.0
    if case == "wall-change":
        wall = np.array([300.0, 300.0, 310.0, 310.0])
    elif case == "adiabatic":
        wall = None
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
                wall_temperature=wall,
                wall_radius=0.05 if case == "axisymmetric" else None,
                transition_distance=0.0001 if case == "transitional" else None,
            )

    analytic = run()
    analytic_picard = sum(s.picard_fallbacks for s in solvers)
    solvers.clear()
    with monkeypatch.context() as patch:
        patch.setattr(TurbulentDiscretization, "_analytic_jacobian_enabled", False)
        finite = run()
    finite_picard = sum(s.picard_fallbacks for s in solvers)
    print(f"{case}: picard_fallbacks FD={finite_picard}, analytic={analytic_picard}")
    assert_results_close(analytic, finite)
    assert analytic_picard == finite_picard


@pytest.mark.parametrize("floor", [1e-30, sa._MINIMUM_VORTICITY])
def test_floor_uses_positive_column_direction_at_ties(floor: float) -> None:
    values = np.array([0.5, 1.0, 2.0]) * floor
    directions = np.array([[-2.0], [0.0], [3.0]]) * floor
    np.testing.assert_array_equal(
        sa._floor_derivative(values, directions, floor),
        np.array([[0.0, 0.0, -2.0], [0.0, 0.0, 0.0], [0.0, 3.0, 3.0]]) * floor,
    )


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


def test_source_join_floor_and_cap_directions() -> None:
    # Positive-vorticity join: both branches meet with the same first derivative.
    n, c, g, y = (np.array([value]) for value in (2.5, 1.0, 1.0, 1.0))
    extra = sa._source_terms(np.ones(1), n, c, g, y)[1]
    vorticity = -extra / sa._CV2
    assert extra == -sa._CV2 * vorticity
    basis = np.eye(5)[:, :, None]
    directions = (basis[0], basis[1], basis[2], basis[3], basis[4])
    modified, fw, dm, dfw = sa._source_derivatives(vorticity, n, c, g, y, directions)
    expected = sa._source_terms(vorticity, n, c, g, y)
    np.testing.assert_array_equal(modified, expected[0])
    np.testing.assert_array_equal(fw, expected[2])
    inputs = np.array([vorticity, n, c, g, y])
    for variable in range(5):
        step = 1e-6 * max(inputs[variable, 0], 1.0)
        plus, minus = inputs.copy(), inputs.copy()
        plus[variable] += step
        minus[variable] -= step
        upper = sa._source_terms(*plus)
        lower = sa._source_terms(*minus)
        np.testing.assert_allclose(
            dm[variable], (upper[0] - lower[0]) / (2.0 * step), rtol=1e-5
        )
        np.testing.assert_allclose(
            dfw[variable], (upper[2] - lower[2]) / (2.0 * step), rtol=1e-5
        )

    # Vorticity floor at an exact tie: positive directions leave the floor,
    # negative directions remain on it. n=0 makes the r denominator floored too.
    zero, one = np.zeros(1), np.ones(1)
    dv = np.array([[1.0], [-1.0]]) * sa._MINIMUM_VORTICITY
    zeros = np.zeros_like(dv)
    _, _, dm, _ = sa._source_derivatives(
        one * sa._MINIMUM_VORTICITY,
        zero,
        one,
        one,
        one * 1e-12,
        (dv, zeros, zeros, zeros, zeros),
    )
    np.testing.assert_array_equal(dm, np.array([[sa._MINIMUM_VORTICITY], [0.0]]))
    # Exercise the limited-branch denominator below its floor as well.
    _, _, dm, _ = sa._source_derivatives(
        one * 1e-32,
        one * 3e-32,
        one * 1e-32,
        one,
        one,
        (dv, zeros, zeros, zeros, zeros),
    )
    np.testing.assert_array_equal(dm, zeros)

    # With the r denominator floored, n=MAX_R*floor is an exact cap tie.
    # fw is numerically saturated here: its tiny inward derivative cannot be
    # recovered by subtracting rounded fw values, so check the analytic formula.
    n = one * (sa._MAX_R * 1e-30)
    dn = np.array([[1.0], [-1.0]]) * n
    _, _, _, dfw = sa._source_derivatives(
        zero,
        n,
        n,
        one,
        one * 1e-12,
        (zeros, dn, zeros, zeros, zeros),
    )
    r = sa._MAX_R
    r6 = r + sa._CW2 * (r**6 - r)
    factor = ((1.0 + sa._CW3**6) / (r6**6 + sa._CW3**6)) ** (1.0 / 6.0)
    inward = (
        -r
        * factor
        * sa._CW3**6
        / (r6**6 + sa._CW3**6)
        * (1.0 + sa._CW2 * (6.0 * r**5 - 1.0))
    )
    assert dfw[0, 0] == 0.0
    assert dfw[1, 0] == pytest.approx(inward, rel=1e-14, abs=0.0)
    _, _, _, dfw = sa._source_derivatives(
        zero,
        2.0 * n,
        n,
        one,
        one * 1e-12,
        (zeros, dn, zeros, zeros, zeros),
    )
    np.testing.assert_array_equal(dfw, zeros)
