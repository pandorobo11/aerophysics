"""Exact discrete linearization and full-march comparisons with colored differences."""

import warnings
from dataclasses import fields, replace

import numpy as np
import pytest
from test_marching_boundary_layer import (
    ETA,
    LinearViscosity,
    howarth_edge,
    plate_edge,
    smooth_edge,
)

from aerophysics import (
    BoundaryLayerRegime,
    MarchingBoundaryLayerResult,
    marching_boundary_layer,
)
from aerophysics._array import FloatArray
from aerophysics._boundary_layer_march import (
    Discretization,
    StationEdge,
    _viscosity_with_derivative,
)
from aerophysics._boundary_layer_sa import TurbulentDiscretization
from aerophysics.exceptions import ApplicabilityWarning, BoundaryLayerSeparationError
from aerophysics.transport import (
    AIR_KEYES_VISCOSITY,
    AIR_VISCOSITY,
    BlottnerModel,
    DynamicViscosityModel,
    KeyesModel,
    SutherlandModel,
    WilkeMixtureViscosityModel,
)

MODELS = [AIR_VISCOSITY, AIR_KEYES_VISCOSITY]


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.parametrize("wall", [None, 1.4])
@pytest.mark.parametrize("smooth", [False, True])
@pytest.mark.parametrize(
    ("xi", "beta", "heating", "coefficient"),
    [(0.0, 0.0, 0.0, 0.0), (0.37, 0.23, 1.8, 1.3), (2.1, -0.15, 4.0, 23.0)],
)
def test_analytic_matches_colored_and_independent_differences(
    model: SutherlandModel | KeyesModel,
    wall: float | None,
    smooth: bool,
    xi: float,
    beta: float,
    heating: float,
    coefficient: float,
) -> None:
    solver = Discretization(ETA, 0.72, model)
    edge = StationEdge(xi, 1.0, beta, heating, 290.0, wall)
    rng = np.random.default_rng(1123)
    state = solver.starting_guess(edge)
    if smooth:
        state += rng.uniform(-0.02, 0.02, (3, 1)) * np.sin(ETA)
        state += rng.uniform(-0.02, 0.02, (3, 1)) * (1.0 - np.exp(-ETA))
    else:
        state += rng.uniform(-0.02, 0.02, state.shape)
    lagged = -0.8 * state + rng.uniform(-0.03, 0.03, state.shape)
    residual = solver.residual(state, edge, coefficient, lagged)
    analytic = solver.jacobian(state, residual, edge, coefficient, lagged)
    solver._analytic_jacobian_enabled = False
    finite = solver.jacobian(state, residual, edge, coefficient, lagged)
    assert analytic.flags.f_contiguous
    assert analytic.shape == (16, state.size)
    # Scale by equation rows, not storage diagonals: near the wall, cancellation
    # in the O(h^-2) fluxes contaminates small f/g couplings in forward differences.
    matrix_rows = np.arange(state.size)[None, :] + np.arange(16)[:, None] - 10
    valid = (matrix_rows >= 0) & (matrix_rows < state.size)
    scale = np.zeros(state.size)
    np.maximum.at(scale, matrix_rows[valid], np.abs(finite[valid]))
    assert np.all(
        np.abs(analytic[valid] - finite[valid]) <= 1e-6 * scale[matrix_rows[valid]]
    )
    np.testing.assert_array_equal(analytic != 0.0, finite != 0.0)

    # Independently perturb individual columns at both boundaries and inside.
    # This checks band indexing and off-stencil zeros without colored assembly.
    for node in (0, 1, len(ETA) // 2, len(ETA) - 2, len(ETA) - 1):
        for variable in range(3):
            step = 1e-5 * max(abs(state[variable, node]), 1.0)
            plus, minus = state.copy(), state.copy()
            plus[variable, node] += step
            minus[variable, node] -= step
            central = (
                solver.residual(plus, edge, coefficient, lagged)
                - solver.residual(minus, edge, coefficient, lagged)
            ).T.ravel() / (2.0 * step)
            column = 3 * node + variable
            expected = np.zeros(state.size)
            rows = np.arange(max(0, column - 5), min(state.size, column + 6))
            expected[rows] = analytic[10 + rows - column, column]
            np.testing.assert_allclose(
                expected, central, rtol=1e-8, atol=1e-8 * np.max(np.abs(central))
            )


@pytest.mark.parametrize("model", MODELS)
def test_viscosity_derivative_retains_values_validation_and_warnings(
    model: SutherlandModel | KeyesModel,
) -> None:
    temperature = np.array([79.0, 290.0, 500.0, 1845.0])
    value, derivative = _viscosity_with_derivative(model, temperature)
    np.testing.assert_array_equal(value, model.dynamic_viscosity(temperature))
    step = temperature * 1e-5
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ApplicabilityWarning)
        central = (
            model.dynamic_viscosity(temperature + step)
            - model.dynamic_viscosity(temperature - step)
        ) / (2.0 * step)
    np.testing.assert_allclose(derivative, central, rtol=1e-9)
    for bad in (0.0, -1.0, np.nan, np.inf):
        with pytest.raises(ValueError):
            _viscosity_with_derivative(model, np.array([bad]))
    if isinstance(model, KeyesModel):
        with pytest.warns(ApplicabilityWarning, match="fitted range"):
            _viscosity_with_derivative(model, np.array([1900.0]))


class SutherlandSubclass(SutherlandModel):
    """A subclass must not silently acquire the built-in analytic derivative."""


class KeyesSubclass(KeyesModel):
    """A subclass must retain the finite-difference dispatch as well."""


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
def test_unsupported_models_keep_finite_differences(
    model: DynamicViscosityModel,
) -> None:
    solver = Discretization(ETA, 0.72, model)
    edge = StationEdge(0.37, 1.0, 0.23, 1.8, 290.0, None)
    state = solver.starting_guess(edge)
    lagged = -0.8 * state
    residual = solver.residual(state, edge, 1.3, lagged)
    np.testing.assert_array_equal(
        solver.jacobian(state, residual, edge, 1.3, lagged),
        solver._finite_difference_jacobian(state, residual, edge, 1.3, lagged),
    )


def test_sa_explicitly_keeps_finite_differences(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    solver = TurbulentDiscretization(ETA, 0.72, 0.9, AIR_VISCOSITY)
    assert not solver._analytic_jacobian_enabled
    edge = StationEdge(0.37, 1.0, 0.23, 1.8, 290.0, None, 1000.0, 0.1, 3.0)
    state = solver.starting_guess(edge)
    lagged = -0.8 * state
    residual = solver.residual(state, edge, 1.3, lagged)
    expected = solver._finite_difference_jacobian(state, residual, edge, 1.3, lagged)
    np.testing.assert_array_equal(
        solver.jacobian(state, residual, edge, 1.3, lagged), expected
    )
    # Even a mistakenly enabled inherited switch cannot route five unknowns
    # through the laminar analytic assembly.
    monkeypatch.setattr(solver, "_analytic_jacobian_enabled", True)
    np.testing.assert_array_equal(
        solver.jacobian(state, residual, edge, 1.3, lagged), expected
    )


def assert_results_close(
    actual: MarchingBoundaryLayerResult, expected: MarchingBoundaryLayerResult
) -> None:
    for field in fields(actual):
        value, reference = getattr(actual, field.name), getattr(expected, field.name)
        if reference is None:
            assert value is None
        else:
            np.testing.assert_allclose(
                value, reference, rtol=1e-8, atol=1e-10, err_msg=field.name
            )


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.parametrize(
    "case", ["plate", "pressure", "wall-change", "adiabatic", "howarth"]
)
def test_complete_marches_match_finite_differences(
    monkeypatch: pytest.MonkeyPatch, model: SutherlandModel | KeyesModel, case: str
) -> None:
    wall: float | FloatArray | None = 300.0
    edge = plate_edge([0.0, 0.5, 0.51, 1.0])
    if case == "pressure":
        edge = smooth_edge(2.0, 3.0, count=41)
    elif case == "wall-change":
        wall = np.array([300.0, 300.0, 400.0, 400.0])
    elif case == "adiabatic":
        wall = None
    elif case == "howarth":
        edge, _ = howarth_edge(0.01, count=41)
        wall = None
    edge = replace(edge, viscosity_model=model)

    def run() -> tuple[MarchingBoundaryLayerResult, float | None]:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ApplicabilityWarning)
            try:
                result = marching_boundary_layer(
                    edge, regime=BoundaryLayerRegime.LAMINAR, wall_temperature=wall
                )
            except BoundaryLayerSeparationError as error:
                return error.result, error.separation_distance
        return result, None

    analytic, separation = run()
    monkeypatch.setattr(Discretization, "_analytic_jacobian_enabled", False)
    finite, reference_separation = run()
    assert_results_close(analytic, finite)
    if case == "howarth":
        assert separation is not None and reference_separation is not None
        assert separation == pytest.approx(reference_separation, rel=1e-8, abs=1e-10)
    else:
        assert separation is reference_separation is None
