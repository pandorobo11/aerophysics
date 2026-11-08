"""Local SA complementarity differs from the projected Picard solve."""

from dataclasses import fields, replace

import numpy as np
import pytest
from test_marching_boundary_layer import smooth_edge
from test_marching_boundary_layer_sa import constant_edge

from aerophysics import (
    BoundaryLayerRegime,
    MarchingBoundaryLayerResult,
    marching_boundary_layer,
)
from aerophysics import _boundary_layer_sa as sa_module
from aerophysics._array import FloatArray
from aerophysics._boundary_layer_march import Discretization, Station, StationEdge
from aerophysics._boundary_layer_sa import TurbulentDiscretization
from aerophysics.transport import AIR_VISCOSITY

StationInput = tuple[
    TurbulentDiscretization, FloatArray, StationEdge, float, FloatArray
]


@pytest.fixture(scope="module")
def clipped_station() -> StationInput:
    """Capture real march arguments where legacy Newton rejects its first step."""
    original = TurbulentDiscretization.newton
    original_damping = TurbulentDiscretization.damping
    calls = 0
    captured: StationInput | None = None

    class FoundStationError(Exception):
        pass

    def damping(
        self: TurbulentDiscretization, state: FloatArray, update: FloatArray
    ) -> float:
        nonlocal calls
        calls += 1
        return original_damping(self, state, update)

    def newton(
        self: TurbulentDiscretization,
        guess: FloatArray,
        edge: StationEdge,
        coefficient: float,
        lagged: FloatArray,
    ) -> Station | None:
        nonlocal calls, captured
        calls = 0
        station = original(self, guess, edge, coefficient, lagged)
        if station is None and calls == 1:
            captured = self, guess.copy(), edge, coefficient, lagged.copy()
            raise FoundStationError
        return station

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(TurbulentDiscretization, "_bound_constrained_newton", False)
        patch.setattr(TurbulentDiscretization, "damping", damping)
        patch.setattr(TurbulentDiscretization, "newton", newton)
        with pytest.raises(FoundStationError):
            marching_boundary_layer(
                smooth_edge(2.0, 3.0),
                regime=BoundaryLayerRegime.TURBULENT,
                wall_temperature=300.0,
            )
    assert captured is not None
    return captured


@pytest.mark.solver_verification
def test_first_update_failure_becomes_local_complementarity(
    clipped_station: StationInput, monkeypatch: pytest.MonkeyPatch
) -> None:
    solver, guess, edge, coefficient, lagged = clipped_station
    monkeypatch.setattr(solver, "_bound_constrained_newton", False)
    assert solver.newton(guess, edge, coefficient, lagged) is None
    monkeypatch.setattr(solver, "_bound_constrained_newton", True)
    actual = solver.newton(guess, edge, coefficient, lagged)
    reference = solver.picard(guess, edge, coefficient, lagged)
    assert actual is not None and reference is not None
    # Observed maxima (f,u,g,n,Y): 2.19e-10, 1.69e-9, 3.16e-10,
    # 1.26e-5, 8.88e-11. The n discrepancy is a different fixed point,
    # not just Picard's 1e-9 stopping tolerance.
    for variable, tolerance in enumerate([1e-9, 4e-9, 1e-9, 3e-5, 5e-10]):
        np.testing.assert_allclose(
            actual.state[variable], reference.state[variable], rtol=0, atol=tolerance
        )
    residual = solver.residual(actual.state, edge, coefficient, lagged)
    bands = solver.jacobian(actual.state, residual, edge, coefficient, lagged)
    diagonal = bands[2 * solver.band, 8:-5:5]
    n = actual.state[3, 1:-1]
    active = n == 0.0
    assert np.any(active)
    assert np.all(n >= 0.0)
    assert np.all(diagonal < 0.0)
    assert np.all(residual[3, 1:-1][active] <= 0.0)
    np.testing.assert_allclose(
        residual[3, 1:-1][~active] / diagonal[~active], 0.0, atol=1e-8
    )
    # Picard leaves a finite residual at positive neighbors of its clipped node.
    picard_residual = solver.residual(reference.state, edge, coefficient, lagged)
    assert np.max(np.abs(picard_residual[3, 1:-1][~active])) > 0.05
    expected = actual.state.copy()
    solver.enforce_boundaries(expected, edge)
    np.testing.assert_array_equal(actual.state, expected)


@pytest.mark.solver_verification
@pytest.mark.parametrize(
    "case", ["zpg", "pressure", "adiabatic", "axisymmetric", "high_reynolds"]
)
def test_complete_marches_and_legacy_switch(
    case: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    edge = (
        smooth_edge(2.0, 3.0)
        if case in ("pressure", "adiabatic")
        else constant_edge(2.0, 220.0, [0.1, 1.0])
    )
    if case == "high_reynolds":
        edge = constant_edge(2.0, 220.0, [100.0], unit_reynolds=1e8)
    wall = None if case == "adiabatic" else 300.0
    radius = 0.1 + 0.01 * edge.distance if case == "axisymmetric" else None
    instances: list[TurbulentDiscretization] = []
    original = TurbulentDiscretization.solve

    def solve(
        self: TurbulentDiscretization,
        guess: FloatArray,
        at: StationEdge,
        coefficient: float,
        lagged: FloatArray,
    ) -> Station | None:
        if not instances or instances[-1] is not self:
            instances.append(self)
        return original(self, guess, at, coefficient, lagged)

    monkeypatch.setattr(TurbulentDiscretization, "solve", solve)

    def run() -> MarchingBoundaryLayerResult:
        return marching_boundary_layer(
            edge,
            regime=BoundaryLayerRegime.TURBULENT,
            wall_temperature=wall,
            wall_radius=radius,
        )

    monkeypatch.setattr(TurbulentDiscretization, "_bound_constrained_newton", False)
    old = run()
    old_fallbacks = instances[-1].picard_fallbacks
    monkeypatch.setattr(TurbulentDiscretization, "_bound_constrained_newton", True)
    new = run()
    new_fallbacks = instances[-1].picard_fallbacks
    for field in fields(old):
        reference, actual = getattr(old, field.name), getattr(new, field.name)
        # Measured maximum |new-old| / max|old|: n = 2.111e-6 (adiabatic),
        # all other fields <= 7.072e-8. Picard stops at a 1e-9 change, but
        # fixed-point differences and step acceptance amplify march differences.
        # Field scaling also covers near-zero SA values without hiding them in rtol.
        tolerance = 4e-6 if field.name == "nu_tilde_ratio" else 2e-7
        if case == "high_reynolds":
            # With the wall guard at Re=1e10: n = 6.031e-5, u = 6.034e-8,
            # g = 8.184e-8, other fields <= 7.222e-7 (delta99), scaled by maxima.
            tolerance = 1e-4 if field.name == "nu_tilde_ratio" else 2e-6
            if field.name in ("velocity_ratio", "temperature_ratio"):
                tolerance = 2e-7
        np.testing.assert_allclose(
            actual, reference, rtol=0.0, atol=tolerance * np.max(np.abs(reference))
        )
    if case == "high_reynolds":
        assert old.nu_tilde_ratio is not None and new.nu_tilde_ratio is not None
        maximum = int(np.argmax(np.abs(new.nu_tilde_ratio[0] - old.nu_tilde_ratio[0])))
        # The guarded trajectory differs most two nodes before clipping (8.099).
        # The last positive node changes 1996.166 -> 2002.957 (6.791), and the
        # following Picard-bound node is effectively zero in both results.
        front = maximum + 1
        assert old.nu_tilde_ratio[0, front] > 0.0
        assert new.nu_tilde_ratio[0, front] > 0.0
        assert old.nu_tilde_ratio[0, front + 1] < 1e-12
        assert new.nu_tilde_ratio[0, front + 1] < 1e-12
        weights = instances[-1].wall_weights
        # Observed relative wall velocity slope difference: 1.537e-7.
        np.testing.assert_allclose(
            new.velocity_ratio[:, :3] @ weights,
            old.velocity_ratio[:, :3] @ weights,
            rtol=4e-7,
            atol=0.0,
        )
    if case in ("pressure", "adiabatic"):
        # Observed old/new counts: pressure 153/0, adiabatic 135/0.
        assert old_fallbacks >= 100
        assert new_fallbacks <= old_fallbacks // 10
    if case == "pressure":
        # Bypass the new override entirely to check the False switch bit for bit,
        # including the old positivity damping and Picard fallback trajectory.
        monkeypatch.setattr(TurbulentDiscretization, "newton", Discretization.newton)
        direct_legacy = run()
        for field in fields(old):
            np.testing.assert_array_equal(
                getattr(old, field.name), getattr(direct_legacy, field.name)
            )
        assert instances[-1].picard_fallbacks == old_fallbacks


@pytest.mark.solver_verification
def test_sa_switch_cannot_change_laminar_marches(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    edge = constant_edge(2.0, 220.0, [0.1, 1.0])
    monkeypatch.setattr(TurbulentDiscretization, "_bound_constrained_newton", False)
    old = marching_boundary_layer(
        edge, regime=BoundaryLayerRegime.LAMINAR, wall_temperature=300.0
    )
    monkeypatch.setattr(TurbulentDiscretization, "_bound_constrained_newton", True)
    new = marching_boundary_layer(
        edge, regime=BoundaryLayerRegime.LAMINAR, wall_temperature=300.0
    )
    for field in fields(old):
        np.testing.assert_array_equal(
            getattr(old, field.name), getattr(new, field.name)
        )


@pytest.fixture
def small_station() -> StationInput:
    """A six-node, positive-SA station; no march is needed for failure-path tests."""
    solver = TurbulentDiscretization(
        np.array([0.0, 0.1, 0.5, 1.0, 3.0, 8.0]), 0.72, 0.9, AIR_VISCOSITY
    )
    edge = StationEdge(0.01, 1.0, 0.0, 0.0, 220.0, 1.0, 1.0, 0.0, 0.1)
    guess = solver.starting_guess(edge)
    solver.enforce_boundaries(guess, edge)
    return solver, guess, edge, 100.0, -100.0 * guess


def test_small_station_and_legacy_switch(
    small_station: StationInput, monkeypatch: pytest.MonkeyPatch
) -> None:
    solver, guess, edge, coefficient, lagged = small_station
    active_set = solver.newton(guess, edge, coefficient, lagged)
    assert active_set is not None
    assert np.all(active_set.state[3, 1:-1] > 0.0)
    np.testing.assert_allclose(active_set.state[2], 1.0, rtol=0.0, atol=1e-14)
    monkeypatch.setattr(solver, "_bound_constrained_newton", False)
    legacy = solver.newton(guess, edge, coefficient, lagged)
    direct = Discretization.newton(solver, guess, edge, coefficient, lagged)
    assert legacy is not None and direct is not None
    np.testing.assert_array_equal(legacy.state, direct.state)
    np.testing.assert_allclose(active_set.state, legacy.state, rtol=1e-12, atol=1e-12)


def test_explicit_zero_seed_preserves_legacy_laminar_limit(
    small_station: StationInput,
) -> None:
    solver, guess, edge, coefficient, lagged = small_station
    edge = replace(edge, sa_ratio=0.0)
    guess[3] = 0.0
    lagged[3] = 0.0
    actual = solver.newton(guess, edge, coefficient, lagged)
    legacy = Discretization.newton(solver, guess, edge, coefficient, lagged)
    assert actual is not None and legacy is not None
    np.testing.assert_array_equal(actual.state, legacy.state)
    np.testing.assert_array_equal(actual.state[3], 0.0)


@pytest.mark.parametrize(("n", "expected"), [(0.0, 0.0), (2.0, 0.45)])
def test_legacy_positivity_damping(
    small_station: StationInput, n: float, expected: float
) -> None:
    solver, guess, _, _, _ = small_station
    guess[3, 2] = n
    update = np.zeros_like(guess)
    update[3, 2] = -4.0
    assert solver.damping(guess, update) == pytest.approx(expected)


@pytest.mark.parametrize("diagonal", [0.0, 1.0, np.nan, np.inf, -np.inf])
def test_invalid_sa_diagonal_requests_picard(
    small_station: StationInput, monkeypatch: pytest.MonkeyPatch, diagonal: float
) -> None:
    solver, guess, edge, coefficient, lagged = small_station
    original = solver.jacobian

    def jacobian(*args: object) -> FloatArray:
        bands = original(
            guess,
            solver.residual(guess, edge, coefficient, lagged),
            edge,
            coefficient,
            lagged,
        )
        bands[2 * solver.band, 8] = diagonal
        return bands

    monkeypatch.setattr(solver, "jacobian", jacobian)
    before = solver.picard_fallbacks
    assert solver.solve(guess, edge, coefficient, lagged) is not None
    assert solver.picard_fallbacks == before + 1


@pytest.mark.parametrize("failure", ["nonfinite", "value", "floating", "singular"])
def test_banded_failure_requests_fallback(
    small_station: StationInput, monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    solver, guess, edge, coefficient, lagged = small_station

    def fail(*args: object, **kwargs: object) -> FloatArray:
        if failure == "nonfinite":
            return np.full(guess.size, np.nan)
        error = {
            "value": ValueError,
            "floating": FloatingPointError,
            "singular": np.linalg.LinAlgError,
        }[failure]
        raise error("invalid trial")

    monkeypatch.setattr(sa_module, "solve_banded", fail)
    assert solver.newton(guess, edge, coefficient, lagged) is None


def test_iteration_exhaustion(
    small_station: StationInput, monkeypatch: pytest.MonkeyPatch
) -> None:
    solver, guess, edge, coefficient, lagged = small_station
    monkeypatch.setattr(solver, "newton_iterations", 0)
    assert solver.newton(guess, edge, coefficient, lagged) is None


@pytest.mark.parametrize("wall_attached", ["none", "active", "projected"])
def test_projection_and_stable_active_set_use_applied_update(
    small_station: StationInput,
    monkeypatch: pytest.MonkeyPatch,
    wall_attached: str,
) -> None:
    solver, guess, edge, coefficient, lagged = small_station
    guess[3, 1:-1] = 0.0
    if wall_attached == "none":
        guess[3, 1] = 1.0
    first_active = 1 if wall_attached == "active" else 2
    calls = 0

    def residual(*args: object) -> FloatArray:
        result = np.zeros_like(guess)
        result[3, first_active:-1] = -1.0
        if wall_attached == "projected":
            result[3, 1] = 1.0  # A positive residual selects the inactive branch.
        return result

    def jacobian(*args: object) -> FloatArray:
        bands = np.ones((3 * solver.band + 1, guess.size), order="F")
        bands[2 * solver.band, 8:-5:5] = -1.0
        return bands

    def update(
        lower: int, upper: int, bands: FloatArray, vector: FloatArray, **kwargs: object
    ) -> FloatArray:
        nonlocal calls
        calls += 1
        # Verify entire active rows, including the first/last interior nodes,
        # in the existing expanded gbsv layout. Inactive entries stay untouched.
        for row in range(5 * first_active + 3, guess.size - 5, 5):
            columns = np.arange(max(0, row - lower), min(guess.size, row + upper + 1))
            np.testing.assert_array_equal(
                bands[2 * lower + row - columns, columns], columns == row
            )
            assert vector[row] == 0.0
        assert bands[2 * lower, 2] == 1.0
        result = np.zeros_like(guess)
        result[3, first_active:-1] = -1.0
        if wall_attached == "projected":
            result[3, 1] = -1.0  # Projection can reach zero on an inactive row.
        return result.T.ravel()

    monkeypatch.setattr(solver, "residual", residual)
    monkeypatch.setattr(solver, "jacobian", jacobian)
    monkeypatch.setattr(sa_module, "solve_banded", update)
    station = solver.newton(guess, edge, coefficient, lagged)
    assert calls == 2  # A zero applied update alone cannot accept the first active set.
    if wall_attached != "none":
        assert station is None
        fallback = solver.station(guess, edge.xi)
        monkeypatch.setattr(solver, "picard", lambda *args: fallback)
        assert solver.solve(guess, edge, coefficient, lagged) is fallback
        assert solver.picard_fallbacks == 1
    else:
        assert station is not None
        np.testing.assert_array_equal(station.state, guess)
