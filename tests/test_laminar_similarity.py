"""Independent similarity numerics, reference provenance, and model accuracy."""

import csv
import runpy
from collections.abc import Callable
from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics.boundary_layer import (
    _BLASIUS_DELTA99,
    BoundaryLayerRegime,
    flat_plate_boundary_layer,
)
from aerophysics.transport import AIR_VISCOSITY

_SOLVER = runpy.run_path(
    str(
        Path(__file__).resolve().parents[1]
        / "docs/scripts/build_laminar_similarity_reference.py"
    )
)
solve_similarity = cast(Callable[..., Any], _SOLVER["solve_similarity"])
reference_solutions = cast(
    Callable[[], list[dict[str, float]]], _SOLVER["reference_solutions"]
)

REFERENCE = (
    Path(__file__).parent
    / "reference_data/viscous_flow/laminar_compressible_similarity.csv"
)


def reference_rows() -> list[dict[str, float]]:
    with REFERENCE.open(newline="", encoding="utf-8") as stream:
        return [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(stream)
        ]


def test_similarity_blasius_limit_and_integral_constants() -> None:
    low_mach = solve_similarity(0.01).quantities()
    assert_allclose(
        [
            low_mach["cf_sqrt_re"],
            low_mach["delta99_sqrt_re_over_x"],
            low_mach["displacement_sqrt_re_over_x"],
        ],
        [0.6641, 4.9100, 1.7208],
        rtol=0.0,
        atol=1e-4,
    )
    blasius = solve_similarity(0.0).quantities()
    assert_allclose(
        [blasius["eta99"], blasius["I1"], blasius["I2"]],
        [_BLASIUS_DELTA99, 3.19409, 2.53484],
        rtol=0.0,
        atol=1e-4,
    )


@pytest.mark.generated_assets
def test_similarity_reference_regenerates_and_obeys_momentum_integral() -> None:
    generated = reference_solutions()
    stored = reference_rows()
    assert len(generated) == len(stored) == 24
    assert {
        (row["mach"], row["wall_to_adiabatic_temperature_ratio"]) for row in stored
    } == {
        (mach, wall) for mach in (0.3, 2, 4, 6, 10, 15) for wall in (0.1, 0.3, 0.5, 1)
    }
    for actual, expected in zip(generated, stored, strict=True):
        for key in expected:
            assert_allclose(
                actual[key],
                expected[key],
                rtol=2e-8,
                atol=1e-10,
                err_msg=f"{REFERENCE}: M={expected['mach']}, {key}",
            )
        assert_allclose(
            actual["momentum_sqrt_re_over_x"], actual["cf_sqrt_re"], rtol=1e-8
        )


@pytest.mark.parametrize(
    "mach,expected",
    [
        (2.0, [0.6356, 6.415, 3.303, 1.676]),
        (10.0, [0.4258, 29.67, 27.67, 17.63]),
    ],
)
def test_similarity_matches_supplied_independent_examples(
    mach: float,
    expected: list[float],
) -> None:
    row = next(
        row
        for row in reference_rows()
        if row["mach"] == mach and row["wall_to_adiabatic_temperature_ratio"] == 1
    )
    assert_allclose(
        [
            row["cf_sqrt_re"],
            row["delta99_sqrt_re_over_x"],
            row["displacement_sqrt_re_over_x"],
            row["adiabatic_temperature_K"] / 240.0,
        ],
        expected,
        rtol=0.005,
    )


@pytest.mark.parametrize("mach,wall_ratio", [(0.3, 0.1), (10.0, 1.0), (15.0, 0.1)])
def test_similarity_domain_and_tolerance_convergence(
    mach: float, wall_ratio: float
) -> None:
    solution = solve_similarity(0.0)
    for value in np.linspace(0, mach, max(2, int(mach * 2) + 1))[1:]:
        solution = solve_similarity(float(value), previous=solution)
    temperature_ratio = solution.quantities()["wall_temperature_K"] / 240.0
    for value in np.linspace(1.0, wall_ratio, 10)[1:]:
        solution = solve_similarity(
            mach, float(value) * temperature_ratio, previous=solution
        )
    refined = solve_similarity(
        mach,
        None if wall_ratio == 1.0 else wall_ratio * temperature_ratio,
        previous=solution,
        extent=20.0,
        tolerance=1e-10,
    )
    first, second = solution.quantities(), refined.quantities()
    for key in (
        "cf_sqrt_re",
        "delta99_sqrt_re_over_x",
        "displacement_sqrt_re_over_x",
        "momentum_sqrt_re_over_x",
        "wall_temperature_K",
    ):
        assert_allclose(first[key], second[key], rtol=2e-7, atol=1e-9)
    assert np.max(refined.bvp.rms_residuals) <= 1e-10


def test_compressible_laminar_accuracy_against_similarity() -> None:
    rows = reference_rows()
    viscosity = float(AIR_VISCOSITY.dynamic_viscosity(240.0))
    for row in rows:
        ratio = row["wall_to_adiabatic_temperature_ratio"]
        # Re=1e8 keeps chi below one over this mathematical comparison.
        # Cooled walls use the reference's solved Taw to set Tw.
        result = flat_plate_boundary_layer(
            1.0,
            1e8 * viscosity,
            1.0,
            viscosity,
            regime=BoundaryLayerRegime.LAMINAR,
            compressible=True,
            mach=row["mach"],
            edge_temperature=240.0,
            wall_temperature=None if ratio == 1 else row["wall_temperature_K"],
        )
        for value, key, tolerance in (
            (result.local_skin_friction_coefficient, "cf_sqrt_re", 1e-3),
            (result.momentum_thickness, "momentum_sqrt_re_over_x", 1e-3),
            (result.boundary_layer_thickness, "delta99_sqrt_re_over_x", 1e-3),
        ):
            assert_allclose(
                float(value) * 1e4,
                row[key],
                rtol=tolerance,
                err_msg=f"M={row['mach']}, Tw/Taw={ratio}, {key}",
            )
        displacement = row["displacement_sqrt_re_over_x"]
        assert_allclose(
            float(result.displacement_thickness) * 1e4,
            displacement,
            rtol=1e-3 if abs(displacement) >= 0.05 else 0,
            atol=0 if abs(displacement) >= 0.05 else 5e-5,
        )
        assert result.recovery_temperature is not None
        assert_allclose(
            result.recovery_temperature, row["adiabatic_temperature_K"], rtol=1e-4
        )
        assert_allclose(
            result.average_skin_friction_coefficient,
            2 * np.asarray(result.local_skin_friction_coefficient),
            rtol=1e-12,
        )
        assert_allclose(
            result.momentum_thickness,
            result.local_skin_friction_coefficient,
            rtol=1e-12,
        )


def test_true_adiabatic_wall_and_prescribed_wall_recovery() -> None:
    mu = float(AIR_VISCOSITY.dynamic_viscosity(240.0))
    kwargs: dict[str, Any] = dict(
        regime=BoundaryLayerRegime.LAMINAR,
        compressible=True,
        mach=10.0,
        edge_temperature=240.0,
    )
    adiabatic = flat_plate_boundary_layer(1.0, 1e8 * mu, 1.0, mu, **kwargs)
    cooled = flat_plate_boundary_layer(
        1.0, 1e8 * mu, 1.0, mu, wall_temperature=300.0, **kwargs
    )
    reference = solve_similarity(10.0)
    assert abs(reference.bvp.y[4, 0]) < 1e-10
    assert adiabatic.recovery_temperature is not None
    taw = float(adiabatic.recovery_temperature)
    approximate = 240 * (1 + np.sqrt(0.72) * 0.2 * 100)
    assert_allclose(taw, reference.quantities()["wall_temperature_K"], rtol=1e-5)
    assert cooled.recovery_temperature is not None
    assert_allclose(cooled.recovery_temperature, taw, rtol=1e-6)
    assert adiabatic.wall_temperature == adiabatic.recovery_temperature
    assert cooled.wall_temperature == 300.0
    assert -0.02 < taw / approximate - 1 < -0.015


@pytest.mark.parametrize("wall", [None, 300.0])
def test_distance_sweep_solves_one_exact_state_per_wall_condition(
    monkeypatch: pytest.MonkeyPatch, wall: float | None
) -> None:
    from aerophysics import _laminar_similarity as solver

    original = solver._richardson
    batches: list[int] = []

    def recording(states: np.ndarray, *args: Any, **kwargs: Any) -> Any:
        batches.append(len(states))
        return original(states, *args, **kwargs)

    monkeypatch.setattr(solver, "_richardson", recording)
    mu = float(AIR_VISCOSITY.dynamic_viscosity(240.0))
    flat_plate_boundary_layer(
        np.linspace(0.1, 1, 10000),
        1e8 * mu,
        1.0,
        mu,
        regime=BoundaryLayerRegime.LAMINAR,
        compressible=True,
        mach=2.0,
        edge_temperature=240.0,
        wall_temperature=wall,
    )
    assert batches == [1 if wall is None else 2]


def test_no_rounding_of_physical_states(monkeypatch: pytest.MonkeyPatch) -> None:
    from aerophysics import _laminar_similarity as solver

    original = solver._richardson
    captured: list[np.ndarray] = []

    def recording(states: np.ndarray, *args: Any, **kwargs: Any) -> Any:
        captured.append(states.copy())
        return original(states, *args, **kwargs)

    monkeypatch.setattr(solver, "_richardson", recording)
    solver.similarity_coefficients(
        np.array(
            [
                [2.0, 240.0, np.nan, 0.72, 1.4],
                [np.nextafter(2.0, 3.0), 240.0, np.nan, 0.72, 1.4],
            ]
        ),
        AIR_VISCOSITY,
    )
    assert len(captured[0]) == 2
    assert captured[0][0, 0] != captured[0][1, 0]


class _PowerViscosity:
    """Array-capable, deterministic mu proportional to T: C is exactly one."""

    def dynamic_viscosity(self, temperature: Any) -> Any:
        return 1e-5 * np.asarray(temperature) / 240


class _DifferentiablePowerViscosity(_PowerViscosity):
    def dynamic_viscosity_derivative(self, temperature: Any) -> Any:
        return np.full_like(temperature, 1e-5 / 240)


@pytest.mark.parametrize("model", [_PowerViscosity(), _DifferentiablePowerViscosity()])
def test_generic_viscosity_and_edge_viscosity_scaling(model: Any) -> None:
    from aerophysics._laminar_similarity import similarity_coefficients

    q = similarity_coefficients(np.array([[6.0, 240.0, np.nan, 1.0, 1.4]]), model)[0]
    assert_allclose(q[0], 0.6641146724, rtol=2e-6)
    assert_allclose(q[4], 1 + 0.2 * 36, rtol=1e-5)
    results = [
        flat_plate_boundary_layer(
            1.0,
            1e8 * mu,
            1.0,
            mu,
            regime=BoundaryLayerRegime.LAMINAR,
            compressible=True,
            mach=6.0,
            edge_temperature=240.0,
            prandtl_number=1.0,
            viscosity_model=model,
        )
        for mu in (1e-5, 4e-5)
    ]
    for field in (
        "local_skin_friction_coefficient",
        "boundary_layer_thickness",
        "displacement_thickness",
        "momentum_thickness",
    ):
        assert_allclose(
            getattr(results[1], field),
            0.5 * np.asarray(getattr(results[0], field)),
            rtol=1e-12,
        )
    assert results[0].recovery_temperature == results[1].recovery_temperature


@pytest.mark.parametrize("wall", [np.nan, 1.25])
def test_collocation_fallback_matches_independent_reference(wall: float) -> None:
    from aerophysics._laminar_similarity import similarity_coefficients

    q = similarity_coefficients(
        np.array([[6.0, 240.0, wall, 0.72, 1.4]]), AIR_VISCOSITY, _max_iterations=0
    )[0]
    reference = solve_similarity(6.0, None if np.isnan(wall) else wall).quantities()
    assert_allclose(
        q[:4],
        [
            reference[key]
            for key in (
                "cf_sqrt_re",
                "delta99_sqrt_re_over_x",
                "displacement_sqrt_re_over_x",
                "momentum_sqrt_re_over_x",
            )
        ],
        rtol=1e-6,
    )
    assert_allclose(q[4], reference["wall_temperature_K"] / 240, rtol=1e-6)


def test_failed_iteration_reports_original_indices_and_state() -> None:
    from aerophysics._laminar_similarity import similarity_coefficients
    from aerophysics.exceptions import ModelRangeError

    states = np.array(
        [
            [10.0, 240.0, np.nan, 0.72, 1.4],
            [2.0, 240.0, 1.25, 0.72, 1.4],
            [2.0, 240.0, 1.25, 0.72, 1.4],
        ]
    )
    with pytest.raises(
        ModelRangeError, match=r"array indices \[1, 2\].*2.0, 240.0, 1.25"
    ):
        similarity_coefficients(
            states, AIR_VISCOSITY, _max_iterations=0, _allow_fallback=False
        )


@pytest.mark.parametrize(
    "column,value", [(0, -1), (1, 0), (2, 0), (2, np.inf), (3, 0), (4, 1)]
)
def test_invalid_state_has_index(column: int, value: float) -> None:
    from aerophysics._laminar_similarity import similarity_coefficients
    from aerophysics.exceptions import ModelRangeError

    states = np.array([[2.0, 240.0, 1.0, 0.72, 1.4]] * 2)
    states[1, column] = value
    with pytest.raises(ModelRangeError, match="array index 1"):
        similarity_coefficients(states, AIR_VISCOSITY)


def test_nonpositive_viscosity_fails_with_state() -> None:
    from aerophysics._laminar_similarity import similarity_coefficients
    from aerophysics.exceptions import ModelRangeError

    class InvalidViscosity(_PowerViscosity):
        def dynamic_viscosity(self, temperature: Any) -> Any:
            return np.where(
                np.asarray(temperature) >= 500,
                -1.0,
                super().dynamic_viscosity(temperature),
            )

    with pytest.raises(
        ModelRangeError, match=r"array indices \[1\].*nonpositive.*viscosity"
    ):
        similarity_coefficients(
            np.array([[0.0, 240.0, 1.0, 0.72, 1.4], [0.0, 600.0, 1.0, 0.72, 1.4]]),
            InvalidViscosity(),
        )


@pytest.mark.generated_assets
def test_random_similarity_states_are_robust_and_match_collocation() -> None:
    from aerophysics._laminar_similarity import _collocation, similarity_coefficients

    rng = np.random.default_rng(20261102)
    mach = rng.uniform(0, 20, 200)
    temperature = np.exp(rng.uniform(np.log(30), np.log(2000), 200))
    ratio = np.exp(rng.uniform(np.log(0.05), np.log(1.5), 200))
    states = np.column_stack(
        (mach, temperature, np.full(200, np.nan), np.full(200, 0.72), np.full(200, 1.4))
    )
    adiabatic = similarity_coefficients(states, AIR_VISCOSITY)
    states[:, 2] = adiabatic[:, 4] * ratio
    actual = similarity_coefficients(states, AIR_VISCOSITY)
    assert np.all(np.isfinite(actual))
    assert_allclose(actual[:, 0], actual[:, 3], rtol=1e-12)
    for index in range(0, 200, 10):
        reference = _collocation(states[index], AIR_VISCOSITY)
        assert_allclose(actual[index, [0, 1, 3, 4]], reference[[0, 1, 3, 4]], rtol=1e-3)
        assert abs(actual[index, 2] - reference[2]) <= max(
            5e-5, 1e-3 * abs(reference[2])
        )


def test_public_failure_preserves_broadcast_index(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from aerophysics import boundary_layer
    from aerophysics._laminar_similarity import similarity_coefficients
    from aerophysics.exceptions import ModelRangeError

    def fail(states: np.ndarray, model: Any, **kwargs: Any) -> np.ndarray:
        return similarity_coefficients(
            states, model, _max_iterations=0, _allow_fallback=False, **kwargs
        )

    monkeypatch.setattr(boundary_layer, "similarity_coefficients", fail)
    with pytest.raises(
        ModelRangeError, match=r"array indices \[\(0, 1\)\].*2.0, 240.0"
    ):
        flat_plate_boundary_layer(
            np.array([[1.0, 2.0]]),
            1000.0,
            0.1,
            1.5e-5,
            regime=BoundaryLayerRegime.LAMINAR,
            compressible=True,
            mach=np.array([[10.0, 2.0]]),
            edge_temperature=240.0,
        )


def test_far_field_gate_routes_to_larger_collocation_domain(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from aerophysics import _laminar_similarity as solver

    state = np.array([[2.0, 240.0, np.nan, 0.72, 1.4]])
    # Leave the cached initial profile unchanged while shortening only FD grids.
    solver._blasius_guess()
    monkeypatch.setattr(solver, "_EXTENT", 4.0)
    grid = solver._finite_difference(state, 321, AIR_VISCOSITY)
    assert not grid.valid[0]
    actual = solver.similarity_coefficients(state, AIR_VISCOSITY)[0]
    reference = solve_similarity(2.0).quantities()
    assert_allclose(actual[0], reference["cf_sqrt_re"], rtol=1e-6)
    assert_allclose(actual[4], reference["wall_temperature_K"] / 240, rtol=1e-6)


def test_viscosity_domain_exception_is_indexed() -> None:
    from aerophysics._laminar_similarity import similarity_coefficients
    from aerophysics.exceptions import ModelRangeError

    class RestrictedViscosity(_PowerViscosity):
        def dynamic_viscosity(self, temperature: Any) -> Any:
            if np.any(np.asarray(temperature) > 500):
                raise ValueError("outside model temperature domain")
            return super().dynamic_viscosity(temperature)

    with pytest.raises(
        ModelRangeError, match=r"array indices \[0\].*outside model temperature domain"
    ):
        similarity_coefficients(
            np.array([[10.0, 240.0, np.nan, 0.72, 1.4]]), RestrictedViscosity()
        )


def test_state_dependent_prandtl_gamma_and_far_field_extension() -> None:
    from aerophysics._laminar_similarity import _collocation, similarity_coefficients

    states = np.array([[2.0, 300.0, np.nan, 0.1, 1.667], [4.0, 100.0, 2.0, 5.0, 1.2]])
    actual = similarity_coefficients(states, AIR_VISCOSITY)
    for state, coefficients in zip(states, actual, strict=True):
        reference = _collocation(state, AIR_VISCOSITY)
        assert_allclose(coefficients, reference, rtol=1e-4)
