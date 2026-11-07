"""Optional native station solver dispatch and numerical equivalence."""

import sys
from dataclasses import replace
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest
from test_marching_boundary_layer import smooth_edge

from aerophysics import BoundaryLayerEdge, BoundaryLayerRegime, marching_boundary_layer
from aerophysics._boundary_layer_march import Discretization, StationEdge, eta_grid
from aerophysics._boundary_layer_sa import TurbulentDiscretization
from aerophysics._native_march import create_adapter, load_extension, viscosity_spec
from aerophysics.exceptions import ApplicabilityWarning
from aerophysics.marching_boundary_layer import _EdgeInterpolant
from aerophysics.transport import (
    AIR_KEYES_VISCOSITY,
    AIR_VISCOSITY,
    BlottnerModel,
    SutherlandModel,
)


class FakeKernel:
    result: Any = None
    out_of_range = False

    def __init__(
        self,
        eta: Any,
        prandtl_number: float,
        kind: int,
        parameters: tuple[float, ...],
        jacobian: str = "analytic",
    ) -> None:
        self.jacobian = jacobian
        self.kind = kind
        self.parameters = parameters

    def newton(self, *args: Any) -> tuple[Any, bool]:
        return self.result, self.out_of_range


@pytest.fixture
def fake_extension(monkeypatch: pytest.MonkeyPatch) -> None:
    FakeKernel.result = None
    FakeKernel.out_of_range = False
    monkeypatch.setitem(
        sys.modules, "aerophysics_native", SimpleNamespace(LaminarKernel=FakeKernel)
    )


@pytest.mark.parametrize("setting", ["", "auto", "1"])
def test_enabled_modes(
    monkeypatch: pytest.MonkeyPatch, fake_extension: None, setting: str
) -> None:
    monkeypatch.setenv("AEROPHYSICS_NATIVE", setting)
    assert load_extension() is sys.modules["aerophysics_native"]


def test_missing_and_invalid_modes(monkeypatch: pytest.MonkeyPatch) -> None:
    def missing(name: str) -> Any:
        raise ImportError(name)

    monkeypatch.setattr("aerophysics._native_march.importlib.import_module", missing)
    monkeypatch.delenv("AEROPHYSICS_NATIVE", raising=False)
    assert load_extension() is None
    monkeypatch.setenv("AEROPHYSICS_NATIVE", "0")
    assert load_extension() is None
    monkeypatch.setenv("AEROPHYSICS_NATIVE", "1")
    with pytest.raises(ImportError, match="requires aerophysics_native"):
        load_extension()
    monkeypatch.setenv("AEROPHYSICS_NATIVE", "other")
    with pytest.raises(ValueError, match="AEROPHYSICS_NATIVE"):
        load_extension()


def test_selection_and_dispatch(
    monkeypatch: pytest.MonkeyPatch, fake_extension: None
) -> None:
    monkeypatch.setenv("AEROPHYSICS_NATIVE", "auto")
    eta = np.array([0.0, 0.5, 1.0])
    edge = StationEdge(0.0, 1.0, 0.0, 0.0, 300.0, 1.0)
    guess = np.ones((3, 3))
    lagged = np.zeros_like(guess)
    assert viscosity_spec(AIR_VISCOSITY) == (
        0,
        (
            AIR_VISCOSITY.reference_viscosity,
            AIR_VISCOSITY.reference_temperature,
            AIR_VISCOSITY.sutherland_temperature,
        ),
    )
    assert viscosity_spec(AIR_KEYES_VISCOSITY) == (
        1,
        (
            AIR_KEYES_VISCOSITY.coefficient,
            AIR_KEYES_VISCOSITY.additive_temperature,
            AIR_KEYES_VISCOSITY.exponential_temperature,
            *AIR_KEYES_VISCOSITY.temperature_range,
        ),
    )
    unsupported = BlottnerModel(0.026, 0.3, -11.3)
    assert viscosity_spec(unsupported) is None

    class DerivedSutherland(SutherlandModel):
        pass

    assert viscosity_spec(DerivedSutherland(1e-5, 300.0, 100.0)) is None
    assert Discretization(eta, 0.72, unsupported)._native is None
    assert TurbulentDiscretization(eta, 0.72, 0.9, AIR_VISCOSITY)._native is None
    disc = Discretization(eta, 0.72, AIR_VISCOSITY)
    assert disc.newton(guess, edge, 0.0, lagged) is None
    FakeKernel.result = guess
    station = disc.newton(guess, edge, 0.0, lagged)
    assert station is not None
    np.testing.assert_array_equal(station.state, guess)
    FakeKernel.out_of_range = True
    assert disc.newton(guess, edge, 0.0, lagged) is not None
    keyes = Discretization(eta, 0.72, AIR_KEYES_VISCOSITY)
    with pytest.warns(
        ApplicabilityWarning, match="outside the fitted range 79--1845 K"
    ):
        assert keyes.newton(guess, edge, 0.0, lagged) is not None
    monkeypatch.setenv("AEROPHYSICS_NATIVE", "0")
    assert Discretization(eta, 0.72, AIR_VISCOSITY)._native is None


@pytest.mark.parametrize("jacobian", ["analytic", "finite-difference"])
@pytest.mark.parametrize("case", ["adiabatic", "isothermal", "keyes"])
def test_real_native_agrees_with_python(
    monkeypatch: pytest.MonkeyPatch, case: str, jacobian: str
) -> None:
    pytest.importorskip("aerophysics_native")
    monkeypatch.setenv("AEROPHYSICS_NATIVE_JACOBIAN", jacobian)
    rtol, atol = (1e-8, 1e-10) if jacobian == "analytic" else (1e-10, 1e-12)
    edge = smooth_edge(2.0, 3.0, count=51)
    if case == "adiabatic":
        edge = BoundaryLayerEdge.from_pressure(
            np.linspace(0.0, 1.0, 51),
            np.full(51, 2e4),
            total_pressure=1.6e5,
            total_temperature=500.0,
        )
    if case == "keyes":
        edge = replace(edge, viscosity_model=AIR_KEYES_VISCOSITY)
    wall = None if case == "adiabatic" else 300.0
    wall_array = None if wall is None else np.full_like(edge.distance, wall)
    station_edge = _EdgeInterpolant(edge, wall_array).state(0.0)
    results = []
    starts = []
    for mode in ("0", "1"):
        monkeypatch.setenv("AEROPHYSICS_NATIVE", mode)
        disc = Discretization(
            eta_grid(1e-3, 1.015, 0.05, 16.0), 0.72, edge.viscosity_model
        )
        start = disc.newton(
            disc.starting_guess(station_edge),
            station_edge,
            0.0,
            np.zeros((3, len(disc.eta))),
        )
        assert start is not None
        starts.append(start.state)
        results.append(
            marching_boundary_layer(
                edge, regime=BoundaryLayerRegime.LAMINAR, wall_temperature=wall
            )
        )
    np.testing.assert_allclose(starts[0], starts[1], rtol=rtol, atol=atol)
    for name in (
        "local_skin_friction_coefficient",
        "wall_temperature" if wall is None else "wall_heat_flux",
        "velocity_ratio",
        "temperature_ratio",
    ):
        np.testing.assert_allclose(
            getattr(results[0], name), getattr(results[1], name), rtol=rtol, atol=atol
        )


@pytest.mark.parametrize("jacobian", ["analytic", "finite-difference"])
def test_real_failure_shapes_and_noncontiguous(
    monkeypatch: pytest.MonkeyPatch, jacobian: str
) -> None:
    monkeypatch.setenv("AEROPHYSICS_NATIVE_JACOBIAN", jacobian)
    native = pytest.importorskip("aerophysics_native")
    eta = np.linspace(0.0, 10.0, 101)
    for invalid in (eta[:2], np.array([0.0, 1.0, 1.0])):
        with pytest.raises(ValueError):
            native.LaminarKernel(invalid, 0.72, 0, (1e-5, 300.0, 100.0))
    for kind, parameters in ((2, ()), (-1, ()), (0, ())):
        with pytest.raises(ValueError):
            native.LaminarKernel(eta, 0.72, kind, parameters)
    with pytest.raises(ValueError, match="jacobian"):
        native.LaminarKernel(eta, 0.72, 0, (1e-5, 300.0, 100.0), "invalid")
    kernel = native.LaminarKernel(eta, 0.72, 0, (1e-5, 300.0, 100.0), jacobian)
    with pytest.raises(ValueError):
        kernel.newton(
            0.0,
            0.0,
            0.0,
            300.0,
            None,
            0.0,
            np.zeros((3, 3)),
            np.zeros((3, 3)),
            20,
            1e-10,
        )
    monkeypatch.setenv("AEROPHYSICS_NATIVE", "0")
    disc = Discretization(eta, 0.72, AIR_VISCOSITY)
    edge = StationEdge(0.0, 1.0, 0.0, 0.0, 300.0, None)
    guess = disc.starting_guess(edge)
    guess[2, 0] = 0.0
    lagged = np.zeros_like(guess)
    assert disc.newton(guess, edge, 0.0, lagged) is None
    monkeypatch.setenv("AEROPHYSICS_NATIVE", "1")
    disc = Discretization(eta, 0.72, AIR_VISCOSITY)
    assert disc.newton(guess, edge, 0.0, lagged) is None
    guess = disc.starting_guess(edge)
    noncontiguous_guess = np.empty((3, 202))[:, ::2]
    noncontiguous_lagged = np.empty((3, 202))[:, ::2]
    noncontiguous_guess[:] = guess
    noncontiguous_lagged[:] = lagged
    assert not noncontiguous_guess.flags.c_contiguous
    args = (0.0, 0.0, 0.0, 300.0, None, 0.0)
    empty_state, empty_out = kernel.newton(*args, lagged, guess, -1, 1e-10)
    assert empty_state is None and not empty_out
    state, out = kernel.newton(*args, lagged, guess, 20, 1e-10)
    strided_state, strided_out = kernel.newton(
        *args, noncontiguous_lagged, noncontiguous_guess, 20, 1e-10
    )
    assert state is not None and strided_state is not None
    assert not out and not strided_out
    np.testing.assert_array_equal(state, strided_state)


@pytest.mark.parametrize("jacobian", ["analytic", "finite-difference"])
def test_real_keyes_range_and_perturbed_temperature(jacobian: str) -> None:
    native = pytest.importorskip("aerophysics_native")
    eta = np.linspace(0.0, 1.0, 5)
    parameters = (1.488e-6, 122.1, 5.0, 79.0, 1845.0)
    kernel = native.LaminarKernel(eta, 0.72, 1, parameters, jacobian)
    state = np.ones((3, 5))
    lagged = np.zeros_like(state)
    # The unperturbed temperature is at the upper limit; the Jacobian's
    # positive temperature perturbations cross it only in finite-difference mode.
    result, out = kernel.newton(
        0.0, 0.0, 0.0, 1845.0, 1.0, 0.0, lagged, state, 1, 1e-10
    )
    assert out == (jacobian == "finite-difference")
    assert result is None
    state[2, 0] = 0.0
    result, out = kernel.newton(
        0.0, 0.0, 0.0, 1900.0, 1.0, 0.0, lagged, state, 1, 1e-10
    )
    assert result is None and out


@pytest.mark.parametrize("jacobian", ["analytic", "finite-difference"])
@pytest.mark.parametrize("temperature", [np.nan, np.inf, -np.inf, 0.0, -1.0])
def test_real_nonfinite_temperatures_fail_both_paths(
    monkeypatch: pytest.MonkeyPatch, temperature: float, jacobian: str
) -> None:
    pytest.importorskip("aerophysics_native")
    monkeypatch.setenv("AEROPHYSICS_NATIVE_JACOBIAN", jacobian)
    eta = np.linspace(0.0, 10.0, 101)
    edge = StationEdge(0.0, 1.0, 0.0, 0.0, 300.0, None)
    for mode in ("0", "1"):
        monkeypatch.setenv("AEROPHYSICS_NATIVE", mode)
        disc = Discretization(eta, 0.72, AIR_VISCOSITY)
        guess = disc.starting_guess(edge)
        guess[2, 10] = temperature
        lagged = np.zeros_like(guess)
        assert disc.newton(guess, edge, 0.0, lagged) is None
        edge_nonfinite = replace(edge, temperature=temperature)
        assert (
            disc.newton(disc.starting_guess(edge), edge_nonfinite, 0.0, lagged) is None
        )


@pytest.mark.parametrize("setting", [None, "", "analytic", "finite-difference"])
def test_jacobian_selection(
    monkeypatch: pytest.MonkeyPatch, fake_extension: None, setting: str | None
) -> None:
    monkeypatch.setenv("AEROPHYSICS_NATIVE", "1")
    if setting is None:
        monkeypatch.delenv("AEROPHYSICS_NATIVE_JACOBIAN", raising=False)
    else:
        monkeypatch.setenv("AEROPHYSICS_NATIVE_JACOBIAN", setting)
    adapter = create_adapter(np.array([0.0, 0.5, 1.0]), 0.72, AIR_VISCOSITY)
    assert adapter is not None
    assert isinstance(adapter.kernel, FakeKernel)
    assert adapter.kernel.jacobian == (setting or "analytic")


def test_invalid_jacobian_selection(
    monkeypatch: pytest.MonkeyPatch, fake_extension: None
) -> None:
    monkeypatch.setenv("AEROPHYSICS_NATIVE", "1")
    monkeypatch.setenv("AEROPHYSICS_NATIVE_JACOBIAN", "other")
    with pytest.raises(ValueError, match="AEROPHYSICS_NATIVE_JACOBIAN"):
        create_adapter(np.array([0.0, 0.5, 1.0]), 0.72, AIR_VISCOSITY)


class FakeEdge:
    def __init__(self, *args: Any) -> None:
        self.args = args


class FakeMarch:
    status = "ok"
    out_of_range = False
    arguments: tuple[Any, ...] = ()

    @classmethod
    def run(cls, *args: Any) -> Any:
        cls.arguments = args
        n = len(args[3])
        return {
            "status": cls.status,
            "candidate": 0.123456789,
            "positions": np.array([0.1]),
            "states": np.ones((1, 3, n)),
            "slopes": np.array([[0.4, 0.2]]),
            "edges": np.ones((1, 7)),
            "accepted_positions": [0.0, 0.1, 0.2],
            "accepted_shear": [0.4, 0.3, 0.2],
            "out_of_range": cls.out_of_range,
        }


@pytest.fixture
def fake_march_extension(monkeypatch: pytest.MonkeyPatch, fake_extension: None) -> None:
    extension = sys.modules["aerophysics_native"]
    monkeypatch.setattr(extension, "NativeEdge", FakeEdge, raising=False)
    monkeypatch.setattr(extension, "march", FakeMarch.run, raising=False)
    FakeMarch.status = "ok"
    FakeMarch.out_of_range = False
    monkeypatch.setenv("AEROPHYSICS_NATIVE", "1")
    monkeypatch.delenv("AEROPHYSICS_NATIVE_MARCH", raising=False)


@pytest.mark.parametrize("setting", [None, "", "1", "0", "invalid"])
def test_full_march_dispatch(
    monkeypatch: pytest.MonkeyPatch, fake_march_extension: None, setting: str | None
) -> None:
    from aerophysics._native_march import native_march

    edge = smooth_edge(2.0, 3.0, count=5)
    disc = Discretization(eta_grid(1e-3, 1.015, 0.05, 16.0), 0.72, AIR_VISCOSITY)
    interpolant = _EdgeInterpolant(edge, None)
    if setting is not None:
        monkeypatch.setenv("AEROPHYSICS_NATIVE_MARCH", setting)
    if setting == "invalid":
        with pytest.raises(ValueError, match="AEROPHYSICS_NATIVE_MARCH"):
            native_march(interpolant, edge.distance[1:], disc, AIR_VISCOSITY)
    else:
        result = native_march(interpolant, edge.distance[1:], disc, AIR_VISCOSITY)
        if setting == "0":
            assert result is None
        else:
            assert result is not None
            solution, values = result
            np.testing.assert_array_equal(solution.velocity_ratio, 1.0)
            np.testing.assert_array_equal(solution.wall_velocity_slope, [0.4])
            assert values.shape == (1, 7)
            assert FakeMarch.arguments[-2:] == (20, 1e-10)


def test_full_march_fallbacks(
    monkeypatch: pytest.MonkeyPatch, fake_march_extension: None
) -> None:
    from aerophysics._native_march import native_march

    edge = smooth_edge(2.0, 3.0, count=5)
    interpolant = _EdgeInterpolant(edge, None)
    eta = np.linspace(0.0, 16.0, 101)
    unsupported = BlottnerModel(0.026, 0.3, -11.3)
    for disc, model in (
        (TurbulentDiscretization(eta, 0.72, 0.9, AIR_VISCOSITY), AIR_VISCOSITY),
        (Discretization(eta, 0.72, unsupported), unsupported),
        (Discretization(eta, 0.72, AIR_VISCOSITY), unsupported),
    ):
        assert native_march(interpolant, edge.distance[1:], disc, model) is None
    disc = Discretization(eta, 0.72, AIR_VISCOSITY)
    monkeypatch.setenv("AEROPHYSICS_NATIVE", "0")
    assert native_march(interpolant, edge.distance[1:], disc, AIR_VISCOSITY) is None
    disc = Discretization(eta, 0.72, AIR_VISCOSITY)
    assert native_march(interpolant, edge.distance[1:], disc, AIR_VISCOSITY) is None


@pytest.mark.parametrize(
    ("status", "message"),
    [
        ("invalid_start", "the march must start at xi = 0"),
        ("leading_edge_failed", "the leading-edge similarity profile did not converge"),
        ("not_converged", "the march did not converge at distance 0.123457"),
        (
            "outer_boundary",
            "the outer boundary is too close to the boundary layer at "
            "distance 0.123457",
        ),
    ],
)
def test_full_march_failure_messages(
    fake_march_extension: None, status: str, message: str
) -> None:
    from aerophysics._native_march import native_march
    from aerophysics.exceptions import BoundaryLayerConvergenceError

    FakeMarch.status = status
    edge = smooth_edge(2.0, 3.0, count=5)
    disc = Discretization(np.linspace(0.0, 16.0, 101), 0.72, AIR_VISCOSITY)
    with pytest.raises(BoundaryLayerConvergenceError) as caught:
        native_march(
            _EdgeInterpolant(edge, None), edge.distance[1:], disc, AIR_VISCOSITY
        )
    assert str(caught.value) == message


def test_full_march_separation_and_warning(fake_march_extension: None) -> None:
    from aerophysics._boundary_layer_march import SeparationError, _separation_position
    from aerophysics._native_march import native_march

    FakeMarch.status = "separated"
    FakeMarch.out_of_range = True
    edge = replace(smooth_edge(2.0, 3.0, count=5), viscosity_model=AIR_KEYES_VISCOSITY)
    disc = Discretization(np.linspace(0.0, 16.0, 101), 0.72, AIR_KEYES_VISCOSITY)
    wall = np.full_like(edge.distance, 300.0)
    with (
        pytest.warns(
            ApplicabilityWarning, match="outside the fitted range"
        ) as recorded,
        pytest.raises(SeparationError) as caught,
    ):
        native_march(
            _EdgeInterpolant(edge, wall), edge.distance[1:], disc, AIR_KEYES_VISCOSITY
        )
    assert len(recorded) == 1
    assert caught.value.distance == _separation_position(
        [0.0, 0.1, 0.2], [0.4, 0.3, 0.2]
    )
    assert caught.value.solution.velocity_ratio.shape == (1, 101)
    FakeMarch.status = "ok"
    # A spurious flag for Sutherland must not generate a Keyes warning.
    disc = Discretization(disc.eta, 0.72, AIR_VISCOSITY)
    assert native_march(
        _EdgeInterpolant(edge, wall), edge.distance[1:], disc, AIR_VISCOSITY
    )


@pytest.mark.parametrize("model", [AIR_VISCOSITY, AIR_KEYES_VISCOSITY])
@pytest.mark.parametrize("isothermal", [False, True])
def test_real_edge_and_ppoly(model: Any, isothermal: bool) -> None:
    from aerophysics._native_march import create_edge

    native = pytest.importorskip("aerophysics_native")
    edge = replace(smooth_edge(2.0, 3.0, count=51), viscosity_model=model)
    wall = 300.0 + 20.0 * np.sin(edge.distance) if isothermal else None
    interpolant = _EdgeInterpolant(edge, wall)
    spec = viscosity_spec(model)
    assert spec is not None
    rust = create_edge(native, interpolant, spec)
    rng = np.random.default_rng(1729)
    # Every knot, both extrapolation intervals, and many random interior points.
    samples = np.concatenate((edge.distance, rng.uniform(-0.01, 2.01, 2000)))
    actual_poly = np.array(rust.polynomials(samples))
    for col, poly in enumerate(
        (
            interpolant.log_pressure,
            interpolant.log_pressure_slope,
            interpolant.wall_temperature,
        )
    ):
        if poly is not None:
            actual = np.asarray(actual_poly[:, col], dtype=np.float64)
            np.testing.assert_array_equal(
                actual.view(np.uint64), poly(samples).view(np.uint64)
            )
    samples = np.concatenate(
        (
            edge.distance,
            (edge.distance[:-1] + edge.distance[1:]) / 2.0,
            rng.uniform(0.0, 2.0, 2000),
        )
    )
    actual_states, out = rust.states(samples)
    expected = []
    for s in samples:
        state = interpolant.state(float(s))
        expected.append(
            [
                state.xi,
                state.velocity,
                state.beta,
                state.temperature,
                float(np.exp(interpolant.log_pressure(s))),
                state.heating,
                np.nan if state.wall_ratio is None else state.wall_ratio,
            ]
        )
    np.testing.assert_allclose(actual_states, expected, rtol=1e-13, atol=0.0)
    assert not out
    differences = np.abs(actual_states - expected) / np.maximum(
        np.abs(expected), 1e-300
    )
    print(
        f"edge {type(model).__name__} wall={isothermal}: "
        f"max rel={np.nanmax(differences):.3e}; PPoly bitwise"
    )


@pytest.mark.parametrize(
    "case", ["adiabatic", "isothermal", "varying-wall", "keyes", "howarth"]
)
def test_real_full_march_three_paths(
    monkeypatch: pytest.MonkeyPatch, case: str
) -> None:
    from dataclasses import fields

    from test_marching_boundary_layer import howarth_edge

    from aerophysics.exceptions import BoundaryLayerSeparationError
    from aerophysics.marching_boundary_layer import MarchingBoundaryLayerResult

    native_extension = pytest.importorskip("aerophysics_native")
    from aerophysics import _boundary_layer_march as reference

    accepted: list[float] = []
    position = 0.0
    original_state = _EdgeInterpolant.state
    original_resolved = reference._step_is_resolved
    original_march = native_extension.march

    def traced_state(self: _EdgeInterpolant, distance: float) -> StationEdge:
        nonlocal position
        position = distance
        return original_state(self, distance)

    def traced_resolved(*args: Any) -> bool:
        resolved = original_resolved(*args)
        if resolved:
            accepted.append(position)
        return resolved

    def traced_march(*args: Any) -> Any:
        raw = original_march(*args)
        accepted.extend(raw["accepted_positions"][1:])
        return raw

    monkeypatch.setattr(_EdgeInterpolant, "state", traced_state)
    monkeypatch.setattr(reference, "_step_is_resolved", traced_resolved)
    monkeypatch.setattr(native_extension, "march", traced_march)
    monkeypatch.setenv("AEROPHYSICS_NATIVE_JACOBIAN", "analytic")
    edge = smooth_edge(2.0, 3.0, count=51)
    wall: Any = 300.0
    if case == "adiabatic":
        wall = None
    elif case == "varying-wall":
        wall = 290.0 + 30.0 * np.tanh((edge.distance - 1.0) / 0.2)
    elif case == "keyes":
        edge = replace(edge, viscosity_model=AIR_KEYES_VISCOSITY)
    elif case == "howarth":
        # The reference's linear viscosity is deliberately unsupported. Use the
        # same low-Mach Howarth pressure/distance data with supported Sutherland.
        edge = replace(howarth_edge(0.01)[0], viscosity_model=AIR_VISCOSITY)
        wall = None
    results = []
    errors = []
    sequences = []
    for mode, march_mode in (("0", "1"), ("1", "0"), ("1", "1")):
        accepted.clear()
        monkeypatch.setenv("AEROPHYSICS_NATIVE", mode)
        monkeypatch.setenv("AEROPHYSICS_NATIVE_MARCH", march_mode)
        try:
            result = marching_boundary_layer(
                edge, regime=BoundaryLayerRegime.LAMINAR, wall_temperature=wall
            )
        except BoundaryLayerSeparationError as error:
            errors.append(error)
            result = error.result
        assert isinstance(result, MarchingBoundaryLayerResult)
        results.append(result)
        sequences.append(np.array(accepted))
    for sequence in sequences[1:]:
        np.testing.assert_array_equal(sequence, sequences[0])
    print(f"march {case}: identical {len(sequences[0])} accepted steps")
    if case == "howarth":
        assert len(errors) == 3
        for subsequent in errors[1:]:
            assert str(subsequent) == str(errors[0])
            assert subsequent.separation_distance == pytest.approx(
                errors[0].separation_distance, rel=1e-8
            )
        print(
            f"Howarth separation distances: {[e.separation_distance for e in errors]}"
        )
    else:
        assert not errors
    maximum = 0.0
    for result in results[1:]:
        for field in fields(result):
            expected, actual = (
                getattr(results[0], field.name),
                getattr(result, field.name),
            )
            if expected is not None:
                np.testing.assert_allclose(
                    actual, expected, rtol=1e-8, atol=1e-10, err_msg=field.name
                )
                # Near-zero values are assessed by the absolute tolerance above.
                maximum = max(
                    maximum,
                    float(
                        np.max(
                            np.abs(actual - expected)
                            / np.maximum(np.abs(expected), 1e-10)
                        )
                    ),
                )
    print(f"march {case}: max scaled difference={maximum:.3e}")


@pytest.mark.parametrize("temperature", [0.0, -1.0, np.nan, np.inf, -np.inf])
def test_real_edge_invalid_temperature(temperature: float) -> None:
    from aerophysics._native_march import create_edge

    native = pytest.importorskip("aerophysics_native")
    edge = smooth_edge(2.0, 3.0, count=5)
    interpolant = _EdgeInterpolant(edge, None)
    interpolant.total_temperature = temperature
    spec = viscosity_spec(AIR_VISCOSITY)
    assert spec is not None
    rust = create_edge(native, interpolant, spec)
    with np.errstate(invalid="ignore", divide="ignore"):
        with pytest.raises(ValueError) as expected:
            interpolant.state(0.0)
        with pytest.raises(ValueError) as actual:
            rust.states(np.array([0.0]))
    assert str(actual.value) == str(expected.value)


def test_real_edge_keyes_range_flag() -> None:
    from aerophysics._native_march import create_edge

    native = pytest.importorskip("aerophysics_native")
    model = replace(AIR_KEYES_VISCOSITY, temperature_range=(290.0, 310.0))
    edge = replace(smooth_edge(2.0, 3.0, count=5), viscosity_model=model)
    with pytest.warns(ApplicabilityWarning):
        interpolant = _EdgeInterpolant(edge, None)
    spec = viscosity_spec(model)
    assert spec is not None
    _, out = create_edge(native, interpolant, spec).states(edge.distance)
    assert out


@pytest.mark.parametrize(
    "failure",
    ["leading_edge_failed", "invalid_start", "outer_boundary", "not_converged"],
)
def test_real_march_failures_match_reference(
    monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    from aerophysics import _boundary_layer_march as reference
    from aerophysics._native_march import native_march
    from aerophysics.exceptions import BoundaryLayerConvergenceError

    pytest.importorskip("aerophysics_native")
    monkeypatch.setenv("AEROPHYSICS_NATIVE", "1")
    monkeypatch.setenv("AEROPHYSICS_NATIVE_MARCH", "1")
    edge = smooth_edge(2.0, 3.0, count=5)
    interpolant = _EdgeInterpolant(edge, None)
    extent = 2.0 if failure == "outer_boundary" else 16.0
    disc = Discretization(eta_grid(1e-3, 1.015, 0.05, extent), 0.72, AIR_VISCOSITY)
    if failure == "leading_edge_failed":
        disc.newton_iterations = 0
    elif failure == "invalid_start":
        interpolant.xi_at_panels[0] = 1.0
    elif failure == "not_converged":
        monkeypatch.setattr(reference, "_MAX_LOG_VELOCITY_CHANGE", -1.0)
    with pytest.raises(BoundaryLayerConvergenceError) as expected:
        reference.march(interpolant.state, edge.distance[1:], disc)
    with pytest.raises(BoundaryLayerConvergenceError) as actual:
        native_march(interpolant, edge.distance[1:], disc, AIR_VISCOSITY)
    assert str(actual.value) == str(expected.value)


def test_real_success_uses_no_python_edge_states(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pytest.importorskip("aerophysics_native")
    monkeypatch.setenv("AEROPHYSICS_NATIVE", "1")
    monkeypatch.setenv("AEROPHYSICS_NATIVE_MARCH", "1")

    def unexpected(*args: Any) -> Any:
        raise AssertionError("Python edge state called in the full native path")

    monkeypatch.setattr(_EdgeInterpolant, "state", unexpected)
    result = marching_boundary_layer(
        smooth_edge(2.0, 3.0, count=5),
        regime=BoundaryLayerRegime.LAMINAR,
        wall_temperature=300.0,
    )
    assert len(result.distance) == 4
