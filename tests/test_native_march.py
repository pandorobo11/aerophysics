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
