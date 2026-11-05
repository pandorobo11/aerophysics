"""Independent SA boundary-layer reference and its comparison with VD2."""

import csv
import runpy
import warnings
from collections.abc import Callable
from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics._laminar_similarity import similarity_coefficients
from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    TurbulentModel,
    flat_plate_boundary_layer,
)
from aerophysics.gas import AIR
from aerophysics.transport import AIR_VISCOSITY

_SOLVER = runpy.run_path(
    str(
        Path(__file__).resolve().parents[1]
        / "docs/scripts/build_turbulent_sa_reference.py"
    )
)
TurbulentSABoundaryLayer = cast(Any, _SOLVER["TurbulentSABoundaryLayer"])
reference_solutions = cast(
    Callable[[], list[dict[str, float | str]]], _SOLVER["reference_solutions"]
)

REFERENCE = (
    Path(__file__).parent
    / "reference_data/viscous_flow/turbulent_sa_boundary_layer.csv"
)
QUANTITIES = (
    "cf",
    "delta99_over_x",
    "displacement99_over_x",
    "momentum_over_x",
    "wall_heat_flux_coefficient",
)


def reference_rows() -> list[dict[str, str]]:
    with REFERENCE.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def test_laminar_limit_matches_compressible_similarity() -> None:
    solver = TurbulentSABoundaryLayer(5.0, 220.0, 220.0, edge_sa_ratio=0.0)
    positions, history = solver.march(0.1, steps=20)
    result = solver.quantities(float(positions[-1]), history[-1])
    root = np.sqrt(solver.unit_reynolds * positions[-1])
    expected = similarity_coefficients(
        np.array([[5.0, 220.0, 1.0, 0.72, 1.4]]), AIR_VISCOSITY
    )[0]
    assert_allclose(
        [result.cf * root, result.delta99_over_x * root, result.momentum_over_x * root],
        [expected[0], expected[1], expected[3]],
        rtol=1e-3,
    )


def test_van_driest_wall_temperature_stays_within_sa_reference_bands() -> None:
    for row in reference_rows():
        temperature = float(row["edge_temperature_K"])
        reynolds = float(row["reynolds_number"])
        mach = float(row["mach"])
        viscosity = float(AIR_VISCOSITY.dynamic_viscosity(temperature))
        velocity = mach * np.sqrt(
            AIR.heat_capacity_ratio * AIR.specific_gas_constant * temperature
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = flat_plate_boundary_layer(
                1.0,
                velocity,
                reynolds * viscosity / velocity,
                viscosity,
                regime=BoundaryLayerRegime.TURBULENT,
                compressible=True,
                turbulent_model=TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE,
                mach=mach,
                edge_temperature=temperature,
                wall_temperature=float(row["wall_temperature_K"]),
            )
        friction = float(result.local_skin_friction_coefficient) / float(row["cf"])
        thickness = float(result.boundary_layer_thickness) / float(
            row["delta99_over_x"]
        )
        label = f"{row['case']} Re_x={reynolds:g}"
        # Fully turbulent VD2 versus SA developing from the leading edge.
        assert abs(friction - 1) < (0.03 if reynolds >= 1e7 else 0.08), label
        # At Re_x = 1e6 SA can still carry the thick laminar run from the edge.
        assert abs(thickness - 1) < (0.12 if reynolds >= 1e7 else 0.16), label
        # Fully developed and away from the near-adiabatic St singularity.
        if reynolds >= 1e7 and float(row["wall_to_recovery_temperature_ratio"]) < 1.0:
            assert result.wall_heat_flux is not None
            coefficient = float(result.wall_heat_flux) / (
                reynolds * viscosity * AIR.cp * temperature
            )
            heat_flux = coefficient / float(row["wall_heat_flux_coefficient"])
            assert abs(heat_flux - 1) < 0.06, label


@pytest.mark.parametrize(
    ("mach", "edge_temperature", "wall_ratio"),
    [(0.3, 288.15, 1.2), (2.0, 220.0, 0.5), (5.0, 220.0, 0.25), (10.0, 60.0, 0.2)],
)
def test_laminar_heat_flux_matches_sa_solver_laminar_limit(
    mach: float, edge_temperature: float, wall_ratio: float
) -> None:
    viscosity = float(AIR_VISCOSITY.dynamic_viscosity(edge_temperature))
    velocity = mach * np.sqrt(
        AIR.heat_capacity_ratio * AIR.specific_gas_constant * edge_temperature
    )
    recovery = flat_plate_boundary_layer(
        1.0,
        velocity,
        1e6 * viscosity / velocity,
        viscosity,
        regime=BoundaryLayerRegime.LAMINAR,
        compressible=True,
        mach=mach,
        edge_temperature=edge_temperature,
    ).recovery_temperature
    assert recovery is not None
    wall = wall_ratio * float(recovery)
    result = flat_plate_boundary_layer(
        1.0,
        velocity,
        1e6 * viscosity / velocity,
        viscosity,
        regime=BoundaryLayerRegime.LAMINAR,
        compressible=True,
        mach=mach,
        edge_temperature=edge_temperature,
        wall_temperature=wall,
    )
    solver = TurbulentSABoundaryLayer(mach, edge_temperature, wall, edge_sa_ratio=0.0)
    positions, history = solver.march(0.1, steps=30)
    expected = solver.quantities(
        float(positions[-1]), history[-1]
    ).wall_heat_flux_coefficient
    assert result.wall_heat_flux is not None
    actual = float(result.wall_heat_flux) / (
        1e6 * viscosity * AIR.cp * edge_temperature
    )
    assert_allclose(actual, expected, rtol=0.01)


def test_sa_diffusion_matches_standard_operator(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # With C = 1, standard SA diffusion in Levy--Lees variables is
    # ((g + n/g) n')'/(sigma g) + cb2 n'^2/(sigma g^2). Other SA terms are removed.
    namespace = TurbulentSABoundaryLayer._residual.__globals__
    monkeypatch.setitem(namespace, "CB1", 0.0)
    monkeypatch.setitem(namespace, "CW1", 0.0)
    errors = []
    for nodes in (401, 801):
        eta = np.linspace(0.0, 4.0, nodes)
        solver = TurbulentSABoundaryLayer(5.0, 220.0, 300.0, eta=eta)
        monkeypatch.setattr(solver, "chapman_rubesin", np.ones_like)
        decay = np.exp(-eta)
        g, dg = 1.0 + 3.0 * decay, -3.0 * decay
        n, dn, ddn = 5.0 * (1.0 - decay), 5.0 * decay, -5.0 * decay
        coefficient = g + n / g
        slope = dg + dn / g - n * dg / g**2
        expected = (slope * dn + coefficient * ddn) / (namespace["SIGMA"] * g) + (
            namespace["CB2"] * dn**2 / (namespace["SIGMA"] * g**2)
        )
        state = np.array([0.0 * eta, 0.0 * eta, g, n, 1.0 + eta])
        residual = solver._residual(state, 0.0, 0.0, np.zeros_like(state))
        errors.append(np.max(np.abs(residual[3, 1:-1] - expected[1:-1])))
    assert errors[0] < 2e-3
    assert errors[0] / errors[1] > 3.5


@pytest.mark.generated_assets
def test_sa_reference_is_insensitive_to_outer_boundary() -> None:
    # Low Mach and high Reynolds number put delta99 closest to the outer boundary.
    label, mach, edge_temperature, ratio = _SOLVER["CASES"][0]
    wall = _SOLVER["wall_temperature"](mach, edge_temperature, ratio)
    results = []
    for eta in (None, _SOLVER["grid"](eta_max=300.0)):
        solver = TurbulentSABoundaryLayer(mach, edge_temperature, wall, eta=eta)
        positions, history = solver.march(10.0)
        results.append(solver.quantities(float(positions[-1]), history[-1]))
    # Wall heat flux is near zero on this wall, so it is not compared.
    for name in QUANTITIES[:4]:
        assert_allclose(
            getattr(results[0], name),
            getattr(results[1], name),
            rtol=1e-4,
            err_msg=f"{label}: {name}",
        )


@pytest.mark.generated_assets
def test_sa_reference_regenerates() -> None:
    generated = reference_solutions()
    stored = reference_rows()
    assert len(generated) == len(stored) == 18
    for actual, expected in zip(generated, stored, strict=True):
        assert actual["case"] == expected["case"]
        for key in ("reynolds_number", "wall_temperature_K", *QUANTITIES):
            assert_allclose(
                float(actual[key]),
                float(expected[key]),
                rtol=1e-6,
                err_msg=f"{REFERENCE}: {expected['case']}, {key}",
            )


def test_newton_matches_picard_profiles() -> None:
    newton = TurbulentSABoundaryLayer(5.0, 220.0, 300.0)
    picard = TurbulentSABoundaryLayer(5.0, 220.0, 300.0, use_newton=False)
    positions, coupled = newton.march(0.1, steps=30)
    old_positions, sequential = picard.march(0.1, steps=30)
    assert_allclose(positions, old_positions, rtol=0.0, atol=0.0)
    for name in ("u", "g"):
        assert_allclose(
            getattr(coupled[-1], name),
            getattr(sequential[-1], name),
            rtol=0.0,
            atol=1e-7,
        )
    assert_allclose(
        coupled[-1].n,
        sequential[-1].n,
        rtol=0.0,
        atol=1e-7 * max(1.0, float(sequential[-1].n.max())),
    )
    actual = newton.quantities(float(positions[-1]), coupled[-1])
    expected = picard.quantities(float(old_positions[-1]), sequential[-1])
    assert_allclose(
        [actual.cf, actual.delta99_over_x],
        [expected.cf, expected.delta99_over_x],
        rtol=1e-7,
        atol=0.0,
    )
    assert 0 < newton.newton_iterations < 10 * len(positions)
    assert newton.picard_fallbacks == picard.picard_fallbacks == 0
    assert picard.newton_iterations == 0


def test_colored_jacobian_matches_single_column_differences() -> None:
    # Include wall, edge, trapezoidal constraints and a nonuniform interior.
    solver = TurbulentSABoundaryLayer(
        5.0, 220.0, 300.0, eta=np.concatenate(([0.0], np.geomspace(0.01, 30.0, 12)))
    )
    initial = solver._initial()
    distance = np.concatenate(
        ([0.0], np.cumsum(np.diff(solver.eta) * (initial.g[1:] + initial.g[:-1]) / 2))
    )
    state = np.array([initial.f, initial.u, initial.g, initial.n, distance])
    c0, old = solver._streamwise([initial], [-1.0], 0.01)
    residual = solver._residual(state, 0.01, c0, old)
    bands = solver._jacobian(state, residual, 0.01, c0, old)
    rows = np.arange(state.size)
    for column in range(state.size):
        node, variable = divmod(column, 5)
        step = np.sqrt(np.finfo(float).eps) * max(abs(state[variable, node]), 1.0)
        perturbed = state.copy()
        perturbed[variable, node] += step
        expected = (
            solver._residual(perturbed, 0.01, c0, old) - residual
        ).T.ravel() / step
        actual = np.zeros(state.size)
        inside = abs(rows - column) <= 9
        actual[inside] = bands[9 + rows[inside] - column, column]
        assert_allclose(actual, expected, rtol=1e-10, atol=1e-10)


def test_failed_newton_restarts_picard(monkeypatch: pytest.MonkeyPatch) -> None:
    solver = TurbulentSABoundaryLayer(5.0, 220.0, 300.0)
    initial = solver._initial()
    expected = solver._solve([initial], [-1.0], 0.0)
    monkeypatch.setattr(solver, "_solve_newton", lambda *_: None)
    positions, history = solver.march(0.001, steps=3)
    assert solver.picard_fallbacks == len(positions)
    for name in ("u", "g", "n", "f"):
        assert_allclose(
            getattr(history[0], name), getattr(expected, name), atol=0.0, rtol=0.0
        )
