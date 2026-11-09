"""Thin axisymmetric verification, without a transverse-curvature model."""

import importlib
import warnings
from collections.abc import Callable
from dataclasses import fields

import numpy as np
import pytest
from numpy.typing import ArrayLike
from test_marching_boundary_layer import smooth_edge
from test_marching_boundary_layer_sa import constant_edge

from aerophysics import (
    BoundaryLayerEdge,
    BoundaryLayerRegime,
    MarchingBoundaryLayerResult,
    marching_boundary_layer,
)
from aerophysics._array import FloatArray
from aerophysics._boundary_layer_march import (
    Discretization,
    MarchSolution,
    SeparationError,
    StationEdge,
    eta_grid,
    march,
)
from aerophysics._boundary_layer_sa import TurbulentDiscretization
from aerophysics.exceptions import (
    ApplicabilityWarning,
    BoundaryLayerSeparationError,
)
from aerophysics.marching_boundary_layer import (
    _EdgeInterpolant,
    _isentropic_state,
    _radius_array,
    _result,
    _turbulent_eta_grid,
)

LAMINAR = BoundaryLayerRegime.LAMINAR
TURBULENT = BoundaryLayerRegime.TURBULENT
QUANTITIES = (
    "local_skin_friction_coefficient",
    "wall_heat_flux",
    "boundary_layer_thickness",
    "momentum_thickness",
    "displacement_thickness",
    "wall_distance",
    "velocity_ratio",
    "temperature_ratio",
    "nu_tilde_ratio",
)


def relative_error(actual: FloatArray, expected: FloatArray) -> float:
    """Maximum relative error on nonzero values, requiring identical zeros."""
    nonzero = expected != 0.0
    np.testing.assert_array_equal(actual[~nonzero], expected[~nonzero])
    return float(np.max(abs(actual[nonzero] / expected[nonzero] - 1.0)))


@pytest.mark.parametrize(
    ("varying", "regime"),
    [
        pytest.param(False, LAMINAR, id="False-laminar"),
        pytest.param(
            True, TURBULENT, id="True-turbulent", marks=pytest.mark.solver_verification
        ),
    ],
)
def test_cylinder_is_planar(regime: BoundaryLayerRegime, varying: bool) -> None:
    edge = smooth_edge(2.0, 3.0) if varying else constant_edge(2.0, 220.0, [0.1, 1.0])
    planar = marching_boundary_layer(edge, regime=regime, wall_temperature=300.0)
    cylinder = marching_boundary_layer(
        edge, regime=regime, wall_temperature=300.0, wall_radius=0.05
    )
    errors = {}
    for name in QUANTITIES:
        expected, actual = getattr(planar, name), getattr(cylinder, name)
        if expected is not None:
            if name in ("velocity_ratio", "temperature_ratio", "nu_tilde_ratio"):
                # Relative infinity norm: pointwise relative error at the exact
                # no-slip zero measures roundoff, not profile equivalence.
                errors[name] = float(
                    np.max(abs(actual - expected)) / np.max(abs(expected))
                )
                assert errors[name] < 1e-9
            else:
                errors[name] = relative_error(actual, expected)
                np.testing.assert_allclose(actual, expected, rtol=1e-9, atol=0.0)


@pytest.mark.parametrize(
    ("mach", "wall", "angle"),
    [(0.3, 300.0, 5.0), (2.0, None, 20.0), (2.0, 300.0, 20.0)],
)
def test_sharp_cone_mangler(mach: float, wall: float | None, angle: float) -> None:
    edge = constant_edge(mach, 220.0, [0.1, 1.0])
    plate = marching_boundary_layer(edge, regime=LAMINAR, wall_temperature=wall)
    cone = marching_boundary_layer(
        edge,
        regime=LAMINAR,
        wall_temperature=wall,
        wall_radius=edge.distance * np.sin(np.radians(angle)),
    )
    errors = {}
    for name, factor in (
        ("local_skin_friction_coefficient", np.sqrt(3.0)),
        ("boundary_layer_thickness", 1.0 / np.sqrt(3.0)),
        ("momentum_thickness", 1.0 / np.sqrt(3.0)),
        ("displacement_thickness", 1.0 / np.sqrt(3.0)),
        ("wall_temperature", 1.0),
    ):
        errors[name] = relative_error(
            getattr(cone, name), getattr(plate, name) * factor
        )
        assert errors[name] < 0.0005
    if wall is not None:
        errors["heat"] = relative_error(
            cone.wall_heat_flux, plate.wall_heat_flux * np.sqrt(3)
        )
        assert errors["heat"] < 0.0005
    else:
        np.testing.assert_array_equal(cone.wall_heat_flux, 0.0)


@pytest.mark.parametrize("kink", [False, True])
@pytest.mark.solver_verification
def test_generalized_mangler(kink: bool) -> None:
    s = np.linspace(0.0, 1.0, 801)
    radius = (
        np.minimum(s, 0.2) * np.sin(np.radians(10))
        if kink
        else 0.15 * s * (1 + 0.3 * s)
    )
    length = 0.1
    # Integral of the square of LINEAR radius, exactly on each input interval.
    mapped = np.r_[
        0.0,
        np.cumsum(
            np.diff(s)
            * (radius[:-1] ** 2 + radius[:-1] * radius[1:] + radius[1:] ** 2)
            / (3 * length**2)
        ),
    ]
    pressure = 2e4 * np.exp(-0.1 * np.sin(np.pi * s) ** 4)
    edge = BoundaryLayerEdge.from_pressure(
        s, pressure, total_pressure=1.6e5, total_temperature=500.0
    )
    equivalent = BoundaryLayerEdge.from_pressure(
        mapped, pressure, total_pressure=1.6e5, total_temperature=500.0
    )
    axisym = marching_boundary_layer(
        edge, regime=LAMINAR, wall_temperature=300.0, wall_radius=radius
    )
    planar = marching_boundary_layer(equivalent, regime=LAMINAR, wall_temperature=300.0)
    # xi_axi = L^2 xi_planar, R_axi = (L/r) R_planar. Thus Cf and q
    # acquire r/L, while all physical thicknesses acquire L/r. The laminar
    # transformed PDE is invariant; the SA source contains R and is not.
    errors = {}
    for name, factor in (
        ("local_skin_friction_coefficient", radius[1:] / length),
        ("wall_heat_flux", radius[1:] / length),
        ("momentum_thickness", length / radius[1:]),
    ):
        errors[name] = relative_error(
            getattr(axisym, name), getattr(planar, name) * factor
        )
        assert errors[name] < 0.001


@pytest.mark.parametrize(
    "radius",
    [
        -1.0,
        0.0,
        np.nan,
        [0, 0.1],
        [-0.1, 0.1, 0.2],
        [0, 0, 0.1],
        [[1, 1, 1]],
        "bad",
    ],
)
def test_invalid_radius(radius: ArrayLike) -> None:
    with pytest.raises(ValueError, match="wall_radius must"):
        marching_boundary_layer(
            constant_edge(2.0, 220.0, [0.1, 1.0]), regime=LAMINAR, wall_radius=radius
        )


def test_radius_ownership_linear_interpolation_and_tip_limits() -> None:
    edge = constant_edge(2.0, 220.0, [0.2, 1.0])
    radius = np.array([0.0, 0.02, 0.02])
    owned = _radius_array(radius, 3)
    assert not np.shares_memory(radius, owned)
    ip = _EdgeInterpolant(edge, None, 3.0, owned)
    before = ip.state(0.1)
    result = marching_boundary_layer(edge, regime=LAMINAR, wall_radius=radius)
    saved = {
        f.name: getattr(result, f.name).copy()
        for f in fields(result)
        if isinstance(getattr(result, f.name), np.ndarray)
    }
    radius[:] = 0.8
    owned[:] = 0.9
    assert ip.state(0.1) == before
    assert ip.radius is not None
    assert float(ip.radius(0.1)) == pytest.approx(0.01)
    assert float(ip.radius(0.6)) == pytest.approx(0.02)
    start = ip.state(0.0)
    assert (start.xi, start.beta, start.viscosity_growth, start.reynolds_root) == (
        0,
        0,
        0,
        0,
    )
    assert before.reynolds_root**2 / 2 == pytest.approx(1e7 * 0.1 / 3, rel=1e-12)
    for name, expected in saved.items():
        actual = getattr(result, name)
        np.testing.assert_array_equal(actual, expected)
        assert not np.shares_memory(actual, radius)
        assert not actual.flags.writeable


@pytest.mark.parametrize(
    "radius, expected", [(None, False), (1.0, False), (np.array([0, 1e-6, 0.01]), True)]
)
def test_curvature_warning(radius: ArrayLike | None, expected: bool) -> None:
    edge = constant_edge(2.0, 220.0, [1e-4, 1.0])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ApplicabilityWarning)
        marching_boundary_layer(edge, regime=LAMINAR, wall_radius=radius)
    curvature = [w for w in caught if "transverse curvature" in str(w.message)]
    assert bool(curvature) is expected
    if curvature:
        assert "maximum at s = 0.0001 m" in str(curvature[0].message)


def axisymmetric_march(
    edge: BoundaryLayerEdge,
    radius: FloatArray,
    regime: BoundaryLayerRegime,
    wall: float | None,
    targets: FloatArray,
    *,
    relative_step: float = 0.03,
) -> MarchingBoundaryLayerResult:
    """V8/V9 methodology: fixed edge interpolation and common output stations."""
    ip = _EdgeInterpolant(
        edge,
        None if wall is None else np.full(edge.distance.size, wall),
        3.0 if regime is TURBULENT else None,
        radius,
    )
    solver: Discretization
    if regime is TURBULENT:
        grid = _turbulent_eta_grid(
            max(ip.state(float(s)).reynolds_root for s in targets) ** 2 / 2
        )
        solver = TurbulentDiscretization(grid, 0.72, 0.9, edge.viscosity_model)
    else:
        grid = eta_grid(0.001, 1.015, 0.05, 16.0)
        solver = Discretization(grid, 0.72, edge.viscosity_model)
    return _result(
        edge, ip, march(ip.state, targets, solver, relative_step=relative_step), 0.72
    )


def balance_metrics(
    edge: BoundaryLayerEdge,
    radius: FloatArray,
    result: MarchingBoundaryLayerResult,
) -> dict[str, float]:
    ip = _EdgeInterpolant(edge, None, wall_radius=radius)
    s = result.distance
    assert ip.radius is not None
    r = ip.radius(s)
    pressure = np.exp(ip.log_pressure(s))
    temperature, velocity = _isentropic_state(
        pressure, ip.total_pressure, ip.total_temperature, edge.gas
    )
    density = pressure / (edge.gas.specific_gas_constant * temperature)
    momentum = (
        np.gradient(r * density * velocity**2 * result.momentum_thickness, s)
        + r
        * density
        * velocity
        * np.gradient(velocity, s)
        * result.displacement_thickness
    ) / (r * result.wall_shear_stress) - 1
    mask = (s > 0.05) & (s < 1.99)
    return {
        "momentum_median_percent": float(np.median(abs(momentum[mask])) * 100),
    }


def flared_radius(s: FloatArray) -> FloatArray:
    """Cone-like tip followed by a smooth flare and waist."""
    return 0.05 * s * (1 + 0.5 * np.exp(-(((s - 0.9) / 0.4) ** 2)))


@pytest.mark.parametrize(
    ("cone", "regime"),
    [
        pytest.param(False, LAMINAR, id="flared-laminar"),
        pytest.param(True, TURBULENT, id="cone-turbulent"),
    ],
)
@pytest.mark.solver_verification
def test_axisymmetric_momentum_balance(cone: bool, regime: BoundaryLayerRegime) -> None:
    edge = smooth_edge(2.0, 3.0)
    radius = (
        edge.distance * np.sin(np.radians(10.0))
        if cone
        else flared_radius(edge.distance)
    )
    result = axisymmetric_march(edge, radius, regime, 300.0, np.linspace(0, 2, 601)[1:])
    metrics = balance_metrics(edge, radius, result)
    assert metrics["momentum_median_percent"] < 0.05


def test_separation_postprocessing_on_a_short_constant_edge(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Check exception postprocessing after one real station and forced separation."""

    def separated(
        edge_at: Callable[[float], StationEdge],
        targets: FloatArray,
        discretization: Discretization,
    ) -> MarchSolution:
        solution = march(edge_at, targets[:1], discretization)
        raise SeparationError(float(targets[0]) * 1.01, solution)

    module = importlib.import_module("aerophysics.marching_boundary_layer")
    monkeypatch.setattr(module, "march", separated)
    edge = constant_edge(2.0, 220.0, [0.1, 1.0])
    with pytest.raises(BoundaryLayerSeparationError) as cylinder:
        marching_boundary_layer(
            edge, regime=LAMINAR, wall_temperature=300.0, wall_radius=0.05
        )
    with pytest.raises(BoundaryLayerSeparationError) as planar:
        marching_boundary_layer(edge, regime=LAMINAR, wall_temperature=300.0)
    assert cylinder.value.separation_distance == pytest.approx(0.1 * 1.01)
    assert isinstance(cylinder.value.__cause__, SeparationError)
    result, reference = cylinder.value.result, planar.value.result
    assert isinstance(result, MarchingBoundaryLayerResult)
    assert isinstance(reference, MarchingBoundaryLayerResult)
    np.testing.assert_array_equal(result.distance, [0.1])
    assert np.all(result.wall_shear_stress > 0.0)
    np.testing.assert_allclose(
        result.momentum_thickness, reference.momentum_thickness, rtol=1e-8
    )


@pytest.mark.solver_verification
def test_sa_cone_ratio_is_informational() -> None:
    edge = constant_edge(2.0, 220.0, [0.1, 1.0, 10.0])
    cone = marching_boundary_layer(
        edge,
        regime=TURBULENT,
        wall_temperature=300.0,
        wall_radius=edge.distance * np.sin(np.radians(10.0)),
    )
    plate = marching_boundary_layer(edge, regime=TURBULENT, wall_temperature=300.0)
    print(
        "V15 SA Re=1e6,1e7,1e8",
        cone.local_skin_friction_coefficient / plate.local_skin_friction_coefficient,
    )
    # Freeze the measured numerical regression, not a physical acceptance criterion.
    np.testing.assert_allclose(
        cone.local_skin_friction_coefficient / plate.local_skin_friction_coefficient,
        [1.17608213, 1.15882533, 1.13586288],
        rtol=1e-4,
        atol=0.0,
    )
