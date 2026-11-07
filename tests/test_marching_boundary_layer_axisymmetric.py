"""Thin axisymmetric verification, without a transverse-curvature model."""

import importlib
import json
import time
import warnings
from dataclasses import fields

import numpy as np
import pytest
from numpy.typing import ArrayLike
from scipy.interpolate import PchipInterpolator
from test_marching_boundary_layer import (
    smooth_edge,
    tighten_streamwise_control,
)
from test_marching_boundary_layer_sa import constant_edge

from aerophysics import (
    BoundaryLayerEdge,
    BoundaryLayerRegime,
    MarchingBoundaryLayerResult,
    marching_boundary_layer,
)
from aerophysics._array import FloatArray
from aerophysics._boundary_layer_march import Discretization, eta_grid, march
from aerophysics._boundary_layer_sa import TurbulentDiscretization
from aerophysics.exceptions import (
    ApplicabilityWarning,
    BoundaryLayerConvergenceError,
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


def test_explicit_none_is_bit_identical() -> None:
    edge = constant_edge(2.0, 220.0, [0.1])
    implicit = marching_boundary_layer(edge, regime=LAMINAR)
    explicit = marching_boundary_layer(edge, regime=LAMINAR, wall_radius=None)
    for field in fields(implicit):
        np.testing.assert_array_equal(
            getattr(implicit, field.name), getattr(explicit, field.name)
        )


@pytest.mark.parametrize("regime", [LAMINAR, TURBULENT])
@pytest.mark.parametrize("varying", [False, True])
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
    print("V13", regime.name, varying, errors)


@pytest.mark.parametrize("mach", [0.3, 2.0])
@pytest.mark.parametrize("wall", [300.0, None])
@pytest.mark.parametrize("angle", [5.0, 20.0])
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
    print("V14", mach, wall, angle, errors)


@pytest.mark.parametrize("kink", [False, True])
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
    # Residual from interpolating the same ln(p) in two distance variables.
    mid_radius = (radius[:-1] + radius[1:]) / 2
    mid_x = mapped[:-1] + np.diff(s) / 2 * (
        radius[:-1] ** 2 + radius[:-1] * mid_radius + mid_radius**2
    ) / (3 * length**2)
    log_difference = PchipInterpolator(mapped, np.log(pressure))(
        mid_x
    ) - PchipInterpolator(s, np.log(pressure))((s[:-1] + s[1:]) / 2)
    errors["pressure_interpolation_percent"] = float(
        np.max(abs(np.expm1(log_difference))) * 100
    )
    print("V15", kink, errors)


@pytest.mark.parametrize(
    "radius",
    [
        -1.0,
        0.0,
        np.nan,
        np.inf,
        [0, 0.1],
        [-0.1, 0.1, 0.2],
        [0, 0, 0.1],
        [0, -0.1, 0.2],
        [0, np.nan, 1],
        [0, 1, np.inf],
        [[1, 1, 1]],
        [[0.0, 0.1], [0.2]],
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
    normal: bool = False,
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
        grid = (
            eta_grid(0.00035, 1.02, 0.5, 150.0)
            if normal
            else _turbulent_eta_grid(
                max(ip.state(float(s)).reynolds_root for s in targets) ** 2 / 2
            )
        )
        solver = TurbulentDiscretization(grid, 0.72, 0.9, edge.viscosity_model)
    else:
        grid = (
            eta_grid(0.0005, 1.0075, 0.025, 16.0)
            if normal
            else eta_grid(0.001, 1.015, 0.05, 16.0)
        )
        solver = Discretization(grid, 0.72, edge.viscosity_model)
    return _result(
        edge, ip, march(ip.state, targets, solver, relative_step=relative_step), 0.72
    )


def balance_metrics(
    edge: BoundaryLayerEdge,
    radius: FloatArray,
    result: MarchingBoundaryLayerResult,
    *,
    adiabatic: bool,
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
    u, g = result.velocity_ratio, result.temperature_ratio
    deficit = np.trapezoid(
        density[:, None]
        / g
        * velocity[:, None]
        * u
        * (
            edge.gas.cp * temperature[:, None] * (g - 1)
            + velocity[:, None] ** 2 / 2 * (u * u - 1)
        ),
        result.wall_distance,
        axis=1,
    )
    derivative = np.gradient(r * deficit, s) / r
    mask = (s > 0.05) & (s < 1.99)
    metrics = {
        "momentum_median_percent": float(np.median(abs(momentum[mask])) * 100),
        "momentum_max_percent": float(np.max(abs(momentum[mask])) * 100),
    }
    if adiabatic:
        normalized = abs(derivative) / (
            density * velocity * edge.gas.cp * ip.total_temperature
        )
        for name, values in (
            ("enthalpy_W_m2", abs(derivative)),
            ("enthalpy_normalized", normalized),
        ):
            metrics[name + "_median"] = float(np.median(values[mask]))
            metrics[name + "_max"] = float(np.max(values[mask]))
    else:
        error = abs(-derivative / result.wall_heat_flux - 1)
        metrics["enthalpy_median_percent"] = float(np.median(error[mask]) * 100)
        metrics["enthalpy_max_percent"] = float(np.max(error[mask]) * 100)
    return metrics


def flared_radius(s: FloatArray) -> FloatArray:
    """Cone-like tip followed by a smooth flare and waist."""
    return 0.05 * s * (1 + 0.5 * np.exp(-(((s - 0.9) / 0.4) ** 2)))


@pytest.mark.parametrize("regime", [LAMINAR, TURBULENT])
def test_integral_balance(regime: BoundaryLayerRegime) -> None:
    edge = smooth_edge(2.0, 3.0)
    radius = flared_radius(edge.distance)
    result = axisymmetric_march(edge, radius, regime, 300.0, np.linspace(0, 2, 601)[1:])
    metrics = balance_metrics(edge, radius, result, adiabatic=False)
    print("V16", regime.name, metrics)
    assert metrics["momentum_median_percent"] < 0.05


@pytest.mark.generated_assets
@pytest.mark.parametrize("regime", [LAMINAR, TURBULENT])
@pytest.mark.parametrize("wall", [300.0, None])
def test_refinement_and_balance_diagnostics(
    monkeypatch: pytest.MonkeyPatch, regime: BoundaryLayerRegime, wall: float | None
) -> None:
    edge = smooth_edge(2.0, 3.0)
    radius = flared_radius(edge.distance)
    targets = np.linspace(0, 2, 601)[1:]
    base = axisymmetric_march(edge, radius, regime, wall, targets)
    normal = axisymmetric_march(edge, radius, regime, wall, targets, normal=True)
    tighten_streamwise_control(monkeypatch)
    # Eight subdivisions of each output interval, preserving all base stations.
    dense = np.concatenate(
        [
            np.linspace(a, b, 9)[1:]
            for a, b in zip(np.r_[0, targets[:-1]], targets, strict=True)
        ]
    )
    stream = axisymmetric_march(edge, radius, regime, wall, dense, relative_step=0.0075)
    indices = np.searchsorted(dense, targets)
    common = {
        f.name: getattr(stream, f.name)[indices]
        if isinstance(getattr(stream, f.name), np.ndarray) and f.name != "eta"
        else getattr(stream, f.name)
        for f in fields(stream)
    }
    stream = MarchingBoundaryLayerResult(**common)
    for label, result in (("base", base), ("normal", normal), ("stream", stream)):
        metrics = balance_metrics(edge, radius, result, adiabatic=wall is None)
        print("V16", regime.name, wall, label, metrics)
        assert metrics["momentum_median_percent"] < 0.05
    for label, result, limit in (("normal", normal, 0.003), ("stream", stream, 0.001)):
        cf = relative_error(
            base.local_skin_friction_coefficient, result.local_skin_friction_coefficient
        )
        thermal = (
            relative_error(base.wall_temperature, result.wall_temperature)
            if wall is None
            else relative_error(base.wall_heat_flux, result.wall_heat_flux)
        )
        print(
            "V17",
            regime.name,
            wall,
            label,
            {"cf_percent": cf * 100, "thermal_percent": thermal * 100},
        )
        assert cf < limit
        if label == "stream":
            assert thermal < (0.001 if wall is None else 0.005)


def test_separation_uses_axisymmetric_postprocessing() -> None:
    edge = smooth_edge(3.0, 1.0)
    with pytest.raises(BoundaryLayerSeparationError) as caught:
        marching_boundary_layer(
            edge, regime=LAMINAR, wall_temperature=300.0, wall_radius=0.05
        )
    with pytest.raises(BoundaryLayerSeparationError) as planar:
        marching_boundary_layer(edge, regime=LAMINAR, wall_temperature=300.0)
    result = caught.value.result
    reference = planar.value.result
    assert isinstance(result, MarchingBoundaryLayerResult)
    assert isinstance(reference, MarchingBoundaryLayerResult)
    np.testing.assert_allclose(
        result.momentum_thickness, reference.momentum_thickness, rtol=1e-8
    )


def outer_oscillation(n: FloatArray) -> dict[str, float]:
    """Section 11: neighbouring extrema, excluding the monotone outer front."""
    derivative = np.diff(n)
    extrema = np.flatnonzero(derivative[:-1] * derivative[1:] < 0) + 1
    outer = extrema[extrema > max(30, int(np.argmax(n)))]
    amplitudes = []
    for i in range(1, len(outer) - 1):
        at = outer[i]
        neighbor = min((outer[i - 1], outer[i + 1]), key=lambda j: abs(n[j] - n[at]))
        a, b = n[at], n[neighbor]
        if a + b > 0:
            amplitudes.append((abs(a - b) / (a + b), abs(a - b) / 2, (a + b) / 2))
    relative, half, mean = max(amplitudes, default=(0.0, 0.0, 0.0))
    return {
        "amplitude_samples": float(len(amplitudes)),
        "reversals": float(np.count_nonzero(extrema >= 30)),
        "clipped": float(np.count_nonzero(n[1:-1] == 0)),
        "relative_amplitude_percent": float(relative * 100),
        "half_amplitude": float(half),
        "local_mean": float(mean),
    }


def _measure_shoulder(
    monkeypatch: pytest.MonkeyPatch,
    width: float,
    regime: BoundaryLayerRegime,
    wall: float | None,
) -> tuple[MarchingBoundaryLayerResult, set[str]]:
    """Record a numerical stress test, NOT a physical shoulder edge condition.

    Failures are printed with their exception and location before propagation.
    Solver grids, tolerances and step controls are untouched.
    """
    from aerophysics import conical_shock, prandtl_meyer_expansion

    cone = conical_shock(3.0, np.radians(10.0))
    expansion = prandtl_meyer_expansion(float(cone.surface_mach), np.radians(10.0))
    cone_pressure = 2e5 / (1 + 0.2 * 3**2) ** 3.5 * float(cone.surface_pressure_ratio)
    after = cone_pressure * float(expansion.static_pressure_ratio)
    s = np.linspace(0.0, 1.0, 801)
    if width <= 0.001:
        # 64 panels across +/-5 widths give spacing width/6.4 < width/5.
        # Round coordinates only to merge floating-point duplicate nodes.
        local = np.linspace(0.3 - 5 * width, 0.3 + 5 * width, 65)
        s = np.unique(np.round(np.r_[s, local], 15))
        inside = s[(s >= local[0]) & (s <= local[-1])]
        assert np.max(np.diff(inside)) <= width / 5
    radius = np.minimum(s, 0.3) * np.sin(np.radians(10.0))
    pressure = (
        cone_pressure + (after - cone_pressure) * (1 + np.tanh((s - 0.3) / width)) / 2
    )
    edge = BoundaryLayerEdge.from_pressure(
        s,
        pressure,
        total_pressure=2e5 * float(cone.total_pressure_ratio),
        total_temperature=450.0,
    )
    solvers = []

    class CountingSolver(TurbulentDiscretization):
        def __init__(
            self,
            eta: FloatArray,
            prandtl_number: float,
            turbulent_prandtl_number: float,
            viscosity_model: object,
        ) -> None:
            # Capture the existing counter without changing iterations or equations.
            super().__init__(
                eta, prandtl_number, turbulent_prandtl_number, edge.viscosity_model
            )
            solvers.append(self)

    module = importlib.import_module("aerophysics.marching_boundary_layer")
    monkeypatch.setattr(module, "TurbulentDiscretization", CountingSolver)
    result = None
    failure = None
    warning_kinds: set[str] = set()
    start = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            result = marching_boundary_layer(
                edge, regime=regime, wall_temperature=wall, wall_radius=radius
            )
        except (BoundaryLayerConvergenceError, BoundaryLayerSeparationError) as error:
            failure = f"{type(error).__name__}: {error}"
            if isinstance(error, BoundaryLayerSeparationError):
                result = error.result
            raise
        finally:
            seconds = time.perf_counter() - start
            for warning in caught:
                message = str(warning.message)
                if issubclass(warning.category, ApplicabilityWarning):
                    if "transverse curvature" in message:
                        kind = "curvature"
                    elif "edge pressure changes" in message:
                        kind = "pressure_length"
                    elif "y+ > 1" in message:
                        kind = "y_plus"
                    elif "viscous interaction" in message:
                        kind = "viscous_interaction"
                    else:
                        kind = message
                else:
                    kind = f"{warning.category.__name__}: {message}"
                warning_kinds.add(kind)
            record: dict[str, object] = {
                "regime": regime.name,
                "wall": wall,
                "width": width,
                "input_nodes": int(s.size),
                "seconds": seconds,
                "picard": sum(solver.picard_fallbacks for solver in solvers),
                "failure": failure,
                "warning_kinds": sorted(warning_kinds),
                "warnings": [str(w.message) for w in caught],
            }
            if result is not None and result.distance.size:
                ip = _EdgeInterpolant(edge, None, wall_radius=radius)
                assert ip.radius is not None
                pressure_length = result.boundary_layer_thickness * abs(
                    ip.log_pressure_slope(result.distance)
                )
                curvature = result.boundary_layer_thickness / ip.radius(result.distance)
                downstream = np.flatnonzero(
                    (result.distance > 0.25) & (result.distance < 1.0)
                )
                for name, values, indices in (
                    (
                        "pressure_length",
                        pressure_length,
                        np.arange(result.distance.size),
                    ),
                    ("delta_radius", curvature, np.arange(result.distance.size)),
                    ("downstream_delta_radius", curvature, downstream),
                ):
                    if indices.size:
                        at = indices[np.argmax(values[indices])]
                        record["max_" + name] = float(values[at])
                        record["max_" + name + "_s"] = float(result.distance[at])
                if result.nu_tilde_ratio is not None:
                    at = int(np.argmin(abs(result.distance - 0.31)))
                    record["oscillation_after_shoulder_s"] = float(result.distance[at])
                    record["oscillation_after_shoulder"] = outer_oscillation(
                        result.nu_tilde_ratio[at]
                    )
                    record["oscillation_final_s"] = float(result.distance[-1])
                    record["oscillation_final"] = outer_oscillation(
                        result.nu_tilde_ratio[-1]
                    )
            print("V20 " + json.dumps(record, sort_keys=True), flush=True)
    assert result is not None
    return result, warning_kinds


@pytest.mark.generated_assets
@pytest.mark.parametrize("width", [0.02, 0.005, 0.001, 0.0002])
@pytest.mark.parametrize("regime", [LAMINAR, TURBULENT])
@pytest.mark.parametrize("wall", [300.0, None])
def test_shoulder_stress_measurement(
    monkeypatch: pytest.MonkeyPatch,
    width: float,
    regime: BoundaryLayerRegime,
    wall: float | None,
) -> None:
    result, warning_kinds = _measure_shoulder(monkeypatch, width, regime, wall)
    assert result.distance[-1] == 1.0
    expected = {"curvature"}
    if width <= 0.0002 or (width <= 0.001 and regime is TURBULENT):
        expected.add("pressure_length")
    assert warning_kinds == expected


@pytest.mark.generated_assets
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


@pytest.mark.parametrize("regime", [LAMINAR, TURBULENT])
@pytest.mark.parametrize("wall", [300.0, None])
def test_cone_integral_balance(regime: BoundaryLayerRegime, wall: float | None) -> None:
    edge = smooth_edge(2.0, 3.0)
    radius = edge.distance * np.sin(np.radians(10.0))
    result = axisymmetric_march(
        edge, radius, regime, wall, np.linspace(0.0, 2.0, 601)[1:]
    )
    metrics = balance_metrics(edge, radius, result, adiabatic=wall is None)
    print("V16 cone", regime.name, wall, metrics)
    assert metrics["momentum_median_percent"] < 0.05
