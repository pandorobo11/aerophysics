"""SA-noft2 verification; CFD agreement is not experimental validation."""

import csv
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from test_marching_boundary_layer import (
    BoundedAirViscosity,
    low_mach_viscosity_edge,
    plate_edge,
    smooth_edge,
)

from aerophysics import (
    AIR,
    BoundaryLayerEdge,
    BoundaryLayerRegime,
    MarchingBoundaryLayerResult,
    marching_boundary_layer,
)
from aerophysics import _boundary_layer_sa as sa_module
from aerophysics._array import FloatArray
from aerophysics._boundary_layer_march import StationEdge, eta_grid, march
from aerophysics._boundary_layer_sa import TurbulentDiscretization
from aerophysics.exceptions import ApplicabilityWarning
from aerophysics.marching_boundary_layer import (
    _SA_ETA_FIRST,
    _SA_ETA_MINIMUM_EXTENT,
    _SA_ETA_RATIO,
    _SA_ETA_SPACING,
    _EdgeInterpolant,
    _turbulent_eta_grid,
)
from aerophysics.transport import AIR_VISCOSITY

REFERENCE = (
    Path(__file__).parent
    / "reference_data/viscous_flow/turbulent_sa_boundary_layer.csv"
)


def discretization() -> TurbulentDiscretization:
    return TurbulentDiscretization(_turbulent_eta_grid(0.0), 0.72, 0.9, AIR_VISCOSITY)


def constant_edge(
    mach: float, temperature: float, targets: list[float], *, unit_reynolds: float = 1e7
) -> BoundaryLayerEdge:
    velocity = mach * np.sqrt(1.4 * AIR.specific_gas_constant * temperature)
    density = unit_reynolds * AIR_VISCOSITY.dynamic_viscosity(temperature) / velocity
    pressure = density * AIR.specific_gas_constant * temperature
    ratio = 1.0 + 0.2 * mach**2
    return BoundaryLayerEdge.from_pressure(
        [0.0, *targets],
        np.full(len(targets) + 1, pressure),
        total_pressure=pressure * ratio**3.5,
        total_temperature=temperature * ratio,
    )


@pytest.mark.parametrize("case", ["M2-AD", "M6-0.25"])
def test_zpg_against_frozen_reference(case: str) -> None:
    with REFERENCE.open() as stream:
        row = next(row for row in csv.DictReader(stream) if row["case"] == case)
    edge = constant_edge(float(row["mach"]), float(row["edge_temperature_K"]), [0.1])
    result = marching_boundary_layer(
        edge,
        regime=BoundaryLayerRegime.TURBULENT,
        wall_temperature=float(row["wall_temperature_K"]),
    )
    for computed, key in (
        (result.local_skin_friction_coefficient[0], "cf"),
        (result.boundary_layer_thickness[0] / 0.1, "delta99_over_x"),
        (result.momentum_thickness[0] / 0.1, "momentum_over_x"),
    ):
        assert computed == pytest.approx(float(row[key]), rel=0.003)
    density = edge.pressure[-1] / (AIR.specific_gas_constant * edge.temperature[-1])
    heat = result.wall_heat_flux[0] / (
        density * edge.velocity[-1] * AIR.cp * edge.temperature[-1]
    )
    reference_heat = float(row["wall_heat_flux_coefficient"])
    scale = float(row["cf"]) if "AD" in case else abs(reference_heat)
    assert abs(heat - reference_heat) < 0.003 * scale


def test_air_viscosity_has_the_reference_chapman_rubesin_factor() -> None:
    g = np.geomspace(0.1, 20.0, 101)
    for temperature in (60.0, 100.0, 220.0, 288.15):
        expected = np.sqrt(g) * (temperature + 110.4) / (temperature * g + 110.4)
        np.testing.assert_allclose(
            discretization().chapman_rubesin(g, temperature), expected, rtol=1e-14
        )


def test_zero_sa_limit_and_read_only_profiles() -> None:
    edge = constant_edge(2.0, 220.0, [0.1, 1.0])
    turbulent = marching_boundary_layer(
        edge,
        regime=BoundaryLayerRegime.TURBULENT,
        wall_temperature=300.0,
        edge_sa_ratio=0.0,
    )
    laminar = marching_boundary_layer(
        edge, regime=BoundaryLayerRegime.LAMINAR, wall_temperature=300.0
    )
    np.testing.assert_allclose(
        turbulent.local_skin_friction_coefficient,
        laminar.local_skin_friction_coefficient,
        rtol=5e-4,
    )
    assert laminar.nu_tilde_ratio is None
    assert turbulent.nu_tilde_ratio is not None
    assert turbulent.nu_tilde_ratio.shape == turbulent.velocity_ratio.shape
    assert np.max(abs(turbulent.nu_tilde_ratio)) == 0.0
    assert not turbulent.nu_tilde_ratio.flags.writeable
    assert not np.shares_memory(turbulent.nu_tilde_ratio, turbulent.velocity_ratio)
    with pytest.raises(ValueError):
        turbulent.nu_tilde_ratio[0, 0] = 1.0


@pytest.mark.parametrize("regime", list(BoundaryLayerRegime))
@pytest.mark.parametrize(
    ("name", "value"),
    [
        ("turbulent_prandtl_number", 0.0),
        ("turbulent_prandtl_number", -1.0),
        ("turbulent_prandtl_number", float("nan")),
        ("edge_sa_ratio", -1.0),
        ("edge_sa_ratio", float("inf")),
    ],
)
def test_new_options_are_validated_for_every_regime(
    regime: BoundaryLayerRegime, name: str, value: float
) -> None:
    with pytest.raises(ValueError, match=name):
        marching_boundary_layer(plate_edge([0.0, 1.0]), regime=regime, **{name: value})


@pytest.mark.solver_verification
def test_edge_viscosity_derivative_and_conservation() -> None:
    edge = smooth_edge(2.0, 3.0)
    interp = _EdgeInterpolant(edge, np.full(edge.distance.size, 300.0), 3.0)
    point, gap = 1.003, 1e-5
    at = interp.state(point)
    left, right = interp.state(point - gap), interp.state(point + gap)
    # n_edge is inversely proportional to nu_e, independently of the analytic slope.
    growth = -2 * at.xi * np.log(right.sa_ratio / left.sa_ratio) / (right.xi - left.xi)
    assert at.viscosity_growth == pytest.approx(growth, rel=1e-6)
    result = marching_boundary_layer(
        edge, regime=BoundaryLayerRegime.TURBULENT, wall_temperature=300.0
    )
    assert result.nu_tilde_ratio is not None
    nu = (
        AIR_VISCOSITY.dynamic_viscosity(edge.temperature)
        * AIR.specific_gas_constant
        * edge.temperature
        / edge.pressure
    )
    conserved = result.nu_tilde_ratio[:, -5:] * nu[1:, None]
    np.testing.assert_allclose(
        conserved, np.full_like(conserved, 3.0 * nu[0]), rtol=0.005
    )
    rho = edge.pressure[1:] / (AIR.specific_gas_constant * edge.temperature[1:])
    velocity = edge.velocity[1:]
    momentum = (
        np.gradient(rho * velocity**2 * result.momentum_thickness, result.distance)
        + rho
        * velocity
        * np.gradient(velocity, result.distance)
        * result.displacement_thickness
    )
    mask = (result.distance > 0.05) & (result.distance < 1.99)
    assert (
        np.median(abs(momentum[mask] / result.wall_shear_stress[mask] - 1.0)) < 0.0005
    )


@pytest.mark.parametrize("sa_ratio", [0.0, 3.0])
def test_turbulent_bounded_viscosity_propagates_range_error(sa_ratio: float) -> None:
    edge = low_mach_viscosity_edge(300.001, BoundedAirViscosity())
    with pytest.raises(ValueError, match="temperature outside 300--1000 K"):
        marching_boundary_layer(
            edge,
            regime=BoundaryLayerRegime.TURBULENT,
            wall_temperature=300.001,
            edge_sa_ratio=sa_ratio,
        )


def test_zero_sa_edge_keeps_viscosity_growth_but_laminar_skips_it() -> None:
    edge = smooth_edge(2.0, 3.0)
    turbulent = _EdgeInterpolant(edge, None, 0.0).state(1.003)
    laminar = _EdgeInterpolant(edge, None).state(1.003)
    assert turbulent.viscosity_growth != 0.0
    assert laminar.viscosity_growth == 0.0
    assert turbulent.sa_ratio == laminar.sa_ratio == 0.0
    assert turbulent.reynolds_root == laminar.reynolds_root


@pytest.mark.parametrize("wall", [1.4, None])
@pytest.mark.parametrize(("beta", "growth"), [(0.1, 0.04), (-0.1, -0.04)])
def test_forced_picard_solves_the_newton_equations(
    monkeypatch: pytest.MonkeyPatch, wall: float | None, beta: float, growth: float
) -> None:
    solver = discretization()
    edge = StationEdge(0.01, 1.0, 0.15, 1.6, 220.0, wall, 2000.0, 0.05, 3.0)
    initial = solver.solve(
        solver.starting_guess(edge), edge, 0.0, np.zeros((5, len(solver.eta)))
    )
    assert initial is not None
    edge = replace(edge, xi=0.0101, beta=beta, viscosity_growth=growth, sa_ratio=2.999)
    coefficient = 10000.0
    lagged = -coefficient * initial.state
    newton = solver.newton(initial.state, edge, coefficient, lagged)
    assert newton is not None
    monkeypatch.setattr(solver, "newton", lambda *_: None)
    picard = solver.solve(initial.state, edge, coefficient, lagged)
    assert picard is not None
    assert solver.picard_fallbacks >= 1
    assert np.min(picard.state[3, 1:-1]) > 0.0
    np.testing.assert_allclose(picard.state, newton.state, atol=1e-6, rtol=1e-8)
    residual = solver.residual(picard.state, edge, coefficient, lagged)
    # Scale transport rows by the local diffusion coefficient of a unit field.
    scaled = residual.copy()
    scaled[1:4, 1:-1] *= solver.width**2
    scaled[3] /= max(1.0, float(np.max(picard.state[3])))
    assert np.max(abs(scaled)) < 1e-5


def test_picard_failure_halves_the_step(monkeypatch: pytest.MonkeyPatch) -> None:
    solver = discretization()
    edge = constant_edge(2.0, 220.0, [0.01])
    interp = _EdgeInterpolant(edge, np.full(2, 300.0), 3.0)
    original_newton, original_picard = solver.newton, solver.picard
    attempted: list[float] = []

    def newton(guess, at, coefficient, lagged):  # type: ignore[no-untyped-def]
        if len(attempted) == 0 and at.xi > 0.0:
            attempted.append(at.xi)
            return None
        if attempted and len(attempted) == 1:
            attempted.append(at.xi)
        return original_newton(guess, at, coefficient, lagged)

    monkeypatch.setattr(solver, "newton", newton)
    monkeypatch.setattr(
        solver,
        "picard",
        lambda *args: None if len(attempted) == 1 else original_picard(*args),
    )
    result = march(interp.state, edge.distance[1:], solver)
    assert result.distance[-1] == 0.01
    assert attempted[1] == pytest.approx(attempted[0] / 2.0)


def test_picard_exhaustion_and_invalid_trials(monkeypatch: pytest.MonkeyPatch) -> None:
    solver = discretization()
    edge = StationEdge(0.0, 1.0, 0.0, 1.6, 220.0, 1.4, 0.0, 0.0, 3.0)
    guess = solver.starting_guess(edge)
    monkeypatch.setattr(sa_module, "_PICARD_ITERATIONS", 0)
    assert solver.picard(guess, edge, 0.0, np.zeros_like(guess)) is None
    monkeypatch.setattr(sa_module, "_PICARD_ITERATIONS", 1)
    guess[2] = -1.0
    assert solver.picard(guess, edge, 0.0, np.zeros_like(guess)) is None


@pytest.mark.solver_verification
def test_adiabatic_wall_matches_a_fixed_recovery_temperature() -> None:
    edge = constant_edge(2.0, 220.0, [1.0])
    adiabatic = marching_boundary_layer(edge, regime=BoundaryLayerRegime.TURBULENT)
    fixed = marching_boundary_layer(
        edge,
        regime=BoundaryLayerRegime.TURBULENT,
        wall_temperature=float(adiabatic.wall_temperature[-1]),
    )
    cold = marching_boundary_layer(
        edge, regime=BoundaryLayerRegime.TURBULENT, wall_temperature=300.0
    )
    assert adiabatic.wall_heat_flux[-1] == 0.0
    assert abs(fixed.wall_heat_flux[-1] / cold.wall_heat_flux[-1]) < 0.005
    assert fixed.local_skin_friction_coefficient[-1] == pytest.approx(
        adiabatic.local_skin_friction_coefficient[-1], rel=0.001
    )


@pytest.mark.solver_verification
def test_strong_compression_separates_with_an_upstream_result() -> None:
    from aerophysics.exceptions import BoundaryLayerSeparationError

    with pytest.raises(BoundaryLayerSeparationError) as caught:
        marching_boundary_layer(
            smooth_edge(3.0, 1.0),
            regime=BoundaryLayerRegime.TURBULENT,
            wall_temperature=300.0,
        )
    error = caught.value
    assert 1.03 < error.separation_distance < 1.05
    assert error.result is not None
    assert error.result.distance[-1] < error.separation_distance
    assert np.all(error.result.wall_shear_stress > 0.0)
    assert error.result.nu_tilde_ratio is not None


def _sa_march(
    edge: BoundaryLayerEdge,
    wall: FloatArray | None,
    targets: FloatArray,
    *,
    eta: FloatArray | None = None,
    relative_step: float = 0.03,
) -> MarchingBoundaryLayerResult:
    from aerophysics.marching_boundary_layer import _result

    interp = _EdgeInterpolant(edge, wall, 3.0)
    if eta is None:
        equivalent = (
            max(interp.state(float(s)).reynolds_root for s in targets) ** 2 / 2.0
        )
        eta = _turbulent_eta_grid(equivalent)
    solver = TurbulentDiscretization(eta, 0.72, 0.9, edge.viscosity_model)
    return _result(
        edge,
        interp,
        march(interp.state, targets, solver, relative_step=relative_step),
        0.72,
    )


@pytest.mark.parametrize(
    ("start", "end", "width"),
    [
        (3.0, 4.0, 0.02),
        (3.0, 2.0, 0.1),
    ],
)
@pytest.mark.solver_verification
def test_sa_streamwise_refinement(
    monkeypatch: pytest.MonkeyPatch, start: float, end: float, width: float
) -> None:
    from test_marching_boundary_layer import tighten_streamwise_control

    edge = smooth_edge(start, end, width=width)
    wall = np.full(edge.distance.size, 300.0)
    base = _sa_march(edge, wall, edge.distance[1:])
    tighten_streamwise_control(monkeypatch)
    targets = np.unique(
        np.concatenate(
            [
                np.linspace(left, right, 9)[1:]
                for left, right in zip(
                    edge.distance[:-1], edge.distance[1:], strict=True
                )
            ]
        )
    )
    fine = _sa_march(edge, wall, targets, relative_step=0.0075)
    ix = np.searchsorted(targets, base.distance)
    assert (
        np.max(
            abs(
                base.local_skin_friction_coefficient
                / fine.local_skin_friction_coefficient[ix]
                - 1.0
            )
        )
        < 0.001
    )
    assert np.max(abs(base.wall_heat_flux / fine.wall_heat_flux[ix] - 1.0)) < 0.005


@pytest.mark.solver_verification
def test_sa_wall_temperature_step_refinement(monkeypatch: pytest.MonkeyPatch) -> None:
    from test_marching_boundary_layer import tighten_streamwise_control

    edge = plate_edge([0.0, 0.5, 0.51, 1.0])
    wall = np.array([300.0, 300.0, 400.0, 400.0])
    base = _sa_march(edge, wall, edge.distance[1:])
    tighten_streamwise_control(monkeypatch)
    targets = np.unique(
        np.concatenate(
            [
                np.linspace(left, right, 81)[1:]
                for left, right in zip(
                    edge.distance[:-1], edge.distance[1:], strict=True
                )
            ]
        )
    )
    fine = _sa_march(edge, wall, targets, relative_step=0.0075)
    ix = np.searchsorted(targets, base.distance)
    assert (
        np.max(
            abs(
                base.local_skin_friction_coefficient
                / fine.local_skin_friction_coefficient[ix]
                - 1.0
            )
        )
        < 0.001
    )
    assert np.max(abs(base.wall_heat_flux / fine.wall_heat_flux[ix] - 1.0)) < 0.005


def test_edge_growth_with_a_non_sutherland_viscosity() -> None:
    from test_marching_boundary_layer import LinearViscosity

    edge = replace(smooth_edge(2.0, 3.0), viscosity_model=LinearViscosity())
    interpolant = _EdgeInterpolant(edge, None, 3.0)
    at = interpolant.state(1.003)
    left, right = interpolant.state(1.003 - 1e-5), interpolant.state(1.003 + 1e-5)
    expected = (
        -2 * at.xi * np.log(right.sa_ratio / left.sa_ratio) / (right.xi - left.xi)
    )
    assert at.viscosity_growth == pytest.approx(expected, rel=1e-6)


def _energy_balance_error(
    edge: BoundaryLayerEdge, result: MarchingBoundaryLayerResult
) -> float:
    from aerophysics.marching_boundary_layer import _isentropic_state

    interpolant = _EdgeInterpolant(edge, None)
    distance = result.distance
    pressure = np.exp(interpolant.log_pressure(distance))
    temperature, velocity = _isentropic_state(
        pressure, interpolant.total_pressure, interpolant.total_temperature, edge.gas
    )
    density = pressure / (AIR.specific_gas_constant * temperature)
    u, g = result.velocity_ratio, result.temperature_ratio
    deficit = np.trapezoid(
        density[:, None]
        / g
        * velocity[:, None]
        * u
        * (
            AIR.cp * temperature[:, None] * (g - 1.0)
            + velocity[:, None] ** 2 / 2.0 * (u**2 - 1.0)
        ),
        result.wall_distance,
        axis=1,
    )
    error = abs(-np.gradient(deficit, distance) / result.wall_heat_flux - 1.0)
    mask = (distance > 0.05) & (distance < 1.99)
    return float(np.median(error[mask]))


@pytest.mark.solver_verification
def test_enthalpy_balance_improves_under_wall_normal_refinement(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from test_marching_boundary_layer import tighten_streamwise_control

    edge = smooth_edge(2.0, 3.0)
    wall = np.full(edge.distance.size, 300.0)
    targets = np.linspace(0.0, 2.0, 601)[1:]
    base = _sa_march(edge, wall, targets)
    normal = _sa_march(
        edge,
        wall,
        targets,
        eta=eta_grid(
            _SA_ETA_FIRST / 2.0,
            1.0 + (_SA_ETA_RATIO - 1.0) / 2.0,
            _SA_ETA_SPACING / 2.0,
            _SA_ETA_MINIMUM_EXTENT,
        ),
    )
    tighten_streamwise_control(monkeypatch)
    streamwise = _sa_march(edge, wall, targets, relative_step=0.0075)
    cf_error = np.max(
        abs(
            base.local_skin_friction_coefficient
            / normal.local_skin_friction_coefficient
            - 1.0
        )
    )
    assert cf_error < 0.003
    base_error = _energy_balance_error(edge, base)
    normal_error = _energy_balance_error(edge, normal)
    streamwise_error = _energy_balance_error(edge, streamwise)
    assert normal_error < 0.5 * base_error
    assert normal_error < streamwise_error


@pytest.mark.parametrize("regime", ["laminar", "turbulent", None])
def test_regime_requires_an_enum_member(regime: str | None) -> None:
    from aerophysics.exceptions import ModelRangeError

    with pytest.raises(ModelRangeError, match=r"LAMINAR.*TURBULENT"):
        marching_boundary_layer(plate_edge([0.0, 1.0]), regime=regime)  # type: ignore[arg-type]


@pytest.mark.parametrize("reynolds", [0.0, 1e6, 1e8, np.nextafter(1e8, np.inf)])
def test_reference_grid_is_unchanged(reynolds: float) -> None:
    expected = eta_grid(
        _SA_ETA_FIRST, _SA_ETA_RATIO, _SA_ETA_SPACING, _SA_ETA_MINIMUM_EXTENT
    )
    np.testing.assert_array_equal(_turbulent_eta_grid(reynolds), expected)


@pytest.mark.solver_verification
@pytest.mark.parametrize(
    ("mach", "temperature", "wall"),
    [(0.3, 288.15, 292.8), (2.0, 220.0, 300.0), (6.0, 100.0, 300.0)],
)
@pytest.mark.parametrize("reynolds", [1e9, 1e10])
def test_high_reynolds_outer_extent(
    mach: float, temperature: float, wall: float, reynolds: float
) -> None:
    edge = constant_edge(mach, temperature, [reynolds / 1e8], unit_reynolds=1e8)
    result = marching_boundary_layer(
        edge, regime=BoundaryLayerRegime.TURBULENT, wall_temperature=wall
    )
    eta99 = float(np.interp(0.99, result.velocity_ratio[-1], result.eta))
    assert result.distance[-1] == reynolds / 1e8
    if reynolds == 1e9:
        refined_grid = eta_grid(
            _SA_ETA_FIRST,
            _SA_ETA_RATIO,
            _SA_ETA_SPACING,
            max(result.eta[-1], 2.0 * eta99),
        )
        refined = _sa_march(
            edge, np.full(edge.distance.size, wall), edge.distance[1:], eta=refined_grid
        )
        assert result.local_skin_friction_coefficient[-1] == pytest.approx(
            refined.local_skin_friction_coefficient[-1], rel=0.001
        )
        assert result.wall_heat_flux[-1] == pytest.approx(
            refined.wall_heat_flux[-1], rel=0.005
        )


@pytest.mark.solver_verification
def test_real_high_reynolds_wall_node_warns() -> None:
    edge = constant_edge(0.3, 288.15, [100.0], unit_reynolds=1e8)
    with pytest.warns(ApplicabilityWarning, match=r"y\+ > 1"):
        result = marching_boundary_layer(
            edge, regime=BoundaryLayerRegime.TURBULENT, wall_temperature=292.8
        )
    wall_density = edge.pressure[-1] / (
        AIR.specific_gas_constant * result.wall_temperature[-1]
    )
    yplus = (
        result.wall_distance[-1, 1]
        * np.sqrt(wall_density * result.wall_shear_stress[-1])
        / AIR_VISCOSITY.dynamic_viscosity(float(result.wall_temperature[-1]))
    )
    assert yplus > 1.0


@pytest.mark.parametrize(("reynolds", "expected"), [(1e3, True), (1e8, False)])
def test_turbulent_leading_edge_viscous_interaction(
    reynolds: float, expected: bool
) -> None:
    import warnings

    edge = constant_edge(10.0, 60.0, [reynolds / 1e7])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ApplicabilityWarning)
        marching_boundary_layer(
            edge, regime=BoundaryLayerRegime.TURBULENT, wall_temperature=300.0
        )
    interaction = [
        item for item in caught if "viscous interaction" in str(item.message)
    ]
    assert bool(interaction) is expected


@pytest.mark.solver_verification
def test_adiabatic_pressure_change_streamwise_refinement(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from test_marching_boundary_layer import tighten_streamwise_control

    edge = smooth_edge(2.0, 3.0)
    base = _sa_march(edge, None, edge.distance[1:])
    tighten_streamwise_control(monkeypatch)
    targets = np.unique(
        np.concatenate(
            [
                np.linspace(left, right, 9)[1:]
                for left, right in zip(
                    edge.distance[:-1], edge.distance[1:], strict=True
                )
            ]
        )
    )
    fine = _sa_march(edge, None, targets, relative_step=0.0075)
    ix = np.searchsorted(targets, base.distance)
    assert (
        np.max(
            abs(
                base.local_skin_friction_coefficient
                / fine.local_skin_friction_coefficient[ix]
                - 1.0
            )
        )
        < 0.001
    )
    assert np.max(abs(base.wall_temperature / fine.wall_temperature[ix] - 1.0)) < 0.001
    assert np.all(base.wall_heat_flux == 0.0)
    assert np.all(fine.wall_heat_flux == 0.0)


@pytest.mark.solver_verification
def test_high_reynolds_compression_keeps_the_outer_boundary_check() -> None:
    from aerophysics.exceptions import BoundaryLayerConvergenceError

    edge = smooth_edge(3.0, 1.5)
    density = edge.pressure / (AIR.specific_gas_constant * edge.temperature)
    local_re = (
        density
        * edge.velocity
        * edge.distance
        / AIR_VISCOSITY.dynamic_viscosity(edge.temperature)
    )
    edge = replace(edge, pressure=edge.pressure * (1e9 / local_re[100]))
    with pytest.raises(
        BoundaryLayerConvergenceError, match="outer boundary is too close"
    ):
        marching_boundary_layer(
            edge, regime=BoundaryLayerRegime.TURBULENT, wall_temperature=300.0
        )
