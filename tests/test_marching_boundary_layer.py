"""Verification of the marching boundary layer against independent solutions."""

import warnings
from dataclasses import FrozenInstanceError
from typing import overload

import numpy as np
import pytest
from numpy.typing import ArrayLike
from scipy.integrate import cumulative_trapezoid, solve_bvp

from aerophysics import (
    AIR,
    BoundaryLayerEdge,
    BoundaryLayerRegime,
    MarchingBoundaryLayerResult,
    flat_plate_boundary_layer,
    marching_boundary_layer,
)
from aerophysics._array import FloatArray
from aerophysics._boundary_layer_march import StationEdge, eta_grid, march
from aerophysics.exceptions import (
    ApplicabilityWarning,
    BoundaryLayerConvergenceError,
    BoundaryLayerSeparationError,
    ModelRangeError,
)
from aerophysics.marching_boundary_layer import _EdgeInterpolant
from aerophysics.transport import AIR_VISCOSITY

GAMMA = AIR.heat_capacity_ratio
GAS_CONSTANT = AIR.specific_gas_constant
ETA = eta_grid(1e-3, 1.015, 0.05, 16.0)
# Howarth (1938) linear retardation: Hartree (ARC R&M 2426, eq. 37) gives
# x_s = 0.9589 for U = 1 - x/8.
HOWARTH_SEPARATION = 0.9589 / 8.0
# NACA Report 1293, Table II, for Pr = 1 and viscosity linear in temperature:
# (beta, S_w, f''(0), S'(0)), with g = 1 + S without dissipation.
COHEN_RESHOTKO_TABLE_II = (
    (0.0, 0.0, 0.4696, 0.0),
    (1.0, 0.0, 1.2326, 0.0),
    (2.0, 0.0, 1.6870, 0.0),
    (0.5, 1.0, 1.2351, -0.5725),
    (2.0, 1.0, 2.4878, -0.6613),
)


class LinearViscosity:
    """``mu`` proportional to ``T``, which makes ``C = 1`` for a thin layer."""

    @overload
    def dynamic_viscosity(self, temperature: float) -> float: ...

    @overload
    def dynamic_viscosity(self, temperature: ArrayLike) -> FloatArray: ...

    def dynamic_viscosity(self, temperature: ArrayLike) -> float | FloatArray:
        value = np.asarray(temperature, dtype=np.float64) * 1e-7
        return float(value) if np.ndim(temperature) == 0 else value


def similarity_reference(beta: float, wall_ratio: float) -> tuple[float, float]:
    """Solve the similar equations with ``C = Pr = 1`` and no dissipation.

    This is the form of Cohen and Reshotko (1956), with ``g = 1 + S``. For
    ``wall_ratio = 1`` it is the Falkner--Skan equation.
    """

    def rhs(_: FloatArray, y: FloatArray) -> FloatArray:
        return np.vstack(
            (y[1], y[2], -y[0] * y[2] - beta * (y[3] - y[1] ** 2), y[4], -y[0] * y[4])
        )

    def boundary(wall: FloatArray, edge: FloatArray) -> FloatArray:
        return np.array(
            (wall[0], wall[1], edge[1] - 1.0, wall[3] - wall_ratio, edge[3] - 1.0)
        )

    x = np.linspace(0.0, 12.0, 400)
    decay = np.exp(-x)
    guess = np.vstack(
        (
            x - 1.0 + decay,
            1.0 - decay,
            decay,
            1.0 + (wall_ratio - 1.0) * decay,
            (1.0 - wall_ratio) * decay,
        )
    )
    solution = solve_bvp(rhs, boundary, x, guess, tol=1e-10, max_nodes=200000)
    assert solution.success
    return float(solution.y[2, 0]), float(solution.y[4, 0])


def wedge_slopes(beta: float, wall_ratio: float) -> tuple[float, float]:
    solution = march(
        lambda xi: StationEdge(xi, 1.0, beta, 0.0, 300.0, wall_ratio),
        np.array([1.0]),
        eta=ETA,
        prandtl_number=1.0,
        viscosity_model=LinearViscosity(),
    )
    return (
        float(solution.wall_velocity_slope[-1]),
        float(solution.wall_temperature_slope[-1]),
    )


@pytest.mark.parametrize("beta", [-0.15, -0.1, 0.0, 0.1, 0.5, 1.0])
def test_falkner_skan_wall_shear_matches_an_independent_solution(beta: float) -> None:
    expected, _ = similarity_reference(beta, 1.0)
    assert wedge_slopes(beta, 1.0)[0] == pytest.approx(expected, rel=5e-4)


@pytest.mark.parametrize("beta", [-0.1, 0.0, 0.5])
@pytest.mark.parametrize("wall_ratio", [0.2, 2.0])
def test_similar_wall_heat_transfer_matches_an_independent_solution(
    beta: float, wall_ratio: float
) -> None:
    shear, heat = similarity_reference(beta, wall_ratio)
    marched_shear, marched_heat = wedge_slopes(beta, wall_ratio)
    assert marched_shear == pytest.approx(shear, rel=5e-4)
    assert marched_heat == pytest.approx(heat, rel=5e-4)


@pytest.mark.parametrize(
    ("beta", "wall_enthalpy", "shear", "heat"), COHEN_RESHOTKO_TABLE_II
)
def test_similar_solutions_match_the_published_table(
    beta: float, wall_enthalpy: float, shear: float, heat: float
) -> None:
    marched_shear, marched_heat = wedge_slopes(beta, 1.0 + wall_enthalpy)
    assert marched_shear == pytest.approx(shear, rel=2e-3)
    assert marched_heat == pytest.approx(heat, rel=2e-3, abs=1e-8)


def howarth_edge(
    mach0: float, count: int = 161
) -> tuple[BoundaryLayerEdge, FloatArray]:
    """Build the Illingworth--Stewartson image of ``M = M0 (1 - X/8)``.

    With ``C = Pr = 1`` and an adiabatic wall, the compressible equations reduce
    to the incompressible ones in ``X`` with ``U proportional to M``, where
    ``d(xi)`` is proportional to ``M dX``.
    """
    transformed = np.linspace(0.0, 1.0, count)
    mach = mach0 * (1.0 - transformed / 8.0)
    total_temperature, total_pressure = 300.0, 1e5
    temperature = total_temperature / (1.0 + (GAMMA - 1.0) / 2.0 * mach**2)
    pressure = total_pressure * (temperature / total_temperature) ** (
        GAMMA / (GAMMA - 1.0)
    )
    density = pressure / (GAS_CONSTANT * temperature)
    velocity = mach * np.sqrt(GAMMA * GAS_CONSTANT * temperature)
    model = LinearViscosity()
    distance = cumulative_trapezoid(
        mach / (density * model.dynamic_viscosity(temperature) * velocity),
        transformed,
        initial=0,
    )
    distance /= distance[-1]
    edge = BoundaryLayerEdge.from_pressure(
        distance,
        pressure,
        total_pressure=total_pressure,
        total_temperature=total_temperature,
        viscosity_model=model,
    )
    return edge, transformed


@pytest.mark.parametrize("mach0", [0.01, 3.0])
def test_howarth_separation_and_its_compressible_image(mach0: float) -> None:
    edge, transformed = howarth_edge(mach0)
    with pytest.raises(BoundaryLayerSeparationError) as caught:
        marching_boundary_layer(
            edge, regime=BoundaryLayerRegime.LAMINAR, prandtl_number=1.0
        )
    separation = np.interp(caught.value.separation_distance, edge.distance, transformed)
    assert separation / 8.0 == pytest.approx(HOWARTH_SEPARATION, rel=2e-3)

    upstream = caught.value.result
    assert isinstance(upstream, MarchingBoundaryLayerResult)
    assert upstream.distance[-1] <= caught.value.separation_distance
    assert np.all(upstream.pressure_gradient_parameter < 0.0)
    assert np.all(np.diff(upstream.local_skin_friction_coefficient[20:]) < 0.0)


def test_adiabatic_unit_prandtl_layer_keeps_the_exact_total_enthalpy() -> None:
    edge, _ = howarth_edge(3.0)
    with pytest.raises(BoundaryLayerSeparationError) as caught:
        marching_boundary_layer(
            edge, regime=BoundaryLayerRegime.LAMINAR, prandtl_number=1.0
        )
    result = caught.value.result
    index = np.arange(result.distance.size)[1:]
    velocity_ratio = result.velocity_ratio[index]
    temperature_ratio = result.temperature_ratio[index]
    temperature = edge.temperature[index + 1]
    m = 300.0 / temperature - 1.0
    exact = 1.0 + m[:, np.newaxis] * (1.0 - velocity_ratio**2)
    assert np.max(np.abs(temperature_ratio / exact - 1.0)) < 5e-4
    assert np.all(result.wall_heat_flux == 0.0)
    assert result.wall_temperature[index] == pytest.approx(300.0, rel=5e-4)


@pytest.mark.parametrize(
    ("mach", "edge_temperature", "wall_temperature"),
    [
        (0.1, 300.0, None),
        (2.0, 220.0, None),
        (6.0, 100.0, None),
        (6.0, 100.0, 400.0),
        (10.0, 60.0, 300.0),
    ],
)
def test_constant_edge_matches_the_compressible_similarity_solution(
    mach: float, edge_temperature: float, wall_temperature: float | None
) -> None:
    pressure = 2e4
    ratio = 1.0 + (GAMMA - 1.0) / 2.0 * mach**2
    edge = BoundaryLayerEdge.from_pressure(
        [0.0, 0.5, 1.0],
        np.full(3, pressure),
        total_pressure=pressure * ratio ** (GAMMA / (GAMMA - 1.0)),
        total_temperature=edge_temperature * ratio,
    )
    density = pressure / (GAS_CONSTANT * edge_temperature)
    velocity = mach * np.sqrt(GAMMA * GAS_CONSTANT * edge_temperature)
    viscosity = float(AIR_VISCOSITY.dynamic_viscosity(edge_temperature))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ApplicabilityWarning)
        result = marching_boundary_layer(
            edge,
            regime=BoundaryLayerRegime.LAMINAR,
            wall_temperature=wall_temperature,
        )
        reference = flat_plate_boundary_layer(
            1.0,
            velocity,
            density,
            viscosity,
            regime=BoundaryLayerRegime.LAMINAR,
            compressible=True,
            mach=mach,
            edge_temperature=edge_temperature,
            wall_temperature=wall_temperature,
        )
    assert result.local_skin_friction_coefficient[-1] == pytest.approx(
        reference.local_skin_friction_coefficient, rel=5e-4
    )
    assert result.momentum_thickness[-1] == pytest.approx(
        reference.momentum_thickness, rel=5e-4
    )
    assert result.displacement_thickness[-1] == pytest.approx(
        reference.displacement_thickness, rel=5e-4
    )
    assert result.boundary_layer_thickness[-1] == pytest.approx(
        reference.boundary_layer_thickness, rel=5e-4
    )
    assert result.reynolds_number[-1] == pytest.approx(reference.reynolds_number)
    assert result.wall_shear_stress[-1] == pytest.approx(
        reference.wall_shear_stress, rel=5e-4
    )
    assert result.pressure_gradient_parameter == pytest.approx(0.0, abs=1e-12)


def smooth_edge(
    mach_start: float, mach_end: float, count: int = 201, width: float = 0.1
) -> BoundaryLayerEdge:
    distance = np.linspace(0.0, 2.0, count)
    mach = (
        mach_start
        + (mach_end - mach_start) * (1.0 + np.tanh((distance - 1.0) / width)) / 2.0
    )
    total_temperature = 450.0
    total_pressure = 2e5
    pressure = total_pressure * (1.0 + (GAMMA - 1.0) / 2.0 * mach**2) ** (
        -GAMMA / (GAMMA - 1.0)
    )
    return BoundaryLayerEdge.from_pressure(
        distance,
        pressure,
        total_pressure=total_pressure,
        total_temperature=total_temperature,
    )


@pytest.mark.parametrize(("start", "end"), [(2.0, 3.0), (3.0, 2.95)])
def test_integral_momentum_and_energy_balances(start: float, end: float) -> None:
    edge = smooth_edge(start, end)
    result = marching_boundary_layer(
        edge, regime=BoundaryLayerRegime.LAMINAR, wall_temperature=300.0
    )
    distance = result.distance
    density = edge.pressure[1:] / (GAS_CONSTANT * edge.temperature[1:])
    velocity = edge.velocity[1:]
    mask = distance > 0.05

    momentum_flux = np.gradient(
        density * velocity**2 * result.momentum_thickness, distance
    )
    pressure_term = (
        density
        * velocity
        * np.gradient(velocity, distance)
        * result.displacement_thickness
    )
    momentum = (momentum_flux + pressure_term) / result.wall_shear_stress
    assert np.median(np.abs(momentum[mask] - 1.0)) < 2e-3

    # d/ds int rho u (H - H_e) dy = -q_w, with q_w positive into the wall.
    enthalpy = AIR.cp * edge.temperature[1:]
    velocity_ratio = result.velocity_ratio
    deficit = np.array(
        [
            np.trapezoid(
                density[i]
                / result.temperature_ratio[i]
                * velocity[i]
                * velocity_ratio[i]
                * (
                    enthalpy[i] * (result.temperature_ratio[i] - 1.0)
                    + velocity[i] ** 2 / 2.0 * (velocity_ratio[i] ** 2 - 1.0)
                ),
                result.wall_distance[i],
            )
            for i in range(distance.size)
        ]
    )
    energy = -np.gradient(deficit, distance) / result.wall_heat_flux
    assert np.median(np.abs(energy[mask] - 1.0)) < 2e-3


def test_streamwise_step_refinement_converges() -> None:
    edge = smooth_edge(2.0, 3.0)
    interpolant = _EdgeInterpolant(edge, np.full(edge.distance.size, 300.0))
    slopes = [
        march(
            interpolant.state,
            edge.distance[1:],
            eta=ETA,
            prandtl_number=0.72,
            viscosity_model=edge.viscosity_model,
            relative_step=step,
        ).wall_velocity_slope
        for step in (0.03, 0.015)
    ]
    assert np.max(np.abs(slopes[0] / slopes[1] - 1.0)) < 1e-3


def test_wall_normal_grid_refinement_converges() -> None:
    edge = smooth_edge(2.0, 3.0)
    interpolant = _EdgeInterpolant(edge, np.full(edge.distance.size, 300.0))
    coarse, fine = (
        march(
            interpolant.state,
            edge.distance[1:],
            eta=grid,
            prandtl_number=0.72,
            viscosity_model=edge.viscosity_model,
        )
        for grid in (ETA, eta_grid(5e-4, 1.0075, 0.025, 16.0))
    )
    assert (
        np.max(np.abs(coarse.wall_velocity_slope / fine.wall_velocity_slope - 1.0))
        < 2e-4
    )
    assert (
        np.max(
            np.abs(coarse.wall_temperature_slope / fine.wall_temperature_slope - 1.0)
        )
        < 1e-3
    )


def test_output_stations_do_not_change_the_solution() -> None:
    def layer(count: int) -> MarchingBoundaryLayerResult:
        distance = np.linspace(0.0, 1.0, count)
        pressure = 2e4 * (1.0 - 0.1 * distance**2)
        return marching_boundary_layer(
            BoundaryLayerEdge.from_pressure(
                distance, pressure, total_pressure=1.6e5, total_temperature=500.0
            ),
            regime=BoundaryLayerRegime.LAMINAR,
        )

    sparse, dense = layer(21), layer(81)
    assert sparse.local_skin_friction_coefficient[-1] == pytest.approx(
        dense.local_skin_friction_coefficient[-1], rel=2e-3
    )


def test_history_matters_for_the_same_local_edge_state() -> None:
    """Downstream of a pressure change a layer differs from the local plate."""
    edge = smooth_edge(2.0, 3.0, width=0.02)
    result = marching_boundary_layer(
        edge, regime=BoundaryLayerRegime.LAMINAR, wall_temperature=300.0
    )
    index = int(np.argmin(np.abs(edge.distance[1:] - 1.2)))
    plate = marching_boundary_layer(
        BoundaryLayerEdge.from_pressure(
            [0.0, edge.distance[index + 1]],
            np.full(2, edge.pressure[index + 1]),
            total_pressure=2e5,
            total_temperature=450.0,
        ),
        regime=BoundaryLayerRegime.LAMINAR,
        wall_temperature=300.0,
    )
    ratio = (
        result.local_skin_friction_coefficient[index]
        / plate.local_skin_friction_coefficient[-1]
    )
    assert abs(ratio - 1.0) > 0.05


def test_wall_temperature_distribution_is_followed() -> None:
    distance = np.linspace(0.0, 1.0, 21)
    wall = 300.0 + 100.0 * distance
    edge = BoundaryLayerEdge.from_pressure(
        distance, np.full(21, 2e4), total_pressure=1.6e5, total_temperature=500.0
    )
    result = marching_boundary_layer(
        edge, regime=BoundaryLayerRegime.LAMINAR, wall_temperature=wall
    )
    assert result.wall_temperature == pytest.approx(wall[1:], rel=1e-10)
    assert result.wall_heat_flux[-1] > 0.0
    assert result.velocity_ratio.shape == (20, ETA.size)
    assert result.wall_distance.shape == result.velocity_ratio.shape
    assert result.eta == pytest.approx(ETA)


def test_results_and_edges_own_read_only_arrays() -> None:
    distance = np.array([0.0, 0.5, 1.0])
    pressure = np.full(3, 2e4)
    edge = BoundaryLayerEdge.from_pressure(
        distance, pressure, total_pressure=1.6e5, total_temperature=500.0
    )
    distance[1] = 0.9
    pressure[1] = 1.0
    assert edge.distance[1] == 0.5
    assert edge.pressure[1] == 2e4
    assert not edge.temperature.flags.writeable
    result = marching_boundary_layer(edge, regime=BoundaryLayerRegime.LAMINAR)
    assert result.distance.tolist() == [0.5, 1.0]
    for values in (result.distance, result.velocity_ratio, result.eta):
        assert not values.flags.writeable
    with pytest.raises(FrozenInstanceError):
        edge.gas = AIR  # type: ignore[misc]


@pytest.mark.parametrize(
    ("distance", "pressure", "message"),
    [
        ([0.0], [2e4], "two or more"),
        (0.0, 2e4, "two or more"),
        ([0.1, 1.0], [2e4, 2e4], "leading edge"),
        ([0.0, 1.0, 0.5], [2e4, 2e4, 2e4], "increase"),
        ([0.0, 1.0], [2e4, 2e4, 2e4], "shape"),
        ([0.0, 1.0], [2e4, -1.0], "greater than zero"),
        ([0.0, 1.0], [2e4, 1.6e5], "below total_pressure"),
    ],
)
def test_invalid_edges_are_rejected(
    distance: ArrayLike, pressure: ArrayLike, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        BoundaryLayerEdge.from_pressure(
            distance, pressure, total_pressure=1.6e5, total_temperature=500.0
        )


def test_invalid_totals_are_rejected() -> None:
    with pytest.raises(ValueError, match="total_pressure"):
        BoundaryLayerEdge.from_pressure(
            [0.0, 1.0], [2e4, 2e4], total_pressure=-1.0, total_temperature=500.0
        )
    with pytest.raises(ValueError, match="total_temperature"):
        BoundaryLayerEdge.from_pressure(
            [0.0, 1.0],
            [2e4, 2e4],
            total_pressure=1.6e5,
            total_temperature=[500.0],  # type: ignore[arg-type]
        )


def test_direct_edges_must_be_isentropic() -> None:
    valid = BoundaryLayerEdge.from_pressure(
        [0.0, 1.0], [2e4, 1.5e4], total_pressure=1.6e5, total_temperature=500.0
    )
    rebuilt = BoundaryLayerEdge(
        valid.distance, valid.pressure, valid.temperature, valid.velocity
    )
    assert rebuilt.velocity == pytest.approx(valid.velocity)
    with pytest.raises(ValueError, match="isentropic"):
        BoundaryLayerEdge(
            valid.distance,
            valid.pressure,
            valid.temperature,
            valid.velocity * np.array([1.0, 1.01]),
        )
    with pytest.raises(ValueError, match="isentropic"):
        BoundaryLayerEdge(
            valid.distance,
            valid.pressure,
            valid.temperature * np.array([1.0, 1.01]),
            valid.velocity,
        )


def test_invalid_march_options_are_rejected() -> None:
    edge = BoundaryLayerEdge.from_pressure(
        [0.0, 1.0], [2e4, 2e4], total_pressure=1.6e5, total_temperature=500.0
    )
    for regime in (BoundaryLayerRegime.TURBULENT, BoundaryLayerRegime.TRANSITIONAL):
        with pytest.raises(ModelRangeError, match="LAMINAR"):
            marching_boundary_layer(edge, regime=regime)
    with pytest.raises(ValueError, match="prandtl_number"):
        marching_boundary_layer(
            edge, regime=BoundaryLayerRegime.LAMINAR, prandtl_number=0.0
        )
    with pytest.raises(ValueError, match="wall_temperature"):
        marching_boundary_layer(
            edge, regime=BoundaryLayerRegime.LAMINAR, wall_temperature=[300.0]
        )


def test_short_pressure_changes_warn() -> None:
    distance = np.linspace(0.0, 1.0, 2001)
    mach = 2.0 + (1.0 + np.tanh((distance - 0.5) / 0.0002)) / 2.0
    pressure = 2e5 * (1.0 + (GAMMA - 1.0) / 2.0 * mach**2) ** (-GAMMA / (GAMMA - 1.0))
    edge = BoundaryLayerEdge.from_pressure(
        distance, pressure, total_pressure=2e5, total_temperature=450.0
    )
    with pytest.warns(ApplicabilityWarning, match="thin-layer"):
        marching_boundary_layer(edge, regime=BoundaryLayerRegime.LAMINAR)


def test_strong_viscous_interaction_warns() -> None:
    mach = 10.0
    ratio = 1.0 + (GAMMA - 1.0) / 2.0 * mach**2
    edge = BoundaryLayerEdge.from_pressure(
        [0.0, 0.01],
        [100.0, 100.0],
        total_pressure=100.0 * ratio ** (GAMMA / (GAMMA - 1.0)),
        total_temperature=60.0 * ratio,
    )
    with pytest.warns(ApplicabilityWarning, match="viscous interaction"):
        marching_boundary_layer(
            edge, regime=BoundaryLayerRegime.LAMINAR, wall_temperature=300.0
        )


def constant_station(xi: float) -> StationEdge:
    return StationEdge(xi, 1.0, 0.0, 0.0, 300.0, 1.0)


def test_march_reports_an_outer_boundary_inside_the_layer() -> None:
    with pytest.raises(BoundaryLayerConvergenceError, match="outer boundary"):
        march(
            constant_station,
            np.array([1.0]),
            eta=eta_grid(1e-3, 1.04, 0.05, 2.0),
            prandtl_number=1.0,
            viscosity_model=LinearViscosity(),
        )


def test_march_rejects_a_start_away_from_the_leading_edge() -> None:
    with pytest.raises(BoundaryLayerConvergenceError, match="xi = 0"):
        march(
            lambda xi: StationEdge(xi + 1.0, 1.0, 0.0, 0.0, 300.0, 1.0),
            np.array([1.0]),
            eta=ETA,
            prandtl_number=1.0,
            viscosity_model=LinearViscosity(),
        )


def test_march_reports_an_unconverged_start() -> None:
    with pytest.raises(BoundaryLayerConvergenceError, match="leading-edge"):
        march(
            lambda xi: StationEdge(xi, 1.0, float("nan"), 0.0, 300.0, 1.0),
            np.array([1.0]),
            eta=ETA,
            prandtl_number=1.0,
            viscosity_model=LinearViscosity(),
        )


def test_march_reports_a_station_that_cannot_be_resolved() -> None:
    def edge_at(xi: float) -> StationEdge:
        beta = float("nan") if xi > 0.5 else 0.0
        return StationEdge(xi, 1.0, beta, 0.0, 300.0, 1.0)

    with pytest.raises(BoundaryLayerConvergenceError, match="did not converge"):
        march(
            edge_at,
            np.array([1.0]),
            eta=ETA,
            prandtl_number=1.0,
            viscosity_model=LinearViscosity(),
        )
