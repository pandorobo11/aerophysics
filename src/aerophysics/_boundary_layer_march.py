"""Downstream march of the transformed laminar boundary-layer equations.

The unknowns at one station are ``f``, ``u = f'`` and ``g = T/T_e`` on a fixed
``eta`` grid, with ``xi`` as the marching variable:

    (C f'')' + f f'' + beta (g - f'^2)        = 2 xi (f' f'_xi - f'' f_xi)
    (C/Pr g')' + f g' + (gamma-1) M_e^2 C f''^2 = 2 xi (f' g_xi - g' f_xi)

Streamwise derivatives use variable-step BDF2 (BDF1 on the first step). Each
station is solved with a damped Newton iteration; the block-tridiagonal Jacobian
is assembled by three-color finite differences and solved with banded LAPACK.
"""

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from scipy.linalg import solve_banded

from aerophysics._array import FloatArray
from aerophysics.exceptions import BoundaryLayerConvergenceError
from aerophysics.transport import DynamicViscosityModel

_UNKNOWNS = 3
_BAND = 2 * _UNKNOWNS - 1
_COLORS = 3
_NEWTON_ITERATIONS = 20
_NEWTON_TOLERANCE = 1e-10
_MAX_VELOCITY_UPDATE = 0.5
_MAX_RELATIVE_TEMPERATURE_UPDATE = 0.7

_FIRST_STEP_FRACTION = 1e-6
_MINIMUM_STEP_FRACTION = 1e-10
_STEP_GROWTH = 1.25
_MAX_STEP_RATIO = 2.0
_MERGE_WITH_TARGET = 1.3
_RELATIVE_STEP = 0.03
_MAX_LOG_VELOCITY_CHANGE = 0.000625
_MAX_LOG_SHEAR_CHANGE = 0.0125
_MAX_HEAT_TRANSFER_CHANGE = 0.02
_HEAT_TRANSFER_FLOOR = 0.1

_SEPARATION_SHEAR_FRACTION = 0.25
_SEPARATION_FIT_POINTS = 6
_EDGE_PROBE_NODES = 5
_EDGE_TOLERANCE = 1e-4


@dataclass(frozen=True, slots=True)
class StationEdge:
    """Edge parameters at one streamwise position.

    ``wall_ratio`` is ``T_w/T_e``, or None for an adiabatic wall.
    """

    xi: float
    velocity: float
    beta: float
    heating: float
    temperature: float
    wall_ratio: float | None


@dataclass(frozen=True, slots=True)
class MarchSolution:
    """Solution at the requested stations, in streamwise order."""

    eta: FloatArray
    distance: FloatArray
    velocity_ratio: FloatArray
    temperature_ratio: FloatArray
    wall_velocity_slope: FloatArray
    wall_temperature_slope: FloatArray


class SeparationError(Exception):
    """The wall shear collapsed; carries the stations solved upstream of it."""

    def __init__(self, distance: float, solution: MarchSolution) -> None:
        super().__init__(f"separation at {distance:.6g}")
        self.distance = distance
        self.solution = solution


@dataclass(frozen=True, slots=True)
class _Station:
    xi: float
    f: FloatArray
    u: FloatArray
    g: FloatArray
    wall_velocity_slope: float
    wall_temperature_slope: float


def eta_grid(
    first: float, ratio: float, largest_spacing: float, maximum: float
) -> FloatArray:
    """Return a geometrically stretched grid with a capped outer spacing."""
    spacing = first
    nodes = [0.0]
    while nodes[-1] < maximum:
        nodes.append(nodes[-1] + spacing)
        spacing = min(spacing * ratio, largest_spacing)
    return np.asarray(nodes, dtype=np.float64)


class _Discretization:
    """Residual and Jacobian of one station on a fixed ``eta`` grid."""

    def __init__(
        self,
        eta: FloatArray,
        prandtl_number: float,
        viscosity_model: DynamicViscosityModel,
    ) -> None:
        self.eta = eta
        self.prandtl_number = prandtl_number
        self.viscosity_model = viscosity_model
        self.h = np.diff(eta)
        left, right = self.h[:-1], self.h[1:]
        self.width = (left + right) / 2.0
        # Second-order first derivative on the nonuniform three-point stencil.
        self.gradient_weights = np.array(
            [
                -right / (left * (left + right)),
                (right - left) / (left * right),
                left / (right * (left + right)),
            ]
        )
        first, second = eta[1], eta[2]
        self.wall_weights = np.array(
            [
                -(first + second) / (first * second),
                second / (first * (second - first)),
                -first / (second * (second - first)),
            ]
        )
        node, offset, equation, variable = np.indices(
            (len(eta), _COLORS, _UNKNOWNS, _UNKNOWNS)
        )
        neighbor = node + offset - 1
        valid = (neighbor >= 0) & (neighbor < len(eta))
        self.band_indices = (
            node[valid],
            neighbor[valid],
            equation[valid],
            variable[valid],
        )

    def chapman_rubesin(self, g: FloatArray, edge_temperature: float) -> FloatArray:
        """Return ``C = rho mu/(rho_e mu_e)`` for a thin layer of perfect gas."""
        model = self.viscosity_model
        edge = float(model.dynamic_viscosity(edge_temperature))
        wall_normal = np.asarray(
            model.dynamic_viscosity(edge_temperature * g), dtype=np.float64
        )
        return np.asarray(wall_normal / edge / g, dtype=np.float64)

    def derivative(self, value: FloatArray) -> FloatArray:
        weights = self.gradient_weights
        return np.asarray(
            weights[0] * value[..., :-2]
            + weights[1] * value[..., 1:-1]
            + weights[2] * value[..., 2:],
            dtype=np.float64,
        )

    def residual(
        self,
        state: FloatArray,
        edge: StationEdge,
        coefficient: float,
        lagged: FloatArray,
    ) -> FloatArray:
        """Return residuals for ``(..., 3, N)`` states ``(f, u, g)``."""
        f, u, g = np.moveaxis(state, -2, 0)
        h, width = self.h, self.width
        c = self.chapman_rubesin(g, edge.temperature)
        scaled_xi = 2.0 * edge.xi
        convection = f[..., 1:-1] + scaled_xi * (
            coefficient * f[..., 1:-1] + lagged[0, 1:-1]
        )

        def transport(value: FloatArray, diffusivity: FloatArray) -> FloatArray:
            flux = (
                (diffusivity[..., 1:] + diffusivity[..., :-1])
                / 2.0
                * (np.diff(value) / h)
            )
            return np.asarray(np.diff(flux) / width, dtype=np.float64)

        velocity_slope = self.derivative(u)
        inner_u, inner_g = u[..., 1:-1], g[..., 1:-1]
        result = np.empty_like(state)
        result[..., 0, 0] = f[..., 0]
        result[..., 0, 1:] = np.diff(f) - h / 2.0 * (u[..., 1:] + u[..., :-1])
        result[..., 1, 0] = u[..., 0]
        result[..., 1, -1] = u[..., -1] - 1.0
        result[..., 1, 1:-1] = (
            transport(u, c)
            + convection * velocity_slope
            - scaled_xi * inner_u * (coefficient * inner_u + lagged[1, 1:-1])
            + edge.beta * (inner_g - inner_u**2)
        )
        conductivity = c / self.prandtl_number
        result[..., 2, 1:-1] = (
            transport(g, conductivity)
            + convection * self.derivative(g)
            + edge.heating * c[..., 1:-1] * velocity_slope**2
            - scaled_xi * inner_u * (coefficient * inner_g + lagged[2, 1:-1])
        )
        result[..., 2, -1] = g[..., -1] - 1.0
        if edge.wall_ratio is None:
            # Half-cell energy balance with zero wall flux: at the wall the
            # convective terms vanish, leaving only viscous dissipation.
            face = (conductivity[..., 0] + conductivity[..., 1]) / 2.0
            wall_slope = (u[..., 1] - u[..., 0]) / h[0]
            result[..., 2, 0] = (
                face * (g[..., 1] - g[..., 0]) / h[0]
                + h[0] / 2.0 * edge.heating * c[..., 0] * wall_slope**2
            )
        else:
            result[..., 2, 0] = g[..., 0] - edge.wall_ratio
        return result

    def jacobian(
        self,
        state: FloatArray,
        residual: FloatArray,
        edge: StationEdge,
        coefficient: float,
        lagged: FloatArray,
    ) -> FloatArray:
        """Assemble a three-color finite-difference Jacobian in LAPACK bands."""
        step = np.sqrt(np.finfo(np.float64).eps) * np.maximum(np.abs(state), 1.0)
        colors = _COLORS * _UNKNOWNS
        perturbed = np.broadcast_to(state, (colors, *state.shape)).copy()
        variables, nodes = np.indices(state.shape)
        perturbed[_UNKNOWNS * (nodes % _COLORS) + variables, variables, nodes] += step
        difference = self.residual(perturbed, edge, coefficient, lagged) - residual
        row, column, equation, variable = self.band_indices
        bands = np.zeros((2 * _BAND + 1, state.size), dtype=np.float64)
        bands[
            _BAND + _UNKNOWNS * (row - column) + equation - variable,
            _UNKNOWNS * column + variable,
        ] = (
            difference[_UNKNOWNS * (column % _COLORS) + variable, equation, row]
            / step[variable, column]
        )
        return bands

    def station(self, state: FloatArray, xi: float) -> _Station:
        f, u, g = state
        return _Station(
            xi=xi,
            f=f.copy(),
            u=u.copy(),
            g=g.copy(),
            wall_velocity_slope=float(np.dot(self.wall_weights, u[:3])),
            wall_temperature_slope=float(np.dot(self.wall_weights, g[:3])),
        )

    def solve(
        self,
        guess: FloatArray,
        edge: StationEdge,
        coefficient: float,
        lagged: FloatArray,
    ) -> _Station | None:
        """Run the damped Newton iteration; None requests a smaller step."""
        state = guess.copy()
        try:
            for _ in range(_NEWTON_ITERATIONS):
                residual = self.residual(state, edge, coefficient, lagged)
                bands = self.jacobian(state, residual, edge, coefficient, lagged)
                update = (
                    solve_banded(
                        (_BAND, _BAND),
                        bands,
                        -residual.T.ravel(),
                        check_finite=False,
                    )
                    .reshape(-1, _UNKNOWNS)
                    .T
                )
                if not np.all(np.isfinite(update)):
                    return None
                velocity_change = float(np.max(np.abs(update[1])))
                temperature_change = float(np.max(np.abs(update[2] / state[2])))
                damping = min(
                    1.0,
                    _MAX_VELOCITY_UPDATE / max(velocity_change, 1e-30),
                    _MAX_RELATIVE_TEMPERATURE_UPDATE / max(temperature_change, 1e-30),
                )
                state += damping * update
                if max(velocity_change, temperature_change) < _NEWTON_TOLERANCE:
                    return self.station(state, edge.xi)
        except (ValueError, FloatingPointError, np.linalg.LinAlgError):
            # Viscosity models reject nonphysical trial temperatures.
            return None
        return None


def _starting_guess(
    discretization: _Discretization, edge: StationEdge, prandtl_number: float
) -> FloatArray:
    eta = discretization.eta
    u = np.tanh(0.5 * eta)
    recovery = 1.0 + np.sqrt(prandtl_number) * edge.heating / 2.0
    wall = recovery if edge.wall_ratio is None else edge.wall_ratio
    g = wall + (recovery - wall) * u + (1.0 - recovery) * u**2
    f = np.concatenate(([0.0], np.cumsum(discretization.h * (u[1:] + u[:-1]) / 2.0)))
    return np.array([f, u, g], dtype=np.float64)


def _streamwise_stencil(xi: float, history: list[_Station]) -> tuple[float, FloatArray]:
    """Return ``c0`` and lagged terms with ``phi_xi = c0 phi + lagged``."""
    last = history[-1]
    latest = np.array([last.f, last.u, last.g])
    first_gap = xi - last.xi
    if len(history) == 1:
        coefficient = 1.0 / first_gap
        return coefficient, -coefficient * latest
    older = history[-2]
    second_gap = xi - older.xi
    coefficient = 1.0 / first_gap + 1.0 / second_gap
    lagged = -second_gap / (
        first_gap * (second_gap - first_gap)
    ) * latest + first_gap / (second_gap * (second_gap - first_gap)) * np.array(
        [older.f, older.u, older.g]
    )
    return coefficient, lagged


def _build_solution(
    eta: FloatArray, distances: list[float], stations: list[_Station]
) -> MarchSolution:
    count = len(stations)
    velocity = np.empty((count, len(eta)), dtype=np.float64)
    temperature = np.empty((count, len(eta)), dtype=np.float64)
    for index, station in enumerate(stations):
        velocity[index] = station.u
        temperature[index] = station.g
    return MarchSolution(
        eta=eta,
        distance=np.asarray(distances, dtype=np.float64),
        velocity_ratio=velocity,
        temperature_ratio=temperature,
        wall_velocity_slope=np.asarray(
            [station.wall_velocity_slope for station in stations], dtype=np.float64
        ),
        wall_temperature_slope=np.asarray(
            [station.wall_temperature_slope for station in stations], dtype=np.float64
        ),
    )


def _step_is_resolved(previous: _Station, station: _Station, edge: StationEdge) -> bool:
    """Reject a converged step that changes the wall shear or heat flux too much."""
    shear = station.wall_velocity_slope
    if shear <= 0.0:
        return False
    if abs(np.log(shear / previous.wall_velocity_slope)) > _MAX_LOG_SHEAR_CHANGE:
        return False
    if edge.wall_ratio is None:
        return True
    # The heat flux can pass through zero, so a relative change needs a floor: the
    # scale of the wall temperature gradient for a wall-to-recovery difference of
    # the order of the edge total temperature.
    slope, earlier = station.wall_temperature_slope, previous.wall_temperature_slope
    floor = _HEAT_TRANSFER_FLOOR * shear * (1.0 + edge.heating / 2.0)
    return abs(slope - earlier) <= _MAX_HEAT_TRANSFER_CHANGE * max(
        abs(slope), abs(earlier), floor
    )


def _separation_position(positions: list[float], shear: list[float]) -> float:
    """Extrapolate the square-root collapse of the wall shear to zero."""
    position = np.asarray(positions[-_SEPARATION_FIT_POINTS:], dtype=np.float64)
    squared = np.asarray(shear[-_SEPARATION_FIT_POINTS:], dtype=np.float64) ** 2
    intercept = np.polyfit(squared, position, 1)[1]
    return max(float(intercept), positions[-1])


def march(
    edge_at: Callable[[float], StationEdge],
    targets: FloatArray,
    *,
    eta: FloatArray,
    prandtl_number: float,
    viscosity_model: DynamicViscosityModel,
    relative_step: float = _RELATIVE_STEP,
) -> MarchSolution:
    """March from the leading edge at zero through the increasing ``targets``."""
    discretization = _Discretization(eta, prandtl_number, viscosity_model)
    end = float(targets[-1])
    start_edge = edge_at(0.0)
    if start_edge.xi != 0.0:
        raise BoundaryLayerConvergenceError("the march must start at xi = 0")
    start = discretization.solve(
        _starting_guess(discretization, start_edge, prandtl_number),
        start_edge,
        0.0,
        np.zeros((_UNKNOWNS, len(eta)), dtype=np.float64),
    )
    if start is None:
        raise BoundaryLayerConvergenceError(
            "the leading-edge similarity profile did not converge"
        )

    history = [start]
    position = 0.0
    velocity = start_edge.velocity
    step = last_step = _FIRST_STEP_FRACTION * min(float(targets[0]), end)
    minimum_step = _MINIMUM_STEP_FRACTION * end
    accepted_positions = [0.0]
    accepted_shear = [start.wall_velocity_slope]
    output_positions: list[float] = []
    output_stations: list[_Station] = []
    for target in (float(value) for value in targets):
        while position < target:
            trial = min(step, _MAX_STEP_RATIO * last_step)
            if position > 0.0:
                trial = min(trial, relative_step * position)
            remaining = target - position
            if remaining <= _MERGE_WITH_TARGET * trial:
                trial = remaining
            elif remaining < 2.0 * trial:
                trial = remaining / 2.0
            candidate = target if trial == remaining else position + trial
            edge = edge_at(candidate)
            station = None
            if abs(np.log(edge.velocity / velocity)) <= _MAX_LOG_VELOCITY_CHANGE:
                coefficient, lagged = _streamwise_stencil(edge.xi, history)
                guess = np.array(
                    [history[-1].f, history[-1].u, history[-1].g], dtype=np.float64
                )
                if edge.wall_ratio is not None:
                    guess[2, 0] = edge.wall_ratio
                station = discretization.solve(guess, edge, coefficient, lagged)
            if station is not None and not _step_is_resolved(
                history[-1], station, edge
            ):
                station = None
            if station is None:
                step = trial / 2.0
                if step >= minimum_step:
                    continue
                peak = max(accepted_shear)
                if (
                    edge.beta < 0.0
                    and accepted_shear[-1] < _SEPARATION_SHEAR_FRACTION * peak
                ):
                    raise SeparationError(
                        _separation_position(accepted_positions, accepted_shear),
                        _build_solution(eta, output_positions, output_stations),
                    )
                raise BoundaryLayerConvergenceError(
                    f"the march did not converge at distance {candidate:.6g}"
                )
            outer = np.abs(station.u[-_EDGE_PROBE_NODES] - 1.0) + np.abs(
                station.g[-_EDGE_PROBE_NODES] - 1.0
            )
            if outer > _EDGE_TOLERANCE:
                raise BoundaryLayerConvergenceError(
                    "the outer boundary is too close to the boundary layer at "
                    f"distance {candidate:.6g}"
                )
            history = [*history[-1:], station]
            position, velocity, last_step = candidate, edge.velocity, trial
            step = trial * _STEP_GROWTH
            accepted_positions.append(position)
            accepted_shear.append(station.wall_velocity_slope)
        output_positions.append(position)
        output_stations.append(history[-1])
    return _build_solution(eta, output_positions, output_stations)
