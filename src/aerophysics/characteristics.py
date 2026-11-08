"""Rotational characteristics and fitted shocks around sharp convex bodies.

SI units; angles are measured from the positive axial direction in radians.
The gas is calorically perfect, the free stream uniform, and incidence zero.
"""

from dataclasses import dataclass, field, fields
from math import asin, atan, atan2, cos, exp, log, sin, sqrt, tan
from numbers import Real
from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

from aerophysics._array import FloatArray, as_float_array
from aerophysics.exceptions import (
    ModelRangeError,
    NoAttachedShockError,
    ShockConvergenceError,
)
from aerophysics.expansion import prandtl_meyer_expansion
from aerophysics.gas import AIR, PerfectGas
from aerophysics.marching_boundary_layer import BoundaryLayerEdge
from aerophysics.shocks import _taylor_maccoll_rhs, conical_shock, oblique_shock
from aerophysics.transport import AIR_VISCOSITY, DynamicViscosityModel

_MAX_ITERATIONS = 30
_MAX_STATIONS = 100000
_MAX_FAN_STEPS = 10000


class BoundaryLayerInputs(NamedTuple):
    """Isentropic edge and corresponding owned, read-only wall radii in metres."""

    edge: BoundaryLayerEdge
    wall_radius: FloatArray


@dataclass(frozen=True, slots=True)
class AxisymmetricBodyFlowResult:
    """Surface and fitted-shock stations, with owned read-only arrays.

    At a convex vertex the surface array retains only the downstream limit.
    Intermediate fan states are internal, at the same physical vertex. The
    sharp-tip state is the limiting cone surface state. ``gas`` is retained
    so that :meth:`boundary_layer_inputs` also works for non-air perfect gases.
    """

    axial_distance: FloatArray
    radius: FloatArray
    surface_distance: FloatArray
    pressure: FloatArray
    temperature: FloatArray
    mach: FloatArray
    flow_angle: FloatArray
    shock_axial_distance: FloatArray
    shock_radius: FloatArray
    shock_angle: FloatArray
    total_pressure: float
    total_temperature: float
    gas: PerfectGas = AIR
    _corner_states: tuple[tuple[int, float], ...] = field(default=(), repr=False)

    def __post_init__(self) -> None:
        for descriptor in fields(self):
            value = getattr(self, descriptor.name)
            if isinstance(value, np.ndarray):
                owned = np.array(value, dtype=np.float64, copy=True)
                owned.flags.writeable = False
                object.__setattr__(self, descriptor.name, owned)

    def boundary_layer_inputs(
        self,
        *,
        corner_length: float | None = None,
        viscosity_model: DynamicViscosityModel = AIR_VISCOSITY,
    ) -> BoundaryLayerInputs:
        """Return the edge and matching radius, resolving each convex corner.

        Insert the upstream limiting pressure at ``s_c - w`` and a logarithmic
        pressure midpoint at ``s_c - w/2``. Here ``w`` is ``corner_length`` in
        metres, or ``0.001 * s_c`` by default. The midpoint exposes the steep
        pressure gradient to the boundary-layer applicability warning sampler.
        This finite ramp does not resolve viscous interaction at a sharp corner.

        Raises ``ValueError`` for a nonpositive/nonfinite width or when the
        width reaches the preceding surface station. Public surface arrays
        retain their original MOC stations and downstream corner limits.
        """
        if corner_length is not None and (
            isinstance(corner_length, bool)
            or not isinstance(corner_length, Real)
            or not np.isfinite(corner_length)
            or corner_length <= 0
        ):
            raise ValueError("corner_length must be a positive finite scalar in metres")
        distances = self.surface_distance.tolist()
        pressures = self.pressure.tolist()
        # Reverse insertion preserves the stored surface indices at later corners.
        for index, upstream_pressure in reversed(self._corner_states):
            sc = float(self.surface_distance[index])
            width = 0.001 * sc if corner_length is None else float(corner_length)
            before, midpoint = sc - width, sc - width / 2
            if not self.surface_distance[index - 1] < before < midpoint < sc:
                raise ValueError(
                    f"corner_length at s={sc:g} m must be smaller than the "
                    "preceding output interval and resolve two distinct points"
                )
            distances[index:index] = [before, midpoint]
            pressures[index:index] = [
                upstream_pressure,
                exp((log(upstream_pressure) + log(self.pressure[index])) / 2),
            ]
        edge = BoundaryLayerEdge.from_pressure(
            distances,
            pressures,
            total_pressure=self.total_pressure,
            total_temperature=self.total_temperature,
            gas=self.gas,
            viscosity_model=viscosity_model,
        )
        radius = np.interp(edge.distance, self.surface_distance, self.radius)
        radius.flags.writeable = False
        return BoundaryLayerInputs(edge, radius)


@dataclass(frozen=True, slots=True)
class _Point:
    x: float
    r: float
    theta: float
    p: float  # log(p / p_inf)
    p0: float  # log(p0 / p_inf); entropy / R differs by -p0


@dataclass(frozen=True, slots=True)
class _Flow:
    mach: float
    gas: PerfectGas
    axisymmetric: bool = True

    @property
    def gamma(self) -> float:
        return self.gas.heat_capacity_ratio

    @property
    def upstream_p0(self) -> float:
        g = self.gamma
        return g / (g - 1.0) * log(1.0 + (g - 1.0) * self.mach**2 / 2.0)

    def mach_at(self, point: _Point) -> float:
        g = self.gamma
        m2 = 2.0 / (g - 1.0) * (exp((point.p0 - point.p) * (g - 1.0) / g) - 1.0)
        if m2 <= 1.0 or not np.isfinite(m2):
            raise _failure(point, "subsonic or non-finite state")
        return sqrt(m2)

    def coefficients(self, point: _Point, sign: int) -> tuple[float, float, float]:
        m = self.mach_at(point)
        mu = asin(1.0 / m)
        direction = point.theta + sign * mu
        if abs(direction) >= np.pi / 2.0 or point.r <= 0.0:
            raise _failure(point, "characteristic cannot march downstream")
        slope = tan(direction)
        k = sqrt(m * m - 1.0) / (self.gamma * m * m)
        source = (
            -sign * sin(mu) * sin(point.theta) / (point.r * cos(direction))
            if self.axisymmetric
            else 0.0
        )
        return slope, sign * k, source

    def shock(self, x: float, r: float, beta: float) -> _Point:
        """Rankine--Hugoniot state at a fitted shock with uniform upstream."""
        g = self.gamma
        mn2 = (self.mach * sin(beta)) ** 2
        if beta <= asin(1.0 / self.mach) or beta >= np.pi / 2:
            raise ModelRangeError(f"invalid shock angle at x={x:g}, r={r:g}")
        density = (g + 1.0) * mn2 / ((g - 1.0) * mn2 + 2.0)
        pressure = 1.0 + 2.0 * g / (g + 1.0) * (mn2 - 1.0)
        theta = beta - atan(tan(beta) / density)
        loss = g / (g - 1) * log(density) + log(
            (g + 1.0) / (2.0 * g * mn2 - (g - 1.0))
        ) / (g - 1.0)
        return _Point(x, r, theta, log(pressure), self.upstream_p0 + loss)


def _failure(point: _Point, reason: str) -> ModelRangeError:
    return ModelRangeError(f"{reason} at x={point.x:.12g} m, r={point.r:.12g} m")


def _intersection(
    a: _Point, b: _Point, la: float, lb: float, *, check_forward: bool = True
) -> tuple[float, float]:
    if abs(la - lb) < 1e-14:
        raise _failure(a, "parallel or coalescing characteristics")
    x = (b.r - a.r + la * a.x - lb * b.x) / (la - lb)
    r = a.r + la * (x - a.x)
    if check_forward and x < max(a.x, b.x) - 1e-10 * max(abs(x), 1.0):
        raise _failure(_Point(x, r, 0, 0, 0), "crossed characteristics")
    return x, r


def _entropy(a: _Point, b: _Point, target: _Point, *, check_cell: bool = True) -> float:
    """Trace a streamline to the known edge a--b, using midpoint direction."""
    theta = 0.5 * (target.theta + 0.5 * (a.theta + b.theta))
    for _ in range(_MAX_ITERATIONS):
        slope = tan(theta)
        denominator = b.r - a.r - slope * (b.x - a.x)
        if abs(denominator) < 1e-15:
            raise _failure(target, "degenerate streamline interpolation")
        fraction = (target.r - a.r - slope * (target.x - a.x)) / denominator
        corrected = 0.5 * (target.theta + a.theta + fraction * (b.theta - a.theta))
        if abs(corrected - theta) < 2e-12:
            break
        theta = corrected
    else:
        raise _failure(target, "streamline corrector did not converge")
    if check_cell and not -1e-7 <= fraction <= 1.0 + 1e-7:
        raise _failure(target, "streamline left its characteristic cell")
    return a.p0 + fraction * (b.p0 - a.p0)


def _interior(a: _Point, b: _Point, flow: _Flow) -> _Point:
    """Intersect C+ from a and C- from b; iterate trapezoidal corrector."""
    ca, cb = flow.coefficients(a, 1), flow.coefficients(b, -1)
    q = _Point(b.x, b.r, (a.theta + b.theta) / 2, (a.p + b.p) / 2, b.p0)
    for iteration in range(_MAX_ITERATIONS):
        da = ca if iteration == 0 else flow.coefficients(q, 1)
        db = cb if iteration == 0 else flow.coefficients(q, -1)
        la, ka, ha = (0.5 * (v + w) for v, w in zip(ca, da, strict=True))
        lb, kb, hb = (0.5 * (v + w) for v, w in zip(cb, db, strict=True))
        # A first-order predictor can temporarily lie upstream of a parent
        # on the curved exact starting line. Enforce causality and streamline
        # containment on the converged trapezoidal cell, not on trial states.
        x, r = _intersection(a, b, la, lb, check_forward=False)
        va = a.theta + ka * a.p + ha * (x - a.x)
        vb = b.theta + kb * b.p + hb * (x - b.x)
        p = (va - vb) / (ka - kb)
        theta = va - ka * p
        trial = _Point(x, r, theta, p, q.p0)
        new = _Point(x, r, theta, p, _entropy(a, b, trial, check_cell=False))
        if _converged(q, new):
            _intersection(a, b, la, lb)
            _entropy(a, b, new)
            flow.mach_at(new)
            return new
        q = new
    raise _failure(q, "interior corrector did not converge")


def _converged(a: _Point, b: _Point) -> bool:
    return (
        max(
            abs(a.x - b.x) / max(abs(b.x), 1e-10),
            abs(a.r - b.r) / max(abs(b.r), 1e-10),
            abs(a.p - b.p),
            abs(a.p0 - b.p0),
            abs(a.theta - b.theta),
        )
        < 2e-10
    )


def _wall(a: _Point, wall: _Point, theta: float, flow: _Flow) -> _Point:
    """C- from a meets the straight wall; wall-streamline p0 is constant."""
    ca = flow.coefficients(a, -1)
    q = _Point(a.x, a.r, theta, a.p, wall.p0)
    for iteration in range(_MAX_ITERATIONS):
        cb = ca if iteration == 0 else flow.coefficients(q, -1)
        slope, k, h = (0.5 * (v + w) for v, w in zip(ca, cb, strict=True))
        x, r = _intersection(wall, a, tan(theta), slope)
        p = a.p + (a.theta - theta + h * (x - a.x)) / k
        new = _Point(x, r, theta, p, wall.p0)
        if _converged(q, new):
            flow.mach_at(new)
            return new
        q = new
    raise _failure(q, "wall corrector did not converge")


class _ShockAngleTrialError(Exception):
    """A numerical root trial, not a physical state, left the shock interval."""


def _shock_point(
    a: _Point,
    previous: _Point,
    beta: float,
    flow: _Flow,
    *,
    x_limit: float = np.inf,
) -> tuple[_Point, float] | None:
    """Fit the shock; optionally abandon an unsolved downstream extrapolation."""
    ca = flow.coefficients(a, 1)
    continuation_x: float | None = None

    def evaluate(angle: float) -> tuple[float, _Point]:
        nonlocal continuation_x
        if not asin(1.0 / flow.mach) < angle < np.pi / 2:
            raise _ShockAngleTrialError
        state = flow.shock(previous.x, previous.r, angle)
        # Source depends on radius; geometric corrector uses the same trapezoid.
        for _ in range(_MAX_ITERATIONS):
            cb = flow.coefficients(state, 1)
            slope, k, h = (0.5 * (v + w) for v, w in zip(ca, cb, strict=True))
            x, r = _intersection(a, previous, slope, 0.5 * (tan(beta) + tan(angle)))
            new = _Point(x, r, state.theta, state.p, state.p0)
            if _converged(state, new):
                state = new
                break
            state = new
        else:
            raise _failure(state, "shock geometry corrector did not converge")
        if angle == beta:
            # Geometric continuation at the last known shock tangent. Retain
            # successful extrapolated shock output for backward compatibility.
            continuation_x = state.x
        residual = state.theta - a.theta + k * (state.p - a.p) - h * (state.x - a.x)
        return residual, state

    try:
        root = root_scalar(
            lambda angle: evaluate(angle)[0],
            x0=beta,
            x1=beta * (1 + 1e-5),
            xtol=1e-12,
        )
        if not root.converged:
            raise _failure(a, "shock corrector did not converge")
        _, point = evaluate(float(root.root))
    except _ShockAngleTrialError as error:
        if continuation_x is not None and continuation_x > x_limit:
            return None
        raise _failure(
            a, "shock angle corrector left its admissible interval"
        ) from error
    return point, float(root.root)


def _initial_line(
    x: float, angle: float, count: int, flow: _Flow
) -> tuple[list[_Point], float, _Point]:
    """Exact leading C+ from the corner through the conical Taylor--Maccoll field."""
    if flow.axisymmetric:
        shock = conical_shock(flow.mach, angle, flow.gas)
        beta = float(shock.shock_angle)
        surface_pressure = float(shock.surface_pressure_ratio)
    else:
        wedge = oblique_shock(flow.mach, angle, gas=flow.gas)
        beta = float(wedge.shock_angle)
        surface_pressure = float(wedge.static_pressure_ratio)
    post = flow.shock(x, x * tan(beta), beta)
    wall = _Point(x, x * tan(angle), angle, log(surface_pressure), post.p0)
    flow.mach_at(wall)
    if not flow.axisymmetric:
        slope = tan(angle + asin(1 / flow.mach_at(wall)))
        shock_x = (wall.r - slope * x) / (tan(beta) - slope)
        xx = np.linspace(x, shock_x, count + 1)
        return (
            [
                _Point(
                    float(v), wall.r + slope * (float(v) - x), angle, wall.p, wall.p0
                )
                for v in xx
            ],
            beta,
            wall,
        )
    g = flow.gamma
    speed = flow.mach / sqrt(flow.mach**2 + 2.0 / (g - 1.0))
    mn2 = (flow.mach * sin(beta)) ** 2
    density = (g + 1) * mn2 / ((g - 1) * mn2 + 2)
    solution = solve_ivp(
        lambda t, y: _taylor_maccoll_rhs(t, y, flow.gas),
        (beta, angle),
        [speed * cos(beta), -speed * sin(beta) / density],
        dense_output=True,
        method="DOP853",
        rtol=1e-11,
        atol=1e-13,
    )
    if not solution.success or solution.sol is None:
        raise _failure(wall, "Taylor--Maccoll starting line failed")

    def state(polar: float, xx: float) -> _Point:
        assert solution.sol is not None
        vr, vt = solution.sol(polar)
        speed2 = float(vr * vr + vt * vt)
        p = post.p0 + g / (g - 1) * log(1 - speed2)
        theta = polar + atan2(float(vt), float(vr))
        return _Point(xx, xx * tan(polar), theta, p, post.p0)

    # An exact C+ starting characteristic avoids interpolating the cone field
    # into an artificial numerical boundary at the corner's leading Mach wave.
    def trajectory(polar: float, _: FloatArray) -> list[float]:
        slope = flow.coefficients(state(polar, 1.0), 1)[0]
        return [1.0 / (slope * cos(polar) ** 2 - sin(polar) * cos(polar))]

    path = solve_ivp(
        trajectory,
        (angle, beta),
        [log(x)],
        t_eval=np.linspace(angle, beta, count + 1),
        rtol=1e-11,
        atol=1e-13,
        method="DOP853",
    )
    if not path.success:
        raise _failure(wall, "starting characteristic integration failed")
    points = [
        state(float(t), exp(float(v))) for t, v in zip(path.t, path.y[0], strict=True)
    ]
    points[0] = wall
    points[-1] = flow.shock(points[-1].x, points[-1].r, beta)
    return points, beta, wall


def _row(
    wall: _Point,
    parents: list[_Point],
    flow: _Flow,
    beta: float,
    x_limit: float,
    has_shock: bool,
) -> tuple[list[_Point], float, bool]:
    points = [wall]
    for parent in parents:
        points.append(_interior(points[-1], parent, flow))
        if points[-1].x >= x_limit and len(points) >= 3:
            return points, beta, False
    if has_shock:
        fitted = _shock_point(
            points[-1],
            parents[-1],
            beta,
            flow,
            x_limit=x_limit if len(points) >= 3 else np.inf,
        )
        if fitted is None:
            return points, beta, False
        shock, beta = fitted
        points.append(shock)
    return points, beta, has_shock


def _body_inputs(
    axial_distance: ArrayLike, radius: ArrayLike
) -> tuple[FloatArray, FloatArray, FloatArray]:
    x, _ = as_float_array(axial_distance, name="axial_distance")
    r, _ = as_float_array(radius, name="radius")
    if x.ndim != 1 or x.size < 2 or x.shape != r.shape:
        raise ValueError(
            "axial_distance and radius must be equal-length 1D arrays "
            "with at least two points"
        )
    if x[0] != 0 or r[0] != 0 or np.any(np.diff(x) <= 0) or np.any(r[1:] <= 0):
        raise ValueError(
            "body must start at (0, 0), have strictly increasing x "
            "and positive downstream radii"
        )
    angles = np.arctan(np.diff(r) / np.diff(x))
    if np.any(np.diff(angles) > 1e-12):
        i = int(np.argmax(np.diff(angles))) + 1
        raise ModelRangeError(f"concave corner at x={x[i]:g} m, r={r[i]:g} m")
    return x.copy(), r.copy(), angles


def axisymmetric_body_flow(
    axial_distance: ArrayLike,
    radius: ArrayLike,
    *,
    mach: float,
    pressure: float,
    temperature: float,
    gas: PerfectGas = AIR,
    characteristic_count: int = 88,
    fan_angle_step: float = np.pi / 180.0,
) -> AxisymmetricBodyFlowResult:
    """Solve sharp convex axisymmetric external flow with a fitted shock.

    ``characteristic_count`` is the number of intervals across the starting
    line; ``fan_angle_step`` is the maximum expansion-fan increment in radians.
    Refinement doubles the count and halves the angular step. The returned
    wall stations include all vertices, and start at the sharp tip.
    Shock stations can extend past the body's last x: the leading wave from
    a corner may first reach the shock farther downstream. No outflow state
    is imposed. An unsolved downstream extrapolation may instead leave the
    shock trace ending before the body's last x; no terminal shock state is
    guaranteed. All input lengths are in metres, pressure in Pa and
    temperature in K; angles and ``fan_angle_step`` are in radians.

    Raises
    ------
    ValueError
        Invalid geometry, scalar state or resolution.
    ModelRangeError
        Concavity, detached shock, subsonic flow, crossed characteristics or
        failure of a unit-process corrector. The location is in the message.
    """
    return _solve(
        axial_distance,
        radius,
        mach=mach,
        pressure=pressure,
        temperature=temperature,
        gas=gas,
        characteristic_count=characteristic_count,
        fan_angle_step=fan_angle_step,
    )[0]


def _coincident_vertex(
    coordinate: float, vertices: FloatArray, tolerances: FloatArray
) -> int | None:
    """Identify roundoff overlap using the vertex's adjacent input intervals."""
    index = min(int(np.searchsorted(vertices, coordinate)), len(vertices) - 1)
    if index > 0 and abs(coordinate - vertices[index - 1]) < abs(
        coordinate - vertices[index]
    ):
        index -= 1
    return index if abs(coordinate - vertices[index]) <= tolerances[index] else None


def _append_surface(
    surface: list[_Point],
    point: _Point,
    vertices: FloatArray,
    radii: FloatArray,
    tolerances: FloatArray,
) -> None:
    """Give input geometry priority over roundoff-close automatic stations."""
    index = _coincident_vertex(point.x, vertices, tolerances)
    if index is not None:
        point = _Point(
            float(vertices[index]), float(radii[index]), point.theta, point.p, point.p0
        )
    if surface[-1].x == point.x:
        surface[-1] = point
    else:
        surface.append(point)


def _solve(
    axial_distance: ArrayLike,
    radius: ArrayLike,
    *,
    mach: float,
    pressure: float,
    temperature: float,
    gas: PerfectGas = AIR,
    characteristic_count: int = 88,
    fan_angle_step: float = np.pi / 180,
    axisymmetric: bool = True,
    starting_distance: float | None = None,
) -> tuple[AxisymmetricBodyFlowResult, list[_Point]]:
    original_x, original_r, original_angles = _body_inputs(axial_distance, radius)
    intervals = np.diff(original_x)
    # A few rounding units relative to the SMALLER adjacent input interval.
    # Never merge two input vertices, or use the full body length as tolerance.
    vertex_tolerances = np.maximum(
        32
        * np.finfo(float).eps
        * np.minimum(np.r_[intervals[0], intervals], np.r_[intervals, intervals[-1]]),
        # Coordinate subtraction cannot resolve below representation spacing,
        # even when an adjacent input interval is extremely short.
        4 * np.spacing(original_x),
    )
    # Collinear vertices are output locations, not new characteristic fronts.
    # Nearly coincident fronts have poor intersection conditioning near a shock.
    keep = np.concatenate(
        (
            [0],
            np.flatnonzero(np.diff(original_angles) < -1e-12) + 1,
            [len(original_x) - 1],
        )
    )
    x, r = original_x[keep], original_r[keep]
    if starting_distance is not None:
        if not np.isfinite(starting_distance) or not 0 < starting_distance < x[1]:
            raise ValueError(
                "starting_distance must lie strictly inside the initial cone"
            )
        # Verification only: a zero-turn virtual vertex starts a genuine MOC
        # march on the remaining cone, using the same exact C+ initial line.
        r = np.insert(r, 1, starting_distance * r[1] / x[1])
        x = np.insert(x, 1, starting_distance)
    angles = np.arctan(np.diff(r) / np.diff(x))
    if not isinstance(gas, PerfectGas):
        raise TypeError("gas must be a PerfectGas")
    for name, value in (
        ("mach", mach),
        ("pressure", pressure),
        ("temperature", temperature),
        ("fan_angle_step", fan_angle_step),
    ):
        data, scalar = as_float_array(value, name=name)
        if not scalar or float(data) <= 0:
            raise ValueError(f"{name} must be a positive scalar")
    if mach <= 1:
        raise ModelRangeError("free stream must be supersonic at x=0 m, r=0 m")
    if (
        isinstance(characteristic_count, bool)
        or not isinstance(characteristic_count, int)
        or characteristic_count < 4
    ):
        raise ValueError("characteristic_count must be an integer of at least four")
    flow = _Flow(mach, gas, axisymmetric)
    corners = np.flatnonzero(np.diff(angles) < -1e-12)
    cone_end = (
        (int(corners[0]) + 1 if corners.size else len(x) - 1)
        if starting_distance is None
        else 1
    )
    try:
        points, beta, wall = _initial_line(
            float(x[cone_end]), float(angles[0]), characteristic_count, flow
        )
    except (NoAttachedShockError, ShockConvergenceError) as error:
        raise ModelRangeError(f"starting cone at x=0 m, r=0 m: {error}") from error
    cone_samples = [
        float(xx)
        for xx in np.linspace(0, x[cone_end], characteristic_count + 1)
        if _coincident_vertex(float(xx), original_x, vertex_tolerances) is None
    ]
    surface = [
        _Point(
            float(xx),
            float(np.interp(xx, original_x, original_r)),
            wall.theta,
            wall.p,
            wall.p0,
        )
        for xx in np.sort(np.r_[cone_samples, original_x[original_x <= x[cone_end]]])
    ]
    shocks = [(0.0, 0.0, beta), (points[-1].x, points[-1].r, beta)]
    corner_states: list[tuple[int, float]] = []
    has_shock = True
    last_complete = points
    for segment in range(cone_end, len(angles)):
        angle = float(angles[segment])
        turn = wall.theta - angle
        if turn > 1e-12:
            corner_states.append((len(surface) - 1, pressure * exp(wall.p)))
            start = wall
            steps = max(1, int(np.ceil(turn / fan_angle_step)))
            if steps > _MAX_FAN_STEPS:
                raise _failure(wall, "expansion-fan point budget exceeded")
            for delta in np.linspace(0, turn, steps + 1)[1:]:
                try:
                    expansion = prandtl_meyer_expansion(
                        flow.mach_at(start), float(delta), gas
                    )
                except ValueError as error:
                    raise _failure(start, f"unresolvable expansion: {error}") from error
                wall = _Point(
                    start.x,
                    start.r,
                    start.theta - float(delta),
                    start.p + log(float(expansion.static_pressure_ratio)),
                    start.p0,
                )
                points, beta, has_shock = _row(
                    wall, points[1:], flow, beta, float(x[-1]), has_shock
                )
                if has_shock:
                    shocks.append((points[-1].x, points[-1].r, beta))
                    last_complete = points
            surface[-1] = wall
        while wall.x < x[segment + 1]:
            new = _wall(points[1], wall, angle, flow)
            target_x, target_r = float(x[segment + 1]), float(r[segment + 1])
            target_index = int(np.searchsorted(original_x, target_x))
            at_end = abs(new.x - target_x) <= vertex_tolerances[target_index]
            output_end = target_x if at_end else min(new.x, target_x)
            for index in np.flatnonzero(
                (original_x > wall.x) & (original_x < output_end)
            ):
                vertex = float(original_x[index])
                # At either endpoint the inverse problem becomes ill-conditioned.
                # Reuse its limiting state; only the output geometry is snapped.
                if abs(vertex - wall.x) <= vertex_tolerances[index]:
                    at_vertex = wall
                elif abs(vertex - new.x) <= vertex_tolerances[index]:
                    at_vertex = new
                else:
                    at_vertex = _wall_at(points, wall, angle, vertex, flow)
                _append_surface(
                    surface,
                    _Point(
                        vertex, float(original_r[index]), angle, at_vertex.p, wall.p0
                    ),
                    original_x,
                    original_r,
                    vertex_tolerances,
                )
            if at_end:
                # A full wall step reaches the segment endpoint to roundoff.
                wall = _Point(target_x, target_r, angle, new.p, new.p0)
                points, beta, has_shock = _row(
                    wall, points[2:], flow, beta, float(x[-1]), has_shock
                )
            elif new.x >= target_x:
                # Inverse wall point: interpolate incoming C- foot on first edge.
                new = _wall_at(points, wall, angle, target_x, flow)
                wall = _Point(target_x, target_r, angle, new.p, new.p0)
                points, beta, has_shock = _row(
                    wall, points[1:], flow, beta, float(x[-1]), has_shock
                )
            else:
                wall = new
                points, beta, has_shock = _row(
                    wall, points[2:], flow, beta, float(x[-1]), has_shock
                )
            _append_surface(surface, wall, original_x, original_r, vertex_tolerances)
            if has_shock:
                shocks.append((points[-1].x, points[-1].r, beta))
                last_complete = points
            if len(surface) > _MAX_STATIONS:
                raise _failure(wall, "station budget exceeded")
    xx = np.array([q.x for q in surface])
    rr = np.array([q.r for q in surface])
    mm = np.array([flow.mach_at(q) for q in surface])
    total_t = temperature * (1 + (flow.gamma - 1) * mach**2 / 2)
    ss = np.concatenate(([0.0], np.cumsum(np.hypot(np.diff(xx), np.diff(rr)))))
    sh = np.array(shocks)
    return AxisymmetricBodyFlowResult(
        xx,
        rr,
        ss,
        pressure * np.exp([q.p for q in surface]),
        total_t / (1 + (flow.gamma - 1) * mm**2 / 2),
        mm,
        np.array([q.theta for q in surface]),
        sh[:, 0],
        sh[:, 1],
        sh[:, 2],
        pressure * exp(wall.p0),
        total_t,
        gas,
        tuple(corner_states),
    ), last_complete


def _wall_at(
    points: list[_Point], wall: _Point, angle: float, target: float, flow: _Flow
) -> _Point:
    """Inverse wall operation with quadratic interpolation on the C+ edge."""
    a, b = points[:2]
    scale = b.x - a.x
    z = (points[2].x - a.x) / scale if len(points) >= 3 else 0.0

    def foot(f: float) -> _Point:
        weights = (
            ((f - 1) * (f - z) / z, f * (f - z) / (1 - z), f * (f - 1) / (z * (z - 1)))
            if len(points) >= 3
            else (1 - f, f)
        )
        return _Point(
            *(
                sum(
                    w * getattr(p, key)
                    for w, p in zip(weights, points[:3], strict=True)
                )
                for key in ("x", "r", "theta", "p", "p0")
            )
        )

    def residual(f: float) -> float:
        return _wall(foot(f), wall, angle, flow).x - target

    try:
        root = root_scalar(residual, bracket=(1e-12, 1.0), xtol=1e-12)
    except ModelRangeError:
        raise
    except (ValueError, RuntimeError) as error:
        raise _failure(
            wall, f"inverse wall root failed for target x={target:g}: {error}"
        ) from error
    if not root.converged:
        raise _failure(
            wall, f"inverse wall root did not converge for target x={target:g}"
        )
    return _wall(foot(float(root.root)), wall, angle, flow)


__all__ = [
    "AxisymmetricBodyFlowResult",
    "BoundaryLayerInputs",
    "axisymmetric_body_flow",
]
