"""Perfect-gas zero-pressure-gradient laminar similarity boundary-value problem.

The Levy--Lees coordinate z uses (C f'')' + f f'' = 0, so physical
lengths contain sqrt(2). C = mu(Te*g)/mu(Te)/g and g = T/Te. State rows
are (Me, Te [K], Tw/Te or NaN for zero wall heat flux, Pr, gamma).
Coefficients are (cf sqrt(Re), delta99 sqrt(Re)/x, delta* sqrt(Re)/x,
theta sqrt(Re)/x, g_wall), initially normalized with model mu(Te).

Conservative second-order differences, damped Newton, and Richardson
extrapolation use 161/321 nodes, then 321/641 only where needed. The final
fallback is residual-controlled collocation. No physical-state cache or
interpolation table is used. The Blasius profile is only an initial guess.
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_bvp
from scipy.interpolate import CubicSpline, PPoly
from scipy.linalg import solve_banded
from scipy.optimize import brentq

from aerophysics._array import FloatArray
from aerophysics.exceptions import ModelRangeError
from aerophysics.transport import DynamicViscosityModel

_EXTENT = 14.0
_BATCH_SIZE = 256
_NEWTON_ITERATIONS = 24
_GRID_TOLERANCE = 2e-4
_GAUSS_POINTS, _GAUSS_WEIGHTS = leggauss(5)


class _ConvergenceError(Exception):
    """Internal numerical failure; dispatcher attaches original array indices."""


def _viscosity(temperature: FloatArray, model: DynamicViscosityModel) -> FloatArray:
    try:
        viscosity = np.asarray(model.dynamic_viscosity(temperature), dtype=np.float64)
    except (ValueError, FloatingPointError) as error:
        raise _ConvergenceError(f"viscosity evaluation: {error}") from error
    if np.any(~np.isfinite(viscosity)) or np.any(viscosity <= 0):
        raise _ConvergenceError("nonpositive or nonfinite viscosity")
    return viscosity


def _properties(
    g: FloatArray, te: FloatArray, model: DynamicViscosityModel
) -> tuple[FloatArray, FloatArray]:
    """Return C and dC/dg without requiring an extended viscosity protocol."""
    if np.any(~np.isfinite(g)) or np.any(g <= 0):
        raise _ConvergenceError("nonpositive or nonfinite temperature")
    temperature = te * g
    mu_e = _viscosity(te, model)
    mu = _viscosity(temperature, model)
    c = mu / mu_e / g
    derivative = getattr(model, "dynamic_viscosity_derivative", None)
    if callable(derivative):
        dc = np.asarray(derivative(temperature)) * te / mu_e / g - c / g
    else:
        epsilon = 1e-5
        plus = _viscosity(temperature * (1 + epsilon), model)
        minus = _viscosity(temperature * (1 - epsilon), model)
        dc = ((plus - minus) / (2 * epsilon) / mu_e / g - c) / g
    if not np.all(np.isfinite(c)) or not np.all(np.isfinite(dc)):
        raise _ConvergenceError("nonfinite viscosity or derivative")
    return c, dc


@lru_cache(maxsize=1)
def _blasius_guess() -> PPoly:
    """One universal numerical initial profile, never a cached physical state."""
    z = np.linspace(0, _EXTENT, 201)
    result = solve_bvp(
        lambda z, y: np.array([y[1], y[2], -y[0] * y[2]]),
        lambda a, b: np.array([a[0], a[1], b[1] - 1]),
        z,
        np.array([np.log(np.cosh(z)), np.tanh(z), 0.47 * np.exp(-z)]),
        tol=1e-10,
        max_nodes=10000,
    )
    if not result.success:
        raise _ConvergenceError("Blasius initialization failed")
    return result.sol


def _initial(
    z: FloatArray, states: FloatArray, model: DynamicViscosityModel
) -> FloatArray:
    mach, te, wall, pr, gamma = states.T
    recovery = 1 + np.sqrt(pr) * (gamma - 1) * mach**2 / 2
    wall = np.where(np.isnan(wall), recovery, wall)
    c, _ = _properties(0.22 * recovery + 0.28 + 0.5 * wall, te, model)
    coordinate = z[:, None] / np.sqrt(c)
    profile = _blasius_guess()(np.minimum(coordinate, _EXTENT))
    u = profile[1]
    f = np.sqrt(c) * profile[0] + np.maximum(z[:, None] - _EXTENT * np.sqrt(c), 0)
    g = wall + (recovery - wall) * u + (1 - recovery) * u**2
    return np.stack((f, u, g), axis=-1)


def _assemble(
    y: FloatArray, z: FloatArray, states: FloatArray, model: DynamicViscosityModel
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    """Residual and block tridiagonal Jacobian for unknowns (f, u=f', g)."""
    nodes, batch, _ = y.shape
    te, wall, pr = states[:, 1], states[:, 2], states[:, 3]
    heating = (states[:, 4] - 1) * states[:, 0] ** 2
    adiabatic = np.isnan(wall)
    f, u, g = y.transpose(2, 0, 1)
    h = np.diff(z)[:, None]
    cm, dc = _properties((g[:-1] + g[1:]) / 2, te, model)
    du, dg = np.diff(u, axis=0) / h, np.diff(g, axis=0) / h
    flux_u, flux_g = cm * du, cm * dg / pr
    lower = np.zeros((nodes, batch, 3, 3))
    diagonal, upper = lower.copy(), lower.copy()
    residual = np.zeros_like(y)
    residual[1:, :, 0] = np.diff(f, axis=0) - h * (u[:-1] + u[1:]) / 2
    lower[1:, :, 0, 0], diagonal[1:, :, 0, 0] = -1, 1
    lower[1:, :, 0, 1] = -h / 2
    diagonal[1:, :, 0, 1] = -h / 2
    width = (h[:-1] + h[1:]) / 2
    uc = (u[2:] - u[:-2]) / (2 * width)
    gc = (g[2:] - g[:-2]) / (2 * width)
    cn, dcn = _properties(g[1:-1], te, model)
    residual[1:-1, :, 1] = np.diff(flux_u, axis=0) / width + f[1:-1] * uc
    residual[1:-1, :, 2] = (
        np.diff(flux_g, axis=0) / width + f[1:-1] * gc + heating * cn * uc**2
    )
    lower[1:-1, :, 1, 1] = cm[:-1] / h[:-1] / width - f[1:-1] / (2 * width)
    upper[1:-1, :, 1, 1] = cm[1:] / h[1:] / width + f[1:-1] / (2 * width)
    diagonal[1:-1, :, 1, 1] = -(cm[:-1] / h[:-1] + cm[1:] / h[1:]) / width
    lower[1:-1, :, 1, 2] = -dc[:-1] * du[:-1] / (2 * width)
    upper[1:-1, :, 1, 2] = dc[1:] * du[1:] / (2 * width)
    diagonal[1:-1, :, 1, 2] = (dc[1:] * du[1:] - dc[:-1] * du[:-1]) / (2 * width)
    diagonal[1:-1, :, 1, 0] = uc
    lower[1:-1, :, 2, 2] = (
        cm[:-1] / pr / h[:-1] / width
        - dc[:-1] * dg[:-1] / (2 * pr * width)
        - f[1:-1] / (2 * width)
    )
    upper[1:-1, :, 2, 2] = (
        cm[1:] / pr / h[1:] / width
        + dc[1:] * dg[1:] / (2 * pr * width)
        + f[1:-1] / (2 * width)
    )
    diagonal[1:-1, :, 2, 2] = (
        -(cm[:-1] / h[:-1] + cm[1:] / h[1:]) / pr / width
        + (dc[1:] * dg[1:] - dc[:-1] * dg[:-1]) / (2 * pr * width)
        + heating * dcn * uc**2
    )
    diagonal[1:-1, :, 2, 0] = gc
    lower[1:-1, :, 2, 1] = -heating * cn * uc / width
    upper[1:-1, :, 2, 1] = heating * cn * uc / width
    residual[0, :, 0], residual[0, :, 1] = f[0], u[0]
    diagonal[0, :, 0, 0], diagonal[0, :, 1, 1] = 1, 1
    # Integrate energy over the wall half-cell with exactly zero boundary flux.
    residual[0, :, 2] = np.where(
        adiabatic, flux_g[0] + h[0] / 2 * heating * cm[0] * du[0] ** 2, g[0] - wall
    )
    diagonal[0, :, 2, 2] = np.where(
        adiabatic,
        (-cm[0] / h[0] + dc[0] * dg[0] / 2) / pr
        + h[0] / 4 * heating * dc[0] * du[0] ** 2,
        1,
    )
    upper[0, :, 2, 2] = np.where(
        adiabatic,
        (cm[0] / h[0] + dc[0] * dg[0] / 2) / pr
        + h[0] / 4 * heating * dc[0] * du[0] ** 2,
        0,
    )
    diagonal[0, :, 2, 1] = np.where(adiabatic, -heating * cm[0] * du[0], 0)
    upper[0, :, 2, 1] = np.where(adiabatic, heating * cm[0] * du[0], 0)
    residual[-1, :, 1], residual[-1, :, 2] = u[-1] - 1, g[-1] - 1
    diagonal[-1, :, 1, 1], diagonal[-1, :, 2, 2] = 1, 1
    return lower, diagonal, upper, residual


def _block_solve(
    lower: FloatArray, diagonal: FloatArray, upper: FloatArray, residual: FloatArray
) -> FloatArray:
    """Banded LAPACK for a singleton; vectorized block Thomas for batches."""
    nodes, batch = diagonal.shape[:2]
    if batch == 1:
        bands = np.zeros((11, 3 * nodes))
        for row in range(3):
            for col in range(3):
                columns = np.arange(nodes) * 3 + col
                bands[5 + row - col, columns] = diagonal[:, 0, row, col]
                bands[8 + row - col, columns[:-1]] = lower[1:, 0, row, col]
                bands[2 + row - col, columns[1:]] = upper[:-1, 0, row, col]
        return solve_banded((5, 5), bands, residual[:, 0].ravel()).reshape(nodes, 1, 3)
    reduced_upper, answer = np.empty_like(upper), np.empty_like(residual)
    for node in range(nodes):
        block = diagonal[node]
        rhs = residual[node]
        if node:
            block = block - lower[node] @ reduced_upper[node - 1]
            rhs = rhs - (lower[node] @ answer[node - 1, ..., None])[..., 0]
        solved = np.linalg.solve(
            block, np.concatenate((upper[node], rhs[..., None]), axis=-1)
        )
        reduced_upper[node], answer[node] = solved[..., :3], solved[..., 3]
    for node in range(nodes - 2, -1, -1):
        answer[node] -= (reduced_upper[node] @ answer[node + 1, ..., None])[..., 0]
    return answer


def _pair(spline: PPoly, positions: FloatArray) -> FloatArray:
    """Evaluate a batched spline at one position per state, without outer products."""
    cells = np.clip(np.searchsorted(spline.x, positions) - 1, 0, len(spline.x) - 2)
    coefficients = spline.c[:, cells, np.arange(len(positions))]
    offset = positions - spline.x[cells]
    result = coefficients[0].copy()
    for coefficient in coefficients[1:]:
        result = result * offset + coefficient
    return np.asarray(result, dtype=np.float64)


def _quantities(
    z: FloatArray, y: FloatArray, states: FloatArray, model: DynamicViscosityModel
) -> FloatArray:
    u, g = y[:, :, 1], y[:, :, 2]
    velocity, temperature = CubicSpline(z, u, axis=0), CubicSpline(z, g, axis=0)
    indices = np.argmax(u >= 0.99, axis=0)
    if np.any(indices == 0):
        raise _ConvergenceError("velocity does not reach 0.99 within domain")
    low, high = z[indices - 1], z[indices]
    root = (low + high) / 2
    derivative = velocity.derivative()
    for _ in range(5):
        root = np.clip(
            root - (_pair(velocity, root) - 0.99) / _pair(derivative, root), low, high
        )
    primitive = temperature.antiderivative()
    delta99 = _pair(primitive, root) - primitive(0)
    displacement = temperature.integrate(0, z[-1]) - velocity.integrate(0, z[-1])
    positions = (z[:-1, None] + z[1:, None]) / 2 + np.diff(z)[
        :, None
    ] / 2 * _GAUSS_POINTS
    uv = velocity(positions)
    momentum = np.sum(
        uv * (1 - uv) * np.diff(z)[:, None, None] / 2 * _GAUSS_WEIGHTS[None, :, None],
        axis=(0, 1),
    )
    cw, _ = _properties(g[0], states[:, 1], model)
    coefficients = np.sqrt(2) * np.stack(
        (cw * velocity(0, 1), delta99, displacement, momentum), axis=1
    )
    return np.column_stack((coefficients, g[0]))


@dataclass(frozen=True)
class _GridSolution:
    z: FloatArray
    y: FloatArray
    quantities: FloatArray
    valid: np.ndarray


def _finite_difference(
    states: FloatArray,
    nodes: int,
    model: DynamicViscosityModel,
    *,
    previous: _GridSolution | None = None,
    max_iterations: int = _NEWTON_ITERATIONS,
) -> _GridSolution:
    z = _EXTENT * np.linspace(0, 1, nodes) ** 1.6
    y = (
        _initial(z, states, model)
        if previous is None
        else CubicSpline(previous.z, previous.y, axis=0)(z)
    )
    scale = np.ones((len(states), 3))
    scale[:, 2] = np.maximum(1, np.max(y[:, :, 2], axis=0))
    converged = np.zeros(len(states), dtype=np.bool_)
    for _ in range(max_iterations):
        lower, diagonal, upper, residual = _assemble(y, z, states, model)
        step = _block_solve(lower, diagonal, upper, -residual)
        converged = np.max(np.abs(step) / scale, axis=(0, 2)) < 2e-9
        if np.all(converged):
            break
        damping = np.minimum(
            1,
            0.7 / np.maximum(np.max(np.abs(step[:, :, 2]) / y[:, :, 2], axis=0), 1e-30),
        )
        damping = np.minimum(
            damping, 0.5 / np.maximum(np.max(np.abs(step[:, :, 1]), axis=0), 1e-30)
        )
        y += damping[None, :, None] * step
    residual = _assemble(y, z, states, model)[3]
    valid = converged & (np.max(np.abs(residual) / scale, axis=(0, 2)) < 2e-6)
    # Values in the last tenth of the domain must already reach the edge state.
    tail = y[z >= 0.9 * _EXTENT]
    valid &= np.max(np.abs(tail[:, :, 1] - 1), axis=0) < 1e-8
    valid &= np.max(np.abs(tail[:, :, 2] - 1), axis=0) / scale[:, 2] < 1e-8
    valid &= np.all(np.isfinite(y), axis=(0, 2)) & np.all(y[:, :, 2] > 0, axis=0)
    return _GridSolution(z, y, _quantities(z, y, states, model), valid)


def _relative_change(first: FloatArray, second: FloatArray) -> FloatArray:
    denominator = np.abs(second).copy()
    # A signed displacement integral crosses zero on cold walls.
    denominator[:, 2] = np.maximum(denominator[:, 2], 0.05)
    return np.max(np.abs(first - second) / denominator, axis=1)


def _consistent(quantities: FloatArray) -> np.ndarray:
    return (
        np.all(np.isfinite(quantities), axis=1)
        & np.all(quantities[:, [0, 1, 3, 4]] > 0, axis=1)
        & (np.abs(quantities[:, 0] / quantities[:, 3] - 1) < _GRID_TOLERANCE)
    )


def _richardson(
    states: FloatArray,
    nodes: int,
    model: DynamicViscosityModel,
    *,
    max_iterations: int,
    check_correction: bool = True,
) -> tuple[FloatArray, np.ndarray]:
    coarse = _finite_difference(states, nodes, model, max_iterations=max_iterations)
    fine = _finite_difference(
        states, 2 * nodes - 1, model, previous=coarse, max_iterations=max_iterations
    )
    extrapolated = (4 * fine.quantities - coarse.quantities) / 3
    valid = coarse.valid & fine.valid & _consistent(extrapolated)
    if check_correction:
        valid &= _relative_change(fine.quantities, extrapolated) < _GRID_TOLERANCE
    return extrapolated, valid


def _collocation(state: FloatArray, model: DynamicViscosityModel) -> FloatArray:
    """Independent flux formulation, checked again on a larger far-field domain."""
    mach, te, wall, pr, gamma = state
    heating = (gamma - 1) * mach**2

    def rhs(z: FloatArray, y: FloatArray) -> FloatArray:
        f, u, stress, g, heat = y
        # A collocation Newton trial may leave the positive-temperature domain.
        # This extension is only for iterations; accepted profiles are checked below.
        c, _ = _properties(np.maximum(g, 1e-12), np.asarray(te), model)
        up, gp = stress / c, pr * heat / c
        return np.array([u, up, -f * up, gp, -f * gp - heating * stress * up])

    def boundary(a: FloatArray, b: FloatArray) -> FloatArray:
        return np.array(
            [a[0], a[1], a[4] if np.isnan(wall) else a[3] - wall, b[1] - 1, b[3] - 1]
        )

    previous: PPoly | None = None
    previous_quantities: FloatArray | None = None
    for extent in (14.0, 21.0, 31.5, 47.25, 70.875):
        z = extent * np.linspace(0, 1, 321) ** 1.6
        if previous is None:
            initial = _initial(z, state[None, :], model)[:, 0].T
            c, _ = _properties(initial[2], np.asarray(te), model)
            derivatives = CubicSpline(z, initial, axis=1)(z, 1)
            guess = np.array(
                [
                    initial[0],
                    initial[1],
                    c * derivatives[1],
                    initial[2],
                    c / pr * derivatives[2],
                ]
            )
        else:
            guess = previous(np.minimum(z, previous.x[-1]))
            beyond = z > previous.x[-1]
            guess[0, beyond] += z[beyond] - previous.x[-1]
        solved = solve_bvp(rhs, boundary, z, guess, tol=1e-8, max_nodes=30000)
        if not solved.success or np.max(solved.rms_residuals) > 1.01e-8:
            raise _ConvergenceError(f"collocation: {solved.message}")
        probe = solved.sol((solved.x[:-1] + solved.x[1:]) / 2)
        if np.min(solved.y[3]) <= 0 or np.min(probe[3]) <= 0:
            raise _ConvergenceError("nonpositive collocation temperature")
        _properties(probe[3], np.asarray(te), model)
        if np.max(np.abs(boundary(solved.y[:, 0], solved.y[:, -1]))) > 1e-8:
            raise _ConvergenceError("collocation boundary residual")
        # Integrate the actual collocation interpolant, retaining its wall stress.
        positions = (solved.x[:-1, None] + solved.x[1:, None]) / 2 + np.diff(solved.x)[
            :, None
        ] / 2 * _GAUSS_POINTS
        values = solved.sol(positions.ravel()).reshape(5, -1, 5)
        weights = np.diff(solved.x)[:, None] / 2 * _GAUSS_WEIGHTS
        z99 = brentq(
            lambda z, profile=solved.sol: float(profile(z)[1]) - 0.99, 0, extent
        )
        knots = np.r_[solved.x[solved.x < z99], z99]
        cut = (knots[:-1, None] + knots[1:, None]) / 2 + np.diff(knots)[
            :, None
        ] / 2 * _GAUSS_POINTS
        g99 = solved.sol(cut.ravel())[3].reshape(-1, 5)
        q = np.sqrt(2) * np.array(
            [
                solved.y[2, 0],
                np.sum(g99 * np.diff(knots)[:, None] / 2 * _GAUSS_WEIGHTS),
                np.sum((values[3] - values[1]) * weights),
                np.sum(values[1] * (1 - values[1]) * weights),
            ]
        )
        quantities = np.r_[q, solved.y[3, 0]]
        tail = solved.sol(np.linspace(0.9 * extent, extent, 10))
        tail_ok = np.max(np.abs(tail[[1, 3]] - 1)) / max(1, solved.y[3].max()) < 1e-9
        flux_ok = (
            np.max(np.abs(tail[[2, 4]]))
            / max(1, abs(solved.y[2, 0]), abs(solved.y[4, 0]))
            < 1e-9
        )
        if (
            previous_quantities is not None
            and tail_ok
            and flux_ok
            and _consistent(quantities[None])[0]
        ):
            if _relative_change(previous_quantities[None], quantities[None])[0] < 1e-6:
                return np.asarray(quantities, dtype=np.float64)
        previous, previous_quantities = solved.sol, quantities
    raise _ConvergenceError("far-field domain or momentum balance did not converge")


def similarity_coefficients(
    states: FloatArray,
    model: DynamicViscosityModel,
    *,
    _max_iterations: int = _NEWTON_ITERATIONS,
    _allow_fallback: bool = True,
    _input_shape: tuple[int, ...] | None = None,
) -> FloatArray:
    """Solve exact unique state rows and return five coefficients per input row.

    Rows are Me, Te, Tw/Te (NaN means adiabatic), Pr, gamma. Errors retain
    original row indices. Callers supplying a different mu_e must multiply
    the first four outputs by sqrt(mu_model(Te)/mu_e). Displacement is signed.
    """
    keys = np.array(states, dtype=np.float64, copy=True)
    if keys.ndim != 2 or keys.shape[1] != 5:
        raise ValueError("similarity states must have shape (n, 5)")
    invalid = (
        ~np.all(np.isfinite(keys[:, [0, 1, 3, 4]]), axis=1)
        | (keys[:, 0] < 0)
        | (keys[:, 1] <= 0)
        | (keys[:, 3] <= 0)
        | (keys[:, 4] <= 1)
        | (~np.isnan(keys[:, 2]) & (~np.isfinite(keys[:, 2]) | (keys[:, 2] <= 0)))
    )
    if np.any(invalid):
        index = int(np.flatnonzero(invalid)[0])
        raise ModelRangeError(
            f"laminar similarity at array index {index}: "
            f"invalid state {keys[index].tolist()}"
        )
    keys[:, 2] = np.where(np.isnan(keys[:, 2]), -1, keys[:, 2])
    unique, inverse = np.unique(keys, axis=0, return_inverse=True)
    unique[:, 2] = np.where(unique[:, 2] == -1, np.nan, unique[:, 2])
    answer = np.empty_like(unique)

    def solve_batch(indices: np.ndarray) -> None:
        batch = unique[indices]
        try:
            quantities, valid = _richardson(
                batch, 161, model, max_iterations=_max_iterations
            )
        except (_ConvergenceError, np.linalg.LinAlgError):
            # Isolate numerical/model failures; unaffected states retain the FD route.
            if len(indices) > 1:
                midpoint = len(indices) // 2
                solve_batch(indices[:midpoint])
                solve_batch(indices[midpoint:])
                return
            quantities, valid = np.full_like(batch, np.nan), np.zeros(1, dtype=np.bool_)
        answer[indices] = quantities
        suspect = np.flatnonzero(~valid)
        if len(suspect):
            try:
                refined, refined_valid = _richardson(
                    batch[suspect],
                    321,
                    model,
                    max_iterations=_max_iterations,
                    check_correction=False,
                )
                # Refinement compares extrapolants, not the larger raw-grid correction.
                refined_valid &= _consistent(refined) & (
                    _relative_change(quantities[suspect], refined) < _GRID_TOLERANCE
                )
                answer[indices[suspect]] = refined
                suspect = suspect[~refined_valid]
            except (_ConvergenceError, np.linalg.LinAlgError):
                pass
        for position in suspect:
            index = indices[position]
            try:
                if not _allow_fallback:
                    raise _ConvergenceError(
                        "finite-difference convergence failed; fallback disabled"
                    )
                answer[index] = _collocation(unique[index], model)
            except (_ConvergenceError, np.linalg.LinAlgError) as error:
                rows = np.flatnonzero(inverse == index)
                original = (
                    rows.tolist()
                    if _input_shape is None
                    else [
                        tuple(
                            int(i)
                            for i in np.unravel_index(
                                int(row) % int(np.prod(_input_shape)), _input_shape
                            )
                        )
                        for row in rows
                    ]
                )
                raise ModelRangeError(
                    f"laminar similarity at array indices {original}: "
                    f"(Me, Te, Tw/Te, Pr, gamma)={unique[index].tolist()}: {error}"
                ) from error

    for start in range(0, len(unique), _BATCH_SIZE):
        solve_batch(np.arange(start, min(start + _BATCH_SIZE, len(unique))))
    # The independent integral was checked above. Enforce the exact ZPG identity
    # in the returned coefficients instead of exposing integration roundoff.
    answer[:, 3] = answer[:, 0]
    return answer[inverse]
