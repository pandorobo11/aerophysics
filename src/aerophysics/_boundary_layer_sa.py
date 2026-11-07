"""Spalart--Allmaras equations on the downstream boundary-layer march.

The model is SA-noft2 in its standard form, in which density does not appear
(Allmaras, Johnson and Spalart 2012, section 2.2). With ``n = nu_tilde/nu_e``,
``Y = int g d(eta)``, ``chi = n/(C g^2)``, and ``eps = chi fv1``, one station solves

    (C(1+eps) f'')' + f f'' + beta (g - f'^2)   = 2 xi (f' f'_xi - f'' f_xi)
    (C/Pr (1 + eps Pr/Prt) g')' + f g'
        + (gamma-1) M_e^2 C (1+eps) f''^2       = 2 xi (f' g_xi - g' f_xi)
    ((C g + n/g) n')'/(sigma g) + cb2 n'^2/(sigma g^2) + f n'
        + cb1 St n - cw1 fw n^2/Y^2 - Lambda f' n = 2 xi (f' n_xi - n' f_xi)

with ``St = R |f''|/g + n fv2/(kappa^2 Y^2)``, ``R = sqrt(2 xi)/(r_w^j mu_e)``,
and ``Lambda = 2 xi d(ln nu_e)/d(xi)``. Here j=0 is planar and j=1 is thin-layer
axisymmetric. The dimensional ``nu_tilde`` is conserved along
the edge streamline, so ``n`` at the outer boundary varies as ``1/nu_e``.

Newton stalls where the positivity limit binds near the outer turbulent front. A
relaxed Picard iteration with the same discrete operators then takes over.
"""

import numpy as np
from scipy.linalg import solve_banded

from aerophysics._array import FloatArray
from aerophysics._boundary_layer_march import Discretization, Station, StationEdge
from aerophysics.transport import DynamicViscosityModel

_CB1 = 0.1355
_CB2 = 0.622
_SIGMA = 2.0 / 3.0
_KAPPA = 0.41
_CV1 = 7.1
_CW2 = 0.3
_CW3 = 2.0
_CW1 = _CB1 / _KAPPA**2 + (1.0 + _CB2) / _SIGMA
# Smooth limiter of the modified vorticity (Allmaras, Johnson and Spalart 2012).
_CV2 = 0.7
_CV3 = 0.9

_MAX_R = 10.0
_MINIMUM_VORTICITY = 1e-12
_POSITIVITY_FRACTION = 0.9
_MINIMUM_DAMPING = 1e-4
_MINIMUM_TEMPERATURE_RATIO = 1e-3

_PICARD_ITERATIONS = 4000
_PICARD_TOLERANCE = 1e-9
_PICARD_UNDER_RELAXATION = 0.5


def _fv1(chi: FloatArray) -> FloatArray:
    return np.asarray(chi**3 / (chi**3 + _CV1**3), dtype=np.float64)


def _source_terms(
    vorticity: FloatArray,
    n: FloatArray,
    c: FloatArray,
    g: FloatArray,
    distance: FloatArray,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Return the modified vorticity, its extra part, and ``fw``."""
    chi = n / (c * g**2)
    fv2 = 1.0 - chi / (1.0 + chi * _fv1(chi))
    extra = n * fv2 / (_KAPPA**2 * distance**2)
    limited = vorticity + vorticity * (_CV2**2 * vorticity + _CV3 * extra) / np.maximum(
        (_CV3 - 2.0 * _CV2) * vorticity - extra, 1e-30
    )
    modified = np.maximum(
        np.where(extra >= -_CV2 * vorticity, vorticity + extra, limited),
        _MINIMUM_VORTICITY,
    )
    r = np.minimum(n / np.maximum(modified * _KAPPA**2 * distance**2, 1e-30), _MAX_R)
    r6 = r + _CW2 * (r**6 - r)
    fw = r6 * ((1.0 + _CW3**6) / (r6**6 + _CW3**6)) ** (1.0 / 6.0)
    return modified, extra, np.asarray(fw, dtype=np.float64)


class TurbulentDiscretization(Discretization):
    """One station of ``(f, u, g, n, Y)`` with a Picard fallback."""

    unknowns = 5
    newton_iterations = 50
    newton_tolerance = 1e-9

    def __init__(
        self,
        eta: FloatArray,
        prandtl_number: float,
        turbulent_prandtl_number: float,
        viscosity_model: DynamicViscosityModel,
    ) -> None:
        super().__init__(eta, prandtl_number, viscosity_model)
        self.turbulent_prandtl_number = turbulent_prandtl_number
        self.picard_fallbacks = 0

    def residual(
        self,
        state: FloatArray,
        edge: StationEdge,
        coefficient: float,
        lagged: FloatArray,
    ) -> FloatArray:
        """Return residuals for ``(..., 5, N)`` states ``(f, u, g, n, Y)``."""
        f, u, g, n, stretched = np.moveaxis(state, -2, 0)
        h, width = self.h, self.width
        c = self.chapman_rubesin(g, edge.temperature)
        chi = n / (c * g**2)
        eddy = chi * _fv1(chi)
        effective = c * (1.0 + eddy)
        conductivity = (
            c
            / self.prandtl_number
            * (1.0 + eddy * self.prandtl_number / self.turbulent_prandtl_number)
        )
        scaled_xi = 2.0 * edge.xi
        convection = f[..., 1:-1] + scaled_xi * (
            coefficient * f[..., 1:-1] + lagged[0, 1:-1]
        )
        inner_u, inner_g = u[..., 1:-1], g[..., 1:-1]
        inner_n = n[..., 1:-1]

        def centered(value: FloatArray) -> FloatArray:
            return np.asarray(
                (value[..., 2:] - value[..., :-2]) / (2.0 * width), dtype=np.float64
            )

        def transport(
            value: FloatArray, diffusivity: FloatArray, history: FloatArray
        ) -> FloatArray:
            flux = (
                (diffusivity[..., 1:] + diffusivity[..., :-1])
                / 2.0
                * (np.diff(value) / h)
            )
            return np.asarray(
                np.diff(flux) / width
                + convection * centered(value)
                - scaled_xi * inner_u * (coefficient * value[..., 1:-1] + history),
                dtype=np.float64,
            )

        velocity_slope = self.derivative(u)
        positive_n = np.maximum(inner_n, 0.0)
        inner_c = c[..., 1:-1]
        distance = stretched[..., 1:-1]
        modified, _, fw = _source_terms(
            edge.reynolds_root * np.abs(velocity_slope) / inner_g,
            positive_n,
            inner_c,
            inner_g,
            distance,
        )

        result = np.empty_like(state)
        result[..., 0, 0] = f[..., 0]
        result[..., 0, 1:] = np.diff(f) - h / 2.0 * (u[..., 1:] + u[..., :-1])
        result[..., 4, 0] = stretched[..., 0]
        result[..., 4, 1:] = np.diff(stretched) - h / 2.0 * (g[..., 1:] + g[..., :-1])

        result[..., 1, 0] = u[..., 0]
        result[..., 1, -1] = u[..., -1] - 1.0
        result[..., 1, 1:-1] = transport(u, effective, lagged[1, 1:-1]) + edge.beta * (
            inner_g - inner_u**2
        )

        result[..., 2, -1] = g[..., -1] - 1.0
        result[..., 2, 1:-1] = (
            transport(g, conductivity, lagged[2, 1:-1])
            + edge.heating * effective[..., 1:-1] * velocity_slope**2
        )
        if edge.wall_ratio is None:
            # Half-cell energy balance with zero wall flux, as in the laminar row.
            face = (conductivity[..., 0] + conductivity[..., 1]) / 2.0
            wall_slope = (u[..., 1] - u[..., 0]) / h[0]
            result[..., 2, 0] = (
                face * (g[..., 1] - g[..., 0]) / h[0]
                + h[0] / 2.0 * edge.heating * effective[..., 0] * wall_slope**2
            )
        else:
            result[..., 2, 0] = g[..., 0] - edge.wall_ratio

        result[..., 3, 0] = n[..., 0]
        result[..., 3, -1] = n[..., -1] - edge.sa_ratio
        result[..., 3, 1:-1] = (
            transport(n, (c + n / g**2) / _SIGMA, lagged[3, 1:-1])
            + (
                _CB2 / _SIGMA * self.derivative(n) / inner_g**2
                + (inner_c + inner_n / inner_g**2)
                / _SIGMA
                * self.derivative(g)
                / inner_g
            )
            * centered(n)
            + _CB1 * modified * positive_n
            - _CW1 * fw * positive_n**2 / distance**2
            - edge.viscosity_growth * inner_u * inner_n
        )
        return result

    def starting_guess(self, edge: StationEdge) -> FloatArray:
        u = np.tanh(0.33 * self.eta)
        recovery = 1.0 + np.sqrt(self.prandtl_number) * edge.heating / 2.0
        wall = recovery if edge.wall_ratio is None else edge.wall_ratio
        g = wall + (recovery - wall) * u + (1.0 - recovery) * u**2
        n = np.full_like(u, edge.sa_ratio)
        n[0] = 0.0
        return np.array(
            [self.cumulative(u), u, g, n, self.cumulative(g)], dtype=np.float64
        )

    def next_guess(self, last: Station, edge: StationEdge) -> FloatArray:
        guess = super().next_guess(last, edge)
        guess[3, -1] = edge.sa_ratio
        return guess

    def enforce_boundaries(self, state: FloatArray, edge: StationEdge) -> None:
        state[[0, 1, 3, 4], 0] = 0.0
        if edge.wall_ratio is not None:
            state[2, 0] = edge.wall_ratio
        state[1:4, -1] = (1.0, 1.0, edge.sa_ratio)

    def damping(self, state: FloatArray, update: FloatArray) -> float:
        damping = super().damping(state, update)
        decreasing = update[3, 1:-1] < 0.0
        if np.any(decreasing):
            damping = min(
                damping,
                float(
                    np.min(
                        -_POSITIVITY_FRACTION
                        * state[3, 1:-1][decreasing]
                        / update[3, 1:-1][decreasing]
                    )
                ),
            )
        # Picard can settle on a clipped n = 0 at the outer turbulent front, where
        # the unclipped equation has no admissible Newton step.
        return damping if damping >= _MINIMUM_DAMPING else 0.0

    def update_size(self, state: FloatArray, update: FloatArray) -> float:
        return max(
            float(np.max(np.abs(update[1]))),
            float(np.max(np.abs(update[2]))) / max(1.0, float(np.max(state[2]))),
            float(np.max(np.abs(update[3]))) / max(1.0, float(np.max(state[3]))),
        )

    def solve(
        self,
        guess: FloatArray,
        edge: StationEdge,
        coefficient: float,
        lagged: FloatArray,
    ) -> Station | None:
        station = self.newton(guess, edge, coefficient, lagged)
        if station is not None:
            return station
        self.picard_fallbacks += 1
        return self.picard(guess, edge, coefficient, lagged)

    def _tridiagonal(
        self,
        lower: FloatArray,
        diagonal: FloatArray,
        upper: FloatArray,
        wall: tuple[float, float, float],
        edge_value: float,
        right_hand_side: FloatArray,
    ) -> FloatArray:
        """Solve interior rows with a wall row ``(a, b, c)``: ``a v0 + b v1 = c``."""
        size = len(diagonal) + 2
        bands = np.zeros((3, size), dtype=np.float64)
        vector = np.zeros(size, dtype=np.float64)
        bands[1, 0], bands[0, 1], vector[0] = wall
        bands[1, -1], vector[-1] = 1.0, edge_value
        bands[1, 1:-1] = diagonal
        bands[0, 2:] = upper
        bands[2, :-2] = lower
        vector[1:-1] = right_hand_side
        return self._solve_tridiagonal(bands, vector)

    @staticmethod
    def _solve_tridiagonal(bands: FloatArray, vector: FloatArray) -> FloatArray:
        return np.asarray(solve_banded((1, 1), bands, vector), dtype=np.float64)

    def picard(
        self,
        guess: FloatArray,
        edge: StationEdge,
        coefficient: float,
        lagged: FloatArray,
    ) -> Station | None:
        """Run the relaxed Picard iteration; None requests a smaller step."""
        try:
            return self._picard(guess, edge, coefficient, lagged)
        except (ValueError, FloatingPointError, np.linalg.LinAlgError):
            return None

    def _picard(
        self,
        guess: FloatArray,
        edge: StationEdge,
        coefficient: float,
        lagged: FloatArray,
    ) -> Station | None:
        h, width = self.h, self.width
        scaled_xi = 2.0 * edge.xi
        beta, growth = edge.beta, edge.viscosity_growth
        implicit_beta, implicit_growth = max(beta, 0.0), max(growth, 0.0)
        u, g, n = guess[1].copy(), guess[2].copy(), guess[3].copy()
        wall_velocity = (1.0, 0.0, 0.0)
        wall_nu = (1.0, 0.0, 0.0)
        for _ in range(_PICARD_ITERATIONS):
            f = self.cumulative(u)
            convection = f[1:-1] + scaled_xi * (coefficient * f[1:-1] + lagged[0, 1:-1])
            c = self.chapman_rubesin(g, edge.temperature)
            chi = n / (c * g**2)
            eddy = chi * _fv1(chi)
            effective = c * (1.0 + eddy)
            inner_u, inner_g = u[1:-1], g[1:-1]

            face = (effective[1:] + effective[:-1]) / 2.0
            time_sink = scaled_xi * (2.0 * coefficient * inner_u + lagged[1, 1:-1])
            implicit_time = np.maximum(time_sink, 0.0)
            u_new = self._tridiagonal(
                face[:-1] / h[:-1] / width - convection / (2.0 * width),
                -(face[:-1] / h[:-1] + face[1:] / h[1:]) / width
                - implicit_time
                - 2.0 * implicit_beta * inner_u,
                face[1:] / h[1:] / width + convection / (2.0 * width),
                wall_velocity,
                1.0,
                -scaled_xi * coefficient * inner_u**2
                + (time_sink - implicit_time) * inner_u
                - beta * inner_g
                - implicit_beta * inner_u**2
                + (beta - implicit_beta) * inner_u**2,
            )

            conductivity = (
                c
                / self.prandtl_number
                * (1.0 + eddy * self.prandtl_number / self.turbulent_prandtl_number)
            )
            face = (conductivity[1:] + conductivity[:-1]) / 2.0
            slope = self.derivative(u_new)
            if edge.wall_ratio is None:
                wall_slope = (u_new[1] - u_new[0]) / h[0]
                wall_energy = (
                    -face[0] / h[0],
                    face[0] / h[0],
                    -h[0] / 2.0 * edge.heating * effective[0] * wall_slope**2,
                )
            else:
                wall_energy = (1.0, 0.0, edge.wall_ratio)
            g_new = self._tridiagonal(
                face[:-1] / h[:-1] / width - convection / (2.0 * width),
                -(face[:-1] / h[:-1] + face[1:] / h[1:]) / width
                - scaled_xi * coefficient * u_new[1:-1],
                face[1:] / h[1:] / width + convection / (2.0 * width),
                wall_energy,
                1.0,
                scaled_xi * u_new[1:-1] * lagged[2, 1:-1]
                - edge.heating * effective[1:-1] * slope**2,
            )
            g_new = np.maximum(g_new, _MINIMUM_TEMPERATURE_RATIO)

            # Production is lagged and destruction implicit, with a pseudo-time
            # relaxation bounded by the source stiffness.
            distance = self.cumulative(g)[1:-1]
            positive_n, inner_c = np.maximum(n[1:-1], 0.0), c[1:-1]
            vorticity = edge.reynolds_root * np.abs(slope) / inner_g
            modified, extra, fw = _source_terms(
                vorticity, positive_n, inner_c, inner_g, distance
            )
            destruction = _CW1 * fw * positive_n / distance**2
            relaxation = (
                4.0 * scaled_xi * coefficient
                + _CB1 * (vorticity + 2.0 * np.abs(extra))
                + destruction
            )
            diffusivity = c + n / g**2
            face = (diffusivity[1:] + diffusivity[:-1]) / 2.0
            sa_convection = (
                convection
                + _CB2 / _SIGMA * self.derivative(n) / inner_g**2
                + diffusivity[1:-1] / _SIGMA * self.derivative(g) / inner_g
            )
            n_new = self._tridiagonal(
                face[:-1] / _SIGMA / h[:-1] / width - sa_convection / (2.0 * width),
                -(face[:-1] / h[:-1] + face[1:] / h[1:]) / _SIGMA / width
                - scaled_xi * coefficient * u_new[1:-1]
                - implicit_growth * u_new[1:-1]
                - destruction
                - relaxation,
                face[1:] / _SIGMA / h[1:] / width + sa_convection / (2.0 * width),
                wall_nu,
                edge.sa_ratio,
                scaled_xi * u_new[1:-1] * lagged[3, 1:-1]
                - _CB1 * modified * positive_n
                - relaxation * positive_n
                + (growth - implicit_growth) * u_new[1:-1] * positive_n,
            )
            n_new = np.maximum(n_new, 0.0)

            change = max(
                float(np.max(np.abs(u_new - u))),
                float(np.max(np.abs(g_new - g))) / max(1.0, float(np.max(g))),
                float(np.max(np.abs(n_new - n))) / max(1.0, float(np.max(n))),
            )
            weight = _PICARD_UNDER_RELAXATION
            u = (1.0 - weight) * u + weight * u_new
            g = (1.0 - weight) * g + weight * g_new
            n = (1.0 - weight) * n + weight * n_new
            if change < _PICARD_TOLERANCE:
                state = np.array(
                    [self.cumulative(u), u, g, n, self.cumulative(g)], dtype=np.float64
                )
                return self.station(state, edge.xi)
        return None
