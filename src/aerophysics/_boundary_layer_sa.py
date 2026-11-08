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

from typing import cast

import numpy as np

from aerophysics._array import FloatArray
from aerophysics._banded import solve_tridiagonal
from aerophysics._boundary_layer_march import (
    Discretization,
    Station,
    StationEdge,
    _viscosity_with_derivative,
)
from aerophysics.transport import DynamicViscosityModel, KeyesModel, SutherlandModel

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
    return chi**3 / (chi**3 + _CV1**3)


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
    return modified, extra, fw


def _floor_derivative(
    value: FloatArray, derivative: FloatArray, floor: float
) -> FloatArray:
    """Differentiate max(value, floor) in each positive column direction."""
    return np.where(
        value > floor,
        derivative,
        np.where(value == floor, np.maximum(derivative, 0.0), 0.0),
    )


def _source_derivatives(
    vorticity: FloatArray,
    n: FloatArray,
    c: FloatArray,
    g: FloatArray,
    distance: FloatArray,
    directions: tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray],
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    """Return modified, fw and their directional derivatives, node-vectorized.

    Keep the residual's _source_terms untouched. At nonsmooth ties each column
    uses the right-sided derivative (a positive perturbation of that unknown),
    including the sign of the incoming direction. Strictly clipped directions
    have zero derivative; at a floor/cap only directions leaving it survive.
    """
    dv, dn, dc, dg, dy = directions
    chi = n / (c * g**2)
    dchi = dn / (c * g**2) - chi * (dc / c + 2.0 * dg / g)
    fv1 = _fv1(chi)
    dfv1 = 3.0 * chi**2 * _CV1**3 / (chi**3 + _CV1**3) ** 2
    fv2 = 1.0 - chi / (1.0 + chi * fv1)
    dfv2 = -(1.0 - chi**2 * dfv1) / (1.0 + chi * fv1) ** 2 * dchi
    extra = n * fv2 / (_KAPPA**2 * distance**2)
    de = (dn * fv2 + n * dfv2) / (_KAPPA**2 * distance**2) - 2.0 * extra * dy / distance
    denominator = (_CV3 - 2.0 * _CV2) * vorticity - extra
    dd = _floor_derivative(denominator, (_CV3 - 2.0 * _CV2) * dv - de, 1e-30)
    denominator = np.maximum(denominator, 1e-30)
    numerator = vorticity * (_CV2**2 * vorticity + _CV3 * extra)
    limited = vorticity + numerator / denominator
    dl = (
        dv
        + ((2.0 * _CV2**2 * vorticity + _CV3 * extra) * dv + _CV3 * vorticity * de)
        / denominator
        - numerator / denominator * dd / denominator
    )
    # The two modified-vorticity branches have matching derivatives at their
    # positive-vorticity join; retain np.where's >= choice, as in the residual.
    regular = extra >= -_CV2 * vorticity
    raw_modified = np.where(regular, vorticity + extra, limited)
    dm = _floor_derivative(
        raw_modified, np.where(regular, dv + de, dl), _MINIMUM_VORTICITY
    )
    modified = np.maximum(raw_modified, _MINIMUM_VORTICITY)
    denominator = modified * _KAPPA**2 * distance**2
    dd = _floor_derivative(
        denominator,
        _KAPPA**2 * (dm * distance**2 + 2.0 * modified * distance * dy),
        1e-30,
    )
    denominator = np.maximum(denominator, 1e-30)
    raw_r = n / denominator
    dr = dn / denominator - raw_r * dd / denominator
    dr = np.where(
        raw_r < _MAX_R, dr, np.where(raw_r == _MAX_R, np.minimum(dr, 0.0), 0.0)
    )
    r = np.minimum(raw_r, _MAX_R)
    r6 = r + _CW2 * (r**6 - r)
    factor = ((1.0 + _CW3**6) / (r6**6 + _CW3**6)) ** (1.0 / 6.0)
    fw = r6 * factor
    dfw = factor * _CW3**6 / (r6**6 + _CW3**6) * (1.0 + _CW2 * (6.0 * r**5 - 1.0)) * dr
    return modified, fw, dm, dfw


class TurbulentDiscretization(Discretization):
    """One station of ``(f, u, g, n, Y)`` with a Picard fallback."""

    unknowns = 5
    _analytic_jacobian_enabled = True
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
            return (value[..., 2:] - value[..., :-2]) / (2.0 * width)

        def transport(
            value: FloatArray, diffusivity: FloatArray, history: FloatArray
        ) -> FloatArray:
            flux = (
                (diffusivity[..., 1:] + diffusivity[..., :-1])
                / 2.0
                * ((value[..., 1:] - value[..., :-1]) / h)
            )
            return cast(
                FloatArray,
                (flux[..., 1:] - flux[..., :-1]) / width
                + convection * centered(value)
                - scaled_xi * inner_u * (coefficient * value[..., 1:-1] + history),
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
        result[..., 0, 1:] = (f[..., 1:] - f[..., :-1]) - h / 2.0 * (
            u[..., 1:] + u[..., :-1]
        )
        result[..., 4, 0] = stretched[..., 0]
        result[..., 4, 1:] = (stretched[..., 1:] - stretched[..., :-1]) - h / 2.0 * (
            g[..., 1:] + g[..., :-1]
        )

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
            + edge.turbulence_production * _CB1 * modified * positive_n
            - _CW1 * fw * positive_n**2 / distance**2
            - edge.viscosity_growth * inner_u * inner_n
        )
        return result

    def _analytic_jacobian(
        self,
        state: FloatArray,
        edge: StationEdge,
        coefficient: float,
        lagged: FloatArray,
    ) -> FloatArray:
        """Differentiate the unchanged five-equation discrete residual."""
        f, u, g, n, distance = state
        model = cast(SutherlandModel | KeyesModel, self.viscosity_model)
        mu_edge = float(model.dynamic_viscosity(edge.temperature))
        mu, dmu = _viscosity_with_derivative(model, edge.temperature * g)
        c = mu / mu_edge / g
        dc = (edge.temperature * dmu / mu_edge - c) / g
        chi = n / (c * g**2)
        fv1 = _fv1(chi)
        eddy = chi * fv1
        de = fv1 + 3.0 * chi**3 * _CV1**3 / (chi**3 + _CV1**3) ** 2
        eddy_g = -de * chi * (dc / c + 2.0 / g)
        eddy_n = de / (c * g**2)
        effective = c * (1.0 + eddy)
        effective_g = dc * (1.0 + eddy) + c * eddy_g
        effective_n = c * eddy_n
        pr, prt = self.prandtl_number, self.turbulent_prandtl_number
        conductivity = c / pr * (1.0 + eddy * pr / prt)
        conductivity_g = dc / pr * (1.0 + eddy * pr / prt) + c * eddy_g / prt
        conductivity_n = c * eddy_n / prt
        diffusivity = np.array([effective, conductivity, (c + n / g**2) / _SIGMA])
        diffusivity_g = np.array(
            [effective_g, conductivity_g, (dc - 2.0 * n / g**3) / _SIGMA]
        )
        diffusivity_n = np.array([effective_n, conductivity_n, 1.0 / g**2 / _SIGMA])
        sx = 2.0 * edge.xi
        convection = f[1:-1] + sx * (coefficient * f[1:-1] + lagged[0, 1:-1])
        values = state[1:4]
        centered = (values[:, 2:] - values[:, :-2]) / (2.0 * self.width)
        slopes = np.diff(values) / self.h / 2.0
        left_slope, right_slope = (
            slopes[:, :-1] / self.width,
            slopes[:, 1:] / self.width,
        )
        face = (diffusivity[:, 1:] + diffusivity[:, :-1]) / 2.0 / self.h
        left, right = face[:, :-1] / self.width, face[:, 1:] / self.width
        # Blocks are indexed by stencil offset, equation, variable, interior node.
        blocks = np.zeros((3, 5, 5, len(g) - 2), dtype=np.float64)
        for k, (linear, flux, weight) in enumerate(
            (
                (left, -left_slope, -1.0 / (2.0 * self.width)),
                (-left - right, right_slope - left_slope, np.zeros_like(self.width)),
                (right, right_slope, 1.0 / (2.0 * self.width)),
            )
        ):
            neighbor = slice(k, len(g) - 2 + k)
            for equation in range(1, 4):
                blocks[k, equation, equation] = (
                    linear[equation - 1] + convection * weight
                )
            blocks[k, 1:4, 2] += diffusivity_g[:, neighbor] * flux
            blocks[k, 1:4, 3] += diffusivity_n[:, neighbor] * flux
        blocks[1, 1:4, 0] = (1.0 + sx * coefficient) * centered
        blocks[1, 1:4, 1] -= sx * (coefficient * values[:, 1:-1] + lagged[1:4, 1:-1])
        for equation in range(1, 4):
            blocks[1, equation, equation] -= sx * u[1:-1] * coefficient
        blocks[1, 1, 1] -= 2.0 * edge.beta * u[1:-1]
        blocks[1, 1, 2] += edge.beta

        du, dg, dn = self.derivative(u), self.derivative(g), self.derivative(n)
        blocks[:, 2, 1] += (
            2.0 * edge.heating * effective[1:-1] * du * self.gradient_weights
        )
        blocks[1, 2, 2] += edge.heating * effective_g[1:-1] * du**2
        blocks[1, 2, 3] += edge.heating * effective_n[1:-1] * du**2
        inner_g, inner_n = g[1:-1], n[1:-1]
        sa_diffusion = diffusivity[2, 1:-1]
        sa_convection = _CB2 / _SIGMA * dn / inner_g**2 + sa_diffusion * dg / inner_g
        for k, weight in enumerate(
            (
                -1.0 / (2.0 * self.width),
                np.zeros_like(self.width),
                1.0 / (2.0 * self.width),
            )
        ):
            blocks[k, 3, 3] += (
                sa_convection * weight
                + _CB2 / _SIGMA / inner_g**2 * self.gradient_weights[k] * centered[2]
            )
            blocks[k, 3, 2] += (
                sa_diffusion / inner_g * self.gradient_weights[k] * centered[2]
            )
        blocks[1, 3, 2] += (
            -2.0 * _CB2 / _SIGMA * dn / inner_g**3
            + (diffusivity_g[2, 1:-1] / inner_g - sa_diffusion / inner_g**2) * dg
        ) * centered[2]
        blocks[1, 3, 3] += diffusivity_n[2, 1:-1] * dg / inner_g * centered[2]

        positive_n = np.maximum(inner_n, 0.0)
        y = distance[1:-1]
        vorticity = edge.reynolds_root * np.abs(du) / inner_g
        # Six directions: u at i-1,i,i+1, then local g,n,Y. At du=0,
        # d|du|/d(u_j)+ = |weight_j|, including negative stencil weights.
        directions = np.zeros((5, 6, len(g) - 2), dtype=np.float64)
        directions[0, :3] = (
            edge.reynolds_root
            / inner_g
            * np.where(
                du == 0.0,
                np.abs(self.gradient_weights),
                np.sign(du) * self.gradient_weights,
            )
        )
        directions[0, 3] = -vorticity / inner_g
        directions[1, 4] = inner_n >= 0.0
        directions[2, 3] = dc[1:-1]
        directions[3, 3] = 1.0
        directions[4, 5] = 1.0
        modified, fw, dm, dfw = _source_derivatives(
            vorticity,
            positive_n,
            c[1:-1],
            inner_g,
            y,
            (directions[0], directions[1], directions[2], directions[3], directions[4]),
        )
        source = (
            edge.turbulence_production
            * _CB1
            * (dm * positive_n + modified * directions[1])
            - _CW1
            / y**2
            * (dfw * positive_n**2 + 2.0 * fw * positive_n * directions[1])
            + 2.0 * _CW1 * fw * positive_n**2 / y**3 * directions[4]
        )
        blocks[:, 3, 1] += source[:3]
        blocks[1, 3, 2:5] += source[3:]
        blocks[1, 3, 1] -= edge.viscosity_growth * inner_n
        blocks[1, 3, 3] -= edge.viscosity_growth * u[1:-1]

        bands = np.zeros((3 * self.band + 1, state.size), dtype=np.float64, order="F")
        diagonal = 2 * self.band
        for k in range(3):
            for equation in range(1, 4):
                for variable in range(5):
                    bands[
                        diagonal - 5 * (k - 1) + equation - variable,
                        5 * k + variable : 5 * (len(g) - 2 + k) + variable : 5,
                    ] = blocks[k, equation, variable]
        # Trapezoidal f and Y integration, including their wall values.
        for equation, variable in ((0, 1), (4, 2)):
            bands[diagonal, equation::5] = 1.0
            bands[diagonal + 5, equation:-5:5] = -1.0
            bands[diagonal + 5 + equation - variable, variable:-5:5] = -self.h / 2.0
            bands[diagonal + equation - variable, 5 + variable :: 5] = -self.h / 2.0
        bands[diagonal, 1] = bands[diagonal, 3] = 1.0
        bands[diagonal, -4:-1] = 1.0
        if edge.wall_ratio is None:
            h = self.h[0]
            slope = (u[1] - u[0]) / h
            wall_face = (conductivity[0] + conductivity[1]) / 2.0
            wall_dg = (g[1] - g[0]) / h
            bands[diagonal + 1, 1] = -edge.heating * effective[0] * slope
            bands[diagonal - 4, 6] = edge.heating * effective[0] * slope
            for variable, derivative, dissipation in (
                (2, conductivity_g, effective_g),
                (3, conductivity_n, effective_n),
            ):
                bands[diagonal + 2 - variable, variable] = (
                    derivative[0] / 2.0 * wall_dg
                    + h / 2.0 * edge.heating * dissipation[0] * slope**2
                )
                bands[diagonal - 3 - variable, 5 + variable] = (
                    derivative[1] / 2.0 * wall_dg
                )
            bands[diagonal, 2] -= wall_face / h
            bands[diagonal - 5, 7] += wall_face / h
        else:
            bands[diagonal, 2] = 1.0
        return bands

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
        return solve_tridiagonal(bands, vector)

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
                + edge.turbulence_production * _CB1 * (vorticity + 2.0 * np.abs(extra))
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
                - edge.turbulence_production * _CB1 * modified * positive_n
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
