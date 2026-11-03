"""Independent compressible ZPG turbulent boundary-layer reference with SA-noft2.

The thin-layer equations are marched downstream from the leading edge for a
calorically perfect gas with constant edge state. Levy--Lees variables are
xi = rho_e mu_e U_e x and eta = U_e/sqrt(2 xi) int rho dy, with u/Ue = f',
g = T/Te, n = nu_tilde/nu_e, C = rho mu/(rho_e mu_e), and Y = int g d(eta),
so that y = x sqrt(2/Re_x) Y:

  (C(1+eps) f'')' + f f''                         = 2x (f' f'_x - f'' f_x)
  (C/Pr (1 + eps Pr/Prt) g')' + f g'
      + (gamma-1) Me^2 C (1+eps) f''^2            = 2x (f' g_x - g' f_x)
  ((C + n/g^2) n')'/sigma + cb2 n'^2/(sigma g^2) + f n'
      + cb1 St n - cw1 fw n^2/Y^2                 = 2x (f' n_x - n' f_x)

with St = sqrt(2 Re_x)|f''|/g + n fv2/(kappa^2 Y^2), chi = n/(C g^2), and
eps = mu_t/mu = chi fv1. The SA transport is the density-weighted form with
rho inside the diffusion term, the trip term ft2 is omitted, and the edge value
is nu_tilde = 3 nu_e. At x = 0 every streamwise term vanishes, so the leading
edge is a regular starting profile in these variables.

No production boundary-layer correlation, fitted thickness coefficient or
package model is used in this solver; the Sutherland law is evaluated here.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import solve_banded

Array = NDArray[np.float64]

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / "tests/reference_data/viscous_flow"
STEM = "turbulent_sa_boundary_layer"

PRANDTL = 0.72
TURBULENT_PRANDTL = 0.9
GAMMA = 1.4
SUTHERLAND_TEMPERATURE = 110.4
EDGE_SA_RATIO = 3.0
REYNOLDS_NUMBERS = (1.0e6, 1.0e7, 1.0e8)
UNIT_REYNOLDS = 1.0e7
STEPS = 300

CB1, CB2, SIGMA, KAPPA, CV1 = 0.1355, 0.622, 2.0 / 3.0, 0.41, 7.1
CW2, CW3 = 0.3, 2.0
CW1 = CB1 / KAPPA**2 + (1.0 + CB2) / SIGMA
# Smooth S-tilde limiter of Allmaras, Johnson and Spalart (2012).
CV2, CV3 = 0.7, 0.9

ITERATIONS = 4000
TOLERANCE = 1e-9

# (label, edge Mach number, edge temperature [K], Tw/Tr or None for Tw = 300 K)
CASES: tuple[tuple[str, float, float, float | None], ...] = (
    ("M0.3-AD", 0.3, 288.15, 1.0),
    ("M2-AD", 2.0, 220.0, 1.0),
    ("M5-TW300", 5.0, 220.0, None),
    ("M6-0.76", 6.0, 100.0, 0.76),
    ("M6-0.25", 6.0, 100.0, 0.25),
    ("M10-0.20", 10.0, 60.0, 0.20),
)


def recovery_temperature(mach: float, edge_temperature: float) -> float:
    """Turbulent recovery temperature with r = Pr**(1/3), as in the VD2 model."""
    return edge_temperature * (
        1.0 + PRANDTL ** (1.0 / 3.0) * (GAMMA - 1.0) / 2.0 * mach**2
    )


def grid(ratio: float = 1.04, eta_max: float = 90.0, first: float = 7e-4) -> Array:
    """Geometric eta grid; the first point lies well inside y+ = 1."""
    count = math.ceil(math.log(eta_max / first) / math.log(ratio)) + 2
    eta = np.concatenate(([0.0], first * ratio ** np.arange(count)))
    return np.asarray(eta[eta <= eta_max], dtype=np.float64)


def _fv1(chi: Array) -> Array:
    return np.asarray(chi**3 / (chi**3 + CV1**3), dtype=np.float64)


def _cumulative(values: Array, eta: Array) -> Array:
    increments = np.diff(eta) * (values[1:] + values[:-1]) / 2.0
    return np.concatenate(([0.0], np.cumsum(increments)))


def _tridiagonal(
    lower: Array, diagonal: Array, upper: Array, wall: float, edge: float, rhs: Array
) -> Array:
    """Solve interior rows with Dirichlet values at both ends."""
    size = len(diagonal) + 2
    bands = np.zeros((3, size))
    vector = np.zeros(size)
    bands[1, 0], vector[0] = 1.0, wall
    bands[1, -1], vector[-1] = 1.0, edge
    bands[1, 1:-1] = diagonal
    bands[0, 2:] = upper
    bands[2, :-2] = lower
    vector[1:-1] = rhs
    return np.asarray(solve_banded((1, 1), bands, vector), dtype=np.float64)


@dataclass(frozen=True)
class Station:
    """Profiles u/Ue, T/Te and nu_tilde/nu_e, and the stream function f."""

    u: Array
    g: Array
    n: Array
    f: Array


@dataclass(frozen=True)
class Quantities:
    """Local results normalized by x and the edge dynamic pressure."""

    cf: float
    delta99_over_x: float
    displacement99_over_x: float
    momentum_over_x: float


class TurbulentSABoundaryLayer:
    """Downstream march for one constant edge state and wall temperature."""

    def __init__(
        self,
        mach: float,
        edge_temperature: float,
        wall_temperature: float,
        *,
        unit_reynolds: float = UNIT_REYNOLDS,
        edge_sa_ratio: float = EDGE_SA_RATIO,
        eta: Array | None = None,
    ) -> None:
        self.mach = mach
        self.edge_temperature = edge_temperature
        self.wall_ratio = wall_temperature / edge_temperature
        self.unit_reynolds = unit_reynolds
        self.edge_sa_ratio = edge_sa_ratio
        self.eta = grid() if eta is None else eta
        self.heating = (GAMMA - 1.0) * mach**2

    def chapman_rubesin(self, g: Array) -> Array:
        """C = rho mu/(rho_e mu_e) for the Sutherland law."""
        te, s = self.edge_temperature, SUTHERLAND_TEMPERATURE
        return np.asarray(np.sqrt(g) * (te + s) / (g * te + s), dtype=np.float64)

    def _initial(self) -> Station:
        u = np.clip(np.tanh(0.33 * self.eta), 0.0, 1.0)
        recovery = 1.0 + math.sqrt(PRANDTL) * self.heating / 2.0
        g = self.wall_ratio + (recovery - self.wall_ratio) * u + (1.0 - recovery) * u**2
        n = np.full_like(u, self.edge_sa_ratio)
        n[0] = 0.0
        return Station(u, g, n, _cumulative(u, self.eta))

    def _solve(
        self, history: list[Station], positions: list[float], x: float
    ) -> Station:
        """Picard iteration at x with BDF2 streamwise derivatives (BDF1 at start)."""
        eta = self.eta
        h = np.diff(eta)
        width = (h[1:] + h[:-1]) / 2.0
        if len(positions) >= 2:
            a1, a2 = x - positions[-1], x - positions[-2]
            c0 = 1.0 / a1 + 1.0 / a2
            c1, c2 = -a2 / (a1 * (a2 - a1)), a1 / (a2 * (a2 - a1))

            def lagged(name: str) -> Array:
                current = getattr(history[-1], name)
                previous = getattr(history[-2], name)
                return np.asarray(c1 * current + c2 * previous, dtype=np.float64)

        else:
            c0 = 1.0 / (x - positions[-1])

            def lagged(name: str) -> Array:
                return np.asarray(-c0 * getattr(history[-1], name), dtype=np.float64)

        old_u, old_g, old_n, old_f = (lagged(name) for name in ("u", "g", "n", "f"))
        u, g, n = history[-1].u.copy(), history[-1].g.copy(), history[-1].n.copy()
        root = math.sqrt(2.0 * self.unit_reynolds * x)
        two_x = 2.0 * x
        change = math.inf
        for _ in range(ITERATIONS):
            f = _cumulative(u, eta)
            convection = f[1:-1] + two_x * (c0 * f[1:-1] + old_f[1:-1])
            c = self.chapman_rubesin(g)
            chi = n / (c * g**2)
            effective = c * (1.0 + chi * _fv1(chi))

            # Momentum: lagged C(1+eps), Newton-linearized u u_x.
            face = (effective[1:] + effective[:-1]) / 2.0
            u_new = _tridiagonal(
                face[:-1] / h[:-1] / width - convection / (2.0 * width),
                -(face[:-1] / h[:-1] + face[1:] / h[1:]) / width
                - two_x * (2.0 * c0 * u[1:-1] + old_u[1:-1]),
                face[1:] / h[1:] / width + convection / (2.0 * width),
                0.0,
                1.0,
                -two_x * c0 * u[1:-1] ** 2,
            )

            # Energy: lagged conductivity and dissipation.
            conductivity = (
                c / PRANDTL * (1.0 + chi * _fv1(chi) * PRANDTL / TURBULENT_PRANDTL)
            )
            face = (conductivity[1:] + conductivity[:-1]) / 2.0
            slope = np.gradient(u_new, eta)[1:-1]
            g_new = _tridiagonal(
                face[:-1] / h[:-1] / width - convection / (2.0 * width),
                -(face[:-1] / h[:-1] + face[1:] / h[1:]) / width
                - two_x * c0 * u_new[1:-1],
                face[1:] / h[1:] / width + convection / (2.0 * width),
                self.wall_ratio,
                1.0,
                two_x * u_new[1:-1] * old_g[1:-1]
                - self.heating * effective[1:-1] * slope**2,
            )
            g_new = np.maximum(g_new, 1e-3)

            # SA: production lagged, destruction implicit, and a pseudo-time
            # relaxation bounded by the source stiffness.
            distance = _cumulative(g, eta)[1:-1]
            nk, gk, ck = np.maximum(n[1:-1], 0.0), g[1:-1], c[1:-1]
            chik = nk / (ck * gk**2)
            fv2 = 1.0 - chik / (1.0 + chik * _fv1(chik))
            vorticity = root * np.abs(slope) / gk
            extra = nk * fv2 / (KAPPA**2 * distance**2)
            limited = vorticity + vorticity * (
                CV2**2 * vorticity + CV3 * extra
            ) / np.maximum((CV3 - 2.0 * CV2) * vorticity - extra, 1e-30)
            modified = np.where(extra >= -CV2 * vorticity, vorticity + extra, limited)
            modified = np.maximum(modified, 1e-12)
            r = np.minimum(
                nk / np.maximum(modified * KAPPA**2 * distance**2, 1e-30), 10.0
            )
            r6 = r + CW2 * (r**6 - r)
            fw = r6 * ((1.0 + CW3**6) / (r6**6 + CW3**6)) ** (1.0 / 6.0)
            destruction = CW1 * fw * nk / distance**2
            relax = (
                4.0 * two_x * c0 + CB1 * (vorticity + 2.0 * np.abs(extra)) + destruction
            )
            diffusivity = c + n / g**2
            face = (diffusivity[1:] + diffusivity[:-1]) / 2.0
            sa_convection = convection + CB2 / SIGMA * np.gradient(n, eta)[1:-1] / gk**2
            n_new = _tridiagonal(
                face[:-1] / SIGMA / h[:-1] / width - sa_convection / (2.0 * width),
                -(face[:-1] / h[:-1] + face[1:] / h[1:]) / SIGMA / width
                - two_x * c0 * u_new[1:-1]
                - destruction
                - relax,
                face[1:] / SIGMA / h[1:] / width + sa_convection / (2.0 * width),
                0.0,
                self.edge_sa_ratio,
                two_x * u_new[1:-1] * old_n[1:-1] - CB1 * modified * nk - relax * nk,
            )
            n_new = np.maximum(n_new, 0.0)

            change = max(
                float(np.max(np.abs(u_new - u))),
                float(np.max(np.abs(g_new - g))) / max(1.0, float(g.max())),
                float(np.max(np.abs(n_new - n))) / max(1.0, float(n.max())),
            )
            u = 0.5 * (u + u_new)
            g = 0.5 * (g + g_new)
            n = 0.5 * (n + n_new)
            if change < TOLERANCE:
                return Station(u, g, n, _cumulative(u, eta))
        raise RuntimeError(f"station x={x:g} did not converge (change {change:.1e})")

    def march(
        self, x_end: float, *, steps: int = STEPS, include: tuple[float, ...] = ()
    ) -> tuple[Array, list[Station]]:
        """Return stations from x = 0 to x_end, with the requested x inserted."""
        start = self._solve([self._initial()], [-1.0], 0.0)
        positions = np.unique(
            np.concatenate(
                ([0.0], np.geomspace(x_end * 1e-8, x_end, steps - 1), include)
            )
        )
        history, done = [start], [0.0]
        for x in positions[1:]:
            history.append(self._solve(history, done, float(x)))
            done.append(float(x))
        return positions, history

    def quantities(self, x: float, station: Station) -> Quantities:
        """cf, delta99, displacement through delta99, and full momentum thickness."""
        eta, u, g = self.eta, station.u, station.g
        reynolds = self.unit_reynolds * x
        h1, h2 = eta[1], eta[2]
        wall_slope = (
            -(h1 + h2) / (h1 * h2) * u[0]
            + h2 / (h1 * (h2 - h1)) * u[1]
            - h1 / (h2 * (h2 - h1)) * u[2]
        )
        wall_c = float(self.chapman_rubesin(g[:1])[0])
        scale = math.sqrt(2.0 / reynolds)
        edge99 = float(np.interp(0.99, u, eta))
        fine = np.linspace(0.0, edge99, 4001)
        inner_u, inner_g = np.interp(fine, eta, u), np.interp(fine, eta, g)
        return Quantities(
            cf=math.sqrt(2.0) * wall_c * float(wall_slope) / math.sqrt(reynolds),
            delta99_over_x=float(np.interp(edge99, eta, _cumulative(g, eta))) * scale,
            displacement99_over_x=float(np.trapezoid(inner_g - inner_u, fine)) * scale,
            momentum_over_x=float(np.trapezoid(u * (1.0 - u), eta)) * scale,
        )


def wall_temperature(
    mach: float, edge_temperature: float, ratio: float | None
) -> float:
    return (
        300.0 if ratio is None else ratio * recovery_temperature(mach, edge_temperature)
    )


def reference_solutions() -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    targets = tuple(reynolds / UNIT_REYNOLDS for reynolds in REYNOLDS_NUMBERS)
    for label, mach, edge_temperature, ratio in CASES:
        wall = wall_temperature(mach, edge_temperature, ratio)
        solver = TurbulentSABoundaryLayer(mach, edge_temperature, wall)
        positions, history = solver.march(targets[-1], include=targets)
        for reynolds, x in zip(REYNOLDS_NUMBERS, targets, strict=True):
            index = int(np.argmin(np.abs(positions - x)))
            result = solver.quantities(float(positions[index]), history[index])
            rows.append(
                {
                    "case": label,
                    "mach": mach,
                    "edge_temperature_K": edge_temperature,
                    "wall_temperature_K": wall,
                    "wall_to_recovery_temperature_ratio": (
                        wall / recovery_temperature(mach, edge_temperature)
                    ),
                    "reynolds_number": reynolds,
                    "cf": result.cf,
                    "delta99_over_x": result.delta99_over_x,
                    "displacement99_over_x": result.displacement99_over_x,
                    "momentum_over_x": result.momentum_over_x,
                }
            )
    return rows


def _format(value: float | str) -> str:
    return value if isinstance(value, str) else f"{value:.8g}"


def main() -> None:
    rows = reference_solutions()
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    with (DIRECTORY / f"{STEM}.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(
            {key: _format(value) for key, value in row.items()} for row in rows
        )
    metadata = {
        "sources": [
            "Spalart and Allmaras (1994), La Recherche Aerospatiale 1, 5-21",
            "Allmaras, Johnson and Spalart (2012), ICCFD7-1902 (S-tilde limiter)",
            "Levy (1954) and Lees (1956) boundary-layer transformation",
        ],
        "role": (
            "independent numerical RANS boundary-layer solution with SA-noft2; "
            "not measurements and not an acceptance target for correlations"
        ),
        "command": "uv run python docs/scripts/build_turbulent_sa_reference.py",
        "flow": "smooth ZPG flat plate, calorically perfect gas, constant edge state",
        "prandtl_number": PRANDTL,
        "turbulent_prandtl_number": TURBULENT_PRANDTL,
        "heat_capacity_ratio": GAMMA,
        "transport": {
            "model": "Sutherland",
            "sutherland_temperature_K": SUTHERLAND_TEMPERATURE,
            "mu_over_mu_e": "g^1.5 (Te+S)/(g Te+S)",
        },
        "turbulence": (
            "SA-noft2, density-weighted with rho inside diffusion; "
            f"edge nu_tilde/nu_e = {EDGE_SA_RATIO:g}; marched from the leading edge"
        ),
        "wall_temperature": (
            "Tw = ratio * Tr with Tr = Te (1 + Pr^(1/3) (gamma-1)/2 Me^2), "
            "or Tw = 300 K for M5-TW300"
        ),
        "normalization": (
            "Re_x = rho_e U_e x/mu_e; cf uses edge dynamic pressure; delta99 at "
            "u/Ue = 0.99; displacement integrated through delta99; momentum "
            "thickness integrated to the edge"
        ),
        "numerics": {
            "eta_grid": "geometric, first point 7e-4, ratio 1.04, eta_max 90",
            "streamwise": (
                f"BDF2, {STEPS} geometric stations from x_end 1e-8 to x_end, "
                "output x inserted"
            ),
            "unit_reynolds_per_m": UNIT_REYNOLDS,
            "picard_tolerance": TOLERANCE,
        },
        "verification": (
            "laminar limit (edge nu_tilde = 0) matches the compressible similarity "
            "solution within 0.01%; halving the streamwise step or refining the "
            "eta grid changes cf and thicknesses by less than 0.1%"
        ),
    }
    (DIRECTORY / f"{STEM}.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
