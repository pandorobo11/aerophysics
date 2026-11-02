"""Independent perfect-gas ZPG laminar similarity reference (SI temperatures).

Illingworth/Levy--Lees coordinate z uses (C f'')' + f f'' = 0, with
u/Ue=f', g=T/Te, and C=(mu/mu_e)/g. Thus the conventional Blasius coordinate
is eta=sqrt(2)*z in the constant-property limit. No production boundary-layer
correlation or fitted thickness coefficient is used in this solver.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import quad, solve_bvp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / "tests/reference_data/viscous_flow"
STEM = "laminar_compressible_similarity"
EDGE_TEMPERATURE = 240.0
PRANDTL = 0.72
GAMMA = 1.4
SUTHERLAND_TEMPERATURE = 110.4
MACH_NUMBERS = (0.3, 2.0, 4.0, 6.0, 10.0, 15.0)
WALL_RATIOS = (0.1, 0.3, 0.5, 1.0)


@dataclass
class SimilaritySolution:
    """Converged BVP and independently integrated edge-normalized quantities."""

    mach: float
    bvp: Any

    def integral(
        self, component: int, power: int = 1, end: float | None = None
    ) -> float:
        upper = float(self.bvp.x[-1]) if end is None else end
        value, _ = quad(
            lambda z: float(self.bvp.sol(z)[component]) ** power,
            0.0,
            upper,
            epsabs=1e-10,
            epsrel=1e-10,
        )
        return float(value)

    def quantities(self) -> dict[str, float]:
        root_two = float(np.sqrt(2.0))
        z99 = brentq(
            lambda z: float(self.bvp.sol(z)[1]) - 0.99,
            0.0,
            float(self.bvp.x[-1]),
            xtol=1e-13,
        )
        # Difference integrands avoid subtraction of large far-field integrals.
        displacement, _ = quad(
            lambda z: float(self.bvp.sol(z)[3] - self.bvp.sol(z)[1]),
            0.0,
            float(self.bvp.x[-1]),
            epsabs=1e-10,
            epsrel=1e-10,
        )
        momentum, _ = quad(
            lambda z: float(self.bvp.sol(z)[1] * (1.0 - self.bvp.sol(z)[1])),
            0.0,
            float(self.bvp.x[-1]),
            epsabs=1e-10,
            epsrel=1e-10,
        )
        return {
            "wall_temperature_K": EDGE_TEMPERATURE * float(self.bvp.y[3, 0]),
            "cf_sqrt_re": root_two * float(self.bvp.y[2, 0]),
            "delta99_sqrt_re_over_x": root_two * self.integral(3, end=z99),
            "displacement_sqrt_re_over_x": root_two * float(displacement),
            "momentum_sqrt_re_over_x": root_two * float(momentum),
            "eta99": root_two * z99,
            "I1": root_two * self.integral(1, end=z99),
            "I2": root_two * self.integral(1, power=2, end=z99),
        }


def solve_similarity(
    mach: float,
    wall_ratio: float | None = None,
    *,
    previous: SimilaritySolution | None = None,
    extent: float = 14.0,
    tolerance: float = 1e-8,
) -> SimilaritySolution:
    """Solve with prescribed Tw/Te, or zero wall heat flux if None.

    State = (f, f', C f'', g, C g'/Pr); f(0)=f'(0)=0,
    f'(infinity)=g(infinity)=1. The fifth wall condition is g(0)=Tw/Te
    or g'(0)=0. Sutherland's viscosity ratio is evaluated independently
    with the package default air constant 110.4 K.
    """
    dissipation = (GAMMA - 1.0) * mach**2

    def equations(
        z: NDArray[np.float64], state: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        del z
        stream, speed, shear, temperature, heat = state
        # Newton iterates may leave the physical domain; the converged
        # temperature is checked below and never clipped in the output.
        positive_temperature = np.maximum(temperature, 1e-8)
        chapman_rubesin = (
            np.sqrt(positive_temperature)
            * (EDGE_TEMPERATURE + SUTHERLAND_TEMPERATURE)
            / (EDGE_TEMPERATURE * positive_temperature + SUTHERLAND_TEMPERATURE)
        )
        speed_gradient = shear / chapman_rubesin
        temperature_gradient = PRANDTL * heat / chapman_rubesin
        return np.vstack(
            (
                speed,
                speed_gradient,
                -stream * speed_gradient,
                temperature_gradient,
                -stream * temperature_gradient - dissipation * shear * speed_gradient,
            )
        )

    def boundaries(
        left: NDArray[np.float64], right: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        thermal = left[4] if wall_ratio is None else left[3] - wall_ratio
        return np.array((left[0], left[1], thermal, right[1] - 1.0, right[3] - 1.0))

    nodes = np.linspace(0.0, extent, 401)
    if previous is None:
        speed = np.tanh(nodes)
        guess = np.vstack(
            (
                np.log(np.cosh(nodes)),
                speed,
                0.47 * np.exp(-nodes),
                np.ones_like(nodes),
                np.zeros_like(nodes),
            )
        )
    else:
        guess = previous.bvp.sol(nodes)
    result = solve_bvp(
        equations, boundaries, nodes, guess, tol=tolerance, max_nodes=30_000
    )
    if not result.success:
        raise RuntimeError(
            f"Similarity BVP failed at M={mach}, Tw/Te={wall_ratio}: {result.message}"
        )
    if np.min(result.y[3]) <= 0.0:
        raise RuntimeError("Similarity solution has nonpositive temperature")
    return SimilaritySolution(mach, result)


def reference_solutions() -> list[dict[str, float]]:
    """Continue in Mach along adiabatic walls, then cool each wall gradually."""
    rows = []
    adiabatic = solve_similarity(0.0)
    for mach in MACH_NUMBERS:
        for intermediate in np.linspace(
            adiabatic.mach, mach, max(2, int(np.ceil(2 * (mach - adiabatic.mach))) + 1)
        )[1:]:
            adiabatic = solve_similarity(float(intermediate), previous=adiabatic)
        taw = adiabatic.quantities()["wall_temperature_K"]
        solution = adiabatic
        last_ratio = 1.0
        by_ratio = {}
        for ratio in reversed(WALL_RATIOS):
            if ratio != 1.0:
                for intermediate in np.linspace(last_ratio, ratio, 5)[1:]:
                    solution = solve_similarity(
                        mach,
                        float(intermediate) * taw / EDGE_TEMPERATURE,
                        previous=solution,
                    )
            quantities = solution.quantities()
            by_ratio[ratio] = {
                "mach": mach,
                "wall_to_adiabatic_temperature_ratio": ratio,
                "wall_temperature_K": quantities["wall_temperature_K"],
                "adiabatic_temperature_K": taw,
                **{
                    key: value
                    for key, value in quantities.items()
                    if key not in ("wall_temperature_K", "eta99", "I1", "I2")
                },
            }
            last_ratio = ratio
        rows.extend(by_ratio[ratio] for ratio in WALL_RATIOS)
    return rows


def main() -> None:
    rows = reference_solutions()
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    with (DIRECTORY / f"{STEM}.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(
            {key: f"{value:.12g}" for key, value in row.items()} for row in rows
        )
    metadata = {
        "sources": [
            "Illingworth (1949), doi:10.1098/rspa.1949.0153",
            "Howarth (1948), doi:10.1098/rspa.1948.0064",
            "Chapman and Rubesin (1949), doi:10.2514/8.11853",
        ],
        "role": (
            "independent numerical perfect-gas similarity solution; not measurements"
        ),
        "command": "uv run python docs/scripts/build_laminar_similarity_reference.py",
        "equations": [
            "(C f'')' + f f'' = 0",
            "(C g'/Pr)' + f g' + (gamma-1) Me^2 C f''^2 = 0",
        ],
        "coordinate": "z: dy sqrt(Re_x)/x = sqrt(2) g dz; Blasius eta=sqrt(2) z",
        "boundary_conditions": "f(0)=f'(0)=0; f'(infinity)=g(infinity)=1; "
        "g'(0)=0 adiabatic or g(0)=Tw/Te isothermal",
        "normalization": "Re_x=rho_e U_e x/mu_e; cf uses edge dynamic pressure; "
        "delta99 at u/Ue=0.99; delta* and theta integrated to infinity",
        "temperature_units": "K",
        "edge_temperature_K": EDGE_TEMPERATURE,
        "prandtl_number": PRANDTL,
        "heat_capacity_ratio": GAMMA,
        "transport": {
            "model": "Sutherland",
            "sutherland_temperature_K": SUTHERLAND_TEMPERATURE,
            "mu_over_mu_e": "g^1.5 (Te+S)/(g Te+S)",
        },
        "wall_temperature": "Tw/Taw uses the solved adiabatic temperature; "
        "ratio 1 uses zero heat flux, not sqrt(Pr) recovery",
        "numerics": {
            "solver": "scipy.integrate.solve_bvp",
            "tolerance": 1e-8,
            "domain": [0.0, 14.0],
            "initial_nodes": 401,
            "max_nodes": 30000,
            "continuation": "Mach increments <=0.5, then decreasing Tw/Taw",
            "quadrature": "adaptive quad, epsabs=epsrel=1e-10; brentq for u/Ue=0.99",
        },
        "limitations": "constant cp and Pr, perfect gas, ZPG, no viscous interaction; "
        "high-temperature rows are mathematical checks, not real-gas validation",
    }
    (DIRECTORY / f"{STEM}.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
