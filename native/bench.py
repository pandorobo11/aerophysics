"""Compare the Python and optional Rust laminar station solvers."""

import os
from time import perf_counter

import numpy as np

from aerophysics import (
    AIR,
    BoundaryLayerEdge,
    BoundaryLayerRegime,
    marching_boundary_layer,
)


def flat_edge():
    return BoundaryLayerEdge.from_pressure(
        np.linspace(0.0, 1.0, 51),
        np.full(51, 2e4),
        total_pressure=1.6e5,
        total_temperature=500.0,
    )


def smooth_edge(count):
    distance = np.linspace(0.0, 2.0, count)
    mach = 2.0 + (1.0 + np.tanh((distance - 1.0) / 0.1)) / 2.0
    gamma = AIR.heat_capacity_ratio
    pressure = 2e5 * (1.0 + (gamma - 1.0) / 2.0 * mach**2) ** (-gamma / (gamma - 1.0))
    return BoundaryLayerEdge.from_pressure(
        distance,
        pressure,
        total_pressure=2e5,
        total_temperature=450.0,
    )


def relative_difference(a, b):
    return float(np.max(np.abs(a - b) / np.maximum(np.abs(a), 1e-300)))


def run(edge, wall, mode):
    os.environ["AEROPHYSICS_NATIVE"] = mode
    return marching_boundary_layer(
        edge, regime=BoundaryLayerRegime.LAMINAR, wall_temperature=wall
    )


def main():
    print("case | Python (s) | native (s) | speedup | max rel Cf | max rel wall")
    print("--- | ---: | ---: | ---: | ---: | ---:")
    for name, edge, wall in (
        ("flat-51", flat_edge(), None),
        ("smooth-201", smooth_edge(201), 300.0),
        ("smooth-1601", smooth_edge(1601), 300.0),
    ):
        timings = []
        results = []
        for mode in ("0", "1"):
            run(edge, wall, mode)
            samples = []
            for _ in range(5):
                start = perf_counter()
                result = run(edge, wall, mode)
                samples.append(perf_counter() - start)
            timings.append(min(samples))
            results.append(result)
        python, native = results
        cf = relative_difference(
            python.local_skin_friction_coefficient,
            native.local_skin_friction_coefficient,
        )
        wall_name = "wall_temperature" if wall is None else "wall_heat_flux"
        thermal = relative_difference(
            getattr(python, wall_name), getattr(native, wall_name)
        )
        print(
            f"{name} | {timings[0]:.3f} | {timings[1]:.3f} | "
            f"{timings[0] / timings[1]:.2f} | {cf:.3e} | {thermal:.3e}",
            flush=True,
        )


if __name__ == "__main__":
    main()
