"""Warm CPU timings for direct laminar similarity; not a CI performance test.

Run from the repository root with ``uv run python
 docs/scripts/benchmark_laminar_similarity.py``. States use seed 20261102,
Me uniform on 0--20, Te log-uniform on 30--2000 K, and Tw/Taw log-uniform
on 0.05--1.5. Solved adiabatic temperatures and imports are prepared before
timing; no physical state solutions are cached. Times include exact grouping,
refinement and any fallback. The public API additionally solves adiabatic
states for recovery temperatures when an isothermal wall is prescribed.
"""

import json
from time import perf_counter

import numpy as np

from aerophysics._laminar_similarity import similarity_coefficients
from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    CompressibilityCorrection,
    flat_plate_boundary_layer,
)
from aerophysics.transport import AIR_VISCOSITY


def main() -> None:
    rng = np.random.default_rng(20261102)
    count = 10000
    mach = rng.uniform(0, 20, count)
    temperature = np.exp(rng.uniform(np.log(30), np.log(2000), count))
    ratio = np.exp(rng.uniform(np.log(0.05), np.log(1.5), count))
    states = np.column_stack(
        (
            mach,
            temperature,
            np.full(count, np.nan),
            np.full(count, 0.72),
            np.full(count, 1.4),
        )
    )
    print("Preparing adiabatic wall temperatures", flush=True)
    states[:, 2] = ratio * similarity_coefficients(states, AIR_VISCOSITY)[:, 4]
    timings = {}
    for size in (1, 100, 10000):
        start = perf_counter()
        result = similarity_coefficients(states[:size], AIR_VISCOSITY)
        timings[f"direct_{size}_states_seconds"] = perf_counter() - start
        assert np.all(np.isfinite(result))
        print(json.dumps(timings), flush=True)
    viscosity = AIR_VISCOSITY.dynamic_viscosity(temperature)
    start = perf_counter()
    flat_plate_boundary_layer(
        1.0,
        1e8 * viscosity,
        1.0,
        viscosity,
        regime=BoundaryLayerRegime.LAMINAR,
        compressibility_correction=CompressibilityCorrection.ECKERT,
        mach=mach,
        edge_temperature=temperature,
        wall_temperature=states[:, 2] * temperature,
    )
    timings["public_10000_isothermal_states_seconds"] = perf_counter() - start
    for wall in (None, 300.0):
        start = perf_counter()
        flat_plate_boundary_layer(
            np.linspace(0.1, 1.0, count),
            1e8 * viscosity[0],
            1.0,
            viscosity[0],
            regime=BoundaryLayerRegime.LAMINAR,
            compressibility_correction=CompressibilityCorrection.ECKERT,
            mach=mach[0],
            edge_temperature=temperature[0],
            wall_temperature=wall,
        )
        label = "adiabatic" if wall is None else "isothermal"
        timings[f"public_10000_x_points_{label}_seconds"] = perf_counter() - start
    print(json.dumps(timings, indent=2), flush=True)


if __name__ == "__main__":
    main()
