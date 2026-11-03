"""Independent SA boundary-layer reference and its comparison with VD2."""

import csv
import runpy
import warnings
from collections.abc import Callable
from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics._laminar_similarity import similarity_coefficients
from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    TurbulentModel,
    flat_plate_boundary_layer,
)
from aerophysics.gas import AIR
from aerophysics.transport import AIR_VISCOSITY

_SOLVER = runpy.run_path(
    str(
        Path(__file__).resolve().parents[1]
        / "docs/scripts/build_turbulent_sa_reference.py"
    )
)
TurbulentSABoundaryLayer = cast(Any, _SOLVER["TurbulentSABoundaryLayer"])
reference_solutions = cast(
    Callable[[], list[dict[str, float | str]]], _SOLVER["reference_solutions"]
)

REFERENCE = (
    Path(__file__).parent
    / "reference_data/viscous_flow/turbulent_sa_boundary_layer.csv"
)
QUANTITIES = ("cf", "delta99_over_x", "displacement99_over_x", "momentum_over_x")


def reference_rows() -> list[dict[str, str]]:
    with REFERENCE.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def test_laminar_limit_matches_compressible_similarity() -> None:
    solver = TurbulentSABoundaryLayer(5.0, 220.0, 220.0, edge_sa_ratio=0.0)
    positions, history = solver.march(0.1, steps=20)
    result = solver.quantities(float(positions[-1]), history[-1])
    root = np.sqrt(solver.unit_reynolds * positions[-1])
    expected = similarity_coefficients(
        np.array([[5.0, 220.0, 1.0, 0.72, 1.4]]), AIR_VISCOSITY
    )[0]
    assert_allclose(
        [result.cf * root, result.delta99_over_x * root, result.momentum_over_x * root],
        [expected[0], expected[1], expected[3]],
        rtol=1e-3,
    )


def test_van_driest_wall_temperature_stays_within_sa_reference_bands() -> None:
    for row in reference_rows():
        temperature = float(row["edge_temperature_K"])
        reynolds = float(row["reynolds_number"])
        mach = float(row["mach"])
        viscosity = float(AIR_VISCOSITY.dynamic_viscosity(temperature))
        velocity = mach * np.sqrt(
            AIR.heat_capacity_ratio * AIR.specific_gas_constant * temperature
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = flat_plate_boundary_layer(
                1.0,
                velocity,
                reynolds * viscosity / velocity,
                viscosity,
                regime=BoundaryLayerRegime.TURBULENT,
                compressible=True,
                turbulent_model=TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE,
                mach=mach,
                edge_temperature=temperature,
                wall_temperature=float(row["wall_temperature_K"]),
            )
        friction = float(result.local_skin_friction_coefficient) / float(row["cf"])
        thickness = float(result.boundary_layer_thickness) / float(
            row["delta99_over_x"]
        )
        label = f"{row['case']} Re_x={reynolds:g}"
        # Fully turbulent VD2 versus SA developing from the leading edge.
        assert abs(friction - 1) < (0.03 if reynolds >= 1e7 else 0.08), label
        assert abs(thickness - 1) < 0.12, label


@pytest.mark.generated_assets
def test_sa_reference_regenerates() -> None:
    generated = reference_solutions()
    stored = reference_rows()
    assert len(generated) == len(stored) == 18
    for actual, expected in zip(generated, stored, strict=True):
        assert actual["case"] == expected["case"]
        for key in ("reynolds_number", "wall_temperature_K", *QUANTITIES):
            assert_allclose(
                float(actual[key]),
                float(expected[key]),
                rtol=1e-6,
                err_msg=f"{REFERENCE}: {expected['case']}, {key}",
            )
