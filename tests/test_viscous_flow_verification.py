"""Source-equation and integration checks for viscous-flow verification."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    CompressibilityCorrection,
    TurbulentCorrelation,
    flat_plate_boundary_layer,
)
from aerophysics.gas import AIR
from aerophysics.protrusion import protrusion_drag
from aerophysics.transport import AIR_VISCOSITY

REFERENCE = Path(__file__).parent / "reference_data/viscous_flow"


def _rows() -> list[dict[str, str]]:
    with (REFERENCE / "flat_plate_source_equations.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        return list(csv.DictReader(stream))


def test_flat_plate_correlations_match_independent_source_equations() -> None:
    for row in _rows():
        reynolds = float(row["reynolds_number"])
        laminar = flat_plate_boundary_layer(
            1.0,
            reynolds,
            1.0,
            1.0,
            regime=BoundaryLayerRegime.LAMINAR,
        )
        power = flat_plate_boundary_layer(
            1.0,
            reynolds,
            1.0,
            1.0,
            regime=BoundaryLayerRegime.TURBULENT,
            turbulent_correlation=TurbulentCorrelation.POWER_LAW,
        )
        schlichting = flat_plate_boundary_layer(
            1.0,
            reynolds,
            1.0,
            1.0,
            regime=BoundaryLayerRegime.TURBULENT,
            turbulent_correlation=TurbulentCorrelation.SCHLICHTING,
        )
        comparisons = (
            (laminar.boundary_layer_thickness, "blasius_delta_over_x"),
            (laminar.displacement_thickness, "blasius_displacement_over_x"),
            (laminar.momentum_thickness, "blasius_momentum_over_x"),
            (laminar.local_skin_friction_coefficient, "blasius_local_cf"),
            (laminar.average_skin_friction_coefficient, "blasius_average_cf"),
            (power.boundary_layer_thickness, "power_delta_over_x"),
            (power.local_skin_friction_coefficient, "power_local_cf"),
            (power.average_skin_friction_coefficient, "power_average_cf"),
            (
                schlichting.average_skin_friction_coefficient,
                "schlichting_average_cf",
            ),
        )
        for actual, column in comparisons:
            expected = float(row[column])
            assert abs(float(actual) / expected - 1.0) <= 1.0e-12


def test_viscous_reference_provenance_records_equation_role() -> None:
    metadata = json.loads(
        (REFERENCE / "flat_plate_source_equations.json").read_text(encoding="utf-8")
    )
    assert len(metadata["sources"]) == 3
    assert metadata["role"].startswith("independent direct evaluation")


def test_van_driest_ii_matches_digitized_nasa_chart_points() -> None:
    with (REFERENCE / "van_driest_ii_chart.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        temperature = float(row["edge_temperature_K"])
        mach = float(row["mach"])
        reynolds = float(row["reynolds_number"])
        viscosity = float(AIR_VISCOSITY.dynamic_viscosity(temperature))
        recovery = temperature * (
            1.0 + 0.72 ** (1.0 / 3.0) * (AIR.heat_capacity_ratio - 1.0) * mach**2 / 2.0
        )
        result = flat_plate_boundary_layer(
            1.0,
            reynolds * viscosity,
            1.0,
            viscosity,
            regime=BoundaryLayerRegime.TURBULENT,
            compressibility_correction=CompressibilityCorrection.VAN_DRIEST_II,
            mach=mach,
            edge_temperature=temperature,
            wall_temperature=(
                float(row["wall_to_adiabatic_temperature_ratio"]) * recovery
            ),
        )
        difference = abs(
            float(result.local_skin_friction_coefficient) - float(row["local_cf"])
        )
        assert difference <= float(row["absolute_tolerance"])

    metadata = json.loads(
        (REFERENCE / "van_driest_ii_chart.json").read_text(encoding="utf-8")
    )
    assert metadata["role"] == "chart-resolution acceptance comparison"


def test_one_seventh_power_protrusion_matches_closed_form_integral() -> None:
    for ratio in (0.01, 0.03, 0.1, 0.3, 0.8):
        result = protrusion_drag(
            1.0,
            ratio,
            1.0,
            1.0,
            1.0,
            1.0,
            integration_points=32_769,
        )
        expected = 7.0 / 9.0 * ratio ** (2.0 / 7.0)
        assert abs(result.shielding_factor / expected - 1.0) <= 1.0e-6


def test_protrusion_grid_and_geometric_limits() -> None:
    coarse = protrusion_drag(1.0, 0.3, 1.0, 1.0, 1.0, 1.0, integration_points=4097)
    fine = protrusion_drag(1.0, 0.3, 1.0, 1.0, 1.0, 1.0, integration_points=32_769)
    assert abs(coarse.shielding_factor / fine.shielding_factor - 1.0) <= 1.0e-5

    vanishing = protrusion_drag(
        1.0, 1.0e-12, 1.0, 1.0, 1.0, 1.0, integration_points=4097
    )
    outside = protrusion_drag(1.0, 1.0e4, 1.0, 1.0, 1.0, 1.0, integration_points=4097)
    assert vanishing.shielding_factor < 1.0e-3
    assert abs(outside.shielding_factor - 1.0) <= 1.0e-4
