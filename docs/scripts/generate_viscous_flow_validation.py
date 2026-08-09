"""Generate boundary-layer and protrusion verification artifacts."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from _verification_common import line_chart_svg, write_or_check

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "tests/reference_data/viscous_flow"
TABLE_PATH = ROOT / "docs/_generated/viscous_flow_validation.rst"
FRICTION_PATH = ROOT / "docs/_static/viscous_skin_friction.svg"
PROTRUSION_PATH = ROOT / "docs/_static/protrusion_shielding.svg"


def _api() -> tuple[object, ...]:
    sys.path.insert(0, str(ROOT / "src"))
    from aerophysics.boundary_layer import (
        BoundaryLayerRegime,
        CompressibilityCorrection,
        TurbulentCorrelation,
        flat_plate_boundary_layer,
    )
    from aerophysics.boundary_layer_profile import (
        CompressibleVelocityTransformation,
        transform_compressible_velocity_profile,
    )
    from aerophysics.gas import AIR
    from aerophysics.protrusion import protrusion_drag
    from aerophysics.transport import AIR_VISCOSITY

    return (
        BoundaryLayerRegime,
        CompressibilityCorrection,
        TurbulentCorrelation,
        CompressibleVelocityTransformation,
        flat_plate_boundary_layer,
        transform_compressible_velocity_profile,
        protrusion_drag,
        AIR_VISCOSITY,
        AIR,
    )


def _rows() -> list[dict[str, str]]:
    with (REFERENCE / "flat_plate_source_equations.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        return list(csv.DictReader(stream))


def _source_comparison(api: tuple[object, ...]) -> dict[str, tuple[float, float]]:
    regime_enum, _, correlation_enum, _, flat_plate, _, _, _, _ = api
    maxima = {
        name: (0.0, 0.0)
        for name in (
            "Blasius delta99",
            "Blasius displacement thickness",
            "Blasius momentum thickness",
            "Blasius local Cf",
            "Blasius average Cf",
            "power-law delta99",
            "power-law local Cf",
            "power-law average Cf",
            "Schlichting average Cf",
        )
    }
    for row in _rows():
        reynolds = float(row["reynolds_number"])
        laminar = flat_plate(
            1.0,
            reynolds,
            1.0,
            1.0,
            regime=regime_enum.LAMINAR,
        )
        power = flat_plate(
            1.0,
            reynolds,
            1.0,
            1.0,
            regime=regime_enum.TURBULENT,
            turbulent_correlation=correlation_enum.POWER_LAW,
        )
        schlichting = flat_plate(
            1.0,
            reynolds,
            1.0,
            1.0,
            regime=regime_enum.TURBULENT,
            turbulent_correlation=correlation_enum.SCHLICHTING,
        )
        comparisons = {
            "Blasius delta99": (
                laminar.boundary_layer_thickness,
                "blasius_delta_over_x",
            ),
            "Blasius displacement thickness": (
                laminar.displacement_thickness,
                "blasius_displacement_over_x",
            ),
            "Blasius momentum thickness": (
                laminar.momentum_thickness,
                "blasius_momentum_over_x",
            ),
            "Blasius local Cf": (
                laminar.local_skin_friction_coefficient,
                "blasius_local_cf",
            ),
            "Blasius average Cf": (
                laminar.average_skin_friction_coefficient,
                "blasius_average_cf",
            ),
            "power-law delta99": (power.boundary_layer_thickness, "power_delta_over_x"),
            "power-law local Cf": (
                power.local_skin_friction_coefficient,
                "power_local_cf",
            ),
            "power-law average Cf": (
                power.average_skin_friction_coefficient,
                "power_average_cf",
            ),
            "Schlichting average Cf": (
                schlichting.average_skin_friction_coefficient,
                "schlichting_average_cf",
            ),
        }
        for label, (actual, column) in comparisons.items():
            expected = float(row[column])
            relative = abs(float(actual) / expected - 1.0)
            if relative > maxima[label][0]:
                maxima[label] = (relative, reynolds)
    return maxima


def _chart_comparison(api: tuple[object, ...]) -> tuple[float, float]:
    regime_enum, correction_enum, _, _, flat_plate, _, _, viscosity_model, gas = api
    maximum = (0.0, 0.0)
    with (REFERENCE / "van_driest_ii_chart.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        for row in csv.DictReader(stream):
            temperature = float(row["edge_temperature_K"])
            mach = float(row["mach"])
            reynolds = float(row["reynolds_number"])
            edge_viscosity = float(viscosity_model.dynamic_viscosity(temperature))
            recovery = temperature * (
                1.0 + np.cbrt(0.72) * (gas.heat_capacity_ratio - 1.0) * mach**2 / 2.0
            )
            result = flat_plate(
                1.0,
                reynolds * edge_viscosity,
                1.0,
                edge_viscosity,
                regime=regime_enum.TURBULENT,
                compressibility_correction=correction_enum.VAN_DRIEST_II,
                mach=mach,
                edge_temperature=temperature,
                wall_temperature=(
                    float(row["wall_to_adiabatic_temperature_ratio"]) * recovery
                ),
            )
            ratio = abs(
                float(result.local_skin_friction_coefficient) - float(row["local_cf"])
            ) / float(row["absolute_tolerance"])
            if ratio > maximum[0]:
                maximum = (ratio, mach)
    return maximum


def _invariants(api: tuple[object, ...]) -> list[tuple[str, str, float, bool]]:
    (
        regime_enum,
        correction_enum,
        _,
        transformation_enum,
        flat_plate,
        transform_profile,
        protrusion_drag,
        viscosity_model,
        _,
    ) = api
    reynolds = np.geomspace(1.0e5, 1.0e9, 300)
    result = flat_plate(
        1.0,
        reynolds,
        1.0,
        1.0,
        regime=regime_enum.TURBULENT,
    )
    average = np.asarray(result.average_skin_friction_coefficient)
    expected_local = average * (1.0 - 2.58 / np.log(reynolds))
    derivative_error = float(
        np.max(
            np.abs(
                np.asarray(result.local_skin_friction_coefficient) / expected_local
                - 1.0
            )
        )
    )
    edge_viscosity = float(viscosity_model.dynamic_viscosity(300.0))
    incompressible = flat_plate(
        1.0,
        100.0,
        1.0,
        edge_viscosity,
        regime=regime_enum.TURBULENT,
        compressibility_correction=correction_enum.VAN_DRIEST_II,
        mach=0.0,
        edge_temperature=300.0,
    )
    low_mach = flat_plate(
        1.0,
        100.0,
        1.0,
        edge_viscosity,
        regime=regime_enum.TURBULENT,
        compressibility_correction=correction_enum.VAN_DRIEST_II,
        mach=1.0e-6,
        edge_temperature=300.0,
    )
    low_mach_error = abs(
        float(low_mach.average_skin_friction_coefficient)
        / float(incompressible.average_skin_friction_coefficient)
        - 1.0
    )
    height = np.linspace(0.0, 0.01, 1001)
    velocity = np.linspace(0.0, 100.0, 1001)
    transformed = transform_profile(
        height,
        velocity,
        np.full_like(height, 1.2),
        np.full_like(height, 1.8e-5),
        12.0,
        transformation=transformation_enum.VAN_DRIEST,
    )
    friction_velocity = np.sqrt(12.0 / 1.2)
    profile_error = max(
        float(
            np.max(
                np.abs(
                    np.asarray(transformed.transformed_wall_coordinate)
                    - height * 1.2 * friction_velocity / 1.8e-5
                )
            )
        ),
        float(
            np.max(
                np.abs(
                    np.asarray(transformed.transformed_velocity_plus)
                    - velocity / friction_velocity
                )
            )
        ),
    )
    height_ratio = np.asarray([0.01, 0.03, 0.1, 0.3, 0.8])
    shielding_error = 0.0
    for ratio in height_ratio:
        drag = protrusion_drag(
            1.0,
            float(ratio),
            1.0,
            1.0,
            1.0,
            1.0,
            integration_points=32_769,
        )
        expected = 7.0 / 9.0 * ratio ** (2.0 / 7.0)
        shielding_error = max(
            shielding_error, abs(float(drag.shielding_factor) / expected - 1.0)
        )
    constant = protrusion_drag(
        1.2,
        0.02,
        0.01,
        100.0,
        1.2,
        0.1,
        profile_height=[0.0, 0.02],
        profile_velocity=[100.0, 100.0],
        profile_density=[1.2, 1.2],
    )
    expected_force = 1.2 * 0.5 * 1.2 * 100.0**2 * 0.02 * 0.01
    constant_error = abs(float(constant.direct_drag) / expected_force - 1.0)
    coarse = protrusion_drag(1.0, 0.3, 1.0, 1.0, 1.0, 1.0, integration_points=4097)
    fine = protrusion_drag(1.0, 0.3, 1.0, 1.0, 1.0, 1.0, integration_points=32_769)
    grid_error = abs(
        float(coarse.shielding_factor) / float(fine.shielding_factor) - 1.0
    )
    vanishing = protrusion_drag(
        1.0, 1.0e-12, 1.0, 1.0, 1.0, 1.0, integration_points=4097
    )
    outside = protrusion_drag(1.0, 1.0e4, 1.0, 1.0, 1.0, 1.0, integration_points=4097)
    return [
        (
            "Schlichting local/average derivative",
            "relative <= 1e-12",
            derivative_error,
            derivative_error <= 1.0e-12,
        ),
        (
            "Van Driest II low-Mach limit",
            "relative <= 1e-10",
            low_mach_error,
            low_mach_error <= 1.0e-10,
        ),
        (
            "Van Driest constant-property transform",
            "absolute <= 1e-10",
            profile_error,
            profile_error <= 1.0e-10,
        ),
        (
            "One-seventh-power protrusion integral",
            "relative <= 1e-6",
            shielding_error,
            shielding_error <= 1.0e-6,
        ),
        (
            "Constant-profile protrusion integral",
            "relative <= 1e-12",
            constant_error,
            constant_error <= 1.0e-12,
        ),
        (
            "Protrusion integration grid convergence",
            "relative <= 1e-5",
            grid_error,
            grid_error <= 1.0e-5,
        ),
        (
            "Protrusion zero-height limit",
            "shielding < 1e-3",
            float(vanishing.shielding_factor),
            float(vanishing.shielding_factor) < 1.0e-3,
        ),
        (
            "Protrusion outside-layer limit",
            "absolute distance from 1 <= 1e-4",
            abs(float(outside.shielding_factor) - 1.0),
            abs(float(outside.shielding_factor) - 1.0) <= 1.0e-4,
        ),
    ]


def _rst(
    source: dict[str, tuple[float, float]],
    chart: tuple[float, float],
    invariants: list[tuple[str, str, float, bool]],
) -> str:
    passed = (
        all(value[0] <= 1.0e-12 for value in source.values())
        and chart[0] <= 1.0
        and all(item[3] for item in invariants)
    )
    status = "Verified" if passed else "Needs revision"
    lines = [
        f"**Overall status: {status}.**",
        "",
        "Published-correlation comparison",
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
        "",
        (
            ".. list-table:: Maximum relative differences from independently "
            "evaluated source equations"
        ),
        "   :header-rows: 1",
        "",
        "   * - Quantity",
        "     - Maximum relative difference",
        "     - Reynolds number",
        "     - Result",
    ]
    for label, (difference, reynolds) in source.items():
        lines.extend(
            [
                f"   * - {label}",
                f"     - {difference:.4g}",
                f"     - {reynolds:.4g}",
                f"     - {'Pass' if difference <= 1.0e-12 else 'Fail'}",
            ]
        )
    lines.extend(
        [
            "",
            "NASA TN D-6945 chart comparison",
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
            "",
            (
                "For representative Figure 3(e) points, the largest absolute "
                f"difference/tolerance ratio is ``{chart[0]:.4g}`` at Mach "
                f"``{chart[1]:g}`` (``{'Pass' if chart[0] <= 1.0 else 'Fail'}``)."
            ),
        ]
    )
    lines.extend(
        [
            "",
            "Physical and numerical invariants",
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
            "",
            ".. list-table:: Profile, limit, and integration checks",
            "   :header-rows: 1",
            "",
            "   * - Check",
            "     - Criterion",
            "     - Maximum difference",
            "     - Result",
        ]
    )
    for label, criterion, value, item_passed in invariants:
        lines.extend(
            [
                f"   * - {label}",
                f"     - {criterion}",
                f"     - {value:.4g}",
                f"     - {'Pass' if item_passed else 'Fail'}",
            ]
        )
    return "\n".join(lines)


def generate(*, check: bool) -> bool:
    api = _api()
    source = _source_comparison(api)
    chart = _chart_comparison(api)
    invariants = _invariants(api)
    reynolds = np.geomspace(1.0e5, 1.0e9, 400)
    laminar = 0.664 / np.sqrt(reynolds)
    power = 0.0592 * reynolds**-0.2
    schlichting_average = 0.455 / np.log10(reynolds) ** 2.58
    schlichting = schlichting_average * (1.0 - 2.58 / np.log(reynolds))
    friction = line_chart_svg(
        title="Smooth flat-plate local skin friction",
        description=(
            "Blasius, one-fifth-power, and Schlichting local skin-friction "
            "correlations over Reynolds number."
        ),
        x_label="Reynolds number [-]",
        y_label="Local Cf [-]",
        x=reynolds,
        series=(
            ("Blasius", "#0072B2", laminar),
            ("Power law", "#D55E00", power),
            ("Schlichting", "#009E73", schlichting),
        ),
        log_x=True,
        log_y=True,
    )
    ratio = np.geomspace(1.0e-3, 1.0, 250)
    shielding = 7.0 / 9.0 * ratio ** (2.0 / 7.0)
    protrusion = line_chart_svg(
        title="Boundary-layer shielding of a protrusion",
        description=(
            "Closed-form effective dynamic-pressure factor for a constant-width "
            "protrusion in a one-seventh-power velocity profile."
        ),
        x_label="Protrusion height / boundary-layer thickness [-]",
        y_label="Shielding factor [-]",
        x=ratio,
        series=(("Analytic 7/9 (h/delta)^(2/7)", "#CC79A7", shielding),),
        log_x=True,
    )
    current = True
    current &= write_or_check(TABLE_PATH, _rst(source, chart, invariants), check=check)
    current &= write_or_check(FRICTION_PATH, friction, check=check)
    current &= write_or_check(PROTRUSION_PATH, protrusion, check=check)
    return current


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if not generate(check=args.check):
        raise SystemExit("viscous-flow verification artifacts are stale")


if __name__ == "__main__":
    main()
