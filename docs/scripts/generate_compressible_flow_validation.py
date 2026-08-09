"""Generate the compressible-flow verification record and SVG figures."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np
from _verification_common import line_chart_svg, write_or_check

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "tests/reference_data/compressible_flow"
TABLE_PATH = ROOT / "docs/_generated/compressible_flow_validation.rst"
PROFILE_PATH = ROOT / "docs/_static/compressible_flow_profiles.svg"
ERROR_PATH = ROOT / "docs/_static/compressible_flow_differences.svg"


def _load_api() -> tuple[object, ...]:
    sys.path.insert(0, str(ROOT / "src"))
    from aerophysics.expansion import (
        mach_from_prandtl_meyer,
        prandtl_meyer_angle,
    )
    from aerophysics.isentropic import (
        MachBranch,
        area_ratio,
        isentropic_ratios,
        mach_from_area_ratio,
        mach_from_total_density_ratio,
        mach_from_total_pressure_ratio,
        mach_from_total_temperature_ratio,
        mass_flow_parameter,
    )
    from aerophysics.shocks import (
        ShockBranch,
        conical_shock,
        maximum_attached_deflection,
        normal_shock,
        oblique_shock,
        shock_angle,
        theta_from_shock_angle,
    )

    return (
        MachBranch,
        ShockBranch,
        area_ratio,
        conical_shock,
        isentropic_ratios,
        mach_from_area_ratio,
        mach_from_prandtl_meyer,
        mach_from_total_density_ratio,
        mach_from_total_pressure_ratio,
        mach_from_total_temperature_ratio,
        mass_flow_parameter,
        maximum_attached_deflection,
        normal_shock,
        oblique_shock,
        prandtl_meyer_angle,
        shock_angle,
        theta_from_shock_angle,
    )


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _maximum_naca_errors(api: tuple[object, ...]) -> dict[str, tuple[float, float]]:
    (
        _,
        _,
        area_ratio,
        _,
        isentropic_ratios,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        normal_shock,
        _,
        prandtl_meyer_angle,
        _,
        _,
    ) = api
    specs = {
        "pressure_over_total": ("pressure_printed_unit", "pressure"),
        "density_over_total": ("density_printed_unit", "density"),
        "temperature_over_total": ("temperature_printed_unit", "temperature"),
        "area_over_critical": ("area_printed_unit", "area"),
        "prandtl_meyer_deg": ("prandtl_meyer_printed_unit_deg", "nu"),
        "downstream_mach": ("downstream_mach_printed_unit", "m2"),
        "pressure_ratio": ("pressure_ratio_printed_unit", "p2"),
        "density_ratio": ("density_ratio_printed_unit", "rho2"),
        "temperature_ratio": ("temperature_ratio_printed_unit", "t2"),
        "total_pressure_ratio": ("total_pressure_ratio_printed_unit", "pt2"),
    }
    maxima = {label: (0.0, 0.0) for _, label in specs.values()}
    for row in _rows(REFERENCE / "naca1135_tables_i_ii.csv"):
        mach = float(row["mach"])
        ratios = isentropic_ratios(mach)
        actual = {
            "pressure_over_total": 1.0 / float(ratios.total_pressure_ratio),
            "density_over_total": 1.0 / float(ratios.total_density_ratio),
            "temperature_over_total": 1.0 / float(ratios.total_temperature_ratio),
        }
        if mach > 0.0:
            actual["area_over_critical"] = float(area_ratio(mach))
        if mach >= 1.0:
            shock = normal_shock(mach)
            actual.update(
                {
                    "prandtl_meyer_deg": math.degrees(float(prandtl_meyer_angle(mach))),
                    "downstream_mach": float(shock.downstream_mach),
                    "pressure_ratio": float(shock.static_pressure_ratio),
                    "density_ratio": float(shock.static_density_ratio),
                    "temperature_ratio": float(shock.static_temperature_ratio),
                    "total_pressure_ratio": float(shock.total_pressure_ratio),
                }
            )
        for name, (unit_name, label) in specs.items():
            if not row[name]:
                continue
            expected = float(row[name])
            tolerance = max(float(row[unit_name]), 1.0e-4 * abs(expected))
            error = abs(actual[name] - expected)
            ratio = error / tolerance
            if ratio > maxima[label][0]:
                maxima[label] = (ratio, mach)
    return maxima


def _maximum_cone_error(api: tuple[object, ...]) -> tuple[float, float, float]:
    conical_shock = api[3]
    maximum = 0.0
    location = (0.0, 0.0)
    for row in _rows(REFERENCE / "nasa_sp3004_conical.csv"):
        mach = float(row["upstream_mach"])
        cone_angle = float(row["cone_half_angle_deg"])
        result = conical_shock(mach, math.radians(cone_angle))
        for name in (
            "shock_angle",
            "surface_mach",
            "post_shock_mach",
            "surface_pressure_ratio",
            "surface_density_ratio",
            "surface_temperature_ratio",
        ):
            reference_name = "shock_angle_rad" if name == "shock_angle" else name
            expected = float(row[reference_name])
            relative = abs(float(getattr(result, name)) - expected) / abs(expected)
            if relative > maximum:
                maximum = relative
                location = (mach, cone_angle)
    return maximum, *location


def _maximum_chart_error(api: tuple[object, ...]) -> tuple[float, float, float]:
    shock_branch_enum = api[1]
    oblique_shock = api[13]
    maximum = 0.0
    location = (0.0, 0.0)
    for row in _rows(REFERENCE / "naca1135_charts_2_4.csv"):
        mach = float(row["upstream_mach"])
        deflection = float(row["deflection_angle_deg"])
        result = oblique_shock(
            mach,
            math.radians(deflection),
            shock_branch_enum(row["branch"]),
        )
        pressure_coefficient = (
            2.0 / (1.4 * mach**2) * (float(result.static_pressure_ratio) - 1.0)
        )
        ratios = (
            abs(math.degrees(float(result.shock_angle)) - float(row["shock_angle_deg"]))
            / float(row["shock_angle_abs_tolerance_deg"]),
            abs(pressure_coefficient - float(row["pressure_coefficient"]))
            / float(row["pressure_coefficient_abs_tolerance"]),
            abs(float(result.downstream_mach) - float(row["downstream_mach"]))
            / float(row["downstream_mach_abs_tolerance"]),
        )
        ratio = max(ratios)
        if ratio > maximum:
            maximum = ratio
            location = (mach, deflection)
    return maximum, *location


def _invariants(api: tuple[object, ...]) -> list[tuple[str, str, float, bool]]:
    (
        mach_branch_enum,
        shock_branch_enum,
        area_ratio,
        conical_shock,
        isentropic_ratios,
        mach_from_area_ratio,
        mach_from_prandtl_meyer,
        mach_from_density,
        mach_from_pressure,
        mach_from_temperature,
        mass_flow_parameter,
        maximum_attached_deflection,
        normal_shock,
        oblique_shock,
        prandtl_meyer_angle,
        shock_angle,
        theta_from_shock_angle,
    ) = api
    mach = np.concatenate(
        (np.linspace(0.01, 0.99, 300), [1.0], np.linspace(1.01, 10.0, 700))
    )
    ratios = isentropic_ratios(mach)
    inverse_error = max(
        float(
            np.max(
                np.abs(
                    np.asarray(mach_from_temperature(ratios.total_temperature_ratio))
                    - mach
                )
                / mach
            )
        ),
        float(
            np.max(
                np.abs(
                    np.asarray(mach_from_pressure(ratios.total_pressure_ratio)) - mach
                )
                / mach
            )
        ),
        float(
            np.max(
                np.abs(np.asarray(mach_from_density(ratios.total_density_ratio)) - mach)
                / mach
            )
        ),
    )
    area_targets = np.asarray(area_ratio(mach))
    subsonic = mach < 1.0
    supersonic = mach > 1.0
    area_error = max(
        float(
            np.max(
                np.abs(
                    np.asarray(
                        mach_from_area_ratio(
                            area_targets[subsonic], mach_branch_enum.SUBSONIC
                        )
                    )
                    - mach[subsonic]
                )
                / mach[subsonic]
            )
        ),
        float(
            np.max(
                np.abs(
                    np.asarray(
                        mach_from_area_ratio(
                            area_targets[supersonic], mach_branch_enum.SUPERSONIC
                        )
                    )
                    - mach[supersonic]
                )
                / mach[supersonic]
            )
        ),
    )
    nu_mach = np.linspace(1.0, 20.0, 800)
    nu = np.asarray(prandtl_meyer_angle(nu_mach))
    nu_error = float(
        np.max(np.abs(np.asarray(mach_from_prandtl_meyer(nu)) - nu_mach) / nu_mach)
    )
    flux_at_one = float(mass_flow_parameter(1.0))
    flux_margin = min(
        flux_at_one - float(mass_flow_parameter(0.999)),
        flux_at_one - float(mass_flow_parameter(1.001)),
    )

    shock_mach = np.linspace(1.001, 10.0, 900)
    shock = normal_shock(shock_mach)
    m2 = np.asarray(shock.downstream_mach)
    p2 = np.asarray(shock.static_pressure_ratio)
    rho2 = np.asarray(shock.static_density_ratio)
    t2 = np.asarray(shock.static_temperature_ratio)
    mass_residual = float(np.max(np.abs(rho2 * m2 * np.sqrt(t2) / shock_mach - 1.0)))
    momentum_left = 1.0 + 1.4 * shock_mach**2
    momentum_right = p2 + 1.4 * rho2 * (m2 * np.sqrt(t2)) ** 2
    momentum_residual = float(np.max(np.abs(momentum_right / momentum_left - 1.0)))
    energy_left = 3.5 + 0.5 * 1.4 * shock_mach**2
    energy_right = 3.5 * t2 + 0.5 * 1.4 * (m2 * np.sqrt(t2)) ** 2
    energy_residual = float(np.max(np.abs(energy_right / energy_left - 1.0)))

    oblique_error = 0.0
    branch_error = 0.0
    for upstream in (1.5, 2.0, 3.0, 5.0, 10.0):
        limit = maximum_attached_deflection(upstream)
        theta = 0.5 * float(limit.deflection_angle)
        weak = oblique_shock(upstream, theta, shock_branch_enum.WEAK)
        normal = normal_shock(upstream * math.sin(float(weak.shock_angle)))
        oblique_error = max(
            oblique_error,
            abs(
                float(weak.static_pressure_ratio) / float(normal.static_pressure_ratio)
                - 1.0
            ),
            abs(
                float(theta_from_shock_angle(upstream, weak.shock_angle)) / theta - 1.0
            ),
        )
        weak_limit = float(
            shock_angle(upstream, limit.deflection_angle, shock_branch_enum.WEAK)
        )
        strong_limit = float(
            shock_angle(upstream, limit.deflection_angle, shock_branch_enum.STRONG)
        )
        branch_error = max(branch_error, abs(weak_limit - strong_limit))
    cone_zero = max(
        abs(float(conical_shock(value, 0.0).surface_pressure_ratio) - 1.0)
        for value in (1.5, 2.0, 3.0, 5.0)
    )
    return [
        (
            "Isentropic ratio inverse round trips",
            "relative <= 1e-10",
            inverse_error,
            inverse_error <= 1.0e-10,
        ),
        (
            "Area-Mach branch round trips",
            "relative <= 1e-10",
            area_error,
            area_error <= 1.0e-10,
        ),
        (
            "Prandtl-Meyer inverse round trip",
            "relative <= 1e-10",
            nu_error,
            nu_error <= 1.0e-10,
        ),
        (
            "Mass-flow maximum at M=1",
            "positive neighbour margin",
            flux_margin,
            flux_margin > 0.0,
        ),
        (
            "Normal-shock mass closure",
            "relative <= 1e-12",
            mass_residual,
            mass_residual <= 1.0e-12,
        ),
        (
            "Normal-shock momentum closure",
            "relative <= 1e-12",
            momentum_residual,
            momentum_residual <= 1.0e-12,
        ),
        (
            "Normal-shock energy closure",
            "relative <= 1e-12",
            energy_residual,
            energy_residual <= 1.0e-12,
        ),
        (
            "Oblique normal-component closure",
            "relative <= 1e-12",
            oblique_error,
            oblique_error <= 1.0e-12,
        ),
        (
            "Weak/strong attached-limit merger",
            "absolute <= 1e-10 rad",
            branch_error,
            branch_error <= 1.0e-10,
        ),
        (
            "Zero-angle cone Mach-wave limit",
            "absolute <= 1e-12",
            cone_zero,
            cone_zero <= 1.0e-12,
        ),
    ]


def _rst(
    naca: dict[str, tuple[float, float]],
    chart: tuple[float, float, float],
    cone: tuple[float, float, float],
    invariants: list[tuple[str, str, float, bool]],
) -> str:
    accepted = all(value[0] <= 1.0 for value in naca.values())
    passed = (
        accepted
        and chart[0] <= 1.0
        and cone[0] <= 1.0e-4
        and all(item[3] for item in invariants)
    )
    status = "Verified" if passed else "Needs revision"
    lines = [
        f"**Overall status: {status}.**",
        "",
        "NACA Report 1135 table comparison",
        "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
        "",
        ".. list-table:: Maximum difference divided by the adopted table tolerance",
        "   :header-rows: 1",
        "",
        "   * - Quantity",
        "     - Maximum ratio",
        "     - Mach",
        "     - Result",
    ]
    for label, (ratio, mach) in naca.items():
        lines.extend(
            [
                f"   * - {label}",
                f"     - {ratio:.4g}",
                f"     - {mach:g}",
                f"     - {'Pass' if ratio <= 1.0 else 'Fail'}",
            ]
        )
    lines.extend(
        [
            "",
            "NACA Charts 2--4 observation",
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
            "",
            f"The largest difference/tolerance ratio is ``{chart[0]:.4g}`` at "
            f"Mach ``{chart[1]:g}``, deflection ``{chart[2]:g} deg`` "
            f"(``{'Pass' if chart[0] <= 1.0 else 'Fail'}``).",
            "",
            "NASA SP-3004 cone-table comparison",
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
            "",
            f"The maximum relative difference is ``{cone[0]:.4g}`` at Mach "
            f"``{cone[1]:g}``, cone half-angle ``{cone[2]:g} deg`` "
            f"(``{'Pass' if cone[0] <= 1.0e-4 else 'Fail'}``).",
            "",
            "Physical and mathematical invariants",
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~",
            "",
            ".. list-table:: Dense-grid invariant checks",
            "   :header-rows: 1",
            "",
            "   * - Check",
            "     - Criterion",
            "     - Maximum or margin",
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
    api = _load_api()
    naca = _maximum_naca_errors(api)
    chart = _maximum_chart_error(api)
    cone = _maximum_cone_error(api)
    invariants = _invariants(api)
    mach = np.geomspace(0.05, 10.0, 320)
    ratios = api[4](mach)
    profile_svg = line_chart_svg(
        title="Calorically perfect-gas isentropic profiles",
        description=(
            "Static-to-total pressure, density, and temperature ratios and "
            "static-to-stagnation state ratios as functions of Mach number."
        ),
        x_label="Mach number [-]",
        y_label="Dimensionless ratio [-]",
        x=mach,
        series=(
            ("p/pt", "#0072B2", 1.0 / np.asarray(ratios.total_pressure_ratio)),
            ("rho/rhot", "#009E73", 1.0 / np.asarray(ratios.total_density_ratio)),
            ("T/Tt", "#D55E00", 1.0 / np.asarray(ratios.total_temperature_ratio)),
        ),
        log_x=True,
    )
    error_mach = np.asarray([value[1] for value in naca.values()], dtype=np.float64)
    error_ratio = np.asarray(
        [max(value[0], 1.0e-12) for value in naca.values()], dtype=np.float64
    )
    order = np.argsort(error_mach)
    error_svg = line_chart_svg(
        title="NACA Report 1135 normalized differences",
        description=(
            "Maximum difference at each controlling Mach number divided by "
            "the adopted printed-table tolerance."
        ),
        x_label="Mach number [-]",
        y_label="Difference / tolerance [-]",
        x=error_mach[order],
        series=(("Table comparison", "#CC79A7", error_ratio[order]),),
        log_x=True,
        log_y=True,
    )
    current = True
    current &= write_or_check(
        TABLE_PATH, _rst(naca, chart, cone, invariants), check=check
    )
    current &= write_or_check(PROFILE_PATH, profile_svg, check=check)
    current &= write_or_check(ERROR_PATH, error_svg, check=check)
    return current


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if not generate(check=args.check):
        raise SystemExit("compressible-flow verification artifacts are stale")


if __name__ == "__main__":
    main()
