"""Offline source and invariant checks for compressible-flow verification."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from aerophysics.expansion import prandtl_meyer_angle
from aerophysics.isentropic import area_ratio, isentropic_ratios
from aerophysics.shocks import ShockBranch, conical_shock, normal_shock, oblique_shock

REFERENCE = Path(__file__).parent / "reference_data/compressible_flow"


def _rows(name: str) -> list[dict[str, str]]:
    with (REFERENCE / name).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def test_complete_naca_1135_table_grid_matches_public_api() -> None:
    rows = _rows("naca1135_tables_i_ii.csv")
    assert len(rows) == 1744
    assert rows[0]["mach"] == "0.00"
    assert rows[-1]["mach"] == "100.00"
    tolerance_columns = {
        "pressure_over_total": "pressure_printed_unit",
        "density_over_total": "density_printed_unit",
        "temperature_over_total": "temperature_printed_unit",
        "area_over_critical": "area_printed_unit",
        "prandtl_meyer_deg": "prandtl_meyer_printed_unit_deg",
        "downstream_mach": "downstream_mach_printed_unit",
        "pressure_ratio": "pressure_ratio_printed_unit",
        "density_ratio": "density_ratio_printed_unit",
        "temperature_ratio": "temperature_ratio_printed_unit",
        "total_pressure_ratio": "total_pressure_ratio_printed_unit",
    }
    for row in rows:
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
        for name, value in actual.items():
            if not row[name]:
                continue
            expected = float(row[name])
            unit_name = tolerance_columns[name]
            tolerance = max(float(row[unit_name]), 1.0e-4 * abs(expected))
            assert abs(value - expected) <= tolerance, (mach, name)


def test_naca_fixture_provenance_is_complete() -> None:
    metadata = json.loads(
        (REFERENCE / "naca1135_tables_i_ii.json").read_text(encoding="utf-8")
    )
    assert metadata["source"] == "NACA Report 1135, Tables I and II"
    assert metadata["printed_pages"] == "633-651"
    assert metadata["gamma"] == 1.4
    assert metadata["row_count"] == 1744


def test_naca_oblique_shock_charts_agree_at_chart_resolution() -> None:
    rows = _rows("naca1135_charts_2_4.csv")
    assert len(rows) == 5
    for row in rows:
        mach = float(row["upstream_mach"])
        result = oblique_shock(
            mach,
            math.radians(float(row["deflection_angle_deg"])),
            ShockBranch(row["branch"]),
        )
        pressure_coefficient = (
            2.0 / (1.4 * mach**2) * (float(result.static_pressure_ratio) - 1.0)
        )
        assert abs(
            math.degrees(float(result.shock_angle)) - float(row["shock_angle_deg"])
        ) <= float(row["shock_angle_abs_tolerance_deg"])
        assert abs(pressure_coefficient - float(row["pressure_coefficient"])) <= float(
            row["pressure_coefficient_abs_tolerance"]
        )
        assert abs(
            float(result.downstream_mach) - float(row["downstream_mach"])
        ) <= float(row["downstream_mach_abs_tolerance"])


def test_nasa_sp3004_conical_table_cells_match_public_api() -> None:
    rows = _rows("nasa_sp3004_conical.csv")
    assert len(rows) == 46
    for row in rows:
        result = conical_shock(
            float(row["upstream_mach"]),
            math.radians(float(row["cone_half_angle_deg"])),
        )
        for field in (
            "surface_mach",
            "post_shock_mach",
            "surface_pressure_ratio",
            "surface_density_ratio",
            "surface_temperature_ratio",
        ):
            expected = float(row[field])
            relative = abs(float(getattr(result, field)) / expected - 1.0)
            assert relative <= float(row["relative_tolerance"]), (
                row["table"],
                field,
            )
        expected_angle = float(row["shock_angle_rad"])
        assert abs(float(result.shock_angle) / expected_angle - 1.0) <= float(
            row["relative_tolerance"]
        )


def test_sp3004_fixture_provenance_records_source_domain() -> None:
    metadata = json.loads(
        (REFERENCE / "nasa_sp3004_conical.json").read_text(encoding="utf-8")
    )
    assert metadata["source"].startswith("NASA SP-3004")
    assert metadata["gamma"] == 1.4
    assert metadata["cone_angles_deg"] == [2.5 * value for value in range(1, 13)]
    assert metadata["mach_numbers"] == [1.5, 2.0, 3.0, 5.0]
