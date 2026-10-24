"""Offline source and invariant checks for compressible-flow verification."""

from __future__ import annotations

import ast
import csv
import json
import math
import runpy
from pathlib import Path
from typing import Any

from aerophysics.expansion import prandtl_meyer_angle
from aerophysics.isentropic import area_ratio, isentropic_ratios
from aerophysics.shocks import ShockBranch, conical_shock, normal_shock, oblique_shock

REFERENCE = Path(__file__).parent / "reference_data/compressible_flow"
SP3004_CAPTURE = (
    Path(__file__).resolve().parents[1] / "docs/scripts/capture_sp3004_reference.py"
)


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


def _sp3004_capture() -> dict[str, Any]:
    return runpy.run_path(str(SP3004_CAPTURE))


def test_sp3004_capture_never_imports_the_solver() -> None:
    tree = ast.parse(SP3004_CAPTURE.read_text(encoding="utf-8"))
    modules = {
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    } | {
        node.module or "" for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)
    }
    assert not any(module.startswith("aerophysics") for module in modules)


def test_sp3004_capture_rejects_undecodable_and_inconsistent_cells() -> None:
    capture = _sp3004_capture()
    cell = capture["_Cell"]
    table = capture["_Table"](
        2.5,
        2.0,
        21,
        35,
        {
            "shock_angle": cell(0.0, 0.0, "+52369890+00"),
            "surface_mstar": cell(0.0, 0.0, "+1624256?+01"),
            "post_shock_mstar": cell(0.0, 0.0, "+16328751+01"),
            "pressure": cell(0.0, 0.0, "+10302014+01"),
            "density": cell(0.0, 0.0, "+10214806+01"),
            "temperature": None,
        },
        2.0,
        0.5236989,
    )
    errors: list[str] = []
    values = capture["_resolve"](table, errors)
    # A damaged or missing cell is never filled from any other source.
    assert "surface_mstar" not in values
    assert "temperature" not in values
    assert len(errors) == 2
    assert all("table 21, PDF pages 35-36" in error for error in errors)

    errors = []
    values.update(surface_mstar=1.6242561, temperature=1.0185375)
    capture["_check"](table, values, errors)
    assert any("p/(rho T)" in error for error in errors)
    assert any("total-enthalpy" in error for error in errors)


def test_sp3004_fixture_records_source_locations_and_transcriptions() -> None:
    capture = _sp3004_capture()
    rows = _rows("nasa_sp3004_conical.csv")
    tables = capture["SOURCE_TABLES"]
    assert [
        (float(row["cone_half_angle_deg"]), float(row["upstream_mach"])) for row in rows
    ] == [(cone, mach) for cone, mach, _, _ in tables]
    for row, (_, _, table, page) in zip(rows, tables, strict=True):
        assert int(row["table"]) == table
        assert int(row["printed_page"]) == page - capture["PRINTED_PAGE_OFFSET"]
    metadata = json.loads(
        (REFERENCE / "nasa_sp3004_conical.json").read_text(encoding="utf-8")
    )
    recorded = {
        (entry["table"], entry["field"]): entry["value"]
        for entry in metadata["manual_cells"]
    }
    assert recorded == {
        key: value for key, (value, _) in capture["MANUAL_CELLS"].items()
    }
