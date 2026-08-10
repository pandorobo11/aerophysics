"""Generate standard-atmosphere verification tables and SVG figures."""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from html import escape
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import numpy as np
from _verification_common import format_bounded_error

if TYPE_CHECKING:
    from collections.abc import Sequence

    from numpy.typing import NDArray

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OFFICIAL_CSV = (
    PROJECT_ROOT / "tests/reference_data/standard_atmosphere/official_1976.csv"
)
FLUIDS_CSV = PROJECT_ROOT / "tests/reference_data/standard_atmosphere/fluids-1.3.1.csv"
TABLE_PATH = PROJECT_ROOT / "docs/_generated/standard_atmosphere_validation.rst"
PROFILE_SVG_PATH = PROJECT_ROOT / "docs/_static/standard_atmosphere_profiles.svg"
COMPARISON_SVG_PATH = PROJECT_ROOT / "docs/_static/standard_atmosphere_comparison.svg"


@dataclass(frozen=True)
class PropertySpec:
    """Mapping between one reference column and an atmosphere-state field."""

    reference_column: str
    tolerance_column: str
    state_field: str
    label: str
    unit: str
    color: str


@dataclass(frozen=True)
class OfficialResult:
    """One comparison against a printed official-table cell."""

    coordinate_type: str
    altitude_m: float
    spec: PropertySpec
    expected: float
    actual: float
    tolerance: float
    verification_tolerance: float
    error: float
    ratio: float
    passed: bool
    verification_passed: bool


@dataclass(frozen=True)
class ExternalResult:
    """Maximum discrepancy for one pinned ``fluids`` snapshot variable."""

    label: str
    criterion: str
    maximum: float
    altitude_m: float
    limit: float
    passed: bool


@dataclass(frozen=True)
class PhysicalResult:
    """Diagnostic result for one physical or mathematical invariant."""

    label: str
    criterion: str
    maximum: float
    passed: bool
    limit: float | None = None


@dataclass(frozen=True)
class Series:
    """One SVG line series."""

    label: str
    color: str
    values: NDArray[np.float64]


PROPERTY_SPECS = (
    PropertySpec(
        "temperature_K",
        "temperature_abs_tolerance_K",
        "temperature",
        "Temperature",
        "K",
        "#0072B2",
    ),
    PropertySpec(
        "pressure_Pa",
        "pressure_abs_tolerance_Pa",
        "pressure",
        "Pressure",
        "Pa",
        "#D55E00",
    ),
    PropertySpec(
        "density_kg_m3",
        "density_abs_tolerance_kg_m3",
        "density",
        "Density",
        "kg/m³",
        "#009E73",
    ),
    PropertySpec(
        "gravity_m_s2",
        "gravity_abs_tolerance_m_s2",
        "gravity",
        "Gravity",
        "m/s²",
        "#CC79A7",
    ),
    PropertySpec(
        "speed_of_sound_m_s",
        "speed_of_sound_abs_tolerance_m_s",
        "speed_of_sound",
        "Speed of sound",
        "m/s",
        "#56B4E9",
    ),
    PropertySpec(
        "dynamic_viscosity_Pa_s",
        "dynamic_viscosity_abs_tolerance_Pa_s",
        "dynamic_viscosity",
        "Dynamic viscosity",
        "Pa·s",
        "#E69F00",
    ),
    PropertySpec(
        "kinematic_viscosity_m2_s",
        "kinematic_viscosity_abs_tolerance_m2_s",
        "kinematic_viscosity",
        "Kinematic viscosity",
        "m²/s",
        "#000000",
    ),
)

EXTERNAL_SPECS = (
    ("temperature_K", "temperature", "Temperature", "absolute", 1.0e-4),
    ("pressure_Pa", "pressure", "Pressure", "relative", 2.0e-5),
    ("density_kg_m3", "density", "Density", "relative", 2.0e-5),
    ("gravity_m_s2", "gravity", "Gravity", "relative", 2.0e-5),
    (
        "speed_of_sound_m_s",
        "speed_of_sound",
        "Speed of sound",
        "relative",
        2.0e-5,
    ),
    (
        "dynamic_viscosity_Pa_s",
        "dynamic_viscosity",
        "Dynamic viscosity",
        "relative",
        2.0e-5,
    ),
    (
        "kinematic_viscosity_m2_s",
        "kinematic_viscosity",
        "Kinematic viscosity",
        "relative",
        2.0e-5,
    ),
    (
        "thermal_conductivity_W_m_K",
        "thermal_conductivity",
        "Thermal conductivity",
        "relative",
        2.0e-5,
    ),
)


def _load_public_api() -> tuple[object, object]:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))
    from aerophysics.atmosphere import (
        geopotential_to_geometric,
        standard_atmosphere,
    )

    return standard_atmosphere, geopotential_to_geometric


def _load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def official_results() -> list[OfficialResult]:
    """Compare public-API values with transcribed official table cells."""
    standard_atmosphere, geopotential_to_geometric = _load_public_api()
    results: list[OfficialResult] = []
    for row in _load_csv(OFFICIAL_CSV):
        altitude = float(row["altitude_m"])
        geometric_altitude = (
            altitude
            if row["coordinate_type"] == "geometric"
            else float(geopotential_to_geometric(altitude))
        )
        state = standard_atmosphere(geometric_altitude)
        for spec in PROPERTY_SPECS:
            if not row[spec.reference_column]:
                continue
            expected = float(row[spec.reference_column])
            tolerance = float(row[spec.tolerance_column])
            actual = float(getattr(state, spec.state_field))
            error = abs(actual - expected)
            guard = 1.0e-14 * max(1.0, abs(expected))
            verification_tolerance = max(
                2.0 * tolerance,
                1.0e-4 * abs(expected),
            )
            results.append(
                OfficialResult(
                    coordinate_type=row["coordinate_type"],
                    altitude_m=altitude,
                    spec=spec,
                    expected=expected,
                    actual=actual,
                    tolerance=tolerance,
                    verification_tolerance=verification_tolerance,
                    error=error,
                    ratio=error / tolerance,
                    passed=error <= tolerance + guard,
                    verification_passed=(error <= verification_tolerance + guard),
                )
            )
    return results


def _external_arrays() -> tuple[NDArray[np.float64], dict[str, NDArray[np.float64]]]:
    rows = _load_csv(FLUIDS_CSV)
    altitude = np.asarray(
        [float(row["geometric_altitude_m"]) for row in rows], dtype=np.float64
    )
    values = {
        column: np.asarray([float(row[column]) for row in rows], dtype=np.float64)
        for column, *_ in EXTERNAL_SPECS
    }
    return altitude, values


def external_results() -> list[ExternalResult]:
    """Compare the public API with the pinned ``fluids`` snapshot."""
    standard_atmosphere, _ = _load_public_api()
    altitude, references = _external_arrays()
    state = standard_atmosphere(altitude)
    results: list[ExternalResult] = []
    for column, state_field, label, kind, threshold in EXTERNAL_SPECS:
        reference = references[column]
        actual = np.asarray(getattr(state, state_field), dtype=np.float64)
        difference = np.abs(actual - reference)
        if kind == "relative":
            difference = difference / np.abs(reference)
            criterion = f"relative ≤ {threshold:.0e}"
        else:
            criterion = f"absolute ≤ {threshold:.0e} K"
        index = int(np.argmax(difference))
        maximum = float(difference[index])
        results.append(
            ExternalResult(
                label=label,
                criterion=criterion,
                maximum=maximum,
                altitude_m=float(altitude[index]),
                limit=threshold,
                passed=maximum <= threshold,
            )
        )
    return results


def physical_results() -> list[PhysicalResult]:
    """Evaluate physical and mathematical invariants on the full range."""
    standard_atmosphere, geopotential_to_geometric = _load_public_api()
    altitude = np.arange(-5_000.0, 86_001.0, 1.0, dtype=np.float64)
    state = standard_atmosphere(altitude)
    results: list[PhysicalResult] = []
    for field in ("pressure", "density", "gravity"):
        step = np.diff(np.asarray(getattr(state, field), dtype=np.float64))
        maximum = float(np.max(step))
        results.append(
            PhysicalResult(
                f"Strict decrease: {field}",
                "largest 1 m step < 0",
                maximum,
                maximum < 0.0,
            )
        )
    viscosity_step = np.diff(np.asarray(state.kinematic_viscosity, dtype=np.float64))
    minimum_viscosity_step = float(np.min(viscosity_step))
    results.append(
        PhysicalResult(
            "Strict increase: kinematic viscosity",
            "smallest 1 m step > 0",
            minimum_viscosity_step,
            minimum_viscosity_step > 0.0,
        )
    )

    boundaries = np.asarray(
        [
            0.0,
            11_000.0,
            20_000.0,
            32_000.0,
            47_000.0,
            51_000.0,
            71_000.0,
            84_852.0,
        ]
    )
    expected_slopes = np.asarray([-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0])
    midpoints = (boundaries[:-1] + boundaries[1:]) / 2.0
    offset = 100.0
    lower = np.asarray(geopotential_to_geometric(midpoints - offset))
    upper = np.asarray(geopotential_to_geometric(midpoints + offset))
    lower_temperature = np.asarray(standard_atmosphere(lower).temperature)
    upper_temperature = np.asarray(standard_atmosphere(upper).temperature)
    slopes = (upper_temperature - lower_temperature) / (2.0 * offset) * 1_000.0
    slope_error = float(np.max(np.abs(slopes - expected_slopes)))
    results.append(
        PhysicalResult(
            "Seven layer temperature gradients",
            "maximum absolute error ≤ 2e-12 K/km",
            slope_error,
            slope_error <= 2.0e-12,
            limit=2.0e-12,
        )
    )

    internal_boundaries = boundaries[1:-1]
    continuity_offset = 1.0e-6
    below = standard_atmosphere(
        geopotential_to_geometric(internal_boundaries - continuity_offset)
    )
    above = standard_atmosphere(
        geopotential_to_geometric(internal_boundaries + continuity_offset)
    )
    continuity_fields = (
        "temperature",
        "pressure",
        "density",
        "speed_of_sound",
        "dynamic_viscosity",
        "kinematic_viscosity",
        "thermal_conductivity",
        "prandtl_number",
    )
    continuity_error = 0.0
    for field in continuity_fields:
        lower_values = np.asarray(getattr(below, field), dtype=np.float64)
        upper_values = np.asarray(getattr(above, field), dtype=np.float64)
        relative = np.abs(lower_values - upper_values) / np.maximum(
            np.abs(lower_values), 1.0e-30
        )
        continuity_error = max(continuity_error, float(np.max(relative)))
    results.append(
        PhysicalResult(
            "Layer-boundary continuity",
            "maximum relative jump ≤ 1e-9",
            continuity_error,
            continuity_error <= 1.0e-9,
            limit=1.0e-9,
        )
    )

    pressure = np.asarray(state.pressure)
    density = np.asarray(state.density)
    temperature = np.asarray(state.temperature)
    gas_constant = np.asarray(state.specific_gas_constant)
    gamma = np.asarray(state.heat_capacity_ratio)
    dynamic_viscosity = np.asarray(state.dynamic_viscosity)
    conductivity = np.asarray(state.thermal_conductivity)
    closure_pairs = (
        (pressure, density * gas_constant * temperature),
        (
            np.asarray(state.speed_of_sound) ** 2,
            gamma * gas_constant * temperature,
        ),
        (np.asarray(state.kinematic_viscosity), dynamic_viscosity / density),
        (
            np.asarray(state.prandtl_number),
            dynamic_viscosity * np.asarray(state.cp) / conductivity,
        ),
        (np.asarray(state.cp) - np.asarray(state.cv), gas_constant),
        (gamma, np.asarray(state.cp) / np.asarray(state.cv)),
    )
    closure_error = max(
        float(np.max(np.abs(actual / expected - 1.0)))
        for actual, expected in closure_pairs
    )
    results.append(
        PhysicalResult(
            "Six state-equation closures",
            "maximum relative error ≤ 1e-12",
            closure_error,
            closure_error <= 1.0e-12,
            limit=1.0e-12,
        )
    )

    geopotential = np.asarray(state.geopotential_altitude)
    derivative = np.gradient(np.log(pressure), geopotential, edge_order=2)
    hydrostatic_expected = -9.80665 / (gas_constant * temperature)
    distance = np.min(np.abs(geopotential[:, np.newaxis] - internal_boundaries), axis=1)
    mask = distance > 2.0
    mask[[0, -1]] = False
    hydrostatic_error = float(
        np.max(np.abs(derivative[mask] / hydrostatic_expected[mask] - 1.0))
    )
    results.append(
        PhysicalResult(
            "Numerical hydrostatic derivative",
            "maximum relative error ≤ 1e-7",
            hydrostatic_error,
            hydrostatic_error <= 1.0e-7,
            limit=1.0e-7,
        )
    )

    geometric_round_trip = np.asarray(
        geopotential_to_geometric(geopotential), dtype=np.float64
    )
    round_trip_error = float(np.max(np.abs(geometric_round_trip - altitude)))
    results.append(
        PhysicalResult(
            "Altitude-coordinate round trip",
            "maximum absolute error ≤ 1e-10 m",
            round_trip_error,
            round_trip_error <= 1.0e-10,
            limit=1.0e-10,
        )
    )
    return results


def _number(value: float) -> str:
    if value == 0.0:
        return "0"
    if abs(value) < 1.0e-3 or abs(value) >= 1.0e5:
        return f"{value:.4e}"
    return f"{value:.6g}"


def _official_summary(results: Sequence[OfficialResult]) -> list[str]:
    lines = [
        ".. list-table:: Official-table comparison summary",
        "   :header-rows: 1",
        "   :widths: 20 12 15 20 18 15",
        "",
        "   * - Quantity",
        "     - Cells",
        "     - Within acceptance",
        "     - Maximum absolute difference",
        "     - Maximum diagnostic ratio",
        "     - Status",
    ]
    for spec in PROPERTY_SPECS:
        selected = [result for result in results if result.spec == spec]
        if not selected:
            continue
        passed = sum(result.verification_passed for result in selected)
        maximum_error = max(result.error for result in selected)
        maximum_ratio = max(result.ratio for result in selected)
        lines.extend(
            [
                f"   * - {spec.label}",
                f"     - {len(selected)}",
                f"     - {passed}",
                f"     - {_number(maximum_error)} {spec.unit}",
                f"     - {maximum_ratio:.3f}",
                f"     - {'Pass' if passed == len(selected) else 'Known deviation'}",
            ]
        )
    return lines


def _external_summary(results: Sequence[ExternalResult]) -> list[str]:
    lines = [
        ".. list-table:: fluids 1.3.1 comparison summary",
        "   :header-rows: 1",
        "   :widths: 24 23 24 16 13",
        "",
        "   * - Quantity",
        "     - Criterion",
        "     - Maximum difference",
        "     - Altitude (m)",
        "     - Status",
    ]
    for result in results:
        bounded = result.maximum < result.limit / 2.0
        maximum = (
            format_bounded_error(result.maximum, result.limit, digits=5)
            if bounded
            else _number(result.maximum)
        )
        lines.extend(
            [
                f"   * - {result.label}",
                f"     - {result.criterion}",
                f"     - {maximum}",
                f"     - {'—' if bounded else f'{result.altitude_m:.0f}'}",
                f"     - {'Pass' if result.passed else 'Fail'}",
            ]
        )
    return lines


def _physical_summary(results: Sequence[PhysicalResult]) -> list[str]:
    lines = [
        ".. list-table:: Physical and mathematical invariant summary",
        "   :header-rows: 1",
        "   :widths: 34 34 20 12",
        "",
        "   * - Check",
        "     - Criterion",
        "     - Diagnostic value",
        "     - Status",
    ]
    for result in results:
        diagnostic = (
            _number(result.maximum)
            if result.limit is None
            else format_bounded_error(result.maximum, result.limit, digits=5)
        )
        lines.extend(
            [
                f"   * - {result.label}",
                f"     - {result.criterion}",
                f"     - {diagnostic}",
                f"     - {'Pass' if result.passed else 'Fail'}",
            ]
        )
    return lines


def _known_deviations(results: Sequence[OfficialResult]) -> list[str]:
    failed = [result for result in results if not result.passed]
    if not failed:
        return ["No strict official-table deviations were found."]
    lines = [
        ".. list-table:: Known deviations from the strict printed-digit criterion",
        "   :header-rows: 1",
        "   :widths: 14 14 22 18 18 14",
        "",
        "   * - Coordinate",
        "     - Altitude (m)",
        "     - Quantity",
        "     - Published",
        "     - Computed",
        "     - Tolerance ratio",
    ]
    for result in failed:
        lines.extend(
            [
                f"   * - {result.coordinate_type}",
                f"     - {result.altitude_m:g}",
                f"     - {result.spec.label}",
                f"     - {_number(result.expected)}",
                f"     - {_number(result.actual)}",
                f"     - {result.ratio:.3f}",
            ]
        )
    return lines


def generate_tables() -> str:
    """Return the generated RST verification fragment."""
    official = official_results()
    external = external_results()
    physical = physical_results()
    verification_pass = (
        all(result.verification_passed for result in official)
        and all(result.passed for result in external)
        and all(result.passed for result in physical)
    )
    has_observations = any(not result.passed for result in official)
    if not verification_pass:
        assessment = "Needs revision"
        assessment_text = (
            "At least one verification acceptance criterion failed.  "
            "No production model code was changed as part of this "
            "verification work."
        )
    elif has_observations:
        assessment = "Verified with observations"
        assessment_text = (
            "All verification acceptance criteria pass.  Differences that "
            "exceed the stricter half-of-the-last-printed-digit diagnostic "
            "are retained as observations below; no production model code "
            "was changed as part of this verification work."
        )
    else:
        assessment = "Verified"
        assessment_text = "All verification and diagnostic criteria pass."
    lines = [
        (
            ".. This file is generated by "
            "docs/scripts/generate_standard_atmosphere_validation.py."
        ),
        ".. Do not edit it by hand.",
        "",
        "Overall assessment",
        "^^^^^^^^^^^^^^^^^^",
        "",
        f"**{assessment}**",
        "",
        assessment_text,
        "",
    ]
    lines.extend(_official_summary(official))
    lines.extend(["", ""])
    lines.extend(_external_summary(external))
    lines.extend(["", ""])
    lines.extend(_physical_summary(physical))
    lines.extend(["", ""])
    lines.extend(_known_deviations(official))
    lines.append("")
    return "\n".join(lines)


def _scale(
    values: NDArray[np.float64],
    domain: tuple[float, float],
    output: tuple[float, float],
    *,
    kind: Literal["linear", "log"] = "linear",
) -> NDArray[np.float64]:
    source = values
    lower, upper = domain
    if kind == "log":
        source = np.log10(source)
        lower, upper = np.log10([lower, upper])
    return output[0] + (source - lower) / (upper - lower) * (output[1] - output[0])


def _panel(
    *,
    x: float,
    y: float,
    width: float,
    height: float,
    title: str,
    x_values: NDArray[np.float64],
    series: Sequence[Series],
    y_domain: tuple[float, float],
    y_ticks: Sequence[float],
    y_label: str,
    y_kind: Literal["linear", "log"] = "linear",
    reference_line: float | None = None,
) -> list[str]:
    left, right = x + 74.0, x + width - 22.0
    top, bottom = y + 42.0, y + height - 54.0
    lines = [
        f'<g aria-label="{escape(title)}">',
        (f'<rect class="panel" x="{x}" y="{y}" width="{width}" height="{height}"/>'),
        (
            f'<text class="panel-title" x="{x + width / 2}" '
            f'y="{y + 25}" text-anchor="middle">{escape(title)}</text>'
        ),
    ]
    for value in y_ticks:
        position = float(
            _scale(np.asarray([value]), y_domain, (bottom, top), kind=y_kind)[0]
        )
        lines.extend(
            [
                (
                    f'<line class="grid" x1="{left}" y1="{position:.2f}" '
                    f'x2="{right}" y2="{position:.2f}"/>'
                ),
                (
                    f'<text class="tick" x="{left - 8}" '
                    f'y="{position + 4:.2f}" text-anchor="end">'
                    f"{_number(value)}</text>"
                ),
            ]
        )
    for value in (-5.0, 0.0, 20.0, 40.0, 60.0, 80.0):
        position = float(_scale(np.asarray([value]), (-5.0, 86.0), (left, right))[0])
        lines.extend(
            [
                (
                    f'<line class="grid" x1="{position:.2f}" y1="{top}" '
                    f'x2="{position:.2f}" y2="{bottom}"/>'
                ),
                (
                    f'<text class="tick" x="{position:.2f}" '
                    f'y="{bottom + 19}" text-anchor="middle">{value:g}</text>'
                ),
            ]
        )
    lines.extend(
        [
            (
                f'<line class="axis" x1="{left}" y1="{bottom}" '
                f'x2="{right}" y2="{bottom}"/>'
            ),
            (f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{bottom}"/>'),
            (
                f'<text class="axis-label" x="{(left + right) / 2}" '
                f'y="{y + height - 13}" text-anchor="middle">'
                "Geometric altitude (km)</text>"
            ),
            (
                f'<text class="axis-label" x="{x + 17}" '
                f'y="{(top + bottom) / 2}" text-anchor="middle" '
                f'transform="rotate(-90 {x + 17} {(top + bottom) / 2})">'
                f"{escape(y_label)}</text>"
            ),
        ]
    )
    if reference_line is not None:
        position = float(
            _scale(
                np.asarray([reference_line]),
                y_domain,
                (bottom, top),
                kind=y_kind,
            )[0]
        )
        lines.append(
            f'<line class="reference" x1="{left}" y1="{position:.2f}" '
            f'x2="{right}" y2="{position:.2f}"/>'
        )
    scaled_x = _scale(x_values, (-5.0, 86.0), (left, right))
    for item in series:
        bounded = np.clip(item.values, y_domain[0], y_domain[1])
        scaled_y = _scale(bounded, y_domain, (bottom, top), kind=y_kind)
        points = " ".join(
            f"{x_value:.2f},{y_value:.2f}"
            for x_value, y_value in zip(scaled_x, scaled_y, strict=True)
        )
        lines.append(
            f'<polyline class="curve" stroke="{item.color}" points="{points}"/>'
        )
    for index, item in enumerate(series):
        row_y = top + 15.0 + index * 17.0
        lines.extend(
            [
                (
                    f'<line class="legend-line" stroke="{item.color}" '
                    f'x1="{left + 7}" y1="{row_y}" '
                    f'x2="{left + 31}" y2="{row_y}"/>'
                ),
                (
                    f'<text class="legend" x="{left + 37}" y="{row_y + 4}">'
                    f"{escape(item.label)}</text>"
                ),
            ]
        )
    lines.append("</g>")
    return lines


def _svg_header(*, width: int, height: int, title: str, description: str) -> list[str]:
    return [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
            f'height="{height}" viewBox="0 0 {width} {height}" role="img" '
            'aria-labelledby="svg-title svg-desc">'
        ),
        f'<title id="svg-title">{escape(title)}</title>',
        f'<desc id="svg-desc">{escape(description)}</desc>',
        "<style>",
        (
            "text{font-family:system-ui,-apple-system,BlinkMacSystemFont,"
            "'Segoe UI',sans-serif;fill:#17202a}"
        ),
        ".background{fill:#fff}.panel{fill:#fbfcfd;stroke:#b9c3cc;stroke-width:1}",
        ".grid{stroke:#dfe5ea;stroke-width:1}.axis{stroke:#4d5b66;stroke-width:1.2}",
        ".curve{fill:none;stroke-width:2.3;stroke-linecap:round;stroke-linejoin:round}",
        ".reference{stroke:#777;stroke-width:1.4;stroke-dasharray:5 4}",
        ".legend-line{stroke-width:2.3}.title{font-size:22px;font-weight:650}",
        ".subtitle{font-size:13px;fill:#52616d}.panel-title{font-size:14px;font-weight:650}",
        ".axis-label{font-size:12px;font-weight:600}.tick,.legend{font-size:10px}",
        "</style>",
        f'<rect class="background" width="{width}" height="{height}"/>',
        (
            f'<text class="title" x="{width / 2}" y="29" '
            f'text-anchor="middle">{escape(title)}</text>'
        ),
    ]


def generate_profile_svg() -> str:
    """Return a six-panel physical-behaviour figure."""
    standard_atmosphere, _ = _load_public_api()
    altitude_m = np.linspace(-5_000.0, 86_000.0, 700, dtype=np.float64)
    altitude_km = altitude_m / 1_000.0
    state = standard_atmosphere(altitude_m)
    panels = (
        (
            "Temperature layer structure",
            "temperature",
            (180.0, 325.0),
            (200.0, 240.0, 280.0, 320.0),
            "Temperature (K)",
            "linear",
        ),
        (
            "Pressure decreases with altitude",
            "pressure",
            (0.3, 2.0e5),
            (1.0, 1.0e2, 1.0e4, 1.0e5),
            "Pressure (Pa, log scale)",
            "log",
        ),
        (
            "Density decreases with altitude",
            "density",
            (5.0e-6, 2.0),
            (1.0e-5, 1.0e-3, 1.0e-1, 1.0),
            "Density (kg/m³, log scale)",
            "log",
        ),
        (
            "Speed of sound follows temperature",
            "speed_of_sound",
            (270.0, 365.0),
            (280.0, 300.0, 320.0, 340.0, 360.0),
            "Speed of sound (m/s)",
            "linear",
        ),
        (
            "Gravity decreases with altitude",
            "gravity",
            (9.53, 9.84),
            (9.55, 9.65, 9.75, 9.85),
            "Gravity (m/s²)",
            "linear",
        ),
        (
            "Kinematic viscosity increases",
            "kinematic_viscosity",
            (8.0e-6, 3.0),
            (1.0e-5, 1.0e-3, 1.0e-1, 1.0),
            "Kinematic viscosity (m²/s, log scale)",
            "log",
        ),
    )
    lines = _svg_header(
        width=1200,
        height=900,
        title="U.S. Standard Atmosphere 1976 physical behaviour",
        description=(
            "Six panels show temperature, pressure, density, speed of sound, "
            "gravity, and kinematic viscosity from minus 5 to 86 kilometres."
        ),
    )
    lines.append(
        '<text class="subtitle" x="600" y="50" text-anchor="middle">'
        "Public API evaluated on a dense geometric-altitude grid; "
        "logarithmic axes are labelled explicitly.</text>"
    )
    profile_colors = {spec.state_field: spec.color for spec in PROPERTY_SPECS}
    for index, (title, field, domain, ticks, label, kind) in enumerate(panels):
        column, row = index % 2, index // 2
        lines.extend(
            _panel(
                x=28.0 + column * 586.0,
                y=66.0 + row * 274.0,
                width=558.0,
                height=252.0,
                title=title,
                x_values=altitude_km,
                series=(
                    Series(
                        title,
                        profile_colors[field],
                        np.asarray(getattr(state, field), dtype=np.float64),
                    ),
                ),
                y_domain=domain,
                y_ticks=ticks,
                y_label=label,
                y_kind=kind,
            )
        )
    lines.append("</svg>")
    return "\n".join(lines) + "\n"


def generate_comparison_svg() -> str:
    """Return official-table and independent-implementation difference plots."""
    standard_atmosphere, _ = _load_public_api()
    official = [
        result for result in official_results() if result.coordinate_type == "geometric"
    ]
    altitude, references = _external_arrays()
    state = standard_atmosphere(altitude)
    external_series: list[Series] = []
    external_colors = (
        "#0072B2",
        "#D55E00",
        "#009E73",
        "#CC79A7",
        "#56B4E9",
        "#E69F00",
        "#000000",
        "#6A3D9A",
    )
    for color, (column, state_field, label, _kind, _threshold) in zip(
        external_colors, EXTERNAL_SPECS, strict=True
    ):
        actual = np.asarray(getattr(state, state_field), dtype=np.float64)
        relative_ppm = np.abs(actual / references[column] - 1.0) * 1.0e6
        external_series.append(Series(label, color, np.maximum(relative_ppm, 1.0e-9)))

    official_series: list[Series] = []
    official_altitude = (
        np.asarray(sorted({result.altitude_m for result in official})) / 1_000.0
    )
    for spec in PROPERTY_SPECS:
        selected = {
            result.altitude_m: result for result in official if result.spec == spec
        }
        if not selected:
            continue
        ratios = np.asarray(
            [
                selected[value * 1_000.0].ratio
                if value * 1_000.0 in selected
                else np.nan
                for value in official_altitude
            ],
            dtype=np.float64,
        )
        valid = np.isfinite(ratios)
        official_series.append(
            Series(
                spec.label,
                spec.color,
                np.where(valid, np.maximum(ratios, 1.0e-3), 1.0e-3),
            )
        )

    lines = _svg_header(
        width=1200,
        height=470,
        title="Standard-atmosphere verification differences",
        description=(
            "The left panel shows absolute relative differences from "
            "fluids 1.3.1 in parts per million. The right panel shows "
            "absolute differences divided by official printed-cell tolerances."
        ),
    )
    lines.append(
        '<text class="subtitle" x="600" y="50" text-anchor="middle">'
        "The official line is a strict diagnostic; the software comparison "
        "uses dimensionless relative differences.</text>"
    )
    lines.extend(
        _panel(
            x=28.0,
            y=66.0,
            width=558.0,
            height=374.0,
            title="Absolute relative difference from fluids 1.3.1",
            x_values=altitude / 1_000.0,
            series=external_series,
            y_domain=(1.0e-9, 30.0),
            y_ticks=(1.0e-6, 1.0e-3, 1.0, 10.0),
            y_label="Absolute relative difference (ppm, log scale)",
            y_kind="log",
        )
    )
    lines.extend(
        _panel(
            x=614.0,
            y=66.0,
            width=558.0,
            height=374.0,
            title="Difference relative to official printed tolerance",
            x_values=official_altitude,
            series=official_series,
            y_domain=(1.0e-3, 5.0),
            y_ticks=(1.0e-2, 1.0e-1, 1.0, 5.0),
            y_label="Absolute difference / tolerance (log scale)",
            y_kind="log",
            reference_line=1.0,
        )
    )
    lines.append("</svg>")
    return "\n".join(lines) + "\n"


def _update(path: Path, expected: str, *, check: bool) -> bool:
    if check:
        if not path.is_file() or path.read_text(encoding="utf-8") != expected:
            print(f"out of date: {path.relative_to(PROJECT_ROOT)}", file=sys.stderr)
            return False
        return True
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(expected, encoding="utf-8")
    print(f"wrote {path.relative_to(PROJECT_ROOT)}")
    return True


def main() -> int:
    """Generate assets or fail when checked-in assets are stale."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if a generated RST or SVG asset is out of date",
    )
    arguments = parser.parse_args()
    results = (
        _update(TABLE_PATH, generate_tables(), check=arguments.check),
        _update(PROFILE_SVG_PATH, generate_profile_svg(), check=arguments.check),
        _update(
            COMPARISON_SVG_PATH,
            generate_comparison_svg(),
            check=arguments.check,
        ),
    )
    return 0 if all(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
