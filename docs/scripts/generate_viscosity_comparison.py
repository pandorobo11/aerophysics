"""Generate the viscosity-model comparison figure and tables."""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Sequence

    from numpy.typing import NDArray

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SVG_PATH = PROJECT_ROOT / "docs/_static/viscosity_model_comparison.svg"
TABLE_PATH = PROJECT_ROOT / "docs/_generated/viscosity_model_comparison.rst"

LOW_TABLE_TEMPERATURES = np.array([79.0, 100.0, 300.0, 1000.0, 1845.0])
HIGH_TABLE_TEMPERATURES = np.array([1000.0, 1845.0, 5000.0, 10000.0, 30000.0])


@dataclass(frozen=True)
class Series:
    """One plotted viscosity or relative-difference series."""

    label: str
    color: str
    values: NDArray[np.float64]


def _models() -> tuple[object, object, object]:
    """Load the models through their documented public API."""
    sys.path.insert(0, str(PROJECT_ROOT / "src"))
    from aerophysics.transport import (
        AIR_BLOTTNER_VISCOSITY,
        AIR_KEYES_VISCOSITY,
        AIR_VISCOSITY,
    )

    return AIR_VISCOSITY, AIR_KEYES_VISCOSITY, AIR_BLOTTNER_VISCOSITY


def _evaluate(model: Any, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.asarray(model.dynamic_viscosity(temperature), dtype=np.float64)


def _relative_difference(
    candidate: NDArray[np.float64], reference: NDArray[np.float64]
) -> NDArray[np.float64]:
    return (candidate / reference - 1.0) * 100.0


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
    fraction = (source - lower) / (upper - lower)
    return output[0] + fraction * (output[1] - output[0])


def _number(value: float) -> str:
    if value == 0.0:
        return "0"
    if abs(value) < 0.001:
        return f"{value:.0e}"
    return f"{value:g}"


def _panel(
    *,
    x: float,
    y: float,
    width: float,
    height: float,
    title: str,
    temperatures: NDArray[np.float64],
    series: Sequence[Series],
    x_domain: tuple[float, float],
    x_ticks: Sequence[float],
    y_domain: tuple[float, float],
    y_ticks: Sequence[float],
    y_label: str,
    x_kind: Literal["linear", "log"] = "linear",
) -> list[str]:
    left = x + 78.0
    right = x + width - 24.0
    top = y + 44.0
    bottom = y + height - 62.0
    lines = [
        f'<g aria-label="{title}">',
        f'<rect class="panel" x="{x}" y="{y}" width="{width}" height="{height}"/>',
        (
            f'<text class="panel-title" x="{x + width / 2}" y="{y + 26}" '
            f'text-anchor="middle">{title}</text>'
        ),
    ]

    y_tick_values = np.asarray(y_ticks, dtype=np.float64)
    y_positions = _scale(y_tick_values, y_domain, (bottom, top))
    for value, position in zip(y_tick_values, y_positions, strict=True):
        lines.extend(
            [
                (
                    f'<line class="grid" x1="{left}" y1="{position:.2f}" '
                    f'x2="{right}" y2="{position:.2f}"/>'
                ),
                (
                    f'<text class="tick" x="{left - 9}" '
                    f'y="{position + 4:.2f}" text-anchor="end">'
                    f"{_number(float(value))}</text>"
                ),
            ]
        )

    x_tick_values = np.asarray(x_ticks, dtype=np.float64)
    x_positions = _scale(x_tick_values, x_domain, (left, right), kind=x_kind)
    for value, position in zip(x_tick_values, x_positions, strict=True):
        lines.extend(
            [
                (
                    f'<line class="grid" x1="{position:.2f}" y1="{top}" '
                    f'x2="{position:.2f}" y2="{bottom}"/>'
                ),
                (
                    f'<text class="tick" x="{position:.2f}" y="{bottom + 20}" '
                    f'text-anchor="middle">{value:g}</text>'
                ),
            ]
        )

    lines.extend(
        [
            (
                f'<line class="axis" x1="{left}" y1="{bottom}" '
                f'x2="{right}" y2="{bottom}"/>'
            ),
            f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{bottom}"/>',
            (
                f'<text class="axis-label" x="{(left + right) / 2}" '
                f'y="{y + height - 17}" text-anchor="middle">'
                "Temperature (K)</text>"
            ),
            (
                f'<text class="axis-label" x="{x + 19}" '
                f'y="{(top + bottom) / 2}" text-anchor="middle" '
                f'transform="rotate(-90 {x + 19} {(top + bottom) / 2})">'
                f"{y_label}</text>"
            ),
        ]
    )

    x_values = _scale(temperatures, x_domain, (left, right), kind=x_kind)
    for item in series:
        y_values = _scale(item.values, y_domain, (bottom, top))
        points = " ".join(
            f"{x_value:.2f},{y_value:.2f}"
            for x_value, y_value in zip(x_values, y_values, strict=True)
        )
        lines.append(
            f'<polyline class="curve" stroke="{item.color}" points="{points}"/>'
        )

    legend_x = left + 8.0
    legend_y = top + 15.0
    for index, item in enumerate(series):
        row_y = legend_y + index * 19.0
        lines.extend(
            [
                (
                    f'<line class="legend-line" stroke="{item.color}" '
                    f'x1="{legend_x}" y1="{row_y}" x2="{legend_x + 28}" '
                    f'y2="{row_y}"/>'
                ),
                (
                    f'<text class="legend" x="{legend_x + 35}" '
                    f'y="{row_y + 4}">{item.label}</text>'
                ),
            ]
        )
    lines.append("</g>")
    return lines


def generate_svg() -> str:
    """Return the deterministic four-panel comparison SVG."""
    sutherland, keyes, blottner = _models()
    low_temperature = np.linspace(79.0, 1845.0, 300, dtype=np.float64)
    high_temperature = np.geomspace(1000.0, 30000.0, 300).astype(np.float64)

    low_sutherland = _evaluate(sutherland, low_temperature)
    low_keyes = _evaluate(keyes, low_temperature)
    high_sutherland = _evaluate(sutherland, high_temperature)
    high_blottner = _evaluate(blottner, high_temperature)

    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            '<svg xmlns="http://www.w3.org/2000/svg" width="1280" height="930" '
            'viewBox="0 0 1280 930" role="img" '
            'aria-labelledby="svg-title svg-desc">'
        ),
        '<title id="svg-title">Comparison of dry-air dynamic-viscosity models</title>',
        (
            '<desc id="svg-desc">Four panels compare Sutherland with Keyes from '
            "79 to 1845 kelvin and Sutherland with the frozen-composition "
            "Blottner and Wilke model from 1000 to 30000 kelvin, showing "
            "viscosity and relative difference.</desc>"
        ),
        "<style>",
        (
            "text{font-family:system-ui,-apple-system,BlinkMacSystemFont,"
            "'Segoe UI',sans-serif;fill:#17202a}"
        ),
        ".background{fill:#fff}.panel{fill:#fbfcfd;stroke:#b9c3cc;stroke-width:1}",
        ".grid{stroke:#dfe5ea;stroke-width:1}.axis{stroke:#4d5b66;stroke-width:1.2}",
        ".curve{fill:none;stroke-width:2.6;stroke-linecap:round;stroke-linejoin:round}",
        ".legend-line{stroke-width:2.6}.title{font-size:22px;font-weight:650}",
        ".subtitle{font-size:13px;fill:#52616d}.panel-title{font-size:15px;font-weight:650}",
        ".axis-label{font-size:12px;font-weight:600}.tick,.legend{font-size:11px}",
        "</style>",
        '<rect class="background" width="1280" height="930"/>',
        (
            '<text class="title" x="640" y="31" text-anchor="middle">'
            "Dry-air dynamic-viscosity model comparison</text>"
        ),
        (
            '<text class="subtitle" x="640" y="53" text-anchor="middle">'
            "Only each fitted model's nominal range is shown; Sutherland is the "
            "comparison baseline.</text>"
        ),
    ]
    lines.extend(
        _panel(
            x=35.0,
            y=68.0,
            width=590.0,
            height=390.0,
            title="Low temperature: dynamic viscosity",
            temperatures=low_temperature,
            series=(
                Series("Sutherland", "#0072B2", low_sutherland),
                Series("Keyes", "#D55E00", low_keyes),
            ),
            x_domain=(79.0, 1845.0),
            x_ticks=(79.0, 300.0, 600.0, 1000.0, 1400.0, 1845.0),
            y_domain=(0.0, 6.5e-5),
            y_ticks=(0.0, 2e-5, 4e-5, 6e-5),
            y_label="Dynamic viscosity (Pa·s)",
        )
    )
    lines.extend(
        _panel(
            x=655.0,
            y=68.0,
            width=590.0,
            height=390.0,
            title="Low temperature: Keyes relative to Sutherland",
            temperatures=low_temperature,
            series=(
                Series(
                    "(Keyes / Sutherland - 1) x 100",
                    "#D55E00",
                    _relative_difference(low_keyes, low_sutherland),
                ),
            ),
            x_domain=(79.0, 1845.0),
            x_ticks=(79.0, 300.0, 600.0, 1000.0, 1400.0, 1845.0),
            y_domain=(0.0, 5.2),
            y_ticks=(0.0, 1.0, 2.0, 3.0, 4.0, 5.0),
            y_label="Relative difference (%)",
        )
    )
    lines.extend(
        _panel(
            x=35.0,
            y=488.0,
            width=590.0,
            height=390.0,
            title="High temperature: dynamic viscosity",
            temperatures=high_temperature,
            series=(
                Series("Sutherland", "#0072B2", high_sutherland),
                Series("Blottner + Wilke", "#009E73", high_blottner),
            ),
            x_domain=(1000.0, 30000.0),
            x_ticks=(1000.0, 2000.0, 5000.0, 10000.0, 20000.0, 30000.0),
            y_domain=(0.0, 6.0e-4),
            y_ticks=(0.0, 2e-4, 4e-4, 6e-4),
            y_label="Dynamic viscosity (Pa·s)",
            x_kind="log",
        )
    )
    lines.extend(
        _panel(
            x=655.0,
            y=488.0,
            width=590.0,
            height=390.0,
            title="High temperature: Blottner/Wilke relative to Sutherland",
            temperatures=high_temperature,
            series=(
                Series(
                    "(Blottner/Wilke / Sutherland - 1) x 100",
                    "#009E73",
                    _relative_difference(high_blottner, high_sutherland),
                ),
            ),
            x_domain=(1000.0, 30000.0),
            x_ticks=(1000.0, 2000.0, 5000.0, 10000.0, 20000.0, 30000.0),
            y_domain=(-10.0, 125.0),
            y_ticks=(0.0, 25.0, 50.0, 75.0, 100.0, 125.0),
            y_label="Relative difference (%)",
            x_kind="log",
        )
    )
    lines.extend(
        [
            (
                '<text class="subtitle" x="640" y="910" text-anchor="middle">'
                "Blottner/Wilke: frozen N₂/O₂/Ar/CO₂ dry-air composition; no "
                "dissociation or chemical reactions.</text>"
            ),
            "</svg>",
        ]
    )
    return "\n".join(lines) + "\n"


def _table(
    *,
    title: str,
    temperatures: NDArray[np.float64],
    reference: NDArray[np.float64],
    candidate: NDArray[np.float64],
    candidate_name: str,
) -> list[str]:
    lines = [
        f".. list-table:: {title}",
        "   :header-rows: 1",
        "   :widths: 14 24 24 25",
        "",
        "   * - Temperature (K)",
        "     - Sutherland (Pa s)",
        f"     - {candidate_name} (Pa s)",
        "     - Relative difference (%)",
    ]
    differences = _relative_difference(candidate, reference)
    for temperature, baseline, value, difference in zip(
        temperatures, reference, candidate, differences, strict=True
    ):
        lines.extend(
            [
                f"   * - {temperature:g}",
                f"     - {baseline:.7e}",
                f"     - {value:.7e}",
                f"     - {difference:+.3f}",
            ]
        )
    return lines


def generate_tables() -> str:
    """Return the deterministic RST table fragment."""
    sutherland, keyes, blottner = _models()
    low_sutherland = _evaluate(sutherland, LOW_TABLE_TEMPERATURES)
    low_keyes = _evaluate(keyes, LOW_TABLE_TEMPERATURES)
    high_sutherland = _evaluate(sutherland, HIGH_TABLE_TEMPERATURES)
    high_blottner = _evaluate(blottner, HIGH_TABLE_TEMPERATURES)
    lines = [
        ".. This file is generated by docs/scripts/generate_viscosity_comparison.py.",
        ".. Do not edit it by hand.",
        "",
    ]
    lines.extend(
        _table(
            title="Low-temperature model values",
            temperatures=LOW_TABLE_TEMPERATURES,
            reference=low_sutherland,
            candidate=low_keyes,
            candidate_name="Keyes",
        )
    )
    lines.extend(["", ""])
    lines.extend(
        _table(
            title="High-temperature model values",
            temperatures=HIGH_TABLE_TEMPERATURES,
            reference=high_sutherland,
            candidate=high_blottner,
            candidate_name="Blottner/Wilke",
        )
    )
    lines.append("")
    return "\n".join(lines)


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
    """Generate assets, or check that committed assets are current."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if the generated SVG or RST fragment is out of date",
    )
    arguments = parser.parse_args()
    results = (
        _update(SVG_PATH, generate_svg(), check=arguments.check),
        _update(TABLE_PATH, generate_tables(), check=arguments.check),
    )
    return 0 if all(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
