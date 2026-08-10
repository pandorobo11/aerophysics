"""Small dependency-free helpers shared by verification generators."""

from __future__ import annotations

import sys
from html import escape
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Sequence

    from numpy.typing import NDArray


def write_or_check(path: Path, content: str, *, check: bool) -> bool:
    """Write ``content`` or return whether an existing file is current."""
    normalized = content.rstrip() + "\n"
    if check:
        current = path.exists() and path.read_text(encoding="utf-8") == normalized
        if not current:
            try:
                display_path = path.relative_to(Path.cwd())
            except ValueError:
                display_path = path
            print(f"out of date or missing: {display_path}", file=sys.stderr)
        return current
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(normalized, encoding="utf-8")
    return True


def line_chart_svg(
    *,
    title: str,
    description: str,
    x_label: str,
    y_label: str,
    x: NDArray[np.float64],
    series: Sequence[tuple[str, str, NDArray[np.float64]]],
    log_x: bool = False,
    log_y: bool = False,
) -> str:
    """Return an accessible SVG line chart with labelled axes and legend."""
    width, height = 920, 520
    left, right, top, bottom = 92, 28, 62, 78
    plot_width = width - left - right
    plot_height = height - top - bottom
    x_values = np.log10(x) if log_x else x
    all_y = np.concatenate([values for _, _, values in series])
    y_values = np.log10(all_y) if log_y else all_y
    x_min, x_max = float(np.min(x_values)), float(np.max(x_values))
    y_min, y_max = float(np.min(y_values)), float(np.max(y_values))
    if y_min == y_max:
        y_min -= 0.5
        y_max += 0.5
    padding = 0.05 * (y_max - y_min)
    y_min -= padding
    y_max += padding

    def point(x_value: float, y_value: float) -> tuple[float, float]:
        transformed_x = np.log10(x_value) if log_x else x_value
        transformed_y = np.log10(y_value) if log_y else y_value
        px = left + (transformed_x - x_min) / (x_max - x_min) * plot_width
        py = top + (y_max - transformed_y) / (y_max - y_min) * plot_height
        return px, py

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        f'<title id="title">{escape(title)}</title>',
        f'<desc id="desc">{escape(description)}</desc>',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2:g}" y="32" text-anchor="middle" '
        f'font-family="sans-serif" font-size="20">{escape(title)}</text>',
        f'<line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" '
        f'y2="{top + plot_height}" stroke="#333"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" '
        'stroke="#333"/>',
    ]
    for index in range(6):
        fraction = index / 5
        px = left + fraction * plot_width
        py = top + (1.0 - fraction) * plot_height
        x_tick = x_min + fraction * (x_max - x_min)
        y_tick = y_min + fraction * (y_max - y_min)
        x_text = f"{10**x_tick:.3g}" if log_x else f"{x_tick:.3g}"
        y_text = f"{10**y_tick:.3g}" if log_y else f"{y_tick:.3g}"
        lines.extend(
            [
                f'<line x1="{px:g}" y1="{top}" x2="{px:g}" '
                f'y2="{top + plot_height}" stroke="#ddd"/>',
                f'<text x="{px:g}" y="{top + plot_height + 24}" text-anchor="middle" '
                f'font-family="sans-serif" font-size="12">{escape(x_text)}</text>',
                f'<line x1="{left}" y1="{py:g}" x2="{left + plot_width}" '
                f'y2="{py:g}" stroke="#ddd"/>',
                f'<text x="{left - 10}" y="{py + 4:g}" text-anchor="end" '
                f'font-family="sans-serif" font-size="12">{escape(y_text)}</text>',
            ]
        )
    for _label, color, values in series:
        points = " ".join(
            f"{px:.2f},{py:.2f}"
            for px, py in (
                point(x_value, y_value)
                for x_value, y_value in zip(x, values, strict=True)
            )
        )
        lines.append(
            f'<polyline points="{points}" fill="none" stroke="{color}" '
            'stroke-width="2"/>'
        )
    lines.extend(
        [
            f'<text x="{left + plot_width / 2:g}" y="{height - 22}" '
            f'text-anchor="middle" font-family="sans-serif" font-size="14">'
            f"{escape(x_label)}</text>",
            f'<text x="24" y="{top + plot_height / 2:g}" text-anchor="middle" '
            f'transform="rotate(-90 24 {top + plot_height / 2:g})" '
            f'font-family="sans-serif" font-size="14">{escape(y_label)}</text>',
        ]
    )
    legend_x = left + 12
    for index, (label, color, _) in enumerate(series):
        y_position = top + 18 + index * 21
        lines.extend(
            [
                f'<line x1="{legend_x}" y1="{y_position}" x2="{legend_x + 28}" '
                f'y2="{y_position}" stroke="{color}" stroke-width="3"/>',
                f'<text x="{legend_x + 36}" y="{y_position + 4}" '
                f'font-family="sans-serif" font-size="12">{escape(label)}</text>',
            ]
        )
    lines.append("</svg>")
    return "\n".join(lines)
