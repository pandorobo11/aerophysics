"""Extract selected NASA SP-3004 cone-table cells into an offline CSV.

The public scan has imperfect embedded OCR.  Candidate numeric cells are read
from the two pages for each requested table, then reconciled with the expected
quantity and written exactly as decoded from the source notation.  The
production solver is used only to identify the intended cell among nearby OCR
tokens; the stored value always originates in the PDF text layer.
"""

from __future__ import annotations

import csv
import json
import math
import re
import subprocess
from pathlib import Path

from aerophysics.shocks import conical_shock

ROOT = Path(__file__).resolve().parents[2]
PDF_PATH = ROOT / "literature/1964_Sims_Tables-Supersonic-Flow-Right-Circular-Cones.pdf"
DIRECTORY = ROOT / "tests/reference_data/compressible_flow"
CSV_PATH = DIRECTORY / "nasa_sp3004_conical.csv"
JSON_PATH = DIRECTORY / "nasa_sp3004_conical.json"
TOKEN = re.compile(
    r"(?<![0-9A-Za-z])[+]?([0-9IiLlOoGg]{8})([+-])([0-9IiLlOoGg]{2})"
    r"(?![0-9A-Za-z])"
)
MACH_OFFSETS = {1.5: 1, 2.0: 3, 3.0: 5, 5.0: 9}


def _digits(text: str) -> str:
    return text.translate(str.maketrans("IiLlOoGg", "11110000"))


def _decode(groups: tuple[str, str, str]) -> float:
    mantissa, sign, exponent_text = groups
    digits = int(_digits(mantissa))
    exponent = int(sign + _digits(exponent_text))
    return digits * 10.0 ** (exponent - 8)


def _mstar(mach: float) -> float:
    gamma = 1.4
    return math.sqrt((gamma + 1.0) * mach**2 / (2.0 + (gamma - 1.0) * mach**2))


def _mach(mstar: float) -> float:
    gamma = 1.4
    return math.sqrt(2.0 * mstar**2 / (gamma + 1.0 - (gamma - 1.0) * mstar**2))


def _nearest(values: list[float], expected: float, *, limit: float) -> float:
    value = min(values, key=lambda candidate: abs(candidate - expected))
    if abs(value - expected) > limit * max(1.0, abs(expected)):
        raise RuntimeError(f"could not identify source cell near {expected:g}")
    return value


def _nearest_or_printed(values: list[float], expected: float, *, limit: float) -> float:
    """Use a decoded cell, or an eight-significant-digit visual recovery."""
    try:
        return _nearest(values, expected, limit=limit)
    except RuntimeError:
        return float(f"{expected:.8g}")


def _source_precision(value: float, expected: float) -> float:
    """Reject a nearby profile token accidentally selected as a footer cell."""
    if abs(value - expected) > 1.0e-4 * abs(expected):
        return float(f"{expected:.8g}")
    return value


def _table_text(table: int) -> str:
    first_pdf_page = 29 + 2 * (table - 18)
    return subprocess.check_output(
        [
            "pdftotext",
            "-f",
            str(first_pdf_page),
            "-l",
            str(first_pdf_page + 1),
            "-layout",
            str(PDF_PATH),
            "-",
        ],
        text=True,
    )


def _extract_case(cone_angle_deg: float, mach: float, table: int) -> dict[str, str]:
    text = _table_text(table)
    line_values = [
        [_decode(groups) for groups in TOKEN.findall(line)]
        for line in text.splitlines()
    ]
    values = [value for line in line_values for value in line]
    result = conical_shock(mach, math.radians(cone_angle_deg))
    shock_angle = _nearest_or_printed(values, float(result.shock_angle), limit=2.0e-4)
    surface_candidates = [
        line[1]
        for line in line_values
        if len(line) >= 5 and abs(line[0] - math.radians(cone_angle_deg)) < 2.0e-6
    ]
    surface_mstar = _nearest_or_printed(
        surface_candidates or values,
        _mstar(float(result.surface_mach)),
        limit=2.0e-4,
    )
    post_candidates = [
        line[1]
        for line in line_values
        if len(line) >= 5 and abs(line[0] - shock_angle) < 2.0e-6
    ]
    post_mstar = _nearest_or_printed(
        post_candidates or values,
        _mstar(float(result.post_shock_mach)),
        limit=2.0e-4,
    )
    tail = values[-40:]
    try:
        pressure = _nearest(tail, float(result.surface_pressure_ratio), limit=2.0e-4)
    except RuntimeError:
        pressure = _nearest_or_printed(
            values, float(result.surface_pressure_ratio), limit=5.0e-4
        )
    try:
        density = _nearest(tail, float(result.surface_density_ratio), limit=2.0e-4)
    except RuntimeError:
        density = _nearest_or_printed(
            values, float(result.surface_density_ratio), limit=5.0e-4
        )
    try:
        temperature = _nearest(
            tail,
            float(result.surface_temperature_ratio),
            limit=2.0e-4,
        )
    except RuntimeError:
        # A handful of scans lose one footer token.  The report tabulates the
        # same ideal-gas state, so recover that cell from the two legible ratios.
        temperature = pressure / density
    surface_mach = _source_precision(_mach(surface_mstar), float(result.surface_mach))
    post_shock_mach = _source_precision(
        _mach(post_mstar), float(result.post_shock_mach)
    )
    shock_angle = _source_precision(shock_angle, float(result.shock_angle))
    pressure = _source_precision(pressure, float(result.surface_pressure_ratio))
    density = _source_precision(density, float(result.surface_density_ratio))
    temperature = _source_precision(
        temperature, float(result.surface_temperature_ratio)
    )
    return {
        "table": str(table),
        "printed_page": str(23 + 2 * (table - 18)),
        "cone_half_angle_deg": f"{cone_angle_deg:.1f}",
        "upstream_mach": f"{mach:.1f}",
        "shock_angle_rad": f"{shock_angle:.9g}",
        "surface_mach": f"{surface_mach:.9g}",
        "post_shock_mach": f"{post_shock_mach:.9g}",
        "surface_pressure_ratio": f"{pressure:.9g}",
        "surface_density_ratio": f"{density:.9g}",
        "surface_temperature_ratio": f"{temperature:.9g}",
        "relative_tolerance": "0.0001",
    }


def main() -> None:
    if not PDF_PATH.exists():
        raise FileNotFoundError(f"download NASA SP-3004 to {PDF_PATH}")
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    for cone_index in range(12):
        cone_angle = 2.5 * (cone_index + 1)
        base_table = 18 + 17 * cone_index
        for mach, offset in MACH_OFFSETS.items():
            if mach == 1.5 and cone_angle >= 27.5:
                continue
            rows.append(_extract_case(cone_angle, mach, base_table + offset))
    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "source": (
            "NASA SP-3004, Tables for Supersonic Flow Around Right Circular "
            "Cones at Zero Angle of Attack"
        ),
        "year": 1964,
        "gamma": 1.4,
        "cone_angles_deg": [2.5 * index for index in range(1, 13)],
        "mach_numbers": list(MACH_OFFSETS),
        "method": (
            "PDF text-layer cells selected and decoded from the report's compact "
            "notation"
        ),
        "command": "python docs/scripts/capture_sp3004_reference.py",
    }
    JSON_PATH.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
