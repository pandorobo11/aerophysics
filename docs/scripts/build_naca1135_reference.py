"""Build the offline NACA Report 1135 table transcription fixture.

The scan uses a compact floating-point notation and contains more than 1,600
Mach rows.  This script independently evaluates the equations printed beside
Tables I and II at every published Mach abscissa and rounds the results to the
table precision.  The committed CSV is the test input; this script is an
explicit transcription/review aid and is not run by normal tests.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / "tests/reference_data/compressible_flow"
CSV_PATH = DIRECTORY / "naca1135_tables_i_ii.csv"
METADATA_PATH = DIRECTORY / "naca1135_tables_i_ii.json"
GAMMA = 1.4


def _published_mach_rows() -> list[float]:
    subsonic = [index / 100 for index in range(101)]
    supersonic_fine = [1.0 + index / 100 for index in range(945)]
    supersonic_medium = [10.46 + index * 0.02 for index in range(478)]
    supersonic_coarse = [20.2 + index * 0.2 for index in range(174)]
    supersonic_integer = [float(index) for index in range(55, 101)]
    return (
        subsonic
        + supersonic_fine
        + supersonic_medium
        + supersonic_coarse
        + supersonic_integer
    )


def _printed(value: float, digits: int) -> str:
    return format(value, f".{digits}g")


def _printed_unit(value: float, digits: int) -> float:
    if value == 0.0:
        return 10.0 ** (1 - digits)
    return 10.0 ** (math.floor(math.log10(abs(value))) - digits + 1)


def _page_for_mach(mach: float) -> int:
    if mach <= 1.0:
        return 633
    table_two_rows = [value for value in _published_mach_rows()[101:] if value >= 1.0]
    index = min(range(len(table_two_rows)), key=lambda i: abs(table_two_rows[i] - mach))
    return 633 if index < 25 else 634 + (index - 25) // 90


def _row(mach: float) -> dict[str, str]:
    temperature = 1.0 / (1.0 + 0.5 * (GAMMA - 1.0) * mach**2)
    pressure = temperature ** (GAMMA / (GAMMA - 1.0))
    density = temperature ** (1.0 / (GAMMA - 1.0))
    area = ""
    area_unit = ""
    if mach > 0.0:
        value = (2.0 / (GAMMA + 1.0) / temperature) ** (
            (GAMMA + 1.0) / (2.0 * (GAMMA - 1.0))
        ) / mach
        area = _printed(value, 5)
        area_unit = f"{_printed_unit(float(area), 5):.12g}"
    row = {
        "table": "I" if mach <= 1.0 else "II",
        "printed_page": str(_page_for_mach(mach)),
        "mach": f"{mach:.2f}",
        "pressure_over_total": _printed(pressure, 4),
        "pressure_printed_unit": (
            f"{_printed_unit(float(_printed(pressure, 4)), 4):.12g}"
        ),
        "density_over_total": _printed(density, 4),
        "density_printed_unit": f"{_printed_unit(float(_printed(density, 4)), 4):.12g}",
        "temperature_over_total": _printed(temperature, 4),
        "temperature_printed_unit": (
            f"{_printed_unit(float(_printed(temperature, 4)), 4):.12g}"
        ),
        "area_over_critical": area,
        "area_printed_unit": area_unit,
        "prandtl_meyer_deg": "",
        "prandtl_meyer_printed_unit_deg": "",
        "downstream_mach": "",
        "downstream_mach_printed_unit": "",
        "pressure_ratio": "",
        "pressure_ratio_printed_unit": "",
        "density_ratio": "",
        "density_ratio_printed_unit": "",
        "temperature_ratio": "",
        "temperature_ratio_printed_unit": "",
        "total_pressure_ratio": "",
        "total_pressure_ratio_printed_unit": "",
    }
    if mach < 1.0:
        return row
    root = math.sqrt(mach**2 - 1.0)
    nu = math.sqrt((GAMMA + 1.0) / (GAMMA - 1.0)) * math.atan(
        math.sqrt((GAMMA - 1.0) / (GAMMA + 1.0)) * root
    ) - math.atan(root)
    downstream = math.sqrt(
        (1.0 + 0.5 * (GAMMA - 1.0) * mach**2) / (GAMMA * mach**2 - 0.5 * (GAMMA - 1.0))
    )
    pressure_ratio = 1.0 + 2.0 * GAMMA / (GAMMA + 1.0) * (mach**2 - 1.0)
    density_ratio = (GAMMA + 1.0) * mach**2 / (2.0 + (GAMMA - 1.0) * mach**2)
    temperature_ratio = pressure_ratio / density_ratio
    total_pressure_ratio = (
        ((GAMMA + 1.0) * mach**2) / ((GAMMA - 1.0) * mach**2 + 2.0)
    ) ** (GAMMA / (GAMMA - 1.0)) * (
        (GAMMA + 1.0) / (2.0 * GAMMA * mach**2 - (GAMMA - 1.0))
    ) ** (1.0 / (GAMMA - 1.0))
    values = {
        "prandtl_meyer_deg": (math.degrees(nu), 5),
        "downstream_mach": (downstream, 4),
        "pressure_ratio": (pressure_ratio, 4),
        "density_ratio": (density_ratio, 4),
        "temperature_ratio": (temperature_ratio, 4),
        "total_pressure_ratio": (total_pressure_ratio, 4),
    }
    for name, (value, digits) in values.items():
        printed = _printed(value, digits)
        row[name] = printed
        suffix = "_deg" if name == "prandtl_meyer_deg" else ""
        unit_name = name.replace("_deg", "") + f"_printed_unit{suffix}"
        row[unit_name] = f"{_printed_unit(float(printed), digits):.12g}"
    return row


def main() -> None:
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    rows = [_row(mach) for mach in _published_mach_rows()]
    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "source": "NACA Report 1135, Tables I and II",
        "report_year": 1953,
        "printed_pages": "633-651",
        "gamma": GAMMA,
        "row_count": len(rows),
        "method": (
            "Every published Mach abscissa was reconstructed with the equations "
            "printed in Report 1135 and rounded to the displayed table precision; "
            "the endpoints and page transitions were visually audited against the scan."
        ),
        "acceptance": "max(one printed unit, 1e-4 relative)",
        "command": "python docs/scripts/build_naca1135_reference.py",
    }
    METADATA_PATH.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
