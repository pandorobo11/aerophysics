"""Extract selected NASA SP-3004 cone-table cells into an offline CSV.

The public scan has imperfect embedded OCR.  Each requested table is located
from an explicit, visually verified PDF page, and cells are selected from the
text layer by their position in the table: the cone-surface row, the shock
row, and the four footer cells.  No aerophysics solver is imported, so every
stored value originates in the report.

A cell whose OCR text cannot be decoded unambiguously must be listed in
``MANUAL_CELLS`` with its value transcribed from the scan image and its
location.  Each table is then checked against relations that hold within the
report itself: the shock-angle header matches the shock row, the footer
satisfies ``p = rho T``, the surface temperature ratio matches total enthalpy
conservation for the tabulated surface velocity, and the post-shock velocity
matches the oblique-shock Rankine--Hugoniot relations at the tabulated shock
angle.  Any unresolved or
inconsistent cell aborts the capture with its table and page context.
"""

from __future__ import annotations

import csv
import html
import json
import math
import re
import subprocess
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import Final

ROOT = Path(__file__).resolve().parents[2]
PDF_PATH = ROOT / "literature/1964_Sims_Tables-Supersonic-Flow-Right-Circular-Cones.pdf"
DIRECTORY = ROOT / "tests/reference_data/compressible_flow"
CSV_PATH = DIRECTORY / "nasa_sp3004_conical.csv"
JSON_PATH = DIRECTORY / "nasa_sp3004_conical.json"
GAMMA: Final = 1.4

# (cone half-angle deg, Mach, table, first PDF page).  Each table occupies the
# listed page and the next one.  Locations were read from the scanned table
# headers.  The scan omits some pages, so tables cannot be located by
# arithmetic.  For cones of 27.5 deg and above there is no Mach 1.5 table,
# which shifts the table numbers.  The printed page is PDF page minus 6.
SOURCE_TABLES: Final = (
    (2.5, 1.5, 19, 31),
    (2.5, 2.0, 21, 35),
    (2.5, 3.0, 23, 39),
    (2.5, 5.0, 27, 47),
    (5.0, 1.5, 36, 65),
    (5.0, 2.0, 38, 69),
    (5.0, 3.0, 40, 73),
    (5.0, 5.0, 44, 81),
    (7.5, 1.5, 53, 99),
    (7.5, 2.0, 55, 103),
    (7.5, 3.0, 57, 107),
    (7.5, 5.0, 61, 115),
    (10.0, 1.5, 70, 133),
    (10.0, 2.0, 72, 137),
    (10.0, 3.0, 74, 141),
    (10.0, 5.0, 78, 149),
    (12.5, 1.5, 87, 167),
    (12.5, 2.0, 89, 171),
    (12.5, 3.0, 91, 175),
    (12.5, 5.0, 95, 183),
    (15.0, 1.5, 104, 200),
    (15.0, 2.0, 106, 204),
    (15.0, 3.0, 108, 208),
    (15.0, 5.0, 112, 216),
    (17.5, 1.5, 121, 233),
    (17.5, 2.0, 123, 237),
    (17.5, 3.0, 125, 241),
    (17.5, 5.0, 129, 249),
    (20.0, 1.5, 138, 266),
    (20.0, 2.0, 140, 270),
    (20.0, 3.0, 142, 274),
    (20.0, 5.0, 146, 282),
    (22.5, 1.5, 155, 299),
    (22.5, 2.0, 157, 303),
    (22.5, 3.0, 159, 307),
    (22.5, 5.0, 163, 315),
    (25.0, 1.5, 172, 332),
    (25.0, 2.0, 174, 336),
    (25.0, 3.0, 176, 340),
    (25.0, 5.0, 180, 348),
    (27.5, 2.0, 190, 368),
    (27.5, 3.0, 192, 372),
    (27.5, 5.0, 196, 380),
    (30.0, 2.0, 206, 400),
    (30.0, 3.0, 208, 404),
    (30.0, 5.0, 212, 412),
)
PRINTED_PAGE_OFFSET: Final = 6

FIELDS: Final = (
    "shock_angle",
    "surface_mstar",
    "post_shock_mstar",
    "pressure",
    "density",
    "temperature",
)

# Cells transcribed from the scan image because their OCR text is damaged or
# missing: (table, field) -> (eight-significant-digit value, location).
MANUAL_CELLS: Final[dict[tuple[int, str], tuple[str, str]]] = {
    (23, "post_shock_mstar"): (
        "1.9636795",
        "PDF p. 40, shock row, M* column; OCR '419636795+01'",
    ),
    (36, "temperature"): (
        "1.0174631",
        "PDF p. 66, footer Ts/Too; OCR reads 10174691",
    ),
    (44, "pressure"): (
        "1.4033150",
        "PDF p. 82, footer Ps/Poo; OCR '+I_033150+01-'",
    ),
    (57, "surface_mstar"): (
        "1.9153173",
        "PDF p. 107, cone-surface row, M* column; no OCR word",
    ),
    (89, "surface_mstar"): (
        "1.5214464",
        "PDF p. 171, cone-surface row, M* column; no OCR word",
    ),
    (95, "surface_mstar"): (
        "2.1463663",
        "PDF p. 183, cone-surface row, M* column; OCR reads 21463669",
    ),
    (104, "shock_angle"): (
        "0.78592364",
        "PDF p. 201, shock row and p. 200 header; OCR reads 78592564 and 78592361",
    ),
    (106, "temperature"): (
        "1.1373037",
        "PDF p. 205, footer Ts/Too; OCR reads 11373097",
    ),
    (108, "density"): (
        "1.6849445",
        "PDF p. 209, footer rho_s/rho_oo; OCR '+168494_+01'",
    ),
    (112, "surface_mstar"): (
        "2.1161444",
        "PDF p. 216, cone-surface row, M* column; OCR '+2116144_+01'",
    ),
    (121, "temperature"): (
        "1.1204385",
        "PDF p. 234, footer Ts/Too; OCR '1204385+01'",
    ),
    (129, "post_shock_mstar"): (
        "2.0941342",
        "PDF p. 250, shock row, M* column; OCR '+20941542+0'",
    ),
    (140, "surface_mstar"): (
        "1.4062005",
        "PDF p. 270, cone-surface row, M* column; OCR '+140620C5+01'",
    ),
    (146, "density"): (
        "3.0369857",
        "PDF p. 283, footer rho_s/rho_oo; OCR '+30369857+01UPs/PoO'",
    ),
    (155, "surface_mstar"): (
        "1.0683100",
        "PDF p. 299, cone-surface row, M* column; no OCR word",
    ),
    (174, "post_shock_mstar"): (
        "1.3857733",
        "PDF p. 337, shock row, M* column; OCR reads 13857753",
    ),
    (190, "post_shock_mstar"): (
        "1.3328842",
        "PDF p. 369, shock row, M* column; OCR '+1.3328842+01'",
    ),
    (190, "shock_angle"): (
        "0.78882392",
        "PDF p. 369, shock row; OCR reads 78692392, p. 368 header reads 78882392",
    ),
    (192, "post_shock_mstar"): (
        "1.6442293",
        "PDF p. 373, shock row, M* column; OCR reads 15442293",
    ),
    (192, "surface_mstar"): (
        "1.6089479",
        "PDF p. 372, cone-surface row, M* column; OCR '+IS089479+01'",
    ),
    (196, "post_shock_mstar"): (
        "1.9088824",
        "PDF p. 381, shock row, M* column; OCR '+190S8824+01'",
    ),
    (196, "pressure"): (
        "9.1261264",
        "PDF p. 381, footer Ps/Poo; OCR '+g1261264+01'",
    ),
    (196, "surface_mstar"): (
        "1.8943970",
        "PDF p. 380, cone-surface row, M* column; OCR '+IS94397g+oi'",
    ),
    (206, "post_shock_mstar"): (
        "1.2749761",
        "PDF p. 401, shock row, M* column; no OCR word",
    ),
    (206, "pressure"): (
        "2.8063347",
        "PDF p. 401, footer Ps/Poo; OCR reads 29063347",
    ),
    (206, "shock_angle"): (
        "0.83913453",
        "PDF p. 401, shock row and p. 400 header; OCR '+8391345.5+00'",
    ),
    (206, "temperature"): (
        "1.3695501",
        "PDF p. 401, footer Ts/Too; OCR '3695501'",
    ),
    (208, "post_shock_mstar"): (
        "1.5879080",
        "PDF p. 405, shock row, M* column; no OCR word",
    ),
    (208, "pressure"): (
        "4.6226939",
        "PDF p. 405, footer Ps/Poo; OCR reads 43226939",
    ),
    (212, "density"): (
        "3.9856925",
        "PDF p. 413, footer rho_s/rho_oo; OCR '+398569_5+01\"'",
    ),
    (212, "post_shock_mstar"): (
        "1.8505089",
        "PDF p. 413, shock row, M* column; OCR reads 19505089",
    ),
    (212, "pressure"): (
        "10.492909",
        "PDF p. 413, footer Ps/Poo; OCR '+10492909+02I'",
    ),
    (212, "surface_mstar"): (
        "1.8350366",
        "PDF p. 412, cone-surface row, M* column; OCR reads 19350366",
    ),
}

# Unambiguous OCR digit confusions.  Other substitutions, such as G/g, S, or
# a leading 4 for "+", are not decoded automatically.
_DIGITS: Final = str.maketrans("IilLOo", "111100")
_TOKEN: Final = re.compile(
    r"(?<![0-9A-Za-z])[+]?([0-9IilLOo]{8})([+-])([0-9IilLOo]{2})(?![0-9A-Za-z])"
)
_WORD: Final = re.compile(
    r'<word xMin="([0-9.]+)" yMin="([0-9.]+)" xMax="[0-9.]+" yMax="[0-9.]+">'
    r"(.*?)</word>"
)
_NUMERIC_LIKE: Final = re.compile(r"[0-9IilLOoSsgGQuUD]{6,}")
_ROW_GAP: Final = 4.0
_COLUMN_HALF_WIDTH: Final = 25.0


@dataclass(frozen=True)
class _Cell:
    x: float
    y: float
    text: str

    @property
    def value(self) -> float | None:
        match = _TOKEN.search(self.text)
        if match is None:
            return None
        mantissa, sign, exponent = match.groups()
        digits = int(mantissa.translate(_DIGITS))
        return digits * 10.0 ** (int(sign + exponent.translate(_DIGITS)) - 8)


@cache
def _page_rows(page: int) -> tuple[tuple[_Cell, ...], ...]:
    """Return the page words grouped into rows ordered from top to bottom."""
    output = subprocess.check_output(
        ["pdftotext", "-f", str(page), "-l", str(page), "-bbox", str(PDF_PATH), "-"],
        text=True,
    )
    cells = sorted(
        (
            _Cell(float(x), float(y), html.unescape(text))
            for x, y, text in _WORD.findall(output)
        ),
        key=lambda cell: (cell.y, cell.x),
    )
    rows: list[list[_Cell]] = []
    for cell in cells:
        if rows and cell.y - rows[-1][-1].y <= _ROW_GAP:
            rows[-1].append(cell)
        else:
            rows.append([cell])
    return tuple(tuple(sorted(row, key=lambda cell: cell.x)) for row in rows)


def _numeric(row: tuple[_Cell, ...]) -> list[_Cell]:
    return [cell for cell in row if _NUMERIC_LIKE.search(cell.text)]


def _data_rows(page: int) -> list[tuple[_Cell, ...]]:
    return [row for row in _page_rows(page) if len(_numeric(row)) >= 3]


def _column_centres(rows: list[tuple[_Cell, ...]]) -> list[float]:
    """Return the median left edges of the five table columns."""
    complete = [_numeric(row) for row in rows if len(_numeric(row)) == 5]
    if not complete:
        raise RuntimeError("could not identify the table columns")
    return [
        sorted(cells[index].x for cells in complete)[len(complete) // 2]
        for index in range(5)
    ]


def _in_column(row: tuple[_Cell, ...], centre: float) -> _Cell | None:
    for cell in _numeric(row):
        if abs(cell.x - centre) <= _COLUMN_HALF_WIDTH:
            return cell
    return None


@dataclass
class _Table:
    cone_deg: float
    mach: float
    table: int
    page: int
    raw: dict[str, _Cell | None]
    header_mach: float
    header_shock_angle: float | None

    def location(self, field: str) -> str:
        return f"table {self.table}, PDF pages {self.page}-{self.page + 1}, {field}"


def _read_table(cone_deg: float, mach: float, table: int, page: int) -> _Table:
    cone = math.radians(cone_deg)
    first_rows = _data_rows(page)
    columns = _column_centres(first_rows)
    velocity_x = columns[1]
    # The header row above the table has no fifth column.  The first row that
    # does is the cone surface, where the flow-angle column repeats the cone
    # half-angle.
    surface_row = next(
        (row for row in first_rows if _in_column(row, columns[4]) is not None), None
    )
    if surface_row is None:
        raise RuntimeError(f"table {table}: no cone-surface row on PDF page {page}")
    angles = [
        cell.value
        for column in (columns[0], columns[2])
        if (cell := _in_column(surface_row, column)) is not None
        and cell.value is not None
    ]
    if angles and not any(abs(angle - cone) <= 1.0e-7 for angle in angles):
        raise RuntimeError(
            f"table {table}: first row on PDF page {page} is not at the cone surface"
        )
    # The header above the columns lists the computed free-stream Mach number
    # in column two and the shock angle in column four.
    header_rows = [row for row in _page_rows(page) if row[0].y < surface_row[0].y]
    mach_cells = [cell for row in header_rows if (cell := _in_column(row, columns[1]))]
    header_mach = mach_cells[-1].value if mach_cells else None
    if header_mach is None or abs(header_mach / mach - 1.0) > 1.0e-5:
        header_mach = mach
    angle_cells = [cell for row in header_rows if (cell := _in_column(row, columns[3]))]
    header_shock_angle = angle_cells[-1].value if angle_cells else None

    second = page + 1
    second_rows = _page_rows(second)
    second_data = _data_rows(second)
    angle_x2, velocity_x2 = _column_centres(second_data)[:2]
    shock_row = second_data[-1]
    shock_index = second_rows.index(shock_row)
    angle_cell = _in_column(shock_row, angle_x2)
    trailing = list(second_rows[shock_index + 1 :])
    if angle_cell is None and trailing:
        # The scan sometimes prints the last angle slightly below its row.
        angle_cell = _in_column(trailing[0], angle_x2)
        if angle_cell is not None:
            trailing = trailing[1:]
    footer = [
        cell
        for row in trailing
        for cell in _numeric(row)
        if abs(cell.x - velocity_x2) <= 3.0 * _COLUMN_HALF_WIDTH
    ]
    # The footer lists an unused integration parameter followed by p, rho
    # and T ratios.  A lost first cell leaves the three ratios in order;
    # misassignment is caught by the p = rho T check.
    ratios: list[_Cell | None] = list(footer[1:4] if len(footer) >= 4 else footer)
    ratios = [None] * (3 - len(ratios)) + ratios
    return _Table(
        cone_deg,
        mach,
        table,
        page,
        {
            "shock_angle": angle_cell,
            "surface_mstar": _in_column(surface_row, velocity_x),
            "post_shock_mstar": _in_column(shock_row, velocity_x2),
            "pressure": ratios[0],
            "density": ratios[1],
            "temperature": ratios[2],
        },
        header_mach,
        header_shock_angle,
    )


def _resolve(table: _Table, errors: list[str]) -> dict[str, float]:
    values: dict[str, float] = {}
    for field in FIELDS:
        cell = table.raw[field]
        decoded = None if cell is None else cell.value
        manual = MANUAL_CELLS.get((table.table, field))
        if manual is not None:
            if decoded is not None and f"{decoded:.8g}" == f"{float(manual[0]):.8g}":
                errors.append(f"{table.location(field)}: stale manual transcription")
            values[field] = float(manual[0])
        elif decoded is None:
            text = "missing" if cell is None else repr(cell.text)
            errors.append(f"{table.location(field)}: OCR cell {text} is undecodable")
        else:
            values[field] = decoded
    return values


def _mach(mstar: float) -> float:
    return math.sqrt(2.0 * mstar**2 / (GAMMA + 1.0 - (GAMMA - 1.0) * mstar**2))


def _check(table: _Table, values: dict[str, float], errors: list[str]) -> None:
    """Check relations that hold between cells of the same source table."""
    if len(values) != len(FIELDS):
        return
    header = table.header_shock_angle
    if (
        header is not None
        and (table.table, "shock_angle") not in MANUAL_CELLS
        and abs(header - values["shock_angle"]) > 2.0e-8
    ):
        errors.append(
            f"{table.location('shock_angle')}: row {values['shock_angle']:.8g} "
            f"differs from header {header:.8g}"
        )
    state = values["pressure"] / (values["density"] * values["temperature"])
    if abs(state - 1.0) > 5.0e-7:
        errors.append(f"{table.location('footer')}: p/(rho T) = {state:.9g}")
    surface_mach = _mach(values["surface_mstar"])
    total = (1.0 + 0.5 * (GAMMA - 1.0) * table.header_mach**2) / (
        1.0 + 0.5 * (GAMMA - 1.0) * surface_mach**2
    )
    if abs(total / values["temperature"] - 1.0) > 2.0e-6:
        errors.append(
            f"{table.location('surface_mstar')}: total-enthalpy temperature "
            f"{total:.8g} differs from footer {values['temperature']:.8g}"
        )
    # Rankine--Hugoniot state behind the tabulated shock angle.
    beta = values["shock_angle"]
    normal = table.header_mach * math.sin(beta)
    downstream_normal = math.sqrt(
        (1.0 + 0.5 * (GAMMA - 1.0) * normal**2)
        / (GAMMA * normal**2 - 0.5 * (GAMMA - 1.0))
    )
    deflection = math.atan(
        2.0
        / math.tan(beta)
        * (normal**2 - 1.0)
        / (table.header_mach**2 * (GAMMA + math.cos(2.0 * beta)) + 2.0)
    )
    post_shock = downstream_normal / math.sin(beta - deflection)
    tabulated = _mach(values["post_shock_mstar"])
    if abs(post_shock / tabulated - 1.0) > 1.0e-6:
        errors.append(
            f"{table.location('post_shock_mstar')}: Rankine-Hugoniot Mach "
            f"{post_shock:.8g} differs from tabulated {tabulated:.8g}"
        )


def _row(table: _Table, values: dict[str, float]) -> dict[str, str]:
    return {
        "table": str(table.table),
        "printed_page": str(table.page - PRINTED_PAGE_OFFSET),
        "cone_half_angle_deg": f"{table.cone_deg:.1f}",
        "upstream_mach": f"{table.mach:.1f}",
        "shock_angle_rad": f"{values['shock_angle']:.8g}",
        "surface_mach": f"{_mach(values['surface_mstar']):.9g}",
        "post_shock_mach": f"{_mach(values['post_shock_mstar']):.9g}",
        "surface_pressure_ratio": f"{values['pressure']:.8g}",
        "surface_density_ratio": f"{values['density']:.8g}",
        "surface_temperature_ratio": f"{values['temperature']:.8g}",
        "relative_tolerance": "0.0001",
    }


def main() -> None:
    if not PDF_PATH.exists():
        raise FileNotFoundError(f"download NASA SP-3004 to {PDF_PATH}")
    errors: list[str] = []
    rows: list[dict[str, str]] = []
    for cone_deg, mach, table_number, page in SOURCE_TABLES:
        table = _read_table(cone_deg, mach, table_number, page)
        values = _resolve(table, errors)
        _check(table, values, errors)
        if len(values) == len(FIELDS):
            rows.append(_row(table, values))
    if errors:
        raise RuntimeError(
            "unresolved NASA SP-3004 cells; transcribe them from the scan into "
            "MANUAL_CELLS:\n" + "\n".join(errors)
        )
    DIRECTORY.mkdir(parents=True, exist_ok=True)
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
        "gamma": GAMMA,
        "cone_angles_deg": sorted({entry[0] for entry in SOURCE_TABLES}),
        "mach_numbers": sorted({entry[1] for entry in SOURCE_TABLES}),
        "method": (
            "Cells selected by table position from the PDF text layer of "
            "visually verified table pages and decoded from the report's "
            "compact notation; no solver values are used. Damaged OCR cells "
            "are transcribed from the scan image. Each table is checked for "
            "header/shock-row agreement, p = rho T, total-enthalpy "
            "consistency of the surface temperature, and Rankine-Hugoniot "
            "consistency of the post-shock velocity."
        ),
        "manual_cells": [
            {"table": table, "field": field, "value": value, "location": location}
            for (table, field), (value, location) in sorted(MANUAL_CELLS.items())
        ],
        "command": "python docs/scripts/capture_sp3004_reference.py",
    }
    JSON_PATH.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
