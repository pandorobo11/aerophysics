"""Build source-equation reference values for smooth flat plates."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / "tests/reference_data/viscous_flow"
CSV_PATH = DIRECTORY / "flat_plate_source_equations.csv"
JSON_PATH = DIRECTORY / "flat_plate_source_equations.json"
REYNOLDS_NUMBERS = (1.0e5, 2.0e5, 5.0e5, 1.0e6, 2.0e6, 5.0e6, 1.0e7, 1.0e8, 1.0e9)


def main() -> None:
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    rows = []
    for reynolds in REYNOLDS_NUMBERS:
        rows.append(
            {
                "reynolds_number": f"{reynolds:.8g}",
                "blasius_delta_over_x": f"{5.0 / math.sqrt(reynolds):.15g}",
                "blasius_displacement_over_x": f"{1.7208 / math.sqrt(reynolds):.15g}",
                "blasius_momentum_over_x": f"{0.664 / math.sqrt(reynolds):.15g}",
                "blasius_local_cf": f"{0.664 / math.sqrt(reynolds):.15g}",
                "blasius_average_cf": f"{1.328 / math.sqrt(reynolds):.15g}",
                "power_delta_over_x": f"{0.37 * reynolds**-0.2:.15g}",
                "power_local_cf": f"{0.0592 * reynolds**-0.2:.15g}",
                "power_average_cf": f"{0.074 * reynolds**-0.2:.15g}",
                "schlichting_average_cf": (
                    f"{0.455 / math.log10(reynolds) ** 2.58:.15g}"
                ),
            }
        )
    with CSV_PATH.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "sources": [
            "Blasius laminar flat-plate similarity solution",
            "one-fifth-power turbulent flat-plate correlation",
            "Schlichting smooth-plate average skin-friction correlation",
        ],
        "role": "independent direct evaluation of published correlation equations",
        "command": "python docs/scripts/build_viscous_reference.py",
    }
    JSON_PATH.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
