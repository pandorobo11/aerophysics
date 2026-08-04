from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from aerophysics.transport import (
    AIR_BLOTTNER_VISCOSITY,
    AIR_KEYES_VISCOSITY,
    AIR_VISCOSITY,
    DynamicViscosityModel,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
GENERATOR = PROJECT_ROOT / "docs/scripts/generate_viscosity_comparison.py"
TABLE_FRAGMENT = PROJECT_ROOT / "docs/_generated/viscosity_model_comparison.rst"
SVG = PROJECT_ROOT / "docs/_static/viscosity_model_comparison.svg"


def _expected_row(
    temperature: float, candidate: DynamicViscosityModel
) -> str:
    baseline_value = AIR_VISCOSITY.dynamic_viscosity(temperature)
    candidate_value = candidate.dynamic_viscosity(temperature)
    relative_difference = (candidate_value / baseline_value - 1.0) * 100.0
    return "\n".join(
        [
            f"   * - {temperature:g}",
            f"     - {baseline_value:.7e}",
            f"     - {candidate_value:.7e}",
            f"     - {relative_difference:+.3f}",
        ]
    )


def test_viscosity_comparison_assets_are_current() -> None:
    subprocess.run(
        [sys.executable, str(GENERATOR), "--check"],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )


def test_viscosity_comparison_tables_match_public_models() -> None:
    fragment = TABLE_FRAGMENT.read_text(encoding="utf-8")
    for temperature in (79.0, 100.0, 300.0, 1000.0, 1845.0):
        assert _expected_row(temperature, AIR_KEYES_VISCOSITY) in fragment
    for temperature in (1000.0, 1845.0, 5000.0, 10000.0, 30000.0):
        assert _expected_row(temperature, AIR_BLOTTNER_VISCOSITY) in fragment


def test_viscosity_comparison_svg_has_accessible_labels() -> None:
    svg = SVG.read_text(encoding="utf-8")
    assert "<title id=\"svg-title\">" in svg
    assert "<desc id=\"svg-desc\">" in svg
    assert "Temperature (K)" in svg
    assert "Dynamic viscosity (Pa·s)" in svg
    assert "Relative difference (%)" in svg
    assert "frozen N₂/O₂/Ar/CO₂ dry-air composition" in svg
