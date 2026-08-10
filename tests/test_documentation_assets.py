from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

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
ATMOSPHERE_FRAGMENT = (
    PROJECT_ROOT / "docs/_generated/standard_atmosphere_validation.rst"
)
ATMOSPHERE_SVGS = (
    PROJECT_ROOT / "docs/_static/standard_atmosphere_profiles.svg",
    PROJECT_ROOT / "docs/_static/standard_atmosphere_comparison.svg",
)
VERIFICATION_GENERATOR = PROJECT_ROOT / "docs/scripts/generate_verification.py"
VERIFICATION_SVGS = (
    PROJECT_ROOT / "docs/_static/compressible_flow_profiles.svg",
    PROJECT_ROOT / "docs/_static/compressible_flow_differences.svg",
    PROJECT_ROOT / "docs/_static/thermophysical_properties.svg",
    PROJECT_ROOT / "docs/_static/thermophysical_transport_accuracy.svg",
    PROJECT_ROOT / "docs/_static/viscous_skin_friction.svg",
    PROJECT_ROOT / "docs/_static/protrusion_shielding.svg",
)


def _expected_row(temperature: float, candidate: DynamicViscosityModel) -> str:
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


def _assert_generator_check_passes(generator: Path) -> None:
    result = subprocess.run(
        [sys.executable, str(generator), "--check"],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
    )
    if result.returncode:
        output = "\n".join(
            part for part in (result.stdout, result.stderr) if part
        ).strip()
        pytest.fail(
            f"{generator.relative_to(PROJECT_ROOT)} reported stale generated assets"
            f"\n{output}"
        )


@pytest.mark.generated_assets
def test_viscosity_comparison_assets_are_current() -> None:
    _assert_generator_check_passes(GENERATOR)


def test_viscosity_comparison_tables_match_public_models() -> None:
    fragment = TABLE_FRAGMENT.read_text(encoding="utf-8")
    for temperature in (79.0, 100.0, 300.0, 1000.0, 1845.0):
        assert _expected_row(temperature, AIR_KEYES_VISCOSITY) in fragment
    for temperature in (1000.0, 1845.0, 5000.0, 10000.0, 30000.0):
        assert _expected_row(temperature, AIR_BLOTTNER_VISCOSITY) in fragment


def test_viscosity_comparison_svg_has_accessible_labels() -> None:
    svg = SVG.read_text(encoding="utf-8")
    assert '<title id="svg-title">' in svg
    assert '<desc id="svg-desc">' in svg
    assert "Temperature (K)" in svg
    assert "Dynamic viscosity (Pa·s)" in svg
    assert "Relative difference (%)" in svg
    assert "frozen N₂/O₂/Ar/CO₂ dry-air composition" in svg


def test_standard_atmosphere_tables_record_the_computed_result() -> None:
    fragment = ATMOSPHERE_FRAGMENT.read_text(encoding="utf-8")
    assert "**Verified with observations**" in fragment
    assert "Official-table comparison summary" in fragment
    assert "fluids 1.3.1 comparison summary" in fragment
    assert "Physical and mathematical invariant summary" in fragment
    assert "9.1690e-05" in fragment


def test_standard_atmosphere_svgs_have_accessible_labels() -> None:
    for path in ATMOSPHERE_SVGS:
        svg = path.read_text(encoding="utf-8")
        assert 'role="img"' in svg
        assert 'aria-labelledby="svg-title svg-desc"' in svg
        assert '<title id="svg-title">' in svg
        assert '<desc id="svg-desc">' in svg
        assert "Geometric altitude (km)" in svg


@pytest.mark.generated_assets
def test_all_verification_assets_are_current() -> None:
    _assert_generator_check_passes(VERIFICATION_GENERATOR)


def test_new_verification_svgs_have_accessible_labels() -> None:
    for path in VERIFICATION_SVGS:
        svg = path.read_text(encoding="utf-8")
        assert 'role="img"' in svg
        assert 'aria-labelledby="title desc"' in svg
        assert '<title id="title">' in svg
        assert '<desc id="desc">' in svg
        assert "[-]" in svg or "[K]" in svg


def test_thermophysical_record_separates_verification_and_accuracy() -> None:
    fragment = (
        PROJECT_ROOT / "docs/_generated/thermophysical_validation.rst"
    ).read_text(encoding="utf-8")
    assert "**Overall status: Verified.**" in fragment
    assert "Primary transport references" in fragment
    assert "NIST physical-accuracy assessment" in fragment
    assert "not an acceptance test" in fragment
