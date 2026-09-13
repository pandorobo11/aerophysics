from __future__ import annotations

import runpy
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path
from typing import cast

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.generated_assets
GENERATOR = PROJECT_ROOT / "docs/scripts/generate_viscosity_comparison.py"
SVG = PROJECT_ROOT / "docs/_static/viscosity_model_comparison.svg"
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


def test_bounded_error_format_hides_insignificant_cpu_noise() -> None:
    namespace = runpy.run_path(
        str(PROJECT_ROOT / "docs/scripts/_verification_common.py")
    )
    formatter = cast(Callable[[float, float], str], namespace["format_bounded_error"])
    epsilon = sys.float_info.epsilon
    assert formatter(2.0 * epsilon, 1.0e-12) == "< 5e-13"
    assert formatter(5.0e-13, 1.0e-12) == "5e-13"


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


def test_viscosity_comparison_assets_are_current() -> None:
    _assert_generator_check_passes(GENERATOR)


def test_viscosity_comparison_svg_has_accessible_labels() -> None:
    svg = SVG.read_text(encoding="utf-8")
    assert '<title id="svg-title">' in svg
    assert '<desc id="svg-desc">' in svg
    assert "Temperature (K)" in svg
    assert "Dynamic viscosity (Pa·s)" in svg
    assert "Relative difference (%)" in svg
    assert "frozen N₂/O₂/Ar/CO₂ dry-air composition" in svg


def test_standard_atmosphere_svgs_have_accessible_labels() -> None:
    for path in ATMOSPHERE_SVGS:
        svg = path.read_text(encoding="utf-8")
        assert 'role="img"' in svg
        assert 'aria-labelledby="svg-title svg-desc"' in svg
        assert '<title id="svg-title">' in svg
        assert '<desc id="svg-desc">' in svg
        assert "Geometric altitude (km)" in svg


def test_deterministic_verification_assets_are_current() -> None:
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
    assert "NIST physical-accuracy assessment" in fragment
    assert "not an acceptance test" in fragment
