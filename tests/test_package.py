"""Package-level smoke tests."""

import os
import shutil
import subprocess
import sys
from importlib import import_module
from pathlib import Path
from textwrap import dedent

import pytest

import aerophysics
from aerophysics import __version__
from aerophysics.exceptions import (
    ApplicabilityWarning,
    ModelRangeError,
    NoAttachedShockError,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_DIRECTORY_ENVIRONMENT_VARIABLE = "AEROPHYSICS_TEST_ARTIFACT_DIRECTORY"


@pytest.fixture(scope="session")
def artifact_directory(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Return the externally built wheel and sdist directory."""
    configured = os.environ.get(ARTIFACT_DIRECTORY_ENVIRONMENT_VARIABLE)
    if configured is not None:
        directory = Path(configured)
        return directory if directory.is_absolute() else PROJECT_ROOT / directory

    directory = tmp_path_factory.mktemp("package-artifacts")
    _run(
        [
            sys.executable,
            "-m",
            "build",
            "--outdir",
            str(directory),
        ],
        cwd=PROJECT_ROOT,
    )
    return directory


def _run(
    command: list[str],
    *,
    cwd: Path,
    environment: dict[str, str] | None = None,
    timeout: float | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run a packaging command and retain actionable output on failure."""
    result = subprocess.run(
        command,
        cwd=cwd,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    if result.returncode != 0:
        pytest.fail(
            f"command failed with exit code {result.returncode}: {command!r}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )
    return result


def _clean_environment() -> dict[str, str]:
    environment = dict(os.environ)
    for name in (
        "AEROPHYSICS_DOCS_DIR",
        "AEROPHYSICS_DOCS_URL",
        "PYTHONHOME",
        "PYTHONPATH",
        "VIRTUAL_ENV",
    ):
        environment.pop(name, None)
    return environment


def _venv_python(venv: Path) -> Path:
    if sys.platform == "win32":
        return venv / "Scripts" / "python.exe"
    return venv / "bin" / "python"


def _venv_entry_point(venv: Path, name: str) -> Path:
    if sys.platform == "win32":
        return venv / "Scripts" / f"{name}.exe"
    return venv / "bin" / name


def _create_clean_venv(venv: Path, *, cwd: Path) -> Path:
    _run(
        [sys.executable, "-m", "venv", "--without-pip", str(venv)],
        cwd=cwd,
    )
    return _venv_python(venv)


def _write_runtime_constraints(path: Path) -> None:
    uv = shutil.which("uv")
    if uv is None:
        pytest.fail("uv executable is required for clean-install package tests")
    _run(
        [
            uv,
            "export",
            "--quiet",
            "--no-dev",
            "--all-extras",
            "--locked",
            "--no-emit-project",
            "--no-hashes",
            "--output-file",
            str(path),
        ],
        cwd=PROJECT_ROOT,
    )


@pytest.fixture(scope="session")
def runtime_constraints(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Export locked versions as constraints without preinstalling dependencies."""
    path = tmp_path_factory.mktemp("package-constraints") / "runtime.txt"
    _write_runtime_constraints(path)
    return path


def _install_artifact(
    artifact: str,
    python: Path,
    constraints: Path,
    *,
    cwd: Path,
) -> None:
    uv = shutil.which("uv")
    if uv is None:
        pytest.fail("uv executable is required for clean-install package tests")
    _run(
        [
            uv,
            "pip",
            "install",
            "--python",
            str(python),
            "--constraint",
            str(constraints),
            artifact,
        ],
        cwd=cwd,
        environment=_clean_environment(),
    )


def test_release_install_instructions_follow_package_version() -> None:
    wheel_url = (
        "https://github.com/pandorobo11/aerophysics/releases/download/"
        f"v{__version__}/aerophysics-{__version__}-py3-none-any.whl"
    )
    documented_paths = (
        "README.md",
        "README.ja.md",
        "docs/getting_started/installation.rst",
        "docs/guides/gui.rst",
    )
    for relative_path in documented_paths:
        content = (PROJECT_ROOT / relative_path).read_text(encoding="utf-8")
        assert wheel_url in content, relative_path


def test_public_diagnostics() -> None:
    assert issubclass(ModelRangeError, ValueError)
    assert issubclass(NoAttachedShockError, ValueError)
    assert issubclass(ApplicabilityWarning, UserWarning)


def test_public_models_are_exported_from_their_owning_modules() -> None:
    exports = {
        "shocks": ("ShockBranch", "normal_shock", "oblique_shock", "conical_shock"),
        "expansion": ("prandtl_meyer_expansion",),
        "detached_shock": (
            "DetachedShockGeometry",
            "DetachedShockModel",
            "shock_standoff_distance",
            "seiff_standoff_distance",
            "seiff_standoff_distance_from_mach",
            "compare_standoff_distances",
            "billig_shock_shape",
        ),
        "thermochemistry": ("AIR_NASA7", "AIR_NASA9", "ThermallyPerfectGas"),
        "real_gas": (
            "AIR_HARMONIC_OSCILLATOR",
            "AIR_BEATTIE_BRIDGEMAN",
            "HarmonicOscillatorGas",
            "BeattieBridgemanGas",
        ),
        "boundary_layer": (
            "BoundaryLayerRegime",
            "TurbulentCorrelation",
            "CompressibilityCorrection",
            "flat_plate_boundary_layer",
        ),
        "boundary_layer_profile": (
            "CompressibleVelocityTransformation",
            "TemperatureVelocityRelation",
            "transform_compressible_velocity_profile",
            "compressible_turbulent_boundary_layer_profile",
        ),
        "protrusion": ("protrusion_drag",),
    }
    for module_name, names in exports.items():
        module = import_module(f"aerophysics.{module_name}")
        for name in names:
            assert getattr(aerophysics, name) is getattr(module, name), name


@pytest.mark.package_artifact
def test_wheel_clean_install_metadata_entry_point_and_docs(
    tmp_path: Path,
    artifact_directory: Path,
    runtime_constraints: Path,
) -> None:
    wheels = sorted(artifact_directory.glob(f"aerophysics-{__version__}-*.whl"))
    assert len(wheels) == 1
    wheel = wheels[0]
    docs_source = PROJECT_ROOT / "docs"
    expected_pages = sorted(
        path.relative_to(docs_source).with_suffix(".html").as_posix()
        for path in docs_source.rglob("*.rst")
        if not any(part.startswith("_") for part in path.relative_to(docs_source).parts)
    )

    outside_checkout = tmp_path / "outside"
    outside_checkout.mkdir()
    venv = tmp_path / "wheel-venv"
    python = _create_clean_venv(venv, cwd=outside_checkout)
    _install_artifact(
        f"{wheel}[gui]",
        python,
        runtime_constraints,
        cwd=outside_checkout,
    )

    installed_check = dedent(
        f"""
        from importlib.metadata import distribution
        from pathlib import Path
        import sys

        import aerophysics
        from aerophysics import normal_shock
        from aerophysics.gui.documentation import DOCUMENTATION_TOPICS
        from aerophysics.gui.launcher import _documentation_directory

        package_path = Path(aerophysics.__file__).resolve()
        assert package_path.is_relative_to(Path(sys.prefix).resolve())
        assert aerophysics.__version__ == {__version__!r}
        assert normal_shock(2.0).downstream_mach < 1.0

        metadata = distribution("aerophysics")
        assert metadata.version == {__version__!r}
        entry_points = [
            entry
            for entry in metadata.entry_points
            if entry.group == "console_scripts"
            and entry.name == "aerophysics-gui"
        ]
        assert len(entry_points) == 1
        assert entry_points[0].value == "aerophysics.gui.launcher:main"

        directory = _documentation_directory()
        assert directory is not None
        assert directory == package_path.parent / "_docs"
        for relative_path in {expected_pages!r}:
            assert (directory / relative_path).is_file(), relative_path
        for relative_path in DOCUMENTATION_TOPICS.values():
            assert (directory / relative_path).is_file(), relative_path
        """
    )
    _run(
        [
            str(python),
            "-c",
            installed_check,
        ],
        cwd=outside_checkout,
        environment=_clean_environment(),
    )

    entry_point = _venv_entry_point(venv, "aerophysics-gui")
    assert entry_point.is_file()
    _run(
        [str(entry_point), "--help"],
        cwd=outside_checkout,
        environment=_clean_environment(),
        timeout=30.0,
    )


@pytest.mark.package_artifact
def test_sdist_installs_into_separate_clean_environment(
    tmp_path: Path,
    artifact_directory: Path,
    runtime_constraints: Path,
) -> None:
    sdists = sorted(artifact_directory.glob(f"aerophysics-{__version__}.tar.gz"))
    assert len(sdists) == 1

    outside_checkout = tmp_path / "outside"
    outside_checkout.mkdir()
    venv = tmp_path / "sdist-venv"
    python = _create_clean_venv(venv, cwd=outside_checkout)
    _install_artifact(
        str(sdists[0]),
        python,
        runtime_constraints,
        cwd=outside_checkout,
    )

    installed_check = dedent(
        f"""
        from importlib.metadata import version
        from pathlib import Path
        import sys

        import aerophysics
        from aerophysics import normal_shock

        package_path = Path(aerophysics.__file__).resolve()
        assert package_path.is_relative_to(Path(sys.prefix).resolve())
        assert version("aerophysics") == {__version__!r}
        assert normal_shock(2.0).downstream_mach < 1.0
        """
    )
    _run(
        [str(python), "-c", installed_check],
        cwd=outside_checkout,
        environment=_clean_environment(),
    )
