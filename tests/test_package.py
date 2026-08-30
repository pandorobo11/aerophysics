"""Package-level smoke tests."""

import os
import shutil
import subprocess
import sys
from pathlib import Path
from textwrap import dedent
from zipfile import ZipFile

import pytest

from aerophysics import (
    AIR_BEATTIE_BRIDGEMAN,
    AIR_HARMONIC_OSCILLATOR,
    AIR_NASA7,
    AIR_NASA9,
    BeattieBridgemanGas,
    BoundaryLayerRegime,
    CompressibilityCorrection,
    CompressibleVelocityTransformation,
    DetachedShockGeometry,
    DetachedShockModel,
    HarmonicOscillatorGas,
    ShockBranch,
    TemperatureVelocityRelation,
    ThermallyPerfectGas,
    TurbulentCorrelation,
    __version__,
    billig_shock_shape,
    compare_standoff_distances,
    compressible_turbulent_boundary_layer_profile,
    conical_shock,
    flat_plate_boundary_layer,
    normal_shock,
    oblique_shock,
    prandtl_meyer_expansion,
    protrusion_drag,
    seiff_standoff_distance,
    seiff_standoff_distance_from_mach,
    shock_standoff_distance,
    transform_compressible_velocity_profile,
)
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


def test_version() -> None:
    assert __version__ == "0.6.0"


def test_public_diagnostics() -> None:
    assert issubclass(ModelRangeError, ValueError)
    assert issubclass(NoAttachedShockError, ValueError)
    assert issubclass(ApplicabilityWarning, UserWarning)


def test_primary_compressible_flow_api_is_exported() -> None:
    assert ShockBranch.WEAK.value == "weak"
    assert normal_shock(2.0).downstream_mach < 1.0
    assert oblique_shock(2.0, 0.1).downstream_mach > 1.0
    assert conical_shock(2.0, 0.1).surface_mach > 1.0
    assert prandtl_meyer_expansion(2.0, 0.1).downstream_mach > 2.0


def test_detached_shock_api_is_exported() -> None:
    geometry = DetachedShockGeometry.AXISYMMETRIC_SPHERE
    aw = shock_standoff_distance(4.0, 0.5, geometry=geometry)
    seiff = seiff_standoff_distance(4.0, 0.5)
    from_mach = seiff_standoff_distance_from_mach(4.0, 0.5)
    comparison = compare_standoff_distances(4.0, 0.5)
    shape = billig_shock_shape(4.0, 0.5, [-0.5, 0.0, 0.5], geometry=geometry)
    assert aw.model is DetachedShockModel.AMBROSIO_WORTMAN
    assert seiff.model is DetachedShockModel.SEIFF
    assert from_mach.density_ratio is not None
    assert from_mach.density_ratio > 1.0
    assert comparison.ambrosio_wortman.model is aw.model
    assert shape.model is DetachedShockModel.BILLIG


def test_thermally_perfect_air_api_is_exported() -> None:
    assert isinstance(AIR_NASA7, ThermallyPerfectGas)
    assert isinstance(AIR_NASA9, ThermallyPerfectGas)
    assert AIR_NASA9.heat_capacity_ratio(300.0) < 1.4
    assert isinstance(AIR_HARMONIC_OSCILLATOR, HarmonicOscillatorGas)
    assert isinstance(AIR_BEATTIE_BRIDGEMAN, BeattieBridgemanGas)


def test_primary_boundary_layer_api_is_exported() -> None:
    result = flat_plate_boundary_layer(
        1.0,
        10.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TURBULENT,
        turbulent_correlation=TurbulentCorrelation.POWER_LAW,
        compressibility_correction=CompressibilityCorrection.NONE,
    )
    assert result.drag_per_unit_width > 0.0


def test_compressible_boundary_layer_profile_api_is_exported() -> None:
    transformed = transform_compressible_velocity_profile(
        [0.0, 1.0],
        [0.0, 1.0],
        [1.0, 1.0],
        [1.0, 1.0],
        1.0,
        transformation=CompressibleVelocityTransformation.VAN_DRIEST,
    )
    assert transformed.transformed_velocity_plus[-1] == 1.0
    predicted = compressible_turbulent_boundary_layer_profile(
        [0.0, 0.05],
        300.0,
        1.0,
        300.0,
        0.05,
        85.0,
        transformation=CompressibleVelocityTransformation.VAN_DRIEST,
        temperature_velocity_relation=(
            TemperatureVelocityRelation.GENERALIZED_REYNOLDS_ANALOGY
        ),
        wall_temperature=250.0,
    )
    assert predicted.edge_velocity_ratio == pytest.approx(0.99)


def test_primary_protrusion_drag_api_is_exported() -> None:
    result = protrusion_drag(1.0, 0.01, 0.005, 10.0, 1.0, 0.02)
    assert result.direct_drag > 0.0


@pytest.mark.package_artifact
def test_wheel_clean_install_metadata_entry_point_and_docs(
    tmp_path: Path,
    artifact_directory: Path,
    runtime_constraints: Path,
) -> None:
    wheels = sorted(artifact_directory.glob(f"aerophysics-{__version__}-*.whl"))
    assert len(wheels) == 1
    wheel = wheels[0]
    with ZipFile(wheel) as archive:
        names = set(archive.namelist())
        required_docs = {
            "aerophysics/_docs/index.html",
            "aerophysics/_docs/getting_started/index.html",
            "aerophysics/_docs/getting_started/installation.html",
            "aerophysics/_docs/getting_started/quickstart.html",
            "aerophysics/_docs/getting_started/conventions.html",
            "aerophysics/_docs/guides/index.html",
            "aerophysics/_docs/guides/atmosphere_flight.html",
            "aerophysics/_docs/guides/compressible_flow.html",
            "aerophysics/_docs/guides/boundary_layers.html",
            "aerophysics/_docs/guides/vectorization_errors.html",
            "aerophysics/_docs/guides/gui.html",
            "aerophysics/_docs/models/index.html",
            "aerophysics/_docs/models/gas_and_atmosphere.html",
            "aerophysics/_docs/models/transport_properties.html",
            "aerophysics/_docs/models/thermochemistry.html",
            "aerophysics/_docs/models/isentropic_flow.html",
            "aerophysics/_docs/models/shock_waves.html",
            "aerophysics/_docs/models/expansion_waves.html",
            "aerophysics/_docs/models/flat_plate_boundary_layer.html",
            "aerophysics/_docs/models/compressible_velocity_transformations.html",
            "aerophysics/_docs/models/protrusion_drag.html",
            "aerophysics/_docs/models/flight_conditions.html",
            "aerophysics/_docs/models/unit_conversions.html",
            "aerophysics/_docs/verification/index.html",
            "aerophysics/_docs/verification/standard_atmosphere.html",
            "aerophysics/_docs/verification/compressible_flow.html",
            "aerophysics/_docs/verification/thermophysical.html",
            "aerophysics/_docs/verification/viscous_flow.html",
            "aerophysics/_docs/api/index.html",
            "aerophysics/_docs/api/thermophysical.html",
            "aerophysics/_docs/api/compressible_flow.html",
            "aerophysics/_docs/api/viscous_flow.html",
            "aerophysics/_docs/api/flight_units.html",
            "aerophysics/_docs/api/errors.html",
            "aerophysics/_docs/references.html",
        }
        legacy_docs = {
            "aerophysics/_docs/quickstart.html",
            "aerophysics/_docs/gas_and_atmosphere.html",
            "aerophysics/_docs/transport_properties.html",
            "aerophysics/_docs/thermochemistry.html",
            "aerophysics/_docs/compressible_flow.html",
            "aerophysics/_docs/boundary_layers.html",
            "aerophysics/_docs/compressible_velocity_transformations.html",
            "aerophysics/_docs/flight_conditions.html",
            "aerophysics/_docs/unit_conversions.html",
            "aerophysics/_docs/verification.html",
            "aerophysics/_docs/verification_compressible_flow.html",
            "aerophysics/_docs/verification_thermophysical.html",
            "aerophysics/_docs/verification_viscous_flow.html",
            "aerophysics/_docs/api.html",
        }
        assert required_docs <= names
        assert legacy_docs.isdisjoint(names)

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
        from aerophysics.gui.launcher import _documentation_directory

        package_path = Path(aerophysics.__file__).resolve()
        assert package_path.is_relative_to(Path(sys.prefix).resolve())
        assert aerophysics.__version__ == {__version__!r}
        assert normal_shock(2.0).downstream_mach < 1.0

        metadata = distribution("aerophysics")
        assert metadata.version == {__version__!r}
        requirements = [
            value.lower().replace('"', "'")
            for value in (metadata.requires or [])
        ]
        assert any(
            value.startswith("numpy") and ">=2.0" in value
            for value in requirements
        )
        assert any(
            value.startswith("scipy") and ">=1.14" in value
            for value in requirements
        )
        assert any(
            value.startswith("plotly")
            and ">=6.9" in value
            and "<7" in value
            and "extra == 'gui'" in value
            for value in requirements
        )
        assert any(
            value.startswith("streamlit")
            and ">=1.59" in value
            and "<2" in value
            and "extra == 'gui'" in value
            for value in requirements
        )

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
        assert (directory / "index.html").is_file()
        assert (directory / "models" / "shock_waves.html").is_file()
        assert (directory / "guides" / "gui.html").is_file()
        assert not (directory / "compressible_flow.html").exists()
        print(package_path)
        """
    )
    result = _run(
        [
            str(python),
            "-c",
            installed_check,
        ],
        cwd=outside_checkout,
        environment=_clean_environment(),
    )
    assert "site-packages" in result.stdout

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
        print(package_path)
        """
    )
    result = _run(
        [str(python), "-c", installed_check],
        cwd=outside_checkout,
        environment=_clean_environment(),
    )
    assert "site-packages" in result.stdout
