"""Package-level smoke tests."""

import os
import subprocess
import sys
from pathlib import Path
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


def test_version() -> None:
    assert __version__ == "0.5.0"


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


def test_wheel_bundles_docs_and_installed_gui_finds_them(tmp_path: Path) -> None:
    wheel_directory = tmp_path / "wheel"
    wheel_directory.mkdir()
    subprocess.run(
        [
            sys.executable,
            "-m",
            "build",
            "--wheel",
            "--outdir",
            str(wheel_directory),
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )

    wheels = sorted(wheel_directory.glob("aerophysics-*.whl"))
    assert len(wheels) == 1
    with ZipFile(wheels[0]) as archive:
        names = set(archive.namelist())
        required_docs = {
            "aerophysics/_docs/index.html",
            "aerophysics/_docs/getting_started/installation.html",
            "aerophysics/_docs/getting_started/quickstart.html",
            "aerophysics/_docs/getting_started/conventions.html",
            "aerophysics/_docs/guides/atmosphere_flight.html",
            "aerophysics/_docs/guides/compressible_flow.html",
            "aerophysics/_docs/guides/boundary_layers.html",
            "aerophysics/_docs/guides/vectorization_errors.html",
            "aerophysics/_docs/guides/gui.html",
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
        archive.extractall(tmp_path / "installed")

    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(tmp_path / "installed")
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "from aerophysics.gui.launcher import _documentation_directory; "
                "directory = _documentation_directory(); "
                "assert directory is not None; "
                "assert (directory / 'index.html').is_file(); "
                "assert (directory / 'models' / 'shock_waves.html').is_file(); "
                "assert (directory / 'guides' / 'gui.html').is_file(); "
                "assert not (directory / 'compressible_flow.html').exists(); "
                "print(directory)"
            ),
        ],
        cwd=tmp_path,
        env=environment,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "_docs" in result.stdout
