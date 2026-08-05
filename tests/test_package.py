"""Package-level smoke tests."""

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
    HarmonicOscillatorGas,
    ShockBranch,
    TemperatureVelocityRelation,
    ThermallyPerfectGas,
    TurbulentCorrelation,
    __version__,
    compressible_turbulent_boundary_layer_profile,
    conical_shock,
    flat_plate_boundary_layer,
    normal_shock,
    oblique_shock,
    prandtl_meyer_expansion,
    protrusion_drag,
    transform_compressible_velocity_profile,
)
from aerophysics.exceptions import (
    ApplicabilityWarning,
    ModelRangeError,
    NoAttachedShockError,
)


def test_version() -> None:
    assert __version__ == "0.3.0"


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
