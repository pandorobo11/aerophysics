"""Traceable engineering models for atmospheric and aerodynamic physics."""

import tomllib
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from aerophysics.atmosphere import AtmosphereState, standard_atmosphere
from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    TurbulentModel,
    flat_plate_boundary_layer,
    flat_plate_boundary_layer_profile,
)
from aerophysics.boundary_layer_profile import (
    CompressibleBoundaryLayerProfileResult,
    CompressibleVelocityTransformation,
    TemperatureVelocityRelation,
    TransformedVelocityProfileResult,
    compressible_turbulent_boundary_layer_profile,
    transform_compressible_velocity_profile,
)
from aerophysics.characteristics import (
    AxisymmetricBodyFlowResult,
    BoundaryLayerInputs,
    CharacteristicNet,
    CharacteristicRow,
    axisymmetric_body_flow,
)
from aerophysics.detached_shock import (
    BilligShockShapeResult,
    DetachedShockComparisonResult,
    DetachedShockGeometry,
    DetachedShockModel,
    DetachedShockStandoffResult,
    billig_shock_shape,
    compare_standoff_distances,
    seiff_standoff_distance,
    seiff_standoff_distance_from_mach,
    shock_standoff_distance,
)
from aerophysics.expansion import prandtl_meyer_expansion
from aerophysics.flight import FlightCondition
from aerophysics.gas import AIR, PerfectGas
from aerophysics.heating import (
    AttachmentLineHeatingModel,
    AttachmentLineHeatingResult,
    StagnationHeatingModel,
    StagnationHeatingResult,
    SweptCylinderHeatingResult,
    sphere_stagnation_heat_flux,
    swept_cylinder_heat_flux,
    swept_cylinder_turbulent_heat_flux,
)
from aerophysics.marching_boundary_layer import (
    BoundaryLayerEdge,
    MarchingBoundaryLayerResult,
    marching_boundary_layer,
)
from aerophysics.protrusion import (
    ProtrusionDragResult,
    ProtrusionProfile,
    protrusion_drag,
)
from aerophysics.real_gas import (
    AIR_BEATTIE_BRIDGEMAN,
    AIR_HARMONIC_OSCILLATOR,
    BeattieBridgemanGas,
    HarmonicOscillatorGas,
    ThermodynamicState,
    VibrationalMode,
)
from aerophysics.shocks import ShockBranch, conical_shock, normal_shock, oblique_shock
from aerophysics.thermochemistry import (
    AIR_NASA7,
    AIR_NASA9,
    IdealGasSpecies,
    NASA7Polynomial,
    NASA9Polynomial,
    ThermallyPerfectGas,
)


def _package_version() -> str:
    try:
        return version("aerophysics")
    except PackageNotFoundError:
        pyproject = Path(__file__).resolve().parents[2] / "pyproject.toml"
        if not pyproject.is_file():
            return "unknown"
        with pyproject.open("rb") as stream:
            project = tomllib.load(stream)["project"]
        return str(project["version"])


__version__ = _package_version()

__all__ = [
    "AIR",
    "AIR_BEATTIE_BRIDGEMAN",
    "AIR_HARMONIC_OSCILLATOR",
    "AIR_NASA7",
    "AIR_NASA9",
    "AtmosphereState",
    "AttachmentLineHeatingModel",
    "AttachmentLineHeatingResult",
    "AxisymmetricBodyFlowResult",
    "BeattieBridgemanGas",
    "BilligShockShapeResult",
    "BoundaryLayerEdge",
    "BoundaryLayerInputs",
    "BoundaryLayerRegime",
    "CharacteristicNet",
    "CharacteristicRow",
    "CompressibleBoundaryLayerProfileResult",
    "CompressibleVelocityTransformation",
    "DetachedShockComparisonResult",
    "DetachedShockGeometry",
    "DetachedShockModel",
    "DetachedShockStandoffResult",
    "FlightCondition",
    "HarmonicOscillatorGas",
    "IdealGasSpecies",
    "MarchingBoundaryLayerResult",
    "NASA7Polynomial",
    "NASA9Polynomial",
    "PerfectGas",
    "ProtrusionDragResult",
    "ProtrusionProfile",
    "ShockBranch",
    "StagnationHeatingModel",
    "StagnationHeatingResult",
    "SweptCylinderHeatingResult",
    "TemperatureVelocityRelation",
    "ThermallyPerfectGas",
    "ThermodynamicState",
    "TransformedVelocityProfileResult",
    "TurbulentModel",
    "VibrationalMode",
    "__version__",
    "axisymmetric_body_flow",
    "billig_shock_shape",
    "compare_standoff_distances",
    "compressible_turbulent_boundary_layer_profile",
    "conical_shock",
    "flat_plate_boundary_layer",
    "flat_plate_boundary_layer_profile",
    "marching_boundary_layer",
    "normal_shock",
    "oblique_shock",
    "prandtl_meyer_expansion",
    "protrusion_drag",
    "seiff_standoff_distance",
    "seiff_standoff_distance_from_mach",
    "shock_standoff_distance",
    "sphere_stagnation_heat_flux",
    "standard_atmosphere",
    "swept_cylinder_heat_flux",
    "swept_cylinder_turbulent_heat_flux",
    "transform_compressible_velocity_profile",
]
