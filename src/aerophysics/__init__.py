"""Traceable engineering models for atmospheric and aerodynamic physics."""

import tomllib
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from aerophysics.atmosphere import AtmosphereState, standard_atmosphere
from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    CompressibilityCorrection,
    TurbulentCorrelation,
    flat_plate_boundary_layer,
)
from aerophysics.boundary_layer_profile import (
    CompressibleBoundaryLayerProfileResult,
    CompressibleVelocityTransformation,
    TemperatureVelocityRelation,
    TransformedVelocityProfileResult,
    compressible_turbulent_boundary_layer_profile,
    transform_compressible_velocity_profile,
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
    "BeattieBridgemanGas",
    "BilligShockShapeResult",
    "BoundaryLayerRegime",
    "CompressibilityCorrection",
    "CompressibleBoundaryLayerProfileResult",
    "CompressibleVelocityTransformation",
    "DetachedShockComparisonResult",
    "DetachedShockGeometry",
    "DetachedShockModel",
    "DetachedShockStandoffResult",
    "FlightCondition",
    "HarmonicOscillatorGas",
    "IdealGasSpecies",
    "NASA7Polynomial",
    "NASA9Polynomial",
    "PerfectGas",
    "ProtrusionDragResult",
    "ProtrusionProfile",
    "ShockBranch",
    "TemperatureVelocityRelation",
    "ThermallyPerfectGas",
    "ThermodynamicState",
    "TransformedVelocityProfileResult",
    "TurbulentCorrelation",
    "VibrationalMode",
    "__version__",
    "billig_shock_shape",
    "compare_standoff_distances",
    "compressible_turbulent_boundary_layer_profile",
    "conical_shock",
    "flat_plate_boundary_layer",
    "normal_shock",
    "oblique_shock",
    "prandtl_meyer_expansion",
    "protrusion_drag",
    "seiff_standoff_distance",
    "seiff_standoff_distance_from_mach",
    "shock_standoff_distance",
    "standard_atmosphere",
    "transform_compressible_velocity_profile",
]
