"""Traceable engineering models for atmospheric and aerodynamic physics."""

from importlib.metadata import version

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
from aerophysics.expansion import prandtl_meyer_expansion
from aerophysics.flight import FlightCondition
from aerophysics.gas import AIR, PerfectGas
from aerophysics.protrusion import (
    ProtrusionDragResult,
    ProtrusionProfile,
    protrusion_drag,
)
from aerophysics.shocks import ShockBranch, normal_shock, oblique_shock

__version__ = version("aerophysics")

__all__ = [
    "AIR",
    "AtmosphereState",
    "BoundaryLayerRegime",
    "CompressibilityCorrection",
    "CompressibleBoundaryLayerProfileResult",
    "CompressibleVelocityTransformation",
    "FlightCondition",
    "PerfectGas",
    "ProtrusionDragResult",
    "ProtrusionProfile",
    "ShockBranch",
    "TemperatureVelocityRelation",
    "TransformedVelocityProfileResult",
    "TurbulentCorrelation",
    "__version__",
    "compressible_turbulent_boundary_layer_profile",
    "flat_plate_boundary_layer",
    "normal_shock",
    "oblique_shock",
    "prandtl_meyer_expansion",
    "protrusion_drag",
    "standard_atmosphere",
    "transform_compressible_velocity_profile",
]
