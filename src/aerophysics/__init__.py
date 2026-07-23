"""Traceable engineering models for atmospheric and aerodynamic physics."""

from importlib.metadata import version

from aerophysics.atmosphere import AtmosphereState, standard_atmosphere
from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    CompressibilityCorrection,
    TurbulentCorrelation,
    flat_plate_boundary_layer,
)
from aerophysics.expansion import prandtl_meyer_expansion
from aerophysics.flight import FlightCondition
from aerophysics.gas import AIR, PerfectGas
from aerophysics.shocks import ShockBranch, normal_shock, oblique_shock

__version__ = version("aerophysics")

__all__ = [
    "AIR",
    "AtmosphereState",
    "BoundaryLayerRegime",
    "CompressibilityCorrection",
    "FlightCondition",
    "PerfectGas",
    "ShockBranch",
    "TurbulentCorrelation",
    "__version__",
    "flat_plate_boundary_layer",
    "normal_shock",
    "oblique_shock",
    "prandtl_meyer_expansion",
    "standard_atmosphere",
]
