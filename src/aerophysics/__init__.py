"""Traceable engineering models for atmospheric and aerodynamic physics."""

from importlib.metadata import version

from aerophysics.atmosphere import AtmosphereState, standard_atmosphere
from aerophysics.flight import FlightCondition
from aerophysics.gas import AIR, PerfectGas

__version__ = version("aerophysics")

__all__ = [
    "AIR",
    "AtmosphereState",
    "FlightCondition",
    "PerfectGas",
    "__version__",
    "standard_atmosphere",
]
