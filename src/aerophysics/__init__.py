"""Traceable engineering models for atmospheric and aerodynamic physics."""

from importlib.metadata import version

from aerophysics.gas import AIR, PerfectGas

__version__ = version("aerophysics")

__all__ = ["AIR", "PerfectGas", "__version__"]
