"""Optional local web GUI for :mod:`aerophysics`."""

from aerophysics.gui.config import (
    CONFIG_SCHEMA_VERSION,
    ConfigurationError,
    dump_configuration,
    load_configuration,
)
from aerophysics.gui.units import UnitPreferences

__all__ = [
    "CONFIG_SCHEMA_VERSION",
    "ConfigurationError",
    "UnitPreferences",
    "dump_configuration",
    "load_configuration",
]
