"""Calorically perfect-gas thermodynamics.

The default dry-air constants follow the U.S. Standard Atmosphere 1976. All
temperatures are in kelvin and all returned properties use SI units. Transport
models remain re-exported here for compatibility; new code should import them
from :mod:`aerophysics.transport`.
"""

from dataclasses import dataclass
from typing import overload

import numpy as np
from numpy.typing import ArrayLike

from aerophysics._array import FloatArray, FloatResult, as_float_array, return_float
from aerophysics.transport import (
    AIR_CONDUCTIVITY,
    AIR_VISCOSITY,
    SutherlandModel,
    USSAConductivityModel,
)


def _require_finite_positive(value: float, *, name: str) -> None:
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and greater than zero")


@dataclass(frozen=True, slots=True)
class PerfectGas:
    """Calorically perfect gas with constant ``R`` and ``gamma``.

    Parameters
    ----------
    specific_gas_constant:
        Specific gas constant in J/(kg K).
    heat_capacity_ratio:
        Ratio of specific heats, ``cp/cv``.
    """

    specific_gas_constant: float
    heat_capacity_ratio: float

    def __post_init__(self) -> None:
        _require_finite_positive(
            self.specific_gas_constant, name="specific_gas_constant"
        )
        if not np.isfinite(self.heat_capacity_ratio) or self.heat_capacity_ratio <= 1.0:
            raise ValueError("heat_capacity_ratio must be finite and greater than one")

    @property
    def cp(self) -> float:
        """Specific heat at constant pressure in J/(kg K)."""
        gamma = self.heat_capacity_ratio
        return gamma * self.specific_gas_constant / (gamma - 1.0)

    @property
    def cv(self) -> float:
        """Specific heat at constant volume in J/(kg K)."""
        gamma = self.heat_capacity_ratio
        return self.specific_gas_constant / (gamma - 1.0)

    @overload
    def speed_of_sound(self, temperature: float) -> float: ...

    @overload
    def speed_of_sound(self, temperature: ArrayLike) -> FloatArray: ...

    def speed_of_sound(self, temperature: ArrayLike) -> FloatResult:
        """Return the ideal-gas speed of sound in m/s."""
        values, scalar = as_float_array(temperature, name="temperature")
        if np.any(values <= 0.0):
            raise ValueError("temperature must be greater than zero")
        result = np.sqrt(self.heat_capacity_ratio * self.specific_gas_constant * values)
        return return_float(result, scalar=scalar)


AIR = PerfectGas(
    specific_gas_constant=8_314.32 / 28.9644,
    heat_capacity_ratio=1.4,
)
"""Calorically perfect dry air using U.S. Standard Atmosphere constants."""

__all__ = [
    "AIR",
    "AIR_CONDUCTIVITY",
    "AIR_VISCOSITY",
    "PerfectGas",
    "SutherlandModel",
    "USSAConductivityModel",
]
