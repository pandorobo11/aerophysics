"""Perfect-gas thermodynamics and air transport-property models.

The default dry-air and transport constants follow the U.S. Standard
Atmosphere 1976. All temperatures are in kelvin and all returned properties
use SI units.
"""

from dataclasses import dataclass
from typing import overload

import numpy as np
from numpy.typing import ArrayLike

from aerophysics._array import FloatArray, FloatResult, as_float_array, return_float


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


@dataclass(frozen=True, slots=True)
class SutherlandModel:
    """Sutherland dynamic-viscosity model for a gas."""

    reference_viscosity: float
    reference_temperature: float
    sutherland_temperature: float

    def __post_init__(self) -> None:
        _require_finite_positive(self.reference_viscosity, name="reference_viscosity")
        _require_finite_positive(
            self.reference_temperature, name="reference_temperature"
        )
        _require_finite_positive(
            self.sutherland_temperature, name="sutherland_temperature"
        )

    @overload
    def dynamic_viscosity(self, temperature: float) -> float: ...

    @overload
    def dynamic_viscosity(self, temperature: ArrayLike) -> FloatArray: ...

    def dynamic_viscosity(self, temperature: ArrayLike) -> FloatResult:
        """Return dynamic viscosity in Pa s."""
        values, scalar = as_float_array(temperature, name="temperature")
        if np.any(values <= 0.0):
            raise ValueError("temperature must be greater than zero")
        ratio = values / self.reference_temperature
        result = (
            self.reference_viscosity
            * ratio**1.5
            * (self.reference_temperature + self.sutherland_temperature)
            / (values + self.sutherland_temperature)
        )
        return return_float(result, scalar=scalar)


@dataclass(frozen=True, slots=True)
class USSAConductivityModel:
    """U.S. Standard Atmosphere 1976 air-conductivity correlation."""

    coefficient: float = 2.64638e-3
    additive_temperature: float = 245.4
    exponential_temperature: float = 12.0

    def __post_init__(self) -> None:
        _require_finite_positive(self.coefficient, name="coefficient")
        _require_finite_positive(self.additive_temperature, name="additive_temperature")
        _require_finite_positive(
            self.exponential_temperature, name="exponential_temperature"
        )

    @overload
    def thermal_conductivity(self, temperature: float) -> float: ...

    @overload
    def thermal_conductivity(self, temperature: ArrayLike) -> FloatArray: ...

    def thermal_conductivity(self, temperature: ArrayLike) -> FloatResult:
        """Return thermal conductivity in W/(m K)."""
        values, scalar = as_float_array(temperature, name="temperature")
        if np.any(values <= 0.0):
            raise ValueError("temperature must be greater than zero")
        denominator = values + self.additive_temperature * 10.0 ** (
            -self.exponential_temperature / values
        )
        result = self.coefficient * values**1.5 / denominator
        return return_float(result, scalar=scalar)


AIR = PerfectGas(specific_gas_constant=287.05287, heat_capacity_ratio=1.4)
"""Calorically perfect dry air using U.S. Standard Atmosphere constants."""

AIR_VISCOSITY = SutherlandModel(
    reference_viscosity=1.7894e-5,
    reference_temperature=288.15,
    sutherland_temperature=110.4,
)
"""Default Sutherland viscosity model for dry air."""

AIR_CONDUCTIVITY = USSAConductivityModel()
"""Default U.S. Standard Atmosphere thermal-conductivity model."""

__all__ = [
    "AIR",
    "AIR_CONDUCTIVITY",
    "AIR_VISCOSITY",
    "PerfectGas",
    "SutherlandModel",
    "USSAConductivityModel",
]
