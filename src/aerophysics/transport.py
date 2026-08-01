"""Temperature-dependent gas transport-property models.

All temperatures are in kelvin. Dynamic viscosity is returned in Pa s and
thermal conductivity is returned in W/(m K). The built-in dry-air constants
use the same frozen composition as the thermochemistry models.
"""

import warnings
from dataclasses import dataclass
from typing import Protocol, overload

import numpy as np
from numpy.typing import ArrayLike

from aerophysics._array import FloatArray, FloatResult, as_float_array, return_float
from aerophysics._nasa_data import DRY_AIR_MOLE_FRACTIONS, NASA9_DATA
from aerophysics.exceptions import ApplicabilityWarning


class DynamicViscosityModel(Protocol):
    """Structural interface for temperature-dependent viscosity models."""

    @overload
    def dynamic_viscosity(self, temperature: float) -> float: ...

    @overload
    def dynamic_viscosity(self, temperature: ArrayLike) -> FloatArray: ...

    def dynamic_viscosity(self, temperature: ArrayLike) -> FloatResult:
        """Return dynamic viscosity in Pa s."""
        ...


def _require_finite_positive(value: float, *, name: str) -> None:
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and greater than zero")


def _normalise_temperature_range(
    temperature_range: tuple[float, float],
) -> tuple[float, float]:
    try:
        minimum, maximum = (float(value) for value in temperature_range)
    except (TypeError, ValueError) as error:
        raise ValueError(
            "temperature_range must contain two real numeric values"
        ) from error
    if not np.isfinite(minimum) or not np.isfinite(maximum):
        raise ValueError("temperature_range must contain only finite values")
    if minimum <= 0.0 or maximum <= minimum:
        raise ValueError(
            "temperature_range must be positive and strictly increasing"
        )
    return minimum, maximum


def _validated_temperature(
    temperature: ArrayLike,
    *,
    temperature_range: tuple[float, float] | None = None,
) -> tuple[FloatArray, bool]:
    values, scalar = as_float_array(temperature, name="temperature")
    if np.any(values <= 0.0):
        raise ValueError("temperature must be greater than zero")
    if temperature_range is not None:
        minimum, maximum = temperature_range
        if np.any((values < minimum) | (values > maximum)):
            warnings.warn(
                (
                    f"temperature is outside the fitted range {minimum:g}--"
                    f"{maximum:g} K; the viscosity correlation was extrapolated"
                ),
                ApplicabilityWarning,
                stacklevel=3,
            )
    return values, scalar


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
        values, scalar = _validated_temperature(temperature)
        ratio = values / self.reference_temperature
        result = (
            self.reference_viscosity
            * ratio**1.5
            * (self.reference_temperature + self.sutherland_temperature)
            / (values + self.sutherland_temperature)
        )
        return return_float(result, scalar=scalar)


@dataclass(frozen=True, slots=True)
class KeyesModel:
    """Keyes dynamic-viscosity correlation for a gas."""

    coefficient: float
    additive_temperature: float
    exponential_temperature: float
    temperature_range: tuple[float, float] = (79.0, 1845.0)

    def __post_init__(self) -> None:
        _require_finite_positive(self.coefficient, name="coefficient")
        _require_finite_positive(
            self.additive_temperature, name="additive_temperature"
        )
        _require_finite_positive(
            self.exponential_temperature, name="exponential_temperature"
        )
        object.__setattr__(
            self,
            "temperature_range",
            _normalise_temperature_range(self.temperature_range),
        )

    @overload
    def dynamic_viscosity(self, temperature: float) -> float: ...

    @overload
    def dynamic_viscosity(self, temperature: ArrayLike) -> FloatArray: ...

    def dynamic_viscosity(self, temperature: ArrayLike) -> FloatResult:
        """Return dynamic viscosity in Pa s."""
        values, scalar = _validated_temperature(
            temperature, temperature_range=self.temperature_range
        )
        denominator = values + self.additive_temperature * 10.0 ** (
            -self.exponential_temperature / values
        )
        result = self.coefficient * values**1.5 / denominator
        return return_float(result, scalar=scalar)


@dataclass(frozen=True, slots=True)
class BlottnerModel:
    """Blottner curve fit for the viscosity of one gas species."""

    a: float
    b: float
    c: float
    temperature_range: tuple[float, float] = (1000.0, 30_000.0)

    def __post_init__(self) -> None:
        if not np.all(np.isfinite((self.a, self.b, self.c))):
            raise ValueError("a, b, and c must be finite")
        object.__setattr__(
            self,
            "temperature_range",
            _normalise_temperature_range(self.temperature_range),
        )

    @overload
    def dynamic_viscosity(self, temperature: float) -> float: ...

    @overload
    def dynamic_viscosity(self, temperature: ArrayLike) -> FloatArray: ...

    def dynamic_viscosity(self, temperature: ArrayLike) -> FloatResult:
        """Return species dynamic viscosity in Pa s."""
        values, scalar = _validated_temperature(
            temperature, temperature_range=self.temperature_range
        )
        logarithm = np.log(values)
        result = 0.1 * np.exp((self.a * logarithm + self.b) * logarithm + self.c)
        return return_float(result, scalar=scalar)


@dataclass(frozen=True, slots=True)
class WilkeMixtureViscosityModel:
    """Frozen-composition gas viscosity using Wilke's mixing rule."""

    component_models: tuple[DynamicViscosityModel, ...]
    molar_masses: tuple[float, ...]
    mole_fractions: tuple[float, ...]

    def __post_init__(self) -> None:
        count = len(self.component_models)
        if count == 0:
            raise ValueError("component_models must contain at least one model")
        if len(self.molar_masses) != count or len(self.mole_fractions) != count:
            raise ValueError(
                "component_models, molar_masses, and mole_fractions must have "
                "the same length"
            )
        if any(
            not callable(getattr(model, "dynamic_viscosity", None))
            for model in self.component_models
        ):
            raise TypeError(
                "component_models must contain dynamic-viscosity models"
            )
        try:
            masses = tuple(float(value) for value in self.molar_masses)
            fractions = tuple(float(value) for value in self.mole_fractions)
        except (TypeError, ValueError) as error:
            raise ValueError(
                "molar_masses and mole_fractions must contain real numeric values"
            ) from error
        if not np.all(np.isfinite(masses)) or any(value <= 0.0 for value in masses):
            raise ValueError("molar_masses must be finite and greater than zero")
        if not np.all(np.isfinite(fractions)) or any(
            value <= 0.0 for value in fractions
        ):
            raise ValueError("mole_fractions must be finite and greater than zero")
        total = sum(fractions)
        if not np.isclose(total, 1.0, rtol=0.0, atol=1.0e-12):
            raise ValueError("mole_fractions must sum to one")
        object.__setattr__(self, "molar_masses", masses)
        object.__setattr__(
            self, "mole_fractions", tuple(value / total for value in fractions)
        )

    @overload
    def dynamic_viscosity(self, temperature: float) -> float: ...

    @overload
    def dynamic_viscosity(self, temperature: ArrayLike) -> FloatArray: ...

    def dynamic_viscosity(self, temperature: ArrayLike) -> FloatResult:
        """Return mixture dynamic viscosity in Pa s."""
        values, scalar = _validated_temperature(temperature)
        component_viscosities = tuple(
            np.asarray(model.dynamic_viscosity(values), dtype=np.float64)
            for model in self.component_models
        )
        result = np.zeros_like(values, dtype=np.float64)
        for i, viscosity_i in enumerate(component_viscosities):
            denominator = np.zeros_like(values, dtype=np.float64)
            for j, viscosity_j in enumerate(component_viscosities):
                phi = (
                    1.0
                    + np.sqrt(viscosity_i / viscosity_j)
                    * (self.molar_masses[j] / self.molar_masses[i]) ** 0.25
                ) ** 2 / np.sqrt(
                    8.0 * (1.0 + self.molar_masses[i] / self.molar_masses[j])
                )
                denominator += self.mole_fractions[j] * phi
            result += self.mole_fractions[i] * viscosity_i / denominator
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
        values, scalar = _validated_temperature(temperature)
        denominator = values + self.additive_temperature * 10.0 ** (
            -self.exponential_temperature / values
        )
        result = self.coefficient * values**1.5 / denominator
        return return_float(result, scalar=scalar)


AIR_VISCOSITY = SutherlandModel(
    reference_viscosity=1.7894e-5,
    reference_temperature=288.15,
    sutherland_temperature=110.4,
)
"""Default Sutherland viscosity model for dry air."""

AIR_KEYES_VISCOSITY = KeyesModel(
    coefficient=1.488e-6,
    additive_temperature=122.1,
    exponential_temperature=5.0,
)
"""Keyes viscosity model for dry air."""

_DRY_AIR_SPECIES = ("N2", "O2", "Ar", "CO2")
_DRY_AIR_BLOTTNER_COEFFICIENTS = (
    (0.0268142, 0.3177838, -11.3155513),
    (0.0449290, -0.0826158, -9.2019475),
    (-0.02201, 1.010, -13.42),
    (-0.041372, 1.3293, -15.016),
)
_raw_air_fractions = tuple(
    DRY_AIR_MOLE_FRACTIONS[name] for name in _DRY_AIR_SPECIES
)
_air_fraction_total = sum(_raw_air_fractions)

AIR_BLOTTNER_VISCOSITY = WilkeMixtureViscosityModel(
    component_models=tuple(
        BlottnerModel(*coefficients)
        for coefficients in _DRY_AIR_BLOTTNER_COEFFICIENTS
    ),
    molar_masses=tuple(NASA9_DATA[name][0] for name in _DRY_AIR_SPECIES),
    mole_fractions=tuple(value / _air_fraction_total for value in _raw_air_fractions),
)
"""Frozen dry-air Blottner viscosities combined with Wilke's rule."""

AIR_CONDUCTIVITY = USSAConductivityModel()
"""Default U.S. Standard Atmosphere thermal-conductivity model."""

__all__ = [
    "AIR_BLOTTNER_VISCOSITY",
    "AIR_CONDUCTIVITY",
    "AIR_KEYES_VISCOSITY",
    "AIR_VISCOSITY",
    "BlottnerModel",
    "DynamicViscosityModel",
    "KeyesModel",
    "SutherlandModel",
    "USSAConductivityModel",
    "WilkeMixtureViscosityModel",
]
