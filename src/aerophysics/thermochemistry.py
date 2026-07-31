"""NASA-polynomial thermochemistry for frozen-composition ideal-gas mixtures.

All public calculation methods use SI units. Temperatures are in kelvin,
pressures are in pascals, molar properties are per mole, and mixture-specific
properties are per kilogram.
"""

import warnings
from dataclasses import dataclass
from itertools import pairwise
from typing import cast

import numpy as np
from numpy.typing import ArrayLike

from aerophysics._array import FloatArray, FloatResult, as_float_array, return_float
from aerophysics._nasa_data import (
    DRY_AIR_MOLE_FRACTIONS,
    NASA7_DATA,
    NASA9_DATA,
    PolynomialData,
)
from aerophysics.exceptions import ApplicabilityWarning, ModelRangeError

UNIVERSAL_GAS_CONSTANT = 8.31446261815324
"""Universal gas constant in J/(mol K)."""

STANDARD_PRESSURE = 100_000.0
"""NASA thermodynamic-data standard pressure in Pa."""

type _CoefficientRows = tuple[tuple[float, ...], ...]


def _normalise_polynomial_data(
    temperature_ranges: tuple[float, ...],
    coefficients: _CoefficientRows,
    *,
    coefficient_count: int,
) -> tuple[tuple[float, ...], _CoefficientRows]:
    try:
        ranges = tuple(float(value) for value in temperature_ranges)
        rows = tuple(tuple(float(value) for value in row) for row in coefficients)
    except (TypeError, ValueError) as error:
        raise ValueError("polynomial data must contain real numeric values") from error

    if not rows:
        raise ValueError("coefficients must contain at least one temperature region")
    if len(ranges) != len(rows) + 1:
        raise ValueError(
            "temperature_ranges must contain one more value than coefficients"
        )
    if not np.all(np.isfinite(ranges)) or any(
        not np.all(np.isfinite(row)) for row in rows
    ):
        raise ValueError("polynomial data must contain only finite values")
    if ranges[0] <= 0.0 or any(upper <= lower for lower, upper in pairwise(ranges)):
        raise ValueError("temperature_ranges must be positive and strictly increasing")
    if any(len(row) != coefficient_count for row in rows):
        raise ValueError(
            f"each coefficient row must contain {coefficient_count} values"
        )
    return ranges, rows


def _validated_temperature(
    temperature: ArrayLike,
    temperature_range: tuple[float, float],
    *,
    allow_extrapolation: bool,
    warn: bool,
) -> tuple[FloatArray, bool]:
    values, scalar = as_float_array(temperature, name="temperature")
    if np.any(values <= 0.0):
        raise ValueError("temperature must be greater than zero")

    minimum, maximum = temperature_range
    outside = np.any((values < minimum) | (values > maximum))
    if outside and not allow_extrapolation:
        raise ModelRangeError(
            f"temperature must be within {minimum:g} K and {maximum:g} K"
        )
    if outside and warn:
        warnings.warn(
            (
                f"temperature is outside the fitted range {minimum:g}--"
                f"{maximum:g} K; the nearest polynomial region was extrapolated"
            ),
            ApplicabilityWarning,
            stacklevel=3,
        )
    return values, scalar


def _region_coefficients(
    temperature: FloatArray,
    temperature_ranges: tuple[float, ...],
    coefficients: _CoefficientRows,
) -> FloatArray:
    # ``side='left'`` makes a shared boundary belong to the lower region.
    indices = np.searchsorted(temperature_ranges[1:], temperature, side="left")
    indices = np.clip(indices, 0, len(coefficients) - 1)
    rows = np.asarray(coefficients, dtype=np.float64)
    return rows[indices]


@dataclass(frozen=True, slots=True)
class NASA7Polynomial:
    """Piecewise NASA seven-coefficient ideal-gas polynomial.

    ``temperature_ranges`` contains the region boundaries in kelvin.
    ``coefficients`` contains one ``(a1, ..., a7)`` row per region.
    """

    temperature_ranges: tuple[float, ...]
    coefficients: tuple[tuple[float, ...], ...]

    def __post_init__(self) -> None:
        ranges, rows = _normalise_polynomial_data(
            self.temperature_ranges,
            self.coefficients,
            coefficient_count=7,
        )
        object.__setattr__(self, "temperature_ranges", ranges)
        object.__setattr__(self, "coefficients", rows)

    @property
    def temperature_range(self) -> tuple[float, float]:
        """Inclusive fitted temperature range in kelvin."""
        return self.temperature_ranges[0], self.temperature_ranges[-1]

    def _dimensionless_properties(
        self,
        temperature: ArrayLike,
        *,
        allow_extrapolation: bool,
        warn: bool,
    ) -> tuple[FloatArray, bool, FloatArray, FloatArray, FloatArray]:
        values, scalar = _validated_temperature(
            temperature,
            self.temperature_range,
            allow_extrapolation=allow_extrapolation,
            warn=warn,
        )
        a = _region_coefficients(values, self.temperature_ranges, self.coefficients)
        cp_over_r = (
            a[..., 0]
            + a[..., 1] * values
            + a[..., 2] * values**2
            + a[..., 3] * values**3
            + a[..., 4] * values**4
        )
        h_over_rt = (
            a[..., 0]
            + a[..., 1] * values / 2.0
            + a[..., 2] * values**2 / 3.0
            + a[..., 3] * values**3 / 4.0
            + a[..., 4] * values**4 / 5.0
            + a[..., 5] / values
        )
        s_over_r = (
            a[..., 0] * np.log(values)
            + a[..., 1] * values
            + a[..., 2] * values**2 / 2.0
            + a[..., 3] * values**3 / 3.0
            + a[..., 4] * values**4 / 4.0
            + a[..., 6]
        )
        return values, scalar, cp_over_r, h_over_rt, s_over_r

    def cp_over_r(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return dimensionless standard-state heat capacity, ``cp/R``."""
        _, scalar, cp_over_r, _, _ = self._dimensionless_properties(
            temperature,
            allow_extrapolation=allow_extrapolation,
            warn=True,
        )
        return return_float(cp_over_r, scalar=scalar)

    def h_over_rt(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return dimensionless standard-state enthalpy, ``h/(R T)``."""
        _, scalar, _, h_over_rt, _ = self._dimensionless_properties(
            temperature,
            allow_extrapolation=allow_extrapolation,
            warn=True,
        )
        return return_float(h_over_rt, scalar=scalar)

    def s_over_r(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return dimensionless standard-state entropy, ``s/R``."""
        _, scalar, _, _, s_over_r = self._dimensionless_properties(
            temperature,
            allow_extrapolation=allow_extrapolation,
            warn=True,
        )
        return return_float(s_over_r, scalar=scalar)


@dataclass(frozen=True, slots=True)
class NASA9Polynomial:
    """Piecewise NASA nine-coefficient ideal-gas polynomial.

    ``temperature_ranges`` contains the region boundaries in kelvin.
    ``coefficients`` contains one ``(a1, ..., a9)`` row per region.
    """

    temperature_ranges: tuple[float, ...]
    coefficients: tuple[tuple[float, ...], ...]

    def __post_init__(self) -> None:
        ranges, rows = _normalise_polynomial_data(
            self.temperature_ranges,
            self.coefficients,
            coefficient_count=9,
        )
        object.__setattr__(self, "temperature_ranges", ranges)
        object.__setattr__(self, "coefficients", rows)

    @property
    def temperature_range(self) -> tuple[float, float]:
        """Inclusive fitted temperature range in kelvin."""
        return self.temperature_ranges[0], self.temperature_ranges[-1]

    def _dimensionless_properties(
        self,
        temperature: ArrayLike,
        *,
        allow_extrapolation: bool,
        warn: bool,
    ) -> tuple[FloatArray, bool, FloatArray, FloatArray, FloatArray]:
        values, scalar = _validated_temperature(
            temperature,
            self.temperature_range,
            allow_extrapolation=allow_extrapolation,
            warn=warn,
        )
        a = _region_coefficients(values, self.temperature_ranges, self.coefficients)
        cp_over_r = (
            a[..., 0] / values**2
            + a[..., 1] / values
            + a[..., 2]
            + a[..., 3] * values
            + a[..., 4] * values**2
            + a[..., 5] * values**3
            + a[..., 6] * values**4
        )
        h_over_rt = (
            -a[..., 0] / values**2
            + a[..., 1] * np.log(values) / values
            + a[..., 2]
            + a[..., 3] * values / 2.0
            + a[..., 4] * values**2 / 3.0
            + a[..., 5] * values**3 / 4.0
            + a[..., 6] * values**4 / 5.0
            + a[..., 7] / values
        )
        s_over_r = (
            -a[..., 0] / (2.0 * values**2)
            - a[..., 1] / values
            + a[..., 2] * np.log(values)
            + a[..., 3] * values
            + a[..., 4] * values**2 / 2.0
            + a[..., 5] * values**3 / 3.0
            + a[..., 6] * values**4 / 4.0
            + a[..., 8]
        )
        return values, scalar, cp_over_r, h_over_rt, s_over_r

    def cp_over_r(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return dimensionless standard-state heat capacity, ``cp/R``."""
        _, scalar, cp_over_r, _, _ = self._dimensionless_properties(
            temperature,
            allow_extrapolation=allow_extrapolation,
            warn=True,
        )
        return return_float(cp_over_r, scalar=scalar)

    def h_over_rt(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return dimensionless standard-state enthalpy, ``h/(R T)``."""
        _, scalar, _, h_over_rt, _ = self._dimensionless_properties(
            temperature,
            allow_extrapolation=allow_extrapolation,
            warn=True,
        )
        return return_float(h_over_rt, scalar=scalar)

    def s_over_r(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return dimensionless standard-state entropy, ``s/R``."""
        _, scalar, _, _, s_over_r = self._dimensionless_properties(
            temperature,
            allow_extrapolation=allow_extrapolation,
            warn=True,
        )
        return return_float(s_over_r, scalar=scalar)


@dataclass(frozen=True, slots=True)
class IdealGasSpecies:
    """Ideal-gas species backed by a NASA polynomial."""

    name: str
    molar_mass: float
    thermo: NASA7Polynomial | NASA9Polynomial
    reference_pressure: float = STANDARD_PRESSURE

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("name must not be empty")
        if not np.isfinite(self.molar_mass) or self.molar_mass <= 0.0:
            raise ValueError("molar_mass must be finite and greater than zero")
        if not isinstance(self.thermo, (NASA7Polynomial, NASA9Polynomial)):
            raise TypeError("thermo must be a NASA7Polynomial or NASA9Polynomial")
        if not np.isfinite(self.reference_pressure) or self.reference_pressure <= 0.0:
            raise ValueError("reference_pressure must be finite and greater than zero")

    @property
    def temperature_range(self) -> tuple[float, float]:
        """Inclusive fitted temperature range in kelvin."""
        return self.thermo.temperature_range

    @property
    def specific_gas_constant(self) -> float:
        """Species-specific gas constant in J/(kg K)."""
        return UNIVERSAL_GAS_CONSTANT / self.molar_mass

    def _standard_molar_properties(
        self,
        temperature: ArrayLike,
        *,
        allow_extrapolation: bool,
        warn: bool,
    ) -> tuple[FloatArray, bool, FloatArray, FloatArray, FloatArray]:
        values, scalar, cp_over_r, h_over_rt, s_over_r = (
            self.thermo._dimensionless_properties(
                temperature,
                allow_extrapolation=allow_extrapolation,
                warn=warn,
            )
        )
        cp = UNIVERSAL_GAS_CONSTANT * cp_over_r
        enthalpy = UNIVERSAL_GAS_CONSTANT * values * h_over_rt
        entropy = UNIVERSAL_GAS_CONSTANT * s_over_r
        return values, scalar, cp, enthalpy, entropy

    def standard_molar_cp(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return standard-state molar ``cp`` in J/(mol K)."""
        _, scalar, cp, _, _ = self._standard_molar_properties(
            temperature,
            allow_extrapolation=allow_extrapolation,
            warn=True,
        )
        return return_float(cp, scalar=scalar)

    def standard_molar_enthalpy(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return standard-state molar enthalpy in J/mol."""
        _, scalar, _, enthalpy, _ = self._standard_molar_properties(
            temperature,
            allow_extrapolation=allow_extrapolation,
            warn=True,
        )
        return return_float(enthalpy, scalar=scalar)

    def standard_molar_entropy(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return standard-state molar entropy in J/(mol K)."""
        _, scalar, _, _, entropy = self._standard_molar_properties(
            temperature,
            allow_extrapolation=allow_extrapolation,
            warn=True,
        )
        return return_float(entropy, scalar=scalar)


@dataclass(frozen=True, slots=True)
class ThermallyPerfectGas:
    """Frozen-composition thermally perfect ideal-gas mixture."""

    species: tuple[IdealGasSpecies, ...]
    mole_fractions: tuple[float, ...]

    def __post_init__(self) -> None:
        if not self.species:
            raise ValueError("species must contain at least one component")
        if len(self.species) != len(self.mole_fractions):
            raise ValueError("species and mole_fractions must have the same length")
        if any(not isinstance(item, IdealGasSpecies) for item in self.species):
            raise TypeError("species must contain only IdealGasSpecies objects")
        if len({item.name for item in self.species}) != len(self.species):
            raise ValueError("species names must be unique")

        try:
            fractions = tuple(float(value) for value in self.mole_fractions)
        except (TypeError, ValueError) as error:
            raise ValueError(
                "mole_fractions must contain real numeric values"
            ) from error
        if not np.all(np.isfinite(fractions)) or any(
            value <= 0.0 for value in fractions
        ):
            raise ValueError("mole_fractions must be finite and greater than zero")
        total = sum(fractions)
        if not np.isclose(total, 1.0, rtol=0.0, atol=1.0e-12):
            raise ValueError("mole_fractions must sum to one")
        fractions = tuple(value / total for value in fractions)
        object.__setattr__(self, "mole_fractions", fractions)

        minimum, maximum = self.temperature_range
        if minimum > maximum:
            raise ValueError("species do not share a common temperature range")

    @property
    def temperature_range(self) -> tuple[float, float]:
        """Inclusive common fitted temperature range in kelvin."""
        return (
            max(item.temperature_range[0] for item in self.species),
            min(item.temperature_range[1] for item in self.species),
        )

    @property
    def molar_mass(self) -> float:
        """Mixture molar mass in kg/mol."""
        return sum(
            fraction * item.molar_mass
            for item, fraction in zip(self.species, self.mole_fractions, strict=True)
        )

    @property
    def mass_fractions(self) -> tuple[float, ...]:
        """Mixture mass fractions in component order."""
        molar_mass = self.molar_mass
        return tuple(
            fraction * item.molar_mass / molar_mass
            for item, fraction in zip(self.species, self.mole_fractions, strict=True)
        )

    @property
    def specific_gas_constant(self) -> float:
        """Mixture-specific gas constant in J/(kg K)."""
        return UNIVERSAL_GAS_CONSTANT / self.molar_mass

    def _validated_temperature(
        self,
        temperature: ArrayLike,
        *,
        allow_extrapolation: bool,
        warn: bool = True,
    ) -> tuple[FloatArray, bool]:
        return _validated_temperature(
            temperature,
            self.temperature_range,
            allow_extrapolation=allow_extrapolation,
            warn=warn,
        )

    def _standard_molar_properties_from_values(
        self, temperature: FloatArray
    ) -> tuple[FloatArray, FloatArray, FloatArray]:
        cp = np.zeros_like(temperature, dtype=np.float64)
        enthalpy = np.zeros_like(temperature, dtype=np.float64)
        entropy = np.zeros_like(temperature, dtype=np.float64)
        for item, fraction in zip(self.species, self.mole_fractions, strict=True):
            _, _, cp_i, enthalpy_i, entropy_i = item._standard_molar_properties(
                temperature,
                allow_extrapolation=True,
                warn=False,
            )
            cp += fraction * cp_i
            enthalpy += fraction * enthalpy_i
            entropy += fraction * entropy_i
        return cp, enthalpy, entropy

    def _cp_array(
        self, temperature: ArrayLike, *, allow_extrapolation: bool
    ) -> tuple[FloatArray, bool, FloatArray]:
        values, scalar = self._validated_temperature(
            temperature,
            allow_extrapolation=allow_extrapolation,
        )
        molar_cp, _, _ = self._standard_molar_properties_from_values(values)
        return values, scalar, molar_cp / self.molar_mass

    def cp(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return mixture ``cp`` in J/(kg K)."""
        _, scalar, cp = self._cp_array(
            temperature, allow_extrapolation=allow_extrapolation
        )
        return return_float(cp, scalar=scalar)

    def cv(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return mixture ``cv`` in J/(kg K)."""
        _, scalar, cp = self._cp_array(
            temperature, allow_extrapolation=allow_extrapolation
        )
        cv = cp - self.specific_gas_constant
        return return_float(cv, scalar=scalar)

    def heat_capacity_ratio(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return temperature-dependent heat-capacity ratio, ``cp/cv``."""
        _, scalar, cp = self._cp_array(
            temperature, allow_extrapolation=allow_extrapolation
        )
        result = cp / (cp - self.specific_gas_constant)
        return return_float(result, scalar=scalar)

    def speed_of_sound(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return frozen-composition ideal-gas speed of sound in m/s."""
        values, scalar, cp = self._cp_array(
            temperature, allow_extrapolation=allow_extrapolation
        )
        gamma = cp / (cp - self.specific_gas_constant)
        result = np.sqrt(gamma * self.specific_gas_constant * values)
        return return_float(result, scalar=scalar)

    def standard_enthalpy(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return standard-state specific enthalpy in J/kg."""
        values, scalar = self._validated_temperature(
            temperature,
            allow_extrapolation=allow_extrapolation,
        )
        _, molar_enthalpy, _ = self._standard_molar_properties_from_values(values)
        return return_float(molar_enthalpy / self.molar_mass, scalar=scalar)

    def standard_internal_energy(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return standard-state specific internal energy in J/kg."""
        values, scalar = self._validated_temperature(
            temperature,
            allow_extrapolation=allow_extrapolation,
        )
        _, molar_enthalpy, _ = self._standard_molar_properties_from_values(values)
        result = molar_enthalpy / self.molar_mass - self.specific_gas_constant * values
        return return_float(result, scalar=scalar)

    def _sensible_property(
        self,
        temperature: ArrayLike,
        reference_temperature: float,
        *,
        internal_energy: bool,
        allow_extrapolation: bool,
    ) -> FloatResult:
        values, scalar = as_float_array(temperature, name="temperature")
        reference = np.asarray(reference_temperature, dtype=np.float64)
        if reference.ndim != 0 or not np.isfinite(reference):
            raise ValueError("reference_temperature must be a finite scalar")
        reference_value = float(reference)
        combined = np.concatenate((values.reshape(-1), [reference_value]))
        _validated_temperature(
            combined,
            self.temperature_range,
            allow_extrapolation=allow_extrapolation,
            warn=True,
        )

        _, enthalpy, _ = self._standard_molar_properties_from_values(values)
        reference_array = np.asarray(reference_value, dtype=np.float64)
        _, reference_enthalpy, _ = self._standard_molar_properties_from_values(
            reference_array
        )
        result = (enthalpy - reference_enthalpy) / self.molar_mass
        if internal_energy:
            result -= self.specific_gas_constant * (values - reference_value)
        return return_float(result, scalar=scalar)

    def sensible_enthalpy(
        self,
        temperature: ArrayLike,
        reference_temperature: float = 298.15,
        *,
        allow_extrapolation: bool = False,
    ) -> FloatResult:
        """Return specific enthalpy relative to a reference temperature."""
        return self._sensible_property(
            temperature,
            reference_temperature,
            internal_energy=False,
            allow_extrapolation=allow_extrapolation,
        )

    def sensible_internal_energy(
        self,
        temperature: ArrayLike,
        reference_temperature: float = 298.15,
        *,
        allow_extrapolation: bool = False,
    ) -> FloatResult:
        """Return specific internal energy relative to a reference temperature."""
        return self._sensible_property(
            temperature,
            reference_temperature,
            internal_energy=True,
            allow_extrapolation=allow_extrapolation,
        )

    def entropy(
        self,
        temperature: ArrayLike,
        pressure: ArrayLike,
        *,
        allow_extrapolation: bool = False,
    ) -> FloatResult:
        """Return ideal-mixture specific entropy in J/(kg K).

        Species standard-state entropy, ideal mixing, and pressure corrections
        are included using each species' declared reference pressure.
        """
        temperature_values, temperature_scalar = self._validated_temperature(
            temperature,
            allow_extrapolation=allow_extrapolation,
        )
        pressure_values, pressure_scalar = as_float_array(pressure, name="pressure")
        if np.any(pressure_values <= 0.0):
            raise ValueError("pressure must be greater than zero")
        try:
            temperatures, pressures = np.broadcast_arrays(
                temperature_values, pressure_values
            )
        except ValueError as error:
            raise ValueError(
                "temperature and pressure must be broadcast-compatible"
            ) from error

        _, _, standard_entropy = self._standard_molar_properties_from_values(
            cast(FloatArray, temperatures)
        )
        pressure_correction = np.zeros_like(temperatures, dtype=np.float64)
        for item, fraction in zip(self.species, self.mole_fractions, strict=True):
            pressure_correction += fraction * np.log(
                fraction * pressures / item.reference_pressure
            )
        result = (
            standard_entropy - UNIVERSAL_GAS_CONSTANT * pressure_correction
        ) / self.molar_mass
        return return_float(
            result,
            scalar=temperature_scalar and pressure_scalar,
        )


def _build_air(data: dict[str, PolynomialData], *, nasa9: bool) -> ThermallyPerfectGas:
    names = ("N2", "O2", "Ar", "CO2")
    raw_fractions = tuple(DRY_AIR_MOLE_FRACTIONS[name] for name in names)
    fraction_total = sum(raw_fractions)
    fractions = tuple(value / fraction_total for value in raw_fractions)
    species: list[IdealGasSpecies] = []
    for name in names:
        molar_mass, ranges, rows = data[name]
        polynomial: NASA7Polynomial | NASA9Polynomial
        if nasa9:
            polynomial = NASA9Polynomial(ranges, rows)
        else:
            polynomial = NASA7Polynomial(ranges, rows)
        species.append(
            IdealGasSpecies(
                name=name,
                molar_mass=molar_mass,
                thermo=polynomial,
                reference_pressure=STANDARD_PRESSURE,
            )
        )
    return ThermallyPerfectGas(tuple(species), fractions)


AIR_NASA7 = _build_air(NASA7_DATA, nasa9=False)
"""Frozen-composition dry air using NASA seven-coefficient polynomials."""

AIR_NASA9 = _build_air(NASA9_DATA, nasa9=True)
"""Frozen-composition dry air using NASA nine-coefficient polynomials."""

__all__ = [
    "AIR_NASA7",
    "AIR_NASA9",
    "STANDARD_PRESSURE",
    "UNIVERSAL_GAS_CONSTANT",
    "IdealGasSpecies",
    "NASA7Polynomial",
    "NASA9Polynomial",
    "ThermallyPerfectGas",
]
