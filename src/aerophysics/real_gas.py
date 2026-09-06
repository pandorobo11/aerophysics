"""Harmonic-oscillator and Beattie--Bridgeman gas models.

The harmonic-oscillator relations follow Kennard's statistical-mechanics
model for vibrational excitation.  The dense-gas corrections follow the
Beattie--Bridgeman equation of state in the form used by Randall and by the
JAXA hypersonic-wind-tunnel implementation.  Chemical reactions, phase
changes, dissociation, and ionisation are not included.
"""

import warnings
from dataclasses import dataclass
from typing import Final, cast

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import brentq

from aerophysics._array import FloatArray, FloatResult, as_float_array, return_float
from aerophysics.exceptions import ApplicabilityWarning, ModelRangeError

_ROOT_XTOL: Final = 1e-12
_ROOT_RTOL: Final = 4.0 * np.finfo(np.float64).eps
_ROOT_MAXITER: Final = 100
_DENSITY_BRACKET_STEPS: Final = 200


def _require_finite_positive(value: float, *, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be finite and greater than zero")
    return result


def _normalise_range(
    value: tuple[float, float] | None, *, name: str
) -> tuple[float, float] | None:
    if value is None:
        return None
    if len(value) != 2:
        raise ValueError(f"{name} must contain two values")
    lower = _require_finite_positive(value[0], name=f"{name}[0]")
    upper = _require_finite_positive(value[1], name=f"{name}[1]")
    if upper <= lower:
        raise ValueError(f"{name} must be strictly increasing")
    return lower, upper


def _warn_or_raise_outside(
    values: FloatArray,
    applicable_range: tuple[float, float] | None,
    *,
    name: str,
    unit: str,
    allow_extrapolation: bool,
    warn: bool,
) -> bool:
    if applicable_range is None:
        return False
    lower, upper = applicable_range
    outside = bool(np.any((values < lower) | (values > upper)))
    if outside and not allow_extrapolation:
        raise ModelRangeError(f"{name} must be within {lower:g}--{upper:g} {unit}")
    if outside and warn:
        warnings.warn(
            f"{name} is outside the documented range {lower:g}--{upper:g} {unit}",
            ApplicabilityWarning,
            stacklevel=3,
        )
    return outside


def _broadcast_temperature_pressure(
    temperature: ArrayLike, pressure: ArrayLike
) -> tuple[FloatArray, FloatArray, bool]:
    temperatures, temperature_scalar = as_float_array(temperature, name="temperature")
    pressures, pressure_scalar = as_float_array(pressure, name="pressure")
    if np.any(temperatures <= 0.0):
        raise ValueError("temperature must be greater than zero")
    if np.any(pressures <= 0.0):
        raise ValueError("pressure must be greater than zero")
    try:
        broadcast_temperature, broadcast_pressure = np.broadcast_arrays(
            temperatures, pressures
        )
    except ValueError as error:
        raise ValueError(
            "temperature and pressure must be broadcast-compatible"
        ) from error
    return (
        np.asarray(broadcast_temperature, dtype=np.float64),
        np.asarray(broadcast_pressure, dtype=np.float64),
        temperature_scalar and pressure_scalar,
    )


def _broadcast_temperature_density(
    temperature: ArrayLike, density: ArrayLike
) -> tuple[FloatArray, FloatArray, bool]:
    temperatures, temperature_scalar = as_float_array(temperature, name="temperature")
    densities, density_scalar = as_float_array(density, name="density")
    if np.any(temperatures <= 0.0):
        raise ValueError("temperature must be greater than zero")
    if np.any(densities <= 0.0):
        raise ValueError("density must be greater than zero")
    try:
        broadcast_temperature, broadcast_density = np.broadcast_arrays(
            temperatures, densities
        )
    except ValueError as error:
        raise ValueError(
            "temperature and density must be broadcast-compatible"
        ) from error
    return (
        np.asarray(broadcast_temperature, dtype=np.float64),
        np.asarray(broadcast_density, dtype=np.float64),
        temperature_scalar and density_scalar,
    )


@dataclass(frozen=True, slots=True)
class VibrationalMode:
    """One harmonic vibrational mode.

    ``weight`` is the mode contribution per mixture mole and
    ``characteristic_temperature`` is the vibrational temperature in kelvin.
    """

    weight: float
    characteristic_temperature: float

    def __post_init__(self) -> None:
        weight = float(self.weight)
        if not np.isfinite(weight) or weight < 0.0:
            raise ValueError("weight must be finite and non-negative")
        object.__setattr__(self, "weight", weight)
        object.__setattr__(
            self,
            "characteristic_temperature",
            _require_finite_positive(
                self.characteristic_temperature,
                name="characteristic_temperature",
            ),
        )


@dataclass(frozen=True, slots=True)
class ThermodynamicState:
    """Thermodynamic properties at one or more temperature-pressure states."""

    temperature: FloatResult
    pressure: FloatResult
    density: FloatResult
    internal_energy: FloatResult
    enthalpy: FloatResult
    entropy: FloatResult
    cp: FloatResult
    cv: FloatResult
    heat_capacity_ratio: FloatResult
    speed_of_sound: FloatResult


@dataclass(frozen=True, slots=True)
class _ScalarThermodynamicState:
    temperature: float
    pressure: float
    density: float
    internal_energy: float
    enthalpy: float
    entropy: float
    cp: float
    cv: float
    heat_capacity_ratio: float
    speed_of_sound: float


def _vibrational_terms(
    temperature: FloatArray, modes: tuple[VibrationalMode, ...]
) -> tuple[FloatArray, FloatArray, FloatArray]:
    heat_capacity = np.zeros_like(temperature)
    energy = np.zeros_like(temperature)
    entropy = np.zeros_like(temperature)
    for mode in modes:
        x = mode.characteristic_temperature / temperature
        exp_minus_x = np.exp(-x)
        denominator = -np.expm1(-x)
        occupation = exp_minus_x / denominator
        energy += mode.weight * x * occupation
        heat_capacity += mode.weight * x**2 * exp_minus_x / denominator**2
        entropy += mode.weight * (x * occupation - np.log1p(-exp_minus_x))
    return heat_capacity, energy, entropy


@dataclass(frozen=True, slots=True)
class HarmonicOscillatorGas:
    """Thermally perfect ideal gas with harmonic vibrational modes."""

    specific_gas_constant: float
    base_heat_capacity_ratio: float
    modes: tuple[VibrationalMode, ...] = ()
    applicable_temperature_range: tuple[float, float] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "specific_gas_constant",
            _require_finite_positive(
                self.specific_gas_constant, name="specific_gas_constant"
            ),
        )
        gamma = float(self.base_heat_capacity_ratio)
        if not np.isfinite(gamma) or gamma <= 1.0:
            raise ValueError(
                "base_heat_capacity_ratio must be finite and greater than one"
            )
        object.__setattr__(self, "base_heat_capacity_ratio", gamma)
        object.__setattr__(self, "modes", tuple(self.modes))
        if any(not isinstance(mode, VibrationalMode) for mode in self.modes):
            raise TypeError("modes must contain VibrationalMode values")
        object.__setattr__(
            self,
            "applicable_temperature_range",
            _normalise_range(
                self.applicable_temperature_range,
                name="applicable_temperature_range",
            ),
        )

    @property
    def temperature_range(self) -> tuple[float, float]:
        """Documented temperature range, or the positive floating-point range."""
        if self.applicable_temperature_range is None:
            return np.finfo(np.float64).tiny, np.finfo(np.float64).max
        return self.applicable_temperature_range

    def _temperature(
        self,
        temperature: ArrayLike,
        *,
        allow_extrapolation: bool,
        warn: bool,
    ) -> tuple[FloatArray, bool, bool]:
        values, scalar = as_float_array(temperature, name="temperature")
        if np.any(values <= 0.0):
            raise ValueError("temperature must be greater than zero")
        outside = _warn_or_raise_outside(
            values,
            self.applicable_temperature_range,
            name="temperature",
            unit="K",
            allow_extrapolation=allow_extrapolation,
            warn=warn,
        )
        return values, scalar, outside

    def _properties_from_values(
        self, temperature: FloatArray
    ) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]:
        vibration_cv, vibration_energy, vibration_entropy = _vibrational_terms(
            temperature, self.modes
        )
        gas_constant = self.specific_gas_constant
        base_cv = gas_constant / (self.base_heat_capacity_ratio - 1.0)
        cv = base_cv + gas_constant * vibration_cv
        cp = cv + gas_constant
        internal_energy = base_cv * temperature + (
            gas_constant * temperature * vibration_energy
        )
        enthalpy = internal_energy + gas_constant * temperature
        entropy_temperature = (base_cv + gas_constant) * np.log(
            temperature
        ) + gas_constant * vibration_entropy
        return cp, cv, internal_energy, enthalpy, entropy_temperature

    def cp(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return constant-pressure specific heat in J/(kg K)."""
        values, scalar, _ = self._temperature(
            temperature, allow_extrapolation=allow_extrapolation, warn=True
        )
        cp, _, _, _, _ = self._properties_from_values(values)
        return return_float(cp, scalar=scalar)

    def cv(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return constant-volume specific heat in J/(kg K)."""
        values, scalar, _ = self._temperature(
            temperature, allow_extrapolation=allow_extrapolation, warn=True
        )
        _, cv, _, _, _ = self._properties_from_values(values)
        return return_float(cv, scalar=scalar)

    def heat_capacity_ratio(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return temperature-dependent ``cp/cv``."""
        values, scalar, _ = self._temperature(
            temperature, allow_extrapolation=allow_extrapolation, warn=True
        )
        cp, cv, _, _, _ = self._properties_from_values(values)
        return return_float(cp / cv, scalar=scalar)

    def standard_enthalpy(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return ideal-gas specific enthalpy in J/kg."""
        values, scalar, _ = self._temperature(
            temperature, allow_extrapolation=allow_extrapolation, warn=True
        )
        _, _, _, enthalpy, _ = self._properties_from_values(values)
        return return_float(enthalpy, scalar=scalar)

    def standard_internal_energy(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return ideal-gas specific internal energy in J/kg."""
        values, scalar, _ = self._temperature(
            temperature, allow_extrapolation=allow_extrapolation, warn=True
        )
        _, _, internal_energy, _, _ = self._properties_from_values(values)
        return return_float(internal_energy, scalar=scalar)

    def entropy(
        self,
        temperature: ArrayLike,
        pressure: ArrayLike,
        *,
        allow_extrapolation: bool = False,
    ) -> FloatResult:
        """Return specific entropy up to an arbitrary additive constant."""
        temperatures, pressures, scalar = _broadcast_temperature_pressure(
            temperature, pressure
        )
        _warn_or_raise_outside(
            temperatures,
            self.applicable_temperature_range,
            name="temperature",
            unit="K",
            allow_extrapolation=allow_extrapolation,
            warn=True,
        )
        _, _, _, _, entropy_temperature = self._properties_from_values(temperatures)
        result = entropy_temperature - self.specific_gas_constant * np.log(pressures)
        return return_float(result, scalar=scalar)

    def speed_of_sound(
        self, temperature: ArrayLike, *, allow_extrapolation: bool = False
    ) -> FloatResult:
        """Return ideal-gas frozen-composition speed of sound in m/s."""
        values, scalar, _ = self._temperature(
            temperature, allow_extrapolation=allow_extrapolation, warn=True
        )
        cp, cv, _, _, _ = self._properties_from_values(values)
        result = np.sqrt(cp / cv * self.specific_gas_constant * values)
        return return_float(result, scalar=scalar)

    def state(
        self,
        temperature: ArrayLike,
        pressure: ArrayLike,
        *,
        allow_extrapolation: bool = False,
    ) -> ThermodynamicState:
        """Return the complete ideal-gas thermodynamic state."""
        temperatures, pressures, scalar = _broadcast_temperature_pressure(
            temperature, pressure
        )
        _warn_or_raise_outside(
            temperatures,
            self.applicable_temperature_range,
            name="temperature",
            unit="K",
            allow_extrapolation=allow_extrapolation,
            warn=True,
        )
        cp, cv, internal_energy, enthalpy, entropy_temperature = (
            self._properties_from_values(temperatures)
        )
        density = pressures / (self.specific_gas_constant * temperatures)
        gamma = cp / cv
        sound_speed = np.sqrt(gamma * self.specific_gas_constant * temperatures)
        entropy = entropy_temperature - self.specific_gas_constant * np.log(pressures)

        def output(value: FloatArray) -> FloatResult:
            return return_float(value, scalar=scalar)

        return ThermodynamicState(
            temperature=output(temperatures),
            pressure=output(pressures),
            density=output(density),
            internal_energy=output(internal_energy),
            enthalpy=output(enthalpy),
            entropy=output(entropy),
            cp=output(cp),
            cv=output(cv),
            heat_capacity_ratio=output(gamma),
            speed_of_sound=output(sound_speed),
        )


@dataclass(frozen=True, slots=True)
class BeattieBridgemanGas:
    """Beattie--Bridgeman gas with harmonic vibrational heat capacity.

    An ``applicable_*_range`` records the conditions tabulated for a
    particular use of the model; it is not an inherent physical-validity
    limit of the equation of state.  In particular, the air preset uses the
    R. E. Randall, *AEDC-TR-57-8* (1957) tabulated range.
    """

    specific_gas_constant: float
    base_heat_capacity_ratio: float
    a0: float
    b0: float
    a: float
    b: float
    c: float
    modes: tuple[VibrationalMode, ...] = ()
    applicable_temperature_range: tuple[float, float] | None = None
    applicable_pressure_range: tuple[float, float] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "specific_gas_constant",
            _require_finite_positive(
                self.specific_gas_constant, name="specific_gas_constant"
            ),
        )
        gamma = float(self.base_heat_capacity_ratio)
        if not np.isfinite(gamma) or gamma <= 1.0:
            raise ValueError(
                "base_heat_capacity_ratio must be finite and greater than one"
            )
        object.__setattr__(self, "base_heat_capacity_ratio", gamma)
        for name in ("a0", "b0", "a", "b", "c"):
            value = float(getattr(self, name))
            if not np.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)
        object.__setattr__(self, "modes", tuple(self.modes))
        if any(not isinstance(mode, VibrationalMode) for mode in self.modes):
            raise TypeError("modes must contain VibrationalMode values")
        object.__setattr__(
            self,
            "applicable_temperature_range",
            _normalise_range(
                self.applicable_temperature_range,
                name="applicable_temperature_range",
            ),
        )
        object.__setattr__(
            self,
            "applicable_pressure_range",
            _normalise_range(
                self.applicable_pressure_range,
                name="applicable_pressure_range",
            ),
        )

    @property
    def temperature_range(self) -> tuple[float, float]:
        """Documented reservoir-temperature range."""
        if self.applicable_temperature_range is None:
            return np.finfo(np.float64).tiny, np.finfo(np.float64).max
        return self.applicable_temperature_range

    def _coefficients(self, temperature: float) -> tuple[float, float, float]:
        gas_constant = self.specific_gas_constant
        e1 = self.b0 - self.a0 / (gas_constant * temperature) - self.c / temperature**3
        e2 = (
            self.a0 * self.a / (gas_constant * temperature)
            - self.b0 * self.b
            - self.b0 * self.c / temperature**3
        )
        e3 = self.b0 * self.b * self.c / temperature**3
        return e1, e2, e3

    def _pressure_scalar(self, temperature: float, density: float) -> float:
        e1, e2, e3 = self._coefficients(temperature)
        return float(
            density
            * self.specific_gas_constant
            * temperature
            * (1.0 + e1 * density + e2 * density**2 + e3 * density**3)
        )

    def _dp_drho_scalar(self, temperature: float, density: float) -> float:
        e1, e2, e3 = self._coefficients(temperature)
        return float(
            self.specific_gas_constant
            * temperature
            * (1.0 + 2.0 * e1 * density + 3.0 * e2 * density**2 + 4.0 * e3 * density**3)
        )

    def _first_spinodal_density(self, temperature: float) -> float | None:
        """Return the lowest positive root of ``(dp/drho)_T``, if present."""
        e1, e2, e3 = self._coefficients(temperature)
        roots = np.roots((4.0 * e3, 3.0 * e2, 2.0 * e1, 1.0))
        positive_roots = [
            float(root.real)
            for root in roots
            if np.isfinite(root.real)
            and np.isfinite(root.imag)
            and root.real > 0.0
            and abs(root.imag) <= 1e-10 * max(1.0, abs(root.real))
        ]
        return min(positive_roots, default=None)

    def _dp_dtemperature_scalar(self, temperature: float, density: float) -> float:
        return float(
            density
            * self.specific_gas_constant
            * (1.0 + 2.0 * self.c * density / temperature**3)
            * (1.0 + self.b0 * density - self.b0 * self.b * density**2)
        )

    def _density_scalar(self, temperature: float, pressure: float) -> float:
        spinodal = self._first_spinodal_density(temperature)
        stable_upper: float | None = None
        if spinodal is not None:
            # The public density inverse is defined on the gas branch connected
            # continuously to rho -> 0.  Restrict the bracket to just below the
            # first isothermal spinodal so a narrow unstable interval cannot be
            # skipped in favour of a separate high-density stable root.
            stable_upper = np.nextafter(spinodal * (1.0 - 1.0e-10), 0.0)
            for _ in range(64):
                if self._dp_drho_scalar(temperature, stable_upper) > 0.0:
                    break
                stable_upper = np.nextafter(stable_upper, 0.0)
            else:
                raise ModelRangeError(
                    "could not establish the mechanically stable "
                    "Beattie--Bridgeman gas branch"
                )

        ideal_density = pressure / (self.specific_gas_constant * temperature)
        lower = 0.0
        upper = max(0.25 * ideal_density, np.finfo(np.float64).tiny)
        if stable_upper is not None:
            upper = min(upper, stable_upper)
        for _ in range(_DENSITY_BRACKET_STEPS):
            residual = self._pressure_scalar(temperature, upper) - pressure
            derivative = self._dp_drho_scalar(temperature, upper)
            if residual >= 0.0:
                density = float(
                    brentq(
                        lambda value: (
                            self._pressure_scalar(temperature, value) - pressure
                        ),
                        lower,
                        upper,
                        xtol=_ROOT_XTOL,
                        rtol=_ROOT_RTOL,
                        maxiter=_ROOT_MAXITER,
                    )
                )
                if self._dp_drho_scalar(temperature, density) <= 0.0:
                    break
                return density
            if derivative <= 0.0:
                break
            lower = upper
            next_upper = upper * 1.5
            if stable_upper is not None:
                next_upper = min(next_upper, stable_upper)
            if next_upper <= upper:
                break
            upper = next_upper
        raise ModelRangeError(
            "no mechanically stable Beattie--Bridgeman gas-phase density root"
        )

    def _scalar_state_from_density(
        self, temperature: float, density: float
    ) -> _ScalarThermodynamicState:
        pressure = self._pressure_scalar(temperature, density)
        if not np.isfinite(pressure) or pressure <= 0.0:
            raise ModelRangeError("Beattie--Bridgeman pressure is non-physical")
        dp_drho = self._dp_drho_scalar(temperature, density)
        dp_dtemperature = self._dp_dtemperature_scalar(temperature, density)
        if not np.isfinite(dp_drho) or dp_drho <= 0.0:
            raise ModelRangeError(
                "Beattie--Bridgeman isothermal compressibility is non-physical"
            )

        temperature_array = np.asarray(temperature, dtype=np.float64)
        vibration_cv, vibration_energy, vibration_entropy = _vibrational_terms(
            temperature_array, self.modes
        )
        gas_constant = self.specific_gas_constant
        base_cv = gas_constant / (self.base_heat_capacity_ratio - 1.0)
        cv = float(
            base_cv
            + gas_constant * vibration_cv
            + 6.0
            * gas_constant
            * self.c
            / temperature**3
            * density
            * (1.0 + 0.5 * self.b0 * density - self.b0 * self.b * density**2 / 3.0)
        )
        cp = cv + temperature * dp_dtemperature**2 / (density**2 * dp_drho)
        gamma = cp / cv
        sound_speed_squared = gamma * dp_drho
        if (
            not np.isfinite(cp)
            or not np.isfinite(cv)
            or cp <= 0.0
            or cv <= 0.0
            or not np.isfinite(gamma)
            or gamma <= 0.0
            or not np.isfinite(sound_speed_squared)
            or sound_speed_squared <= 0.0
        ):
            raise ModelRangeError(
                "Beattie--Bridgeman heat capacity or sound speed is non-physical"
            )

        correction = (
            -density
            * (self.a0 / (gas_constant * temperature) + 3.0 * self.c / temperature**3)
            - density**2
            * (
                3.0 * self.b0 * self.c / (2.0 * temperature**3)
                - self.a0 * self.a / (2.0 * gas_constant * temperature)
            )
            + density**3 * self.b0 * self.b * self.c / temperature**3
        )
        internal_energy = float(
            base_cv * temperature
            + gas_constant * temperature * vibration_energy
            + gas_constant * temperature * correction
        )
        enthalpy = internal_energy + pressure / density
        log_f = (
            density * (self.b0 + 2.0 * self.c / temperature**3)
            + density**2 * (self.b0 * self.c / temperature**3 - 0.5 * self.b0 * self.b)
            - density**3 * 2.0 * self.b0 * self.b * self.c / (3.0 * temperature**3)
        )
        entropy = float(
            gas_constant
            * (
                np.log(temperature) / (self.base_heat_capacity_ratio - 1.0)
                - np.log(density)
                + vibration_entropy
                - log_f
            )
        )
        return _ScalarThermodynamicState(
            temperature=temperature,
            pressure=pressure,
            density=density,
            internal_energy=internal_energy,
            enthalpy=enthalpy,
            entropy=entropy,
            cp=cp,
            cv=cv,
            heat_capacity_ratio=gamma,
            speed_of_sound=float(np.sqrt(sound_speed_squared)),
        )

    def _scalar_state(
        self, temperature: float, pressure: float
    ) -> _ScalarThermodynamicState:
        density = self._density_scalar(temperature, pressure)
        return self._scalar_state_from_density(temperature, density)

    def _check_applicability(
        self,
        temperature: FloatArray,
        pressure: FloatArray,
        *,
        allow_extrapolation: bool,
        warn: bool,
    ) -> bool:
        temperature_outside = _warn_or_raise_outside(
            temperature,
            self.applicable_temperature_range,
            name="temperature",
            unit="K",
            allow_extrapolation=allow_extrapolation,
            warn=False,
        )
        pressure_outside = _warn_or_raise_outside(
            pressure,
            self.applicable_pressure_range,
            name="pressure",
            unit="Pa",
            allow_extrapolation=allow_extrapolation,
            warn=False,
        )
        outside = temperature_outside or pressure_outside
        if outside and warn:
            warnings.warn(
                f"state is outside the {self._applicability_description()}",
                ApplicabilityWarning,
                stacklevel=3,
            )
        return outside

    def _applicability_description(self) -> str:
        """Describe this model's documented table without misattribution."""
        if self is AIR_BEATTIE_BRIDGEMAN:
            return "Randall, AEDC-TR-57-8 Beattie--Bridgeman air tabulated range"
        return "configured Beattie--Bridgeman tabulated range"

    def pressure(
        self,
        temperature: ArrayLike,
        density: ArrayLike,
        *,
        allow_extrapolation: bool = False,
    ) -> FloatResult:
        """Return pressure in Pa from temperature and density."""
        temperatures, densities, scalar = _broadcast_temperature_density(
            temperature, density
        )
        _warn_or_raise_outside(
            temperatures,
            self.applicable_temperature_range,
            name="temperature",
            unit="K",
            allow_extrapolation=allow_extrapolation,
            warn=False,
        )
        result = np.empty_like(temperatures)
        for index in np.ndindex(result.shape):
            result[index] = self._pressure_scalar(
                float(temperatures[index]), float(densities[index])
            )
        if np.any(result <= 0.0) or not np.all(np.isfinite(result)):
            raise ModelRangeError("Beattie--Bridgeman pressure is non-physical")
        self._check_applicability(
            temperatures,
            result,
            allow_extrapolation=allow_extrapolation,
            warn=True,
        )
        return return_float(result, scalar=scalar)

    def density(
        self,
        temperature: ArrayLike,
        pressure: ArrayLike,
        *,
        allow_extrapolation: bool = False,
    ) -> FloatResult:
        """Return the mechanically stable gas-phase density in kg/m³."""
        temperatures, pressures, scalar = _broadcast_temperature_pressure(
            temperature, pressure
        )
        self._check_applicability(
            temperatures,
            pressures,
            allow_extrapolation=allow_extrapolation,
            warn=True,
        )
        result = np.empty_like(temperatures)
        for index in np.ndindex(result.shape):
            result[index] = self._density_scalar(
                float(temperatures[index]), float(pressures[index])
            )
        return return_float(result, scalar=scalar)

    def state(
        self,
        temperature: ArrayLike,
        pressure: ArrayLike,
        *,
        allow_extrapolation: bool = False,
    ) -> ThermodynamicState:
        """Return the complete Beattie--Bridgeman thermodynamic state."""
        temperatures, pressures, scalar = _broadcast_temperature_pressure(
            temperature, pressure
        )
        self._check_applicability(
            temperatures,
            pressures,
            allow_extrapolation=allow_extrapolation,
            warn=True,
        )
        fields = {
            name: np.empty_like(temperatures)
            for name in (
                "density",
                "internal_energy",
                "enthalpy",
                "entropy",
                "cp",
                "cv",
                "heat_capacity_ratio",
                "speed_of_sound",
            )
        }
        for index in np.ndindex(temperatures.shape):
            state = self._scalar_state(
                float(temperatures[index]), float(pressures[index])
            )
            for name, values in fields.items():
                values[index] = getattr(state, name)

        def output(value: FloatArray) -> FloatResult:
            return return_float(value, scalar=scalar)

        return ThermodynamicState(
            temperature=output(temperatures),
            pressure=output(pressures),
            density=output(fields["density"]),
            internal_energy=output(fields["internal_energy"]),
            enthalpy=output(fields["enthalpy"]),
            entropy=output(fields["entropy"]),
            cp=output(fields["cp"]),
            cv=output(fields["cv"]),
            heat_capacity_ratio=output(fields["heat_capacity_ratio"]),
            speed_of_sound=output(fields["speed_of_sound"]),
        )

    def _property(
        self,
        name: str,
        temperature: ArrayLike,
        pressure: ArrayLike,
        *,
        allow_extrapolation: bool,
    ) -> FloatResult:
        return cast(
            FloatResult,
            getattr(
                self.state(
                    temperature,
                    pressure,
                    allow_extrapolation=allow_extrapolation,
                ),
                name,
            ),
        )

    def cp(
        self,
        temperature: ArrayLike,
        pressure: ArrayLike,
        *,
        allow_extrapolation: bool = False,
    ) -> FloatResult:
        """Return constant-pressure specific heat in J/(kg K)."""
        return self._property(
            "cp", temperature, pressure, allow_extrapolation=allow_extrapolation
        )

    def cv(
        self,
        temperature: ArrayLike,
        pressure: ArrayLike,
        *,
        allow_extrapolation: bool = False,
    ) -> FloatResult:
        """Return constant-volume specific heat in J/(kg K)."""
        return self._property(
            "cv", temperature, pressure, allow_extrapolation=allow_extrapolation
        )

    def heat_capacity_ratio(
        self,
        temperature: ArrayLike,
        pressure: ArrayLike,
        *,
        allow_extrapolation: bool = False,
    ) -> FloatResult:
        """Return ``cp/cv``."""
        return self._property(
            "heat_capacity_ratio",
            temperature,
            pressure,
            allow_extrapolation=allow_extrapolation,
        )

    def internal_energy(
        self,
        temperature: ArrayLike,
        pressure: ArrayLike,
        *,
        allow_extrapolation: bool = False,
    ) -> FloatResult:
        """Return specific internal energy in J/kg."""
        return self._property(
            "internal_energy",
            temperature,
            pressure,
            allow_extrapolation=allow_extrapolation,
        )

    def enthalpy(
        self,
        temperature: ArrayLike,
        pressure: ArrayLike,
        *,
        allow_extrapolation: bool = False,
    ) -> FloatResult:
        """Return specific enthalpy in J/kg."""
        return self._property(
            "enthalpy",
            temperature,
            pressure,
            allow_extrapolation=allow_extrapolation,
        )

    def entropy(
        self,
        temperature: ArrayLike,
        pressure: ArrayLike,
        *,
        allow_extrapolation: bool = False,
    ) -> FloatResult:
        """Return specific entropy up to an arbitrary additive constant."""
        return self._property(
            "entropy",
            temperature,
            pressure,
            allow_extrapolation=allow_extrapolation,
        )

    def speed_of_sound(
        self,
        temperature: ArrayLike,
        pressure: ArrayLike,
        *,
        allow_extrapolation: bool = False,
    ) -> FloatResult:
        """Return equilibrium thermodynamic speed of sound in m/s."""
        return self._property(
            "speed_of_sound",
            temperature,
            pressure,
            allow_extrapolation=allow_extrapolation,
        )


AIR_HARMONIC_OSCILLATOR = HarmonicOscillatorGas(
    specific_gas_constant=287.05287,
    base_heat_capacity_ratio=1.4,
    modes=(VibrationalMode(1.0, 3055.56),),
    applicable_temperature_range=(400.0, 2000.0),
)
"""Air represented by one effective harmonic vibrational mode."""

AIR_BEATTIE_BRIDGEMAN = BeattieBridgemanGas(
    specific_gas_constant=287.05287,
    base_heat_capacity_ratio=1.4,
    a0=1.57222e2,
    b0=1.59228e-3,
    a=6.66814e-4,
    b=-3.80199e-4,
    c=1.49868e3,
    modes=(
        VibrationalMode(0.78088, 3394.3),
        VibrationalMode(0.20950, 2273.4),
    ),
    applicable_temperature_range=(38.888888888888886, 1222.2222222222222),
    applicable_pressure_range=(172.3689323292, 27579029.172672),
)
"""Dry-air model tabulated by R. E. Randall, AEDC-TR-57-8 (1957).

The 70--2200 degR and 0.025--4000 psia ranges identify Randall's tabulated
air-property range, not a strict physical-validity range of the equation of
state.  Set ``allow_extrapolation=True`` to evaluate beyond those tables.
"""


__all__ = [
    "AIR_BEATTIE_BRIDGEMAN",
    "AIR_HARMONIC_OSCILLATOR",
    "BeattieBridgemanGas",
    "HarmonicOscillatorGas",
    "ThermodynamicState",
    "VibrationalMode",
]
