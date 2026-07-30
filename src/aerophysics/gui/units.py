"""Display-unit conversion isolated from the SI calculation API."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike

QuantityKind = Literal[
    "angle",
    "density",
    "length",
    "pressure",
    "speed",
    "temperature",
]

_LENGTH_FACTORS = {"m": 1.0, "ft": 0.3048}
_SPEED_FACTORS = {"m/s": 1.0, "kt": 1852.0 / 3600.0}
_PRESSURE_FACTORS = {"Pa": 1.0, "psi": 6894.757293168}
_DENSITY_FACTORS = {"kg/m³": 1.0, "slug/ft³": 515.3788183931961}
_ANGLE_FACTORS = {"rad": 1.0, "deg": np.pi / 180.0}


@dataclass(frozen=True, slots=True)
class UnitPreferences:
    """Units selected for GUI input and output."""

    length: str = "m"
    speed: str = "m/s"
    pressure: str = "Pa"
    temperature: str = "K"
    density: str = "kg/m³"
    angle: str = "deg"

    def __post_init__(self) -> None:
        allowed = {
            "length": _LENGTH_FACTORS,
            "speed": _SPEED_FACTORS,
            "pressure": _PRESSURE_FACTORS,
            "temperature": {"K": 1.0, "°F": 1.0},
            "density": _DENSITY_FACTORS,
            "angle": _ANGLE_FACTORS,
        }
        for name, choices in allowed.items():
            value = getattr(self, name)
            if value not in choices:
                raise ValueError(f"unsupported {name} unit: {value}")

    def to_dict(self) -> dict[str, str]:
        """Return a JSON-serializable mapping."""
        return asdict(self)

    @classmethod
    def from_dict(cls, values: object) -> UnitPreferences:
        """Build preferences from an untrusted mapping."""
        if not isinstance(values, dict):
            raise ValueError("units must be an object")
        expected = {
            "length",
            "speed",
            "pressure",
            "temperature",
            "density",
            "angle",
        }
        if set(values) != expected or not all(
            isinstance(value, str) for value in values.values()
        ):
            raise ValueError("units contain missing or unsupported fields")
        return cls(**values)


def selected_unit(kind: QuantityKind, preferences: UnitPreferences) -> str:
    """Return the selected display unit for a quantity kind."""
    return str(getattr(preferences, kind))


def _factor(kind: QuantityKind, unit: str) -> float:
    tables = {
        "length": _LENGTH_FACTORS,
        "speed": _SPEED_FACTORS,
        "pressure": _PRESSURE_FACTORS,
        "density": _DENSITY_FACTORS,
        "angle": _ANGLE_FACTORS,
    }
    if kind == "temperature":
        if unit not in {"K", "°F"}:
            raise ValueError(f"unsupported temperature unit: {unit}")
        return 1.0
    try:
        return tables[kind][unit]
    except KeyError as error:
        raise ValueError(f"unsupported {kind} unit: {unit}") from error


def to_si(value: ArrayLike, kind: QuantityKind, unit: str) -> float | np.ndarray:
    """Convert a display value to the SI convention used by the core API."""
    values = np.asarray(value, dtype=np.float64)
    if np.any(~np.isfinite(values)):
        raise ValueError("display value must be finite")
    if kind == "temperature" and unit == "°F":
        converted = (values - 32.0) * (5.0 / 9.0) + 273.15
    else:
        converted = values * _factor(kind, unit)
    return float(converted) if converted.ndim == 0 else converted


def from_si(value: ArrayLike, kind: QuantityKind, unit: str) -> float | np.ndarray:
    """Convert an SI value to a selected display unit."""
    values = np.asarray(value, dtype=np.float64)
    if kind == "temperature" and unit == "°F":
        converted = (values - 273.15) * (9.0 / 5.0) + 32.0
    else:
        converted = values / _factor(kind, unit)
    return float(converted) if converted.ndim == 0 else converted
