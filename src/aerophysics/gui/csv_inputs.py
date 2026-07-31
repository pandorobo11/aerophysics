"""Validated CSV boundaries for GUI profile and frontal-shape inputs."""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass

import numpy as np

from aerophysics.gui.units import UnitPreferences, to_si


@dataclass(frozen=True, slots=True)
class ProfileCSV:
    wall_distance: np.ndarray
    velocity: np.ndarray
    density: np.ndarray


@dataclass(frozen=True, slots=True)
class ShapeCSV:
    height: np.ndarray
    width: np.ndarray


def _columns(data: bytes | str, required: tuple[str, ...]) -> dict[str, np.ndarray]:
    try:
        text = data.decode("utf-8-sig") if isinstance(data, bytes) else data
        reader = csv.DictReader(io.StringIO(text))
        if reader.fieldnames is None or not set(required) <= set(reader.fieldnames):
            raise ValueError(f"CSV must contain columns: {', '.join(required)}")
        rows = list(reader)
        if len(rows) < 2:
            raise ValueError("CSV must contain at least two data rows")
        result = {
            name: np.asarray([float(row[name]) for row in rows], dtype=np.float64)
            for name in required
        }
    except (UnicodeDecodeError, TypeError, ValueError) as error:
        if isinstance(error, ValueError) and str(error).startswith("CSV"):
            raise
        raise ValueError("CSV contains invalid numeric data") from error
    if any(np.any(~np.isfinite(values)) for values in result.values()):
        raise ValueError("CSV values must be finite")
    return result


def _wall_grid(values: np.ndarray, name: str) -> None:
    if values[0] != 0.0 or np.any(np.diff(values) <= 0.0):
        raise ValueError(f"{name} must start at zero and be strictly increasing")


def parse_profile_csv(data: bytes | str, units: UnitPreferences) -> ProfileCSV:
    """Parse profile values expressed in the currently selected units."""
    values = _columns(data, ("wall_distance", "velocity", "density"))
    height = np.asarray(
        to_si(values["wall_distance"], "length", units.length), dtype=np.float64
    )
    velocity = np.asarray(
        to_si(values["velocity"], "speed", units.speed), dtype=np.float64
    )
    density = np.asarray(
        to_si(values["density"], "density", units.density), dtype=np.float64
    )
    _wall_grid(height, "wall_distance")
    if np.any(velocity < 0.0):
        raise ValueError("velocity must be non-negative")
    if np.any(density <= 0.0):
        raise ValueError("density must be greater than zero")
    return ProfileCSV(height, velocity, density)


def parse_shape_csv(data: bytes | str, units: UnitPreferences) -> ShapeCSV:
    """Parse frontal width values expressed in the selected length unit."""
    values = _columns(data, ("height", "width"))
    height = np.asarray(
        to_si(values["height"], "length", units.length), dtype=np.float64
    )
    width = np.asarray(to_si(values["width"], "length", units.length), dtype=np.float64)
    _wall_grid(height, "height")
    if np.any(width < 0.0):
        raise ValueError("width must be non-negative")
    return ShapeCSV(height, width)


def profile_csv_template() -> str:
    """Return a minimal profile template in display units."""
    return "\ufeffwall_distance,velocity,density\n0,0,1.2\n0.01,100,1.2\n"


def shape_csv_template() -> str:
    """Return a minimal frontal-shape template in display units."""
    return "\ufeffheight,width\n0,0.01\n0.01,0\n"
