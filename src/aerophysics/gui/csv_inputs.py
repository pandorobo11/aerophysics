"""Validated CSV boundaries for GUI profile and frontal-shape inputs."""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass

import numpy as np

from aerophysics.gui.units import UnitPreferences, from_si, to_si


@dataclass(frozen=True, slots=True)
class ProfileCSV:
    wall_distance: np.ndarray
    velocity: np.ndarray
    density: np.ndarray


@dataclass(frozen=True, slots=True)
class ShapeCSV:
    height: np.ndarray
    width: np.ndarray


def _columns(
    data: bytes | str, required: tuple[str, ...], optional: tuple[str, ...] = ()
) -> dict[str, np.ndarray]:
    try:
        text = (
            data.decode("utf-8-sig")
            if isinstance(data, bytes)
            else data.removeprefix("\ufeff")
        )
        reader = csv.DictReader(io.StringIO(text))
        if reader.fieldnames is None or not set(required) <= set(reader.fieldnames):
            raise ValueError(f"CSV must contain columns: {', '.join(required)}")
        rows = list(reader)
        if len(rows) < 2:
            raise ValueError("CSV must contain at least two data rows")
        result = {
            name: np.asarray([float(row[name]) for row in rows], dtype=np.float64)
            for name in (
                *required,
                *(name for name in optional if name in reader.fieldnames),
            )
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


@dataclass(frozen=True, slots=True)
class MeridianCSV:
    """Body meridian vertices in SI metres."""

    axial_distance: np.ndarray
    radius: np.ndarray


def parse_meridian_csv(data: bytes | str, units: UnitPreferences) -> MeridianCSV:
    """Parse x,r in display length units, leaving convexity to the solver."""
    values = _columns(data, ("x", "r"))
    x = np.asarray(to_si(values["x"], "length", units.length), dtype=np.float64)
    r = np.asarray(to_si(values["r"], "length", units.length), dtype=np.float64)
    if x[0] != 0 or r[0] != 0:
        raise ValueError("meridian must start at (0, 0)")
    if np.any(np.diff(x) <= 0):
        raise ValueError("x must be strictly increasing")
    if np.any(r[1:] <= 0):
        raise ValueError("r must be positive after the first point")
    return MeridianCSV(x, r)


def meridian_csv_template() -> str:
    """Return a cone-cylinder template in display length units with a BOM."""
    return "\ufeffx,r\n0,0\n0.3,0.05\n1.3,0.05\n"


@dataclass(frozen=True, slots=True)
class EdgeCSV:
    """Boundary-layer edge pressure and optional wall radius in SI units."""

    distance: np.ndarray
    pressure: np.ndarray
    radius: np.ndarray | None


def parse_edge_csv(data: bytes | str, units: UnitPreferences) -> EdgeCSV:
    """Parse s,p and an optional r column in display length and pressure units."""
    values = _columns(data, ("s", "p"), ("r",))
    distance = np.asarray(to_si(values["s"], "length", units.length), dtype=np.float64)
    pressure = np.asarray(
        to_si(values["p"], "pressure", units.pressure), dtype=np.float64
    )
    _wall_grid(distance, "s")
    if np.any(pressure <= 0.0):
        raise ValueError("p must be greater than zero")
    radius = None
    if "r" in values:
        radius = np.asarray(
            to_si(values["r"], "length", units.length), dtype=np.float64
        )
        if radius[0] < 0.0 or np.any(radius[1:] <= 0.0):
            raise ValueError("r must be non-negative at s = 0 and positive after it")
    return EdgeCSV(distance, pressure, radius)


_EDGE_TEMPLATE_SI = (
    (0.0, 10000.0),
    (0.25, 9800.0),
    (0.5, 9600.0),
    (0.75, 9400.0),
    (1.0, 9200.0),
)


def edge_csv_template(units: UnitPreferences) -> str:
    """Return a planar, mildly accelerating edge template in display units.

    The values are written in ``units`` because :func:`parse_edge_csv` reads
    an upload in the display units active when it arrives.
    """
    lines = (
        f"{float(from_si(s, 'length', units.length)):.12g},"
        f"{float(from_si(p, 'pressure', units.pressure)):.12g}"
        for s, p in _EDGE_TEMPLATE_SI
    )
    return "\ufeffs,p\n" + "".join(f"{line}\n" for line in lines)
