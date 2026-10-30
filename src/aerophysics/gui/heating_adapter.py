"""SI adapters for aerodynamic heating calculations."""

from __future__ import annotations

import warnings

import numpy as np

from aerophysics.atmosphere import standard_atmosphere
from aerophysics.exceptions import ApplicabilityWarning
from aerophysics.gui.adapters import CalculationResult, Row, sweep_values
from aerophysics.heating import (
    StagnationHeatingModel,
    sphere_stagnation_heat_flux,
    swept_cylinder_heat_flux,
)


def heating_condition(
    *,
    upstream_mach: float | np.ndarray,
    temperature: float | np.ndarray,
    pressure: float | np.ndarray,
    radius: float,
    wall_temperature: float,
    geometry: str,
    model: str,
    sweep_angle: float | np.ndarray = 0.0,
    prandtl_number: float = 0.72,
    geometric_altitude: float | np.ndarray | None = None,
) -> CalculationResult:
    """Calculate one or broadcast many heating conditions."""
    if (
        geometry not in {"sphere", "swept_cylinder"}
        or model not in {"fay_riddell", "tauber", "comparison"}
        or (geometry == "swept_cylinder" and model != "tauber")
    ):
        raise ValueError("invalid heating geometry and model combination")
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always", ApplicabilityWarning)
        fr = tauber = cylinder = None
        if geometry == "sphere":
            if model in {"fay_riddell", "comparison"}:
                fr = sphere_stagnation_heat_flux(
                    upstream_mach,
                    temperature,
                    pressure,
                    radius,
                    wall_temperature=wall_temperature,
                    model=StagnationHeatingModel.FAY_RIDDELL,
                    prandtl_number=prandtl_number,
                )
            if model in {"tauber", "comparison"}:
                tauber = sphere_stagnation_heat_flux(
                    upstream_mach,
                    temperature,
                    pressure,
                    radius,
                    wall_temperature=wall_temperature,
                    model=StagnationHeatingModel.TAUBER,
                )
        else:
            cylinder = swept_cylinder_heat_flux(
                upstream_mach,
                temperature,
                pressure,
                radius,
                sweep_angle,
                wall_temperature=wall_temperature,
            )
    result = fr or tauber or cylinder
    assert result is not None
    arrays = np.broadcast_arrays(
        np.asarray(upstream_mach),
        np.asarray(temperature),
        np.asarray(pressure),
        np.asarray(result.heat_flux),
        np.asarray(sweep_angle),
        np.asarray(geometric_altitude if geometric_altitude is not None else 0.0),
    )
    shape = arrays[0].shape

    def at(value: object, index: int) -> float | None:
        if value is None:
            return None
        return float(np.broadcast_to(np.asarray(value), shape).reshape(-1)[index])

    rows: list[Row] = []
    for i in range(arrays[0].size):
        f = at(fr.heat_flux, i) if fr else None
        t = (
            at(tauber.heat_flux, i)
            if tauber
            else at(cylinder.heat_flux, i)
            if cylinder
            else None
        )
        sphere = fr or tauber
        rows.append(
            {
                "geometric_altitude": at(geometric_altitude, i),
                "upstream_mach": at(upstream_mach, i),
                "temperature": at(temperature, i),
                "pressure": at(pressure, i),
                "freestream_density": at(result.freestream_density, i),
                "freestream_velocity": at(result.freestream_velocity, i),
                "radius": radius,
                "sweep_angle": at(sweep_angle, i) if cylinder else None,
                "wall_temperature": wall_temperature,
                "geometry": geometry,
                "model": model,
                "heat_flux": None if model == "comparison" else f if fr else t,
                "fay_riddell_heat_flux": f,
                "tauber_heat_flux": t,
                "heat_flux_ratio": f / t
                if f is not None and t not in (None, 0.0) and model == "comparison"
                else None,
                "stagnation_temperature": at(sphere.stagnation_temperature, i)
                if sphere
                else None,
                "stagnation_pressure": at(sphere.stagnation_pressure, i)
                if sphere
                else None,
                "total_enthalpy": at(sphere.stagnation_enthalpy, i) if sphere else None,
                "recovery_enthalpy": at(cylinder.recovery_enthalpy, i)
                if cylinder
                else None,
                "velocity_gradient": at(fr.velocity_gradient, i) if fr else None,
                "status": "ok",
                "message": "",
            }
        )
    messages = tuple(
        dict.fromkeys(
            str(w.message)
            for w in captured
            if issubclass(w.category, ApplicabilityWarning)
        )
    )
    return CalculationResult(tuple(rows), messages)


def heating_sweep(
    *,
    field: str,
    start: float,
    stop: float,
    points: int,
    upstream_mach: float | None,
    temperature: float | None,
    pressure: float | None,
    geometric_altitude: float | None,
    radius: float,
    wall_temperature: float,
    geometry: str,
    model: str,
    sweep_angle: float = 0.0,
    prandtl_number: float = 0.72,
) -> CalculationResult:
    """Calculate a Mach, geometric-altitude, or cylinder sweep-angle sweep."""
    values = sweep_values(start, stop, points)
    mach: float | np.ndarray
    temp: float | np.ndarray
    press: float | np.ndarray
    altitude: float | np.ndarray | None
    angle: float | np.ndarray = sweep_angle
    if field == "upstream_mach":
        if temperature is None or pressure is None:
            raise ValueError("Mach sweep requires temperature and pressure")
        mach = values
        temp, press, altitude = temperature, pressure, geometric_altitude
    elif field == "geometric_altitude":
        if upstream_mach is None:
            raise ValueError("altitude sweep requires Mach")
        atmosphere = standard_atmosphere(values)
        mach, temp, press, altitude = (
            upstream_mach,
            atmosphere.temperature,
            atmosphere.pressure,
            values,
        )
    elif field == "sweep_angle":
        if geometry != "swept_cylinder":
            raise ValueError("sweep-angle sweep requires a swept cylinder")
        if upstream_mach is None or temperature is None or pressure is None:
            raise ValueError(
                "sweep-angle sweep requires Mach, temperature, and pressure"
            )
        mach, temp, press, altitude = (
            upstream_mach,
            temperature,
            pressure,
            geometric_altitude,
        )
        angle = values
    else:
        raise ValueError("invalid heating sweep field")
    return heating_condition(
        upstream_mach=mach,
        temperature=temp,
        pressure=press,
        geometric_altitude=altitude,
        radius=radius,
        wall_temperature=wall_temperature,
        geometry=geometry,
        model=model,
        sweep_angle=angle,
        prandtl_number=prandtl_number,
    )
