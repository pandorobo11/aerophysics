"""Heating GUI adapters preserve the core calculation results."""

import warnings
from typing import Any, cast

import numpy as np
import pytest

from aerophysics.atmosphere import standard_atmosphere
from aerophysics.gui.figures import heating_trends
from aerophysics.gui.heating_adapter import heating_condition, heating_sweep
from aerophysics.gui.tables import display_rows
from aerophysics.gui.units import UnitPreferences
from aerophysics.heating import (
    StagnationHeatingModel,
    sphere_stagnation_heat_flux,
    swept_cylinder_heat_flux,
)


@pytest.mark.parametrize("model", ["fay_riddell", "tauber", "comparison"])
def test_sphere(model: str) -> None:
    kwargs: dict[str, Any] = dict(
        upstream_mach=10.0,
        temperature=226.5,
        pressure=1197.0,
        radius=0.1,
        wall_temperature=300.0,
        geometry="sphere",
        model=model,
    )
    result = heating_condition(**kwargs)
    row = result.rows[0]
    for name in ("fay_riddell", "tauber"):
        if model in (name, "comparison"):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                core = sphere_stagnation_heat_flux(
                    10,
                    226.5,
                    1197,
                    0.1,
                    wall_temperature=300,
                    model=StagnationHeatingModel(name),
                )
            assert row[name + "_heat_flux"] == pytest.approx(core.heat_flux, rel=1e-12)
    if model == "comparison":
        assert row["heat_flux"] is None
        assert row["heat_flux_ratio"] == pytest.approx(
            cast(float, row["fay_riddell_heat_flux"])
            / cast(float, row["tauber_heat_flux"]),
            rel=1e-12,
        )
    assert len(result.warnings) == len(set(result.warnings))


def test_cylinder_and_invalid_models() -> None:
    kwargs: dict[str, Any] = dict(
        upstream_mach=5.0,
        temperature=226.5,
        pressure=1197.0,
        radius=0.1,
        wall_temperature=300.0,
        geometry="swept_cylinder",
        sweep_angle=0.4,
    )
    result = heating_condition(**kwargs, model="tauber")
    with pytest.warns(Warning):
        core = swept_cylinder_heat_flux(5, 226.5, 1197, 0.1, 0.4, wall_temperature=300)
    assert result.rows[0]["heat_flux"] == pytest.approx(core.heat_flux, rel=1e-12)
    assert result.warnings
    for model in ("fay_riddell", "comparison", "invalid"):
        with pytest.raises(ValueError):
            heating_condition(**kwargs, model=model)


def test_sweeps_and_display() -> None:
    common: dict[str, Any] = dict(
        upstream_mach=10.0,
        temperature=226.5,
        pressure=1197.0,
        geometric_altitude=30000.0,
        radius=0.1,
        wall_temperature=300.0,
        geometry="sphere",
        model="comparison",
    )
    result = heating_sweep(field="upstream_mach", start=5, stop=10, points=3, **common)
    assert [row["upstream_mach"] for row in result.rows] == [5, 7.5, 10]
    assert len(result.warnings) == len(set(result.warnings))
    rows = display_rows(
        "aerodynamic_heating", result.rows, UnitPreferences(heat_flux="MW/m²")
    )
    assert rows[0]["Tauber q_w [MW/m²]"] == pytest.approx(
        cast(float, result.rows[0]["tauber_heat_flux"]) * 1e-6
    )
    figure = heating_trends(
        result.rows, UnitPreferences(heat_flux="MW/m²"), field="upstream_mach"
    )["熱流束"]
    assert len(figure.data) == 2
    assert "MW/m²" in figure.layout.yaxis.title.text
    altitude = heating_sweep(
        field="geometric_altitude", start=20000, stop=60000, points=3, **common
    )
    state = standard_atmosphere(np.array([20000.0, 40000.0, 60000.0]))
    assert [r["temperature"] for r in altitude.rows] == pytest.approx(state.temperature)
    assert [r["pressure"] for r in altitude.rows] == pytest.approx(state.pressure)
    cylinder = {**common, "geometry": "swept_cylinder", "model": "tauber"}
    angles = heating_sweep(
        field="sweep_angle", start=0.0, stop=np.pi / 2, points=3, **cylinder
    )
    expected = swept_cylinder_heat_flux(
        10.0,
        226.5,
        1197.0,
        0.1,
        np.array([0.0, np.pi / 4, np.pi / 2]),
        wall_temperature=300.0,
    )
    assert [r["sweep_angle"] for r in angles.rows] == pytest.approx(
        [0.0, np.pi / 4, np.pi / 2]
    )
    assert [r["heat_flux"] for r in angles.rows] == pytest.approx(expected.heat_flux)
    assert {r["upstream_mach"] for r in angles.rows} == {10.0}


@pytest.mark.parametrize(
    ("field", "missing", "message"),
    [
        (
            "upstream_mach",
            "temperature",
            "Mach sweep requires temperature and pressure",
        ),
        ("upstream_mach", "pressure", "Mach sweep requires temperature and pressure"),
        ("geometric_altitude", "upstream_mach", "altitude sweep requires Mach"),
        ("sweep_angle", None, "sweep-angle sweep requires a swept cylinder"),
        ("invalid", None, "invalid heating sweep field"),
    ],
)
def test_sweep_invalid_inputs(field: str, missing: str | None, message: str) -> None:
    inputs: dict[str, Any] = {
        "upstream_mach": 10.0,
        "temperature": 226.5,
        "pressure": 1197.0,
        "geometric_altitude": 30000.0,
        "radius": 0.1,
        "wall_temperature": 300.0,
        "geometry": "sphere",
        "model": "tauber",
    }
    if missing is not None:
        inputs[missing] = None
    with pytest.raises(ValueError, match=message):
        heating_sweep(field=field, start=5.0, stop=10.0, points=3, **inputs)


def test_old_nine_unit_configuration() -> None:
    values = UnitPreferences().to_dict()
    del values["heat_flux"]
    assert UnitPreferences.from_dict(values).heat_flux == "W/m²"


@pytest.mark.parametrize("missing", ["upstream_mach", "temperature", "pressure"])
def test_sweep_angle_requires_fixed_state(missing: str) -> None:
    inputs: dict[str, Any] = {
        "upstream_mach": 10.0,
        "temperature": 226.5,
        "pressure": 1197.0,
        "geometric_altitude": None,
        "radius": 0.1,
        "wall_temperature": 300.0,
        "geometry": "swept_cylinder",
        "model": "tauber",
        missing: None,
    }
    with pytest.raises(ValueError, match="requires Mach, temperature, and pressure"):
        heating_sweep(field="sweep_angle", start=0.0, stop=1.0, points=3, **inputs)
