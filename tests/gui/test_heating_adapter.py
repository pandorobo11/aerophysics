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
    swept_cylinder_turbulent_heat_flux,
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
    for model in ("fay_riddell", "invalid"):
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


TURBULENT_INPUTS: dict[str, Any] = dict(
    upstream_mach=5.0,
    temperature=226.5,
    pressure=1197.0,
    radius=0.1,
    wall_temperature=300.0,
    geometry="swept_cylinder",
    sweep_angle=0.4,
    prandtl_number=0.8,
)
POLL_FIELDS = {
    "poll_turbulent_heat_flux": "heat_flux",
    "edge_temperature": "edge_temperature",
    "edge_pressure": "edge_pressure",
    "edge_mach": "edge_mach",
    "recovery_temperature": "recovery_temperature",
    "reference_temperature": "reference_temperature",
    "attachment_line_reynolds_number": "reynolds_number",
    "stanton_number": "stanton_number",
    "velocity_gradient": "velocity_gradient",
}
LEGACY_KEYS = {
    "geometric_altitude",
    "upstream_mach",
    "temperature",
    "pressure",
    "freestream_density",
    "freestream_velocity",
    "radius",
    "sweep_angle",
    "wall_temperature",
    "geometry",
    "model",
    "heat_flux",
    "fay_riddell_heat_flux",
    "tauber_heat_flux",
    "heat_flux_ratio",
    "stagnation_temperature",
    "stagnation_pressure",
    "total_enthalpy",
    "recovery_enthalpy",
    "velocity_gradient",
    "status",
    "message",
}


@pytest.mark.parametrize("model", ["poll_turbulent", "comparison"])
def test_heating_poll_core_and_comparison(model: str) -> None:
    row = heating_condition(**TURBULENT_INPUTS, model=model).rows[0]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        poll = swept_cylinder_turbulent_heat_flux(
            5, 226.5, 1197, 0.1, 0.4, wall_temperature=300, prandtl_number=0.8
        )
        tauber = swept_cylinder_heat_flux(
            5, 226.5, 1197, 0.1, 0.4, wall_temperature=300
        )
    for key, field in POLL_FIELDS.items():
        assert row[key] == pytest.approx(getattr(poll, field), rel=1e-12)
    assert row["freestream_density"] == pytest.approx(
        poll.freestream_density, rel=1e-12
    )
    assert row["freestream_velocity"] == pytest.approx(
        poll.freestream_velocity, rel=1e-12
    )
    if model == "comparison":
        assert row["heat_flux"] is None
        assert row["tauber_heat_flux"] == pytest.approx(tauber.heat_flux, rel=1e-12)
        assert row["turbulent_laminar_ratio"] == pytest.approx(
            poll.heat_flux / tauber.heat_flux, rel=1e-12
        )
        assert row["recovery_enthalpy"] == pytest.approx(tauber.recovery_enthalpy)
    else:
        assert row["heat_flux"] == pytest.approx(poll.heat_flux, rel=1e-12)
        assert row["turbulent_laminar_ratio"] is None
        assert row["tauber_heat_flux"] is None
        assert row["recovery_enthalpy"] is None


@pytest.mark.parametrize(
    "geometry,model",
    [
        ("sphere", "fay_riddell"),
        ("sphere", "tauber"),
        ("sphere", "comparison"),
        ("swept_cylinder", "tauber"),
    ],
)
def test_heating_legacy_row_and_column_contract(geometry: str, model: str) -> None:
    row = heating_condition(
        **{**TURBULENT_INPUTS, "geometry": geometry}, model=model
    ).rows[0]
    assert set(row) == LEGACY_KEYS
    from aerophysics.gui.tables import HEATING_COLUMNS

    expected = {
        column.heading(UnitPreferences())
        for column in HEATING_COLUMNS
        if column.key in LEGACY_KEYS
    }
    displayed = display_rows("aerodynamic_heating", (row,), UnitPreferences())[0]
    assert set(displayed) == expected
    figure = heating_trends((row,), UnitPreferences(), field="upstream_mach")["熱流束"]
    assert [trace.name for trace in figure.data] == (
        ["Fay\u2013Riddell", "Tauber"]
        if model == "comparison"
        else ["Fay\u2013Riddell"]
        if model == "fay_riddell"
        else ["Tauber"]
    )


@pytest.mark.parametrize("model", ["poll_turbulent", "comparison"])
@pytest.mark.parametrize("field", ["sweep_angle", "upstream_mach"])
def test_heating_poll_partial_domain_and_plots(model: str, field: str) -> None:
    inputs = {**TURBULENT_INPUTS, "geometric_altitude": None}
    if field == "sweep_angle":
        start, stop = 0.0, np.pi / 2
        invalid = [0, 4]
    else:
        inputs["sweep_angle"] = np.pi / 3
        start, stop = 1.1, 5.0
        invalid = [0]
    result = heating_sweep(
        **inputs, model=model, field=field, start=start, stop=stop, points=5
    )
    for i, row in enumerate(result.rows):
        assert cast(float, row["freestream_density"]) > 0
        assert cast(float, row["freestream_velocity"]) > 0
        if i in invalid:
            assert row["status"] == "out_of_domain"
            assert row["message"] == (
                "turbulent attachment line requires 0 < sweep angle < 90 deg "
                "and M cos(sweep) > 1"
            )
            assert all(row[key] is None for key in POLL_FIELDS)
            assert row["turbulent_laminar_ratio"] is None
        else:
            assert row["status"] == "ok" and row["message"] == ""
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                core = swept_cylinder_turbulent_heat_flux(
                    cast(float, row["upstream_mach"]),
                    226.5,
                    1197,
                    0.1,
                    cast(float, row["sweep_angle"]),
                    wall_temperature=300,
                    prandtl_number=0.8,
                )
            for key, core_field in POLL_FIELDS.items():
                assert row[key] == pytest.approx(getattr(core, core_field), rel=1e-12)
        if model == "comparison":
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                laminar = swept_cylinder_heat_flux(
                    cast(float, row["upstream_mach"]),
                    226.5,
                    1197,
                    0.1,
                    cast(float, row["sweep_angle"]),
                    wall_temperature=300,
                )
            assert row["tauber_heat_flux"] == pytest.approx(
                laminar.heat_flux, rel=1e-12
            )
            assert row["heat_flux"] is None
        else:
            assert row["heat_flux"] == row["poll_turbulent_heat_flux"]
    figure = heating_trends(result.rows, UnitPreferences(), field=field)["熱流束"]
    poll_trace = next(trace for trace in figure.data if trace.name == "Poll 乱流")
    assert np.isnan(np.asarray(poll_trace.y)[invalid]).all()
    assert np.isfinite(np.delete(np.asarray(poll_trace.y), invalid)).all()
    assert not poll_trace.connectgaps


@pytest.mark.parametrize("angle", [0.0, np.pi / 2, 1.5])
def test_heating_poll_entirely_invalid(angle: float) -> None:
    inputs = {**TURBULENT_INPUTS, "sweep_angle": angle}
    with pytest.raises(ValueError):
        heating_condition(**inputs, model="poll_turbulent")
    result = heating_condition(**inputs, model="comparison")
    assert result.rows[0]["status"] == "out_of_domain"
    assert all(result.rows[0][key] is None for key in POLL_FIELDS)
    assert result.rows[0]["tauber_heat_flux"] is not None


@pytest.mark.parametrize(
    "geometry,model",
    [
        ("sphere", "poll_turbulent"),
        ("swept_cylinder", "fay_riddell"),
        ("sphere", "unknown"),
        ("unknown", "poll_turbulent"),
    ],
)
def test_heating_invalid_geometry_model(geometry: str, model: str) -> None:
    with pytest.raises(
        ValueError, match=r"^invalid heating geometry and model combination$"
    ):
        heating_condition(**{**TURBULENT_INPUTS, "geometry": geometry}, model=model)


def test_heating_poll_warnings_and_units() -> None:
    result = heating_condition(
        **{**TURBULENT_INPUTS, "upstream_mach": np.array([10.0, 10.0])},
        model="poll_turbulent",
    )
    assert any("2000 K total temperature" in message for message in result.warnings)
    assert len(result.warnings) == len(set(result.warnings))
    row = result.rows[0]
    displayed = display_rows(
        "aerodynamic_heating",
        result.rows,
        UnitPreferences(heat_flux="MW/m²", temperature="°C", pressure="kPa"),
    )[0]
    for heading, key, factor in [
        ("乱流 Poll q_w [MW/m²]", "poll_turbulent_heat_flux", 1e-6),
        ("付着線外縁圧力 p_e [kPa]", "edge_pressure", 1e-3),
    ]:
        assert displayed[heading] == pytest.approx(cast(float, row[key]) * factor)
    for label, key in [
        ("付着線外縁温度 T_e", "edge_temperature"),
        ("回復温度 T_aw", "recovery_temperature"),
        ("基準温度 T*", "reference_temperature"),
    ]:
        assert displayed[label + " [°C]"] == pytest.approx(
            cast(float, row[key]) - 273.15
        )
    for label, key in [
        ("付着線外縁 Mach M_e", "edge_mach"),
        ("付着線 Reynolds R̄*", "attachment_line_reynolds_number"),
        ("Stanton St*", "stanton_number"),
    ]:
        assert displayed[label] == row[key]


def test_heating_poll_broadcast_shape() -> None:
    result = heating_condition(
        **{
            **TURBULENT_INPUTS,
            "upstream_mach": np.array([[5.0], [10.0]]),
            "sweep_angle": np.array([0.0, 0.4, np.pi / 2]),
        },
        model="poll_turbulent",
    )
    assert len(result.rows) == 6
    assert [row["status"] for row in result.rows] == [
        "out_of_domain",
        "ok",
        "out_of_domain",
    ] * 2
    for index, mach in [(1, 5.0), (4, 10.0)]:
        expected = heating_condition(
            **{**TURBULENT_INPUTS, "upstream_mach": mach}, model="poll_turbulent"
        ).rows[0]
        assert result.rows[index] == expected
