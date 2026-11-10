"""Axisymmetric body-flow adapter, replay and headless page regressions."""

from __future__ import annotations

import math
from typing import Any, cast
from unittest.mock import patch

import numpy as np
import pytest
import streamlit as st
from streamlit.testing.v1 import AppTest

from aerophysics.atmosphere import geopotential_to_geometric, standard_atmosphere
from aerophysics.characteristics import axisymmetric_body_flow
from aerophysics.exceptions import ShockConvergenceError
from aerophysics.gas import AIR
from aerophysics.gui.body_flow_adapter import (
    BodyFlowCalculation,
    body_flow_inputs_match,
    body_flow_net_rows,
    build_meridian,
    calculate_body_flow,
    resolve_freestream,
)
from aerophysics.gui.config import (
    ConfigurationError,
    dump_configuration,
    load_configuration,
    make_configuration,
)
from aerophysics.gui.csv_inputs import meridian_csv_template, parse_meridian_csv
from aerophysics.gui.figures import body_flow_figures
from aerophysics.gui.tables import display_rows, rows_to_csv
from aerophysics.gui.units import UnitPreferences

_SCRIPT = """
import streamlit as st
from aerophysics.gui.body_flow_page import render_axisymmetric_body
from aerophysics.gui.units import UnitPreferences
render_axisymmetric_body(UnitPreferences(**st.session_state.get("test_units", {})))
"""


@pytest.mark.parametrize(
    ("current", "saved", "expected"),
    [
        (0.3, 0.3 * (1 + 5e-13), True),
        (0.3, 0.3 * (1 + 2e-12), False),
        (0.0, 1e-30, False),
        ([0.0, 0.3, 1.3], [0.0, 0.3 * (1 + 5e-13), 1.3], True),
        ([0.0, 0.3], [0.0, 0.3, 1.3], False),
        ([0.0, 0.3], [0.0, 0.31], False),
        (88, 88, True),
        (88, 89, False),
        (88, 88.0, False),
        (None, None, True),
        (None, 0.0, False),
        ("cone", "cone", True),
        ("cone", "tangent_ogive", False),
        (True, True, True),
        (True, False, False),
        (True, 1, False),
    ],
)
def test_input_comparison(current: object, saved: object, expected: bool) -> None:
    assert body_flow_inputs_match({"value": current}, {"value": saved}) is expected
    assert body_flow_inputs_match({"value": saved}, {"value": current}) is expected


def test_input_comparison_checks_field_names() -> None:
    assert not body_flow_inputs_match({"mach": 3.0}, {"radius": 3.0})
    assert not body_flow_inputs_match({"mach": 3.0}, {"mach": 3.0, "altitude": None})


@pytest.mark.parametrize(
    ("nose", "boattail"),
    [("cone", False), ("cone", True), ("tangent_ogive", True)],
)
def test_display_unit_changes_do_not_mark_inputs_changed(
    nose: str, boattail: bool
) -> None:
    app = AppTest.from_string(_SCRIPT, default_timeout=30).run()
    if nose != "cone":
        app.selectbox(key="body_flow_nose").set_value(nose).run()
    if boattail:
        app.checkbox(key="body_flow_boattail").check().run()
    app.button(key="body_flow_calculate").click().run()
    assert not app.exception and not app.error and not app.warning
    payload = app.session_state["body_flow_payload"]
    with patch(
        "aerophysics.gui.body_flow_page.calculate_body_flow",
        side_effect=AssertionError("must not recompute"),
    ):
        for units in (
            {"length": "ft"},
            {"length": "in", "pressure": "psf", "temperature": "°R", "angle": "rad"},
        ):
            app.session_state["test_units"] = units
            app.run()
            assert not app.exception and not app.error and not app.warning
            assert app.session_state["body_flow_payload"] is payload
        app.number_input(key="body_flow_mach").set_value(4.0).run()
        assert not app.exception and not app.error
        assert any("入力が変更されています" in warning.value for warning in app.warning)
        assert app.session_state["body_flow_payload"] is payload


def _configuration(**overrides: Any) -> dict[str, Any]:
    inputs = {
        "mach": 3.0,
        "pressure": 101325.0,
        "temperature": 288.15,
        "altitude": None,
        "nose_length": 0.3,
        "base_radius": 0.05,
        "ogive_segments": 32,
        "cylinder_length": 1.0,
        "boattail_length": None,
        "boattail_radius": None,
        "meridian_x": None,
        "meridian_r": None,
        "characteristic_count": 88,
        "fan_angle_step": math.pi / 180,
    }
    inputs.update(overrides)
    return make_configuration(
        calculator="axisymmetric_body",
        mode="single",
        inputs_si=inputs,
        models={
            "geometry_source": "parametric",
            "nose": "cone",
            "boattail": False,
            "freestream_basis": "static",
            "altitude_basis": "geometric",
            "pressure_quantity": "cp",
        },
        units=UnitPreferences(),
    )


@pytest.fixture(scope="module")
def calculation() -> BodyFlowCalculation:
    return calculate_body_flow(build_meridian())


def test_default_adapter_matches_public_solver_exactly(
    calculation: BodyFlowCalculation,
) -> None:
    direct = axisymmetric_body_flow(
        [0, 0.3, 1.3], [0, 0.05, 0.05], mach=3, pressure=101325, temperature=288.15
    )
    for field in (
        "axial_distance",
        "radius",
        "surface_distance",
        "pressure",
        "temperature",
        "mach",
        "flow_angle",
        "shock_axial_distance",
        "shock_radius",
        "shock_angle",
    ):
        np.testing.assert_array_equal(
            getattr(calculation.raw, field), getattr(direct, field)
        )
    assert calculation.raw.total_pressure == direct.total_pressure
    assert calculation.raw.total_temperature == direct.total_temperature
    assert calculation.raw.gas == AIR
    for key in (
        "axial_distance",
        "radius",
        "surface_distance",
        "pressure",
        "temperature",
        "mach",
        "flow_angle",
    ):
        np.testing.assert_array_equal(
            [row[key] for row in calculation.surface.rows], getattr(direct, key)
        )
    for key in ("shock_axial_distance", "shock_radius", "shock_angle"):
        np.testing.assert_array_equal(
            [row[key] for row in calculation.shock_rows], getattr(direct, key)
        )


def test_pressure_formulas(calculation: BodyFlowCalculation) -> None:
    for i in (0, len(calculation.surface.rows) // 2, -1):
        row = calculation.surface.rows[i]
        pressure = calculation.raw.pressure[i]
        assert row["pressure_ratio"] == pressure / 101325.0
        assert row["pressure_rise"] == (pressure - 101325.0) / 101325.0
        assert row["cp"] == (pressure - 101325.0) / (
            0.5 * AIR.heat_capacity_ratio * 101325.0 * 3.0**2
        )


def test_ogive_arc_and_appended_segments() -> None:
    meridian = build_meridian(
        nose="tangent_ogive",
        cylinder_length=1,
        boattail=True,
        boattail_length=0.2,
        boattail_radius=0.025,
    )
    x, r = meridian.axial_distance, meridian.radius
    assert (x[0], r[0]) == (0, 0)
    assert (x[32], r[32]) == (0.3, 0.05)
    rho = (0.05**2 + 0.3**2) / (2 * 0.05)
    np.testing.assert_allclose(
        (x[:33] - 0.3) ** 2 + (r[:33] - 0.05 + rho) ** 2, rho**2, rtol=0, atol=5e-16
    )
    angles = np.arctan2(np.diff(r), np.diff(x))
    assert np.all(np.diff(angles) < 0)
    np.testing.assert_array_equal(x[-2:], [1.3, 1.5])
    np.testing.assert_array_equal(r[-2:], [0.05, 0.025])
    no_cylinder = build_meridian(nose="tangent_ogive", cylinder_length=0)
    assert len(no_cylinder.radius) == 33
    assert (no_cylinder.axial_distance[-1], no_cylinder.radius[-1]) == (0.3, 0.05)
    tail = build_meridian(
        cylinder_length=0, boattail=True, boattail_length=0.2, boattail_radius=0.025
    )
    np.testing.assert_array_equal(tail.axial_distance, [0, 0.3, 0.5])


@pytest.mark.parametrize("radius", [0, -0.1, 0.05, 0.06, math.nan, None])
def test_invalid_boattail_radius(radius: float | None) -> None:
    with pytest.raises(ValueError, match="0 < boattail_radius < base_radius"):
        build_meridian(boattail=True, boattail_length=0.2, boattail_radius=radius)


@pytest.mark.parametrize("basis", ["geometric", "geopotential"])
def test_atmosphere_freestream(basis: str) -> None:
    pressure, temperature = resolve_freestream(
        freestream_basis="atmosphere", altitude=10000, altitude_basis=basis
    )
    state = standard_atmosphere(
        10000 if basis == "geometric" else geopotential_to_geometric(10000)
    )
    assert pressure == state.pressure
    assert temperature == state.temperature
    assert resolve_freestream(pressure=12345, temperature=250) == (12345, 250)


@pytest.mark.parametrize("csv", [False, True])
def test_settings_roundtrip(csv: bool) -> None:
    configuration = _configuration()
    if csv:
        configuration["models"]["geometry_source"] = "csv"
        configuration["inputs_si"]["meridian_x"] = [0, 0.3, 1.3]
        configuration["inputs_si"]["meridian_r"] = [0, 0.05, 0.05]
    assert load_configuration(dump_configuration(configuration)) == configuration


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("characteristic_count", 3),
        ("characteristic_count", 1001),
        ("characteristic_count", 88.5),
        ("ogive_segments", 1),
        ("ogive_segments", 201),
        ("fan_angle_step", math.pi / 7200),
        ("fan_angle_step", math.pi / 18),
    ],
)
def test_settings_resolution_limits(field: str, value: float) -> None:
    configuration = _configuration()
    configuration["inputs_si"][field] = value
    with pytest.raises(ConfigurationError):
        dump_configuration(configuration)


def test_settings_single_only() -> None:
    configuration = _configuration()
    configuration["mode"] = "sweep"
    with pytest.raises(ConfigurationError):
        dump_configuration(configuration)


def test_figures_tables_and_display_units(calculation: BodyFlowCalculation) -> None:
    units = UnitPreferences(length="ft", pressure="kPa", temperature="°C", angle="rad")
    rows = display_rows("axisymmetric_body", calculation.surface.rows, units)
    assert rows[0]["静圧 p [kPa]"] == calculation.raw.pressure[0] / 1000
    assert rows[-1]["軸方向距離 x [ft]"] == 1.3 / 0.3048
    assert rows[0]["静温 T [°C]"] == calculation.raw.temperature[0] - 273.15
    assert "Cp" in rows_to_csv(rows)
    shocks = display_rows("axisymmetric_body_shock", calculation.shock_rows, units)
    assert shocks[0]["衝撃波角 β [rad]"] == calculation.raw.shock_angle[0]
    figures = body_flow_figures(calculation.surface.rows, calculation.shock_rows, units)
    assert figures["形状と衝撃波"].layout.yaxis.scaleanchor == "x"
    body, shock = figures["形状と衝撃波"].data
    assert body.fill == "toself" and shock.name == "衝撃波"
    radius = np.array([r["radius"] for r in calculation.surface.rows]) / 0.3048
    np.testing.assert_array_equal(body.y, np.concatenate((radius, -radius[::-1])))
    assert figures["形状と衝撃波"].layout.shapes[0].line.dash == "dash"
    np.testing.assert_array_equal(
        figures["表面圧"].data[0].y, [r["cp"] for r in calculation.surface.rows]
    )
    figures = body_flow_figures(
        calculation.surface.rows,
        calculation.shock_rows,
        units,
        pressure_quantity="pressure_rise",
    )
    np.testing.assert_array_equal(
        figures["表面圧"].data[0].y,
        [r["pressure_rise"] for r in calculation.surface.rows],
    )


def test_page_calculate_switch_save_replay_and_reset() -> None:
    app = AppTest.from_string(_SCRIPT, default_timeout=30).run()
    assert not app.exception and not app.error
    assert app.title[0].value == "軸対称物体流れ\uff08特性曲線法\uff09"
    assert not app.metric
    assert app.number_input(key="body_flow_characteristic_count").value == 88
    assert app.number_input(key="body_flow_fan_angle_step").value == 1.0
    app.button(key="body_flow_calculate").click().run()
    assert not app.exception and not app.error
    assert len(app.metric) == 6 and len(app.dataframe) == 2
    assert len(app.get("plotly_chart")) == 4
    assert "特性曲線網" in [tab.label for tab in app.tabs]
    calculation, config = app.session_state["body_flow_payload"]
    with patch(
        "aerophysics.gui.body_flow_page.calculate_body_flow",
        side_effect=AssertionError("must not recompute"),
    ):
        app.selectbox(key="body_flow_pressure_quantity").set_value(
            "pressure_rise"
        ).run()
        assert not app.exception and not app.warning
        assert app.session_state["body_flow_payload"][0] is calculation
        for coloring in ("mach", "cp", "none"):
            app.selectbox(key="body_flow_net_coloring").set_value(coloring).run()
            assert not app.exception and not app.warning
            assert app.session_state["body_flow_payload"][0] is calculation
        app.button(key="body_flow_save_case").click().run()
        assert app.session_state["current_axisymmetric_body_case"] is calculation.raw
        assert app.success
        app.number_input(key="body_flow_mach").set_value(4).run()
        assert app.warning and app.session_state["body_flow_payload"][0] is calculation
    replay = AppTest.from_string(_SCRIPT, default_timeout=30)
    replay.session_state["pending_axisymmetric_body_configuration"] = (
        load_configuration(dump_configuration(config))
    )
    replay.run()
    replay.button(key="body_flow_calculate").click().run()
    assert not replay.exception and not replay.error
    np.testing.assert_array_equal(
        replay.session_state["body_flow_payload"][0].raw.pressure,
        calculation.raw.pressure,
    )
    app.button(key="body_flow_reset").click().run()
    assert not app.exception and not app.metric
    assert app.number_input(key="body_flow_mach").value == 3


def test_csv_settings_restore_without_upload() -> None:
    configuration = _configuration()
    configuration["models"]["geometry_source"] = "csv"
    vertices = parse_meridian_csv(meridian_csv_template(), UnitPreferences())
    configuration["inputs_si"]["meridian_x"] = vertices.axial_distance.tolist()
    configuration["inputs_si"]["meridian_r"] = vertices.radius.tolist()
    app = AppTest.from_string(_SCRIPT, default_timeout=30)
    app.session_state["pending_axisymmetric_body_configuration"] = load_configuration(
        dump_configuration(configuration)
    )
    app.run()
    app.run()  # Embedded arrays survive another rerun without the uploader.
    app.button(key="body_flow_calculate").click().run()
    assert not app.exception and not app.error
    assert app.session_state["body_flow_payload"][1]["inputs_si"]["meridian_x"] == [
        0,
        0.3,
        1.3,
    ]


@pytest.mark.parametrize("concave", [True, False])
def test_solver_errors_clear_stale_payload(
    concave: bool, calculation: BodyFlowCalculation
) -> None:
    configuration = _configuration(mach=1.2, nose_length=0.1, base_radius=0.1)
    if concave:
        configuration = _configuration()
        configuration["models"]["geometry_source"] = "csv"
        configuration["inputs_si"].update(
            meridian_x=[0, 0.3, 1.3], meridian_r=[0, 0.05, 0.3]
        )
    app = AppTest.from_string(_SCRIPT, default_timeout=30)
    app.session_state["pending_axisymmetric_body_configuration"] = load_configuration(
        dump_configuration(configuration)
    )
    app.run()
    app.session_state["body_flow_payload"] = (calculation, _configuration())
    app.button(key="body_flow_calculate").click().run()
    assert not app.exception and app.error
    if concave:
        assert "concave corner" in app.error[0].value
    else:
        assert "shock" in app.error[0].value.lower()
    assert "body_flow_payload" not in app.session_state
    assert not app.dataframe


def test_convergence_error_preserves_message() -> None:
    app = AppTest.from_string(_SCRIPT, default_timeout=30).run()
    with patch(
        "aerophysics.gui.body_flow_page.calculate_body_flow",
        side_effect=ShockConvergenceError("failed at x=0.42"),
    ):
        app.button(key="body_flow_calculate").click().run()
    assert not app.exception
    assert app.error[0].value == "failed at x=0.42"
    assert "body_flow_payload" not in app.session_state


def test_page_ogive_atmosphere_and_unit_changes() -> None:
    app = AppTest.from_string(_SCRIPT, default_timeout=30).run()
    app.selectbox(key="body_flow_nose").set_value("tangent_ogive").run()
    assert app.number_input(key="body_flow_ogive_segments").value == 32
    assert "正確な先端半頂角" in app.info[0].value
    app.selectbox(key="body_flow_freestream_basis").set_value("atmosphere").run()
    app.selectbox(key="body_flow_altitude_basis").set_value("geopotential").run()
    app.session_state["test_units"] = {
        "length": "ft",
        "pressure": "kPa",
        "temperature": "°C",
        "angle": "rad",
    }
    app.run()
    assert not app.exception
    assert app.number_input(key="body_flow_nose_length").value == pytest.approx(
        0.3 / 0.3048
    )
    assert app.number_input(key="body_flow_fan_angle_step").value == math.pi / 180
    assert app.number_input(key="body_flow_altitude").value == pytest.approx(
        10000 / 0.3048
    )
    app.button(key="body_flow_calculate").click().run()
    assert not app.exception and not app.error
    configuration = app.session_state["body_flow_payload"][1]
    pressure, temperature = resolve_freestream(
        freestream_basis="atmosphere", altitude=10000, altitude_basis="geopotential"
    )
    assert configuration["inputs_si"]["pressure"] == pressure
    assert configuration["inputs_si"]["temperature"] == temperature


@pytest.mark.parametrize(("x", "r"), [([], []), ([0], [0]), ([0, 0.3], [0])])
def test_malformed_embedded_vertices_show_error(x: list[float], r: list[float]) -> None:
    configuration = _configuration(meridian_x=x, meridian_r=r)
    configuration["models"]["geometry_source"] = "csv"
    app = AppTest.from_string(_SCRIPT, default_timeout=30)
    app.session_state["pending_axisymmetric_body_configuration"] = load_configuration(
        dump_configuration(configuration)
    )
    app.run()
    assert not app.exception
    app.button(key="body_flow_calculate").click().run()
    assert not app.exception and app.error
    assert "body_flow_payload" not in app.session_state


def _upload(app: AppTest, name: str, data: bytes) -> None:
    uploader: Any = next(
        item
        for item in app.get("file_uploader")
        if item.key == "body_flow_meridian_upload"
    )
    uploader.upload(name, data, "text/csv").run()


def test_uploaded_csv_keeps_its_si_shape_across_unit_changes() -> None:
    app = AppTest.from_string(_SCRIPT, default_timeout=30).run()
    app.selectbox(key="body_flow_geometry_source").set_value("csv").run()
    _upload(app, "meridian.csv", meridian_csv_template().encode("utf-8"))
    app.button(key="body_flow_calculate").click().run()
    assert not app.exception and not app.error and not app.warning
    payload = app.session_state["body_flow_payload"]
    expected_x, expected_r = [0.0, 0.3, 1.3], [0.0, 0.05, 0.05]
    with patch(
        "aerophysics.gui.body_flow_page.calculate_body_flow",
        side_effect=AssertionError("must not recompute"),
    ):
        for units in ({"length": "mm"}, {"length": "ft"}):
            app.session_state["test_units"] = units
            app.run()
            assert not app.exception and not app.error and not app.warning
            assert app.session_state["body_flow_payload"] is payload
            cached = app.session_state["body_flow_meridian_upload_si"]
            assert (cached[1], cached[2]) == (expected_x, expected_r)
    app.button(key="body_flow_calculate").click().run()
    assert not app.exception and not app.error and not app.warning
    inputs = app.session_state["body_flow_payload"][1]["inputs_si"]
    assert (inputs["meridian_x"], inputs["meridian_r"]) == (expected_x, expected_r)
    # A replacement file is read in the length unit active when it arrives.
    _upload(app, "feet.csv", b"x,r\n0,0\n1,0.1\n4,0.1\n")
    app.button(key="body_flow_calculate").click().run()
    assert not app.exception and not app.error
    inputs = app.session_state["body_flow_payload"][1]["inputs_si"]
    assert inputs["meridian_x"] == pytest.approx([0.0, 0.3048, 1.2192])
    assert inputs["meridian_r"] == pytest.approx([0.0, 0.03048, 0.03048])
    app.button(key="body_flow_reset").click().run()
    assert "body_flow_meridian_upload_si" not in app.session_state


def test_csv_upload_and_invalid_boattail() -> None:
    app = AppTest.from_string(_SCRIPT, default_timeout=30).run()
    app.checkbox(key="body_flow_boattail").check().run()
    app.number_input(key="body_flow_boattail_radius").set_value(0.06).run()
    app.button(key="body_flow_calculate").click().run()
    assert not app.exception
    assert "0 < boattail_radius < base_radius" in app.error[0].value
    app.selectbox(key="body_flow_geometry_source").set_value("csv").run()
    # Streamlit exposes file upload through the testing element's upload method.
    uploader: Any = next(
        item
        for item in app.get("file_uploader")
        if item.key == "body_flow_meridian_upload"
    )
    uploader.upload(
        "meridian.csv", meridian_csv_template().encode("utf-8"), "text/csv"
    ).run()
    app.button(key="body_flow_calculate").click().run()
    assert not app.exception and not app.error
    assert app.session_state["body_flow_payload"][1]["inputs_si"]["meridian_x"] == [
        0,
        0.3,
        1.3,
    ]
    uploader = next(
        item
        for item in app.get("file_uploader")
        if item.key == "body_flow_meridian_upload"
    )
    uploader.upload("bad.csv", b"x,r\n0,0\n1,0\n", "text/csv").run()
    app.button(key="body_flow_calculate").click().run()
    assert not app.exception and app.error
    assert "positive" in app.error[0].value
    assert "body_flow_payload" not in app.session_state


def test_unrelated_runtime_errors_propagate() -> None:
    app = AppTest.from_string(_SCRIPT, default_timeout=30).run()
    with patch(
        "aerophysics.gui.body_flow_page.calculate_body_flow",
        side_effect=RuntimeError("unexpected failure"),
    ):
        app.button(key="body_flow_calculate").click().run()
    assert app.exception[0].message == "unexpected failure"
    assert not app.error


@pytest.mark.parametrize("coloring", ["none", "mach", "cp"])
def test_net_settings_roundtrip_and_legacy(coloring: str) -> None:
    configuration = _configuration()
    configuration["models"]["net_coloring"] = coloring
    assert load_configuration(dump_configuration(configuration)) == configuration
    del configuration["models"]["net_coloring"]
    assert (
        cast(
            dict[str, Any],
            load_configuration(dump_configuration(configuration))["models"],
        )["net_coloring"]
        == "none"
    )
    configuration["models"]["net_coloring"] = "temperature"
    with pytest.raises(ConfigurationError, match="net_coloring"):
        dump_configuration(configuration)


def test_net_csv_all_nodes_and_download(calculation: BodyFlowCalculation) -> None:
    units = UnitPreferences(length="ft", pressure="kPa", angle="deg")
    net = calculation.raw.characteristic_net
    assert net is not None
    rows = body_flow_net_rows(calculation.raw, pressure=101325, mach=3)
    assert len(rows) == sum(len(row.radius) for row in net.rows)
    displayed = display_rows("axisymmetric_body_net", rows, units)
    csv = rows_to_csv(displayed)
    assert (
        csv.lstrip("\ufeff").splitlines()[0] == "行番号,節点番号,行種別,衝撃波で終端,"
        "軸方向距離 x [ft],半径 r [ft],流れ角 θ [deg],"
        "静圧 p [kPa],全圧 p₀ [kPa],Mach M,Cp"
    )
    assert csv.splitlines()[1].split(",")[:4] == ["0", "0", "initial", "True"]
    assert displayed[0]["軸方向距離 x [ft]"] == net.rows[0].axial_distance[0] / 0.3048
    assert displayed[0]["流れ角 θ [deg]"] == pytest.approx(
        math.degrees(net.rows[0].flow_angle[0])
    )
    assert displayed[0]["静圧 p [kPa]"] == net.rows[0].pressure[0] / 1000
    assert displayed[0]["全圧 p₀ [kPa]"] == net.rows[0].total_pressure[0] / 1000
    assert rows[0]["cp"] == (net.rows[0].pressure[0] - 101325) / (
        0.5 * 1.4 * 101325 * 9
    )
    app = AppTest.from_string(_SCRIPT, default_timeout=30)
    app.session_state["test_units"] = units.to_dict()
    app.run()
    with (
        patch(
            "aerophysics.gui.body_flow_page.calculate_body_flow",
            return_value=calculation,
        ),
        patch(
            "aerophysics.gui.body_flow_page.st.download_button",
            wraps=st.download_button,
        ) as download,
    ):
        app.button(key="body_flow_calculate").click().run()
    assert not app.exception and not app.error
    call = next(
        call
        for call in download.call_args_list
        if call.kwargs.get("key") == "body_flow_net_csv"
    )
    assert call.args[0] == "特性曲線網CSVをダウンロード"
    assert call.args[1] == csv


def test_fine_net_display_limits_and_caption() -> None:
    calculation = calculate_body_flow(build_meridian(), characteristic_count=200)
    net = calculation.raw.characteristic_net
    assert net is not None
    units = UnitPreferences(length="ft")
    for coloring in ("none", "mach", "cp"):
        figure = body_flow_figures(
            calculation.surface.rows,
            calculation.shock_rows,
            units,
            characteristic_net=net,
            net_coloring=coloring,
        )["特性曲線網"]
        assert figure.layout.yaxis.scaleanchor == "x"
        assert figure.layout.meta["downsampled"]
        plus, minus = figure.data[2:4]
        assert plus.name == "C+" and minus.name == "C\N{MINUS SIGN}"
        assert plus.line.color != minus.line.color
        assert sum(value is None for value in plus.x) == min(150, len(net.rows))
        assert sum(value is None for value in minus.x) == min(
            80, len(net.minus_characteristics())
        )
        assert plus.x[0] == net.rows[0].axial_distance[0] / 0.3048
        assert plus.x[-2] == net.rows[-1].axial_distance[-1] / 0.3048
        if coloring == "none":
            assert len(figure.data) == 4
        else:
            markers = figure.data[4]
            assert markers.type == "scattergl" and markers.marker.showscale
            assert len(markers.x) <= 20000
            rows = [
                net.rows[i]
                for i in np.linspace(
                    0, len(net.rows) - 1, min(150, len(net.rows)), dtype=int
                )
            ]
            values = np.concatenate(
                [
                    row.mach
                    if coloring == "mach"
                    else (row.pressure - 101325) / (0.5 * 1.4 * 101325 * 9)
                    for row in rows
                ]
            )
            indices = np.linspace(
                0, len(values) - 1, min(20000, len(values)), dtype=int
            )
            np.testing.assert_array_equal(markers.marker.color, values[indices])
    app = AppTest.from_string(_SCRIPT, default_timeout=30).run()
    app.number_input(key="body_flow_characteristic_count").set_value(200).run()
    with patch(
        "aerophysics.gui.body_flow_page.calculate_body_flow", return_value=calculation
    ):
        app.button(key="body_flow_calculate").click().run()
    assert not app.exception and not app.error
    assert len(app.get("plotly_chart")) == 4
    assert any("表示用に間引いています" in caption.value for caption in app.caption)


def test_net_marker_cap_is_exercised() -> None:
    calculation = calculate_body_flow(build_meridian(), characteristic_count=400)
    net = calculation.raw.characteristic_net
    assert net is not None
    figure = body_flow_figures(
        calculation.surface.rows,
        calculation.shock_rows,
        UnitPreferences(),
        characteristic_net=net,
        net_coloring="mach",
    )["特性曲線網"]
    assert sum(value is not None for value in figure.data[2].x) > 20000
    assert len(figure.data[4].x) == 20000
    assert figure.layout.meta["downsampled"]


def test_net_csv_cache_reuses_and_rebuilds() -> None:
    app = AppTest.from_string(_SCRIPT, default_timeout=30).run()
    with patch(
        "aerophysics.gui.body_flow_page.body_flow_net_rows", wraps=body_flow_net_rows
    ) as build_csv_rows:
        app.button(key="body_flow_calculate").click().run()
        assert not app.exception and not app.error
        assert build_csv_rows.call_count == 1
        cached = app.session_state["body_flow_net_csv_cache"]
        calculation = app.session_state["body_flow_payload"][0]
        assert cached[0] is calculation
        with patch(
            "aerophysics.gui.body_flow_page.calculate_body_flow",
            side_effect=AssertionError("must not recompute"),
        ):
            for key, value in (
                ("body_flow_net_coloring", "mach"),
                ("body_flow_net_coloring", "cp"),
                ("body_flow_pressure_quantity", "pressure_rise"),
            ):
                app.selectbox(key=key).set_value(value).run()
                assert not app.exception and not app.error and not app.warning
                assert build_csv_rows.call_count == 1
                assert app.session_state["body_flow_net_csv_cache"] is cached
            app.number_input(key="body_flow_mach").set_value(4.0).run()
            assert not app.exception and not app.error and app.warning
            assert build_csv_rows.call_count == 1
            assert app.session_state["body_flow_net_csv_cache"] is cached
            app.session_state["test_units"] = {"length": "ft", "pressure": "kPa"}
            app.run()
            assert not app.exception and not app.error
            assert build_csv_rows.call_count == 2
            converted = app.session_state["body_flow_net_csv_cache"]
            assert converted[0] is calculation
            assert converted[2] != cached[2]
            assert "軸方向距離 x [ft]" in converted[2]
            app.selectbox(key="body_flow_net_coloring").set_value("none").run()
            assert not app.exception and not app.error
            assert build_csv_rows.call_count == 2
            assert app.session_state["body_flow_net_csv_cache"] is converted
        app.button(key="body_flow_calculate").click().run()
        assert not app.exception and not app.error and not app.warning
        assert build_csv_rows.call_count == 3
        updated = app.session_state["body_flow_net_csv_cache"]
        assert updated[0] is app.session_state["body_flow_payload"][0]
        assert updated[0] is not calculation
        assert updated[2] != converted[2]
        with patch(
            "aerophysics.gui.body_flow_page.calculate_body_flow",
            side_effect=ValueError("calculation failed"),
        ):
            app.button(key="body_flow_calculate").click().run()
        assert not app.exception and app.error
        assert "body_flow_payload" not in app.session_state
        assert "body_flow_net_csv_cache" not in app.session_state
        assert build_csv_rows.call_count == 3


@pytest.mark.parametrize("action", ["reset", "import", "remove_payload"])
def test_net_csv_cache_cleared_with_payload(action: str) -> None:
    app = AppTest.from_string(_SCRIPT, default_timeout=30).run()
    app.button(key="body_flow_calculate").click().run()
    assert not app.exception and not app.error
    assert "body_flow_net_csv_cache" in app.session_state
    with patch(
        "aerophysics.gui.body_flow_page.body_flow_net_rows",
        side_effect=AssertionError("no payload to export"),
    ):
        if action == "reset":
            app.button(key="body_flow_reset").click().run()
        elif action == "import":
            app.session_state["pending_axisymmetric_body_configuration"] = (
                _configuration()
            )
            app.run()
        else:
            del app.session_state["body_flow_payload"]
            app.run()
    assert not app.exception and not app.error
    assert "body_flow_payload" not in app.session_state
    assert "body_flow_net_csv_cache" not in app.session_state
