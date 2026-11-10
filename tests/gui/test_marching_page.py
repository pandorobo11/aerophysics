"""Marching boundary-layer adapter, replay and headless page regressions."""

from __future__ import annotations

import math
from types import SimpleNamespace
from typing import Any, cast
from unittest.mock import patch

import numpy as np
import pytest
from streamlit.testing.v1 import AppTest

from aerophysics.boundary_layer import BoundaryLayerRegime
from aerophysics.exceptions import (
    BoundaryLayerConvergenceError,
    BoundaryLayerSeparationError,
)
from aerophysics.gas import AIR
from aerophysics.gui.body_flow_adapter import build_meridian, calculate_body_flow
from aerophysics.gui.config import (
    ConfigurationError,
    dump_configuration,
    load_configuration,
    make_configuration,
)
from aerophysics.gui.csv_inputs import edge_csv_template, parse_edge_csv
from aerophysics.gui.figures import (
    marching_boundary_layer_figures,
    marching_boundary_layer_profile_figure,
)
from aerophysics.gui.marching_adapter import (
    EdgeDefinition,
    MarchingCalculation,
    body_case_edge,
    calculate_marching_boundary_layer,
    profile_rows,
    total_conditions,
    uniform_edge,
)
from aerophysics.gui.tables import display_rows, rows_to_csv
from aerophysics.gui.units import UnitPreferences
from aerophysics.marching_boundary_layer import (
    BoundaryLayerEdge,
    MarchingBoundaryLayerResult,
    marching_boundary_layer,
)

_SCRIPT = """
import streamlit as st
from aerophysics.gui.marching_page import render_marching_boundary_layer
from aerophysics.gui.units import UnitPreferences
render_marching_boundary_layer(
    UnitPreferences(**st.session_state.get("test_units", {}))
)
"""
_PENDING = "pending_marching_boundary_layer_configuration"
_TOTAL_PRESSURE, _TOTAL_TEMPERATURE = total_conditions(2.0, 10000.0, 220.0)


def _app() -> AppTest:
    return AppTest.from_string(_SCRIPT, default_timeout=60)


def _configuration(
    models: dict[str, str] | None = None, **overrides: Any
) -> dict[str, Any]:
    inputs = {
        "mach": 2.0,
        "pressure": 10000.0,
        "temperature": 220.0,
        "length": 1.0,
        "station_count": 50,
        "cone_half_angle": None,
        "total_pressure": _TOTAL_PRESSURE,
        "total_temperature": _TOTAL_TEMPERATURE,
        "edge_distance": None,
        "edge_pressure": None,
        "edge_radius": None,
        "wall_temperature": None,
        "transition_distance": None,
        "prandtl_number": 0.72,
        "turbulent_prandtl_number": 0.9,
        "edge_sa_ratio": 3.0,
    }
    inputs.update(overrides)
    return make_configuration(
        calculator="marching_boundary_layer",
        mode="single",
        inputs_si=inputs,
        models={
            "edge_source": "uniform",
            "geometry": "planar",
            "regime": "laminar",
            "wall": "adiabatic",
            **(models or {}),
        },
        units=UnitPreferences(),
    )


def _separating_edge() -> EdgeDefinition:
    distance = np.linspace(0.0, 1.0, 41)
    return EdgeDefinition(
        distance,
        10000.0 * (1.0 + distance**2),
        None,
        _TOTAL_PRESSURE,
        _TOTAL_TEMPERATURE,
    )


@pytest.fixture(scope="module")
def calculation() -> MarchingCalculation:
    return calculate_marching_boundary_layer(uniform_edge())


def test_total_conditions_are_isentropic() -> None:
    assert _TOTAL_TEMPERATURE == pytest.approx(220.0 * 1.8, rel=1e-14)
    assert _TOTAL_PRESSURE == pytest.approx(10000.0 * 1.8**3.5, rel=1e-14)


def test_default_adapter_matches_public_solver(
    calculation: MarchingCalculation,
) -> None:
    edge = BoundaryLayerEdge.from_pressure(
        np.linspace(0.0, 1.0, 51),
        np.full(51, 10000.0),
        total_pressure=_TOTAL_PRESSURE,
        total_temperature=_TOTAL_TEMPERATURE,
    )
    direct = marching_boundary_layer(edge, regime=BoundaryLayerRegime.LAMINAR)
    rows = calculation.stations.rows
    assert len(rows) == 50 and not calculation.stations.warnings
    assert calculation.separation_distance is None
    for key, field in (
        ("surface_distance", "distance"),
        ("reynolds_number", "reynolds_number"),
        ("boundary_layer_thickness", "boundary_layer_thickness"),
        ("displacement_thickness", "displacement_thickness"),
        ("momentum_thickness", "momentum_thickness"),
        ("local_skin_friction_coefficient", "local_skin_friction_coefficient"),
        ("wall_shear_stress", "wall_shear_stress"),
        ("wall_temperature", "wall_temperature"),
        ("wall_heat_flux", "wall_heat_flux"),
        ("pressure_gradient_parameter", "pressure_gradient_parameter"),
    ):
        np.testing.assert_array_equal([r[key] for r in rows], getattr(direct, field))
    last = rows[-1]
    assert "radius" not in last
    assert last["edge_pressure"] == 10000.0
    assert last["edge_temperature"] == pytest.approx(220.0, rel=1e-13)
    assert last["edge_mach"] == pytest.approx(2.0, rel=1e-13)
    assert last["edge_velocity"] == pytest.approx(
        2.0 * math.sqrt(1.4 * AIR.specific_gas_constant * 220.0), rel=1e-13
    )
    assert last["shape_factor"] == (
        direct.displacement_thickness[-1] / direct.momentum_thickness[-1]
    )


def test_cone_edge_reproduces_the_mangler_factor(
    calculation: MarchingCalculation,
) -> None:
    angle = math.radians(10.0)
    cone = calculate_marching_boundary_layer(
        uniform_edge(geometry="cone", cone_half_angle=angle)
    )
    assert cone.stations.rows[-1]["radius"] == pytest.approx(math.sin(angle))
    assert cone.stations.rows[0]["radius"] == pytest.approx(0.02 * math.sin(angle))
    np.testing.assert_allclose(
        cone.raw.local_skin_friction_coefficient,
        math.sqrt(3.0) * calculation.raw.local_skin_friction_coefficient,
        rtol=1e-10,
    )


@pytest.mark.parametrize(
    "arguments",
    [
        {"mach": 0.0},
        {"pressure": -1.0},
        {"temperature": math.nan},
        {"length": math.inf},
        {"station_count": 1},
        {"station_count": 1001},
        {"station_count": 50.0},
        {"station_count": True},
        {"geometry": "cone"},
        {"geometry": "cone", "cone_half_angle": 0.0},
        {"geometry": "cone", "cone_half_angle": math.pi / 2},
        {"geometry": "wedge"},
    ],
)
def test_uniform_edge_rejects_invalid_inputs(arguments: dict[str, Any]) -> None:
    with pytest.raises(ValueError):
        uniform_edge(**arguments)


def test_wall_temperature_regime_and_warnings_are_forwarded() -> None:
    cooled = calculate_marching_boundary_layer(
        uniform_edge(station_count=4),
        regime="transitional",
        wall_temperature=300.0,
        transition_distance=0.3,
        prandtl_number=0.7,
        turbulent_prandtl_number=0.85,
        edge_sa_ratio=4.0,
    )
    edge = cooled.edge
    direct = marching_boundary_layer(
        edge,
        regime=BoundaryLayerRegime.TRANSITIONAL,
        wall_temperature=300.0,
        transition_distance=0.3,
        prandtl_number=0.7,
        turbulent_prandtl_number=0.85,
        edge_sa_ratio=4.0,
    )
    np.testing.assert_array_equal(cooled.raw.wall_heat_flux, direct.wall_heat_flux)
    np.testing.assert_array_equal(cooled.raw.distance, edge.distance[1:])
    assert cooled.raw.nu_tilde_ratio is not None
    assert np.all(cooled.raw.wall_heat_flux > 0)
    # A thick low-Reynolds layer near a fast pressure drop triggers the sampler.
    distance = np.linspace(0.0, 0.01, 6)
    warned = calculate_marching_boundary_layer(
        EdgeDefinition(
            distance,
            10.0 * np.exp(-200.0 * distance),
            None,
            10.0 * 1.8**3.5,
            _TOTAL_TEMPERATURE,
        )
    )
    assert any("thin-layer" in message for message in warned.stations.warnings)


def test_separation_keeps_the_upstream_stations() -> None:
    separated = calculate_marching_boundary_layer(_separating_edge())
    assert separated.separation_distance == pytest.approx(0.5905, abs=1e-3)
    assert len(separated.stations.rows) == separated.raw.distance.size == 23
    assert separated.raw.distance[-1] < separated.separation_distance
    assert separated.raw.pressure_gradient_parameter[-1] < 0
    assert len(profile_rows(separated, 22)) == separated.raw.eta.size
    empty = BoundaryLayerSeparationError(
        "separated before the first station",
        separation_distance=0.01,
        result=cast(MarchingBoundaryLayerResult, SimpleNamespace(distance=np.empty(0))),
    )
    with (
        patch(
            "aerophysics.gui.marching_adapter.marching_boundary_layer",
            side_effect=empty,
        ),
        pytest.raises(BoundaryLayerSeparationError),
    ):
        calculate_marching_boundary_layer(uniform_edge())


def test_profile_rows(calculation: MarchingCalculation) -> None:
    rows = profile_rows(calculation, 49)
    raw = calculation.raw
    assert len(rows) == raw.eta.size
    assert rows[0]["wall_distance"] == 0.0
    assert rows[0]["velocity"] == pytest.approx(0.0, abs=1e-12)
    assert "nu_tilde_ratio" not in rows[0]
    assert rows[-1]["velocity"] == pytest.approx(
        calculation.stations.rows[-1]["edge_velocity"]
    )
    assert rows[0]["temperature"] == pytest.approx(
        calculation.stations.rows[-1]["wall_temperature"], rel=1e-12
    )
    np.testing.assert_array_equal(
        [row["velocity_ratio"] for row in rows], raw.velocity_ratio[49]
    )
    for station in (-1, 50):
        with pytest.raises(ValueError, match="outside the solved range"):
            profile_rows(calculation, station)
    turbulent = calculate_marching_boundary_layer(
        uniform_edge(station_count=2), regime="turbulent"
    )
    assert turbulent.raw.nu_tilde_ratio is not None
    np.testing.assert_array_equal(
        [row["nu_tilde_ratio"] for row in profile_rows(turbulent, 0)],
        turbulent.raw.nu_tilde_ratio[0],
    )


def test_body_case_edge_matches_boundary_layer_inputs() -> None:
    case = calculate_body_flow(build_meridian()).raw
    definition = body_case_edge(case)
    edge, radius = case.boundary_layer_inputs()
    np.testing.assert_array_equal(definition.distance, edge.distance)
    np.testing.assert_array_equal(definition.pressure, edge.pressure)
    assert definition.radius is not None
    np.testing.assert_array_equal(definition.radius, radius)
    assert definition.total_pressure == case.total_pressure
    assert definition.total_temperature == case.total_temperature
    result = calculate_marching_boundary_layer(definition)
    direct = marching_boundary_layer(
        edge, regime=BoundaryLayerRegime.LAMINAR, wall_radius=radius
    )
    np.testing.assert_array_equal(
        result.raw.local_skin_friction_coefficient,
        direct.local_skin_friction_coefficient,
    )
    np.testing.assert_array_equal(
        [row["radius"] for row in result.stations.rows], radius[1:]
    )


def test_edge_csv_parsing() -> None:
    units = UnitPreferences(length="mm", pressure="kPa")
    planar = parse_edge_csv(edge_csv_template(), UnitPreferences())
    assert planar.radius is None
    np.testing.assert_array_equal(planar.distance, [0, 0.25, 0.5, 0.75, 1])
    np.testing.assert_array_equal(planar.pressure, [10000, 9800, 9600, 9400, 9200])
    body = parse_edge_csv(b"s,p,r,note\n0,10,0,a\n500,9,50,b\n", units)
    np.testing.assert_allclose(body.distance, [0, 0.5])
    np.testing.assert_allclose(body.pressure, [10000, 9000])
    assert body.radius is not None
    np.testing.assert_allclose(body.radius, [0, 0.05])
    for data, message in (
        ("s,q\n0,1\n1,1\n", "CSV must contain columns: s, p"),
        ("s,p\n0.1,1\n1,1\n", "s must start at zero"),
        ("s,p\n0,1\n0,1\n", "s must start at zero"),
        ("s,p\n0,1\n1,0\n", "p must be greater than zero"),
        ("s,p,r\n0,1,-1\n1,1,1\n", "r must be non-negative"),
        ("s,p,r\n0,1,0\n1,1,0\n", "r must be non-negative"),
        ("s,p,r\n0,1,0\n1,1,x\n", "invalid numeric"),
    ):
        with pytest.raises(ValueError, match=message):
            parse_edge_csv(data, UnitPreferences())


def test_settings_roundtrip_and_limits() -> None:
    configuration = _configuration()
    assert load_configuration(dump_configuration(configuration)) == configuration
    embedded = _configuration(
        {"edge_source": "csv"},
        edge_distance=[0, 0.5, 1],
        edge_pressure=[10000, 9500, 9000],
        edge_radius=[0, 0.05, 0.05],
        wall_temperature=300.0,
    )
    assert load_configuration(dump_configuration(embedded)) == embedded
    invalid: tuple[tuple[str, object], ...] = (
        ("mach", 0.0),
        ("station_count", 1),
        ("station_count", 1001),
        ("station_count", 50.5),
        ("cone_half_angle", 0.0),
        ("cone_half_angle", math.pi / 2),
        ("wall_temperature", 0.0),
        ("transition_distance", -0.1),
        ("prandtl_number", 0.0),
        ("edge_sa_ratio", -1.0),
        ("edge_distance", [0.0, "a"]),
    )
    for field, value in invalid:
        broken = _configuration()
        broken["inputs_si"][field] = value
        with pytest.raises(ConfigurationError):
            dump_configuration(broken)
    for field, value in (("regime", "roughness"), ("edge_source", "flight")):
        broken = _configuration()
        broken["models"][field] = value
        with pytest.raises(ConfigurationError):
            dump_configuration(broken)
    sweep = _configuration()
    sweep["mode"] = "sweep"
    with pytest.raises(ConfigurationError):
        dump_configuration(sweep)


def test_figures_tables_and_display_units(calculation: MarchingCalculation) -> None:
    units = UnitPreferences(
        length="mm", pressure="kPa", temperature="°C", heat_flux="kW/m²"
    )
    raw = calculation.raw
    rows = display_rows("marching_boundary_layer", calculation.stations.rows, units)
    assert rows[-1]["表面距離 s [mm]"] == 1000.0
    assert rows[-1]["境界層厚さ δ99 [mm]"] == raw.boundary_layer_thickness[-1] / 1e-3
    assert rows[-1]["壁面せん断応力 τw [kPa]"] == raw.wall_shear_stress[-1] / 1000
    assert rows[-1]["壁温 Tw [°C]"] == raw.wall_temperature[-1] - 273.15
    assert "壁面半径 r_w [mm]" not in rows[0]
    assert "局所摩擦係数 Cf" in rows_to_csv(rows)
    figures = marching_boundary_layer_figures(calculation.stations.rows, units)
    assert list(figures) == ["壁面摩擦・熱", "境界層厚さ", "エッジ条件"]
    wall = figures["壁面摩擦・熱"]
    np.testing.assert_array_equal(wall.data[0].y, raw.local_skin_friction_coefficient)
    np.testing.assert_array_equal(wall.data[2].y, raw.wall_temperature - 273.15)
    np.testing.assert_array_equal(wall.data[0].x, raw.distance / 1e-3)
    assert not wall.layout.shapes
    np.testing.assert_array_equal(
        figures["境界層厚さ"].data[1].y, raw.displacement_thickness / 1e-3
    )
    marked = marching_boundary_layer_figures(
        calculation.stations.rows,
        units,
        adiabatic=False,
        transition_distance=0.3,
        separation_distance=0.9,
    )
    wall = marked["壁面摩擦・熱"]
    np.testing.assert_array_equal(wall.data[2].y, raw.wall_heat_flux / 1000)
    assert sorted(shape.x0 for shape in wall.layout.shapes) == [300.0] * 3 + [900.0] * 3
    profile = profile_rows(calculation, 49)
    table = display_rows("marching_boundary_layer_profile", profile, units)
    assert table[1]["壁面距離 y [mm]"] == raw.wall_distance[49, 1] / 1e-3
    assert len(table[0]) == 6  # No SA column for laminar flow.
    figure = marching_boundary_layer_profile_figure(
        profile, units, boundary_layer_thickness=raw.boundary_layer_thickness[49]
    )
    assert len(figure.data) == 2
    assert 0 < max(figure.data[0].y) <= 2000 * raw.boundary_layer_thickness[49]
    assert max(figure.data[0].x) > 0.99


def test_page_calculate_profile_replay_and_reset() -> None:
    app = _app().run()
    assert not app.exception and not app.error
    assert app.title[0].value == "境界層マーチング解析"
    assert not app.metric
    assert app.number_input(key="march_bl_station_count").value == 50
    app.button(key="march_bl_calculate").click().run()
    assert not app.exception and not app.error and not app.warning
    assert len(app.metric) == 6 and len(app.dataframe) == 2
    assert len(app.get("plotly_chart")) == 4
    assert app.metric[5].label == "断熱壁温 [K]"
    calculation, configuration = app.session_state["march_bl_payload"]
    assert configuration["inputs_si"]["total_pressure"] == _TOTAL_PRESSURE
    assert configuration["inputs_si"]["edge_distance"] is None
    assert app.selectbox(key="march_bl_profile_station").value == 49
    with patch(
        "aerophysics.gui.marching_page.calculate_marching_boundary_layer",
        side_effect=AssertionError("must not recompute"),
    ):
        app.selectbox(key="march_bl_profile_station").set_value(0).run()
        assert not app.exception and not app.warning
        for units in (
            {"length": "mm"},
            {"length": "ft", "pressure": "psf", "temperature": "°R"},
        ):
            app.session_state["test_units"] = units
            app.run()
            assert not app.exception and not app.error and not app.warning
        assert app.number_input(key="march_bl_length").value == pytest.approx(
            1 / 0.3048
        )
        assert app.session_state["march_bl_payload"][0] is calculation
        app.number_input(key="march_bl_mach").set_value(3.0).run()
        assert any("入力が変更されています" in item.value for item in app.warning)
        assert app.session_state["march_bl_payload"][0] is calculation
    replay = _app()
    replay.session_state[_PENDING] = load_configuration(
        dump_configuration(configuration)
    )
    replay.run()
    replay.button(key="march_bl_calculate").click().run()
    assert not replay.exception and not replay.error
    np.testing.assert_array_equal(
        replay.session_state["march_bl_payload"][0].raw.wall_shear_stress,
        calculation.raw.wall_shear_stress,
    )
    app.button(key="march_bl_reset").click().run()
    assert not app.exception and not app.metric
    assert app.number_input(key="march_bl_mach").value == 2


def test_page_cone_isothermal_transitional() -> None:
    app = _app().run()
    app.selectbox(key="march_bl_geometry").set_value("cone").run()
    assert app.number_input(key="march_bl_cone_half_angle").value == pytest.approx(10)
    app.number_input(key="march_bl_station_count").set_value(4).run()
    app.selectbox(key="march_bl_regime").set_value("transitional").run()
    assert app.number_input(key="march_bl_transition_distance").value == 0.5
    app.selectbox(key="march_bl_wall").set_value("isothermal").run()
    assert app.number_input(key="march_bl_wall_temperature").value == 300
    app.button(key="march_bl_calculate").click().run()
    assert not app.exception and not app.error
    calculation, configuration = app.session_state["march_bl_payload"]
    assert configuration["inputs_si"]["cone_half_angle"] == pytest.approx(
        math.radians(10)
    )
    assert configuration["inputs_si"]["transition_distance"] == 0.5
    assert configuration["inputs_si"]["wall_temperature"] == 300
    assert calculation.raw.nu_tilde_ratio is not None
    assert calculation.stations.rows[-1]["radius"] == pytest.approx(
        math.sin(math.radians(10))
    )
    assert app.metric[5].label == "qw [W/m²]"
    # An invalid transition position is the solver's error and clears the result.
    app.number_input(key="march_bl_transition_distance").set_value(1.0).run()
    app.button(key="march_bl_calculate").click().run()
    assert not app.exception
    assert "transition_distance must be less than" in app.error[0].value
    assert "march_bl_payload" not in app.session_state
    app.selectbox(key="march_bl_regime").set_value("laminar").run()
    app.selectbox(key="march_bl_wall").set_value("adiabatic").run()
    app.selectbox(key="march_bl_geometry").set_value("planar").run()
    draft = app.session_state["march_bl_draft"]["inputs_si"]
    assert draft["transition_distance"] is None
    assert draft["wall_temperature"] is None and draft["cone_half_angle"] is None


def _upload(app: AppTest, name: str, data: bytes) -> None:
    uploader: Any = next(
        item for item in app.get("file_uploader") if item.key == "march_bl_edge_upload"
    )
    uploader.upload(name, data, "text/csv").run()


def test_page_csv_upload_keeps_si_values_across_unit_changes() -> None:
    app = _app().run()
    app.selectbox(key="march_bl_edge_source").set_value("csv").run()
    assert app.number_input(key="march_bl_total_pressure").value == pytest.approx(
        _TOTAL_PRESSURE
    )
    app.button(key="march_bl_calculate").click().run()
    assert "エッジ圧力CSVを選択してください" in app.error[0].value
    _upload(app, "edge.csv", edge_csv_template().encode("utf-8"))
    app.button(key="march_bl_calculate").click().run()
    assert not app.exception and not app.error and not app.warning
    payload = app.session_state["march_bl_payload"]
    assert payload[1]["inputs_si"]["edge_pressure"] == [10000, 9800, 9600, 9400, 9200]
    assert payload[1]["inputs_si"]["edge_radius"] is None
    assert len(payload[0].stations.rows) == 4
    assert payload[0].stations.rows[-1]["pressure_gradient_parameter"] > 0
    with patch(
        "aerophysics.gui.marching_page.calculate_marching_boundary_layer",
        side_effect=AssertionError("must not recompute"),
    ):
        app.session_state["test_units"] = {"length": "mm", "pressure": "kPa"}
        app.run()
        assert not app.exception and not app.error and not app.warning
        assert app.session_state["march_bl_payload"] is payload
    # A replacement file is read in the units active when it arrives.
    _upload(app, "body.csv", b"s,p,r\n0,10,0\n500,9.5,50\n1000,9,50\n")
    app.button(key="march_bl_calculate").click().run()
    assert not app.exception and not app.error
    inputs = app.session_state["march_bl_payload"][1]["inputs_si"]
    assert inputs["edge_distance"] == pytest.approx([0, 0.5, 1])
    assert inputs["edge_pressure"] == pytest.approx([10000, 9500, 9000])
    assert inputs["edge_radius"] == pytest.approx([0, 0.05, 0.05])
    _upload(app, "bad.csv", b"s,p\n0,1\n1,-1\n")
    app.button(key="march_bl_calculate").click().run()
    assert not app.exception
    assert "p must be greater than zero" in app.error[0].value
    assert "march_bl_payload" not in app.session_state
    # Arrays of one source are dropped when another source is selected.
    replay = _app()
    replay.session_state[_PENDING] = load_configuration(
        dump_configuration(
            _configuration(
                {"edge_source": "csv"},
                edge_distance=[0, 0.5, 1],
                edge_pressure=[10000, 9500, 9000],
            )
        )
    )
    replay.run()
    replay.run()  # Embedded arrays survive another rerun without the uploader.
    replay.button(key="march_bl_calculate").click().run()
    assert not replay.exception and not replay.error
    assert len(replay.session_state["march_bl_payload"][0].stations.rows) == 2
    replay.selectbox(key="march_bl_edge_source").set_value("body").run()
    assert replay.session_state["march_bl_draft"]["inputs_si"]["edge_distance"] is None
    assert any("保存" in item.value for item in replay.info)


def test_page_body_case_and_embedded_edge() -> None:
    case = calculate_body_flow(build_meridian()).raw
    other = calculate_body_flow(build_meridian(nose_length=0.4)).raw
    app = _app().run()
    app.selectbox(key="march_bl_edge_source").set_value("body").run()
    app.button(key="march_bl_calculate").click().run()
    assert not app.exception and "保存" in app.error[0].value
    app.session_state["current_axisymmetric_body_case"] = case
    app.run()
    assert any("軸対称" in item.value for item in app.info)
    app.button(key="march_bl_calculate").click().run()
    assert not app.exception and not app.error and not app.warning
    calculation, configuration = app.session_state["march_bl_payload"]
    edge, radius = case.boundary_layer_inputs()
    assert configuration["inputs_si"]["edge_distance"] == edge.distance.tolist()
    assert configuration["inputs_si"]["edge_radius"] == radius.tolist()
    assert configuration["inputs_si"]["total_pressure"] == case.total_pressure
    assert len(calculation.stations.rows) == edge.distance.size - 1
    # A newly saved case changes the inputs and flags the retained result.
    app.session_state["current_axisymmetric_body_case"] = other
    app.run()
    assert any("入力が変更されています" in item.value for item in app.warning)
    # A loaded settings file keeps its embedded edge until the user switches.
    replay = _app()
    replay.session_state["current_axisymmetric_body_case"] = other
    replay.session_state[_PENDING] = load_configuration(
        dump_configuration(configuration)
    )
    replay.run()
    draft = replay.session_state["march_bl_draft"]["inputs_si"]
    assert draft["edge_distance"] == edge.distance.tolist()
    replay.button(key="march_bl_use_case").click().run()
    assert not replay.exception
    draft = replay.session_state["march_bl_draft"]["inputs_si"]
    assert (
        draft["edge_distance"] == other.boundary_layer_inputs().edge.distance.tolist()
    )
    assert draft["total_pressure"] == other.total_pressure
    alone = _app()
    alone.session_state[_PENDING] = load_configuration(
        dump_configuration(configuration)
    )
    alone.run()
    assert not [button for button in alone.button if button.key == "march_bl_use_case"]
    alone.button(key="march_bl_calculate").click().run()
    assert not alone.exception and not alone.error
    np.testing.assert_array_equal(
        alone.session_state["march_bl_payload"][0].raw.wall_shear_stress,
        calculation.raw.wall_shear_stress,
    )


def test_page_separation_shows_upstream_result() -> None:
    edge = _separating_edge()
    app = _app()
    app.session_state[_PENDING] = load_configuration(
        dump_configuration(
            _configuration(
                {"edge_source": "csv"},
                edge_distance=edge.distance.tolist(),
                edge_pressure=edge.pressure.tolist(),
            )
        )
    )
    app.run()
    app.button(key="march_bl_calculate").click().run()
    assert not app.exception
    assert "剥離しました" in app.error[0].value and "0.59" in app.error[0].value
    assert len(app.metric) == 6
    assert len(app.session_state["march_bl_payload"][0].stations.rows) == 23
    assert app.selectbox(key="march_bl_profile_station").value == 22


@pytest.mark.parametrize(
    "error",
    [
        BoundaryLayerConvergenceError("march failed at s=0.42"),
        ValueError("edge is not isentropic"),
    ],
)
def test_page_solver_errors_clear_stale_payload(
    error: Exception, calculation: MarchingCalculation
) -> None:
    app = _app().run()
    app.session_state["march_bl_payload"] = (calculation, _configuration())
    with patch(
        "aerophysics.gui.marching_page.calculate_marching_boundary_layer",
        side_effect=error,
    ):
        app.button(key="march_bl_calculate").click().run()
    assert not app.exception
    assert app.error[0].value == str(error)
    assert "march_bl_payload" not in app.session_state
    assert not app.dataframe


def test_page_separation_without_stations_and_runtime_errors(
    calculation: MarchingCalculation,
) -> None:
    app = _app().run()
    with patch(
        "aerophysics.gui.marching_page.calculate_marching_boundary_layer",
        side_effect=BoundaryLayerSeparationError(
            "the wall shear collapses", separation_distance=0.01, result=calculation.raw
        ),
    ):
        app.button(key="march_bl_calculate").click().run()
    assert not app.exception
    assert "最初の出力点より上流で剥離" in app.error[0].value
    assert "march_bl_payload" not in app.session_state
    with patch(
        "aerophysics.gui.marching_page.calculate_marching_boundary_layer",
        side_effect=RuntimeError("unexpected failure"),
    ):
        app.button(key="march_bl_calculate").click().run()
    assert app.exception[0].message == "unexpected failure"
    assert not app.error


@pytest.mark.parametrize(
    ("distance", "pressure", "radius"),
    [([], [], None), ([0.0], [1.0], None), ([0.0, 1.0], [1e4, 1e4], [0.0])],
)
def test_page_malformed_embedded_edge_shows_error(
    distance: list[float], pressure: list[float], radius: list[float] | None
) -> None:
    app = _app()
    app.session_state[_PENDING] = load_configuration(
        dump_configuration(
            _configuration(
                {"edge_source": "csv"},
                edge_distance=distance,
                edge_pressure=pressure,
                edge_radius=radius,
            )
        )
    )
    app.run()
    assert not app.exception
    app.button(key="march_bl_calculate").click().run()
    assert not app.exception and app.error
    assert "march_bl_payload" not in app.session_state
