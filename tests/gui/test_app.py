"""Headless Streamlit GUI tests."""

from pathlib import Path
from typing import Any, cast

import pytest
from streamlit.testing.v1 import AppTest

from aerophysics.gui.adapters import FlightCase
from aerophysics.gui.advanced_adapters import BoundaryLayerCase, BoundaryProfileCase
from aerophysics.gui.config import make_configuration
from aerophysics.gui.units import UnitPreferences

APP = Path("src/aerophysics/gui/app.py")


def test_main_app_default_page_and_calculation() -> None:
    app = AppTest.from_file(APP, default_timeout=15).run()
    assert not app.exception
    assert app.title[0].value == "大気・飛行条件"
    app.button(key="FormSubmitter:flight_form-計算").click().run()
    assert not app.exception
    assert not app.error
    assert len(app.metric) == 4
    assert len(app.dataframe) == 1
    assert len(app.get("plotly_chart")) == 3
    assert len(app.download_button) == 2
    app.button(key="flight_save_case").click().run()
    assert isinstance(app.session_state["current_flight_case"], FlightCase)


def test_flight_page_sweep_mode() -> None:
    app = AppTest.from_file(APP, default_timeout=15).run()
    app.radio(key="flight_mode").set_value("1変数スイープ").run()
    app.button(key="FormSubmitter:flight_form-計算").click().run()
    assert not app.exception
    assert not app.error
    assert len(app.dataframe[0].value) == 101
    assert len(app.get("plotly_chart")) == 3


def test_shock_page_calculation() -> None:
    script = """
from aerophysics.gui.pages import render_shock
from aerophysics.gui.units import UnitPreferences
render_shock(UnitPreferences())
"""
    app = AppTest.from_string(script, default_timeout=15).run()
    app.button(key="FormSubmitter:shock_form-計算").click().run()
    assert not app.exception
    assert not app.error
    assert app.metric[0].label == "衝撃波角 β"
    assert len(app.get("plotly_chart")) == 3


def test_shock_page_sweep_marks_non_attached_rows() -> None:
    script = """
from aerophysics.gui.pages import render_shock
from aerophysics.gui.units import UnitPreferences
render_shock(UnitPreferences())
"""
    app = AppTest.from_string(script, default_timeout=15).run()
    app.radio(key="shock_mode").set_value("1変数スイープ").run()
    app.button(key="FormSubmitter:shock_form-計算").click().run()
    assert not app.exception
    assert not app.error
    assert any(
        status == "no_attached_shock"
        for status in app.dataframe[0].value["status"].tolist()
    )
    assert app.warning


def test_conical_shock_page_calculation_and_sweep() -> None:
    script = """
from aerophysics.gui.pages import render_conical_shock
from aerophysics.gui.units import UnitPreferences
render_conical_shock(UnitPreferences())
"""
    app = AppTest.from_string(script, default_timeout=30).run()
    app.button(key="FormSubmitter:cone_shock_form-計算").click().run()
    assert not app.exception
    assert not app.error
    assert app.metric[0].label == "衝撃波角 β"
    assert len(app.get("plotly_chart")) == 3

    app.radio(key="cone_shock_mode").set_value("1変数スイープ").run()
    app.number_input(key="cone_shock_sweep_points").set_value(3).run()
    app.button(key="FormSubmitter:cone_shock_form-計算").click().run()
    assert not app.exception
    assert any(
        status == "no_attached_shock"
        for status in app.dataframe[0].value["status"].tolist()
    )
    assert app.warning


def test_boundary_layer_page_calculation_and_case_source() -> None:
    script = """
from aerophysics.gui.pages import render_boundary_layer
from aerophysics.gui.units import UnitPreferences
render_boundary_layer(UnitPreferences())
"""
    app = AppTest.from_string(script, default_timeout=15)
    app.session_state["current_flight_case"] = FlightCase(
        geometric_altitude=10000.0,
        mach=0.8,
        velocity=240.0,
        density=0.41,
        dynamic_viscosity=1.46e-5,
        temperature=223.0,
    )
    app.run()
    assert app.radio(key="boundary_source").options == ["手入力", "現在の飛行ケース"]
    app.radio(key="boundary_source").set_value("現在の飛行ケース").run()
    app.button(key="FormSubmitter:boundary_form-計算").click().run()
    assert not app.exception
    assert not app.error
    assert len(app.metric) == 4
    assert len(app.get("plotly_chart")) == 3
    app.button(key="boundary_save_case").click().run()
    assert isinstance(
        app.session_state["current_boundary_layer_case"], BoundaryLayerCase
    )


def test_additional_compressible_flow_pages() -> None:
    pages = (
        ("render_isentropic", "isentropic_form", "等エントロピー流れ", 2),
        ("render_normal_shock", "normal_form", "垂直衝撃波", 2),
        ("render_expansion", "expansion_form", "Prandtl\u2013Meyer膨張", 3),
    )
    for function, form, title, plot_count in pages:
        script = f"""
from aerophysics.gui.flow_pages import {function}
from aerophysics.gui.units import UnitPreferences
{function}(UnitPreferences())
"""
        app = AppTest.from_string(script, default_timeout=15).run()
        assert app.title[0].value == title
        app.button(key=f"FormSubmitter:{form}-計算").click().run()
        assert not app.exception
        assert not app.error
        assert len(app.metric) == 4
        assert len(app.get("plotly_chart")) == plot_count
        assert len(app.download_button) == 2


def test_isentropic_page_supports_thermally_perfect_air() -> None:
    script = """
from aerophysics.gui.flow_pages import render_isentropic
from aerophysics.gui.units import UnitPreferences
render_isentropic(UnitPreferences())
"""
    app = AppTest.from_string(script, default_timeout=30).run()
    app.selectbox(key="isentropic_gas_model").set_value("NASA9").run()
    app.number_input(key="isentropic_input").set_value(2.0).run()
    app.number_input(key="isentropic_total_temperature").set_value(1000.0).run()
    app.button(key="FormSubmitter:isentropic_form-計算").click().run()
    assert not app.exception
    assert not app.error
    assert app.dataframe[0].value["気体モデル"].tolist() == ["NASA9"]
    assert app.dataframe[0].value["静温 T [K]"].iloc[0] == pytest.approx(
        580.6729799
    )


def test_isentropic_page_loads_legacy_configuration_as_air() -> None:
    script = """
from aerophysics.gui.flow_pages import render_isentropic
from aerophysics.gui.units import UnitPreferences
render_isentropic(UnitPreferences())
"""
    configuration = make_configuration(
        calculator="isentropic",
        mode="single",
        inputs_si={
            "input_value": 2.0,
            "total_pressure": None,
            "total_temperature": None,
        },
        models={
            "input_basis": "mach",
            "branch": "subsonic",
            "with_mass_flux": False,
        },
        units=UnitPreferences(),
    )
    app = AppTest.from_string(script, default_timeout=30)
    app.session_state["pending_isentropic_configuration"] = configuration
    app.run()
    assert app.selectbox(key="isentropic_gas_model").value == "AIR"
    assert app.number_input(key="isentropic_total_temperature").value == 300.0
    app.button(key="FormSubmitter:isentropic_form-計算").click().run()
    assert not app.exception
    assert not app.error


def test_expansion_sweep_marks_limit_rows() -> None:
    script = """
from aerophysics.gui.flow_pages import render_expansion
from aerophysics.gui.units import UnitPreferences
render_expansion(UnitPreferences())
"""
    app = AppTest.from_string(script, default_timeout=15).run()
    app.radio(key="expansion_mode").set_value("1変数スイープ").run()
    app.number_input(key="expansion_sweep_stop").set_value(130.0).run()
    app.button(key="FormSubmitter:expansion_form-計算").click().run()
    assert not app.exception
    assert any(
        status == "limit_exceeded"
        for status in app.dataframe[0].value["status"].tolist()
    )
    assert app.warning


def test_advanced_analysis_pages_default_calculations() -> None:
    pages = (
        ("render_boundary_layer_profile", "profile_form", 4),
        ("render_protrusion_drag", "protrusion_form", 1),
        ("render_thermochemistry", "thermo_form", 4),
    )
    for function, form, plot_count in pages:
        script = f"""
from aerophysics.gui.analysis_pages import {function}
from aerophysics.gui.units import UnitPreferences
{function}(UnitPreferences())
"""
        app = AppTest.from_string(script, default_timeout=30).run()
        app.button(key=f"FormSubmitter:{form}-計算").click().run()
        assert not app.exception
        assert not app.error
        assert len(app.metric) == 4
        assert len(app.dataframe) == 1
        assert len(app.get("plotly_chart")) == plot_count


def test_profile_to_protrusion_session_transfer() -> None:
    profile_script = """
from aerophysics.gui.analysis_pages import render_boundary_layer_profile
from aerophysics.gui.advanced_adapters import BoundaryLayerCase
from aerophysics.gui.units import UnitPreferences
import streamlit as st
st.session_state["current_boundary_layer_case"] = BoundaryLayerCase(
    300.0, 1.0, 300.0, 0.05, 85.0
)
render_boundary_layer_profile(UnitPreferences())
"""
    app = AppTest.from_string(profile_script, default_timeout=30).run()
    app.radio(key="profile_source").set_value("現在の乱流平板境界層ケース").run()
    app.button(key="FormSubmitter:profile_form-計算").click().run()
    app.button(key="profile_save_case").click().run()
    saved = app.session_state["current_boundary_profile"]
    assert isinstance(saved, BoundaryProfileCase)

    protrusion_script = """
from aerophysics.gui.analysis_pages import render_protrusion_drag
from aerophysics.gui.units import UnitPreferences
render_protrusion_drag(UnitPreferences())
"""
    drag = AppTest.from_string(protrusion_script, default_timeout=30)
    drag.session_state["current_boundary_profile"] = saved
    drag.run()
    assert "現在の境界層プロファイル" in drag.radio(key="protrusion_source").options
    drag.radio(key="protrusion_source").set_value("現在の境界層プロファイル").run()
    drag.button(key="FormSubmitter:protrusion_form-計算").click().run()
    assert not drag.exception
    assert not drag.error
    assert drag.dataframe[0].value["プロファイル"].iloc[0] == "provided"


def test_protrusion_csv_uploads() -> None:
    script = """
from aerophysics.gui.analysis_pages import render_protrusion_drag
from aerophysics.gui.units import UnitPreferences
render_protrusion_drag(UnitPreferences())
"""
    app = AppTest.from_string(script, default_timeout=30).run()
    app.radio(key="protrusion_source").set_value("プロファイルCSV").run()
    profile = cast(
        Any,
        next(
            item
            for item in app.get("file_uploader")
            if item.key == "protrusion_profile_upload"
        ),
    )
    profile.upload(
        "profile.csv",
        b"wall_distance,velocity,density\n0,0,1\n0.05,300,1\n",
        "text/csv",
    ).run()
    app.selectbox(key="protrusion_shape").set_value("形状CSV").run()
    shape = cast(
        Any,
        next(
            item
            for item in app.get("file_uploader")
            if item.key == "protrusion_shape_upload"
        ),
    )
    shape.upload(
        "shape.csv",
        b"height,width\n0,0.005\n0.02,0\n",
        "text/csv",
    ).run()
    app.button(key="FormSubmitter:protrusion_form-計算").click().run()
    assert not app.exception
    assert not app.error
    assert app.dataframe[0].value["投影形状"].iloc[0] == "csv"


def test_protrusion_import_uses_embedded_profile() -> None:
    script = """
from aerophysics.gui.analysis_pages import render_protrusion_drag
from aerophysics.gui.units import UnitPreferences
render_protrusion_drag(UnitPreferences())
"""
    configuration = make_configuration(
        calculator="protrusion_drag",
        mode="single",
        inputs_si={
            "drag_coefficient": 1.0,
            "height": 0.01,
            "base_width": 0.005,
            "edge_velocity": 100.0,
            "edge_density": 1.0,
            "boundary_layer_thickness": 0.02,
            "mach": None,
            "edge_temperature": None,
            "wall_temperature": None,
            "profile_height": [0.0, 0.02],
            "profile_velocity": [0.0, 100.0],
            "profile_density": [1.0, 1.0],
            "shape_height": None,
            "shape_width": None,
        },
        models={
            "profile_source": "saved",
            "shape": "rectangle",
            "compressible": False,
        },
        units=UnitPreferences(),
    )
    app = AppTest.from_string(script, default_timeout=30)
    app.session_state["pending_protrusion_drag_configuration"] = configuration
    app.run()
    assert app.radio(key="protrusion_source").value == "csv"
    app.button(key="FormSubmitter:protrusion_form-計算").click().run()
    assert not app.exception
    assert not app.error
    assert app.dataframe[0].value["プロファイル"].iloc[0] == "provided"
