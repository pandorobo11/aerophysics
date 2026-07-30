"""Headless Streamlit GUI tests."""

from pathlib import Path

from streamlit.testing.v1 import AppTest

from aerophysics.gui.adapters import FlightCase

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
