"""GUI state/output controls, unit conversion and settings replay."""

import pytest
from streamlit.testing.v1 import AppTest

from aerophysics.gui.config import dump_configuration, load_configuration


@pytest.mark.parametrize(
    "module,function,prefix,calculator",
    [
        ("flow_pages", "render_normal_shock", "normal", "normal_shock"),
        ("pages", "render_shock", "shock", "oblique_shock"),
        ("pages", "render_conical_shock", "cone_shock", "conical_shock"),
        ("flow_pages", "render_isentropic", "isentropic", "isentropic"),
    ],
)
def test_flow_output_settings_units_and_sweep_replay(
    module: str, function: str, prefix: str, calculator: str
) -> None:
    script = f"""
from aerophysics.gui.{module} import {function}
from aerophysics.gui.components import render_unit_sidebar
{function}(render_unit_sidebar())
"""
    app = AppTest.from_string(script, default_timeout=30).run()
    app.selectbox(key=f"{prefix}_gas_model").set_value("NASA9").run()
    app.checkbox(key=f"{prefix}_with_heat_capacities").check().run()
    if prefix == "isentropic":
        app.checkbox(key="isentropic_with_flux").check().run()
        app.number_input(key="isentropic_total_temperature").set_value(1000.0).run()
        pressure_key = "isentropic_total_pressure"
    else:
        app.checkbox(key=f"{prefix}_with_pressure").check().run()
        app.number_input(key=f"{prefix}_mach").set_value(3.0).run()
        app.number_input(key=f"{prefix}_upstream_temperature").set_value(500.0).run()
        pressure_key = f"{prefix}_upstream_pressure"
    app.number_input(key=pressure_key).set_value(30_000.0).run()
    app.selectbox(key=f"{prefix}_viscosity_model").set_value("Keyes").run()
    app.checkbox(key=f"{prefix}_with_length").check().run()
    app.number_input(key=f"{prefix}_characteristic_length").set_value(0.2).run()
    app.selectbox(key="unit_length").set_value("mm").run()
    app.selectbox(key="unit_pressure").set_value("kPa").run()
    app.selectbox(key="unit_temperature").set_value("°C").run()
    assert app.number_input(
        key=f"{prefix}_characteristic_length"
    ).value == pytest.approx(200)
    assert app.number_input(key=pressure_key).value == pytest.approx(30)
    assert not app.exception
    assert f"{prefix}_payload" not in app.session_state
    app.button(key=f"FormSubmitter:{prefix}_form-計算").click().run()
    assert not app.exception and not app.error
    table = app.dataframe[0].value
    assert any("単位Re" in heading and "1/mm" in heading for heading in table.columns)
    assert any("定圧比熱" in heading for heading in table.columns)
    assert table["粘性モデル"].tolist() == ["Keyes"]
    app.radio(key=f"{prefix}_mode").set_value("sweep").run()
    app.number_input(key=f"{prefix}_sweep_start").set_value(
        1.0 if prefix == "isentropic" else 3.0 if prefix == "normal" else 5.0
    ).run()
    app.number_input(key=f"{prefix}_sweep_stop").set_value(
        2.0 if prefix == "isentropic" else 5.4 if prefix == "normal" else 10.0
    ).run()
    app.number_input(key=f"{prefix}_sweep_points").set_value(2).run()
    app.button(key=f"FormSubmitter:{prefix}_form-計算").click().run()
    assert not app.exception and not app.error
    result, configuration = app.session_state[f"{prefix}_payload"]
    assert len(result.rows) == 2
    assert configuration["inputs_si"]["characteristic_length"] == pytest.approx(0.2)
    assert configuration["inputs_si"][
        "total_pressure" if prefix == "isentropic" else "upstream_pressure"
    ] == pytest.approx(30_000)
    assert configuration["models"]["viscosity_model"] == "Keyes"
    assert configuration["models"]["with_heat_capacities"] is True
    replay = AppTest.from_string(script, default_timeout=30)
    replay.session_state[f"pending_{calculator}_configuration"] = load_configuration(
        dump_configuration(configuration)
    )
    # Configuration import sets the sidebar units before queuing this payload.
    for name, unit in configuration["display_units"].items():
        replay.session_state[f"unit_{name}"] = unit
    replay.session_state["_gui_display_units"] = configuration["display_units"]
    replay.run()
    assert replay.checkbox(key=f"{prefix}_with_heat_capacities").value
    assert replay.checkbox(key=f"{prefix}_with_length").value
    assert replay.selectbox(key=f"{prefix}_viscosity_model").value == "Keyes"
    replay.button(key=f"FormSubmitter:{prefix}_form-計算").click().run()
    assert not replay.exception and not replay.error
    assert replay.session_state[f"{prefix}_payload"][0].rows == result.rows
    assert replay.session_state[f"{prefix}_payload"][1] == configuration


def test_air_temperature_and_transport_warning_can_be_disabled() -> None:
    script = """
from aerophysics.gui.flow_pages import render_normal_shock
from aerophysics.gui.units import UnitPreferences
render_normal_shock(UnitPreferences())
"""
    app = AppTest.from_string(script).run()
    app.checkbox(key="normal_with_temperature").check().run()
    app.button(key="FormSubmitter:normal_form-計算").click().run()
    assert "速度 V₁ [m/s]" in app.dataframe[0].value.columns
    assert "単位Re₁ [1/m]" not in app.dataframe[0].value.columns
    app.checkbox(key="normal_with_pressure").check().run()
    app.selectbox(key="normal_viscosity_model").set_value("Blottner/Wilke").run()
    app.button(key="FormSubmitter:normal_form-計算").click().run()
    assert not app.exception and not app.error
    assert app.warning
    assert app.session_state["normal_payload"][0].rows[0]["status"] == "ok"
    app.checkbox(key="normal_with_pressure").uncheck().run()
    app.checkbox(key="normal_with_temperature").uncheck().run()
    app.button(key="FormSubmitter:normal_form-計算").click().run()
    assert not app.warning
    assert "速度 V₁ [m/s]" not in app.dataframe[0].value.columns
    assert "単位Re₁ [1/m]" not in app.dataframe[0].value.columns
