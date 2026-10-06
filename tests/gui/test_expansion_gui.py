"""Thermal expansion GUI units, replay, ranges, plots, and failure isolation."""

import numpy as np
import pytest
from streamlit.testing.v1 import AppTest

from aerophysics.exceptions import ExpansionConvergenceError, ModelRangeError
from aerophysics.expansion import prandtl_meyer_expansion
from aerophysics.gui.adapters import expansion_condition, expansion_sweep
from aerophysics.gui.config import (
    ConfigurationError,
    dump_configuration,
    load_configuration,
    make_configuration,
)
from aerophysics.gui.figures import expansion_figures
from aerophysics.gui.tables import display_rows, rows_to_csv
from aerophysics.gui.units import UnitPreferences


@pytest.mark.parametrize("gas_model", ["NASA7", "NASA9", "HARMONIC_OSCILLATOR"])
def test_thermal_page_units_and_settings_replay(gas_model: str) -> None:
    script = """
from aerophysics.gui.flow_pages import render_expansion
from aerophysics.gui.components import render_unit_sidebar
render_expansion(render_unit_sidebar())
"""
    app = AppTest.from_string(script, default_timeout=30).run()
    app.selectbox(key="expansion_gas_model").set_value(gas_model).run()
    app.selectbox(key="unit_temperature").set_value("°C").run()
    assert app.number_input(key="expansion_temperature").value == pytest.approx(726.85)
    app.button(key="FormSubmitter:expansion_form-計算").click().run()
    assert not app.exception and not app.error
    table = app.dataframe[0].value
    assert table["気体モデル"].tolist() == [gas_model]
    assert table["下流静温 T₂ [°C]"].iloc[0] < 726.85
    app.radio(key="expansion_mode").set_value("1変数スイープ").run()
    app.number_input(key="expansion_sweep_stop").set_value(80.0).run()
    app.number_input(key="expansion_sweep_points").set_value(3).run()
    app.button(key="FormSubmitter:expansion_form-計算").click().run()
    assert not app.exception and not app.error and app.warning
    result, configuration = app.session_state["expansion_payload"]
    assert len(result.rows) == 3
    assert (
        result.rows[0]["status"] == "ok" and result.rows[-1]["status"] == "out_of_range"
    )
    assert configuration["models"]["gas_model"] == gas_model
    assert configuration["inputs_si"]["upstream_temperature"] == 1000.0
    replay = AppTest.from_string(
        """
from aerophysics.gui.flow_pages import render_expansion
from aerophysics.gui.units import UnitPreferences
render_expansion(UnitPreferences(temperature='°C'))
""",
        default_timeout=30,
    )
    replay.session_state["pending_expansion_configuration"] = load_configuration(
        dump_configuration(configuration)
    )
    replay.run()
    assert replay.selectbox(key="expansion_gas_model").value == gas_model
    assert replay.number_input(key="expansion_temperature").value == pytest.approx(
        726.85
    )
    replay.button(key="FormSubmitter:expansion_form-計算").click().run()
    assert not replay.exception and not replay.error
    assert replay.session_state["expansion_payload"][0].rows == result.rows
    assert replay.session_state["expansion_payload"][1] == configuration


def test_legacy_configuration_range_errors_and_missing_nu() -> None:
    script = """
from aerophysics.gui.flow_pages import render_expansion
from aerophysics.gui.units import UnitPreferences
render_expansion(UnitPreferences())
"""
    app = AppTest.from_string(script, default_timeout=30)
    app.session_state["pending_expansion_configuration"] = make_configuration(
        calculator="expansion",
        mode="single",
        inputs_si={"upstream_mach": 3.0, "turn_angle": float(np.deg2rad(10.0))},
        models={},
        units=UnitPreferences(),
    )
    app.run()
    assert app.selectbox(key="expansion_gas_model").value == "AIR"
    assert "expansion_temperature" not in {widget.key for widget in app.number_input}
    app.selectbox(key="expansion_gas_model").set_value("HARMONIC_OSCILLATOR").run()
    app.button(key="FormSubmitter:expansion_form-計算").click().run()
    assert not app.exception and not app.error and app.info
    result = app.session_state["expansion_payload"][0]
    assert result.rows[0]["upstream_prandtl_meyer_angle"] is None
    assert result.rows[0]["downstream_mach"] == pytest.approx(3.4950150316)
    app.number_input(key="expansion_temperature").set_value(399.0).run()
    app.button(key="FormSubmitter:expansion_form-計算").click().run()
    assert not app.exception and app.error
    assert "expansion_payload" not in app.session_state


def test_thermal_sweep_keeps_range_rows_and_exports_temperature_units() -> None:
    result = expansion_sweep(
        fixed_mach=2.0,
        fixed_turn_angle=0.1,
        sweep_field="turn_angle",
        start=0.0,
        stop=float(np.deg2rad(30.0)),
        points=3,
        gas_model="HARMONIC_OSCILLATOR",
        upstream_temperature=500.0,
    )
    assert [row["status"] for row in result.rows] == [
        "ok",
        "out_of_range",
        "out_of_range",
    ]
    assert result.rows[0]["maximum_turn_angle"] is None
    assert result.rows[0]["temperature_limited_turn_angle"] is not None
    units = UnitPreferences(temperature="°C", angle="rad")
    table = display_rows("expansion", result.rows, units)
    assert table[0]["上流静温 T₁ [°C]"] == pytest.approx(226.85)
    csv = rows_to_csv(table)
    assert "下流静温 T₂ [°C]" in csv and "out_of_range" in csv
    plots = expansion_figures(result.rows, units, sweep_field="turn_angle")
    assert np.isnan(plots["Mach数"].data[1].y[1])
    assert "温度範囲内の最大膨張角" in [trace.name for trace in plots["角度"].data]
    missing = expansion_condition(
        upstream_mach=3.0,
        turn_angle=0.1,
        gas_model="HARMONIC_OSCILLATOR",
        upstream_temperature=1000.0,
    )
    missing_plots = expansion_figures(missing.rows, units, sweep_field="mach")
    assert [trace.name for trace in missing_plots["角度"].data] == [
        "温度範囲内の最大膨張角"
    ]
    mach_sweep = expansion_sweep(
        fixed_mach=2.0,
        fixed_turn_angle=0.1,
        sweep_field="mach",
        start=1.0,
        stop=3.0,
        points=3,
        gas_model="NASA9",
        upstream_temperature=1000.0,
    )
    assert all(row["status"] == "ok" for row in mach_sweep.rows)
    with pytest.raises(ValueError, match="gas_model"):
        expansion_condition(upstream_mach=2.0, turn_angle=0.1, gas_model="unknown")
    with pytest.raises(ModelRangeError):
        expansion_condition(
            upstream_mach=2.0,
            turn_angle=0.1,
            gas_model="NASA7",
            upstream_temperature=100.0,
        )


def test_numerical_failure_is_displayed_and_sweep_continues(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import aerophysics.gui.adapters as adapters

    original = prandtl_meyer_expansion

    def fail_one(upstream_mach: float, *args: object, **kwargs: object) -> object:
        if upstream_mach == 2.0:
            raise ExpansionConvergenceError("could not resolve expansion")
        return original(upstream_mach, *args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(adapters, "prandtl_meyer_expansion", fail_one)
    sweep = expansion_sweep(
        fixed_mach=2.0,
        fixed_turn_angle=0.1,
        sweep_field="mach",
        start=2.0,
        stop=3.0,
        points=2,
        gas_model="NASA9",
        upstream_temperature=1000.0,
    )
    assert [row["status"] for row in sweep.rows] == ["error", "ok"]
    app = AppTest.from_string(
        """
from aerophysics.gui.flow_pages import render_expansion
from aerophysics.gui.units import UnitPreferences
render_expansion(UnitPreferences())
""",
        default_timeout=30,
    ).run()
    app.selectbox(key="expansion_gas_model").set_value("NASA9").run()
    app.button(key="FormSubmitter:expansion_form-計算").click().run()
    assert not app.exception and app.error
    assert "could not resolve expansion" in app.error[0].value


@pytest.mark.parametrize("temperature", [0.0, "500", np.nan])
def test_thermal_configuration_validates_fields(temperature: object) -> None:
    with pytest.raises(ConfigurationError):
        make_configuration(
            calculator="expansion",
            mode="single",
            inputs_si={
                "upstream_mach": 2.0,
                "turn_angle": 0.1,
                "upstream_temperature": temperature,
            },
            models={"gas_model": "NASA9"},
            units=UnitPreferences(),
        )
