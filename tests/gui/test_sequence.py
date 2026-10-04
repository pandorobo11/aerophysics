"""Sequence GUI replay, editing and schema validation."""

import copy
import math
from typing import Any, cast

import pytest
from streamlit.testing.v1 import AppTest
from streamlit.testing.v1.element_tree import Button

from aerophysics.gui.config import (
    ConfigurationError,
    dump_configuration,
    load_configuration,
    make_configuration,
)
from aerophysics.gui.sequence_adapter import (
    calculate_sequence,
    decode_step,
    sequence_rows,
    validate_sequence_payload,
)
from aerophysics.gui.sequence_page import _default_step, _display_rows
from aerophysics.gui.units import UnitPreferences

SCRIPT = """
from aerophysics.gui.sequence_page import render_sequence
from aerophysics.gui.components import render_unit_sidebar
render_sequence(render_unit_sidebar())
"""


def configuration(basis: str = "static_mach", gas: str = "AIR") -> dict[str, Any]:
    values = {
        "static_mach": {"pressure": 10000, "temperature": 500, "mach": 2},
        "static_velocity": {"pressure": 10000, "temperature": 500, "velocity": 900},
        "total": {"pressure": 100000, "temperature": 1000, "mach": 2},
        "atmosphere_mach": {"altitude": 10000, "mach": 2},
        "atmosphere_velocity": {"altitude": 10000, "velocity": 900},
    }
    return make_configuration(
        calculator="flow_sequence",
        mode="single",
        inputs_si={
            "initial": {"basis": basis, **values[basis]},
            "steps": [_default_step("oblique"), _default_step("expansion")],
        },
        models={"gas_model": gas},
        units=UnitPreferences(),
    )


def button(app: AppTest, label: str) -> Button:
    return next(item for item in app.button if item.label == label)


@pytest.mark.parametrize(
    "basis",
    [
        "static_mach",
        "static_velocity",
        "total",
        "atmosphere_mach",
        "atmosphere_velocity",
    ],
)
def test_replay_all_initial_modes(basis: str) -> None:
    config = configuration(basis)
    replay = cast(dict[str, Any], load_configuration(dump_configuration(config)))
    a = calculate_sequence(config["inputs_si"], config["models"])
    b = calculate_sequence(replay["inputs_si"], replay["models"])
    assert a == b
    assert a.complete
    assert len(sequence_rows(a)) == 3


@pytest.mark.parametrize("gas", ["AIR", "NASA7", "NASA9", "HARMONIC_OSCILLATOR"])
def test_gui_saved_config_and_units(gas: str) -> None:
    config = configuration(gas=gas)
    app = AppTest.from_string(SCRIPT, default_timeout=30)
    app.session_state["sequence_draft"] = config
    app.session_state["sequence_steps"] = copy.deepcopy(config["inputs_si"]["steps"])
    app.run()
    assert not app.exception
    button(app, "計算").click().run()
    assert not app.exception
    assert not app.error
    previous, result = app.session_state["sequence_result"]
    assert result.complete
    assert len(app.dataframe[0].value) == 3
    app.selectbox(key="unit_temperature").set_value("°C").run()
    assert app.number_input(
        key="sequence_initial_1_static_mach_temperature_°C"
    ).value == pytest.approx(226.85)
    assert app.warning
    assert (
        cast(dict[str, Any], load_configuration(previous))["models"]["gas_model"] == gas
    )
    button(app, "計算").click().run()
    assert app.session_state["sequence_result"][
        1
    ].initial.state.temperature == pytest.approx(500)


def test_editor_reorders_clones_deletes_and_keeps_initial() -> None:
    app = AppTest.from_string(SCRIPT, default_timeout=30).run()
    app.number_input(key="sequence_initial_0_static_mach_mach").set_value(4).run()
    button(app, "段を追加").click().run()
    button(app, "複製").click().run()
    assert len(app.session_state["sequence_steps"]) == 2
    button(app, "↓").click().run()
    app.button(key="sequence_2_1_up").click().run()
    assert app.number_input(key="sequence_initial_3_static_mach_mach").value == 4
    button(app, "削除").click().run()
    assert len(app.session_state["sequence_steps"]) == 1
    button(app, "計算").click().run()
    assert app.session_state["sequence_result"][1].initial.state.mach == 4
    assert not app.exception


@pytest.mark.parametrize("kind", ["conical", "isentropic", "expansion"])
def test_stage_controls(kind: str) -> None:
    app = AppTest.from_string(SCRIPT, default_timeout=30)
    app.session_state["sequence_steps"] = [_default_step(kind)]
    app.run()
    if kind == "isentropic":
        app.selectbox(key="sequence_0_0_basis").set_value("area_ratio").run()
        app.selectbox(key="sequence_0_0_exit").set_value("subsonic").run()
    button(app, "計算").click().run()
    assert not app.exception
    assert not app.error
    assert app.session_state["sequence_result"][1].complete


def test_failed_and_invalid_initial_gui() -> None:
    app = AppTest.from_string(SCRIPT, default_timeout=30)
    app.session_state["sequence_steps"] = [
        _default_step("normal"),
        _default_step("normal"),
        _default_step("normal"),
    ]
    app.run()
    button(app, "計算").click().run()
    assert app.error
    rows = app.dataframe[0].value
    assert rows["status"].tolist() == ["ok", "ok", "invalid_input", "not_computed"]
    app.number_input(key="sequence_initial_0_static_mach_pressure_Pa").set_value(
        -1
    ).run()
    button(app, "計算").click().run()
    assert app.error
    assert not app.exception
    assert "sequence_result" not in app.session_state


@pytest.mark.parametrize(
    "basis", ["total", "atmosphere_mach", "atmosphere_velocity", "static_velocity"]
)
def test_initial_controls(basis: str) -> None:
    app = AppTest.from_string(SCRIPT, default_timeout=30).run()
    app.selectbox(key="sequence_initial_basis_0").set_value(basis).run()
    button(app, "計算").click().run()
    assert not app.exception
    assert not app.error


def test_warning_and_display_units() -> None:
    config = configuration(gas="HARMONIC_OSCILLATOR")
    config["inputs_si"]["initial"].update(temperature=1500, mach=3)
    config["inputs_si"]["steps"] = []
    app = AppTest.from_string(SCRIPT, default_timeout=30)
    app.session_state["sequence_draft"] = config
    app.run()
    button(app, "計算").click().run()
    assert app.warning
    assert not app.exception
    rows = sequence_rows(app.session_state["sequence_result"][1])
    shown = _display_rows(rows, UnitPreferences(temperature="°C", pressure="kPa"))
    assert shown[0]["total_temperature [°C]"] is None
    assert shown[0]["pressure [kPa]"] == 10


@pytest.mark.parametrize(
    "data",
    [
        {},
        {"kind": "bad", "name": "x"},
        {"kind": "normal", "name": "x", "extra": 1},
        {"kind": "oblique", "name": "x", "angle": math.nan, "branch": "weak"},
        {"kind": "expansion", "name": "x", "angle": True},
        {"kind": "isentropic", "name": "x", "basis": "bad", "value": 1, "branch": None},
        {
            "kind": "isentropic",
            "name": "x",
            "basis": "area_ratio",
            "value": 1,
            "branch": None,
        },
    ],
)
def test_reject_bad_step(data: dict[str, Any]) -> None:
    with pytest.raises(ValueError):
        decode_step(data)


@pytest.mark.parametrize(
    "inputs,models",
    [
        ({}, {"gas_model": "AIR"}),
        ({"initial": {}, "steps": []}, {"gas_model": "AIR"}),
        ({"initial": {"basis": "static_mach"}, "steps": []}, {"gas_model": "AIR"}),
        (
            {
                "initial": {
                    "basis": "static_mach",
                    "pressure": 1,
                    "temperature": 300,
                    "mach": 2,
                },
                "steps": {},
            },
            {"gas_model": "AIR"},
        ),
        ({}, {}),
        ({}, {"gas_model": "bad"}),
    ],
)
def test_reject_bad_payload(inputs: dict[str, Any], models: dict[str, Any]) -> None:
    with pytest.raises(ValueError):
        validate_sequence_payload(inputs, models)


def test_envelope_rejects_sweep_and_bad_payload() -> None:
    config = configuration()
    config["mode"] = "sweep"
    with pytest.raises(ConfigurationError):
        load_configuration(__import__("json").dumps(config))
    config["mode"] = "single"
    config["inputs_si"] = {}
    with pytest.raises(ConfigurationError):
        load_configuration(__import__("json").dumps(config))


def test_pending_import_and_repeated_unit_changes() -> None:
    config = configuration()
    app = AppTest.from_string(SCRIPT, default_timeout=30)
    app.session_state["pending_flow_sequence_configuration"] = config
    app.run()
    assert len(app.session_state["sequence_steps"]) == 2
    button(app, "設定を適用").click().run()
    assert app.error
    app.selectbox(key="unit_temperature").set_value("°C").run()
    app.number_input(key="sequence_initial_2_static_mach_temperature_°C").set_value(
        300
    ).run()
    app.selectbox(key="unit_temperature").set_value("K").run()
    assert app.number_input(
        key="sequence_initial_3_static_mach_temperature_K"
    ).value == pytest.approx(573.15)
    button(app, "計算").click().run()
    assert not app.exception
    assert app.session_state["sequence_result"][
        1
    ].initial.state.temperature == pytest.approx(573.15)


@pytest.mark.parametrize("basis", ["atmosphere_mach", "atmosphere_velocity"])
def test_sequence_geopotential_replay(basis: str) -> None:
    config = configuration(basis)
    config["inputs_si"]["initial"].update(
        altitude=11000.0, altitude_basis="geopotential"
    )
    replay = load_configuration(dump_configuration(config))
    result = calculate_sequence(
        cast(dict[str, Any], replay["inputs_si"]),
        cast(dict[str, Any], replay["models"]),
    )
    rows = sequence_rows(result)
    assert rows[0]["temperature"] == pytest.approx(216.65)
    app = AppTest.from_string(SCRIPT, default_timeout=30)
    app.session_state["pending_flow_sequence_configuration"] = replay
    app.run()
    assert not app.exception
    assert app.radio(key="sequence_altitude_basis_1").value == "geopotential"
    app.button[-1].click().run()
    assert not app.exception
    assert not app.error


def test_sequence_unknown_altitude_coordinate() -> None:
    config = configuration("atmosphere_mach")
    config["inputs_si"]["initial"]["altitude_basis"] = "unknown"
    with pytest.raises(ValueError, match="altitude basis"):
        validate_sequence_payload(config["inputs_si"], config["models"])
