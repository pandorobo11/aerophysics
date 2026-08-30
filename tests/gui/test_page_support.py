"""Tests for the typed Streamlit page-support boundary."""

from __future__ import annotations

import numpy as np
import pytest
from streamlit.testing.v1 import AppTest

from aerophysics.gui import page_support
from aerophysics.gui.adapters import CalculationResult
from aerophysics.gui.page_support import PageDataError


class _FakeStreamlit:
    def __init__(self) -> None:
        self.session_state: dict[str, object] = {}
        self.metrics: list[tuple[str, str]] = []

    def metric(self, label: str, value: str) -> None:
        self.metrics.append((label, value))


def test_configuration_defaults_are_typed_copies_with_explicit_fallbacks() -> None:
    source = {
        "inputs_si": {"temperature": 300, "invalid": "hot"},
        "models": {"gas": "AIR"},
        "sweep_si": "not an object",
    }

    defaults = page_support.configuration_defaults(source)

    assert defaults.inputs == {"temperature": 300, "invalid": "hot"}
    assert defaults.models == {"gas": "AIR"}
    assert defaults.sweep == {}
    assert page_support.numeric_default(defaults.inputs, "temperature", 200.0) == 300.0
    assert page_support.numeric_default(defaults.inputs, "invalid", 200.0) == 200.0
    assert page_support.numeric_default(defaults.inputs, "missing", 250.0) == 250.0
    defaults.inputs["temperature"] = 400.0
    assert source["inputs_si"] == {"temperature": 300, "invalid": "hot"}


def test_page_unit_and_array_defaults() -> None:
    assert page_support.display_value(1.0, "length", "mm") == pytest.approx(1000.0)
    assert page_support.si_value(1000.0, "length", "mm") == pytest.approx(1.0)
    assert page_support.array_default({"values": None}, "values") is None
    stored = page_support.array_default({"values": [1.0, 2.0]}, "values")
    assert stored is not None
    np.testing.assert_allclose(stored, np.asarray([1.0, 2.0]))
    with pytest.raises(PageDataError, match="one-dimensional"):
        page_support.array_default({"values": [[1.0], [2.0]]}, "values")
    with pytest.raises(PageDataError, match="numeric array"):
        page_support.array_default({"values": ["not numeric"]}, "values")


def test_session_payload_validates_calculate_to_render_contract(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake = _FakeStreamlit()
    monkeypatch.setattr(page_support, "st", fake)
    result = CalculationResult(({"value": 1.0},))

    assert page_support.session_payload("result", CalculationResult) is None
    fake.session_state["result"] = (result, {"mode": "single"})
    assert page_support.session_payload("result", CalculationResult) == (
        result,
        {"mode": "single"},
    )

    fake.session_state["result"] = ("wrong result", {"mode": "single"})
    with pytest.raises(PageDataError, match="invalid result"):
        page_support.session_payload("result", CalculationResult)
    fake.session_state["result"] = (result,)
    with pytest.raises(PageDataError, match="result/configuration pair"):
        page_support.session_payload("result", CalculationResult)


def test_metric_fields_are_unique_numeric_and_rendered(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake = _FakeStreamlit()
    monkeypatch.setattr(page_support, "st", fake)

    row = {"下流 Mach M₂": 1.25, "status": "ok"}
    assert page_support.metric_value(row, "下流 Mach") == pytest.approx(1.25)
    page_support.render_metric(row, "M₂", "下流 Mach")
    assert fake.metrics == [("M₂", "1.25")]

    with pytest.raises(PageDataError, match="found 0"):
        page_support.metric_value(row, "圧力")
    with pytest.raises(PageDataError, match="found 2"):
        page_support.metric_value({"value a": 1.0, "value b": 2.0}, "value")
    assert page_support.metric_value({"pressure": None}, "pressure") is None
    page_support.render_metric({"pressure": None}, "p", "pressure")
    assert fake.metrics[-1] == ("p", "—")
    with pytest.raises(PageDataError, match="must be numeric"):
        page_support.metric_value({"pressure": "unknown"}, "pressure")


@pytest.mark.parametrize(
    ("module", "renderer", "form", "payload_key", "calculator", "metric"),
    (
        (
            "aerophysics.gui.flow_pages",
            "render_normal_shock",
            "normal_form",
            "normal_payload",
            "normal_shock",
            "M₂",
        ),
        (
            "aerophysics.gui.analysis_pages",
            "render_thermochemistry",
            "thermo_form",
            "thermo_payload",
            "thermochemistry",
            "c_p",
        ),
    ),
)
def test_migrated_pages_preserve_calculate_render_and_export_pipeline(
    module: str,
    renderer: str,
    form: str,
    payload_key: str,
    calculator: str,
    metric: str,
) -> None:
    script = f"""
from {module} import {renderer}
from aerophysics.gui.units import UnitPreferences
{renderer}(UnitPreferences())
"""
    app = AppTest.from_string(script, default_timeout=30).run()
    if calculator == "thermochemistry":
        app.radio(key="thermo_mode").set_value("single").run()
    app.button(key=f"FormSubmitter:{form}-計算").click().run()

    assert not app.exception
    assert metric in {item.label for item in app.metric}
    assert len(app.download_button) == 2
    payload = app.session_state[payload_key]
    assert isinstance(payload, tuple)
    assert isinstance(payload[1], dict)
    assert payload[1]["calculator"] == calculator


def test_protrusion_sweep_renders_expected_unavailable_first_row() -> None:
    script = """
from aerophysics.gui.analysis_pages import render_protrusion_drag
from aerophysics.gui.units import UnitPreferences
render_protrusion_drag(UnitPreferences())
"""
    app = AppTest.from_string(script, default_timeout=30).run()
    app.radio(key="protrusion_mode").set_value("1変数スイープ").run()
    app.number_input(key="protrusion_sweep_start").set_value(0.0).run()
    app.number_input(key="protrusion_sweep_points").set_value(3).run()
    app.button(key="FormSubmitter:protrusion_form-計算").click().run()

    assert not app.exception
    assert app.warning
    assert app.metric[0].value == "—"
    assert len(app.dataframe[0].value) == 3
    download_keys = {item.key for item in app.download_button}
    assert "aerophysics-protrusion-drag_csv" in download_keys
    assert "aerophysics-protrusion-drag_json" in download_keys
