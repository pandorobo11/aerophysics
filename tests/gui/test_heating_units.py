"""Heating sweep input-state regression tests."""

from types import SimpleNamespace

import pytest

from aerophysics.gui import components
from aerophysics.gui.units import UnitPreferences, from_si


@pytest.mark.parametrize("old_unit,new_unit", [("deg", "rad"), ("rad", "deg")])
@pytest.mark.parametrize(
    "field", ["sweep_angle", "upstream_mach", "geometric_altitude"]
)
def test_heating_sweep_angle_unit_changes(
    monkeypatch: pytest.MonkeyPatch, old_unit: str, new_unit: str, field: str
) -> None:
    old = UnitPreferences(angle=old_unit)
    new = UnitPreferences(angle=new_unit)
    state: dict[str, object] = {
        f"unit_{name}": value for name, value in new.to_dict().items()
    }
    state["_gui_display_units"] = old.to_dict()
    state["heating_sweep_field"] = field
    expected = {}
    for key, value in (("heating_sweep_start", 0.2), ("heating_sweep_stop", 1.2)):
        state[key] = (
            float(from_si(value, "angle", old_unit))
            if field == "sweep_angle"
            else value
        )
        expected[key] = (
            float(from_si(value, "angle", new_unit))
            if field == "sweep_angle"
            else value
        )
    monkeypatch.setattr(components, "st", SimpleNamespace(session_state=state))

    components._convert_display_input_state()

    for key, value in expected.items():
        assert state[key] == pytest.approx(value)
    assert state["_gui_display_units"] == new.to_dict()


@pytest.mark.parametrize("field", ["radius", "geometric_altitude"])
def test_heating_length_sweep_unit_changes(
    monkeypatch: pytest.MonkeyPatch, field: str
) -> None:
    old = UnitPreferences(length="m")
    new = UnitPreferences(length="mm")
    state: dict[str, object] = {
        f"unit_{name}": value for name, value in new.to_dict().items()
    }
    state["_gui_display_units"] = old.to_dict()
    state["heating_sweep_field"] = field
    state["heating_sweep_start"] = 0.05
    state["heating_sweep_stop"] = 0.5
    monkeypatch.setattr(components, "st", SimpleNamespace(session_state=state))

    components._convert_display_input_state()

    assert state["heating_sweep_start"] == pytest.approx(50.0)
    assert state["heating_sweep_stop"] == pytest.approx(500.0)
