"""Regression tests for laminar-profile input unit changes."""

from typing import cast

import pytest

from aerophysics.gui import components
from aerophysics.gui.units import QuantityKind, UnitPreferences, to_si


@pytest.mark.parametrize(
    ("kind", "unit"),
    [
        ("length", "mm"),
        ("speed", "ft/s"),
        ("density", "lbm/ft³"),
        ("temperature", "°C"),
        ("temperature", "°F"),
    ],
)
def test_laminar_profile_input_units_preserve_physical_values(
    monkeypatch: pytest.MonkeyPatch, kind: QuantityKind, unit: str
) -> None:
    original = UnitPreferences()
    inputs: dict[str, tuple[QuantityKind, float]] = {
        "laminar_profile_distance": ("length", 0.2),
        "laminar_profile_velocity": ("speed", 340.0),
        "laminar_profile_density": ("density", 0.91),
        "laminar_profile_temperature": ("temperature", 275.5),
        "laminar_profile_wall_temperature": ("temperature", 300.0),
    }
    state: dict[str, object] = {
        "_gui_display_units": original.to_dict(),
        f"unit_{kind}": unit,
        "laminar_profile_points": 51,
        **{key: value for key, (_, value) in inputs.items()},
    }
    monkeypatch.setattr(components.st, "session_state", state)

    components._convert_display_input_state()
    converted = UnitPreferences.from_dict(state["_gui_display_units"])
    for key, (quantity, value) in inputs.items():
        assert float(
            to_si(cast(float, state[key]), quantity, getattr(converted, quantity))
        ) == pytest.approx(value)
    assert state["laminar_profile_points"] == 51

    state[f"unit_{kind}"] = getattr(original, kind)
    components._convert_display_input_state()
    for key, (_, value) in inputs.items():
        assert state[key] == pytest.approx(value)
    assert state["_gui_display_units"] == original.to_dict()
