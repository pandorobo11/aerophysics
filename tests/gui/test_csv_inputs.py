"""Tests for validated GUI CSV inputs."""

import pytest

from aerophysics.gui.csv_inputs import (
    parse_profile_csv,
    parse_shape_csv,
    profile_csv_template,
    shape_csv_template,
)
from aerophysics.gui.units import UnitPreferences


def test_profile_csv_uses_selected_display_units() -> None:
    profile = parse_profile_csv(
        "wall_distance,velocity,density\n0,0,0.002\n1,100,0.002\n",
        UnitPreferences(length="ft", speed="kt", density="slug/ft³"),
    )
    assert profile.wall_distance[-1] == pytest.approx(0.3048)
    assert profile.velocity[-1] == pytest.approx(51.4444, rel=1.0e-5)
    assert profile.density[0] == pytest.approx(1.03076, rel=1.0e-5)
    assert "wall_distance" in profile_csv_template()


def test_shape_csv_uses_selected_length_unit() -> None:
    shape = parse_shape_csv(b"height,width\n0,1\n2,0\n", UnitPreferences(length="ft"))
    assert shape.height[-1] == pytest.approx(0.6096)
    assert shape.width[0] == pytest.approx(0.3048)
    assert "height,width" in shape_csv_template()


@pytest.mark.parametrize(
    ("text", "message"),
    [
        ("x,velocity,density\n0,0,1\n1,1,1\n", "columns"),
        ("wall_distance,velocity,density\n0,0,1\n", "two"),
        ("wall_distance,velocity,density\n0,0,1\n0,1,1\n", "increasing"),
        ("wall_distance,velocity,density\n0,0,1\n1,-1,1\n", "non-negative"),
        ("wall_distance,velocity,density\n0,0,1\n1,1,0\n", "greater"),
        ("wall_distance,velocity,density\n0,0,1\n1,nan,1\n", "finite"),
        ("wall_distance,velocity,density\n0,0,1\n1,bad,1\n", "numeric"),
    ],
)
def test_profile_csv_rejects_invalid_data(text: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        parse_profile_csv(text, UnitPreferences())


def test_shape_csv_rejects_negative_width() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        parse_shape_csv("height,width\n0,1\n1,-1\n", UnitPreferences())
