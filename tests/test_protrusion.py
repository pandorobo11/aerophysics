"""Tests for boundary-layer-immersed protrusion drag."""

import numpy as np
import pytest

from aerophysics.exceptions import ApplicabilityWarning
from aerophysics.protrusion import ProtrusionProfile, protrusion_drag


def test_constant_width_one_seventh_power_matches_analytic_solution() -> None:
    height = 0.02
    thickness = 0.1
    result = protrusion_drag(
        1.2,
        height,
        0.01,
        edge_velocity=100.0,
        edge_density=1.2,
        boundary_layer_thickness=thickness,
    )
    expected_shielding = 7.0 / 9.0 * (height / thickness) ** (2.0 / 7.0)
    expected_area = height * 0.01
    expected_edge_pressure = 0.5 * 1.2 * 100.0**2
    assert result.frontal_area == pytest.approx(expected_area)
    assert result.edge_dynamic_pressure == pytest.approx(expected_edge_pressure)
    assert result.shielding_factor == pytest.approx(expected_shielding, rel=2e-4)
    assert result.effective_dynamic_pressure == pytest.approx(
        expected_edge_pressure * expected_shielding, rel=2e-4
    )
    assert result.direct_drag == pytest.approx(
        1.2 * expected_edge_pressure * expected_area * expected_shielding,
        rel=2e-4,
    )
    assert result.height_to_boundary_layer_thickness == pytest.approx(0.2)
    assert result.profile is ProtrusionProfile.TURBULENT_ONE_SEVENTH_POWER
    assert result.compressibility_applied is False


def test_protrusion_above_boundary_layer_approaches_edge_dynamic_pressure() -> None:
    result = protrusion_drag(
        1.0,
        100.0,
        1.0,
        edge_velocity=10.0,
        edge_density=1.0,
        boundary_layer_thickness=0.1,
    )
    assert result.shielding_factor == pytest.approx(1.0, abs=3e-4)


def test_height_dependent_frontal_width_is_integrated() -> None:
    result = protrusion_drag(
        1.0,
        0.2,
        lambda y: 0.1 * (1.0 - y / 0.2),
        edge_velocity=20.0,
        edge_density=1.0,
        boundary_layer_thickness=0.5,
    )
    assert result.frontal_area == pytest.approx(0.01)
    assert 0.0 < result.shielding_factor < 1.0


def test_provided_profile_is_interpolated_and_integrated() -> None:
    result = protrusion_drag(
        2.0,
        1.0,
        0.5,
        edge_velocity=10.0,
        edge_density=2.0,
        boundary_layer_thickness=1.0,
        profile_height=[0.0, 0.5, 1.0],
        profile_velocity=[0.0, 5.0, 10.0],
        profile_density=[2.0, 2.0, 2.0],
    )
    assert result.effective_dynamic_pressure == pytest.approx(100.0 / 3.0)
    assert result.shielding_factor == pytest.approx(1.0 / 3.0)
    assert result.frontal_area == pytest.approx(0.5)
    assert result.direct_drag == pytest.approx(100.0 / 3.0)
    assert result.profile is ProtrusionProfile.PROVIDED
    assert result.compressibility_applied is False


def test_provided_profile_uses_edge_conditions_above_boundary_layer() -> None:
    result = protrusion_drag(
        1.0,
        2.0,
        1.0,
        edge_velocity=10.0,
        edge_density=2.0,
        boundary_layer_thickness=1.0,
        profile_height=[0.0, 1.0],
        profile_velocity=[0.0, 10.0],
        profile_density=[2.0, 2.0],
    )
    assert result.effective_dynamic_pressure == pytest.approx(200.0 / 3.0)
    assert result.shielding_factor == pytest.approx(2.0 / 3.0)


def test_walz_compressibility_changes_density_and_drag() -> None:
    incompressible = protrusion_drag(1.0, 0.05, 0.01, 500.0, 0.4, 0.1)
    compressible = protrusion_drag(
        1.0,
        0.05,
        0.01,
        500.0,
        0.4,
        0.1,
        mach=2.0,
        edge_temperature=220.0,
    )
    assert compressible.compressibility_applied is True
    assert compressible.direct_drag < incompressible.direct_drag
    assert 0.0 < compressible.shielding_factor < 1.0


def test_zero_mach_compressible_model_reduces_to_incompressible_model() -> None:
    incompressible = protrusion_drag(1.0, 0.05, 0.01, 100.0, 1.0, 0.1)
    compressible = protrusion_drag(
        1.0,
        0.05,
        0.01,
        100.0,
        1.0,
        0.1,
        mach=0.0,
        edge_temperature=300.0,
    )
    assert compressible.shielding_factor == pytest.approx(
        incompressible.shielding_factor
    )
    assert compressible.direct_drag == pytest.approx(incompressible.direct_drag)


def test_transonic_scalar_drag_coefficient_warns() -> None:
    with pytest.warns(ApplicabilityWarning, match="transonic"):
        protrusion_drag(
            1.0,
            0.05,
            0.01,
            300.0,
            1.0,
            0.1,
            mach=1.0,
            edge_temperature=250.0,
        )


@pytest.mark.parametrize(
    ("keyword", "value"),
    [
        ("drag_coefficient", -1.0),
        ("height", 0.0),
        ("edge_velocity", 0.0),
        ("edge_density", -1.0),
        ("boundary_layer_thickness", 0.0),
        ("prandtl_number", np.nan),
    ],
)
def test_invalid_scalar_inputs(keyword: str, value: float) -> None:
    arguments = {
        "drag_coefficient": 1.0,
        "height": 0.1,
        "frontal_width": 0.01,
        "edge_velocity": 10.0,
        "edge_density": 1.0,
        "boundary_layer_thickness": 0.2,
        "prandtl_number": 0.72,
    }
    arguments[keyword] = value
    with pytest.raises(ValueError, match=keyword):
        protrusion_drag(**arguments)  # type: ignore[arg-type]


def test_scalar_inputs_reject_arrays() -> None:
    with pytest.raises(ValueError, match="drag_coefficient must be a scalar"):
        protrusion_drag([1.0], 0.1, 0.01, 10.0, 1.0, 0.2)  # type: ignore[arg-type]


@pytest.mark.parametrize("integration_points", [True, 31, 32.0])
def test_invalid_integration_points(integration_points: object) -> None:
    with pytest.raises(ValueError, match="integration_points"):
        protrusion_drag(
            1.0,
            0.1,
            0.01,
            10.0,
            1.0,
            0.2,
            integration_points=integration_points,  # type: ignore[arg-type]
        )


def test_profile_inputs_must_be_complete() -> None:
    with pytest.raises(ValueError, match="supplied together"):
        protrusion_drag(
            1.0,
            0.1,
            0.01,
            10.0,
            1.0,
            0.2,
            profile_height=[0.0, 0.1],
        )


@pytest.mark.parametrize(
    ("profile_height", "profile_velocity", "profile_density", "match"),
    [
        ([[0.0, 0.1]], [0.0, 10.0], [1.0, 1.0], "one-dimensional"),
        ([0.0], [0.0], [1.0], "length"),
        ([0.01, 0.1], [0.0, 10.0], [1.0, 1.0], "start at zero"),
        ([0.0, 0.1, 0.05], [0.0, 10.0, 5.0], [1.0, 1.0, 1.0], "increasing"),
        ([0.0, 0.05], [0.0, 5.0], [1.0, 1.0], "cover"),
        ([0.0, 0.1], [0.0, -1.0], [1.0, 1.0], "profile_velocity"),
        ([0.0, 0.1], [0.0, 10.0], [1.0, 0.0], "profile_density"),
    ],
)
def test_invalid_provided_profiles(
    profile_height: object,
    profile_velocity: object,
    profile_density: object,
    match: str,
) -> None:
    with pytest.raises(ValueError, match=match):
        protrusion_drag(
            1.0,
            0.1,
            0.01,
            10.0,
            1.0,
            0.2,
            profile_height=profile_height,  # type: ignore[arg-type]
            profile_velocity=profile_velocity,  # type: ignore[arg-type]
            profile_density=profile_density,  # type: ignore[arg-type]
        )


def test_thermal_inputs_are_exclusive_and_complete() -> None:
    with pytest.raises(ValueError, match="supplied together"):
        protrusion_drag(1.0, 0.1, 0.01, 10.0, 1.0, 0.2, mach=2.0)
    with pytest.raises(ValueError, match="wall_temperature requires"):
        protrusion_drag(1.0, 0.1, 0.01, 10.0, 1.0, 0.2, wall_temperature=300.0)
    with pytest.raises(ValueError, match="cannot be combined"):
        protrusion_drag(
            1.0,
            0.1,
            0.01,
            10.0,
            1.0,
            0.2,
            profile_height=[0.0, 0.1],
            profile_velocity=[0.0, 10.0],
            profile_density=[1.0, 1.0],
            mach=2.0,
            edge_temperature=250.0,
        )


@pytest.mark.parametrize(
    ("keyword", "value"),
    [
        ("mach", -1.0),
        ("edge_temperature", 0.0),
        ("wall_temperature", -1.0),
    ],
)
def test_invalid_thermal_values(keyword: str, value: float) -> None:
    arguments: dict[str, float] = {
        "mach": 2.0,
        "edge_temperature": 250.0,
        "wall_temperature": 300.0,
    }
    arguments[keyword] = value
    with pytest.raises(ValueError, match=keyword):
        protrusion_drag(
            1.0,
            0.1,
            0.01,
            10.0,
            1.0,
            0.2,
            **arguments,  # type: ignore[arg-type]
        )


def test_invalid_frontal_width() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        protrusion_drag(1.0, 0.1, -0.01, 10.0, 1.0, 0.2)
    with pytest.raises(ValueError, match="positive frontal area"):
        protrusion_drag(1.0, 0.1, 0.0, 10.0, 1.0, 0.2)
    with pytest.raises(ValueError, match="matching"):
        protrusion_drag(
            1.0,
            0.1,
            lambda _: [0.1, 0.2],
            10.0,
            1.0,
            0.2,
        )
    with pytest.raises(ValueError, match="could not be evaluated"):
        protrusion_drag(
            1.0,
            0.1,
            lambda _: (_ for _ in ()).throw(RuntimeError("failure")),
            10.0,
            1.0,
            0.2,
        )
