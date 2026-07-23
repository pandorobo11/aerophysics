"""Tests for smooth flat-plate boundary-layer correlations."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    TurbulentCorrelation,
    flat_plate_boundary_layer,
)
from aerophysics.exceptions import ApplicabilityWarning, ModelRangeError


def test_blasius_laminar_reference_values() -> None:
    result = flat_plate_boundary_layer(
        1.0,
        10.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.LAMINAR,
    )
    assert result.reynolds_number == pytest.approx(1e6)
    assert result.boundary_layer_thickness == pytest.approx(0.005)
    assert result.displacement_thickness == pytest.approx(0.0017208)
    assert result.momentum_thickness == pytest.approx(0.000664)
    assert result.local_skin_friction_coefficient == pytest.approx(0.000664)
    assert result.average_skin_friction_coefficient == pytest.approx(0.001328)
    assert result.wall_shear_stress == pytest.approx(0.0332)
    assert result.drag_per_unit_width == pytest.approx(0.0664)


def test_one_fifth_power_turbulent_reference_values() -> None:
    result = flat_plate_boundary_layer(
        1.0,
        10.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TURBULENT,
        turbulent_correlation=TurbulentCorrelation.POWER_LAW,
    )
    factor = 1e6**-0.2
    assert result.boundary_layer_thickness == pytest.approx(0.37 * factor)
    assert result.displacement_thickness == pytest.approx(0.37 * factor / 8.0)
    assert result.momentum_thickness == pytest.approx(0.37 * factor * 7.0 / 72.0)
    assert result.local_skin_friction_coefficient == pytest.approx(0.0592 * factor)
    assert result.average_skin_friction_coefficient == pytest.approx(0.074 * factor)


def test_schlichting_is_default_and_local_is_drag_derivative() -> None:
    reynolds = 1e7
    result = flat_plate_boundary_layer(
        1.0,
        100.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TURBULENT,
    )
    average = 0.455 / np.log10(reynolds) ** 2.58
    assert result.average_skin_friction_coefficient == pytest.approx(average)
    assert result.local_skin_friction_coefficient == pytest.approx(
        average * (1.0 - 2.58 / np.log(reynolds))
    )


def test_explicit_transition_offsets_accumulated_drag() -> None:
    result = flat_plate_boundary_layer(
        [0.25, 0.5, 1.0],
        10.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TRANSITIONAL,
        turbulent_correlation=TurbulentCorrelation.POWER_LAW,
        transition_reynolds=5e5,
    )
    reynolds = np.array([2.5e5, 5e5, 1e6])
    laminar_average_at_transition = 1.328 / np.sqrt(5e5)
    turbulent_average = 0.074 * reynolds[-1] ** -0.2
    turbulent_average_at_transition = 0.074 * (5e5) ** -0.2
    expected_last = turbulent_average + 0.5 * (
        laminar_average_at_transition - turbulent_average_at_transition
    )
    assert_allclose(result.reynolds_number, reynolds)
    assert np.asarray(result.local_skin_friction_coefficient)[0] == pytest.approx(
        0.664 / np.sqrt(reynolds[0])
    )
    assert np.asarray(result.average_skin_friction_coefficient)[1] == pytest.approx(
        laminar_average_at_transition
    )
    assert np.asarray(result.average_skin_friction_coefficient)[2] == pytest.approx(
        expected_last
    )


def test_array_inputs_broadcast_and_return_float64() -> None:
    result = flat_plate_boundary_layer(
        [[0.5], [1.0]],
        [10.0, 20.0, 30.0],
        1.2,
        1.8e-5,
        regime=BoundaryLayerRegime.LAMINAR,
    )
    for value in (
        result.distance,
        result.reynolds_number,
        result.boundary_layer_thickness,
        result.displacement_thickness,
        result.momentum_thickness,
        result.local_skin_friction_coefficient,
        result.average_skin_friction_coefficient,
        result.wall_shear_stress,
        result.drag_per_unit_width,
    ):
        assert isinstance(value, np.ndarray)
        assert value.shape == (2, 3)
        assert value.dtype == np.float64


def test_turbulent_range_warning() -> None:
    with pytest.warns(ApplicabilityWarning):
        flat_plate_boundary_layer(
            1.0,
            1.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.TURBULENT,
        )


@pytest.mark.parametrize(
    ("name", "values"),
    [
        ("distance", (0.0, 1.0, 1.0, 1e-5)),
        ("edge_velocity", (1.0, 0.0, 1.0, 1e-5)),
        ("edge_density", (1.0, 1.0, -1.0, 1e-5)),
        ("edge_dynamic_viscosity", (1.0, 1.0, 1.0, -1.0)),
    ],
)
def test_invalid_physical_inputs(name: str, values: tuple[float, ...]) -> None:
    with pytest.raises(ValueError, match=name):
        flat_plate_boundary_layer(
            *values,
            regime=BoundaryLayerRegime.LAMINAR,
        )


def test_inputs_must_broadcast() -> None:
    with pytest.raises(ValueError, match="broadcastable"):
        flat_plate_boundary_layer(
            [1.0, 2.0],
            [1.0, 2.0, 3.0],
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.LAMINAR,
        )


def test_transition_reynolds_is_explicit_and_exclusive() -> None:
    with pytest.raises(ValueError, match="required"):
        flat_plate_boundary_layer(
            1.0,
            10.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.TRANSITIONAL,
        )
    with pytest.raises(ValueError, match="finite"):
        flat_plate_boundary_layer(
            1.0,
            10.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.TRANSITIONAL,
            transition_reynolds=-1.0,
        )
    with pytest.raises(ValueError, match="only valid"):
        flat_plate_boundary_layer(
            1.0,
            10.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.LAMINAR,
            transition_reynolds=5e5,
        )


def test_invalid_enum_values_are_rejected() -> None:
    with pytest.raises(ValueError, match="regime"):
        flat_plate_boundary_layer(
            1.0,
            10.0,
            1.0,
            1e-5,
            regime="laminar",  # type: ignore[arg-type]
        )
    with pytest.raises(ValueError, match="turbulent_correlation"):
        flat_plate_boundary_layer(
            1.0,
            10.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.LAMINAR,
            turbulent_correlation="power_law",  # type: ignore[arg-type]
        )


def test_schlichting_rejects_undefined_reynolds_number() -> None:
    with pytest.raises(ModelRangeError, match="greater than one"):
        flat_plate_boundary_layer(
            1e-6,
            1.0,
            1.0,
            1.0,
            regime=BoundaryLayerRegime.TURBULENT,
        )
