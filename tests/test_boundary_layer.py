"""Tests for smooth flat-plate boundary-layer correlations."""

from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    CompressibilityCorrection,
    TurbulentCorrelation,
    _van_driest_ii_state,
    flat_plate_boundary_layer,
)
from aerophysics.exceptions import ApplicabilityWarning, ModelRangeError
from aerophysics.gas import AIR, AIR_VISCOSITY


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
    assert result.effective_reynolds_number == result.reynolds_number
    assert result.recovery_temperature is None
    assert result.wall_temperature is None


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
        result.effective_reynolds_number,
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


def test_eckert_laminar_adiabatic_wall() -> None:
    edge_temperature = 250.0
    mach = 2.0
    viscosity = float(AIR_VISCOSITY.dynamic_viscosity(edge_temperature))
    result = flat_plate_boundary_layer(
        1.0,
        500.0,
        0.5,
        viscosity,
        regime=BoundaryLayerRegime.LAMINAR,
        compressibility_correction=CompressibilityCorrection.ECKERT,
        mach=mach,
        edge_temperature=edge_temperature,
    )
    recovery = edge_temperature * (
        1.0 + np.sqrt(0.72) * 0.5 * (AIR.heat_capacity_ratio - 1.0) * mach**2
    )
    reference = 0.22 * recovery + 0.28 * edge_temperature + 0.50 * recovery
    edge_reynolds = 0.5 * 500.0 / viscosity
    effective_reynolds = (
        edge_reynolds
        * edge_temperature
        / reference
        * viscosity
        / float(AIR_VISCOSITY.dynamic_viscosity(reference))
    )
    assert result.recovery_temperature == pytest.approx(recovery)
    assert result.wall_temperature == pytest.approx(recovery)
    assert result.effective_reynolds_number == pytest.approx(effective_reynolds)
    assert result.local_skin_friction_coefficient == pytest.approx(
        0.664 / np.sqrt(effective_reynolds)
    )


def test_eckert_turbulent_uses_specified_wall_temperature() -> None:
    result = flat_plate_boundary_layer(
        1.0,
        100.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TURBULENT,
        turbulent_correlation=TurbulentCorrelation.POWER_LAW,
        compressibility_correction=CompressibilityCorrection.ECKERT,
        mach=3.0,
        edge_temperature=220.0,
        wall_temperature=300.0,
    )
    recovery = 220.0 * (
        1.0 + np.cbrt(0.72) * 0.5 * (AIR.heat_capacity_ratio - 1.0) * 3.0**2
    )
    reference = 0.22 * recovery + 0.28 * 220.0 + 0.50 * 300.0
    expected_reynolds = (
        1e7
        * 220.0
        / reference
        * float(AIR_VISCOSITY.dynamic_viscosity(220.0))
        / float(AIR_VISCOSITY.dynamic_viscosity(reference))
    )
    assert result.recovery_temperature == pytest.approx(recovery)
    assert result.wall_temperature == 300.0
    assert result.effective_reynolds_number == pytest.approx(expected_reynolds)
    assert result.average_skin_friction_coefficient == pytest.approx(
        0.074 * expected_reynolds**-0.2
    )


def _van_driest_factors(
    mach: float,
    edge_temperature: float,
    wall_temperature: float | None = None,
) -> tuple[float, float, float, float]:
    recovery_factor = np.cbrt(0.72)
    m = (
        recovery_factor
        * (AIR.heat_capacity_ratio - 1.0)
        * mach**2
        / 2.0
    )
    recovery = edge_temperature * (1.0 + m)
    wall = recovery if wall_temperature is None else wall_temperature
    if m <= 1e-12:
        return 1.0, 1.0, 1.0, recovery
    temperature_factor = wall / edge_temperature
    denominator = np.sqrt(
        (m + 1.0 + temperature_factor) ** 2 - 4.0 * temperature_factor
    )
    alpha = np.clip(
        (m - 1.0 + temperature_factor) / denominator, -1.0, 1.0
    )
    beta = np.clip(
        (m + 1.0 - temperature_factor) / denominator, -1.0, 1.0
    )
    friction_factor = m / (np.arcsin(alpha) + np.arcsin(beta)) ** 2
    momentum_factor = float(
        AIR_VISCOSITY.dynamic_viscosity(edge_temperature)
    ) / float(AIR_VISCOSITY.dynamic_viscosity(wall))
    return (
        friction_factor,
        momentum_factor,
        momentum_factor / friction_factor,
        recovery,
    )


def test_van_driest_factors_and_effective_reynolds() -> None:
    edge_temperature = 250.0
    mach = 2.0
    reynolds = np.asarray(1e7)
    state = _van_driest_ii_state(
        reynolds,
        np.asarray(edge_temperature),
        np.asarray(mach),
        None,
        prandtl_number=0.72,
        gas=AIR,
        viscosity_model=AIR_VISCOSITY,
    )
    friction_factor, momentum_factor, reynolds_factor, recovery = (
        _van_driest_factors(mach, edge_temperature)
    )
    assert state.friction_factor == pytest.approx(friction_factor)
    assert state.momentum_factor == pytest.approx(momentum_factor)
    assert state.reynolds_factor == pytest.approx(
        state.momentum_factor / state.friction_factor
    )
    assert state.reynolds_factor == pytest.approx(reynolds_factor)
    assert state.friction_reynolds == pytest.approx(reynolds * reynolds_factor)
    assert state.thickness_reynolds == pytest.approx(reynolds * momentum_factor)

    result = flat_plate_boundary_layer(
        1.0,
        100.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TURBULENT,
        turbulent_correlation=TurbulentCorrelation.POWER_LAW,
        compressibility_correction=CompressibilityCorrection.VAN_DRIEST_II,
        mach=mach,
        edge_temperature=edge_temperature,
    )
    assert result.effective_reynolds_number == pytest.approx(
        1e7 * reynolds_factor
    )
    assert result.wall_temperature == pytest.approx(recovery)
    assert result.boundary_layer_thickness == pytest.approx(
        0.37 * (1e7 * momentum_factor) ** -0.2
    )


@pytest.mark.parametrize(
    ("mach", "wall_temperature"),
    [(2.0, None), (5.0, None), (5.0, 300.0), (8.0, 300.0)],
)
def test_van_driest_local_and_average_implicit_residuals(
    mach: float,
    wall_temperature: float | None,
) -> None:
    edge_temperature = 220.0
    result = flat_plate_boundary_layer(
        1.0,
        100.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TURBULENT,
        compressibility_correction=CompressibilityCorrection.VAN_DRIEST_II,
        mach=mach,
        edge_temperature=edge_temperature,
        wall_temperature=wall_temperature,
    )
    friction_factor, momentum_factor, reynolds_factor, _ = _van_driest_factors(
        mach, edge_temperature, wall_temperature
    )
    assert reynolds_factor == pytest.approx(momentum_factor / friction_factor)
    assert result.effective_reynolds_number == pytest.approx(
        result.reynolds_number * reynolds_factor
    )

    local = float(result.local_skin_friction_coefficient)
    average = float(result.average_skin_friction_coefficient)
    local_residual = 0.242 / np.sqrt(local * friction_factor) - (
        0.41
        + np.log10(
            float(result.reynolds_number)
            * reynolds_factor
            * local
            * friction_factor
        )
    )
    average_residual = 0.242 / np.sqrt(average * friction_factor) - np.log10(
        float(result.reynolds_number)
        * reynolds_factor
        * average
        * friction_factor
    )
    assert local_residual == pytest.approx(0.0, abs=2e-14)
    assert average_residual == pytest.approx(0.0, abs=2e-14)
    assert result.drag_per_unit_width == pytest.approx(
        0.5 * 1.0 * 100.0**2 * 1.0 * average
    )


def test_willems_equation_7_direct_local_residual() -> None:
    friction_factor, _, reynolds_factor, _ = _van_driest_factors(5.0, 220.0)
    result = flat_plate_boundary_layer(
        1.0,
        100.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TURBULENT,
        compressibility_correction=CompressibilityCorrection.VAN_DRIEST_II,
        mach=5.0,
        edge_temperature=220.0,
    )
    local_i = (
        float(result.local_skin_friction_coefficient) * friction_factor
    )
    reynolds_i = float(result.reynolds_number) * reynolds_factor
    residual = (
        0.242 / np.sqrt(local_i)
        - 0.41
        - np.log10(reynolds_i * local_i)
    )
    assert residual == pytest.approx(0.0, abs=2e-14)


def test_van_driest_ignores_turbulent_correlation_selection() -> None:
    arguments: dict[str, Any] = {
        "distance": 1.0,
        "edge_velocity": 100.0,
        "edge_density": 1.0,
        "edge_dynamic_viscosity": 1e-5,
        "regime": BoundaryLayerRegime.TURBULENT,
        "compressibility_correction": CompressibilityCorrection.VAN_DRIEST_II,
        "mach": 5.0,
        "edge_temperature": 220.0,
    }
    power = flat_plate_boundary_layer(
        **arguments,
        turbulent_correlation=TurbulentCorrelation.POWER_LAW,
    )
    schlichting = flat_plate_boundary_layer(
        **arguments,
        turbulent_correlation=TurbulentCorrelation.SCHLICHTING,
    )
    assert power.local_skin_friction_coefficient == pytest.approx(
        schlichting.local_skin_friction_coefficient
    )
    assert power.average_skin_friction_coefficient == pytest.approx(
        schlichting.average_skin_friction_coefficient
    )


def test_van_driest_transition_uses_eckert_then_van_driest() -> None:
    mixed = flat_plate_boundary_layer(
        [0.1, 1.0],
        100.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TRANSITIONAL,
        transition_reynolds=2e6,
        turbulent_correlation=TurbulentCorrelation.POWER_LAW,
        compressibility_correction=CompressibilityCorrection.VAN_DRIEST_II,
        mach=2.0,
        edge_temperature=250.0,
    )
    laminar = flat_plate_boundary_layer(
        0.1,
        100.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.LAMINAR,
        compressibility_correction=CompressibilityCorrection.ECKERT,
        mach=2.0,
        edge_temperature=250.0,
    )
    assert np.asarray(mixed.local_skin_friction_coefficient)[0] == pytest.approx(
        laminar.local_skin_friction_coefficient
    )
    assert np.asarray(mixed.recovery_temperature)[0] == pytest.approx(
        laminar.recovery_temperature
    )
    assert (
        np.asarray(mixed.recovery_temperature)[1]
        > np.asarray(mixed.recovery_temperature)[0]
    )


def test_van_driest_transition_accumulated_drag_is_continuous() -> None:
    transition = 2e6
    x_transition = transition * 1e-5 / 100.0
    distances = np.array([x_transition, x_transition * (1.0 + 1e-9)])
    result = flat_plate_boundary_layer(
        distances,
        100.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TRANSITIONAL,
        transition_reynolds=transition,
        compressibility_correction=CompressibilityCorrection.VAN_DRIEST_II,
        mach=5.0,
        edge_temperature=220.0,
        wall_temperature=300.0,
    )
    drag = np.asarray(result.drag_per_unit_width)
    assert drag[1] == pytest.approx(drag[0], rel=2e-9)


def test_compressibility_inputs_broadcast() -> None:
    result = flat_plate_boundary_layer(
        [[0.5], [1.0]],
        100.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.LAMINAR,
        compressibility_correction=CompressibilityCorrection.ECKERT,
        mach=[1.0, 2.0, 3.0],
        edge_temperature=250.0,
    )
    assert isinstance(result.wall_temperature, np.ndarray)
    assert result.wall_temperature.shape == (2, 3)
    assert isinstance(result.effective_reynolds_number, np.ndarray)
    assert result.effective_reynolds_number.shape == (2, 3)


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


def test_compressibility_selection_and_required_inputs() -> None:
    with pytest.raises(ValueError, match="thermal inputs"):
        flat_plate_boundary_layer(
            1.0,
            10.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.LAMINAR,
            mach=1.0,
        )
    with pytest.raises(ValueError, match="required"):
        flat_plate_boundary_layer(
            1.0,
            10.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.LAMINAR,
            compressibility_correction=CompressibilityCorrection.ECKERT,
            mach=1.0,
        )
    with pytest.raises(ValueError, match="compressibility_correction"):
        flat_plate_boundary_layer(
            1.0,
            10.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.LAMINAR,
            compressibility_correction="eckert",  # type: ignore[arg-type]
        )


@pytest.mark.parametrize(
    ("keyword", "value", "match"),
    [
        ("mach", -1.0, "mach"),
        ("edge_temperature", 0.0, "edge_temperature"),
        ("wall_temperature", -1.0, "wall_temperature"),
    ],
)
def test_invalid_thermal_inputs(keyword: str, value: float, match: str) -> None:
    arguments: dict[str, object] = {
        "mach": 1.0,
        "edge_temperature": 250.0,
        "wall_temperature": 300.0,
    }
    arguments[keyword] = value
    with pytest.raises(ValueError, match=match):
        flat_plate_boundary_layer(
            1.0,
            10.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.LAMINAR,
            compressibility_correction=CompressibilityCorrection.ECKERT,
            **arguments,  # type: ignore[arg-type]
        )


@pytest.mark.parametrize("prandtl_number", [0.0, np.nan])
def test_invalid_prandtl_number(prandtl_number: float) -> None:
    with pytest.raises(ValueError, match="prandtl_number"):
        flat_plate_boundary_layer(
            1.0,
            10.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.LAMINAR,
            prandtl_number=prandtl_number,
        )


def test_van_driest_noncompressible_limit_is_stable() -> None:
    reynolds = np.asarray(1e6)
    state = _van_driest_ii_state(
        reynolds,
        np.asarray(250.0),
        np.asarray(0.0),
        None,
        prandtl_number=0.72,
        gas=AIR,
        viscosity_model=AIR_VISCOSITY,
    )
    assert state.friction_factor == 1.0
    assert state.momentum_factor == 1.0
    assert state.reynolds_factor == 1.0
    result = flat_plate_boundary_layer(
        1.0,
        10.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TURBULENT,
        compressibility_correction=CompressibilityCorrection.VAN_DRIEST_II,
        mach=0.0,
        edge_temperature=250.0,
    )
    assert result.effective_reynolds_number == result.reynolds_number
    assert result.recovery_temperature == 250.0
    assert result.wall_temperature == 250.0
    local = float(result.local_skin_friction_coefficient)
    residual = (
        0.242 / np.sqrt(local)
        - 0.41
        - np.log10(float(result.reynolds_number) * local)
    )
    assert residual == pytest.approx(0.0, abs=2e-14)


def test_van_driest_array_inputs_broadcast() -> None:
    result = flat_plate_boundary_layer(
        [[0.5], [1.0]],
        100.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TURBULENT,
        compressibility_correction=CompressibilityCorrection.VAN_DRIEST_II,
        mach=[2.0, 5.0, 8.0],
        edge_temperature=220.0,
        wall_temperature=300.0,
    )
    for value in (
        result.effective_reynolds_number,
        result.local_skin_friction_coefficient,
        result.average_skin_friction_coefficient,
        result.drag_per_unit_width,
    ):
        assert isinstance(value, np.ndarray)
        assert value.shape == (2, 3)
        assert value.dtype == np.float64


def test_van_driest_selection_uses_eckert_for_laminar_flow() -> None:
    result = flat_plate_boundary_layer(
        1.0,
        10.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.LAMINAR,
        compressibility_correction=CompressibilityCorrection.VAN_DRIEST_II,
        mach=0.0,
        edge_temperature=250.0,
    )
    assert result.effective_reynolds_number == result.reynolds_number
    assert result.wall_temperature == 250.0


def test_van_driest_rejects_degenerate_thermal_state() -> None:
    with pytest.raises(ModelRangeError, match="undefined"):
        flat_plate_boundary_layer(
            1.0,
            100.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.TURBULENT,
            compressibility_correction=CompressibilityCorrection.VAN_DRIEST_II,
            mach=2.0,
            edge_temperature=250.0,
            wall_temperature=1e100,
        )


def test_van_driest_reports_unbracketed_willems_solution() -> None:
    with pytest.raises(ModelRangeError, match="bracket"):
        flat_plate_boundary_layer(
            0.1,
            1.0,
            1.0,
            1.0,
            regime=BoundaryLayerRegime.TURBULENT,
            compressibility_correction=CompressibilityCorrection.VAN_DRIEST_II,
            mach=0.0,
            edge_temperature=250.0,
        )


def test_thermal_inputs_must_broadcast() -> None:
    with pytest.raises(ValueError, match="broadcastable"):
        flat_plate_boundary_layer(
            [1.0, 2.0],
            10.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.LAMINAR,
            compressibility_correction=CompressibilityCorrection.ECKERT,
            mach=[1.0, 2.0, 3.0],
            edge_temperature=250.0,
        )
