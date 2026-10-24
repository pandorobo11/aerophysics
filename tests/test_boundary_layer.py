"""Tests for smooth flat-plate boundary-layer correlations."""

from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    CompressibilityCorrection,
    TurbulentCorrelation,
    flat_plate_boundary_layer,
)
from aerophysics.exceptions import ApplicabilityWarning, ModelRangeError
from aerophysics.gas import AIR, AIR_VISCOSITY
from aerophysics.transport import AIR_BLOTTNER_VISCOSITY


def test_blasius_dimensional_outputs() -> None:
    # Nondimensional Blasius values are covered by the source-equation verification.
    result = flat_plate_boundary_layer(
        1.0,
        10.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.LAMINAR,
    )
    assert result.reynolds_number == pytest.approx(1e6)
    assert result.wall_shear_stress == pytest.approx(0.0332)
    assert result.drag_per_unit_width == pytest.approx(0.0664)
    assert result.effective_reynolds_number == result.reynolds_number
    assert result.recovery_temperature is None
    assert result.wall_temperature is None


def test_one_fifth_power_turbulent_integral_thicknesses() -> None:
    result = flat_plate_boundary_layer(
        1.0,
        10.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TURBULENT,
        turbulent_correlation=TurbulentCorrelation.POWER_LAW,
    )
    # delta and C_f are covered by the source-equation verification.
    delta = float(result.boundary_layer_thickness)
    assert result.displacement_thickness == pytest.approx(delta / 8.0)
    assert result.momentum_thickness == pytest.approx(delta * 7.0 / 72.0)


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


def test_eckert_accepts_blottner_mixture_viscosity_model() -> None:
    edge_temperature = 1000.0
    edge_viscosity = float(AIR_BLOTTNER_VISCOSITY.dynamic_viscosity(edge_temperature))
    result = flat_plate_boundary_layer(
        1.0,
        1000.0,
        0.1,
        edge_viscosity,
        regime=BoundaryLayerRegime.LAMINAR,
        compressibility_correction=CompressibilityCorrection.ECKERT,
        mach=1.0,
        edge_temperature=edge_temperature,
        wall_temperature=1000.0,
        viscosity_model=AIR_BLOTTNER_VISCOSITY,
    )
    sutherland = flat_plate_boundary_layer(
        1.0,
        1000.0,
        0.1,
        edge_viscosity,
        regime=BoundaryLayerRegime.LAMINAR,
        compressibility_correction=CompressibilityCorrection.ECKERT,
        mach=1.0,
        edge_temperature=edge_temperature,
        wall_temperature=1000.0,
    )
    assert result.effective_reynolds_number != pytest.approx(
        sutherland.effective_reynolds_number, rel=1e-3
    )


def test_eckert_uses_supplied_edge_viscosity_in_correction() -> None:
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
        * 1e-5
        / float(AIR_VISCOSITY.dynamic_viscosity(reference))
    )
    assert result.recovery_temperature == pytest.approx(recovery)
    assert result.wall_temperature == 300.0
    assert result.effective_reynolds_number == pytest.approx(expected_reynolds)
    assert result.average_skin_friction_coefficient == pytest.approx(
        220.0 / reference * 0.074 * expected_reynolds**-0.2
    )


@pytest.mark.parametrize(
    ("correction", "regime", "transition_reynolds"),
    [
        (CompressibilityCorrection.ECKERT, BoundaryLayerRegime.LAMINAR, None),
        (CompressibilityCorrection.ECKERT, BoundaryLayerRegime.TURBULENT, None),
        (CompressibilityCorrection.ECKERT, BoundaryLayerRegime.TRANSITIONAL, 2e6),
        (CompressibilityCorrection.VAN_DRIEST_II, BoundaryLayerRegime.LAMINAR, None),
    ],
)
def test_eckert_reduces_to_none_in_incompressible_limit(
    correction: CompressibilityCorrection,
    regime: BoundaryLayerRegime,
    transition_reynolds: float | None,
) -> None:
    # M = 0 and T_w = T_e give T* = T_r = T_e and mu* = mu_e, so the
    # reference-density conversion T_e/T* and the Reynolds scaling are unity.
    # Van Driest II laminar portions use Eckert; its turbulent portions use
    # the Willems equations rather than the selected correlation.
    edge_temperature = 250.0
    viscosity = float(AIR_VISCOSITY.dynamic_viscosity(edge_temperature))
    arguments: dict[str, Any] = {
        "distance": [0.05, 0.2, 0.8],
        "edge_velocity": 120.0,
        "edge_density": 1.2,
        "edge_dynamic_viscosity": viscosity,
        "regime": regime,
        "transition_reynolds": transition_reynolds,
    }
    reference = flat_plate_boundary_layer(**arguments)
    corrected = flat_plate_boundary_layer(
        **arguments,
        compressibility_correction=correction,
        mach=0.0,
        edge_temperature=edge_temperature,
        wall_temperature=edge_temperature,
    )
    for name in (
        "effective_reynolds_number",
        "boundary_layer_thickness",
        "displacement_thickness",
        "momentum_thickness",
        "local_skin_friction_coefficient",
        "average_skin_friction_coefficient",
        "wall_shear_stress",
        "drag_per_unit_width",
    ):
        assert_allclose(
            getattr(corrected, name),
            getattr(reference, name),
            rtol=1e-12,
            err_msg=name,
        )


def _hypersonic_cold_wall_state() -> tuple[float, float, float]:
    """Return U_e, rho_e and mu_e for M=6, T_e=240 K and Re'=2e7 1/m."""
    edge_viscosity = float(AIR_VISCOSITY.dynamic_viscosity(240.0))
    edge_velocity = 6.0 * np.sqrt(
        AIR.heat_capacity_ratio * AIR.specific_gas_constant * 240.0
    )
    edge_density = 2e7 * edge_viscosity / edge_velocity
    return float(edge_velocity), float(edge_density), edge_viscosity


def test_eckert_hypersonic_cold_wall_uses_reference_dynamic_pressure() -> None:
    edge_velocity, edge_density, edge_viscosity = _hypersonic_cold_wall_state()
    distance = np.array([0.05, 0.5, 1.0])
    result = flat_plate_boundary_layer(
        distance,
        edge_velocity,
        edge_density,
        edge_viscosity,
        regime=BoundaryLayerRegime.TURBULENT,
        compressibility_correction=CompressibilityCorrection.ECKERT,
        mach=6.0,
        edge_temperature=240.0,
        wall_temperature=300.0,
    )

    recovery = 240.0 * (1.0 + np.cbrt(0.72) * 0.2 * 36.0)
    reference_temperature = 0.22 * recovery + 0.28 * 240.0 + 0.50 * 300.0
    assert reference_temperature == pytest.approx(610.73, abs=0.01)
    reference_density = edge_density * 240.0 / reference_temperature
    reference_viscosity = float(AIR_VISCOSITY.dynamic_viscosity(reference_temperature))
    reference_reynolds = (
        reference_density * edge_velocity * distance / reference_viscosity
    )
    average_reference = 0.455 / np.log10(reference_reynolds) ** 2.58
    local_reference = average_reference * (1.0 - 2.58 / np.log(reference_reynolds))

    assert_allclose(result.effective_reynolds_number, reference_reynolds, rtol=1e-12)
    assert_allclose(
        result.wall_shear_stress,
        0.5 * reference_density * edge_velocity**2 * local_reference,
        rtol=1e-12,
    )
    assert_allclose(
        result.drag_per_unit_width,
        0.5 * reference_density * edge_velocity**2 * distance * average_reference,
        rtol=1e-12,
    )
    edge_dynamic_pressure = 0.5 * edge_density * edge_velocity**2
    assert_allclose(
        result.local_skin_friction_coefficient,
        np.asarray(result.wall_shear_stress) / edge_dynamic_pressure,
        rtol=1e-12,
    )
    # Regression values for Re_x = 1e6, 1e7 and 2e7.
    assert_allclose(
        result.wall_shear_stress, [550.21382, 367.23657, 328.78010], rtol=1e-6
    )


def test_eckert_hypersonic_cold_wall_laminar_uses_reference_dynamic_pressure() -> None:
    edge_velocity, edge_density, edge_viscosity = _hypersonic_cold_wall_state()
    distance = np.array([0.05, 0.5])
    arguments: dict[str, Any] = {
        "regime": BoundaryLayerRegime.LAMINAR,
        "mach": 6.0,
        "edge_temperature": 240.0,
        "wall_temperature": 300.0,
    }
    eckert = flat_plate_boundary_layer(
        distance,
        edge_velocity,
        edge_density,
        edge_viscosity,
        compressibility_correction=CompressibilityCorrection.ECKERT,
        **arguments,
    )
    van_driest = flat_plate_boundary_layer(
        distance,
        edge_velocity,
        edge_density,
        edge_viscosity,
        compressibility_correction=CompressibilityCorrection.VAN_DRIEST_II,
        **arguments,
    )

    recovery = 240.0 * (1.0 + np.sqrt(0.72) * 0.2 * 36.0)
    reference_temperature = 0.22 * recovery + 0.28 * 240.0 + 0.50 * 300.0
    reference_density = edge_density * 240.0 / reference_temperature
    reference_viscosity = float(AIR_VISCOSITY.dynamic_viscosity(reference_temperature))
    reference_reynolds = (
        reference_density * edge_velocity * distance / reference_viscosity
    )
    expected_shear = (
        0.5 * reference_density * edge_velocity**2 * 0.664 / np.sqrt(reference_reynolds)
    )
    # Chapman--Rubesin form: C_f sqrt(Re_x) = 0.664 sqrt(rho* mu* / rho_e mu_e).
    chapman_rubesin = (
        reference_density * reference_viscosity / (edge_density * edge_viscosity)
    )
    edge_reynolds = edge_density * edge_velocity * distance / edge_viscosity

    for result in (eckert, van_driest):
        assert_allclose(result.wall_shear_stress, expected_shear, rtol=1e-12)
        assert_allclose(
            result.local_skin_friction_coefficient,
            0.664 * np.sqrt(chapman_rubesin / edge_reynolds),
            rtol=1e-12,
        )
        assert_allclose(result.drag_per_unit_width, 2.0 * distance * expected_shear)
    # The thickness path is unchanged: Blasius lengths evaluated at Re*.
    assert_allclose(
        eckert.boundary_layer_thickness,
        5.0 * distance / np.sqrt(reference_reynolds),
        rtol=1e-12,
    )


def test_eckert_transition_mixes_edge_basis_average_friction() -> None:
    edge_velocity, edge_density, edge_viscosity = _hypersonic_cold_wall_state()
    transition = 2e6
    x_transition = transition / 2e7
    distances = np.array([x_transition, x_transition * (1.0 + 1e-9), 1.0])
    result = flat_plate_boundary_layer(
        distances,
        edge_velocity,
        edge_density,
        edge_viscosity,
        regime=BoundaryLayerRegime.TRANSITIONAL,
        transition_reynolds=transition,
        turbulent_correlation=TurbulentCorrelation.POWER_LAW,
        compressibility_correction=CompressibilityCorrection.ECKERT,
        mach=6.0,
        edge_temperature=240.0,
        wall_temperature=300.0,
    )

    def reference_state(recovery_factor: float) -> tuple[float, float]:
        recovery = 240.0 * (1.0 + recovery_factor * 0.2 * 36.0)
        temperature = 0.22 * recovery + 0.28 * 240.0 + 0.50 * 300.0
        reynolds_scale = (
            240.0
            / temperature
            * edge_viscosity
            / float(AIR_VISCOSITY.dynamic_viscosity(temperature))
        )
        return 240.0 / temperature, reynolds_scale

    laminar_ratio, laminar_scale = reference_state(np.sqrt(0.72))
    turbulent_ratio, turbulent_scale = reference_state(np.cbrt(0.72))
    laminar_average = laminar_ratio * 1.328 / np.sqrt(laminar_scale * transition)
    expected_average = turbulent_ratio * 0.074 * (turbulent_scale * 2e7) ** -0.2 + (
        transition / 2e7
    ) * (
        laminar_average
        - turbulent_ratio * 0.074 * (turbulent_scale * transition) ** -0.2
    )
    average = np.asarray(result.average_skin_friction_coefficient)
    drag = np.asarray(result.drag_per_unit_width)
    assert average[2] == pytest.approx(expected_average, rel=1e-12)
    # Accumulated drag is continuous: across the infinitesimal step it grows
    # only by the local turbulent wall shear times the step length.  The
    # power-law local coefficient is the exact derivative of its average.
    shear = np.asarray(result.wall_shear_stress)
    assert drag[1] - drag[0] == pytest.approx(
        shear[1] * (distances[1] - distances[0]), rel=1e-5
    )


def _van_driest_factors(
    mach: float,
    edge_temperature: float,
    wall_temperature: float | None = None,
    edge_viscosity: float = 1e-5,
) -> tuple[float, float, float, float]:
    recovery_factor = np.cbrt(0.72)
    m = recovery_factor * (AIR.heat_capacity_ratio - 1.0) * mach**2 / 2.0
    recovery = edge_temperature * (1.0 + m)
    wall = recovery if wall_temperature is None else wall_temperature
    if m <= 1e-12:
        return 1.0, 1.0, 1.0, recovery
    temperature_factor = wall / edge_temperature
    denominator = np.sqrt(
        (m + 1.0 + temperature_factor) ** 2 - 4.0 * temperature_factor
    )
    alpha = np.clip((m - 1.0 + temperature_factor) / denominator, -1.0, 1.0)
    beta = np.clip((m + 1.0 - temperature_factor) / denominator, -1.0, 1.0)
    friction_factor = m / (np.arcsin(alpha) + np.arcsin(beta)) ** 2
    momentum_factor = edge_viscosity / float(AIR_VISCOSITY.dynamic_viscosity(wall))
    return (
        friction_factor,
        momentum_factor,
        momentum_factor / friction_factor,
        recovery,
    )


@pytest.mark.parametrize(
    ("mach", "wall_temperature"),
    [(2.0, None), (5.0, 300.0)],
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
    friction_factor, momentum_factor, reynolds_factor, recovery = _van_driest_factors(
        mach, edge_temperature, wall_temperature
    )
    expected_wall = recovery if wall_temperature is None else wall_temperature
    assert result.wall_temperature == pytest.approx(expected_wall)
    assert result.effective_reynolds_number == pytest.approx(
        result.reynolds_number * reynolds_factor
    )
    # Thickness uses Re_x*F_theta, not the friction Reynolds number.
    assert result.boundary_layer_thickness == pytest.approx(
        0.37 * (float(result.reynolds_number) * momentum_factor) ** -0.2
    )

    local = float(result.local_skin_friction_coefficient)
    average = float(result.average_skin_friction_coefficient)
    local_residual = 0.242 / np.sqrt(local * friction_factor) - (
        0.41
        + np.log10(
            float(result.reynolds_number) * reynolds_factor * local * friction_factor
        )
    )
    average_residual = 0.242 / np.sqrt(average * friction_factor) - np.log10(
        float(result.reynolds_number) * reynolds_factor * average * friction_factor
    )
    assert local_residual == pytest.approx(0.0, abs=2e-14)
    assert average_residual == pytest.approx(0.0, abs=2e-14)
    assert result.drag_per_unit_width == pytest.approx(
        0.5 * 1.0 * 100.0**2 * 1.0 * average
    )


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
    # Accumulated drag is continuous: across the infinitesimal step it grows
    # only by about the local turbulent wall shear times the step length.  The
    # Willems local and average equations are separate fits, so the local
    # coefficient is not the exact derivative of the average.
    shear = np.asarray(result.wall_shear_stress)
    assert drag[1] - drag[0] == pytest.approx(
        shear[1] * (distances[1] - distances[0]), rel=0.1
    )


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
    # Without an explicit wall temperature the wall is adiabatic.
    assert result.recovery_temperature is not None
    assert_allclose(result.wall_temperature, np.asarray(result.recovery_temperature))
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
        0.242 / np.sqrt(local) - 0.41 - np.log10(float(result.reynolds_number) * local)
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
