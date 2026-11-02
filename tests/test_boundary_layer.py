"""Tests for smooth flat-plate boundary-layer correlations."""

import warnings
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.optimize import brentq

from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    TurbulentModel,
    flat_plate_boundary_layer,
)
from aerophysics.exceptions import ApplicabilityWarning, ModelRangeError
from aerophysics.gas import AIR, AIR_VISCOSITY
from aerophysics.transport import AIR_BLOTTNER_VISCOSITY, SutherlandModel


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
        turbulent_model=TurbulentModel.SCHLICHTING,
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


def _schlichting_average(reynolds: Any) -> Any:
    return 0.455 / np.log10(reynolds) ** 2.58


def _schlichting_local(reynolds: Any) -> Any:
    return _schlichting_average(reynolds) * (1.0 - 2.58 / np.log(reynolds))


def _schlichting_start(drag_reynolds: float) -> float:
    """Return R with fully turbulent Schlichting R * Cbar_f = drag_reynolds."""
    return float(
        brentq(
            lambda value: value * _schlichting_average(value) - drag_reynolds,
            20.0,
            1e15,
            rtol=1e-15,
        )
    )


def test_transition_continues_equivalent_turbulent_plate() -> None:
    transition = 5e5
    reynolds = np.array([2.5e5, transition, transition * (1.0 + 1e-9), 1e6, 1e7])
    x = reynolds * 1e-5 / 10.0
    arguments: dict[str, Any] = {
        "edge_velocity": 10.0,
        "edge_density": 1.0,
        "edge_dynamic_viscosity": 1e-5,
        "regime": BoundaryLayerRegime.TRANSITIONAL,
        "transition_reynolds": transition,
    }
    result = flat_plate_boundary_layer(
        x, turbulent_model=TurbulentModel.SCHLICHTING, **arguments
    )
    # The Schlichting plate with the Blasius drag at transition starts
    # Re_x - Re_tr upstream of the equivalent fully turbulent state.
    start = _schlichting_start(1.328 * np.sqrt(transition))
    assert start == pytest.approx(1.41726e5, rel=1e-5)
    equivalent = start + reynolds[2:] - transition
    local = np.asarray(result.local_skin_friction_coefficient)
    average = np.asarray(result.average_skin_friction_coefficient)
    assert local[0] == pytest.approx(0.664 / np.sqrt(reynolds[0]))
    assert average[1] == pytest.approx(1.328 / np.sqrt(transition))
    assert_allclose(local[2:], _schlichting_local(equivalent), rtol=1e-12)
    assert_allclose(
        average[2:],
        _schlichting_average(equivalent) * equivalent / reynolds[2:],
        rtol=1e-12,
    )
    assert_allclose(
        np.asarray(result.effective_reynolds_number),
        np.concatenate((reynolds[:2], equivalent)),
        rtol=1e-12,
    )
    # Accumulated drag is continuous and grows by the local turbulent shear.
    drag = np.asarray(result.drag_per_unit_width)
    shear = np.asarray(result.wall_shear_stress)
    assert drag[2] - drag[1] == pytest.approx(shear[2] * (x[2] - x[1]), rel=1e-5)
    # Downstream, the transitional plate approaches the fully turbulent one.
    turbulent = flat_plate_boundary_layer(
        x[-1],
        regime=BoundaryLayerRegime.TURBULENT,
        edge_velocity=10.0,
        edge_density=1.0,
        edge_dynamic_viscosity=1e-5,
    )
    assert average[-1] < float(turbulent.average_skin_friction_coefficient)
    assert average[-1] == pytest.approx(
        turbulent.average_skin_friction_coefficient, rel=0.04
    )


def test_power_law_thickness_uses_transition_virtual_origin() -> None:
    transition = 5e5
    reynolds = np.array([transition * (1.0 + 1e-12), 1e6, 5e6])
    x = reynolds * 1e-5 / 10.0
    result = flat_plate_boundary_layer(
        x,
        10.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TRANSITIONAL,
        transition_reynolds=transition,
    )
    # The friction-law origin also places the one-fifth-power thicknesses.
    effective = _schlichting_start(1.328 * np.sqrt(transition)) + reynolds - transition
    delta = 0.37 * (x * effective / reynolds) * effective**-0.2
    assert_allclose(result.boundary_layer_thickness, delta, rtol=1e-12)
    assert_allclose(result.displacement_thickness, delta / 8.0, rtol=1e-12)
    assert_allclose(result.momentum_thickness, delta * 7.0 / 72.0, rtol=1e-12)
    # The drag-based theta is exactly continuous with Blasius. The returned
    # power-law theta is not: it carries the one-fifth-power law's own
    # mismatch with Schlichting friction, +1.268% at Re_tr = 5e5.
    laminar_theta = 0.664 * x[0] / np.sqrt(transition)
    assert x[0] * np.asarray(result.average_skin_friction_coefficient)[
        0
    ] / 2.0 == pytest.approx(laminar_theta, rel=1e-9)
    assert np.asarray(result.momentum_thickness)[0] / laminar_theta - 1.0 == (
        pytest.approx(0.012684, abs=1e-6)
    )
    turbulent = flat_plate_boundary_layer(
        x[-1], 10.0, 1.0, 1e-5, regime=BoundaryLayerRegime.TURBULENT
    )
    assert float(turbulent.boundary_layer_thickness) == pytest.approx(
        0.37 * x[-1] * reynolds[-1] ** -0.2, rel=1e-12
    )
    assert np.asarray(result.boundary_layer_thickness)[-1] < float(
        turbulent.boundary_layer_thickness
    )


def test_transition_below_turbulent_correlation_domain_raises() -> None:
    # Schlichting R * Cbar_f has its minimum about 4.5, so a transition
    # momentum thickness below that has no equivalent turbulent plate.
    with pytest.raises(ModelRangeError, match="increase transition_reynolds"):
        flat_plate_boundary_layer(
            1.0,
            10.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.TRANSITIONAL,
            transition_reynolds=10.0,
        )


def _thickness_range_warnings(**arguments: Any) -> list[str]:
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        # Re_x = 2e7 lies above the one-fifth-power thickness range.
        flat_plate_boundary_layer(1.0, 200.0, 1.0, 1e-5, **arguments)
    return [str(item.message) for item in captured if "one-fifth" in str(item.message)]


def test_power_law_thickness_warns_above_its_validity_range() -> None:
    turbulent = BoundaryLayerRegime.TURBULENT
    assert _thickness_range_warnings(regime=turbulent)
    assert not _thickness_range_warnings(
        regime=turbulent, turbulent_model=TurbulentModel.SCHLICHTING_LOG_WAKE
    )
    assert not _thickness_range_warnings(
        regime=turbulent,
        compressible=True,
        turbulent_model=TurbulentModel.ECKERT,
        mach=0.3,
        edge_temperature=288.0,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        flat_plate_boundary_layer(
            1.0, 10.0, 1.0, 1e-5, regime=BoundaryLayerRegime.TURBULENT
        )


def test_log_wake_thickness_closes_the_momentum_integral() -> None:
    from aerophysics.boundary_layer_profile import (
        CompressibleVelocityTransformation,
        TemperatureVelocityRelation,
        compressible_turbulent_boundary_layer_profile,
    )

    temperature = 288.15
    viscosity = float(AIR_VISCOSITY.dynamic_viscosity(temperature))
    velocity, density, x = 10.0, 1.225, 10.0
    result = flat_plate_boundary_layer(
        x,
        velocity,
        density,
        viscosity,
        regime=BoundaryLayerRegime.TURBULENT,
        turbulent_model=TurbulentModel.SCHLICHTING_LOG_WAKE,
    )
    power = flat_plate_boundary_layer(
        x, velocity, density, viscosity, regime=BoundaryLayerRegime.TURBULENT
    )
    for field in (
        "effective_reynolds_number",
        "local_skin_friction_coefficient",
        "average_skin_friction_coefficient",
        "wall_shear_stress",
        "drag_per_unit_width",
    ):
        assert getattr(result, field) == getattr(power, field)
    assert result.momentum_thickness == pytest.approx(
        x * float(result.average_skin_friction_coefficient) / 2.0, rel=1e-14
    )
    delta = float(result.boundary_layer_thickness)
    # Mach 0.03 makes the public compressible profile almost incompressible.
    public = compressible_turbulent_boundary_layer_profile(
        np.array([0.0, delta]),
        velocity,
        density,
        temperature,
        delta,
        float(result.wall_shear_stress),
        transformation=CompressibleVelocityTransformation.VAN_DRIEST,
        temperature_velocity_relation=TemperatureVelocityRelation.WALZ,
        wall_temperature=temperature,
    )
    assert 0.0 <= public.wake_parameter <= 1.0
    assert public.momentum_thickness == pytest.approx(
        0.978 * float(result.momentum_thickness), rel=1e-4
    )
    assert result.displacement_thickness == pytest.approx(
        public.displacement_thickness, rel=1e-4
    )
    compressible = flat_plate_boundary_layer(
        x,
        velocity,
        density,
        viscosity,
        regime=BoundaryLayerRegime.TURBULENT,
        compressible=True,
        turbulent_model=TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE,
        mach=velocity / float(AIR.speed_of_sound(temperature)),
        edge_temperature=temperature,
    )
    assert result.boundary_layer_thickness == pytest.approx(
        compressible.boundary_layer_thickness, rel=1e-3
    )


def test_log_wake_transition_is_admissible_after_late_transition() -> None:
    # After transition (cf, theta) are a fully turbulent state, so the wake
    # closure is admissible immediately downstream of a late transition.
    transition = 3e6
    reynolds = np.array([transition * 0.5, transition * 1.01, transition * 3.0])
    x = reynolds * 1e-5 / 10.0
    with warnings.catch_warnings():
        warnings.simplefilter("error", ApplicabilityWarning)
        result = flat_plate_boundary_layer(
            x,
            10.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.TRANSITIONAL,
            transition_reynolds=transition,
            turbulent_model=TurbulentModel.SCHLICHTING_LOG_WAKE,
        )
    thickness = np.asarray(result.boundary_layer_thickness)
    assert thickness[0] == pytest.approx(4.91 * x[0] / np.sqrt(reynolds[0]))
    assert np.all(np.isfinite(thickness))
    assert np.all(np.isfinite(result.displacement_thickness))
    assert_allclose(
        np.asarray(result.momentum_thickness)[1:],
        x[1:] * np.asarray(result.average_skin_friction_coefficient)[1:] / 2.0,
        rtol=1e-14,
    )
    # The turbulent state equals the fully turbulent plate at Re_eff.
    effective = float(np.asarray(result.effective_reynolds_number)[1])
    turbulent = flat_plate_boundary_layer(
        x[1] * effective / reynolds[1],
        10.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TURBULENT,
        turbulent_model=TurbulentModel.SCHLICHTING_LOG_WAKE,
    )
    assert np.asarray(result.wall_shear_stress)[1] == pytest.approx(
        turbulent.wall_shear_stress, rel=1e-12
    )
    assert thickness[1] == pytest.approx(turbulent.boundary_layer_thickness, rel=1e-5)


def test_log_wake_returns_nan_where_wake_is_inadmissible() -> None:
    reynolds = np.array([1e6, 1e14])
    x = reynolds * 1e-5 / 10.0
    with pytest.warns(ApplicabilityWarning, match="at 1 turbulent point"):
        result = flat_plate_boundary_layer(
            x,
            10.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.TURBULENT,
            turbulent_model=TurbulentModel.SCHLICHTING_LOG_WAKE,
        )
    thickness = np.asarray(result.boundary_layer_thickness)
    displacement = np.asarray(result.displacement_thickness)
    assert np.isfinite(thickness[0]) and np.isfinite(displacement[0])
    assert np.isnan(thickness[1]) and np.isnan(displacement[1])
    assert np.isfinite(np.asarray(result.momentum_thickness)[1])


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
        compressible=True,
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
        compressible=True,
        mach=1.0,
        edge_temperature=edge_temperature,
        wall_temperature=1000.0,
    )
    assert result.effective_reynolds_number == result.reynolds_number
    assert result.local_skin_friction_coefficient != pytest.approx(
        sutherland.local_skin_friction_coefficient, rel=1e-4
    )


def test_eckert_uses_supplied_edge_viscosity_in_correction() -> None:
    result = flat_plate_boundary_layer(
        1.0,
        100.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TURBULENT,
        compressible=True,
        turbulent_model=TurbulentModel.ECKERT,
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
        220.0 / reference * 0.455 / np.log10(expected_reynolds) ** 2.58
    )


@pytest.mark.parametrize(
    ("correction", "regime", "transition_reynolds"),
    [
        (TurbulentModel.ECKERT, BoundaryLayerRegime.LAMINAR, None),
        (TurbulentModel.ECKERT, BoundaryLayerRegime.TURBULENT, None),
        (TurbulentModel.ECKERT, BoundaryLayerRegime.TRANSITIONAL, 2e6),
        (TurbulentModel.VAN_DRIEST_II, BoundaryLayerRegime.LAMINAR, None),
    ],
)
def test_eckert_reduces_to_none_in_incompressible_limit(
    correction: TurbulentModel,
    regime: BoundaryLayerRegime,
    transition_reynolds: float | None,
) -> None:
    # M = 0 and T_w = T_e give T* = T_r = T_e and mu* = mu_e, so the
    # reference-density conversion T_e/T* and the Reynolds scaling are unity.
    # Compressible laminar portions approach the unrounded Blasius solution.
    # VD2 turbulent portions use Willems rather than the selected correlation.
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
        compressible=True,
        turbulent_model=None if regime is BoundaryLayerRegime.LAMINAR else correction,
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
        expected = np.asarray(getattr(reference, name))
        if name in (
            "boundary_layer_thickness",
            "displacement_thickness",
            "momentum_thickness",
        ):
            turbulent = np.zeros_like(expected, dtype=bool)
            if regime is BoundaryLayerRegime.TURBULENT:
                turbulent[:] = True
            elif regime is BoundaryLayerRegime.TRANSITIONAL:
                turbulent = np.asarray(reference.reynolds_number) > transition_reynolds
            expected = np.where(turbulent, np.nan, expected)
        assert_allclose(
            getattr(corrected, name),
            expected,
            # Transitional turbulent points depend on the unrounded laminar
            # momentum thickness at transition.
            rtol=(
                1e-12
                if regime is BoundaryLayerRegime.TURBULENT
                or (
                    name == "effective_reynolds_number"
                    and regime is BoundaryLayerRegime.LAMINAR
                )
                else 2e-4
            ),
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
        compressible=True,
        turbulent_model=TurbulentModel.ECKERT,
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


def test_hypersonic_cold_wall_laminar_uses_similarity_edge_dynamic_pressure() -> None:
    velocity, density, viscosity = _hypersonic_cold_wall_state()
    distance = np.array([0.05, 0.5])
    result = flat_plate_boundary_layer(
        distance,
        velocity,
        density,
        viscosity,
        regime=BoundaryLayerRegime.LAMINAR,
        compressible=True,
        mach=6.0,
        edge_temperature=240.0,
        wall_temperature=300.0,
    )
    assert_allclose(result.effective_reynolds_number, result.reynolds_number)
    assert_allclose(
        result.wall_shear_stress,
        0.5
        * density
        * velocity**2
        * np.asarray(result.local_skin_friction_coefficient),
    )
    assert_allclose(
        result.drag_per_unit_width,
        2 * distance * np.asarray(result.wall_shear_stress),
    )
    assert_allclose(
        result.momentum_thickness,
        distance * np.asarray(result.local_skin_friction_coefficient),
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
        compressible=True,
        turbulent_model=TurbulentModel.ECKERT,
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

    turbulent_ratio, turbulent_scale = reference_state(np.cbrt(0.72))
    laminar_average = float(
        flat_plate_boundary_layer(
            x_transition,
            edge_velocity,
            edge_density,
            edge_viscosity,
            regime=BoundaryLayerRegime.LAMINAR,
            compressible=True,
            mach=6.0,
            edge_temperature=240.0,
            wall_temperature=300.0,
        ).average_skin_friction_coefficient
    )
    # Edge-basis drag Re * Cbar_f = ratio * Re * Cbar_i(scale * Re) of the
    # equivalent fully turbulent plate equals the laminar value at transition.
    start = (
        _schlichting_start(
            transition * laminar_average * turbulent_scale / turbulent_ratio
        )
        / turbulent_scale
    )
    equivalent = start + 2e7 - transition
    expected_average = (
        turbulent_ratio
        * _schlichting_average(equivalent * turbulent_scale)
        * equivalent
        / 2e7
    )
    average = np.asarray(result.average_skin_friction_coefficient)
    drag = np.asarray(result.drag_per_unit_width)
    assert average[2] == pytest.approx(expected_average, rel=1e-12)
    # Accumulated drag is continuous: across the infinitesimal step it grows
    # only by the local turbulent wall shear times the step length.  The
    # Schlichting local coefficient is the exact derivative of its average.
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
        compressible=True,
        turbulent_model=TurbulentModel.VAN_DRIEST_II,
        mach=mach,
        edge_temperature=edge_temperature,
        wall_temperature=wall_temperature,
    )
    friction_factor, _, reynolds_factor, recovery = _van_driest_factors(
        mach, edge_temperature, wall_temperature
    )
    expected_wall = recovery if wall_temperature is None else wall_temperature
    assert result.wall_temperature == pytest.approx(expected_wall)
    assert result.effective_reynolds_number == pytest.approx(
        result.reynolds_number * reynolds_factor
    )
    assert np.isnan(result.boundary_layer_thickness)

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


def test_van_driest_transition_uses_eckert_then_van_driest() -> None:
    mixed = flat_plate_boundary_layer(
        [0.1, 1.0],
        100.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TRANSITIONAL,
        transition_reynolds=2e6,
        compressible=True,
        turbulent_model=TurbulentModel.VAN_DRIEST_II,
        mach=2.0,
        edge_temperature=250.0,
    )
    laminar = flat_plate_boundary_layer(
        0.1,
        100.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.LAMINAR,
        compressible=True,
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
        compressible=True,
        turbulent_model=TurbulentModel.VAN_DRIEST_II,
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
        compressible=True,
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
            compressible=True,
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
    with pytest.raises(ValueError, match="turbulent_model"):
        flat_plate_boundary_layer(
            1.0,
            10.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.LAMINAR,
            turbulent_model="power_law",  # type: ignore[arg-type]
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
            compressible=True,
            mach=1.0,
        )
    with pytest.raises(ValueError, match="compressible"):
        flat_plate_boundary_layer(
            1.0,
            10.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.LAMINAR,
            compressible="eckert",  # type: ignore[arg-type]
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
            compressible=True,
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
        compressible=True,
        turbulent_model=TurbulentModel.VAN_DRIEST_II,
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
        compressible=True,
        turbulent_model=TurbulentModel.VAN_DRIEST_II,
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
            compressible=True,
            turbulent_model=TurbulentModel.VAN_DRIEST_II,
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
            compressible=True,
            turbulent_model=TurbulentModel.VAN_DRIEST_II,
            mach=0.0,
            edge_temperature=250.0,
        )


# The independent fixture uses explicit local edge states and wall temperatures,
# including AD rows. SU2 delta99 is deliberately not a test target.
def _wall_temperature_reference_rows() -> list[dict[str, str]]:
    import csv
    from pathlib import Path

    path = (
        Path(__file__).parent
        / "reference_data/viscous_flow/van_driest_ii_wall_temperature_reference.csv"
    )
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


@pytest.mark.parametrize(
    "row", _wall_temperature_reference_rows(), ids=lambda row: row["case"]
)
def test_wall_temperature_reference_and_public_profile(row: dict[str, str]) -> None:
    from aerophysics.boundary_layer_profile import (
        CompressibleVelocityTransformation,
        TemperatureVelocityRelation,
        _van_driest_walz_profile,
        compressible_turbulent_boundary_layer_profile,
    )
    from aerophysics.gas import PerfectGas
    from aerophysics.transport import SutherlandModel

    r = {key: float(value) for key, value in row.items() if key != "case"}
    gas = PerfectGas(287.053, 1.4)
    viscosity = SutherlandModel(1.7894e-5, 288.15, 110.4)
    result = flat_plate_boundary_layer(
        r["x_m"],
        r["edge_velocity"],
        r["edge_density"],
        r["edge_viscosity"],
        regime=BoundaryLayerRegime.TURBULENT,
        compressible=True,
        turbulent_model=TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE,
        mach=r["mach"],
        edge_temperature=r["edge_temperature"],
        wall_temperature=r["wall_temperature"],
        prandtl_number=0.72,
        gas=gas,
        viscosity_model=viscosity,
    )
    for field, reference, tolerance in (
        ("local_skin_friction_coefficient", "cf_local", 2e-3),
        ("momentum_thickness", "theta", 2e-3),
        ("boundary_layer_thickness", "delta99", 5e-3),
        ("displacement_thickness", "displacement_thickness", 5e-3),
    ):
        assert isinstance(getattr(result, field), float)
        assert getattr(result, field) == pytest.approx(r[reference], rel=tolerance)
    assert result.momentum_thickness == pytest.approx(
        r["x_m"] * float(result.average_skin_friction_coefficient) / 2.0, rel=1e-14
    )
    delta = float(result.boundary_layer_thickness)
    shear = float(result.wall_shear_stress)
    public = compressible_turbulent_boundary_layer_profile(
        np.array([0.0, delta]),
        r["edge_velocity"],
        r["edge_density"],
        r["edge_temperature"],
        delta,
        shear,
        transformation=CompressibleVelocityTransformation.VAN_DRIEST,
        temperature_velocity_relation=TemperatureVelocityRelation.WALZ,
        wall_temperature=r["wall_temperature"],
        gas=gas,
        viscosity_model=viscosity,
        prandtl_number=0.72,
    )
    fast = _van_driest_walz_profile(
        r["edge_velocity"],
        r["edge_density"],
        r["edge_temperature"],
        r["wall_temperature"],
        shear,
        gas=gas,
        viscosity_model=viscosity,
        prandtl_number=0.72,
    ).integrals(delta)
    assert_allclose(
        fast,
        [
            public.momentum_thickness,
            public.displacement_thickness,
            public.wake_parameter,
        ],
        rtol=1e-6,
    )
    assert public.momentum_thickness == pytest.approx(
        0.978 * float(result.momentum_thickness), rel=1e-4
    )
    assert public.edge_velocity_ratio == pytest.approx(0.99, rel=1e-8)
    assert result.displacement_thickness == pytest.approx(
        public.displacement_thickness, rel=1e-6
    )


@pytest.mark.parametrize("wall", [None, 206.0])
@pytest.mark.parametrize(
    "regime", [BoundaryLayerRegime.TURBULENT, BoundaryLayerRegime.TRANSITIONAL]
)
def test_wall_temperature_friction_and_momentum_integral(
    wall: float | None, regime: BoundaryLayerRegime
) -> None:
    args: dict[str, Any] = dict(
        edge_velocity=620.0,
        edge_density=0.5,
        edge_dynamic_viscosity=1.55e-5,
        regime=regime,
        transition_reynolds=5e5 if regime is BoundaryLayerRegime.TRANSITIONAL else None,
        compressible=True,
        turbulent_model=TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE,
        mach=2.0,
        edge_temperature=240.0,
        wall_temperature=wall,
    )
    x = 0.5
    result = flat_plate_boundary_layer(x, **args)
    fc, _, fx, recovery = _van_driest_factors(2.0, 240.0, wall, 1.55e-5)
    weight = 1.0 if wall is None else (wall / recovery) ** 0.064
    re = 0.5 * 620.0 * x / 1.55e-5
    friction_re = re
    if regime is BoundaryLayerRegime.TRANSITIONAL:
        laminar = flat_plate_boundary_layer(
            5e5 * 1.55e-5 / (0.5 * 620.0),
            620.0,
            0.5,
            1.55e-5,
            regime=BoundaryLayerRegime.LAMINAR,
            mach=2.0,
            edge_temperature=240.0,
            wall_temperature=wall,
            compressible=True,
        )
        # VD2 maps edge Re * Cbar_f to incompressible F_x Re * Cbar_i * F_C/w.
        laminar_drag = 5e5 * float(laminar.average_skin_friction_coefficient)
        start = _schlichting_start(laminar_drag * fx * fc / weight) / fx
        friction_re = start + re - 5e5
    average = weight / fc * _schlichting_average(fx * friction_re) * friction_re / re
    local = weight / fc * _schlichting_local(fx * friction_re)
    assert result.local_skin_friction_coefficient == pytest.approx(local, rel=1e-12)
    assert result.average_skin_friction_coefficient == pytest.approx(average, rel=1e-12)
    assert result.momentum_thickness == pytest.approx(x * average / 2.0, rel=1e-12)
    assert np.isfinite(result.boundary_layer_thickness)
    assert np.isfinite(result.displacement_thickness)
    h = x * 1e-5
    below = flat_plate_boundary_layer(x - h, **args)
    above = flat_plate_boundary_layer(x + h, **args)
    derivative = (
        (x + h) * float(above.average_skin_friction_coefficient)
        - (x - h) * float(below.average_skin_friction_coefficient)
    ) / (2 * h)
    assert result.local_skin_friction_coefficient == pytest.approx(derivative, rel=1e-4)


@pytest.mark.parametrize(
    "correction",
    [TurbulentModel.ECKERT, TurbulentModel.VAN_DRIEST_II],
)
@pytest.mark.parametrize(
    "regime", [BoundaryLayerRegime.TURBULENT, BoundaryLayerRegime.TRANSITIONAL]
)
def test_compressible_turbulent_thicknesses_are_nan(
    correction: TurbulentModel, regime: BoundaryLayerRegime
) -> None:
    arguments: dict[str, Any] = dict(
        edge_velocity=620.0,
        edge_density=0.5,
        edge_dynamic_viscosity=1.55e-5,
        compressible=True,
        turbulent_model=correction,
        mach=2.0,
        edge_temperature=240.0,
        wall_temperature=206.0,
    )
    x = np.array([0.01, 0.025, 0.5])
    result = flat_plate_boundary_layer(
        x,
        **arguments,
        regime=regime,
        transition_reynolds=5e5 if regime is BoundaryLayerRegime.TRANSITIONAL else None,
    )
    arguments["turbulent_model"] = None
    laminar = flat_plate_boundary_layer(
        x, **arguments, regime=BoundaryLayerRegime.LAMINAR
    )
    mask = (
        np.ones(3, dtype=bool) if regime is BoundaryLayerRegime.TURBULENT else x > 0.025
    )
    for field in (
        "boundary_layer_thickness",
        "displacement_thickness",
        "momentum_thickness",
    ):
        values = np.asarray(getattr(result, field))
        assert np.all(np.isnan(values[mask]))
        assert_allclose(
            values[~mask], np.asarray(getattr(laminar, field))[~mask], rtol=1e-14
        )
    for field in (
        "local_skin_friction_coefficient",
        "average_skin_friction_coefficient",
        "wall_shear_stress",
        "drag_per_unit_width",
    ):
        assert np.all(np.isfinite(getattr(result, field)))


def test_wall_temperature_broadcast_transition_and_laminar() -> None:
    args: dict[str, Any] = dict(
        edge_velocity=[620.0, 621.0],
        edge_density=0.5,
        edge_dynamic_viscosity=1.55e-5,
        mach=2.0,
        edge_temperature=240.0,
        wall_temperature=[206.0, 250.0],
    )
    mixed = flat_plate_boundary_layer(
        [[0.01], [0.5]],
        **args,
        regime=BoundaryLayerRegime.TRANSITIONAL,
        transition_reynolds=5e5,
        compressible=True,
        turbulent_model=TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE,
    )
    laminar = flat_plate_boundary_layer(
        [[0.01], [0.5]],
        **args,
        regime=BoundaryLayerRegime.LAMINAR,
        compressible=True,
    )
    new_laminar = flat_plate_boundary_layer(
        [[0.01], [0.5]],
        **args,
        regime=BoundaryLayerRegime.LAMINAR,
        compressible=True,
    )
    for field in (
        "boundary_layer_thickness",
        "displacement_thickness",
        "momentum_thickness",
        "local_skin_friction_coefficient",
        "average_skin_friction_coefficient",
        "wall_shear_stress",
        "drag_per_unit_width",
        "recovery_temperature",
        "wall_temperature",
    ):
        values = np.asarray(getattr(mixed, field))
        assert values.shape == (2, 2)
        assert values.dtype == np.float64
        assert np.all(np.isfinite(values))
        assert_allclose(values[0], np.asarray(getattr(laminar, field))[0])
        assert_allclose(getattr(new_laminar, field), getattr(laminar, field))
    assert_allclose(
        np.asarray(mixed.momentum_thickness)[1],
        0.5 * np.asarray(mixed.average_skin_friction_coefficient)[1] / 2.0,
    )


def test_wall_temperature_unresolvable_wake_reports_array_index() -> None:
    with (
        pytest.warns(ApplicabilityWarning, match="turbulent flat-plate correlations"),
        pytest.raises(ModelRangeError, match=r"index \(0, 1\).*Pi > 1"),
    ):
        flat_plate_boundary_layer(
            [[0.5, 1e-6]],
            620.0,
            0.5,
            1.55e-5,
            regime=BoundaryLayerRegime.TURBULENT,
            compressible=True,
            turbulent_model=TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE,
            mach=2.0,
            edge_temperature=240.0,
            wall_temperature=206.0,
        )


def test_wall_temperature_m10_sweep_has_no_low_re_tau_warning() -> None:
    from aerophysics.gas import PerfectGas
    from aerophysics.transport import SutherlandModel

    row = next(
        row for row in _wall_temperature_reference_rows() if row["case"] == "M10-AD"
    )
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        result = flat_plate_boundary_layer(
            np.linspace(0.01, 8.0, 200),
            float(row["edge_velocity"]),
            float(row["edge_density"]),
            float(row["edge_viscosity"]),
            regime=BoundaryLayerRegime.TURBULENT,
            compressible=True,
            turbulent_model=TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE,
            mach=float(row["mach"]),
            edge_temperature=float(row["edge_temperature"]),
            wall_temperature=float(row["wall_temperature"]),
            gas=PerfectGas(287.053, 1.4),
            viscosity_model=SutherlandModel(1.7894e-5, 288.15, 110.4),
        )

    assert np.all(np.isfinite(result.boundary_layer_thickness))
    # The first station is below the existing Re_x correlation range. The
    # validated low-Re_tau stations must not introduce additional warnings.
    assert len(captured) == 1
    assert captured[0].category is ApplicabilityWarning
    assert "turbulent flat-plate correlations" in str(captured[0].message)


_COMPRESSIBLE_MODES = tuple(mode for mode in TurbulentModel if mode.compressible)


@pytest.mark.parametrize("wall_temperature", [None, 300.0])
def test_compressible_laminar_momentum_integral(
    wall_temperature: float | None,
) -> None:
    viscosity = float(AIR_VISCOSITY.dynamic_viscosity(240.0))
    result = flat_plate_boundary_layer(
        [[0.2], [1.0]],
        1e8 * viscosity,
        1.0,
        viscosity,
        regime=BoundaryLayerRegime.LAMINAR,
        compressible=True,
        mach=[0.3, 2.0, 6.0, 15.0],
        edge_temperature=240.0,
        wall_temperature=wall_temperature,
    )
    assert_allclose(
        result.momentum_thickness,
        np.asarray(result.distance) * result.local_skin_friction_coefficient,
        rtol=1e-12,
    )


@pytest.mark.parametrize("mach", [0.0, 1e-7])
def test_compressible_laminar_blasius_limit(
    mach: float,
) -> None:
    viscosity = float(AIR_VISCOSITY.dynamic_viscosity(240.0))
    inputs = (np.array([0.1, 1.0]), 1e6 * viscosity, 1.0, viscosity)
    plain = flat_plate_boundary_layer(*inputs, regime=BoundaryLayerRegime.LAMINAR)
    corrected = flat_plate_boundary_layer(
        *inputs,
        regime=BoundaryLayerRegime.LAMINAR,
        compressible=True,
        mach=mach,
        edge_temperature=240.0,
        wall_temperature=240.0,
    )
    scale = np.asarray(plain.distance) / np.sqrt(plain.reynolds_number)
    for field, coefficient in (
        ("boundary_layer_thickness", 4.9100),
        ("displacement_thickness", 1.7208),
        ("momentum_thickness", 0.664),
    ):
        assert_allclose(getattr(plain, field), coefficient * scale, rtol=1e-12)
        assert_allclose(getattr(corrected, field), getattr(plain, field), rtol=2e-4)


@pytest.mark.parametrize("interaction", [0.999, 1.0, 1.001])
def test_viscous_interaction_threshold(
    interaction: float,
) -> None:
    # Tw=Te and mu_e=mu(Tw) give Cw=1 exactly; M=10 gives Re=1e6/chi^2.
    viscosity = float(AIR_VISCOSITY.dynamic_viscosity(240.0))
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        flat_plate_boundary_layer(
            1.0,
            1e6 / interaction**2,
            viscosity,
            viscosity,
            regime=BoundaryLayerRegime.LAMINAR,
            compressible=True,
            mach=10.0,
            edge_temperature=240.0,
            wall_temperature=240.0,
        )
    assert len(captured) == int(interaction > 1.0)
    if captured:
        assert captured[0].category is ApplicabilityWarning
        assert "chi_bar > 1 at 1 point(s)" in str(captured[0].message)
        assert "maximum chi_bar = 1.001" in str(captured[0].message)
        assert captured[0].filename == __file__


@pytest.mark.parametrize("correction", _COMPRESSIBLE_MODES)
@pytest.mark.parametrize(
    "viscosity_model", [AIR_VISCOSITY, SutherlandModel(2e-5, 300.0, 180.0)]
)
def test_viscous_interaction_uses_wall_viscosity_and_laminar_mask(
    correction: TurbulentModel,
    viscosity_model: Any,
) -> None:
    temperature = 240.0
    viscosity = float(viscosity_model.dynamic_viscosity(temperature))
    mach = np.array([[10.0], [12.0]])
    reynolds = np.array([1e5, 2e5, 1e7])
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        result = flat_plate_boundary_layer(
            1.0,
            mach * float(AIR.speed_of_sound(temperature)),
            reynolds * viscosity / (mach * float(AIR.speed_of_sound(temperature))),
            viscosity,
            regime=BoundaryLayerRegime.TRANSITIONAL,
            transition_reynolds=2e5,
            compressible=True,
            turbulent_model=correction,
            mach=mach,
            edge_temperature=temperature,
            viscosity_model=viscosity_model,
        )
    assert len(captured) == 1
    assert captured[0].category is ApplicabilityWarning
    assert "chi_bar > 1 at 4 point(s)" in str(captured[0].message)
    wall = np.asarray(result.wall_temperature)[:, :2]
    wall_c = temperature / wall * viscosity_model.dynamic_viscosity(wall) / viscosity
    chi = mach**3 * np.sqrt(wall_c / reynolds[:2])
    assert f"maximum chi_bar = {np.max(chi):.6g}" in str(captured[0].message)
    assert_allclose(wall, np.asarray(result.recovery_temperature)[:, :2])
    assert_allclose(
        np.asarray(result.wall_temperature)[:, :2], np.broadcast_to(wall, (2, 2))
    )


@pytest.mark.parametrize("correction", tuple(TurbulentModel))
def test_viscous_interaction_excludes_none_and_turbulent(
    correction: TurbulentModel,
) -> None:
    viscosity = float(AIR_VISCOSITY.dynamic_viscosity(240.0))
    thermal: dict[str, Any] = {}
    regime = BoundaryLayerRegime.LAMINAR
    if correction.compressible:
        thermal = {"mach": 15.0, "edge_temperature": 240.0}
        regime = BoundaryLayerRegime.TURBULENT
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        flat_plate_boundary_layer(
            1.0,
            15.0 * float(AIR.speed_of_sound(240.0)),
            1e6 * viscosity / (15.0 * float(AIR.speed_of_sound(240.0))),
            viscosity,
            regime=regime,
            compressible=correction.compressible,
            turbulent_model=correction if correction.compressible else None,
            **thermal,
        )
    assert not captured


@pytest.mark.parametrize("correction", _COMPRESSIBLE_MODES)
def test_compressible_transition_uses_new_laminar_thickness(
    correction: TurbulentModel,
) -> None:
    viscosity = float(AIR_VISCOSITY.dynamic_viscosity(240.0))
    x = np.array([0.5, 1.0, 10.0])
    thermal: dict[str, Any] = {
        "compressible": True,
        "turbulent_model": correction,
        "mach": 2.0,
        "edge_temperature": 240.0,
    }
    mixed = flat_plate_boundary_layer(
        x,
        1e6 * viscosity,
        1.0,
        viscosity,
        regime=BoundaryLayerRegime.TRANSITIONAL,
        transition_reynolds=1e6,
        **thermal,
    )
    thermal["turbulent_model"] = None
    laminar = flat_plate_boundary_layer(
        x[:2],
        1e6 * viscosity,
        1.0,
        viscosity,
        regime=BoundaryLayerRegime.LAMINAR,
        **thermal,
    )
    for field in (
        "boundary_layer_thickness",
        "displacement_thickness",
        "momentum_thickness",
    ):
        assert_allclose(
            np.asarray(getattr(mixed, field))[:2], getattr(laminar, field), rtol=1e-12
        )
        if correction is not TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE:
            assert np.isnan(np.asarray(getattr(mixed, field))[2])
    assert_allclose(
        np.asarray(mixed.momentum_thickness)[:2],
        x[:2] * np.asarray(mixed.local_skin_friction_coefficient)[:2],
        rtol=1e-12,
    )
    if correction is TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE:
        assert_allclose(
            np.asarray(mixed.momentum_thickness)[2],
            x[2] * np.asarray(mixed.average_skin_friction_coefficient)[2] / 2,
            rtol=1e-12,
        )


@pytest.mark.parametrize("regime", list(BoundaryLayerRegime))
@pytest.mark.parametrize("compressible", [False, True])
@pytest.mark.parametrize("model", [None, *TurbulentModel])
def test_model_selection_matrix(
    regime: BoundaryLayerRegime, compressible: bool, model: TurbulentModel | None
) -> None:
    valid = (
        model is None
        if regime is BoundaryLayerRegime.LAMINAR
        else (model is not None and model.compressible)
        if compressible
        else (model is None or not model.compressible)
    )
    arguments: dict[str, Any] = {
        "regime": regime,
        "compressible": compressible,
        "turbulent_model": model,
        "transition_reynolds": 5e5
        if regime is BoundaryLayerRegime.TRANSITIONAL
        else None,
    }
    if compressible:
        arguments.update(mach=2.0, edge_temperature=240.0)
    if valid:
        result = flat_plate_boundary_layer(1.0, 100.0, 1.0, 1e-5, **arguments)
        assert np.isfinite(result.local_skin_friction_coefficient)
    else:
        with pytest.raises(ValueError, match="turbulent_model"):
            flat_plate_boundary_layer(1.0, 100.0, 1.0, 1e-5, **arguments)


@pytest.mark.parametrize("value", [0, 1, "true", None, np.array([True])])
def test_compressible_requires_boolean(value: Any) -> None:
    with pytest.raises(ValueError, match="compressible must be a bool"):
        flat_plate_boundary_layer(
            1.0,
            100.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.LAMINAR,
            compressible=value,
        )


@pytest.mark.parametrize("compressible", [np.bool_(False), np.bool_(True)])
def test_numpy_boolean_compressibility(compressible: Any) -> None:
    thermal: dict[str, Any] = (
        {"mach": 2.0, "edge_temperature": 240.0} if compressible else {}
    )
    actual = flat_plate_boundary_layer(
        1.0,
        100.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.LAMINAR,
        compressible=compressible,
        **thermal,
    )
    expected = flat_plate_boundary_layer(
        1.0,
        100.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.LAMINAR,
        compressible=bool(compressible),
        **thermal,
    )
    assert actual == expected


@pytest.mark.parametrize(
    "regime", [BoundaryLayerRegime.TURBULENT, BoundaryLayerRegime.TRANSITIONAL]
)
def test_incompressible_default_is_schlichting(regime: BoundaryLayerRegime) -> None:
    arguments: dict[str, Any] = {
        "regime": regime,
        "transition_reynolds": 5e5
        if regime is BoundaryLayerRegime.TRANSITIONAL
        else None,
    }
    default = flat_plate_boundary_layer(1.0, 100.0, 1.0, 1e-5, **arguments)
    explicit = flat_plate_boundary_layer(
        1.0,
        100.0,
        1.0,
        1e-5,
        **arguments,
        turbulent_model=TurbulentModel.SCHLICHTING,
    )
    assert default == explicit


@pytest.mark.parametrize("field", ["mach", "edge_temperature", "wall_temperature"])
def test_incompressible_rejects_each_thermal_input(field: str) -> None:
    with pytest.raises(ValueError, match="thermal inputs require compressible=True"):
        flat_plate_boundary_layer(
            1.0,
            100.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.LAMINAR,
            **{field: 300.0},  # type: ignore[arg-type]
        )


@pytest.mark.parametrize("missing", ["mach", "edge_temperature"])
def test_compressible_requires_both_edge_thermal_inputs(missing: str) -> None:
    thermal: dict[str, Any] = {"mach": 2.0, "edge_temperature": 240.0}
    del thermal[missing]
    with pytest.raises(ValueError, match="mach and edge_temperature are required"):
        flat_plate_boundary_layer(
            1.0,
            100.0,
            1.0,
            1e-5,
            regime=BoundaryLayerRegime.LAMINAR,
            compressible=True,
            **thermal,
        )


def test_turbulent_models_and_compressibility() -> None:
    assert {model.value for model in TurbulentModel} == {
        "schlichting",
        "schlichting_log_wake",
        "eckert",
        "van_driest_ii",
        "van_driest_ii_wall_temperature",
    }
    assert not TurbulentModel.SCHLICHTING.compressible
    assert not TurbulentModel.SCHLICHTING_LOG_WAKE.compressible
    assert all(
        model.compressible
        for model in (
            TurbulentModel.ECKERT,
            TurbulentModel.VAN_DRIEST_II,
            TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE,
        )
    )
    for removed in ("POWER_LAW", "ECKERT_POWER_LAW", "ECKERT_SCHLICHTING"):
        assert not hasattr(TurbulentModel, removed)
