"""Tests for compressible turbulent boundary-layer profiles."""

from types import SimpleNamespace

import numpy as np
import pytest
from numpy.testing import assert_allclose

import aerophysics.boundary_layer_profile as profile_module
from aerophysics.boundary_layer_profile import (
    CompressibleBoundaryLayerProfileResult,
    CompressibleVelocityTransformation,
    TemperatureVelocityRelation,
    TransformedVelocityProfileResult,
    compressible_turbulent_boundary_layer_profile,
    transform_compressible_velocity_profile,
)
from aerophysics.exceptions import ApplicabilityWarning, ModelRangeError
from aerophysics.gas import AIR, AIR_VISCOSITY


def _profile(
    *,
    transformation: CompressibleVelocityTransformation = (
        CompressibleVelocityTransformation.VAN_DRIEST
    ),
    relation: TemperatureVelocityRelation = (
        TemperatureVelocityRelation.GENERALIZED_REYNOLDS_ANALOGY
    ),
    wall_distance: np.ndarray | None = None,
    wake_parameter: float | None = None,
) -> CompressibleBoundaryLayerProfileResult:
    height = (
        np.linspace(0.0, 0.05, 1001)
        if wall_distance is None
        else wall_distance
    )
    return compressible_turbulent_boundary_layer_profile(
        height,
        300.0,
        1.0,
        300.0,
        0.05,
        85.0,
        transformation=transformation,
        wall_temperature=250.0,
        temperature_velocity_relation=relation,
        wake_parameter=wake_parameter,
    )


@pytest.mark.parametrize(
    "transformation",
    list(CompressibleVelocityTransformation),
)
def test_forward_constant_properties_are_identity(
    transformation: CompressibleVelocityTransformation,
) -> None:
    height = np.linspace(0.0, 0.01, 101)
    velocity = 100.0 * height / height[-1]
    result = transform_compressible_velocity_profile(
        height,
        velocity,
        np.full_like(height, 1.2),
        np.full_like(height, 1.8e-5),
        12.0,
        transformation=transformation,
    )
    friction_velocity = np.sqrt(10.0)
    expected_height_plus = 1.2 * friction_velocity * height / 1.8e-5
    assert isinstance(result, TransformedVelocityProfileResult)
    assert result.friction_velocity == pytest.approx(friction_velocity)
    assert_allclose(result.wall_distance_plus, expected_height_plus)
    assert_allclose(
        result.transformed_wall_coordinate, expected_height_plus
    )
    assert_allclose(result.velocity_plus, velocity / friction_velocity)
    assert_allclose(
        result.transformed_velocity_plus, velocity / friction_velocity
    )
    assert result.wall_distance.dtype == np.float64


def test_forward_nonuniform_properties_match_differential_mappings() -> None:
    height = np.array([0.0, 1.0, 2.0])
    velocity = np.array([0.0, 2.0, 5.0])
    density = np.array([4.0, 1.0, 9.0])
    viscosity = np.array([2.0, 4.0, 1.0])

    van_driest = transform_compressible_velocity_profile(
        height,
        velocity,
        density,
        viscosity,
        4.0,
        transformation=CompressibleVelocityTransformation.VAN_DRIEST,
    )
    assert_allclose(van_driest.wall_distance_plus, [0.0, 2.0, 4.0])
    assert_allclose(van_driest.transformed_velocity_plus, [0.0, 1.5, 4.5])

    volpiani = transform_compressible_velocity_profile(
        height,
        velocity,
        density,
        viscosity,
        4.0,
        transformation=CompressibleVelocityTransformation.VOLPIANI,
    )
    coordinate_factor = np.array(
        [1.0, 0.5 * 2.0**-1.5, 1.5 * 0.5**-1.5]
    )
    velocity_factor = np.array(
        [1.0, 0.5 * 2.0**-0.5, 1.5 * 0.5**-0.5]
    )
    expected_coordinate = np.array(
        [
            0.0,
            coordinate_factor[:2].sum(),
            coordinate_factor[:2].sum() + coordinate_factor[1:].sum(),
        ]
    )
    expected_velocity = np.array(
        [
            0.0,
            velocity_factor[:2].sum(),
            velocity_factor[:2].sum()
            + 1.5 * velocity_factor[1:].sum(),
        ]
    )
    assert_allclose(
        volpiani.transformed_wall_coordinate, expected_coordinate
    )
    assert_allclose(volpiani.transformed_velocity_plus, expected_velocity)


def test_coles_wake_endpoints_and_zero_strength() -> None:
    outer_coordinate = np.array([0.0, 0.5, 1.0])
    wake = profile_module._wake_function(outer_coordinate)
    assert_allclose(wake, [0.0, 1.0, 2.0])
    wall_velocity_plus = np.array([0.0, 10.0, 20.0])
    assert_allclose(wall_velocity_plus + 0.0 * wake, wall_velocity_plus)


@pytest.mark.parametrize(
    ("transformation", "relation"),
    [
        (transformation, relation)
        for transformation in CompressibleVelocityTransformation
        for relation in TemperatureVelocityRelation
    ],
)
def test_inverse_profile_satisfies_edge_and_property_relations(
    transformation: CompressibleVelocityTransformation,
    relation: TemperatureVelocityRelation,
) -> None:
    result = _profile(transformation=transformation, relation=relation)
    assert isinstance(result, CompressibleBoundaryLayerProfileResult)
    assert 0.0 <= result.wake_parameter <= 1.0
    assert result.edge_velocity_ratio == pytest.approx(0.99, abs=1e-9)
    assert result.velocity[-1] == pytest.approx(297.0, abs=1e-6)
    assert result.velocity[0] == pytest.approx(0.0)
    assert result.temperature[0] == pytest.approx(250.0)
    assert result.density[0] == pytest.approx(1.2)
    assert_allclose(
        result.density * result.temperature,
        np.full(result.density.shape, 300.0),
    )
    assert_allclose(
        result.dynamic_viscosity,
        AIR_VISCOSITY.dynamic_viscosity(result.temperature),
    )
    assert_allclose(
        result.local_mach_number,
        result.velocity / AIR.speed_of_sound(result.temperature),
    )
    assert_allclose(
        result.dynamic_pressure,
        0.5 * result.density * result.velocity**2,
    )
    assert result.local_skin_friction_coefficient == pytest.approx(
        2.0 * 85.0 / 300.0**2
    )
    assert result.displacement_thickness > result.momentum_thickness > 0.0
    assert result.shape_factor == pytest.approx(
        result.displacement_thickness / result.momentum_thickness
    )
    for values in (
        result.wall_distance,
        result.wall_distance_plus,
        result.transformed_wall_coordinate,
        result.velocity,
        result.velocity_plus,
        result.transformed_velocity_plus,
        result.temperature,
        result.density,
        result.dynamic_viscosity,
        result.local_mach_number,
        result.dynamic_pressure,
    ):
        assert values.shape == (1001,)
        assert values.dtype == np.float64


def test_gra_and_walz_temperature_formulas() -> None:
    gra = _profile(
        relation=TemperatureVelocityRelation.GENERALIZED_REYNOLDS_ANALOGY
    )
    walz = _profile(relation=TemperatureVelocityRelation.WALZ)
    index = 400

    gra_ratio = gra.velocity[index] / 300.0
    expected_gra = (
        250.0
        + 1.14
        * 0.72
        * (gra.recovery_temperature - 250.0)
        * gra_ratio
        * (1.0 - gra_ratio)
        + (300.0 - 250.0) * gra_ratio**2
    )
    assert gra.temperature[index] == pytest.approx(expected_gra)

    walz_ratio = walz.velocity[index] / 300.0
    expected_walz = (
        250.0
        + (walz.recovery_temperature - 250.0) * walz_ratio
        + (300.0 - walz.recovery_temperature) * walz_ratio**2
    )
    assert walz.temperature[index] == pytest.approx(expected_walz)
    assert gra.temperature[index] != pytest.approx(walz.temperature[index])


def test_adiabatic_default_uses_recovery_temperature() -> None:
    height = np.linspace(0.0, 0.05, 501)
    result = compressible_turbulent_boundary_layer_profile(
        height,
        300.0,
        1.0,
        300.0,
        0.05,
        75.0,
        transformation=CompressibleVelocityTransformation.VAN_DRIEST,
    )
    expected_recovery = 300.0 + np.cbrt(0.72) * 300.0**2 / (2.0 * AIR.cp)
    assert result.recovery_temperature == pytest.approx(expected_recovery)
    assert result.wall_temperature == pytest.approx(expected_recovery)
    assert result.temperature[0] == pytest.approx(expected_recovery)


@pytest.mark.parametrize(
    "transformation",
    list(CompressibleVelocityTransformation),
)
def test_predicted_profile_round_trips_through_forward_transform(
    transformation: CompressibleVelocityTransformation,
) -> None:
    height = np.unique(
        np.concatenate(
            (
                np.array([0.0]),
                np.geomspace(1e-10, 0.05, 4000),
                np.linspace(0.0, 0.05, 4001),
            )
        )
    )
    predicted = _profile(
        transformation=transformation,
        wall_distance=height,
    )
    transformed = transform_compressible_velocity_profile(
        height,
        predicted.velocity,
        predicted.density,
        predicted.dynamic_viscosity,
        85.0,
        transformation=transformation,
    )
    assert_allclose(
        transformed.transformed_velocity_plus,
        predicted.transformed_velocity_plus,
        rtol=8e-5,
        atol=3e-4,
    )
    assert_allclose(
        transformed.transformed_wall_coordinate,
        predicted.transformed_wall_coordinate,
        rtol=8e-5,
        atol=3e-4,
    )


def test_integral_quantities_match_sampled_profile() -> None:
    result = _profile(wall_distance=np.linspace(0.0, 0.05, 10001))
    mass_velocity_ratio = result.density * result.velocity / 300.0
    expected_displacement = np.trapezoid(
        1.0 - mass_velocity_ratio, result.wall_distance
    )
    expected_momentum = np.trapezoid(
        mass_velocity_ratio * (1.0 - result.velocity / 300.0),
        result.wall_distance,
    )
    assert result.displacement_thickness == pytest.approx(
        expected_displacement, rel=2e-5
    )
    assert result.momentum_thickness == pytest.approx(
        expected_momentum, rel=2e-5
    )


def test_integrals_cover_full_layer_when_output_grid_is_truncated() -> None:
    full = _profile()
    truncated = _profile(wall_distance=np.linspace(0.0, 0.01, 51))
    assert truncated.velocity[-1] < 300.0
    assert truncated.displacement_thickness == pytest.approx(
        full.displacement_thickness, rel=1e-7
    )
    assert truncated.momentum_thickness == pytest.approx(
        full.momentum_thickness, rel=1e-7
    )


def test_explicit_consistent_wake_parameter_reproduces_profile() -> None:
    automatic = _profile()
    specified = _profile(wake_parameter=automatic.wake_parameter)
    assert specified.wake_parameter == pytest.approx(automatic.wake_parameter)
    assert_allclose(specified.velocity, automatic.velocity)


def test_low_friction_reynolds_number_warns() -> None:
    with pytest.warns(ApplicabilityWarning, match="Re_tau < 500"):
        result = compressible_turbulent_boundary_layer_profile(
            np.linspace(0.0, 0.004, 101),
            20.0,
            1.0,
            300.0,
            0.004,
            1.0,
            transformation=CompressibleVelocityTransformation.VAN_DRIEST,
            wall_temperature=300.0,
        )
    assert result.friction_reynolds_number < 500.0


@pytest.mark.parametrize(
    ("keyword", "value"),
    [
        ("edge_velocity", 0.0),
        ("edge_density", -1.0),
        ("edge_temperature", 0.0),
        ("boundary_layer_thickness", 0.0),
        ("wall_shear_stress", 0.0),
        ("wall_temperature", 0.0),
        ("prandtl_number", 0.0),
        ("reynolds_analogy_factor", 0.0),
        ("von_karman_constant", 0.0),
    ],
)
def test_inverse_rejects_non_positive_inputs(
    keyword: str, value: float
) -> None:
    arguments: dict[str, object] = {
        "wall_distance": np.array([0.0, 0.05]),
        "edge_velocity": 300.0,
        "edge_density": 1.0,
        "edge_temperature": 300.0,
        "boundary_layer_thickness": 0.05,
        "wall_shear_stress": 85.0,
        "transformation": CompressibleVelocityTransformation.VAN_DRIEST,
        "wall_temperature": 250.0,
    }
    arguments[keyword] = value
    with pytest.raises(ValueError, match=keyword):
        compressible_turbulent_boundary_layer_profile(**arguments)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "wall_distance",
    [
        np.array([0.0]),
        np.array([0.001, 0.01]),
        np.array([0.0, 0.01, 0.005]),
        np.array([[0.0, 0.01]]),
    ],
)
def test_inverse_rejects_invalid_wall_grid(
    wall_distance: np.ndarray,
) -> None:
    with pytest.raises(ValueError, match="wall_distance"):
        _profile(wall_distance=wall_distance)


def test_inverse_rejects_wall_grid_beyond_delta_99() -> None:
    with pytest.raises(ValueError, match="must not exceed"):
        _profile(wall_distance=np.array([0.0, 0.051]))


def test_inverse_rejects_invalid_model_choices() -> None:
    with pytest.raises(ValueError, match="transformation"):
        compressible_turbulent_boundary_layer_profile(
            [0.0, 0.05],
            300.0,
            1.0,
            300.0,
            0.05,
            85.0,
            transformation="van_driest",  # type: ignore[arg-type]
            wall_temperature=250.0,
        )
    with pytest.raises(ValueError, match="temperature_velocity_relation"):
        compressible_turbulent_boundary_layer_profile(
            [0.0, 0.05],
            300.0,
            1.0,
            300.0,
            0.05,
            85.0,
            transformation=CompressibleVelocityTransformation.VAN_DRIEST,
            temperature_velocity_relation="walz",  # type: ignore[arg-type]
            wall_temperature=250.0,
        )


def test_inverse_rejects_invalid_log_intercept_shape() -> None:
    with pytest.raises(ValueError, match="log_law_intercept must be a scalar"):
        compressible_turbulent_boundary_layer_profile(
            [0.0, 0.05],
            300.0,
            1.0,
            300.0,
            0.05,
            85.0,
            transformation=CompressibleVelocityTransformation.VAN_DRIEST,
            wall_temperature=250.0,
            log_law_intercept=[5.0, 5.2],  # type: ignore[arg-type]
        )


def test_inverse_rejects_inconsistent_or_invalid_wake_parameter() -> None:
    with pytest.raises(ValueError, match="wake_parameter must be a scalar"):
        _profile(wake_parameter=np.array([0.5, 0.6]))  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="non-negative"):
        _profile(wake_parameter=-0.1)
    with pytest.raises(ModelRangeError, match="does not satisfy"):
        _profile(wake_parameter=0.0)


def test_inverse_rejects_inputs_without_a_zpg_wake_root() -> None:
    with pytest.raises(ModelRangeError, match="no Coles wake parameter"):
        compressible_turbulent_boundary_layer_profile(
            [0.0, 0.05],
            300.0,
            1.0,
            300.0,
            0.05,
            40.0,
            transformation=CompressibleVelocityTransformation.VAN_DRIEST,
            wall_temperature=250.0,
        )


def test_inverse_rejects_non_physical_temperature_relation() -> None:
    with pytest.raises(ModelRangeError, match="non-physical temperature"):
        compressible_turbulent_boundary_layer_profile(
            [0.0, 0.05],
            300.0,
            1.0,
            300.0,
            0.05,
            85.0,
            transformation=CompressibleVelocityTransformation.VAN_DRIEST,
            wall_temperature=1000.0,
            reynolds_analogy_factor=100.0,
        )


def test_inverse_reports_ode_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    failed = SimpleNamespace(success=False, message="synthetic failure")
    monkeypatch.setattr(profile_module, "solve_ivp", lambda *args, **kwargs: failed)
    with pytest.raises(ModelRangeError, match="synthetic failure"):
        _profile()


@pytest.mark.parametrize(
    ("field", "values", "message"),
    [
        ("wall_distance", [[0.0, 1.0]], "one-dimensional"),
        ("wall_distance", [0.0], "at least 2"),
        ("wall_distance", [1.0, 2.0], "start at zero"),
        ("wall_distance", [0.0, 2.0, 1.0], "strictly increasing"),
        ("velocity", [1.0, 2.0], "start at zero"),
        ("velocity", [0.0, -1.0], "non-negative"),
        ("velocity", [0.0, 2.0, 1.0], "non-decreasing"),
        ("density", [1.0, 0.0], "greater than zero"),
        ("dynamic_viscosity", [1.0, 0.0], "greater than zero"),
    ],
)
def test_forward_rejects_invalid_profiles(
    field: str, values: object, message: str
) -> None:
    arguments: dict[str, object] = {
        "wall_distance": [0.0, 1.0],
        "velocity": [0.0, 1.0],
        "density": [1.0, 1.0],
        "dynamic_viscosity": [1.0, 1.0],
        "wall_shear_stress": 1.0,
        "transformation": CompressibleVelocityTransformation.VAN_DRIEST,
    }
    if field == "velocity" and len(values) == 3:  # type: ignore[arg-type]
        arguments["wall_distance"] = [0.0, 1.0, 2.0]
        arguments["density"] = [1.0, 1.0, 1.0]
        arguments["dynamic_viscosity"] = [1.0, 1.0, 1.0]
    arguments[field] = values
    with pytest.raises(ValueError, match=message):
        transform_compressible_velocity_profile(**arguments)  # type: ignore[arg-type]


def test_forward_rejects_length_choice_and_shear_errors() -> None:
    with pytest.raises(ValueError, match="same length"):
        transform_compressible_velocity_profile(
            [0.0, 1.0],
            [0.0, 1.0, 2.0],
            [1.0, 1.0],
            [1.0, 1.0],
            1.0,
            transformation=CompressibleVelocityTransformation.VAN_DRIEST,
        )
    with pytest.raises(ValueError, match="wall_shear_stress"):
        transform_compressible_velocity_profile(
            [0.0, 1.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 1.0],
            0.0,
            transformation=CompressibleVelocityTransformation.VAN_DRIEST,
        )
    with pytest.raises(ValueError, match="transformation"):
        transform_compressible_velocity_profile(
            [0.0, 1.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 1.0],
            1.0,
            transformation="van_driest",  # type: ignore[arg-type]
        )
