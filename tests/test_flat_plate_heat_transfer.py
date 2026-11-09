"""Flat-plate wall heat transfer and combined mean-property profiles."""

from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics import flat_plate_boundary_layer_profile
from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    TurbulentModel,
    flat_plate_boundary_layer,
)
from aerophysics.boundary_layer_profile import (
    CompressibleVelocityTransformation,
    compressible_turbulent_boundary_layer_profile,
)
from aerophysics.gas import AIR
from aerophysics.transport import AIR_VISCOSITY

MACH = 5.0
EDGE_TEMPERATURE = 220.0
VISCOSITY = float(AIR_VISCOSITY.dynamic_viscosity(EDGE_TEMPERATURE))
VELOCITY = MACH * float(
    np.sqrt(AIR.heat_capacity_ratio * AIR.specific_gas_constant * EDGE_TEMPERATURE)
)
DENSITY = 1e7 * VISCOSITY / VELOCITY


def plate(**kwargs: Any) -> Any:
    settings: dict[str, Any] = {
        "regime": BoundaryLayerRegime.TURBULENT,
        "compressible": True,
        "turbulent_model": TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE,
        "mach": MACH,
        "edge_temperature": EDGE_TEMPERATURE,
        "wall_temperature": 300.0,
    }
    settings.update(kwargs)
    return flat_plate_boundary_layer(
        1.0,
        VELOCITY,
        DENSITY,
        VISCOSITY,
        **settings,
    )


def test_turbulent_heat_transfer_follows_reynolds_analogy() -> None:
    result = plate(reynolds_analogy_factor=1.2)
    stanton = 1.2 * result.local_skin_friction_coefficient / 2.0
    coefficient = stanton * DENSITY * VELOCITY * AIR.cp
    assert_allclose(result.stanton_number, stanton, rtol=1e-14)
    assert_allclose(result.heat_transfer_coefficient, coefficient, rtol=1e-14)
    assert_allclose(
        result.wall_heat_flux,
        coefficient * (result.recovery_temperature - 300.0),
        rtol=1e-14,
    )


def test_laminar_heat_transfer_uses_prandtl_analogy() -> None:
    result = plate(
        regime=BoundaryLayerRegime.LAMINAR,
        turbulent_model=None,
        reynolds_analogy_factor=2.0,
        prandtl_number=0.8,
    )
    assert_allclose(
        result.stanton_number,
        0.8 ** (-2.0 / 3.0) * result.local_skin_friction_coefficient / 2.0,
        rtol=1e-14,
    )


def test_wall_heat_flux_sign_and_adiabatic_wall() -> None:
    assert plate().wall_heat_flux > 0.0
    hot = plate(wall_temperature=2000.0)
    assert hot.wall_heat_flux < 0.0
    assert hot.heat_transfer_coefficient > 0.0
    adiabatic = plate(wall_temperature=None)
    assert adiabatic.wall_heat_flux == 0.0
    assert adiabatic.stanton_number > 0.0


def test_transitional_heat_transfer_switches_analogy_at_transition() -> None:
    result: Any = flat_plate_boundary_layer(
        np.array([0.01, 1.0]),
        VELOCITY,
        DENSITY,
        VISCOSITY,
        regime=BoundaryLayerRegime.TRANSITIONAL,
        compressible=True,
        turbulent_model=TurbulentModel.VAN_DRIEST_II_WALL_TEMPERATURE,
        transition_reynolds=5e5,
        mach=MACH,
        edge_temperature=EDGE_TEMPERATURE,
        wall_temperature=300.0,
    )
    assert_allclose(
        2.0 * result.stanton_number / result.local_skin_friction_coefficient,
        [0.72 ** (-2.0 / 3.0), 1.14],
        rtol=1e-14,
    )


def test_incompressible_result_has_no_heat_transfer() -> None:
    result = flat_plate_boundary_layer(
        1.0, 10.0, 1.2, 1.8e-5, regime=BoundaryLayerRegime.LAMINAR
    )
    assert result.stanton_number is None
    assert result.heat_transfer_coefficient is None
    assert result.wall_heat_flux is None


@pytest.mark.parametrize("factor", [0.0, -1.0, float("nan"), float("inf")])
def test_invalid_reynolds_analogy_factor(factor: float) -> None:
    with pytest.raises(ValueError, match="reynolds_analogy_factor"):
        plate(reynolds_analogy_factor=factor)


def test_low_mach_laminar_profile_matches_blasius() -> None:
    mach = 0.01
    velocity = mach * float(
        np.sqrt(AIR.heat_capacity_ratio * AIR.specific_gas_constant * EDGE_TEMPERATURE)
    )
    # Howarth's Blasius values of f'(eta), eta = y sqrt(U/(nu x)).
    eta = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    expected = np.array([0.0, 0.32979, 0.62977, 0.84605, 0.95552])
    result = flat_plate_boundary_layer_profile(
        1.0,
        velocity,
        1e6 * VISCOSITY / velocity,
        EDGE_TEMPERATURE,
        regime=BoundaryLayerRegime.LAMINAR,
        wall_distance_ratio=eta / 4.91,
    )
    assert result.regime is BoundaryLayerRegime.LAMINAR
    assert_allclose(result.velocity / velocity, expected, atol=2e-3)
    assert_allclose(result.temperature, EDGE_TEMPERATURE, rtol=1e-3)


def test_laminar_profile_edge_states_and_wall_flux() -> None:
    ratio = np.array([0.0, 5e-4, 1e-3, 1.0])
    result = flat_plate_boundary_layer_profile(
        1.0,
        VELOCITY,
        DENSITY,
        EDGE_TEMPERATURE,
        regime=BoundaryLayerRegime.LAMINAR,
        wall_temperature=300.0,
        wall_distance_ratio=ratio,
    )
    layer: Any = result.boundary_layer
    assert_allclose(result.wall_distance, ratio * layer.boundary_layer_thickness)
    assert_allclose(result.velocity[0], 0.0, atol=1e-9)
    assert_allclose(result.temperature[0], 300.0, rtol=1e-8)
    assert_allclose(result.velocity[-1] / VELOCITY, 0.99, rtol=1e-8)
    assert_allclose(result.density * result.temperature, DENSITY * EDGE_TEMPERATURE)
    assert_allclose(
        result.local_mach_number,
        result.velocity
        / np.sqrt(
            AIR.heat_capacity_ratio * AIR.specific_gas_constant * result.temperature
        ),
    )
    gradient = (result.temperature[1] - result.temperature[0]) / result.wall_distance[1]
    conductivity = result.dynamic_viscosity[0] * AIR.cp / 0.72
    # The profile's own wall gradient supports the analogy-based flux.
    assert_allclose(conductivity * gradient, layer.wall_heat_flux, rtol=0.01)
    wall_shear = (
        result.dynamic_viscosity[0] * result.velocity[1] / result.wall_distance[1]
    )
    assert_allclose(wall_shear, layer.wall_shear_stress, rtol=0.01)


def test_adiabatic_laminar_profile_uses_recovery_temperature() -> None:
    result = flat_plate_boundary_layer_profile(
        1.0, VELOCITY, DENSITY, EDGE_TEMPERATURE, regime=BoundaryLayerRegime.LAMINAR
    )
    layer: Any = result.boundary_layer
    assert_allclose(result.temperature[0], layer.recovery_temperature, rtol=1e-5)
    assert result.wall_distance.size == 201


def test_turbulent_profile_uses_plate_shear_and_thickness() -> None:
    ratio = np.array([0.0, 1e-5, 0.1, 1.0])
    result = flat_plate_boundary_layer_profile(
        1.0,
        VELOCITY,
        DENSITY,
        EDGE_TEMPERATURE,
        regime=BoundaryLayerRegime.TURBULENT,
        wall_temperature=300.0,
        wall_distance_ratio=ratio,
    )
    layer: Any = result.boundary_layer
    assert result.regime is BoundaryLayerRegime.TURBULENT
    direct = compressible_turbulent_boundary_layer_profile(
        ratio * layer.boundary_layer_thickness,
        VELOCITY,
        DENSITY,
        EDGE_TEMPERATURE,
        layer.boundary_layer_thickness,
        layer.wall_shear_stress,
        transformation=CompressibleVelocityTransformation.VAN_DRIEST,
        wall_temperature=300.0,
    )
    assert_allclose(result.velocity, direct.velocity)
    assert_allclose(result.temperature, direct.temperature)
    assert_allclose(
        result.dynamic_viscosity[0]
        * AIR.cp
        / 0.72
        * (result.temperature[1] - result.temperature[0])
        / result.wall_distance[1],
        layer.wall_heat_flux,
        rtol=0.01,
    )


def test_transitional_profile_selects_local_regime() -> None:
    arguments = (1.0, VELOCITY, DENSITY, EDGE_TEMPERATURE)
    before = flat_plate_boundary_layer_profile(
        *arguments,
        regime=BoundaryLayerRegime.TRANSITIONAL,
        transition_reynolds=2e7,
        wall_temperature=300.0,
    )
    after = flat_plate_boundary_layer_profile(
        *arguments,
        regime=BoundaryLayerRegime.TRANSITIONAL,
        transition_reynolds=5e5,
        wall_temperature=300.0,
    )
    assert before.regime is BoundaryLayerRegime.LAMINAR
    assert after.regime is BoundaryLayerRegime.TURBULENT


@pytest.mark.parametrize(
    ("ratio", "message"),
    [
        ([0.0, 1.5], "wall_distance_ratio"),
        ([0.1, 0.5], "start|zero"),
        ([0.0, 0.5, 0.5], "increasing"),
    ],
)
def test_invalid_wall_distance_ratio(ratio: list[float], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        flat_plate_boundary_layer_profile(
            1.0,
            VELOCITY,
            DENSITY,
            EDGE_TEMPERATURE,
            regime=BoundaryLayerRegime.LAMINAR,
            wall_distance_ratio=ratio,
        )


def test_profile_requires_positive_distance() -> None:
    with pytest.raises(ValueError, match="distance"):
        flat_plate_boundary_layer_profile(
            0.0, VELOCITY, DENSITY, EDGE_TEMPERATURE, regime=BoundaryLayerRegime.LAMINAR
        )
