"""Tests for advanced GUI calculation adapters."""

import numpy as np
import pytest

from aerophysics import AIR_NASA9, protrusion_drag
from aerophysics.boundary_layer_profile import (
    CompressibleVelocityTransformation,
    TemperatureVelocityRelation,
    compressible_turbulent_boundary_layer_profile,
)
from aerophysics.gui.advanced_adapters import (
    boundary_layer_profiles,
    protrusion_condition,
    protrusion_sweep,
    thermochemistry_condition,
    thermochemistry_sweep,
    wall_normal_grid,
)


def test_profile_adapter_matches_core_and_compares_models() -> None:
    calculation = boundary_layer_profiles(
        edge_velocity=300.0,
        edge_density=1.0,
        edge_temperature=300.0,
        boundary_layer_thickness=0.05,
        wall_shear_stress=85.0,
        transformations=tuple(CompressibleVelocityTransformation),
        temperature_velocity_relation=(
            TemperatureVelocityRelation.GENERALIZED_REYNOLDS_ANALOGY
        ),
        wall_temperature=250.0,
        wake_parameter=None,
        points=101,
    )
    grid = wall_normal_grid(0.05, 101)
    direct = compressible_turbulent_boundary_layer_profile(
        grid,
        300.0,
        1.0,
        300.0,
        0.05,
        85.0,
        transformation=CompressibleVelocityTransformation.VAN_DRIEST,
        wall_temperature=250.0,
    )
    assert len(calculation.result.rows) == 202
    assert len(calculation.profiles) == 2
    assert calculation.result.rows[100]["velocity"] == pytest.approx(
        direct.velocity[-1]
    )
    assert calculation.profiles[0].density[-1] == pytest.approx(direct.density[-1])


def test_profile_grid_and_model_validation() -> None:
    grid = wall_normal_grid(0.1, 51)
    assert grid[0] == 0.0
    assert grid[-1] == pytest.approx(0.1)
    assert np.all(np.diff(grid) > 0.0)
    with pytest.raises(ValueError, match="positive"):
        wall_normal_grid(0.0, 51)
    with pytest.raises(ValueError, match="between"):
        wall_normal_grid(0.1, 50)
    with pytest.raises(ValueError, match="unique"):
        boundary_layer_profiles(
            edge_velocity=300.0,
            edge_density=1.0,
            edge_temperature=300.0,
            boundary_layer_thickness=0.05,
            wall_shear_stress=85.0,
            transformations=(),
            temperature_velocity_relation=TemperatureVelocityRelation.WALZ,
            wall_temperature=None,
            wake_parameter=None,
            points=51,
        )


@pytest.mark.parametrize("shape", ["rectangle", "triangle", "ellipse"])
def test_protrusion_shapes_match_direct_model(shape: str) -> None:
    adapted = protrusion_condition(
        drag_coefficient=1.2,
        height=0.01,
        base_width=0.005,
        shape=shape,
        edge_velocity=100.0,
        edge_density=1.2,
        boundary_layer_thickness=0.05,
    )
    assert float(adapted.rows[0]["direct_drag"]) > 0.0  # type: ignore[arg-type]
    if shape == "rectangle":
        direct = protrusion_drag(1.2, 0.01, 0.005, 100.0, 1.2, 0.05)
        assert adapted.rows[0]["direct_drag"] == pytest.approx(direct.direct_drag)


def test_protrusion_provided_profile_compressibility_and_sweep_gaps() -> None:
    height = np.array([0.0, 0.01, 0.02])
    provided = protrusion_condition(
        drag_coefficient=1.0,
        height=0.01,
        base_width=0.005,
        shape="csv",
        edge_velocity=100.0,
        edge_density=1.0,
        boundary_layer_thickness=0.02,
        profile_height=height,
        profile_velocity=np.array([0.0, 50.0, 100.0]),
        profile_density=np.ones(3),
        shape_height=height,
        shape_width=np.array([0.005, 0.003, 0.0]),
    )
    assert provided.rows[0]["profile"] == "provided"
    compressed = protrusion_condition(
        drag_coefficient=1.0,
        height=0.01,
        base_width=0.005,
        shape="rectangle",
        edge_velocity=500.0,
        edge_density=0.4,
        boundary_layer_thickness=0.02,
        mach=2.0,
        edge_temperature=220.0,
    )
    assert compressed.rows[0]["compressibility_applied"] == "True"
    swept = protrusion_sweep(
        sweep_field="height",
        start=0.005,
        stop=0.03,
        points=3,
        drag_coefficient=1.0,
        height=0.01,
        base_width=0.005,
        shape="csv",
        edge_velocity=100.0,
        edge_density=1.0,
        boundary_layer_thickness=0.02,
        shape_height=height,
        shape_width=np.array([0.005, 0.003, 0.0]),
    )
    assert swept.rows[-1]["status"] == "invalid"
    with pytest.raises(ValueError, match="Mach sweep"):
        protrusion_sweep(
            sweep_field="mach",
            start=0.0,
            stop=2.0,
            points=3,
            drag_coefficient=1.0,
            height=0.01,
            base_width=0.005,
            shape="rectangle",
            edge_velocity=100.0,
            edge_density=1.0,
            boundary_layer_thickness=0.02,
        )
    with pytest.raises(ValueError, match="shape must"):
        protrusion_condition(
            drag_coefficient=1.0,
            height=0.01,
            base_width=0.005,
            shape="other",
            edge_velocity=100.0,
            edge_density=1.0,
            boundary_layer_thickness=0.02,
        )
    with pytest.raises(ValueError, match="CSV shape requires"):
        protrusion_condition(
            drag_coefficient=1.0,
            height=0.01,
            base_width=0.005,
            shape="csv",
            edge_velocity=100.0,
            edge_density=1.0,
            boundary_layer_thickness=0.02,
        )
    with pytest.raises(ValueError, match="unsupported"):
        protrusion_sweep(
            sweep_field="other",
            start=0.0,
            stop=1.0,
            points=3,
            drag_coefficient=1.0,
            height=0.01,
            base_width=0.005,
            shape="rectangle",
            edge_velocity=100.0,
            edge_density=1.0,
            boundary_layer_thickness=0.02,
        )


def test_thermochemistry_adapter_matches_core_and_extrapolates() -> None:
    single = thermochemistry_condition(
        temperature=1000.0,
        pressure=101_325.0,
        reference_temperature=298.15,
        models=("NASA9",),
        allow_extrapolation=False,
    )
    assert single.rows[0]["cp"] == pytest.approx(AIR_NASA9.cp(1000.0))
    compared = thermochemistry_sweep(
        start=200.0,
        stop=6000.0,
        points=3,
        pressure=101_325.0,
        reference_temperature=298.15,
        models=("NASA7", "NASA9"),
        allow_extrapolation=False,
    )
    assert len(compared.rows) == 6
    extrapolated = thermochemistry_condition(
        temperature=7000.0,
        pressure=101_325.0,
        reference_temperature=298.15,
        models=("NASA7",),
        allow_extrapolation=True,
    )
    assert extrapolated.warnings
    with pytest.raises(ValueError, match="unique"):
        thermochemistry_condition(
            temperature=300.0,
            pressure=101_325.0,
            reference_temperature=298.15,
            models=(),
            allow_extrapolation=False,
        )
    with pytest.raises(ValueError, match="NASA7 or NASA9"):
        thermochemistry_condition(
            temperature=300.0,
            pressure=101_325.0,
            reference_temperature=298.15,
            models=("other",),
            allow_extrapolation=False,
        )
