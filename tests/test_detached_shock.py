"""Tests for detached-shock engineering correlations."""

import json
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

from aerophysics.detached_shock import (
    DetachedShockGeometry,
    DetachedShockModel,
    billig_shock_shape,
    compare_standoff_distances,
    seiff_standoff_distance,
    seiff_standoff_distance_from_mach,
    shock_standoff_distance,
)
from aerophysics.shocks import normal_shock

REFERENCE_PATH = (
    Path(__file__).parent
    / "reference_data"
    / "compressible_flow"
    / "detached_shock_sources.json"
)


def test_ambrosio_wortman_reference_values() -> None:
    reference = json.loads(REFERENCE_PATH.read_text())
    for case in reference["representative_values"]:
        mach = case["mach"]
        sphere = shock_standoff_distance(
            mach,
            1.0,
            geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
        )
        cylinder = shock_standoff_distance(
            mach,
            1.0,
            geometry=DetachedShockGeometry.CYLINDRICAL_NOSE_2D,
        )
        assert sphere.normalized_standoff_distance == pytest.approx(
            case["sphere_delta_over_rn"]
        )
        assert cylinder.normalized_standoff_distance == pytest.approx(
            case["cylinder_delta_over_rn"]
        )
        sphere_shape = billig_shock_shape(
            mach,
            1.0,
            [0.0],
            geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
        )
        cylinder_shape = billig_shock_shape(
            mach,
            1.0,
            [0.0],
            geometry=DetachedShockGeometry.CYLINDRICAL_NOSE_2D,
        )
        assert sphere_shape.vertex_curvature_radius == pytest.approx(
            case["sphere_curvature_over_rn"]
        )
        assert cylinder_shape.vertex_curvature_radius == pytest.approx(
            case["cylinder_curvature_over_rn"]
        )
        assert seiff_standoff_distance_from_mach(
            mach, 1.0
        ).normalized_standoff_distance == pytest.approx(
            case["seiff_delta_over_rn_gamma_1_4"]
        )


def test_standoff_vectorization_and_broadcasting() -> None:
    result = shock_standoff_distance(
        [[2.0], [4.0]],
        [1.0, 2.0, 3.0],
        geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
    )
    for value in (
        result.upstream_mach,
        result.nose_radius,
        result.normalized_standoff_distance,
        result.standoff_distance,
    ):
        assert isinstance(value, np.ndarray)
        assert value.shape == (2, 3)
        assert value.dtype == np.float64
    assert_allclose(
        result.standoff_distance,
        np.asarray(result.normalized_standoff_distance)
        * np.asarray(result.nose_radius),
    )
    scalar = shock_standoff_distance(
        4.0,
        2.0,
        geometry=DetachedShockGeometry.CYLINDRICAL_NOSE_2D,
    )
    assert isinstance(scalar.standoff_distance, float)


def test_seiff_low_level_and_normal_shock_api_agree() -> None:
    ratio = normal_shock(4.0).static_density_ratio
    direct = seiff_standoff_distance(ratio, 2.0)
    from_mach = seiff_standoff_distance_from_mach(4.0, 2.0)
    assert direct.upstream_mach is None
    assert from_mach.upstream_mach == 4.0
    assert direct.density_ratio == pytest.approx(ratio)
    assert direct.normalized_standoff_distance == pytest.approx(0.78 / ratio)
    assert direct.standoff_distance == pytest.approx(from_mach.standoff_distance)
    assert direct.geometry is DetachedShockGeometry.AXISYMMETRIC_SPHERE
    assert direct.model is DetachedShockModel.SEIFF


def test_seiff_vectorizes_and_broadcasts() -> None:
    result = seiff_standoff_distance([[2.0], [4.0]], [1.0, 2.0])
    assert isinstance(result.density_ratio, np.ndarray)
    assert result.density_ratio.shape == (2, 2)
    assert isinstance(result.standoff_distance, np.ndarray)
    assert result.standoff_distance.shape == (2, 2)


def test_billig_shape_metadata_vertex_curvature_and_symmetry() -> None:
    radius = 2.0
    transverse = np.array([-0.1, 0.0, 0.1])
    result = billig_shock_shape(
        4.0,
        radius,
        transverse,
        geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
    )
    assert result.model is DetachedShockModel.BILLIG
    assert result.standoff_model is DetachedShockModel.AMBROSIO_WORTMAN
    assert result.shock_x.shape == (3,)
    assert_allclose(result.shock_y, transverse)
    assert result.shock_x[0] == pytest.approx(result.shock_x[2])
    assert result.shock_x[1] == pytest.approx(radius + float(result.standoff_distance))
    curvature_from_vertex = 2.0 * (result.shock_x[1] - result.shock_x[2]) / 0.1**2
    assert curvature_from_vertex == pytest.approx(
        1.0 / float(result.vertex_curvature_radius), rel=2e-4
    )


def test_billig_shape_case_axis_and_asymptote() -> None:
    result = billig_shock_shape(
        [2.0, 4.0],
        1.0,
        [0.0, 1.0, 1.0e8],
        geometry=DetachedShockGeometry.CYLINDRICAL_NOSE_2D,
    )
    assert result.shock_x.shape == (2, 3)
    assert result.shock_y.shape == (2, 3)
    beta = np.arcsin(1.0 / 4.0)
    far_slope = (result.shock_x[1, 2] - result.shock_x[1, 1]) / (
        result.shock_y[1, 2] - result.shock_y[1, 1]
    )
    assert far_slope == pytest.approx(-1.0 / np.tan(beta), rel=2e-7)


def test_comparison_differences() -> None:
    result = compare_standoff_distances([2.0, 4.0], 3.0)
    aw = np.asarray(result.ambrosio_wortman.normalized_standoff_distance)
    seiff = np.asarray(result.seiff.normalized_standoff_distance)
    assert_allclose(result.normalized_standoff_difference, seiff - aw)
    assert_allclose(result.relative_difference, (seiff - aw) / aw)
    assert_allclose(result.standoff_distance_difference, 3.0 * (seiff - aw))


@pytest.mark.parametrize(
    ("function", "arguments"),
    [
        (
            shock_standoff_distance,
            (1.0, 1.0),
        ),
        (seiff_standoff_distance, (1.0, 1.0)),
        (seiff_standoff_distance, (2.0, 0.0)),
        (seiff_standoff_distance_from_mach, (np.inf, 1.0)),
    ],
)
def test_invalid_physical_inputs(
    function: object, arguments: tuple[object, ...]
) -> None:
    if function is shock_standoff_distance:
        with pytest.raises(ValueError):
            shock_standoff_distance(
                *arguments,  # type: ignore[arg-type]
                geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
            )
    else:
        with pytest.raises(ValueError):
            function(*arguments)  # type: ignore[operator]


def test_invalid_geometry_and_shape_coordinates() -> None:
    with pytest.raises(ValueError, match="DetachedShockGeometry"):
        shock_standoff_distance(2.0, 1.0, geometry="sphere")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="one-dimensional"):
        billig_shock_shape(
            2.0,
            1.0,
            [[0.0, 1.0]],
            geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
        )
    with pytest.raises(ValueError, match="must not be empty"):
        billig_shock_shape(
            2.0,
            1.0,
            [],
            geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
        )


def test_case_inputs_must_broadcast() -> None:
    with pytest.raises(ValueError, match="broadcastable"):
        shock_standoff_distance(
            [2.0, 3.0],
            [1.0, 2.0, 3.0],
            geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
        )
