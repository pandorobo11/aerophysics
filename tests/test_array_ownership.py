"""Ownership and immutability checks for public vectorized results."""

from collections.abc import Callable

import numpy as np
import pytest
from numpy.typing import NDArray

from aerophysics import AIR_HARMONIC_OSCILLATOR, standard_atmosphere
from aerophysics.boundary_layer import BoundaryLayerRegime, flat_plate_boundary_layer
from aerophysics.detached_shock import DetachedShockGeometry, shock_standoff_distance
from aerophysics.expansion import prandtl_meyer_expansion
from aerophysics.isentropic import isentropic_ratios
from aerophysics.shocks import normal_shock

type FloatArray = NDArray[np.float64]
type VectorCalculation = Callable[[FloatArray], float | FloatArray]


@pytest.mark.parametrize(
    "values, calculate",
    [
        (
            np.asarray([0.0, 1000.0]),
            lambda values: standard_atmosphere(values).geometric_altitude,
        ),
        (
            np.asarray([1.0, 2.0]),
            lambda values: isentropic_ratios(values).mach,
        ),
        (
            np.asarray([2.0, 3.0]),
            lambda values: normal_shock(values).upstream_mach,
        ),
        (
            np.asarray([2.0, 3.0]),
            lambda values: prandtl_meyer_expansion(values, 0.01).upstream_mach,
        ),
        (
            np.asarray([1.0, 2.0]),
            lambda values: (
                flat_plate_boundary_layer(
                    values,
                    100.0,
                    1.0,
                    1.0e-5,
                    regime=BoundaryLayerRegime.LAMINAR,
                ).distance
            ),
        ),
        (
            np.asarray([500.0, 700.0]),
            lambda values: AIR_HARMONIC_OSCILLATOR.state(values, 101_325.0).temperature,
        ),
        (
            np.asarray([2.0, 3.0]),
            lambda values: (
                shock_standoff_distance(
                    values,
                    0.1,
                    geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
                ).upstream_mach
            ),
        ),
    ],
    ids=(
        "atmosphere",
        "isentropic",
        "normal-shock",
        "expansion",
        "boundary-layer",
        "real-gas",
        "detached-shock",
    ),
)
def test_vectorized_results_are_owned_read_only_snapshots(
    values: FloatArray, calculate: VectorCalculation
) -> None:
    source = values.copy()
    result = calculate(source)

    assert isinstance(result, np.ndarray)
    assert result.dtype == np.dtype(np.float64)
    assert result.shape == source.shape
    assert result.flags.owndata
    assert not result.flags.writeable
    assert not np.shares_memory(source, result)

    expected = result.copy()
    source[0] *= 1.1
    np.testing.assert_array_equal(result, expected)
    with pytest.raises(ValueError, match="read-only"):
        result[0] = 0.0


def test_scalar_results_remain_python_floats() -> None:
    assert isinstance(standard_atmosphere(0.0).pressure, float)
