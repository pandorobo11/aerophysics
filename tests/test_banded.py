"""Exact regression checks against SciPy's LAPACK banded-solve paths."""

import numpy as np
import pytest
from scipy.linalg import solve_banded as scipy_solve_banded

from aerophysics._array import FloatArray
from aerophysics._banded import solve_banded, solve_tridiagonal


def random_bands(lower: int, upper: int, size: int, seed: int) -> FloatArray:
    """Make a strictly diagonally dominant float64 banded system."""
    rng = np.random.default_rng(seed)
    bands = rng.uniform(-1.0, 1.0, (lower + upper + 1, size))
    bands[upper] = lower + upper + 2.0
    return bands


@pytest.mark.parametrize("size", [2, 17, 64])
def test_tridiagonal_is_bit_identical_and_preserves_inputs(size: int) -> None:
    bands = random_bands(1, 1, size, size)
    vector = np.random.default_rng(size + 1).normal(size=size)
    saved_bands, saved_vector = bands.copy(), vector.copy()
    expected = scipy_solve_banded((1, 1), bands, vector)
    actual = solve_tridiagonal(bands, vector)
    assert np.array_equal(actual, expected)
    assert actual.tobytes() == expected.tobytes()
    assert np.array_equal(bands, saved_bands)
    assert np.array_equal(vector, saved_vector)


@pytest.mark.parametrize(("lower", "upper"), [(5, 5), (9, 9), (2, 4)])
@pytest.mark.parametrize("size", [23, 64])
@pytest.mark.parametrize("overwrite", [False, True])
def test_expanded_bands_are_bit_identical(
    lower: int, upper: int, size: int, overwrite: bool
) -> None:
    compact = random_bands(lower, upper, size, size + lower)
    expanded = np.zeros((2 * lower + upper + 1, size), dtype=np.float64, order="F")
    expanded[lower:] = compact
    vector = np.random.default_rng(size).normal(size=size)
    saved_bands, saved_vector = expanded.copy(), vector.copy()
    expected = scipy_solve_banded((lower, upper), compact, vector, check_finite=False)
    actual = solve_banded(
        lower, upper, expanded, vector, overwrite_ab=overwrite, overwrite_b=overwrite
    )
    assert np.array_equal(actual, expected)
    assert actual.tobytes() == expected.tobytes()
    if not overwrite:
        assert np.array_equal(expanded, saved_bands)
        assert np.array_equal(vector, saved_vector)
    else:
        # The Newton path really uses its disposable workspace and RHS in place.
        assert not np.array_equal(expanded, saved_bands)
        assert np.shares_memory(actual, vector)


@pytest.mark.parametrize("band", [1, 5, 9])
def test_pivoting_is_bit_identical(band: int) -> None:
    size = 31
    compact = random_bands(band, band, size, band)
    compact[band, 0] = 0.0  # Force a row interchange at the first pivot.
    vector = np.linspace(-1.0, 1.0, size)
    expected = scipy_solve_banded((band, band), compact, vector)
    if band == 1:
        actual = solve_tridiagonal(compact, vector)
    else:
        expanded = np.zeros((3 * band + 1, size), dtype=np.float64)
        expanded[band:] = compact
        actual = solve_banded(band, band, expanded, vector)
    assert np.array_equal(actual, expected)
    assert actual.tobytes() == expected.tobytes()


@pytest.mark.parametrize("band", [1, 5, 9])
def test_singular_system_matches_scipy(band: int) -> None:
    compact = np.zeros((2 * band + 1, 31), dtype=np.float64)
    vector = np.ones(31)
    with pytest.raises(np.linalg.LinAlgError, match=r"^singular matrix$"):
        scipy_solve_banded((band, band), compact, vector)
    with pytest.raises(np.linalg.LinAlgError, match=r"^singular matrix$"):
        if band == 1:
            solve_tridiagonal(compact, vector)
        else:
            solve_banded(band, band, np.zeros((3 * band + 1, 31)), vector)


@pytest.mark.parametrize("value", [np.nan])
@pytest.mark.parametrize("location", ["diagonal", "unused_corner", "rhs"])
def test_nonfinite_input_matches_scipy(value: float, location: str) -> None:
    bands = random_bands(1, 1, 17, 1)
    vector = np.ones(17)
    if location == "rhs":
        vector[3] = value
    elif location == "diagonal":
        bands[1, 3] = value
    else:
        bands[0, 0] = value
    with pytest.raises(ValueError) as expected:
        scipy_solve_banded((1, 1), bands, vector, check_finite=True)
    with pytest.raises(ValueError) as actual:
        solve_tridiagonal(bands, vector, check_finite=True)
    assert str(actual.value) == str(expected.value)


def test_disabled_finite_check_ignores_unused_nonfinite_corner() -> None:
    bands = random_bands(1, 1, 17, 2)
    bands[0, 0] = np.nan
    vector = np.ones(17)
    assert np.array_equal(
        solve_tridiagonal(bands, vector, check_finite=False),
        scipy_solve_banded((1, 1), bands, vector, check_finite=False),
    )


def test_tridiagonal_accepts_readonly_strided_inputs() -> None:
    bands = random_bands(1, 1, 34, 3)[:, ::2]
    vector = np.linspace(0.0, 1.0, 34)[::2]
    bands.flags.writeable = vector.flags.writeable = False
    assert np.array_equal(
        solve_tridiagonal(bands, vector), scipy_solve_banded((1, 1), bands, vector)
    )
