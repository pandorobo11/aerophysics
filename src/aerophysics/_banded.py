"""Float64 LAPACK solves for the marching solver's prevalidated systems.

Callers provide correctly shaped arrays with at least two unknowns; the march
always has at least three grid nodes. These are not general ``solve_banded``
replacements: dtype, shape, batching, and empty/singleton handling are omitted.
"""

from collections.abc import Callable
from typing import cast

import numpy as np
from scipy.linalg import get_lapack_funcs

from aerophysics._array import FloatArray

_gtsv = cast(
    Callable[..., tuple[FloatArray, FloatArray, FloatArray, FloatArray, int]],
    get_lapack_funcs("gtsv", dtype=np.float64),
)
_gbsv = cast(
    Callable[..., tuple[FloatArray, object, FloatArray, int]],
    get_lapack_funcs("gbsv", dtype=np.float64),
)


def _check_info(info: int) -> None:
    """Use the same LAPACK failure translation as SciPy's ``solve_banded``."""
    if info > 0:
        raise np.linalg.LinAlgError("singular matrix")
    if info < 0:
        raise ValueError(f"illegal value in {-info}-th argument of internal gbsv/gtsv")


def solve_tridiagonal(
    bands: FloatArray, vector: FloatArray, *, check_finite: bool = True
) -> FloatArray:
    """Solve SciPy's (1, 1) band layout, preserving both input arrays."""
    # Check the entire arrays, including unused band corners, just as SciPy does.
    if check_finite and (not np.isfinite(bands).all() or not np.isfinite(vector).all()):
        raise ValueError("array must not contain infs or NaNs")
    _, _, _, solution, info = _gtsv(
        bands[2, :-1], bands[1, :], bands[0, 1:], vector, False, False, False, False
    )
    _check_info(info)
    return solution


def solve_banded(
    lower: int,
    upper: int,
    bands: FloatArray,
    vector: FloatArray,
    *,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
) -> FloatArray:
    """Solve expanded ``gbsv`` bands without finite checks.

    ``bands`` has ``2*lower + upper + 1`` rows; its bottom ``lower + upper + 1``
    rows contain SciPy's compact bands. Fortran order and ``overwrite_ab=True``
    avoid the LAPACK input copy. Only enable overwrites for disposable arrays.
    """
    _, _, solution, info = _gbsv(
        lower, upper, bands, vector, overwrite_ab=overwrite_ab, overwrite_b=overwrite_b
    )
    _check_info(info)
    return solution
