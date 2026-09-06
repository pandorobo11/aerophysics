"""Internal array conversion and validation helpers."""

import numpy as np
from numpy.typing import ArrayLike, NDArray

type FloatArray = NDArray[np.float64]
type FloatResult = float | FloatArray


def as_float_array(value: ArrayLike, *, name: str) -> tuple[FloatArray, bool]:
    """Convert a scalar or array-like value to a finite float64 array."""
    scalar = np.ndim(value) == 0
    try:
        array = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must contain real numeric values") from error
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array, scalar


def return_float(array: FloatArray, *, scalar: bool) -> FloatResult:
    """Return a scalar or an owned, read-only float64 result array."""
    if scalar:
        return float(array)
    result = np.array(array, dtype=np.float64, copy=True)
    result.flags.writeable = False
    return result
