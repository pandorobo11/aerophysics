"""Public exceptions and warnings used by :mod:`aerophysics`."""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from aerophysics.marching_boundary_layer import MarchingBoundaryLayerResult


class ModelRangeError(ValueError):
    """Raised when an input is outside a model's implemented range."""


class ApplicabilityWarning(UserWarning):
    """Warn when a computable result is outside a model's validated range."""


class NoAttachedShockError(ValueError):
    """Raised when no attached oblique- or conical-shock solution exists."""


class ShockConvergenceError(RuntimeError):
    """Raised when a shock solver cannot resolve the requested numerical solution."""


class ExpansionConvergenceError(RuntimeError):
    """Raised when a centered expansion cannot be resolved numerically."""


class BoundaryLayerSeparationError(ModelRangeError):
    """Raised when a marched boundary layer separates.

    Parameters
    ----------
    message:
        Description of the failure.
    separation_distance:
        Extrapolated separation distance along the wall in m.
    result:
        Result at the requested stations solved upstream of separation.
    """

    def __init__(
        self,
        message: str,
        *,
        separation_distance: float,
        result: "MarchingBoundaryLayerResult",
    ) -> None:
        super().__init__(message)
        self.separation_distance = separation_distance
        self.result = result


class BoundaryLayerConvergenceError(RuntimeError):
    """Raised when a boundary-layer march cannot resolve the solution."""


__all__ = [
    "ApplicabilityWarning",
    "BoundaryLayerConvergenceError",
    "BoundaryLayerSeparationError",
    "ExpansionConvergenceError",
    "ModelRangeError",
    "NoAttachedShockError",
    "ShockConvergenceError",
]
