"""Public exceptions and warnings used by :mod:`aerophysics`."""


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


__all__ = [
    "ApplicabilityWarning",
    "ExpansionConvergenceError",
    "ModelRangeError",
    "NoAttachedShockError",
    "ShockConvergenceError",
]
