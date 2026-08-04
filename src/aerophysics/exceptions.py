"""Public exceptions and warnings used by :mod:`aerophysics`."""


class ModelRangeError(ValueError):
    """Raised when an input is outside a model's implemented range."""


class ApplicabilityWarning(UserWarning):
    """Warn when a computable result is outside a model's validated range."""


class NoAttachedShockError(ValueError):
    """Raised when no attached oblique- or conical-shock solution exists."""


__all__ = ["ApplicabilityWarning", "ModelRangeError", "NoAttachedShockError"]
