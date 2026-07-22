"""Package-level smoke tests."""

from aerophysics import __version__
from aerophysics.exceptions import ApplicabilityWarning, ModelRangeError


def test_version() -> None:
    assert __version__ == "0.1.0"


def test_public_diagnostics() -> None:
    assert issubclass(ModelRangeError, ValueError)
    assert issubclass(ApplicabilityWarning, UserWarning)
