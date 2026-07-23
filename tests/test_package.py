"""Package-level smoke tests."""

from aerophysics import (
    ShockBranch,
    __version__,
    normal_shock,
    oblique_shock,
    prandtl_meyer_expansion,
)
from aerophysics.exceptions import (
    ApplicabilityWarning,
    ModelRangeError,
    NoAttachedShockError,
)


def test_version() -> None:
    assert __version__ == "0.2.0"


def test_public_diagnostics() -> None:
    assert issubclass(ModelRangeError, ValueError)
    assert issubclass(NoAttachedShockError, ValueError)
    assert issubclass(ApplicabilityWarning, UserWarning)


def test_primary_compressible_flow_api_is_exported() -> None:
    assert ShockBranch.WEAK.value == "weak"
    assert normal_shock(2.0).downstream_mach < 1.0
    assert oblique_shock(2.0, 0.1).downstream_mach > 1.0
    assert prandtl_meyer_expansion(2.0, 0.1).downstream_mach > 2.0
