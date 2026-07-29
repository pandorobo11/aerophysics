"""Package-level smoke tests."""

from aerophysics import (
    BoundaryLayerRegime,
    CompressibilityCorrection,
    ShockBranch,
    TurbulentCorrelation,
    __version__,
    flat_plate_boundary_layer,
    normal_shock,
    oblique_shock,
    prandtl_meyer_expansion,
    protrusion_drag,
)
from aerophysics.exceptions import (
    ApplicabilityWarning,
    ModelRangeError,
    NoAttachedShockError,
)


def test_version() -> None:
    assert __version__ == "0.3.0"


def test_public_diagnostics() -> None:
    assert issubclass(ModelRangeError, ValueError)
    assert issubclass(NoAttachedShockError, ValueError)
    assert issubclass(ApplicabilityWarning, UserWarning)


def test_primary_compressible_flow_api_is_exported() -> None:
    assert ShockBranch.WEAK.value == "weak"
    assert normal_shock(2.0).downstream_mach < 1.0
    assert oblique_shock(2.0, 0.1).downstream_mach > 1.0
    assert prandtl_meyer_expansion(2.0, 0.1).downstream_mach > 2.0


def test_primary_boundary_layer_api_is_exported() -> None:
    result = flat_plate_boundary_layer(
        1.0,
        10.0,
        1.0,
        1e-5,
        regime=BoundaryLayerRegime.TURBULENT,
        turbulent_correlation=TurbulentCorrelation.POWER_LAW,
        compressibility_correction=CompressibilityCorrection.NONE,
    )
    assert result.drag_per_unit_width > 0.0


def test_primary_protrusion_drag_api_is_exported() -> None:
    result = protrusion_drag(1.0, 0.01, 0.005, 10.0, 1.0, 0.02)
    assert result.direct_drag > 0.0
