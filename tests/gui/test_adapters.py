"""Tests for pure GUI calculation adapters."""

import numpy as np
import pytest

from aerophysics import FlightCondition
from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    CompressibilityCorrection,
    TurbulentCorrelation,
)
from aerophysics.gui.adapters import (
    FlightCase,
    conical_shock_condition,
    conical_shock_sweep,
    expansion_condition,
    expansion_sweep,
    flat_plate,
    flat_plate_sweep,
    flight_condition,
    flight_sweep,
    isentropic_condition,
    isentropic_sweep,
    normal_shock_condition,
    normal_shock_sweep,
    oblique_shock_condition,
    oblique_shock_sweep,
    sweep_values,
)
from aerophysics.isentropic import MachBranch, isentropic_ratios
from aerophysics.shocks import ShockBranch, conical_shock, oblique_shock


def test_sweep_values_validation_and_spacing() -> None:
    assert sweep_values(0.0, 1.0, 3).tolist() == [0.0, 0.5, 1.0]
    assert sweep_values(0.01, 1.0, 3, log=True) == pytest.approx([0.01, 0.1, 1.0])
    for arguments in (
        (1.0, 0.0, 3),
        (0.0, 1.0, 1),
        (0.0, 1.0, 502),
        (np.nan, 1.0, 3),
    ):
        with pytest.raises(ValueError):
            sweep_values(*arguments)
    with pytest.raises(ValueError, match="positive"):
        sweep_values(0.0, 1.0, 3, log=True)


def test_flight_adapter_matches_core_and_builds_case() -> None:
    adapted = flight_condition(
        geometric_altitude=10_000.0,
        motion=0.8,
        motion_basis="mach",
        characteristic_length=2.0,
    )
    direct = FlightCondition.from_mach(10_000.0, 0.8, 2.0)
    row = adapted.rows[0]
    assert row["velocity"] == pytest.approx(direct.velocity)
    assert row["dynamic_pressure"] == pytest.approx(direct.dynamic_pressure)
    case = FlightCase.from_row(row)
    assert case.mach == 0.8
    assert case.density == pytest.approx(direct.atmosphere.density)
    with pytest.raises(ValueError, match="numeric"):
        FlightCase.from_row({**row, "mach": None})


def test_flight_velocity_and_sweeps() -> None:
    velocity = flight_condition(
        geometric_altitude=0.0,
        motion=100.0,
        motion_basis="velocity",
        characteristic_length=None,
    )
    assert velocity.rows[0]["velocity"] == pytest.approx(100.0)
    assert velocity.rows[0]["reynolds_number"] is None
    altitude = flight_sweep(
        fixed_altitude=0.0,
        fixed_motion=0.5,
        motion_basis="mach",
        sweep_field="altitude",
        start=0.0,
        stop=2000.0,
        points=3,
        characteristic_length=1.0,
    )
    assert len(altitude.rows) == 3
    motion = flight_sweep(
        fixed_altitude=1000.0,
        fixed_motion=0.5,
        motion_basis="mach",
        sweep_field="motion",
        start=0.2,
        stop=1.0,
        points=3,
        characteristic_length=None,
    )
    assert motion.rows[-1]["mach"] == pytest.approx(1.0)
    with pytest.raises(ValueError, match="motion_basis"):
        flight_condition(
            geometric_altitude=0.0,
            motion=1.0,
            motion_basis="other",
            characteristic_length=None,
        )


def test_isentropic_adapter_forward_inverse_and_mass_flux() -> None:
    forward = isentropic_condition(
        input_value=2.0,
        input_basis="mach",
        total_pressure=101_325.0,
        total_temperature=300.0,
    )
    direct = isentropic_ratios(2.0)
    assert forward.rows[0]["total_pressure_ratio"] == pytest.approx(
        direct.total_pressure_ratio
    )
    assert forward.rows[0]["mass_flux"] is not None
    assert forward.rows[0]["choked_mass_flux"] is not None
    area_inverse = isentropic_condition(
        input_value=2.0,
        input_basis="area_ratio",
        branch=MachBranch.SUPERSONIC,
    )
    assert float(area_inverse.rows[0]["mach"]) > 1.0  # type: ignore[arg-type]
    pressure_inverse = isentropic_condition(
        input_value=float(direct.total_pressure_ratio),
        input_basis="pressure_ratio",
    )
    assert pressure_inverse.rows[0]["mach"] == pytest.approx(2.0)
    zero = isentropic_condition(input_value=0.0, input_basis="mach")
    assert zero.rows[0]["area_ratio"] is None


def test_isentropic_sweep_and_validation() -> None:
    result = isentropic_sweep(
        input_basis="temperature_ratio",
        branch=MachBranch.SUBSONIC,
        start=1.0,
        stop=2.0,
        points=3,
    )
    assert len(result.rows) == 3
    with pytest.raises(ValueError, match="specified together"):
        isentropic_condition(
            input_value=1.0,
            input_basis="mach",
            total_pressure=101_325.0,
        )
    with pytest.raises(ValueError, match="unsupported"):
        isentropic_condition(input_value=1.0, input_basis="other")


def test_normal_shock_adapter_and_sweep() -> None:
    single = normal_shock_condition(upstream_mach=2.0)
    assert single.rows[0]["downstream_mach"] == pytest.approx(0.57735, rel=1e-5)
    assert float(single.rows[0]["pitot_pressure_ratio"]) > 1.0  # type: ignore[arg-type]
    sweep = normal_shock_sweep(start=1.0, stop=3.0, points=3)
    assert len(sweep.rows) == 3
    first_ratio = sweep.rows[0]["total_pressure_ratio"]
    last_ratio = sweep.rows[-1]["total_pressure_ratio"]
    assert isinstance(first_ratio, float)
    assert isinstance(last_ratio, float)
    assert last_ratio < first_ratio


def test_expansion_adapter_and_limit_rows() -> None:
    single = expansion_condition(upstream_mach=2.0, turn_angle=float(np.deg2rad(10.0)))
    assert float(single.rows[0]["downstream_mach"]) > 2.0  # type: ignore[arg-type]
    sweep = expansion_sweep(
        fixed_mach=2.0,
        fixed_turn_angle=0.1,
        sweep_field="turn_angle",
        start=0.0,
        stop=float(np.deg2rad(130.0)),
        points=3,
    )
    assert sweep.rows[-1]["status"] == "limit_exceeded"
    mach_sweep = expansion_sweep(
        fixed_mach=2.0,
        fixed_turn_angle=0.01,
        sweep_field="mach",
        start=1.0,
        stop=3.0,
        points=3,
    )
    assert len(mach_sweep.rows) == 3
    with pytest.raises(ValueError, match="sweep_field"):
        expansion_sweep(
            fixed_mach=2.0,
            fixed_turn_angle=0.1,
            sweep_field="other",
            start=1.0,
            stop=2.0,
            points=2,
        )
    with pytest.raises(ValueError, match="sweep_field"):
        flight_sweep(
            fixed_altitude=0.0,
            fixed_motion=1.0,
            motion_basis="mach",
            sweep_field="other",
            start=0.0,
            stop=1.0,
            points=2,
            characteristic_length=None,
        )


def test_shock_adapter_matches_core() -> None:
    adapted = oblique_shock_condition(
        upstream_mach=2.0,
        deflection_angle=np.deg2rad(10.0),
        branch=ShockBranch.WEAK,
    )
    direct = oblique_shock(2.0, np.deg2rad(10.0), ShockBranch.WEAK)
    assert adapted.rows[0]["downstream_mach"] == pytest.approx(direct.downstream_mach)
    assert adapted.rows[0]["status"] == "ok"


def test_shock_sweep_preserves_non_attached_points() -> None:
    result = oblique_shock_sweep(
        fixed_mach=2.0,
        fixed_deflection=np.deg2rad(10.0),
        branch=ShockBranch.WEAK,
        sweep_field="deflection",
        start=0.0,
        stop=np.deg2rad(30.0),
        points=4,
    )
    assert result.rows[0]["status"] == "ok"
    assert result.rows[-1]["status"] == "no_attached_shock"
    assert result.rows[-1]["shock_angle"] is None
    mach = oblique_shock_sweep(
        fixed_mach=2.0,
        fixed_deflection=np.deg2rad(5.0),
        branch=ShockBranch.STRONG,
        sweep_field="mach",
        start=1.1,
        stop=2.0,
        points=3,
    )
    assert len(mach.rows) == 3
    with pytest.raises(ValueError, match="sweep_field"):
        oblique_shock_sweep(
            fixed_mach=2.0,
            fixed_deflection=0.1,
            branch=ShockBranch.WEAK,
            sweep_field="other",
            start=0.0,
            stop=1.0,
            points=2,
        )


def test_conical_shock_adapter_matches_core() -> None:
    angle = np.deg2rad(10.0)
    adapted = conical_shock_condition(upstream_mach=2.0, cone_half_angle=angle)
    direct = conical_shock(2.0, angle)
    row = adapted.rows[0]
    assert row["shock_angle"] == pytest.approx(direct.shock_angle)
    assert row["surface_mach"] == pytest.approx(direct.surface_mach)
    assert row["maximum_cone_half_angle"] > angle


def test_conical_shock_sweeps_preserve_non_attached_points() -> None:
    angles = conical_shock_sweep(
        fixed_mach=2.0,
        fixed_cone_half_angle=np.deg2rad(10.0),
        sweep_field="cone_half_angle",
        start=0.0,
        stop=np.deg2rad(50.0),
        points=3,
    )
    assert angles.rows[0]["status"] == "ok"
    assert angles.rows[-1]["status"] == "no_attached_shock"
    mach = conical_shock_sweep(
        fixed_mach=2.0,
        fixed_cone_half_angle=np.deg2rad(10.0),
        sweep_field="mach",
        start=1.1,
        stop=3.0,
        points=3,
    )
    assert mach.rows[-1]["status"] == "ok"
    with pytest.raises(ValueError, match="sweep_field"):
        conical_shock_sweep(
            fixed_mach=2.0,
            fixed_cone_half_angle=0.1,
            sweep_field="other",
            start=0.0,
            stop=1.0,
            points=3,
        )


def test_flat_plate_adapter_and_warning_capture() -> None:
    laminar = flat_plate(
        distance=1.0,
        edge_velocity=10.0,
        edge_density=1.0,
        edge_dynamic_viscosity=1e-5,
        regime=BoundaryLayerRegime.LAMINAR,
        turbulent_correlation=TurbulentCorrelation.SCHLICHTING,
        transition_reynolds=None,
        compressibility_correction=CompressibilityCorrection.NONE,
        mach=None,
        edge_temperature=None,
        wall_temperature=None,
    )
    assert laminar.rows[0]["boundary_layer_thickness"] == pytest.approx(0.005)
    assert laminar.rows[0]["recovery_temperature"] is None
    warned = flat_plate(
        distance=1.0,
        edge_velocity=1.0,
        edge_density=1.0,
        edge_dynamic_viscosity=1e-5,
        regime=BoundaryLayerRegime.TURBULENT,
        turbulent_correlation=TurbulentCorrelation.SCHLICHTING,
        transition_reynolds=None,
        compressibility_correction=CompressibilityCorrection.NONE,
        mach=None,
        edge_temperature=None,
        wall_temperature=None,
    )
    assert warned.warnings


def test_flat_plate_compressible_distance_sweep() -> None:
    result = flat_plate_sweep(
        start=0.1,
        stop=1.0,
        points=3,
        logarithmic=True,
        edge_velocity=300.0,
        edge_density=0.5,
        edge_dynamic_viscosity=1.6e-5,
        regime=BoundaryLayerRegime.TRANSITIONAL,
        turbulent_correlation=TurbulentCorrelation.POWER_LAW,
        transition_reynolds=2e6,
        compressibility_correction=CompressibilityCorrection.ECKERT,
        mach=2.0,
        edge_temperature=250.0,
        wall_temperature=None,
    )
    assert len(result.rows) == 3
    assert result.rows[-1]["wall_temperature"] is not None
    assert result.rows[-1]["distance"] == pytest.approx(1.0)
