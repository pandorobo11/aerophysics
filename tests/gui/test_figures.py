"""Tests for GUI Plotly figures."""

import numpy as np
import pytest

from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    CompressibilityCorrection,
    TurbulentCorrelation,
)
from aerophysics.gui.adapters import (
    expansion_sweep,
    flat_plate_sweep,
    flight_sweep,
    isentropic_sweep,
    normal_shock_sweep,
    oblique_shock_condition,
    oblique_shock_sweep,
)
from aerophysics.gui.figures import (
    boundary_layer_figures,
    expansion_figures,
    flight_figures,
    isentropic_figures,
    normal_shock_figures,
    shock_geometry,
    shock_trends,
)
from aerophysics.gui.units import UnitPreferences
from aerophysics.isentropic import MachBranch
from aerophysics.shocks import ShockBranch


def test_flight_figures_have_expected_panels() -> None:
    result = flight_sweep(
        fixed_altitude=0.0,
        fixed_motion=0.8,
        motion_basis="mach",
        sweep_field="altitude",
        start=0.0,
        stop=10_000.0,
        points=3,
        characteristic_length=1.0,
    )
    figures = flight_figures(result.rows, UnitPreferences(length="ft"))
    assert set(figures) == {"大気状態", "飛行状態", "全状態量"}
    assert len(figures["大気状態"].data) == 4
    assert "ft" in str(figures["大気状態"].layout.xaxis.title.text)
    motion_figures = flight_figures(
        result.rows,
        UnitPreferences(),
        sweep_field="motion",
        motion_basis="velocity",
    )
    assert "速度" in str(motion_figures["飛行状態"].layout.xaxis.title.text)


def test_shock_geometry_and_both_sweep_axes() -> None:
    single = oblique_shock_condition(
        upstream_mach=2.0,
        deflection_angle=np.deg2rad(10.0),
        branch=ShockBranch.WEAK,
    )
    geometry = shock_geometry(single.rows[0], UnitPreferences())
    assert len(geometry.data) == 2
    assert "deg" in str(geometry.layout.title.text)
    bad_row = {**single.rows[0], "shock_angle": None}
    with pytest.raises(ValueError, match="successful"):
        shock_geometry(bad_row, UnitPreferences())

    theta = oblique_shock_sweep(
        fixed_mach=2.0,
        fixed_deflection=0.1,
        branch=ShockBranch.WEAK,
        sweep_field="deflection",
        start=0.0,
        stop=0.2,
        points=3,
    )
    theta_figures = shock_trends(theta.rows, UnitPreferences())
    assert "偏向角" in str(theta_figures["状態量"].layout.xaxis.title.text)
    mach = oblique_shock_sweep(
        fixed_mach=2.0,
        fixed_deflection=0.05,
        branch=ShockBranch.WEAK,
        sweep_field="mach",
        start=1.5,
        stop=3.0,
        points=3,
    )
    mach_figures = shock_trends(mach.rows, UnitPreferences())
    assert "Mach" in str(mach_figures["状態量"].layout.xaxis.title.text)


def test_boundary_layer_figures_include_transition_and_thermal() -> None:
    result = flat_plate_sweep(
        start=0.1,
        stop=1.0,
        points=3,
        logarithmic=False,
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
    figures = boundary_layer_figures(
        result.rows, UnitPreferences(), transition_distance=0.3
    )
    assert set(figures) == {"厚さ", "摩擦係数", "せん断・抗力", "温度"}
    assert len(figures["厚さ"].layout.shapes) == 1
    assert len(figures["温度"].data) == 2


def test_additional_compressible_flow_figures() -> None:
    isentropic = isentropic_sweep(
        input_basis="mach",
        branch=MachBranch.SUBSONIC,
        start=0.1,
        stop=3.0,
        points=3,
        total_pressure=101_325.0,
        total_temperature=300.0,
    )
    isentropic_plots = isentropic_figures(isentropic.rows, input_label="Mach M")
    assert set(isentropic_plots) == {"状態量比", "面積・流量", "質量流束"}
    assert len(isentropic_plots["状態量比"].data) == 3

    shock = normal_shock_sweep(start=1.0, stop=3.0, points=3)
    shock_plots = normal_shock_figures(shock.rows)
    assert set(shock_plots) == {"状態量", "全圧・ピトー"}
    assert len(shock_plots["状態量"].data) == 4

    expansion = expansion_sweep(
        fixed_mach=2.0,
        fixed_turn_angle=0.05,
        sweep_field="turn_angle",
        start=0.0,
        stop=0.2,
        points=3,
    )
    expansion_plots = expansion_figures(
        expansion.rows, UnitPreferences(), sweep_field="turn_angle"
    )
    assert set(expansion_plots) == {"Mach数", "角度", "状態量比"}
    assert "膨張角" in str(expansion_plots["Mach数"].layout.xaxis.title.text)
