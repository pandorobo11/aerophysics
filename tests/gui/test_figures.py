"""Tests for GUI Plotly figures."""

import numpy as np
import pytest

from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    CompressibilityCorrection,
    TurbulentCorrelation,
)
from aerophysics.boundary_layer_profile import (
    CompressibleVelocityTransformation,
    TemperatureVelocityRelation,
)
from aerophysics.detached_shock import DetachedShockGeometry
from aerophysics.gui.adapters import (
    conical_shock_condition,
    conical_shock_sweep,
    detached_shock_condition,
    detached_shock_shape,
    detached_shock_sweep,
    expansion_sweep,
    flat_plate_sweep,
    flight_sweep,
    isentropic_sweep,
    normal_shock_sweep,
    oblique_shock_condition,
    oblique_shock_sweep,
)
from aerophysics.gui.advanced_adapters import (
    boundary_layer_profiles,
    protrusion_sweep,
    thermochemistry_sweep,
    viscosity_sweep,
)
from aerophysics.gui.figures import (
    boundary_layer_figures,
    boundary_layer_profile_figures,
    conical_shock_geometry,
    conical_shock_trends,
    detached_shock_geometry,
    detached_shock_trends,
    expansion_figures,
    flight_figures,
    isentropic_figures,
    normal_shock_figures,
    protrusion_figures,
    protrusion_shape_figure,
    shock_geometry,
    shock_trends,
    thermochemistry_figures,
    viscosity_figures,
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
    figures = flight_figures(
        result.rows, UnitPreferences(length="ft", inverse_length="1/ft")
    )
    assert set(figures) == {"大気状態", "飛行状態", "全状態量"}
    assert len(figures["大気状態"].data) == 4
    assert "ft" in str(figures["大気状態"].layout.xaxis.title.text)
    assert "1/ft" in str(figures["飛行状態"].layout.yaxis4.title.text)
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


def test_conical_shock_geometry_and_sweep_axes() -> None:
    single = conical_shock_condition(
        upstream_mach=2.0, cone_half_angle=np.deg2rad(10.0)
    )
    geometry = conical_shock_geometry(single.rows[0], UnitPreferences())
    assert len(geometry.data) == 2
    assert "deg" in str(geometry.layout.title.text)
    with pytest.raises(ValueError, match="successful"):
        conical_shock_geometry(
            {**single.rows[0], "shock_angle": None}, UnitPreferences()
        )
    angles = conical_shock_sweep(
        fixed_mach=2.0,
        fixed_cone_half_angle=np.deg2rad(10.0),
        sweep_field="cone_half_angle",
        start=0.0,
        stop=np.deg2rad(20.0),
        points=3,
    )
    assert len(conical_shock_trends(angles.rows, UnitPreferences())["角度"].data) == 3
    mach = conical_shock_sweep(
        fixed_mach=2.0,
        fixed_cone_half_angle=np.deg2rad(10.0),
        sweep_field="mach",
        start=2.0,
        stop=3.0,
        points=3,
    )
    figures = conical_shock_trends(mach.rows, UnitPreferences())
    assert "Mach" in str(figures["状態量"].layout.xaxis.title.text)


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


def test_detached_shock_geometry_and_trends() -> None:
    shape = detached_shock_shape(
        upstream_mach=4.0,
        nose_radius=0.1,
        geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
    )
    geometry = detached_shock_geometry(shape, UnitPreferences(length="ft"))
    assert len(geometry.data) == 2
    assert geometry.layout.xaxis.scaleanchor == "y"
    assert "ft" in str(geometry.layout.xaxis.title.text)

    comparison = detached_shock_sweep(
        start=2.0,
        stop=4.0,
        points=3,
        nose_radius=0.1,
        geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
        selection="comparison",
    )
    figures = detached_shock_trends(comparison.rows, UnitPreferences())
    assert set(figures) == {"無次元離脱距離", "寸法・曲率", "モデル差"}
    assert len(figures["無次元離脱距離"].data) == 2

    cylinder = detached_shock_condition(
        upstream_mach=np.array([2.0, 4.0]),
        nose_radius=0.1,
        geometry=DetachedShockGeometry.CYLINDRICAL_NOSE_2D,
        selection="ambrosio_wortman",
    )
    cylinder_figures = detached_shock_trends(cylinder.rows, UnitPreferences())
    assert set(cylinder_figures) == {"無次元離脱距離", "寸法・曲率"}
    assert len(cylinder_figures["無次元離脱距離"].data) == 1

    seiff = detached_shock_sweep(
        start=2.0,
        stop=4.0,
        points=3,
        nose_radius=0.1,
        geometry=DetachedShockGeometry.AXISYMMETRIC_SPHERE,
        selection="seiff",
    )
    seiff_figures = detached_shock_trends(seiff.rows, UnitPreferences())
    assert len(seiff_figures["無次元離脱距離"].data) == 1
    assert seiff_figures["無次元離脱距離"].data[0].name == "Seiff"


def test_boundary_profile_and_protrusion_figures() -> None:
    profile = boundary_layer_profiles(
        edge_velocity=300.0,
        edge_density=1.0,
        edge_temperature=300.0,
        boundary_layer_thickness=0.05,
        wall_shear_stress=85.0,
        transformations=tuple(CompressibleVelocityTransformation),
        temperature_velocity_relation=(
            TemperatureVelocityRelation.GENERALIZED_REYNOLDS_ANALOGY
        ),
        wall_temperature=250.0,
        wake_parameter=None,
        points=51,
    )
    profile_plots = boundary_layer_profile_figures(
        profile.result.rows, UnitPreferences(length="ft")
    )
    assert set(profile_plots) == {"速度分布", "壁法則", "熱物性", "局所流れ"}
    assert len(profile_plots["速度分布"].data) == 2
    assert profile_plots["壁法則"].layout.xaxis.type == "log"

    sweep = protrusion_sweep(
        sweep_field="height",
        start=0.005,
        stop=0.02,
        points=3,
        drag_coefficient=1.0,
        height=0.01,
        base_width=0.005,
        shape="rectangle",
        edge_velocity=100.0,
        edge_density=1.0,
        boundary_layer_thickness=0.05,
    )
    trends = protrusion_figures(
        sweep.rows,
        UnitPreferences(length="ft", force="lbf"),
        sweep_field="height",
    )
    assert set(trends) == {"抗力・動圧", "遮蔽"}
    assert "ft" in str(trends["抗力・動圧"].layout.xaxis.title.text)
    assert "lbf" in str(trends["抗力・動圧"].layout.yaxis.title.text)
    expected_drag: list[float] = []
    for row in sweep.rows:
        direct_drag = row["direct_drag"]
        assert isinstance(direct_drag, float)
        expected_drag.append(direct_drag / 4.4482216152605)
    assert list(trends["抗力・動圧"].data[0].y) == pytest.approx(expected_drag)
    coefficient = protrusion_figures(
        sweep.rows, UnitPreferences(), sweep_field="drag_coefficient"
    )
    assert "抗力係数" in str(coefficient["遮蔽"].layout.xaxis.title.text)


@pytest.mark.parametrize("shape", ["rectangle", "triangle", "ellipse"])
def test_representative_protrusion_shape_figures(shape: str) -> None:
    figure = protrusion_shape_figure(
        height=0.01,
        base_width=0.005,
        boundary_layer_thickness=0.02,
        shape=shape,
        preferences=UnitPreferences(),
    )
    assert len(figure.data) == 1
    assert len(figure.layout.shapes) == 1


def test_csv_shape_and_thermochemistry_figures() -> None:
    shape = protrusion_shape_figure(
        height=0.01,
        base_width=0.005,
        boundary_layer_thickness=0.02,
        shape="csv",
        preferences=UnitPreferences(),
        shape_height=np.array([0.0, 0.01]),
        shape_width=np.array([0.005, 0.0]),
    )
    assert list(shape.data[0].x) == pytest.approx([0.005, 0.0])
    thermo = thermochemistry_sweep(
        start=200.0,
        stop=6000.0,
        points=3,
        pressure=101_325.0,
        reference_temperature=298.15,
        models=("NASA7", "NASA9"),
        allow_extrapolation=False,
    )
    figures = thermochemistry_figures(thermo.rows, UnitPreferences(temperature="°F"))
    assert set(figures) == {"比熱", "比熱比・音速", "エネルギー", "エントロピー"}
    assert len(figures["比熱"].data) == 4
    assert len(figures["比熱"].layout.shapes) == 2


def test_viscosity_figures_show_gaps_relative_difference_and_axis_scale() -> None:
    result = viscosity_sweep(
        start=79.0,
        stop=30_000.0,
        points=5,
        models=("Sutherland", "Keyes", "Blottner/Wilke"),
        allow_extrapolation=False,
        log_temperature=True,
    )
    figures = viscosity_figures(
        result.rows, UnitPreferences(temperature="°F"), log_temperature=True
    )
    assert set(figures) == {"粘性係数", "Sutherland基準相対差"}
    assert len(figures["粘性係数"].data) == 3
    assert figures["粘性係数"].layout.xaxis.type == "log"
    assert "K" in str(figures["粘性係数"].layout.xaxis.title.text)
    assert np.isnan(np.asarray(figures["粘性係数"].data[1].y, dtype=float)).any()
    assert np.isnan(np.asarray(figures["粘性係数"].data[2].y, dtype=float)).any()
    assert np.asarray(
        figures["Sutherland基準相対差"].data[0].y, dtype=float
    ) == pytest.approx(np.zeros(5))

    linear = viscosity_figures(
        result.rows, UnitPreferences(temperature="°F"), log_temperature=False
    )
    assert linear["粘性係数"].layout.xaxis.type == "linear"
    assert "°F" in str(linear["粘性係数"].layout.xaxis.title.text)
