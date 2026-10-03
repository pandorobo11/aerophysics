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
    Row,
    detached_shock_condition,
    detached_shock_shape,
    detached_shock_sweep,
    expansion_condition,
    expansion_sweep,
    flat_plate_sweep,
    flight_sweep,
    isentropic_sweep,
    normal_shock_condition,
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
    expansion_geometry,
    flight_figures,
    isentropic_figures,
    normal_shock_figures,
    normal_shock_geometry,
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


@pytest.mark.parametrize(
    ("gas_model", "mach"),
    [
        ("AIR", 1.0),
        ("AIR", 3.0),
        ("NASA7", 3.0),
        ("NASA9", 3.0),
        ("HARMONIC_OSCILLATOR", 3.0),
        ("HARMONIC_OSCILLATOR", 4.2),
    ],
)
def test_normal_geometry_uses_solved_mach_and_collinear_flow(
    gas_model: str, mach: float
) -> None:
    result = normal_shock_condition(
        upstream_mach=mach,
        gas_model=gas_model,
        upstream_temperature=None if gas_model == "AIR" else 500.0,
    )
    row = result.rows[0]
    figure = normal_shock_geometry(row)
    assert figure.layout.yaxis.scaleanchor == "x"
    assert figure.layout.yaxis.scaleratio == 1
    shock = figure.data[0]
    assert list(shock.x) == [0.0, 0.0]
    assert shock.y[0] < 0.0 < shock.y[1]
    assert shock.line.color == "#d62728"
    arrows = [item for item in figure.layout.annotations if item.showarrow]
    assert len(arrows) == 2
    for arrow in arrows:
        assert arrow.x > arrow.ax
        assert arrow.y == arrow.ay == 0.0
        assert arrow.xref == arrow.axref == "x"
        assert arrow.yref == arrow.ayref == "y"
    assert arrows[0].x < 0.0 < arrows[1].ax
    text = " ".join(str(item.text) for item in figure.layout.annotations)
    assert f"M₁ = {mach:.3g}" in text
    assert f"M₂ = {row['downstream_mach']:.3g}" in text
    if mach == 1.0:
        assert "衝撃波強度ゼロ" in text
    if gas_model == "HARMONIC_OSCILLATOR" and mach == 4.2:
        assert row["pitot_pressure_ratio"] is None


@pytest.mark.parametrize("unit", ["deg", "rad"])
@pytest.mark.parametrize(
    ("mach", "theta_degrees"),
    [(1.0, 15.0), (2.0, 0.0), (2.0, 15.0), (2.0, 80.0), (2.0, 100.0), (1.05, 110.0)],
)
def test_expansion_geometry_fan_and_local_mach_angles(
    mach: float, theta_degrees: float, unit: str
) -> None:
    theta = float(np.deg2rad(theta_degrees))
    row = expansion_condition(upstream_mach=mach, turn_angle=theta).rows[0]
    figure = expansion_geometry(row, UnitPreferences(angle=unit))
    assert figure.layout.yaxis.scaleanchor == "x"
    assert figure.layout.yaxis.scaleratio == 1
    mu1 = float(np.arcsin(1.0 / mach))
    mach2 = row["downstream_mach"]
    assert isinstance(mach2, float)
    mu2 = float(np.arcsin(1.0 / mach2))
    wall = figure.data[0]
    assert np.arctan2(wall.y[-1], wall.x[-1]) == pytest.approx(-theta)
    fan = [trace for trace in figure.data if "マッハ線" in trace.name]
    assert len(fan) == (9 if theta > 0.0 else 1)
    assert np.arctan2(fan[0].y[-1], fan[0].x[-1]) == pytest.approx(mu1)
    assert np.arctan2(fan[-1].y[-1], fan[-1].x[-1]) == pytest.approx(mu2 - theta)
    for ray in fan:
        assert ray.x[0] == ray.y[0] == 0.0
        assert np.hypot(ray.x[-1], ray.y[-1]) == pytest.approx(1.6)
        assert ray.line.color == "#16836c"
    arcs = figure.data[-3:]
    for arc, start, end in zip(
        arcs, [0.0, 0.0, -theta], [-theta, mu1, mu2 - theta], strict=True
    ):
        assert np.arctan2(arc.y[0], arc.x[0]) == pytest.approx(start)
        assert np.arctan2(arc.y[-1], arc.x[-1]) == pytest.approx(end)
        value = abs(end - start) if unit == "rad" else np.rad2deg(abs(end - start))
        assert f"{value:.3g} {unit}" in arc.name
    arrows = [item for item in figure.layout.annotations if item.showarrow]
    assert len(arrows) == 2
    assert arrows[0].x > arrows[0].ax and arrows[0].y == arrows[0].ay
    downstream = arrows[1]
    assert np.arctan2(
        downstream.y - downstream.ay, downstream.x - downstream.ax
    ) == pytest.approx(-theta)
    assert np.cos(theta) * downstream.ay + np.sin(theta) * downstream.ax > 0.0
    # The segment must lie in uniform downstream flow, above the wall and
    # below the final Mach line; its direction alone does not establish that.
    points = np.linspace(
        [downstream.ax, downstream.ay], [downstream.x, downstream.y], 17
    )
    wall_normal = np.array([np.sin(theta), np.cos(theta)])
    final_angle = np.arctan2(fan[-1].y[-1], fan[-1].x[-1])
    final_normal = np.array([-np.sin(final_angle), np.cos(final_angle)])
    assert np.all(points @ wall_normal > 0.0)
    assert np.all(points @ final_normal < 0.0)
    for trace in figure.data:
        assert np.all(np.isfinite(trace.x)) and np.all(np.isfinite(trace.y))
    if theta == 0.0:
        assert "膨張なし" in fan[0].name


@pytest.mark.parametrize(
    ("kind", "key"),
    [
        ("normal", "upstream_mach"),
        ("normal", "downstream_mach"),
        ("expansion", "upstream_mach"),
        ("expansion", "downstream_mach"),
        ("expansion", "turn_angle"),
    ],
)
def test_new_flow_geometry_requires_successful_scalar_states(
    kind: str, key: str
) -> None:
    row: Row = {"upstream_mach": 2.0, "downstream_mach": 2.6, "turn_angle": 0.2}
    row[key] = None
    with pytest.raises(ValueError, match="successful"):
        if kind == "normal":
            normal_shock_geometry(row)
        else:
            expansion_geometry(row, UnitPreferences())


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
    assert len(geometry.data) == 4
    assert "deg" in str(geometry.layout.title.text)
    text = " ".join(str(item.text) for item in geometry.layout.annotations)
    assert "上流 M₁ = 2" in text
    assert f"下流 M₂ = {single.rows[0]['downstream_mach']:.3g}" in text
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
    # Deliberate plotting inputs: physical cone solutions are exercised by AppTest.
    rows: tuple[Row, ...] = (
        {
            "upstream_mach": 2.0,
            "cone_half_angle": np.deg2rad(10.0),
            "shock_angle": np.deg2rad(31.0),
            "maximum_cone_half_angle": np.deg2rad(40.0),
            "post_shock_mach": 1.9,
            "surface_mach": 1.8,
        },
        {
            "upstream_mach": 2.0,
            "cone_half_angle": np.deg2rad(20.0),
            "shock_angle": np.deg2rad(38.0),
            "maximum_cone_half_angle": np.deg2rad(40.0),
            "surface_mach": 1.6,
        },
    )
    geometry = conical_shock_geometry(rows[0], UnitPreferences())
    assert len(geometry.data) == 4
    assert "deg" in str(geometry.layout.title.text)
    text = " ".join(str(item.text) for item in geometry.layout.annotations)
    assert "上流 M∞ = 2" in text
    assert "表面 Mₛ = 1.8" in text
    assert "表面 Mₛ = 1.9" not in text
    with pytest.raises(ValueError, match="successful"):
        conical_shock_geometry({**rows[0], "shock_angle": None}, UnitPreferences())
    angles = conical_shock_trends(rows, UnitPreferences())
    assert list(angles["角度"].data[0].x) == pytest.approx([10.0, 20.0])
    assert list(angles["角度"].data[0].y) == pytest.approx([31.0, 38.0])
    mach_rows = (
        rows[0],
        {
            **rows[1],
            "upstream_mach": 3.0,
            "cone_half_angle": rows[0]["cone_half_angle"],
        },
    )
    figures = conical_shock_trends(mach_rows, UnitPreferences())
    assert "Mach" in str(figures["状態量"].layout.xaxis.title.text)
    assert list(figures["状態量"].data[0].x) == pytest.approx([2.0, 3.0])


@pytest.mark.parametrize("unit", ["deg", "rad"])
@pytest.mark.parametrize(
    ("conical", "theta_degrees", "beta_degrees"),
    [
        (False, 0.0, 90.0),
        (False, 10.0, 27.3),
        (False, 10.0, 86.8),
        (True, 10.0, 21.7),
        (True, 40.0, 65.0),
    ],
)
def test_attached_geometry_preserves_angles_and_flow_directions(
    conical: bool, theta_degrees: float, beta_degrees: float, unit: str
) -> None:
    theta, beta = np.deg2rad([theta_degrees, beta_degrees])
    angle_key = "cone_half_angle" if conical else "deflection_angle"
    plot = conical_shock_geometry if conical else shock_geometry
    figure = plot(
        {
            angle_key: float(theta),
            "shock_angle": float(beta),
            "upstream_mach": 3.0,
            "surface_mach" if conical else "downstream_mach": 1.81234,
        },
        UnitPreferences(angle=unit),
    )
    assert figure.layout.yaxis.scaleanchor == "x"
    assert figure.layout.yaxis.scaleratio == 1
    if conical:
        assert "円錐軸" in str(figure.layout.annotations[-1].text)

    # Wall/shock rays and the two angle arcs must use the same physical angles,
    # independently of display units, including the vertical strong-shock limit.
    for trace, angle in zip(figure.data, (theta, beta, theta, beta), strict=True):
        x, y = np.asarray(trace.x), np.asarray(trace.y)
        assert np.all(np.isfinite(x)) and np.all(np.isfinite(y))
        assert max(np.max(np.abs(x)), np.max(np.abs(y))) <= 2.0
        assert np.arctan2(y[-1], x[-1]) == pytest.approx(angle, abs=1e-14)
    theta_display = theta_degrees if unit == "deg" else theta
    assert f"{theta_display:.3g} {unit}" in figure.data[2].name

    arrows = [item for item in figure.layout.annotations if item.showarrow]
    assert len(arrows) == 2
    for arrow, angle in zip(arrows, (0.0, theta), strict=True):
        assert arrow.axref == arrow.xref == "x"
        assert arrow.ayref == arrow.yref == "y"
        assert arrow.arrowhead > 0
        assert np.arctan2(arrow.y - arrow.ay, arrow.x - arrow.ax) == pytest.approx(
            angle, abs=1e-14
        )
    text = " ".join(str(item.text) for item in figure.layout.annotations)
    assert ("上流 M∞ = 3" if conical else "上流 M₁ = 3") in text
    assert ("表面 Mₛ = 1.81" if conical else "下流 M₂ = 1.81") in text
    assert "角度の基準" in text or "円錐軸" in text


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


@pytest.mark.parametrize("geometry", list(DetachedShockGeometry))
@pytest.mark.parametrize(
    ("unit", "factor"), [("m", 1.0), ("mm", 1000.0), ("ft", 1.0 / 0.3048)]
)
def test_detached_geometry_dimensions_and_upstream_coordinates(
    geometry: DetachedShockGeometry, unit: str, factor: float
) -> None:
    shape = detached_shock_shape(upstream_mach=4.0, nose_radius=0.1, geometry=geometry)
    figure = detached_shock_geometry(shape, UnitPreferences(length=unit))
    radius = 0.1 * factor
    distance = float(shape.standoff_distance) * factor
    assert np.asarray(figure.data[1].x) == pytest.approx(shape.shock_x * factor)
    assert np.asarray(figure.data[1].y) == pytest.approx(shape.shock_y * factor)
    assert figure.layout.xaxis.scaleratio == 1
    assert "Billig" in figure.layout.title.text
    assert "Ambrosio" in figure.layout.title.text
    assert "+x" in figure.layout.xaxis.title.text
    assert "上流" in figure.layout.xaxis.title.text
    assert "左側" in figure.layout.xaxis.title.text
    x_left, x_right = figure.layout.xaxis.range
    assert x_left > x_right

    bracket = next(
        line
        for line in figure.layout.shapes
        if line.type == "line" and line.y0 == line.y1 and line.y0 < 0.0
    )
    assert bracket.x0 == pytest.approx(radius)
    assert bracket.x1 == pytest.approx(radius + distance)
    assert bracket.x1 - bracket.x0 == pytest.approx(distance)
    extensions = [
        line
        for line in figure.layout.shapes
        if line.type == "line" and line.x0 == line.x1
    ]
    assert [line.x0 for line in extensions] == pytest.approx(
        [radius, radius + distance]
    )
    assert all(line.y0 == 0.0 for line in extensions)

    arrows = [item for item in figure.layout.annotations if item.showarrow]
    assert len(arrows) == 2
    radius_arrow = next(item for item in arrows if item.arrowcolor == "#555")
    assert radius_arrow.ax == radius_arrow.ay == 0.0
    assert np.hypot(radius_arrow.x, radius_arrow.y) == pytest.approx(radius)
    flow_arrow = next(item for item in arrows if item.arrowcolor == "#2463a5")
    assert flow_arrow.x < flow_arrow.ax  # Freestream velocity points toward -x.
    assert flow_arrow.y == flow_arrow.ay
    assert flow_arrow.x > radius + distance  # Arrow stays upstream of the shock.
    assert flow_arrow.axref == flow_arrow.xref == "x"
    assert flow_arrow.ayref == flow_arrow.yref == "y"
    assert flow_arrow.arrowhead > 0
    # Rendering reverses x only: the arrow runs left-to-right, and both the
    # arrow and shock vertex stay upstream (left) of the body vertex.
    tail_fraction = (flow_arrow.ax - x_left) / (x_right - x_left)
    tip_fraction = (flow_arrow.x - x_left) / (x_right - x_left)
    shock_fraction = (radius + distance - x_left) / (x_right - x_left)
    body_fraction = (radius - x_left) / (x_right - x_left)
    assert 0.0 <= tail_fraction < tip_fraction < shock_fraction < body_fraction <= 1.0
    text = " ".join(str(item.text) for item in figure.layout.annotations)
    assert f"Rₙ = {radius:.3g} {unit}" in text
    assert f"Δ = {distance:.3g} {unit}" in text
    assert "頭部曲率中心" in text and "流れ方向" in text
    assert (
        "半球頭部"
        if geometry is DetachedShockGeometry.AXISYMMETRIC_SPHERE
        else "円柱頭部"
    ) in text


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


@pytest.mark.parametrize(
    ("shape", "expected_widths"),
    [
        ("rectangle", [0.005, 0.005, 0.005]),
        ("triangle", [0.005, 0.0025, 0.0]),
        ("ellipse", [0.005, 0.004330127018922193, 0.0]),
    ],
)
def test_representative_protrusion_shapes_and_boundary_edge(
    shape: str, expected_widths: list[float]
) -> None:
    figure = protrusion_shape_figure(
        height=0.01,
        base_width=0.005,
        boundary_layer_thickness=0.02,
        shape=shape,
        preferences=UnitPreferences(),
    )
    width = np.asarray(figure.data[0].x, dtype=float)
    height = np.asarray(figure.data[0].y, dtype=float)
    assert np.interp([0.0, 0.005, 0.01], height, width) == pytest.approx(
        expected_widths
    )
    assert figure.layout.shapes[0].y0 == pytest.approx(0.02)
    assert figure.layout.shapes[0].y1 == pytest.approx(0.02)


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
