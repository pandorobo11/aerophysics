"""Plotly figures for GUI calculation results."""

from __future__ import annotations

import math

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from aerophysics.detached_shock import (
    BilligShockShapeResult,
    DetachedShockGeometry,
)
from aerophysics.gui.adapters import Row
from aerophysics.gui.units import UnitPreferences, from_si, selected_unit


def _numeric(rows: tuple[Row, ...], key: str) -> np.ndarray:
    return np.asarray(
        [np.nan if not isinstance(row.get(key), float) else row[key] for row in rows],
        dtype=np.float64,
    )


def _converted(
    rows: tuple[Row, ...], key: str, kind: str, preferences: UnitPreferences
) -> np.ndarray:
    unit = selected_unit(kind, preferences)  # type: ignore[arg-type]
    return np.asarray(from_si(_numeric(rows, key), kind, unit), dtype=np.float64)  # type: ignore[arg-type]


def _style(figure: go.Figure, title: str) -> go.Figure:
    figure.update_layout(
        title=title,
        template="plotly_white",
        height=560,
        margin={"l": 60, "r": 30, "t": 70, "b": 50},
        hovermode="x unified",
        legend={"orientation": "h", "y": 1.08, "x": 0.0},
    )
    return figure


def flight_figures(
    rows: tuple[Row, ...],
    preferences: UnitPreferences,
    *,
    sweep_field: str = "altitude",
    motion_basis: str = "mach",
) -> dict[str, go.Figure]:
    """Create atmosphere, flight, and total-state figures."""
    if sweep_field == "motion" and motion_basis == "mach":
        x = _numeric(rows, "mach")
        x_title = "Mach M"
    elif sweep_field == "motion" and motion_basis == "velocity":
        x = _converted(rows, "velocity", "speed", preferences)
        x_title = f"速度 V [{preferences.speed}]"
    else:
        x = _converted(rows, "geometric_altitude", "length", preferences)
        x_title = f"h [{preferences.length}]"

    atmosphere = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=("静温 T", "静圧 p", "密度 ρ", "音速 a"),
    )
    quantities = (
        ("temperature", "temperature", preferences.temperature),
        ("pressure", "pressure", preferences.pressure),
        ("density", "density", preferences.density),
        ("speed_of_sound", "speed", preferences.speed),
    )
    for index, (key, kind, unit) in enumerate(quantities):
        row, col = divmod(index, 2)
        atmosphere.add_trace(
            go.Scatter(x=x, y=_converted(rows, key, kind, preferences), name=key),
            row=row + 1,
            col=col + 1,
        )
        atmosphere.update_yaxes(title_text=unit, row=row + 1, col=col + 1)
        atmosphere.update_xaxes(title_text=x_title, row=row + 1, col=col + 1)

    flight = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=("Mach M", "速度 V", "動圧 q", "Reynolds数/長さ"),
    )
    flight_values = (
        ("mach", _numeric(rows, "mach"), "–"),
        (
            "velocity",
            _converted(rows, "velocity", "speed", preferences),
            preferences.speed,
        ),
        (
            "dynamic_pressure",
            _converted(rows, "dynamic_pressure", "pressure", preferences),
            preferences.pressure,
        ),
        (
            "reynolds_number_per_length",
            _numeric(rows, "reynolds_number_per_length"),
            "1/m",
        ),
    )
    for index, (name, values, unit) in enumerate(flight_values):
        row, col = divmod(index, 2)
        flight.add_trace(go.Scatter(x=x, y=values, name=name), row=row + 1, col=col + 1)
        flight.update_yaxes(title_text=unit, row=row + 1, col=col + 1)
        flight.update_xaxes(title_text=x_title, row=row + 1, col=col + 1)

    total = make_subplots(
        rows=1, cols=3, subplot_titles=("全温 T₀", "全圧 p₀", "全密度 ρ₀")
    )
    totals = (
        ("total_temperature", "temperature", preferences.temperature),
        ("total_pressure", "pressure", preferences.pressure),
        ("total_density", "density", preferences.density),
    )
    for index, (key, kind, unit) in enumerate(totals):
        total.add_trace(
            go.Scatter(x=x, y=_converted(rows, key, kind, preferences), name=key),
            row=1,
            col=index + 1,
        )
        total.update_xaxes(title_text=x_title, row=1, col=index + 1)
        total.update_yaxes(title_text=unit, row=1, col=index + 1)
    return {
        "大気状態": _style(atmosphere, "標準大気"),
        "飛行状態": _style(flight, "飛行条件"),
        "全状態量": _style(total, "全状態量"),
    }


def shock_geometry(row: Row, preferences: UnitPreferences) -> go.Figure:
    """Create a schematic wedge and attached-shock diagram."""
    theta = row.get("deflection_angle")
    beta = row.get("shock_angle")
    if not isinstance(theta, float) or not isinstance(beta, float):
        raise ValueError("geometry requires a successful shock result")
    length = 1.0
    wedge_y = math.tan(theta)
    shock_y = math.tan(beta)
    figure = go.Figure()
    figure.add_trace(
        go.Scatter(
            x=[0.0, length],
            y=[0.0, wedge_y],
            mode="lines",
            line={"width": 8, "color": "#555"},
            name="くさび面",
        )
    )
    figure.add_trace(
        go.Scatter(
            x=[0.0, length],
            y=[0.0, shock_y],
            mode="lines",
            line={"width": 4, "color": "#d62728"},
            name="衝撃波",
        )
    )
    figure.add_annotation(x=0.45, y=-0.08, text="M₁", showarrow=True, ax=-70, ay=0)
    figure.add_annotation(
        x=0.62,
        y=0.62 * math.tan(0.5 * (theta + beta)),
        text="M₂",
        showarrow=True,
        ax=-55,
        ay=30,
    )
    unit = preferences.angle
    theta_display = float(from_si(theta, "angle", unit))
    beta_display = float(from_si(beta, "angle", unit))
    figure.update_layout(
        title=(
            f"付着衝撃波模式図 — θ={theta_display:.3g} {unit}, "
            f"β={beta_display:.3g} {unit}"
        ),
        template="plotly_white",
        height=430,
        xaxis={"visible": False, "range": [-0.15, 1.1]},
        yaxis={
            "visible": False,
            "scaleanchor": "x",
            "scaleratio": 1,
            "range": [-0.15, min(max(shock_y * 1.1, 0.5), 5.0)],
        },
        margin={"l": 20, "r": 20, "t": 70, "b": 20},
    )
    return figure


def shock_trends(
    rows: tuple[Row, ...], preferences: UnitPreferences
) -> dict[str, go.Figure]:
    """Create shock-angle/Mach and state-ratio trend figures."""
    changing_theta = np.nanmax(_numeric(rows, "deflection_angle")) != np.nanmin(
        _numeric(rows, "deflection_angle")
    )
    if changing_theta:
        x = _converted(rows, "deflection_angle", "angle", preferences)
        x_title = f"偏向角 θ [{preferences.angle}]"
    else:
        x = _numeric(rows, "upstream_mach")
        x_title = "上流 Mach M₁"
    state = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=("下流 Mach M₂", "p₂/p₁", "ρ₂/ρ₁", "T₂/T₁・p₀₂/p₀₁"),
    )
    series = (
        ("downstream_mach", "M₂", 1, 1),
        ("static_pressure_ratio", "p₂/p₁", 1, 2),
        ("static_density_ratio", "ρ₂/ρ₁", 2, 1),
        ("static_temperature_ratio", "T₂/T₁", 2, 2),
        ("total_pressure_ratio", "p₀₂/p₀₁", 2, 2),
    )
    for key, name, row, col in series:
        state.add_trace(
            go.Scatter(x=x, y=_numeric(rows, key), name=name), row=row, col=col
        )
    state.update_xaxes(title_text=x_title)

    angles = go.Figure()
    angles.add_trace(
        go.Scatter(
            x=x,
            y=_converted(rows, "shock_angle", "angle", preferences),
            name="衝撃波角 β",
        )
    )
    angles.add_trace(
        go.Scatter(
            x=x,
            y=_converted(rows, "maximum_deflection_angle", "angle", preferences),
            name="最大付着偏向角 θmax",
            line={"dash": "dash"},
        )
    )
    angles.update_xaxes(title_text=x_title)
    angles.update_yaxes(title_text=f"角度 [{preferences.angle}]")
    return {
        "状態量": _style(state, "斜め衝撃波の状態量"),
        "角度": _style(angles, "付着衝撃波の角度"),
    }


def conical_shock_geometry(row: Row, preferences: UnitPreferences) -> go.Figure:
    """Create a meridional schematic of a cone and its attached shock."""
    cone_angle = row.get("cone_half_angle")
    shock_angle = row.get("shock_angle")
    if not isinstance(cone_angle, float) or not isinstance(shock_angle, float):
        raise ValueError("geometry requires a successful conical-shock result")
    length = 1.0
    cone_y = math.tan(cone_angle)
    shock_y = math.tan(shock_angle)
    figure = go.Figure()
    figure.add_trace(
        go.Scatter(
            x=[0.0, length],
            y=[0.0, cone_y],
            mode="lines",
            line={"width": 8, "color": "#555"},
            name="円錐面",
        )
    )
    figure.add_trace(
        go.Scatter(
            x=[0.0, length],
            y=[0.0, shock_y],
            mode="lines",
            line={"width": 4, "color": "#d62728"},
            name="円錐衝撃波",
        )
    )
    figure.add_annotation(x=0.45, y=-0.08, text="M∞", showarrow=True, ax=-70, ay=0)
    figure.add_annotation(
        x=0.68,
        y=0.68 * cone_y,
        text="Mₛ",
        showarrow=True,
        ax=-45,
        ay=30,
    )
    unit = preferences.angle
    cone_display = float(from_si(cone_angle, "angle", unit))
    shock_display = float(from_si(shock_angle, "angle", unit))
    figure.update_layout(
        title=(
            f"円錐衝撃波模式図 — θc={cone_display:.3g} {unit}, "
            f"β={shock_display:.3g} {unit}"
        ),
        template="plotly_white",
        height=430,
        xaxis={"visible": False, "range": [-0.15, 1.1]},
        yaxis={
            "visible": False,
            "scaleanchor": "x",
            "scaleratio": 1,
            "range": [-0.15, min(max(shock_y * 1.1, 0.5), 5.0)],
        },
        margin={"l": 20, "r": 20, "t": 70, "b": 20},
    )
    return figure


def conical_shock_trends(
    rows: tuple[Row, ...], preferences: UnitPreferences
) -> dict[str, go.Figure]:
    """Create cone-shock angle, Mach, and surface-state trend figures."""
    cone_angles = _numeric(rows, "cone_half_angle")
    changing_angle = np.nanmax(cone_angles) != np.nanmin(cone_angles)
    if changing_angle:
        x = _converted(rows, "cone_half_angle", "angle", preferences)
        x_title = f"円錐半頂角 θc [{preferences.angle}]"
    else:
        x = _numeric(rows, "upstream_mach")
        x_title = "上流 Mach M∞"

    state = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=(
            "表面・衝撃波直後Mach",
            "pₛ/p∞",
            "ρₛ/ρ∞",
            "Tₛ/T∞・p₀₂/p₀∞",
        ),
    )
    for key, name in (
        ("surface_mach", "表面 Mach Mₛ"),
        ("post_shock_mach", "衝撃波直後 Mach M₂"),
    ):
        state.add_trace(go.Scatter(x=x, y=_numeric(rows, key), name=name), row=1, col=1)
    for key, name, row, col in (
        ("surface_pressure_ratio", "pₛ/p∞", 1, 2),
        ("surface_density_ratio", "ρₛ/ρ∞", 2, 1),
        ("surface_temperature_ratio", "Tₛ/T∞", 2, 2),
        ("total_pressure_ratio", "p₀₂/p₀∞", 2, 2),
    ):
        state.add_trace(
            go.Scatter(x=x, y=_numeric(rows, key), name=name), row=row, col=col
        )
    state.update_xaxes(title_text=x_title)

    angles = go.Figure()
    for key, name, dash in (
        ("shock_angle", "衝撃波角 β", None),
        ("cone_half_angle", "円錐半頂角 θc", "dot"),
        ("maximum_cone_half_angle", "最大付着半頂角 θc,max", "dash"),
    ):
        line = {} if dash is None else {"dash": dash}
        angles.add_trace(
            go.Scatter(
                x=x,
                y=_converted(rows, key, "angle", preferences),
                name=name,
                line=line,
            )
        )
    angles.update_xaxes(title_text=x_title)
    angles.update_yaxes(title_text=f"角度 [{preferences.angle}]")
    return {
        "状態量": _style(state, "円錐表面の状態量"),
        "角度": _style(angles, "円錐衝撃波の角度"),
    }


def detached_shock_geometry(
    shape: BilligShockShapeResult,
    preferences: UnitPreferences,
) -> go.Figure:
    """Plot a blunt nose and its Billig detached-shock shape."""
    if not isinstance(shape.nose_radius, float):
        raise ValueError("geometry figure requires a scalar Billig result")
    radius = shape.nose_radius
    body_y = np.linspace(radius, -radius, 181, dtype=np.float64)
    body_x = np.sqrt(np.maximum(radius**2 - body_y**2, 0.0))
    afterbody_x = -2.0 * radius
    outline_x = np.concatenate(([-2.0 * radius, 0.0], body_x, [afterbody_x]))
    outline_y = np.concatenate(([radius, radius], body_y, [-radius]))
    unit = preferences.length

    figure = go.Figure()
    figure.add_trace(
        go.Scatter(
            x=from_si(outline_x, "length", unit),
            y=from_si(outline_y, "length", unit),
            mode="lines",
            line={"width": 5, "color": "#555"},
            fill="toself",
            fillcolor="rgba(100,100,100,0.12)",
            name=(
                "半球頭部"
                if shape.geometry is DetachedShockGeometry.AXISYMMETRIC_SPHERE
                else "2D円柱頭部"
            ),
        )
    )
    figure.add_trace(
        go.Scatter(
            x=from_si(shape.shock_x, "length", unit),
            y=from_si(shape.shock_y, "length", unit),
            mode="lines",
            line={"width": 4, "color": "#d62728"},
            name="Billig衝撃波",
        )
    )
    figure.add_annotation(
        x=float(from_si(1.6 * radius, "length", unit)),
        y=float(from_si(1.7 * radius, "length", unit)),
        text="M∞",
        showarrow=True,
        ax=60,
        ay=0,
    )
    figure.update_layout(
        title="離脱衝撃波形状（Billig、離脱距離はAmbrosio–Wortman）",
        template="plotly_white",
        height=560,
        xaxis={"title": f"x [{unit}]", "scaleanchor": "y", "scaleratio": 1},
        yaxis={"title": f"y [{unit}]"},
        margin={"l": 60, "r": 30, "t": 70, "b": 55},
        legend={"orientation": "h", "y": 1.08, "x": 0.0},
    )
    return figure


def detached_shock_trends(
    rows: tuple[Row, ...], preferences: UnitPreferences
) -> dict[str, go.Figure]:
    """Create Mach-sweep standoff, curvature, and comparison figures."""
    mach = _numeric(rows, "upstream_mach")
    selection = str(rows[0].get("selection", "comparison"))
    show_aw = selection in {"ambrosio_wortman", "comparison"}
    show_seiff = selection in {"seiff", "comparison"}
    normalized = go.Figure()
    if show_aw:
        normalized.add_trace(
            go.Scatter(
                x=mach,
                y=_numeric(rows, "aw_normalized_standoff_distance"),
                name="Ambrosio–Wortman",
            )
        )
    seiff_normalized = _numeric(rows, "seiff_normalized_standoff_distance")
    if show_seiff and np.any(np.isfinite(seiff_normalized)):
        normalized.add_trace(go.Scatter(x=mach, y=seiff_normalized, name="Seiff"))
    normalized.update_xaxes(title_text="Mach M∞")
    normalized.update_yaxes(title_text="Δ/Rn")

    dimensional = go.Figure()
    if show_aw:
        dimensional.add_trace(
            go.Scatter(
                x=mach,
                y=_converted(rows, "aw_standoff_distance", "length", preferences),
                name="AW Δ",
            )
        )
    seiff_distance = _numeric(rows, "seiff_standoff_distance")
    if show_seiff and np.any(np.isfinite(seiff_distance)):
        dimensional.add_trace(
            go.Scatter(
                x=mach,
                y=_converted(rows, "seiff_standoff_distance", "length", preferences),
                name="Seiff Δ",
            )
        )
    dimensional.add_trace(
        go.Scatter(
            x=mach,
            y=_converted(
                rows,
                "billig_vertex_curvature_radius",
                "length",
                preferences,
            ),
            name="Billig Rc",
            line={"dash": "dash"},
        )
    )
    dimensional.update_xaxes(title_text="Mach M∞")
    dimensional.update_yaxes(title_text=f"距離 [{preferences.length}]")

    figures = {
        "無次元離脱距離": _style(normalized, "離脱距離のMach依存性"),
        "寸法・曲率": _style(dimensional, "離脱距離と衝撃波頂点曲率半径"),
    }
    relative = _numeric(rows, "relative_difference")
    if selection == "comparison" and np.any(np.isfinite(relative)):
        comparison = go.Figure()
        comparison.add_trace(
            go.Scatter(x=mach, y=100.0 * relative, name="(Seiff−AW)/AW")
        )
        comparison.update_xaxes(title_text="Mach M∞")
        comparison.update_yaxes(title_text="相対差 [%]")
        figures["モデル差"] = _style(comparison, "SeiffとAmbrosio–Wortmanの差")
    return figures


def isentropic_figures(
    rows: tuple[Row, ...], *, input_label: str
) -> dict[str, go.Figure]:
    """Create isentropic ratio and mass-flow figures."""
    x = _numeric(rows, "input_value")
    ratios = go.Figure()
    for key, name in (
        ("total_temperature_ratio", "T₀/T"),
        ("total_pressure_ratio", "p₀/p"),
        ("total_density_ratio", "ρ₀/ρ"),
    ):
        ratios.add_trace(go.Scatter(x=x, y=_numeric(rows, key), name=name))
    ratios.update_xaxes(title_text=input_label)
    ratios.update_yaxes(title_text="全量/静量")

    flow = make_subplots(specs=[[{"secondary_y": True}]])
    flow.add_trace(
        go.Scatter(x=x, y=_numeric(rows, "area_ratio"), name="A/A*"),
        secondary_y=False,
    )
    flow.add_trace(
        go.Scatter(
            x=x,
            y=_numeric(rows, "mass_flow_parameter"),
            name="質量流量パラメータ",
        ),
        secondary_y=True,
    )
    flow.update_xaxes(title_text=input_label)
    flow.update_yaxes(title_text="A/A*", secondary_y=False)
    flow.update_yaxes(title_text="質量流量パラメータ", secondary_y=True)
    figures = {
        "状態量比": _style(ratios, "等エントロピー状態量比"),
        "面積・流量": _style(flow, "面積-Mach関係と質量流量"),
    }
    if np.any(np.isfinite(_numeric(rows, "mass_flux"))):
        flux = go.Figure()
        flux.add_trace(go.Scatter(x=x, y=_numeric(rows, "mass_flux"), name="質量流束"))
        flux.add_trace(
            go.Scatter(
                x=x,
                y=_numeric(rows, "choked_mass_flux"),
                name="チョーク質量流束",
                line={"dash": "dash"},
            )
        )
        flux.update_xaxes(title_text=input_label)
        flux.update_yaxes(title_text="kg/(m²·s)")
        figures["質量流束"] = _style(flux, "質量流束")
    return figures


def normal_shock_figures(rows: tuple[Row, ...]) -> dict[str, go.Figure]:
    """Create normal-shock Mach and pressure-ratio figures."""
    x = _numeric(rows, "upstream_mach")
    states = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=("下流 Mach M₂", "p₂/p₁", "ρ₂/ρ₁", "T₂/T₁"),
    )
    for index, (key, name) in enumerate(
        (
            ("downstream_mach", "M₂"),
            ("static_pressure_ratio", "p₂/p₁"),
            ("static_density_ratio", "ρ₂/ρ₁"),
            ("static_temperature_ratio", "T₂/T₁"),
        )
    ):
        row, col = divmod(index, 2)
        states.add_trace(
            go.Scatter(x=x, y=_numeric(rows, key), name=name),
            row=row + 1,
            col=col + 1,
        )
    states.update_xaxes(title_text="上流 Mach M₁")

    pressure = go.Figure()
    pressure.add_trace(
        go.Scatter(x=x, y=_numeric(rows, "total_pressure_ratio"), name="p₀₂/p₀₁")
    )
    pressure.add_trace(
        go.Scatter(x=x, y=_numeric(rows, "pitot_pressure_ratio"), name="p₀₂/p₁")
    )
    pressure.update_xaxes(title_text="上流 Mach M₁")
    pressure.update_yaxes(title_text="圧力比")
    return {
        "状態量": _style(states, "垂直衝撃波の状態量"),
        "全圧・ピトー": _style(pressure, "全圧損失と超音速ピトー圧力比"),
    }


def expansion_figures(
    rows: tuple[Row, ...],
    preferences: UnitPreferences,
    *,
    sweep_field: str,
) -> dict[str, go.Figure]:
    """Create Prandtl-Meyer expansion trend figures."""
    if sweep_field == "turn_angle":
        x = _converted(rows, "turn_angle", "angle", preferences)
        x_title = f"膨張角 θ [{preferences.angle}]"
    else:
        x = _numeric(rows, "upstream_mach")
        x_title = "上流 Mach M₁"
    mach = go.Figure()
    mach.add_trace(go.Scatter(x=x, y=_numeric(rows, "upstream_mach"), name="M₁"))
    mach.add_trace(go.Scatter(x=x, y=_numeric(rows, "downstream_mach"), name="M₂"))
    mach.update_xaxes(title_text=x_title)
    mach.update_yaxes(title_text="Mach数")

    angles = go.Figure()
    for key, name in (
        ("upstream_prandtl_meyer_angle", "ν₁"),
        ("downstream_prandtl_meyer_angle", "ν₂"),
        ("maximum_turn_angle", "最大膨張角"),
    ):
        angles.add_trace(
            go.Scatter(
                x=x,
                y=_converted(rows, key, "angle", preferences),
                name=name,
            )
        )
    angles.update_xaxes(title_text=x_title)
    angles.update_yaxes(title_text=f"角度 [{preferences.angle}]")

    ratios = go.Figure()
    for key, name in (
        ("static_temperature_ratio", "T₂/T₁"),
        ("static_pressure_ratio", "p₂/p₁"),
        ("static_density_ratio", "ρ₂/ρ₁"),
    ):
        ratios.add_trace(go.Scatter(x=x, y=_numeric(rows, key), name=name))
    ratios.update_xaxes(title_text=x_title)
    ratios.update_yaxes(title_text="下流/上流")
    return {
        "Mach数": _style(mach, "Prandtl-Meyer膨張のMach数"),
        "角度": _style(angles, "Prandtl-Meyer角"),
        "状態量比": _style(ratios, "膨張による静的状態量変化"),
    }


def boundary_layer_figures(
    rows: tuple[Row, ...],
    preferences: UnitPreferences,
    *,
    transition_distance: float | None = None,
) -> dict[str, go.Figure]:
    """Create boundary-layer thickness, friction, load, and thermal figures."""
    x = _converted(rows, "distance", "length", preferences)
    x_title = f"x [{preferences.length}]"
    thickness = go.Figure()
    for key, name in (
        ("boundary_layer_thickness", "δ₉₉"),
        ("displacement_thickness", "δ*"),
        ("momentum_thickness", "θ"),
    ):
        thickness.add_trace(
            go.Scatter(
                x=x,
                y=_converted(rows, key, "length", preferences),
                name=name,
            )
        )
    thickness.update_xaxes(title_text=x_title)
    thickness.update_yaxes(title_text=f"厚さ [{preferences.length}]")

    friction = go.Figure()
    friction.add_trace(
        go.Scatter(
            x=x,
            y=_numeric(rows, "local_skin_friction_coefficient"),
            name="局所 C_f",
        )
    )
    friction.add_trace(
        go.Scatter(
            x=x,
            y=_numeric(rows, "average_skin_friction_coefficient"),
            name="平均 C̄_f",
        )
    )
    friction.update_xaxes(title_text=x_title)
    friction.update_yaxes(title_text="摩擦係数")

    loads = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=("壁面せん断応力 τ_w", "単位幅抗力 D′"),
    )
    loads.add_trace(
        go.Scatter(
            x=x,
            y=_converted(rows, "wall_shear_stress", "pressure", preferences),
            name="τ_w",
        ),
        row=1,
        col=1,
    )
    loads.add_trace(
        go.Scatter(x=x, y=_numeric(rows, "drag_per_unit_width"), name="D′"),
        row=1,
        col=2,
    )
    loads.update_xaxes(title_text=x_title)
    loads.update_yaxes(title_text=preferences.pressure, row=1, col=1)
    loads.update_yaxes(title_text="N/m", row=1, col=2)

    figures = {
        "厚さ": _style(thickness, "境界層厚さ"),
        "摩擦係数": _style(friction, "表面摩擦"),
        "せん断・抗力": _style(loads, "壁面荷重"),
    }
    if transition_distance is not None:
        transition = float(from_si(transition_distance, "length", preferences.length))
        for figure in figures.values():
            figure.add_vline(
                x=transition,
                line_dash="dash",
                line_color="#d62728",
                annotation_text="指定遷移",
            )
    if np.any(np.isfinite(_numeric(rows, "recovery_temperature"))):
        thermal = go.Figure()
        for key, name in (
            ("recovery_temperature", "回復温度 T_r"),
            ("wall_temperature", "壁温 T_w"),
        ):
            thermal.add_trace(
                go.Scatter(
                    x=x,
                    y=_converted(rows, key, "temperature", preferences),
                    name=name,
                )
            )
        thermal.update_xaxes(title_text=x_title)
        thermal.update_yaxes(title_text=preferences.temperature)
        figures["温度"] = _style(thermal, "境界層の熱状態")
    return figures


def _model_rows(rows: tuple[Row, ...]) -> dict[str, tuple[Row, ...]]:
    models = {str(row.get("model", "")) for row in rows}
    return {
        model: tuple(row for row in rows if row.get("model") == model)
        for model in sorted(models)
    }


def boundary_layer_profile_figures(
    rows: tuple[Row, ...], preferences: UnitPreferences
) -> dict[str, go.Figure]:
    """Create profile, wall-law, property, and local-flow figures."""
    velocity = go.Figure()
    wall_law = go.Figure()
    properties = make_subplots(
        rows=1,
        cols=3,
        subplot_titles=("温度 T", "密度 ρ", "粘性係数 μ"),
    )
    local = make_subplots(rows=1, cols=2, subplot_titles=("局所 Mach", "動圧 q"))
    for model, selected in _model_rows(rows).items():
        outer = _numeric(selected, "outer_coordinate")
        velocity.add_trace(
            go.Scatter(x=_numeric(selected, "velocity_ratio"), y=outer, name=model)
        )
        wall_law.add_trace(
            go.Scatter(
                x=_numeric(selected, "transformed_wall_coordinate"),
                y=_numeric(selected, "transformed_velocity_plus"),
                name=model,
            )
        )
        property_values = (
            _converted(selected, "temperature", "temperature", preferences),
            _converted(selected, "density", "density", preferences),
            _numeric(selected, "dynamic_viscosity"),
        )
        for index, values in enumerate(property_values):
            properties.add_trace(
                go.Scatter(x=outer, y=values, name=model, legendgroup=model),
                row=1,
                col=index + 1,
            )
        local.add_trace(
            go.Scatter(x=outer, y=_numeric(selected, "local_mach_number"), name=model),
            row=1,
            col=1,
        )
        local.add_trace(
            go.Scatter(
                x=outer,
                y=_converted(selected, "dynamic_pressure", "pressure", preferences),
                name=model,
                legendgroup=model,
            ),
            row=1,
            col=2,
        )
    velocity.update_xaxes(title_text="U/U_e")
    velocity.update_yaxes(title_text="y/δ₉₉")
    wall_law.update_xaxes(title_text="変換壁座標", type="log")
    wall_law.update_yaxes(title_text="変換速度 U⁺")
    for column in range(1, 4):
        properties.update_xaxes(title_text="y/δ₉₉", row=1, col=column)
    properties.update_yaxes(title_text=preferences.temperature, row=1, col=1)
    properties.update_yaxes(title_text=preferences.density, row=1, col=2)
    properties.update_yaxes(title_text="Pa·s", row=1, col=3)
    local.update_xaxes(title_text="y/δ₉₉")
    local.update_yaxes(title_text="Mach", row=1, col=1)
    local.update_yaxes(title_text=preferences.pressure, row=1, col=2)
    return {
        "速度分布": _style(velocity, "圧縮性境界層速度分布"),
        "壁法則": _style(wall_law, "変換壁座標と速度"),
        "熱物性": _style(properties, "温度・密度・粘性係数"),
        "局所流れ": _style(local, "Mach数と動圧"),
    }


def protrusion_shape_figure(
    *,
    height: float,
    base_width: float,
    boundary_layer_thickness: float,
    shape: str,
    preferences: UnitPreferences,
    shape_height: np.ndarray | None = None,
    shape_width: np.ndarray | None = None,
) -> go.Figure:
    """Draw the projected protrusion width and boundary-layer edge."""
    if shape == "csv" and shape_height is not None and shape_width is not None:
        selected = shape_height < height
        y = np.append(shape_height[selected], height)
        width = np.append(
            shape_width[selected], np.interp(height, shape_height, shape_width)
        )
    else:
        y = np.linspace(0.0, height, 201)
        if shape == "rectangle":
            width = np.full_like(y, base_width)
        elif shape == "triangle":
            width = base_width * (1.0 - y / height)
        else:
            width = base_width * np.sqrt(np.maximum(1.0 - (y / height) ** 2, 0.0))
    display_y = np.asarray(from_si(y, "length", preferences.length))
    display_width = np.asarray(from_si(width, "length", preferences.length))
    figure = go.Figure()
    figure.add_trace(
        go.Scatter(
            x=display_width,
            y=display_y,
            fill="tozerox",
            name="投影幅 b(y)",
        )
    )
    delta = float(from_si(boundary_layer_thickness, "length", preferences.length))
    figure.add_hline(
        y=delta,
        line_dash="dash",
        line_color="#d62728",
        annotation_text="δ₉₉",
    )
    figure.update_xaxes(title_text=f"幅 [{preferences.length}]")
    figure.update_yaxes(title_text=f"壁面高さ [{preferences.length}]")
    return _style(figure, "突起の投影形状")


def protrusion_figures(
    rows: tuple[Row, ...],
    preferences: UnitPreferences,
    *,
    sweep_field: str,
) -> dict[str, go.Figure]:
    """Create protrusion drag and shielding trends."""
    if sweep_field in {"height", "base_width", "boundary_layer_thickness"}:
        x = _converted(rows, sweep_field, "length", preferences)
        x_title = f"{sweep_field} [{preferences.length}]"
    else:
        x = _numeric(rows, sweep_field)
        x_title = "Mach M" if sweep_field == "mach" else "抗力係数 C_D"
    loads = make_subplots(rows=1, cols=2, subplot_titles=("直接抗力", "実効動圧"))
    loads.add_trace(
        go.Scatter(x=x, y=_numeric(rows, "direct_drag"), name="D"), row=1, col=1
    )
    loads.add_trace(
        go.Scatter(
            x=x,
            y=_converted(rows, "effective_dynamic_pressure", "pressure", preferences),
            name="q_eff",
        ),
        row=1,
        col=2,
    )
    loads.update_xaxes(title_text=x_title)
    loads.update_yaxes(title_text="N", row=1, col=1)
    loads.update_yaxes(title_text=preferences.pressure, row=1, col=2)
    shielding = go.Figure()
    shielding.add_trace(
        go.Scatter(x=x, y=_numeric(rows, "shielding_factor"), name="遮蔽係数")
    )
    shielding.add_trace(
        go.Scatter(
            x=x,
            y=_numeric(rows, "height_to_boundary_layer_thickness"),
            name="h/δ",
        )
    )
    shielding.update_xaxes(title_text=x_title)
    return {
        "抗力・動圧": _style(loads, "突起抗力"),
        "遮蔽": _style(shielding, "境界層による遮蔽"),
    }


def thermochemistry_figures(
    rows: tuple[Row, ...], preferences: UnitPreferences
) -> dict[str, go.Figure]:
    """Create frozen-air heat-capacity, acoustic, energy, and entropy figures."""
    heat = go.Figure()
    acoustic = make_subplots(rows=1, cols=2, subplot_titles=("比熱比 γ", "音速 a"))
    energy = make_subplots(
        rows=1, cols=2, subplot_titles=("エンタルピー", "内部エネルギー")
    )
    entropy = go.Figure()
    for model, selected in _model_rows(rows).items():
        temperature = _converted(selected, "temperature", "temperature", preferences)
        heat.add_trace(
            go.Scatter(x=temperature, y=_numeric(selected, "cp"), name=f"{model} c_p")
        )
        heat.add_trace(
            go.Scatter(x=temperature, y=_numeric(selected, "cv"), name=f"{model} c_v")
        )
        acoustic.add_trace(
            go.Scatter(
                x=temperature, y=_numeric(selected, "heat_capacity_ratio"), name=model
            ),
            row=1,
            col=1,
        )
        acoustic.add_trace(
            go.Scatter(
                x=temperature,
                y=_converted(selected, "speed_of_sound", "speed", preferences),
                name=model,
                legendgroup=model,
            ),
            row=1,
            col=2,
        )
        for column, standard, sensible in (
            (1, "standard_enthalpy", "sensible_enthalpy"),
            (2, "standard_internal_energy", "sensible_internal_energy"),
        ):
            energy.add_trace(
                go.Scatter(
                    x=temperature,
                    y=_numeric(selected, standard),
                    name=f"{model} standard",
                ),
                row=1,
                col=column,
            )
            energy.add_trace(
                go.Scatter(
                    x=temperature,
                    y=_numeric(selected, sensible),
                    name=f"{model} sensible",
                ),
                row=1,
                col=column,
            )
        entropy.add_trace(
            go.Scatter(x=temperature, y=_numeric(selected, "entropy"), name=model)
        )
    x_title = f"温度 [{preferences.temperature}]"
    heat.update_xaxes(title_text=x_title)
    heat.update_yaxes(title_text="J/(kg·K)")
    acoustic.update_xaxes(title_text=x_title)
    acoustic.update_yaxes(title_text="γ", row=1, col=1)
    acoustic.update_yaxes(title_text=preferences.speed, row=1, col=2)
    energy.update_xaxes(title_text=x_title)
    energy.update_yaxes(title_text="J/kg")
    entropy.update_xaxes(title_text=x_title)
    entropy.update_yaxes(title_text="J/(kg·K)")
    figures = {
        "比熱": _style(heat, "温度依存比熱"),
        "比熱比・音速": _style(acoustic, "比熱比と音速"),
        "エネルギー": _style(energy, "標準値と顕熱値"),
        "エントロピー": _style(entropy, "理想混合気体エントロピー"),
    }
    for limit, label in ((200.0, "適用下限"), (6000.0, "適用上限")):
        display_limit = float(from_si(limit, "temperature", preferences.temperature))
        for figure in figures.values():
            figure.add_vline(
                x=display_limit,
                line_dash="dot",
                line_color="#777777",
                annotation_text=label,
            )
    return figures


def viscosity_figures(
    rows: tuple[Row, ...],
    preferences: UnitPreferences,
    *,
    log_temperature: bool,
) -> dict[str, go.Figure]:
    """Create dynamic-viscosity and Sutherland-relative comparison figures."""
    absolute = go.Figure()
    relative = go.Figure()
    model_names = tuple(
        dict.fromkeys(
            str(row["model"]) for row in rows if isinstance(row.get("model"), str)
        )
    )
    for model in model_names:
        selected = tuple(row for row in rows if row.get("model") == model)
        temperature = (
            _numeric(selected, "temperature")
            if log_temperature
            else _converted(selected, "temperature", "temperature", preferences)
        )
        absolute.add_trace(
            go.Scatter(
                x=temperature,
                y=_numeric(selected, "dynamic_viscosity"),
                name=model,
                connectgaps=False,
            )
        )
        relative.add_trace(
            go.Scatter(
                x=temperature,
                y=_numeric(selected, "relative_difference"),
                name=model,
                connectgaps=False,
            )
        )

    x_title = "温度 [K]" if log_temperature else f"温度 [{preferences.temperature}]"
    axis_type = "log" if log_temperature else "linear"
    for figure in (absolute, relative):
        figure.update_xaxes(title_text=x_title, type=axis_type)
    absolute.update_yaxes(title_text="Pa·s")
    relative.update_yaxes(title_text="%")
    relative.add_hline(line_dash="dot", line_color="#777777", y=0.0)
    return {
        "粘性係数": _style(absolute, "乾燥空気の動的粘性係数"),
        "Sutherland基準相対差": _style(relative, "Sutherland基準の相対差"),
    }
