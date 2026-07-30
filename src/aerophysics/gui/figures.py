"""Plotly figures for GUI calculation results."""

from __future__ import annotations

import math

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

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
