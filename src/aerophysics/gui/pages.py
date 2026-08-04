"""Streamlit pages for the GUI prototype."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np
import streamlit as st

from aerophysics.boundary_layer import (
    BoundaryLayerRegime,
    CompressibilityCorrection,
    TurbulentCorrelation,
)
from aerophysics.gui.adapters import (
    CalculationResult,
    FlightCase,
    conical_shock_condition,
    conical_shock_sweep,
    flat_plate,
    flat_plate_sweep,
    flight_condition,
    flight_sweep,
    oblique_shock_condition,
    oblique_shock_sweep,
)
from aerophysics.gui.advanced_adapters import BoundaryLayerCase
from aerophysics.gui.components import (
    finite_number,
    pop_pending_configuration,
    render_configuration_import,
    render_reset_button,
    render_result_bundle,
)
from aerophysics.gui.config import make_configuration
from aerophysics.gui.figures import (
    boundary_layer_figures,
    conical_shock_geometry,
    conical_shock_trends,
    flight_figures,
    shock_geometry,
    shock_trends,
)
from aerophysics.gui.units import UnitPreferences, from_si, to_si
from aerophysics.shocks import (
    ShockBranch,
    maximum_attached_cone_angle,
    maximum_attached_deflection,
)


def _configuration_defaults(
    configuration: dict[str, object] | None,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    if configuration is None:
        return {}, {}, {}
    inputs = configuration.get("inputs_si", {})
    models = configuration.get("models", {})
    sweep = configuration.get("sweep_si", {})
    return (
        dict(inputs) if isinstance(inputs, dict) else {},
        dict(models) if isinstance(models, dict) else {},
        dict(sweep) if isinstance(sweep, dict) else {},
    )


def _display(value: float, kind: str, unit: str) -> float:
    return float(from_si(value, kind, unit))  # type: ignore[arg-type]


def _si(value: float, kind: str, unit: str) -> float:
    return float(to_si(value, kind, unit))  # type: ignore[arg-type]


def _metric(
    label: str, row: Mapping[str, object], heading: str, fmt: str = ".5g"
) -> None:
    value = row.get(heading)
    if isinstance(value, (int, float)):
        st.metric(label, format(float(value), fmt))
    else:
        st.metric(label, "—")


def _result_payload(key: str) -> tuple[CalculationResult, dict[str, object]] | None:
    value = st.session_state.get(key)
    if (
        isinstance(value, tuple)
        and len(value) == 2
        and isinstance(value[0], CalculationResult)
        and isinstance(value[1], dict)
    ):
        return value
    return None


def render_flight(preferences: UnitPreferences) -> None:
    """Render standard-atmosphere and flight-condition calculations."""
    st.title("大気・飛行条件")
    st.caption("U.S. Standard Atmosphere 1976と統合飛行条件を計算します。")
    imported = pop_pending_configuration("flight")
    inputs, models, sweep = _configuration_defaults(imported)
    render_configuration_import("flight", "flight")
    render_reset_button("flight", "flight_payload")

    default_mode = str(imported.get("mode", "single")) if imported else "single"
    default_basis = str(models.get("motion_basis", "mach"))
    altitude_si = float(inputs.get("geometric_altitude", 10_000.0))
    motion_si = float(inputs.get("motion", 0.8))
    length_value = inputs.get("characteristic_length", 1.0)
    length_si = float(length_value) if isinstance(length_value, (int, float)) else 1.0
    has_length = length_value is not None

    with st.form("flight_form"):
        mode = st.radio(
            "計算モード",
            ("single", "sweep"),
            index=0 if default_mode == "single" else 1,
            format_func=lambda value: "単点" if value == "single" else "1変数スイープ",
            horizontal=True,
            key="flight_mode",
        )
        basis = st.radio(
            "運動条件",
            ("mach", "velocity"),
            index=0 if default_basis == "mach" else 1,
            format_func=lambda value: "Mach数" if value == "mach" else "速度",
            horizontal=True,
            key="flight_basis",
        )
        altitude = finite_number(
            f"幾何高度 h [{preferences.length}]",
            _display(altitude_si, "length", preferences.length),
            key="flight_altitude",
        )
        if basis == "mach":
            motion = finite_number(
                "Mach M", motion_si, key="flight_motion", min_value=0.0
            )
        else:
            default_motion = (
                _display(motion_si, "speed", preferences.speed)
                if default_basis == "velocity"
                else _display(250.0, "speed", preferences.speed)
            )
            motion = finite_number(
                f"速度 V [{preferences.speed}]",
                default_motion,
                key="flight_motion",
                min_value=0.0,
            )
        use_length = st.checkbox(
            "代表長さを指定",
            value=has_length,
            key="flight_use_length",
        )
        length = (
            finite_number(
                f"代表長さ L [{preferences.length}]",
                _display(length_si, "length", preferences.length),
                key="flight_length",
                min_value=1e-12,
            )
            if use_length
            else None
        )
        sweep_field = "altitude"
        sweep_start = sweep_stop = 0.0
        points = 101
        if mode == "sweep":
            sweep_field = st.selectbox(
                "スイープ変数",
                ("altitude", "motion"),
                index=0 if sweep.get("field", "altitude") == "altitude" else 1,
                format_func=lambda value: (
                    "幾何高度" if value == "altitude" else "運動条件"
                ),
                key="flight_sweep_field",
            )
            kind = (
                "length"
                if sweep_field == "altitude"
                else ("speed" if basis == "velocity" else None)
            )
            unit = (
                preferences.length
                if sweep_field == "altitude"
                else (preferences.speed if basis == "velocity" else "–")
            )
            default_start_si = float(
                sweep.get("start", 0.0 if sweep_field == "altitude" else 0.2)
            )
            default_stop_si = float(
                sweep.get("stop", 20_000.0 if sweep_field == "altitude" else 3.0)
            )
            default_start = (
                _display(default_start_si, kind, unit) if kind else default_start_si
            )
            default_stop = (
                _display(default_stop_si, kind, unit) if kind else default_stop_si
            )
            left, right, count = st.columns(3)
            with left:
                sweep_start = finite_number(
                    f"開始 [{unit}]",
                    default_start,
                    key="flight_sweep_start",
                )
            with right:
                sweep_stop = finite_number(
                    f"終了 [{unit}]",
                    default_stop,
                    key="flight_sweep_stop",
                )
            with count:
                points = int(
                    st.number_input(
                        "点数",
                        2,
                        501,
                        int(sweep.get("points", 101)),
                        1,
                        key="flight_sweep_points",
                    )
                )
        submitted = st.form_submit_button("計算", type="primary")

    if submitted:
        st.session_state.pop("flight_payload", None)
        try:
            altitude_value_si = _si(altitude, "length", preferences.length)
            motion_value_si = (
                motion if basis == "mach" else _si(motion, "speed", preferences.speed)
            )
            length_result_si = (
                _si(length, "length", preferences.length)
                if length is not None
                else None
            )
            sweep_config: dict[str, object] | None = None
            if mode == "single":
                result = flight_condition(
                    geometric_altitude=altitude_value_si,
                    motion=motion_value_si,
                    motion_basis=basis,
                    characteristic_length=length_result_si,
                )
            else:
                sweep_kind = (
                    "length"
                    if sweep_field == "altitude"
                    else ("speed" if basis == "velocity" else None)
                )
                sweep_unit = (
                    preferences.length
                    if sweep_field == "altitude"
                    else preferences.speed
                )
                start_si = (
                    _si(sweep_start, sweep_kind, sweep_unit)
                    if sweep_kind
                    else sweep_start
                )
                stop_si = (
                    _si(sweep_stop, sweep_kind, sweep_unit)
                    if sweep_kind
                    else sweep_stop
                )
                result = flight_sweep(
                    fixed_altitude=altitude_value_si,
                    fixed_motion=motion_value_si,
                    motion_basis=basis,
                    sweep_field=sweep_field,
                    start=start_si,
                    stop=stop_si,
                    points=points,
                    characteristic_length=length_result_si,
                )
                sweep_config = {
                    "field": sweep_field,
                    "start": start_si,
                    "stop": stop_si,
                    "points": points,
                }
            configuration = make_configuration(
                calculator="flight",
                mode=mode,
                inputs_si={
                    "geometric_altitude": altitude_value_si,
                    "motion": motion_value_si,
                    "characteristic_length": length_result_si,
                },
                models={"motion_basis": basis},
                units=preferences,
                sweep_si=sweep_config,
            )
        except ValueError as error:
            st.error(str(error), icon="🚫")
        else:
            st.session_state["flight_payload"] = (result, configuration)

    payload = _result_payload("flight_payload")
    if payload is None:
        with st.expander("モデルの前提・適用範囲"):
            st.write("幾何高度 -5 km〜86 kmの標準大気です。気象予報ではありません。")
        return
    result, configuration = payload
    if configuration["mode"] == "single":
        if st.button("現在の飛行ケースとして保存", key="flight_save_case"):
            st.session_state["current_flight_case"] = FlightCase.from_row(
                result.rows[0]
            )
            st.success("境界層画面へ引き継ぐ飛行ケースを保存しました。")

    def metrics(row: Mapping[str, object]) -> None:
        columns = st.columns(4)
        headings = list(row)
        with columns[0]:
            _metric("Mach M", row, next(x for x in headings if "Mach M" in x))
        with columns[1]:
            _metric("速度 V", row, next(x for x in headings if "速度 V" in x))
        with columns[2]:
            _metric("動圧 q", row, next(x for x in headings if "動圧 q" in x))
        with columns[3]:
            _metric(
                "Reynolds数",
                row,
                next(x for x in headings if x.startswith("Reynolds数 Re")),
            )

    configuration_sweep = configuration.get("sweep_si")
    figure_sweep_field = (
        str(configuration_sweep.get("field", "altitude"))
        if isinstance(configuration_sweep, dict)
        else "altitude"
    )
    configuration_models = configuration.get("models")
    figure_motion_basis = (
        str(configuration_models.get("motion_basis", "mach"))
        if isinstance(configuration_models, dict)
        else "mach"
    )
    render_result_bundle(
        calculator="flight",
        result=result,
        configuration=configuration,
        preferences=preferences,
        figures=flight_figures(
            result.rows,
            preferences,
            sweep_field=figure_sweep_field,
            motion_basis=figure_motion_basis,
        ),
        filename_prefix="aerophysics-flight",
        metrics=metrics,
    )
    with st.expander("モデルの前提・適用範囲"):
        st.write("U.S. Standard Atmosphere 1976、完全気体AIR、SI計算を使用します。")
        st.page_link(
            "https://github.com/pandorobo11/aerophysics",
            label="モデル文書と数式",
        )


def render_shock(preferences: UnitPreferences) -> None:
    """Render attached oblique-shock calculations."""
    st.title("斜め衝撃波")
    st.caption("theta–beta–Mach関係の弱解・強解を明示的に選択します。")
    imported = pop_pending_configuration("oblique_shock")
    inputs, models, sweep = _configuration_defaults(imported)
    render_configuration_import("oblique_shock", "shock")
    render_reset_button("shock", "shock_payload")
    default_mode = str(imported.get("mode", "single")) if imported else "single"
    default_branch = str(models.get("branch", ShockBranch.WEAK.value))

    with st.form("shock_form"):
        mode = st.radio(
            "計算モード",
            ("single", "sweep"),
            index=0 if default_mode == "single" else 1,
            format_func=lambda value: "単点" if value == "single" else "1変数スイープ",
            horizontal=True,
            key="shock_mode",
        )
        mach = finite_number(
            "上流 Mach M₁",
            float(inputs.get("upstream_mach", 2.0)),
            key="shock_mach",
            min_value=1.0000001,
        )
        theta_default = inputs.get("deflection_angle", np.deg2rad(10.0))
        if not isinstance(theta_default, (int, float)):
            raise ValueError("deflection_angle must be numeric")
        theta_si = float(theta_default)
        theta = finite_number(
            f"偏向角 θ [{preferences.angle}]",
            _display(theta_si, "angle", preferences.angle),
            key="shock_theta",
            min_value=0.0,
        )
        branch = st.selectbox(
            "解の分岐",
            tuple(ShockBranch),
            index=0 if default_branch == ShockBranch.WEAK.value else 1,
            format_func=lambda value: "弱解" if value is ShockBranch.WEAK else "強解",
            key="shock_branch",
        )
        assert branch is not None
        sweep_field = "deflection"
        sweep_start = sweep_stop = 0.0
        points = 101
        if mode == "sweep":
            sweep_field = st.selectbox(
                "スイープ変数",
                ("deflection", "mach"),
                index=0 if sweep.get("field", "deflection") == "deflection" else 1,
                format_func=lambda value: (
                    "偏向角 θ" if value == "deflection" else "Mach M₁"
                ),
                key="shock_sweep_field",
            )
            if sweep_field == "deflection":
                default_limit = float(
                    maximum_attached_deflection(mach).deflection_angle
                )
                start_si = float(sweep.get("start", 0.0))
                stop_si = float(sweep.get("stop", default_limit * 1.05))
                start_default = _display(start_si, "angle", preferences.angle)
                stop_default = _display(stop_si, "angle", preferences.angle)
                unit = preferences.angle
            else:
                start_default = float(sweep.get("start", 1.1))
                stop_default = float(sweep.get("stop", 5.0))
                unit = "–"
            left, right, count = st.columns(3)
            with left:
                sweep_start = finite_number(
                    f"開始 [{unit}]", start_default, key="shock_sweep_start"
                )
            with right:
                sweep_stop = finite_number(
                    f"終了 [{unit}]", stop_default, key="shock_sweep_stop"
                )
            with count:
                points = int(
                    st.number_input(
                        "点数",
                        2,
                        501,
                        int(sweep.get("points", 101)),
                        1,
                        key="shock_sweep_points",
                    )
                )
        submitted = st.form_submit_button("計算", type="primary")

    if submitted:
        st.session_state.pop("shock_payload", None)
        try:
            theta_value_si = _si(theta, "angle", preferences.angle)
            sweep_config: dict[str, object] | None = None
            if mode == "single":
                result = oblique_shock_condition(
                    upstream_mach=mach,
                    deflection_angle=theta_value_si,
                    branch=branch,
                )
            else:
                start_si = (
                    _si(sweep_start, "angle", preferences.angle)
                    if sweep_field == "deflection"
                    else sweep_start
                )
                stop_si = (
                    _si(sweep_stop, "angle", preferences.angle)
                    if sweep_field == "deflection"
                    else sweep_stop
                )
                result = oblique_shock_sweep(
                    fixed_mach=mach,
                    fixed_deflection=theta_value_si,
                    branch=branch,
                    sweep_field=sweep_field,
                    start=start_si,
                    stop=stop_si,
                    points=points,
                )
                sweep_config = {
                    "field": sweep_field,
                    "start": start_si,
                    "stop": stop_si,
                    "points": points,
                }
            configuration = make_configuration(
                calculator="oblique_shock",
                mode=mode,
                inputs_si={
                    "upstream_mach": mach,
                    "deflection_angle": theta_value_si,
                },
                models={"branch": branch.value},
                units=preferences,
                sweep_si=sweep_config,
            )
        except ValueError as error:
            st.error(str(error), icon="🚫")
        else:
            st.session_state["shock_payload"] = (result, configuration)

    payload = _result_payload("shock_payload")
    if payload is None:
        with st.expander("モデルの前提・適用範囲"):
            st.write("定常・熱量的完全気体・付着衝撃波を仮定します。")
        return
    result, configuration = payload

    def metrics(row: Mapping[str, object]) -> None:
        headings = list(row)
        columns = st.columns(4)
        wanted = ("衝撃波角 β", "下流 Mach", "p₂/p₁", "p₀₂/p₀₁")
        for column, label in zip(columns, wanted, strict=True):
            with column:
                _metric(label, row, next(x for x in headings if label in x))

    figures = shock_trends(result.rows, preferences)
    if configuration["mode"] == "single":
        figures = {"模式図": shock_geometry(result.rows[0], preferences), **figures}
    render_result_bundle(
        calculator="oblique_shock",
        result=result,
        configuration=configuration,
        preferences=preferences,
        figures=figures,
        filename_prefix="aerophysics-oblique-shock",
        metrics=metrics,
    )
    invalid = sum(row["status"] != "ok" for row in result.rows)
    if invalid:
        st.warning(
            f"{invalid}点は付着解がないため、垂直衝撃波へ置換せず欠損値としました。"
        )
    with st.expander("モデルの前提・適用範囲"):
        st.write(
            "角度はGUI境界でradianへ変換し、NACA Report 1135の関係式を使用します。"
        )


def render_conical_shock(preferences: UnitPreferences) -> None:
    """Render axisymmetric attached conical-shock calculations."""
    st.title("円錐衝撃波")
    st.caption("Taylor–Maccoll理論による軸対称・付着弱解を計算します。")
    imported = pop_pending_configuration("conical_shock")
    inputs, _, sweep = _configuration_defaults(imported)
    render_configuration_import("conical_shock", "cone_shock")
    render_reset_button("cone_shock", "cone_shock_payload")
    default_mode = str(imported.get("mode", "single")) if imported else "single"

    with st.form("cone_shock_form"):
        mode = st.radio(
            "計算モード",
            ("single", "sweep"),
            index=0 if default_mode == "single" else 1,
            format_func=lambda value: "単点" if value == "single" else "1変数スイープ",
            horizontal=True,
            key="cone_shock_mode",
        )
        mach = finite_number(
            "上流 Mach M∞",
            float(inputs.get("upstream_mach", 2.0)),
            key="cone_shock_mach",
            min_value=1.0000001,
        )
        angle_default = inputs.get("cone_half_angle", np.deg2rad(10.0))
        if not isinstance(angle_default, (int, float)):
            raise ValueError("cone_half_angle must be numeric")
        cone_half_angle = finite_number(
            f"円錐半頂角 θc [{preferences.angle}]",
            _display(float(angle_default), "angle", preferences.angle),
            key="cone_shock_angle",
            min_value=0.0,
        )
        sweep_field = "cone_half_angle"
        sweep_start = sweep_stop = 0.0
        points = 31
        if mode == "sweep":
            sweep_field = st.selectbox(
                "スイープ変数",
                ("cone_half_angle", "mach"),
                index=(
                    0
                    if sweep.get("field", "cone_half_angle") == "cone_half_angle"
                    else 1
                ),
                format_func=lambda value: (
                    "円錐半頂角 θc" if value == "cone_half_angle" else "Mach M∞"
                ),
                key="cone_shock_sweep_field",
            )
            if sweep_field == "cone_half_angle":
                default_limit = float(maximum_attached_cone_angle(mach).cone_half_angle)
                start_si = float(sweep.get("start", 0.0))
                stop_si = float(sweep.get("stop", default_limit * 1.05))
                start_default = _display(start_si, "angle", preferences.angle)
                stop_default = _display(stop_si, "angle", preferences.angle)
                unit = preferences.angle
            else:
                start_default = float(sweep.get("start", 1.1))
                stop_default = float(sweep.get("stop", 5.0))
                unit = "–"
            left, right, count = st.columns(3)
            with left:
                sweep_start = finite_number(
                    f"開始 [{unit}]", start_default, key="cone_shock_sweep_start"
                )
            with right:
                sweep_stop = finite_number(
                    f"終了 [{unit}]", stop_default, key="cone_shock_sweep_stop"
                )
            with count:
                points = int(
                    st.number_input(
                        "点数",
                        2,
                        201,
                        int(sweep.get("points", 31)),
                        1,
                        key="cone_shock_sweep_points",
                    )
                )
        submitted = st.form_submit_button("計算", type="primary")

    if submitted:
        st.session_state.pop("cone_shock_payload", None)
        try:
            angle_si = _si(cone_half_angle, "angle", preferences.angle)
            sweep_config: dict[str, object] | None = None
            if mode == "single":
                result = conical_shock_condition(
                    upstream_mach=mach, cone_half_angle=angle_si
                )
            else:
                start_si = (
                    _si(sweep_start, "angle", preferences.angle)
                    if sweep_field == "cone_half_angle"
                    else sweep_start
                )
                stop_si = (
                    _si(sweep_stop, "angle", preferences.angle)
                    if sweep_field == "cone_half_angle"
                    else sweep_stop
                )
                result = conical_shock_sweep(
                    fixed_mach=mach,
                    fixed_cone_half_angle=angle_si,
                    sweep_field=sweep_field,
                    start=start_si,
                    stop=stop_si,
                    points=points,
                )
                sweep_config = {
                    "field": sweep_field,
                    "start": start_si,
                    "stop": stop_si,
                    "points": points,
                }
            configuration = make_configuration(
                calculator="conical_shock",
                mode=mode,
                inputs_si={
                    "upstream_mach": mach,
                    "cone_half_angle": angle_si,
                },
                models={},
                units=preferences,
                sweep_si=sweep_config,
            )
        except ValueError as error:
            st.error(str(error), icon="🚫")
        else:
            st.session_state["cone_shock_payload"] = (result, configuration)

    payload = _result_payload("cone_shock_payload")
    if payload is None:
        with st.expander("モデルの前提・適用範囲"):
            st.write(
                "迎角0°の鋭い円錐、完全気体、軸対称・非粘性の付着弱解を仮定します。"
            )
        return
    result, configuration = payload

    def metrics(row: Mapping[str, object]) -> None:
        headings = list(row)
        columns = st.columns(4)
        wanted = ("衝撃波角 β", "表面 Mach Mₛ", "pₛ/p∞", "p₀₂/p₀∞")
        for column, label in zip(columns, wanted, strict=True):
            with column:
                _metric(label, row, next(x for x in headings if label in x))

    figures = conical_shock_trends(result.rows, preferences)
    if configuration["mode"] == "single":
        figures = {
            "模式図": conical_shock_geometry(result.rows[0], preferences),
            **figures,
        }
    render_result_bundle(
        calculator="conical_shock",
        result=result,
        configuration=configuration,
        preferences=preferences,
        figures=figures,
        filename_prefix="aerophysics-conical-shock",
        metrics=metrics,
    )
    invalid = sum(row["status"] != "ok" for row in result.rows)
    if invalid:
        st.warning(f"{invalid}点は付着弱解がないため欠損値としました。")
    with st.expander("モデルの前提・適用範囲"):
        st.write(
            "角度はGUI境界でradianへ変換し、Taylor–Maccoll方程式を数値積分します。"
        )


def render_boundary_layer(preferences: UnitPreferences) -> None:
    """Render smooth flat-plate boundary-layer calculations."""
    st.title("平板境界層")
    st.caption("滑面・ゼロ圧力勾配・一定外縁条件の片面平板を計算します。")
    imported = pop_pending_configuration("boundary_layer")
    inputs, models, sweep = _configuration_defaults(imported)
    render_configuration_import("boundary_layer", "boundary")
    render_reset_button("boundary", "boundary_payload")
    case = st.session_state.get("current_flight_case")
    has_case = isinstance(case, FlightCase)
    default_mode = str(imported.get("mode", "single")) if imported else "single"

    with st.form("boundary_form"):
        source_choices = ("manual", "flight") if has_case else ("manual",)
        source = st.radio(
            "外縁条件の入力元",
            source_choices,
            format_func=lambda value: (
                "手入力" if value == "manual" else "現在の飛行ケース"
            ),
            horizontal=True,
            key="boundary_source",
        )
        mode = st.radio(
            "計算モード",
            ("single", "sweep"),
            index=0 if default_mode == "single" else 1,
            format_func=lambda value: "単点" if value == "single" else "距離スイープ",
            horizontal=True,
            key="boundary_mode",
        )
        linked = case if source == "flight" and isinstance(case, FlightCase) else None
        velocity_si = (
            linked.velocity if linked else float(inputs.get("edge_velocity", 100.0))
        )
        density_si = (
            linked.density if linked else float(inputs.get("edge_density", 1.225))
        )
        viscosity_si = (
            linked.dynamic_viscosity
            if linked
            else float(inputs.get("edge_dynamic_viscosity", 1.7894e-5))
        )
        mach_si = linked.mach if linked else float(inputs.get("mach", 0.3))
        temperature_si = (
            linked.temperature
            if linked
            else float(inputs.get("edge_temperature", 288.15))
        )
        distance_si = float(inputs.get("distance", 1.0))
        distance = finite_number(
            f"前縁からの距離 x [{preferences.length}]",
            _display(distance_si, "length", preferences.length),
            key="boundary_distance",
            min_value=1e-12,
            disabled=mode == "sweep",
        )
        edge_velocity = finite_number(
            f"外縁速度 V_e [{preferences.speed}]",
            _display(velocity_si, "speed", preferences.speed),
            key="boundary_velocity",
            min_value=1e-12,
            disabled=linked is not None,
        )
        edge_density = finite_number(
            f"外縁密度 ρ_e [{preferences.density}]",
            _display(density_si, "density", preferences.density),
            key="boundary_density",
            min_value=1e-12,
            disabled=linked is not None,
        )
        edge_viscosity = finite_number(
            "外縁粘性係数 μ_e [Pa·s]",
            viscosity_si,
            key="boundary_viscosity",
            min_value=1e-16,
            disabled=linked is not None,
            format="%.8g",
        )
        regime = st.selectbox(
            "境界層状態",
            tuple(BoundaryLayerRegime),
            index=list(BoundaryLayerRegime).index(
                BoundaryLayerRegime(str(models.get("regime", "turbulent")))
            ),
            format_func=lambda value: {
                BoundaryLayerRegime.LAMINAR: "層流",
                BoundaryLayerRegime.TURBULENT: "乱流",
                BoundaryLayerRegime.TRANSITIONAL: "指定遷移",
            }[value],
            key="boundary_regime",
        )
        assert regime is not None
        correlation = st.selectbox(
            "乱流摩擦相関",
            tuple(TurbulentCorrelation),
            index=list(TurbulentCorrelation).index(
                TurbulentCorrelation(
                    str(models.get("turbulent_correlation", "schlichting"))
                )
            ),
            format_func=lambda value: (
                "Schlichting"
                if value is TurbulentCorrelation.SCHLICHTING
                else "1/5乗則"
            ),
            key="boundary_correlation",
            disabled=regime is BoundaryLayerRegime.LAMINAR,
        )
        assert correlation is not None
        transition_reynolds = (
            finite_number(
                "指定遷移 Reynolds数 Re_tr",
                float(inputs.get("transition_reynolds", 5e5)),
                key="boundary_transition_re",
                min_value=1e-12,
            )
            if regime is BoundaryLayerRegime.TRANSITIONAL
            else None
        )
        correction = st.selectbox(
            "圧縮性補正",
            tuple(CompressibilityCorrection),
            index=list(CompressibilityCorrection).index(
                CompressibilityCorrection(
                    str(models.get("compressibility_correction", "none"))
                )
            ),
            format_func=lambda value: {
                CompressibilityCorrection.NONE: "なし",
                CompressibilityCorrection.ECKERT: "Eckert基準温度法",
                CompressibilityCorrection.VAN_DRIEST_II: "Van Driest II",
            }[value],
            key="boundary_correction",
        )
        assert correction is not None
        mach: float | None = None
        edge_temperature: float | None = None
        wall_temperature: float | None = None
        if correction is not CompressibilityCorrection.NONE:
            mach = finite_number(
                "Mach M_e",
                mach_si,
                key="boundary_mach",
                min_value=0.0,
                disabled=linked is not None,
            )
            edge_temperature_display = finite_number(
                f"外縁温度 T_e [{preferences.temperature}]",
                _display(temperature_si, "temperature", preferences.temperature),
                key="boundary_temperature",
                disabled=linked is not None,
            )
            edge_temperature = _si(
                edge_temperature_display, "temperature", preferences.temperature
            )
            adiabatic = st.checkbox(
                "断熱壁（壁温を自動計算）",
                value=inputs.get("wall_temperature") is None,
                key="boundary_adiabatic",
            )
            if not adiabatic:
                wall_default = float(inputs.get("wall_temperature", 300.0))
                wall_display = finite_number(
                    f"壁温 T_w [{preferences.temperature}]",
                    _display(wall_default, "temperature", preferences.temperature),
                    key="boundary_wall_temperature",
                )
                wall_temperature = _si(
                    wall_display, "temperature", preferences.temperature
                )
        start_display = stop_display = 0.0
        points = 201
        logarithmic = True
        if mode == "sweep":
            start_si = float(sweep.get("start", 0.01))
            stop_si = float(sweep.get("stop", 2.0))
            left, right, count = st.columns(3)
            with left:
                start_display = finite_number(
                    f"開始距離 [{preferences.length}]",
                    _display(start_si, "length", preferences.length),
                    key="boundary_sweep_start",
                    min_value=1e-12,
                )
            with right:
                stop_display = finite_number(
                    f"終了距離 [{preferences.length}]",
                    _display(stop_si, "length", preferences.length),
                    key="boundary_sweep_stop",
                    min_value=1e-12,
                )
            with count:
                points = int(
                    st.number_input(
                        "点数",
                        2,
                        501,
                        int(sweep.get("points", 201)),
                        1,
                        key="boundary_sweep_points",
                    )
                )
            logarithmic = st.checkbox(
                "対数間隔",
                value=bool(sweep.get("logarithmic", True)),
                key="boundary_sweep_log",
            )
        submitted = st.form_submit_button("計算", type="primary")

    if submitted:
        st.session_state.pop("boundary_payload", None)
        try:
            velocity_value_si = _si(edge_velocity, "speed", preferences.speed)
            density_value_si = _si(edge_density, "density", preferences.density)
            sweep_config: dict[str, object] | None = None
            if mode == "single":
                distance_value_si = _si(distance, "length", preferences.length)
                result = flat_plate(
                    distance=distance_value_si,
                    edge_velocity=velocity_value_si,
                    edge_density=density_value_si,
                    edge_dynamic_viscosity=edge_viscosity,
                    regime=regime,
                    turbulent_correlation=correlation,
                    transition_reynolds=transition_reynolds,
                    compressibility_correction=correction,
                    mach=mach,
                    edge_temperature=edge_temperature,
                    wall_temperature=wall_temperature,
                )
            else:
                start_value_si = _si(start_display, "length", preferences.length)
                stop_value_si = _si(stop_display, "length", preferences.length)
                distance_value_si = start_value_si
                result = flat_plate_sweep(
                    start=start_value_si,
                    stop=stop_value_si,
                    points=points,
                    logarithmic=logarithmic,
                    edge_velocity=velocity_value_si,
                    edge_density=density_value_si,
                    edge_dynamic_viscosity=edge_viscosity,
                    regime=regime,
                    turbulent_correlation=correlation,
                    transition_reynolds=transition_reynolds,
                    compressibility_correction=correction,
                    mach=mach,
                    edge_temperature=edge_temperature,
                    wall_temperature=wall_temperature,
                )
                sweep_config = {
                    "field": "distance",
                    "start": start_value_si,
                    "stop": stop_value_si,
                    "points": points,
                    "logarithmic": logarithmic,
                }
            configuration = make_configuration(
                calculator="boundary_layer",
                mode=mode,
                inputs_si={
                    "distance": distance_value_si,
                    "edge_velocity": velocity_value_si,
                    "edge_density": density_value_si,
                    "edge_dynamic_viscosity": edge_viscosity,
                    "transition_reynolds": transition_reynolds,
                    "mach": mach,
                    "edge_temperature": (
                        edge_temperature
                        if edge_temperature is not None
                        else temperature_si
                    ),
                    "wall_temperature": wall_temperature,
                },
                models={
                    "source": source,
                    "regime": regime.value,
                    "turbulent_correlation": correlation.value,
                    "compressibility_correction": correction.value,
                },
                units=preferences,
                sweep_si=sweep_config,
            )
        except ValueError as error:
            st.error(str(error), icon="🚫")
        else:
            st.session_state["boundary_payload"] = (result, configuration)

    payload = _result_payload("boundary_payload")
    if payload is None:
        if not has_case:
            st.info(
                "大気・飛行条件画面で単点ケースを保存すると外縁条件を引き継げます。"
            )
        with st.expander("モデルの前提・適用範囲"):
            st.write("滑面、鋭い前縁、ゼロ圧力勾配、一定外縁条件を仮定します。")
        return
    result, configuration = payload
    input_values = configuration["inputs_si"]
    transition_distance = None
    if isinstance(input_values, dict):
        transition = input_values.get("transition_reynolds")
        velocity = input_values.get("edge_velocity")
        density = input_values.get("edge_density")
        viscosity = input_values.get("edge_dynamic_viscosity")
        if (
            isinstance(transition, (int, float))
            and isinstance(velocity, (int, float))
            and isinstance(density, (int, float))
            and isinstance(viscosity, (int, float))
        ):
            transition_distance = (
                float(transition)
                * float(viscosity)
                / (float(density) * float(velocity))
            )

    def metrics(row: Mapping[str, object]) -> None:
        headings = list(row)
        columns = st.columns(4)
        wanted = ("Reynolds数 Re_x", "境界層厚さ", "局所摩擦係数", "単位幅抗力")
        for column, label in zip(columns, wanted, strict=True):
            with column:
                _metric(label, row, next(x for x in headings if label in x))

    render_result_bundle(
        calculator="boundary_layer",
        result=result,
        configuration=configuration,
        preferences=preferences,
        figures=boundary_layer_figures(
            result.rows,
            preferences,
            transition_distance=transition_distance,
        ),
        filename_prefix="aerophysics-boundary-layer",
        metrics=metrics,
    )
    config_models = configuration.get("models")
    if (
        configuration.get("mode") == "single"
        and isinstance(config_models, dict)
        and config_models.get("regime") == BoundaryLayerRegime.TURBULENT.value
    ):
        if st.button("現在の境界層ケースとして保存", key="boundary_save_case"):
            row = result.rows[0]
            input_values = configuration.get("inputs_si")
            if isinstance(input_values, dict):
                st.session_state["current_boundary_layer_case"] = BoundaryLayerCase(
                    edge_velocity=float(input_values["edge_velocity"]),
                    edge_density=float(input_values["edge_density"]),
                    edge_temperature=float(input_values["edge_temperature"]),
                    boundary_layer_thickness=float(row["boundary_layer_thickness"]),
                    wall_shear_stress=float(row["wall_shear_stress"]),
                )
                st.success("乱流境界層ケースを保存しました。")
    with st.expander("モデルの前提・適用範囲"):
        st.write(
            "乱流相関の公称範囲は 5e5 ≤ Re_x ≤ 1e9です。"
            "表面粗さ、圧力勾配、剥離、自然遷移予測は含みません。"
        )
