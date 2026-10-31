"""Aerodynamic heating page."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, cast

import numpy as np
import streamlit as st

from aerophysics.atmosphere import standard_atmosphere
from aerophysics.gui.components import (
    calculation_button,
    clear_widget_state,
    finite_number,
    pop_pending_configuration,
    render_configuration_import,
    render_reset_button,
    render_result_bundle,
)
from aerophysics.gui.config import make_configuration
from aerophysics.gui.documentation import render_documentation_link
from aerophysics.gui.figures import heating_trends
from aerophysics.gui.flow_pages import _defaults, _display, _number, _payload, _si
from aerophysics.gui.heating_adapter import heating_condition, heating_sweep
from aerophysics.gui.units import UnitPreferences


def _on_geometry_change() -> None:
    """Select a model valid for the new geometry and drop stale results."""
    state = st.session_state
    sphere = state.get("heating_geometry") == "sphere"
    state["heating_model"] = "fay_riddell" if sphere else "tauber"
    if sphere and state.get("heating_sweep_field") == "sweep_angle":
        state["heating_sweep_field"] = "upstream_mach"
        clear_widget_state(("heating_sweep_start", "heating_sweep_stop"))
    state.pop("heating_payload", None)


def _on_freestream_change() -> None:
    """Leave an altitude sweep when the freestream is entered manually."""
    state = st.session_state
    if (
        state.get("heating_freestream") == "manual"
        and state.get("heating_sweep_field") == "geometric_altitude"
    ):
        state["heating_sweep_field"] = "upstream_mach"
        clear_widget_state(("heating_sweep_start", "heating_sweep_stop"))


def render_heating(preferences: UnitPreferences) -> None:
    """Render sphere and swept-cylinder aerodynamic heating."""
    st.title("空力加熱")
    render_documentation_link("models/aerodynamic_heating.html")
    imported = pop_pending_configuration("aerodynamic_heating")
    inputs, models, sweep = _defaults(imported)
    render_configuration_import("aerodynamic_heating", "heating")
    render_reset_button("heating", "heating_payload")
    mode = st.radio(
        "計算モード",
        ("single", "sweep"),
        index=1 if imported and imported.get("mode") == "sweep" else 0,
        format_func=lambda x: {"single": "単点", "sweep": "スイープ"}[x],
        horizontal=True,
        key="heating_mode",
    )
    geometry = st.selectbox(
        "形状",
        ("sphere", "swept_cylinder"),
        index=1 if models.get("geometry") == "swept_cylinder" else 0,
        format_func=lambda x: {
            "sphere": "球よどみ点",
            "swept_cylinder": "後退円柱よどみ線",
        }[x],
        key="heating_geometry",
        on_change=_on_geometry_change,
    )
    choices = (
        ("fay_riddell", "tauber", "comparison") if geometry == "sphere" else ("tauber",)
    )
    default_model = models.get("model", choices[0])
    model = st.selectbox(
        "モデル",
        choices,
        index=choices.index(default_model) if default_model in choices else 0,
        format_func=lambda x: {
            "fay_riddell": "Fay–Riddell",
            "tauber": "Tauber",
            "comparison": "Fay–Riddell / Tauber 比較",
        }[x],
        key="heating_model",
        on_change=clear_widget_state,
        args=(("heating_payload",),),
    )
    freestream = st.radio(
        "一様流の与え方",
        ("atmosphere", "manual"),
        index=1 if models.get("freestream") == "manual" else 0,
        format_func=lambda x: {
            "atmosphere": "標準大気（幾何高度）",
            "manual": "静温・静圧を直接入力",
        }[x],
        key="heating_freestream",
        on_change=_on_freestream_change,
    )
    field = "upstream_mach"
    if mode == "sweep":
        fields = (
            "upstream_mach",
            *(("geometric_altitude",) if freestream == "atmosphere" else ()),
            *(("sweep_angle",) if geometry == "swept_cylinder" else ()),
        )
        default_field = sweep.get("field", "upstream_mach")
        field = st.radio(
            "スイープ変数",
            fields,
            index=fields.index(default_field) if default_field in fields else 0,
            format_func=lambda x: {
                "upstream_mach": "Mach",
                "geometric_altitude": "幾何高度",
                "sweep_angle": "後退角",
            }[x],
            key="heating_sweep_field",
            on_change=clear_widget_state,
            args=(("heating_sweep_start", "heating_sweep_stop"),),
        )
    mach = _number(inputs, "upstream_mach", 10.0)
    if mode == "single" or field != "upstream_mach":
        mach = finite_number(
            "上流 Mach M∞", mach, key="heating_mach", min_value=1.0000001
        )
    altitude = None
    temperature = pressure = None
    if freestream == "atmosphere":
        if mode == "single" or field != "geometric_altitude":
            altitude = _si(
                finite_number(
                    f"幾何高度 [{preferences.length}]",
                    _display(
                        float(inputs["geometric_altitude"])
                        if inputs.get("geometric_altitude") is not None
                        else 30000.0,
                        "length",
                        preferences.length,
                    ),
                    key="heating_altitude",
                    min_value=0.0,
                ),
                "length",
                preferences.length,
            )
    else:
        temperature = _si(
            finite_number(
                f"静温 T∞ [{preferences.temperature}]",
                _display(
                    float(inputs.get("temperature") or 226.5),
                    "temperature",
                    preferences.temperature,
                ),
                key="heating_temperature",
            ),
            "temperature",
            preferences.temperature,
        )
        pressure = _si(
            finite_number(
                f"静圧 p∞ [{preferences.pressure}]",
                _display(
                    float(inputs.get("pressure") or 1197.0),
                    "pressure",
                    preferences.pressure,
                ),
                key="heating_pressure",
                min_value=1e-12,
            ),
            "pressure",
            preferences.pressure,
        )
    radius_label = "nose radius Rn" if geometry == "sphere" else "cylinder radius r"
    radius = _si(
        finite_number(
            f"{radius_label} [{preferences.length}]",
            _display(_number(inputs, "radius", 0.1), "length", preferences.length),
            key="heating_radius",
            min_value=1e-12,
        ),
        "length",
        preferences.length,
    )
    angle = _number(inputs, "sweep_angle", np.pi / 6)
    if geometry == "swept_cylinder" and (mode == "single" or field != "sweep_angle"):
        angle = _si(
            st.number_input(
                f"後退角 Λ [{preferences.angle}]",
                min_value=0.0,
                max_value=_display(np.pi / 2, "angle", preferences.angle),
                value=_display(angle, "angle", preferences.angle),
                key="heating_sweep_angle",
            ),
            "angle",
            preferences.angle,
        )
    wall = _si(
        finite_number(
            f"壁温 T_w [{preferences.temperature}]",
            _display(
                _number(inputs, "wall_temperature", 300.0),
                "temperature",
                preferences.temperature,
            ),
            key="heating_wall_temperature",
        ),
        "temperature",
        preferences.temperature,
    )
    prandtl = _number(inputs, "prandtl_number", 0.72)
    if model in {"fay_riddell", "comparison"}:
        prandtl = finite_number(
            "Prandtl 数 Pr", prandtl, key="heating_prandtl", min_value=1e-12
        )
    start = stop = 0.0
    points = 101
    if mode == "sweep":
        kind, unit, name, default_start, default_stop, minimum = {
            "upstream_mach": (None, "", "Mach", 5.0, 20.0, 1.0000001),
            "geometric_altitude": (
                "length",
                preferences.length,
                "高度",
                20000.0,
                60000.0,
                0.0,
            ),
            "sweep_angle": ("angle", preferences.angle, "後退角", 0.0, np.pi / 2, 0.0),
        }[field]
        start_default = _number(sweep, "start", default_start)
        stop_default = _number(sweep, "stop", default_stop)
        if kind is not None:
            start_default = _display(start_default, kind, unit)
            stop_default = _display(stop_default, kind, unit)
            name = f"{name} [{unit}]"
        maximum = _display(np.pi / 2, "angle", unit) if field == "sweep_angle" else None
        a, b, c = st.columns(3)
        with a:
            start = finite_number(
                f"開始 {name}",
                start_default,
                key="heating_sweep_start",
                min_value=minimum,
                max_value=maximum,
            )
        with b:
            stop = finite_number(
                f"終了 {name}",
                stop_default,
                key="heating_sweep_stop",
                min_value=minimum,
                max_value=maximum,
            )
        with c:
            points = int(
                st.number_input(
                    "点数",
                    2,
                    501,
                    int(sweep.get("points", 101)),
                    key="heating_sweep_points",
                )
            )
        if kind is not None:
            start, stop = _si(start, kind, unit), _si(stop, kind, unit)
    if calculation_button("heating_form"):
        st.session_state.pop("heating_payload", None)
        try:
            if altitude is not None:
                atmosphere = standard_atmosphere(altitude)
                temperature, pressure = (
                    float(atmosphere.temperature),
                    float(atmosphere.pressure),
                )
            common: dict[str, Any] = dict(
                radius=radius,
                wall_temperature=wall,
                geometry=geometry,
                model=model,
                sweep_angle=angle,
                prandtl_number=prandtl,
            )
            if mode == "single":
                assert temperature is not None and pressure is not None
                result = heating_condition(
                    upstream_mach=mach,
                    temperature=temperature,
                    pressure=pressure,
                    geometric_altitude=altitude,
                    **common,
                )
                sweep_config = None
            else:
                result = heating_sweep(
                    field=field,
                    start=start,
                    stop=stop,
                    points=points,
                    upstream_mach=mach,
                    temperature=temperature,
                    pressure=pressure,
                    geometric_altitude=altitude,
                    **common,
                )
                sweep_config = {
                    "field": field,
                    "start": start,
                    "stop": stop,
                    "points": points,
                }
            configuration = make_configuration(
                calculator="aerodynamic_heating",
                mode=mode,
                inputs_si={
                    "upstream_mach": mach,
                    "radius": radius,
                    "wall_temperature": wall,
                    "prandtl_number": prandtl,
                    "sweep_angle": angle,
                    "geometric_altitude": altitude,
                    "temperature": temperature if freestream == "manual" else None,
                    "pressure": pressure if freestream == "manual" else None,
                },
                models={"geometry": geometry, "model": model, "freestream": freestream},
                units=preferences,
                sweep_si=sweep_config,
            )
        except ValueError as error:
            st.error(str(error), icon="🚫")
        else:
            st.session_state["heating_payload"] = (result, configuration)
    with st.expander("モデルの前提・適用範囲"):
        st.write(
            "層流、完全気体の空気、超音速一様流を仮定します。"
            "熱流束は壁に入る向きを正とし、高温壁では負になります。"
            "Tauber は 14 km/s、壁温 4000 K、h∞/(V²/2) > 0.1 で警告し、"
            "Fay–Riddell は T₀ > 2000 K で警告します。"
            "円柱は無限円柱のよどみ線で、アタッチメントライン遷移や"
            "翼前縁の合成（TP-2914 式49）は含みません。"
        )
    payload = _payload("heating_payload")
    if payload is None:
        return
    result, configuration = payload
    config_models = configuration["models"]
    assert isinstance(config_models, dict)
    selected_model = config_models["model"]
    selected_geometry = config_models["geometry"]

    def metrics(row: Mapping[str, object]) -> None:
        names = (
            ["Fay–Riddell q_w", "Tauber q_w", "FR/Tauber 比"]
            if selected_model == "comparison"
            else ["熱流束 q_w"]
        )
        names += [
            "よどみ点温度 T₀"
            if selected_geometry == "sphere"
            else "回復エンタルピー h_aw",
            "速度 V∞",
        ]
        columns = st.columns(len(names))
        for column, name in zip(columns, names, strict=True):
            heading = next(
                key for key in row if key == name or key.startswith(name + " [")
            )
            value = row[heading]
            with column:
                st.metric(
                    heading,
                    f"{float(value):.5g}" if isinstance(value, (float, int)) else "—",
                )

    figures = {}
    if configuration["mode"] == "sweep":
        saved_sweep = configuration["sweep_si"]
        assert isinstance(saved_sweep, dict)
        figures = heating_trends(
            result.rows, preferences, field=str(saved_sweep["field"])
        )
    elif selected_geometry == "swept_cylinder":
        row = result.rows[0]
        angle_rows = heating_sweep(
            field="sweep_angle",
            start=0.0,
            stop=np.pi / 2,
            points=91,
            upstream_mach=cast(float, row["upstream_mach"]),
            temperature=cast(float, row["temperature"]),
            pressure=cast(float, row["pressure"]),
            geometric_altitude=None,
            radius=cast(float, row["radius"]),
            wall_temperature=cast(float, row["wall_temperature"]),
            geometry="swept_cylinder",
            model="tauber",
        ).rows
        figures = heating_trends(angle_rows, preferences, field="sweep_angle")
    render_result_bundle(
        calculator="aerodynamic_heating",
        result=result,
        configuration=configuration,
        preferences=preferences,
        figures=figures,
        filename_prefix="aerophysics-heating",
        metrics=metrics,
    )
