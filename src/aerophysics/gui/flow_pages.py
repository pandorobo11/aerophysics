"""Streamlit pages for additional compressible-flow calculators."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import streamlit as st

from aerophysics.detached_shock import (
    BilligShockShapeResult,
    DetachedShockGeometry,
)
from aerophysics.gui.adapters import (
    CalculationResult,
    detached_shock_condition,
    detached_shock_shape,
    detached_shock_sweep,
    expansion_condition,
    expansion_sweep,
    isentropic_condition,
    isentropic_sweep,
    normal_shock_condition,
    normal_shock_sweep,
)
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
from aerophysics.gui.figures import (
    detached_shock_geometry,
    detached_shock_trends,
    expansion_figures,
    isentropic_figures,
    normal_shock_figures,
)
from aerophysics.gui.page_support import (
    configuration_defaults,
    display_value,
    numeric_default,
    render_metric,
    session_payload,
    si_value,
)
from aerophysics.gui.tables import detached_shock_shape_csv
from aerophysics.gui.units import UnitPreferences
from aerophysics.isentropic import MachBranch

_ISENTROPIC_LABELS = {
    "mach": "Mach M",
    "temperature_ratio": "全温/静温 T₀/T",
    "pressure_ratio": "全圧/静圧 p₀/p",
    "density_ratio": "全密度/静密度 ρ₀/ρ",
    "area_ratio": "面積比 A/A*",
}


def render_isentropic(preferences: UnitPreferences) -> None:
    """Render isentropic forward/inverse and mass-flow calculations."""
    st.title("等エントロピー流れ")
    st.caption("状態量比、面積-Mach関係、逆計算、質量流束を計算します。")
    imported = pop_pending_configuration("isentropic")
    inputs, models, sweep = configuration_defaults(imported)
    render_configuration_import("isentropic", "isentropic")
    render_reset_button("isentropic", "isentropic_payload")
    default_mode = str(imported.get("mode", "single")) if imported else "single"
    default_basis = str(models.get("input_basis", "mach"))
    default_branch = str(models.get("branch", MachBranch.SUBSONIC.value))
    default_gas_model = str(models.get("gas_model", "AIR"))
    stored_total_temperature = inputs.get("total_temperature", 300.0)
    default_total_temperature = (
        float(stored_total_temperature)
        if isinstance(stored_total_temperature, (int, float))
        else 300.0
    )

    with st.container():
        mode = st.radio(
            "計算モード",
            ("single", "sweep"),
            index=0 if default_mode == "single" else 1,
            format_func=lambda value: "単点" if value == "single" else "入力値スイープ",
            horizontal=True,
            key="isentropic_mode",
        )
        gas_models = (
            "AIR",
            "NASA7",
            "NASA9",
            "HARMONIC_OSCILLATOR",
            "BEATTIE_BRIDGEMAN",
        )
        gas_model = st.selectbox(
            "気体モデル",
            gas_models,
            index=(
                gas_models.index(default_gas_model)
                if default_gas_model in gas_models
                else 0
            ),
            format_func=lambda value: {
                "AIR": "定比熱 AIR",
                "NASA7": "熱的完全気体 NASA7",
                "NASA9": "熱的完全気体 NASA9",
                "HARMONIC_OSCILLATOR": "Harmonic-oscillator air",
                "BEATTIE_BRIDGEMAN": "Beattie–Bridgeman air",
            }[value],
            key="isentropic_gas_model",
        )
        assert gas_model is not None
        bases = tuple(_ISENTROPIC_LABELS)
        basis = st.selectbox(
            "入力基準",
            bases,
            index=bases.index(default_basis) if default_basis in bases else 0,
            format_func=_ISENTROPIC_LABELS.__getitem__,
            key="isentropic_basis",
            on_change=clear_widget_state,
            args=(
                (
                    "isentropic_input",
                    "isentropic_sweep_start",
                    "isentropic_sweep_stop",
                ),
            ),
        )
        assert basis is not None
        minimum = 0.0 if basis == "mach" else 1.0
        input_value = finite_number(
            _ISENTROPIC_LABELS[basis],
            float(inputs.get("input_value", 2.0)),
            key="isentropic_input",
            min_value=minimum,
        )
        branch = MachBranch.SUBSONIC
        if basis == "area_ratio":
            branch = st.selectbox(
                "面積-Mach逆計算の分岐",
                tuple(MachBranch),
                index=0 if default_branch == MachBranch.SUBSONIC.value else 1,
                format_func=lambda value: (
                    "亜音速枝" if value is MachBranch.SUBSONIC else "超音速枝"
                ),
                key="isentropic_branch",
            )
            assert branch is not None
        requires_pressure = gas_model == "BEATTIE_BRIDGEMAN"
        with_mass_flux_selection = st.checkbox(
            "全圧を指定して質量流束を計算",
            value=requires_pressure or bool(models.get("with_mass_flux", False)),
            disabled=requires_pressure,
            key="isentropic_with_flux",
        )
        with_mass_flux = requires_pressure or with_mass_flux_selection
        total_temperature_display = finite_number(
            f"全温 T₀ [{preferences.temperature}]",
            display_value(
                default_total_temperature,
                "temperature",
                preferences.temperature,
            ),
            key="isentropic_total_temperature",
        )
        total_temperature = si_value(
            total_temperature_display, "temperature", preferences.temperature
        )
        allow_extrapolation = st.checkbox(
            "モデルの文書化された適用範囲外への外挿を許可",
            value=bool(models.get("allow_extrapolation", True)),
            disabled=gas_model == "AIR",
            key="isentropic_allow_extrapolation",
        )
        total_pressure = None
        if with_mass_flux:
            pressure_display = finite_number(
                f"全圧 p₀ [{preferences.pressure}]",
                display_value(
                    float(inputs.get("total_pressure", 101_325.0)),
                    "pressure",
                    preferences.pressure,
                ),
                key="isentropic_total_pressure",
                min_value=1e-12,
            )
            total_pressure = si_value(
                pressure_display, "pressure", preferences.pressure
            )
        start = stop = 0.0
        points = 101
        if mode == "sweep":
            start_default = float(sweep.get("start", 0.1 if basis == "mach" else 1.0))
            stop_default = float(sweep.get("stop", 5.0))
            left, right, count = st.columns(3)
            with left:
                start = finite_number(
                    "開始",
                    start_default,
                    key="isentropic_sweep_start",
                    min_value=minimum,
                )
            with right:
                stop = finite_number(
                    "終了", stop_default, key="isentropic_sweep_stop", min_value=minimum
                )
            with count:
                points = int(
                    st.number_input(
                        "点数",
                        2,
                        501,
                        int(sweep.get("points", 101)),
                        1,
                        key="isentropic_sweep_points",
                    )
                )
        submitted = calculation_button("isentropic_form")

    if submitted:
        st.session_state.pop("isentropic_payload", None)
        try:
            sweep_configuration: dict[str, object] | None = None
            if mode == "single":
                result = isentropic_condition(
                    input_value=input_value,
                    input_basis=basis,
                    branch=branch,
                    gas_model=gas_model,
                    total_pressure=total_pressure,
                    total_temperature=total_temperature,
                    allow_extrapolation=allow_extrapolation,
                )
            else:
                result = isentropic_sweep(
                    input_basis=basis,
                    branch=branch,
                    start=start,
                    stop=stop,
                    points=points,
                    gas_model=gas_model,
                    total_pressure=total_pressure,
                    total_temperature=total_temperature,
                    allow_extrapolation=allow_extrapolation,
                )
                sweep_configuration = {
                    "field": "input_value",
                    "start": start,
                    "stop": stop,
                    "points": points,
                }
            configuration = make_configuration(
                calculator="isentropic",
                mode=mode,
                inputs_si={
                    "input_value": input_value,
                    "total_pressure": total_pressure,
                    "total_temperature": total_temperature,
                },
                models={
                    "input_basis": basis,
                    "branch": branch.value,
                    "gas_model": gas_model,
                    "with_mass_flux": with_mass_flux,
                    "allow_extrapolation": allow_extrapolation,
                },
                units=preferences,
                sweep_si=sweep_configuration,
            )
        except ValueError as error:
            st.error(str(error), icon="🚫")
        else:
            st.session_state["isentropic_payload"] = (result, configuration)

    payload = session_payload("isentropic_payload", CalculationResult)
    if payload is None:
        with st.expander("モデルの前提・適用範囲"):
            st.write(
                "定常・断熱・可逆な流れを仮定します。NASA7/NASA9と調和振動子"
                "モデルは凍結組成の熱的完全気体です。Beattie–Bridgemanモデルは"
                "圧力依存の高密度補正を含みます。"
            )
        return
    result, configuration = payload

    def metrics(row: Mapping[str, object]) -> None:
        columns = st.columns(4)
        wanted = (
            ("Mach M", "Mach M"),
            ("T₀/T", "T₀/T"),
            ("p₀/p", "p₀/p"),
            ("A/A*", "A/A*"),
        )
        for column, (label, contains) in zip(columns, wanted, strict=True):
            with column:
                render_metric(row, label, contains)

    render_result_bundle(
        calculator="isentropic",
        result=result,
        configuration=configuration,
        preferences=preferences,
        figures=isentropic_figures(result.rows, input_label=_ISENTROPIC_LABELS[basis]),
        filename_prefix="aerophysics-isentropic",
        metrics=metrics,
    )
    with st.expander("モデルの前提・適用範囲"):
        st.write(
            "状態量比は全量/静量です。A/A* > 1の逆計算では亜音速枝または"
            "超音速枝を必ず選択します。NASA7/NASA9と調和振動子モデルは凍結"
            "組成です。Beattie–Bridgemanを含め、解離・反応・電離・相変化は"
            "扱いません。"
        )


def render_normal_shock(preferences: UnitPreferences) -> None:
    """Render normal-shock and supersonic-pitot calculations."""
    st.title("垂直衝撃波")
    st.caption("衝撃波前後の状態量比、全圧損失、超音速ピトー圧力比を計算します。")
    imported = pop_pending_configuration("normal_shock")
    inputs, _, sweep = configuration_defaults(imported)
    render_configuration_import("normal_shock", "normal")
    render_reset_button("normal", "normal_payload")
    default_mode = str(imported.get("mode", "single")) if imported else "single"

    with st.container():
        mode = st.radio(
            "計算モード",
            ("single", "sweep"),
            index=0 if default_mode == "single" else 1,
            format_func=lambda value: "単点" if value == "single" else "Machスイープ",
            horizontal=True,
            key="normal_mode",
        )
        mach = finite_number(
            "上流 Mach M₁",
            float(inputs.get("upstream_mach", 2.0)),
            key="normal_mach",
            min_value=1.0,
        )
        start = stop = 0.0
        points = 101
        if mode == "sweep":
            left, right, count = st.columns(3)
            with left:
                start = finite_number(
                    "開始 Mach",
                    float(sweep.get("start", 1.0)),
                    key="normal_sweep_start",
                    min_value=1.0,
                )
            with right:
                stop = finite_number(
                    "終了 Mach",
                    float(sweep.get("stop", 8.0)),
                    key="normal_sweep_stop",
                    min_value=1.0,
                )
            with count:
                points = int(
                    st.number_input(
                        "点数",
                        2,
                        501,
                        int(sweep.get("points", 101)),
                        1,
                        key="normal_sweep_points",
                    )
                )
        submitted = calculation_button("normal_form")

    if submitted:
        st.session_state.pop("normal_payload", None)
        try:
            sweep_configuration: dict[str, object] | None = None
            if mode == "single":
                result = normal_shock_condition(upstream_mach=mach)
            else:
                result = normal_shock_sweep(start=start, stop=stop, points=points)
                sweep_configuration = {
                    "field": "upstream_mach",
                    "start": start,
                    "stop": stop,
                    "points": points,
                }
            configuration = make_configuration(
                calculator="normal_shock",
                mode=mode,
                inputs_si={"upstream_mach": mach},
                models={},
                units=preferences,
                sweep_si=sweep_configuration,
            )
        except ValueError as error:
            st.error(str(error), icon="🚫")
        else:
            st.session_state["normal_payload"] = (result, configuration)

    payload = session_payload("normal_payload", CalculationResult)
    if payload is None:
        with st.expander("モデルの前提・適用範囲"):
            st.write("定常・断熱な垂直衝撃波と熱量的完全気体AIRを仮定します。")
        return
    result, configuration = payload

    def metrics(row: Mapping[str, object]) -> None:
        columns = st.columns(4)
        wanted = (
            ("M₂", "下流 Mach"),
            ("p₂/p₁", "p₂/p₁"),
            ("T₂/T₁", "T₂/T₁"),
            ("p₀₂/p₀₁", "p₀₂/p₀₁"),
        )
        for column, (label, contains) in zip(columns, wanted, strict=True):
            with column:
                render_metric(row, label, contains)

    render_result_bundle(
        calculator="normal_shock",
        result=result,
        configuration=configuration,
        preferences=preferences,
        figures=normal_shock_figures(result.rows),
        filename_prefix="aerophysics-normal-shock",
        metrics=metrics,
    )
    with st.expander("モデルの前提・適用範囲"):
        st.write("状態量比は下流/上流、全圧比はp₀₂/p₀₁です。")


def render_expansion(preferences: UnitPreferences) -> None:
    """Render centered Prandtl-Meyer expansion calculations."""
    st.title("Prandtl–Meyer膨張")
    st.caption("超音速流の中心膨張におけるMach数、角度、静的状態量比を計算します。")
    imported = pop_pending_configuration("expansion")
    inputs, _, sweep = configuration_defaults(imported)
    render_configuration_import("expansion", "expansion")
    render_reset_button("expansion", "expansion_payload")
    default_mode = str(imported.get("mode", "single")) if imported else "single"

    with st.container():
        mode = st.radio(
            "計算モード",
            ("single", "sweep"),
            index=0 if default_mode == "single" else 1,
            format_func=lambda value: "単点" if value == "single" else "1変数スイープ",
            horizontal=True,
            key="expansion_mode",
        )
        mach = finite_number(
            "上流 Mach M₁",
            float(inputs.get("upstream_mach", 2.0)),
            key="expansion_mach",
            min_value=1.0,
        )
        turn_si = numeric_default(inputs, "turn_angle", float(np.deg2rad(15.0)))
        turn_display = finite_number(
            f"膨張角 θ [{preferences.angle}]",
            display_value(turn_si, "angle", preferences.angle),
            key="expansion_turn",
            min_value=0.0,
        )
        sweep_field = "turn_angle"
        start = stop = 0.0
        points = 101
        if mode == "sweep":
            sweep_field = st.selectbox(
                "スイープ変数",
                ("turn_angle", "mach"),
                index=0 if sweep.get("field", "turn_angle") == "turn_angle" else 1,
                format_func=lambda value: (
                    "膨張角 θ" if value == "turn_angle" else "上流 Mach M₁"
                ),
                key="expansion_sweep_field",
                on_change=clear_widget_state,
                args=(("expansion_sweep_start", "expansion_sweep_stop"),),
            )
            assert sweep_field is not None
            if sweep_field == "turn_angle":
                start_default = display_value(
                    float(sweep.get("start", 0.0)), "angle", preferences.angle
                )
                stop_default = display_value(
                    numeric_default(sweep, "stop", float(np.deg2rad(80.0))),
                    "angle",
                    preferences.angle,
                )
                minimum = 0.0
                unit = preferences.angle
            else:
                start_default = float(sweep.get("start", 1.0))
                stop_default = float(sweep.get("stop", 6.0))
                minimum = 1.0
                unit = "–"
            left, right, count = st.columns(3)
            with left:
                start = finite_number(
                    f"開始 [{unit}]",
                    start_default,
                    key="expansion_sweep_start",
                    min_value=minimum,
                )
            with right:
                stop = finite_number(
                    f"終了 [{unit}]",
                    stop_default,
                    key="expansion_sweep_stop",
                    min_value=minimum,
                )
            with count:
                points = int(
                    st.number_input(
                        "点数",
                        2,
                        501,
                        int(sweep.get("points", 101)),
                        1,
                        key="expansion_sweep_points",
                    )
                )
        submitted = calculation_button("expansion_form")

    if submitted:
        st.session_state.pop("expansion_payload", None)
        try:
            turn_angle = si_value(turn_display, "angle", preferences.angle)
            sweep_configuration: dict[str, object] | None = None
            if mode == "single":
                result = expansion_condition(upstream_mach=mach, turn_angle=turn_angle)
            else:
                start_si = (
                    si_value(start, "angle", preferences.angle)
                    if sweep_field == "turn_angle"
                    else start
                )
                stop_si = (
                    si_value(stop, "angle", preferences.angle)
                    if sweep_field == "turn_angle"
                    else stop
                )
                result = expansion_sweep(
                    fixed_mach=mach,
                    fixed_turn_angle=turn_angle,
                    sweep_field=sweep_field,
                    start=start_si,
                    stop=stop_si,
                    points=points,
                )
                sweep_configuration = {
                    "field": sweep_field,
                    "start": start_si,
                    "stop": stop_si,
                    "points": points,
                }
            configuration = make_configuration(
                calculator="expansion",
                mode=mode,
                inputs_si={"upstream_mach": mach, "turn_angle": turn_angle},
                models={},
                units=preferences,
                sweep_si=sweep_configuration,
            )
        except ValueError as error:
            st.error(str(error), icon="🚫")
        else:
            st.session_state["expansion_payload"] = (result, configuration)

    payload = session_payload("expansion_payload", CalculationResult)
    if payload is None:
        with st.expander("モデルの前提・適用範囲"):
            st.write("定常・等エントロピーな中心膨張と熱量的完全気体AIRを仮定します。")
        return
    result, configuration = payload

    def metrics(row: Mapping[str, object]) -> None:
        columns = st.columns(4)
        wanted = (
            ("M₂", "下流 Mach"),
            ("ν₂", "ν₂"),
            ("p₂/p₁", "p₂/p₁"),
            ("T₂/T₁", "T₂/T₁"),
        )
        for column, (label, contains) in zip(columns, wanted, strict=True):
            with column:
                render_metric(row, label, contains)

    config_sweep = configuration.get("sweep_si")
    figure_field = (
        str(config_sweep.get("field", "mach"))
        if isinstance(config_sweep, dict)
        else "mach"
    )
    render_result_bundle(
        calculator="expansion",
        result=result,
        configuration=configuration,
        preferences=preferences,
        figures=expansion_figures(result.rows, preferences, sweep_field=figure_field),
        filename_prefix="aerophysics-expansion",
        metrics=metrics,
    )
    invalid = sum(row["status"] != "ok" for row in result.rows)
    if invalid:
        st.warning(f"{invalid}点はPrandtl–Meyer角の極限を超えるため欠損値としました。")
    with st.expander("モデルの前提・適用範囲"):
        st.write("膨張前後で全温・全圧は一定です。角度はGUI境界でradianへ変換します。")


def render_detached_shock(preferences: UnitPreferences) -> None:
    """Render detached-shock standoff and Billig-shape calculations."""
    st.title("離脱衝撃波")
    st.caption("鈍頭物体の離脱距離を推算し、Billigの双曲線衝撃波形状を表示します。")
    imported = pop_pending_configuration("detached_shock")
    inputs, models, sweep = configuration_defaults(imported)
    render_configuration_import("detached_shock", "detached_shock")
    render_reset_button("detached_shock", "detached_shock_payload")
    default_mode = str(imported.get("mode", "single")) if imported else "single"
    default_geometry = str(
        models.get("geometry", DetachedShockGeometry.AXISYMMETRIC_SPHERE.value)
    )
    default_selection = str(
        models.get("model", models.get("selection", "ambrosio_wortman"))
    )

    with st.container():
        mode = st.radio(
            "計算モード",
            ("single", "sweep"),
            index=0 if default_mode == "single" else 1,
            format_func=lambda value: "単点" if value == "single" else "Machスイープ",
            horizontal=True,
            key="detached_shock_mode",
        )
        geometries = tuple(DetachedShockGeometry)
        geometry = st.selectbox(
            "geometry",
            geometries,
            index=(
                next(
                    (
                        index
                        for index, item in enumerate(geometries)
                        if item.value == default_geometry
                    ),
                    0,
                )
            ),
            format_func=lambda value: (
                "axisymmetric sphere / hemispherical nose"
                if value is DetachedShockGeometry.AXISYMMETRIC_SPHERE
                else "2D cylindrical nose"
            ),
            key="detached_shock_geometry",
        )
        assert geometry is not None
        selections = (
            ("ambrosio_wortman", "seiff", "comparison")
            if geometry is DetachedShockGeometry.AXISYMMETRIC_SPHERE
            else ("ambrosio_wortman",)
        )
        selection = st.selectbox(
            "離脱距離モデル",
            selections,
            index=(
                selections.index(default_selection)
                if default_selection in selections
                else 0
            ),
            format_func=lambda value: {
                "ambrosio_wortman": "Ambrosio–Wortman",
                "seiff": "Seiff（定比熱AIRの垂直衝撃波密度比）",
                "comparison": "Ambrosio–Wortman / Seiff 比較",
            }[value],
            key="detached_shock_selection",
        )
        assert selection is not None
        mach = finite_number(
            "上流 Mach M∞",
            numeric_default(inputs, "upstream_mach", 4.0),
            key="detached_shock_mach",
            min_value=1.0000001,
        )
        radius_si = numeric_default(inputs, "nose_radius", 0.1)
        radius_display = finite_number(
            f"nose radius Rn [{preferences.length}]",
            display_value(radius_si, "length", preferences.length),
            key="detached_shock_radius",
            min_value=1.0e-12,
        )
        radius = si_value(radius_display, "length", preferences.length)
        start = stop = 0.0
        points = 101
        if mode == "sweep":
            left, right, count = st.columns(3)
            with left:
                start = finite_number(
                    "開始 Mach",
                    numeric_default(sweep, "start", 1.5),
                    key="detached_shock_sweep_start",
                    min_value=1.0000001,
                )
            with right:
                stop = finite_number(
                    "終了 Mach",
                    numeric_default(sweep, "stop", 10.0),
                    key="detached_shock_sweep_stop",
                    min_value=1.0000001,
                )
            with count:
                points = int(
                    st.number_input(
                        "点数",
                        2,
                        501,
                        int(sweep.get("points", 101)),
                        1,
                        key="detached_shock_sweep_points",
                    )
                )
        submitted = calculation_button("detached_shock_form")

    if submitted:
        st.session_state.pop("detached_shock_payload", None)
        st.session_state.pop("detached_shock_shape_result", None)
        try:
            sweep_configuration: dict[str, object] | None = None
            if mode == "single":
                result = detached_shock_condition(
                    upstream_mach=mach,
                    nose_radius=radius,
                    geometry=geometry,
                    selection=selection,
                )
                shape = detached_shock_shape(
                    upstream_mach=mach,
                    nose_radius=radius,
                    geometry=geometry,
                )
            else:
                result = detached_shock_sweep(
                    start=start,
                    stop=stop,
                    points=points,
                    nose_radius=radius,
                    geometry=geometry,
                    selection=selection,
                )
                shape = None
                sweep_configuration = {
                    "field": "upstream_mach",
                    "start": start,
                    "stop": stop,
                    "points": points,
                }
            configuration = make_configuration(
                calculator="detached_shock",
                mode=mode,
                inputs_si={"upstream_mach": mach, "nose_radius": radius},
                models={"geometry": geometry.value, "model": selection},
                units=preferences,
                sweep_si=sweep_configuration,
            )
        except ValueError as error:
            st.error(str(error), icon="🚫")
        else:
            st.session_state["detached_shock_payload"] = (result, configuration)
            if shape is not None:
                st.session_state["detached_shock_shape_result"] = shape

    payload = session_payload("detached_shock_payload", CalculationResult)
    if payload is None:
        with st.expander("モデルの前提・適用範囲"):
            st.write(
                "連続流・低温完全気体の経験相関です。Seiffは球形のみ、Billigの"
                "離脱距離は常にAmbrosio–Wortmanを使用します。"
            )
        return
    result, configuration = payload
    stored_shape = st.session_state.get("detached_shock_shape_result")
    shape = stored_shape if isinstance(stored_shape, BilligShockShapeResult) else None
    configured_models = configuration.get("models")
    configured_model = (
        configured_models.get("model") if isinstance(configured_models, dict) else None
    )

    def metrics(row: Mapping[str, object]) -> None:
        wanted: tuple[tuple[str, str], ...]
        if configured_model == "comparison":
            wanted = (
                ("AW Δ/Rn", "Ambrosio–Wortman Δ/Rn"),
                ("Seiff Δ/Rn", "Seiff Δ/Rn"),
                ("AW Δ", "Ambrosio–Wortman Δ ["),
                ("Seiff Δ", "Seiff Δ ["),
                ("Billig Rc", "Billig頂点曲率半径"),
            )
        else:
            wanted = (
                ("Δ/Rn", "選択モデル Δ/Rn"),
                ("Δ", "選択モデル Δ ["),
                ("AW Δ/Rn", "Ambrosio–Wortman Δ/Rn"),
                ("Billig Rc", "Billig頂点曲率半径"),
            )
        columns = st.columns(len(wanted))
        for column, (label, contains) in zip(columns, wanted, strict=True):
            with column:
                render_metric(row, label, contains)

    figures = (
        {"衝撃波形状": detached_shock_geometry(shape, preferences)}
        if shape is not None
        else detached_shock_trends(result.rows, preferences)
    )
    render_result_bundle(
        calculator="detached_shock",
        result=result,
        configuration=configuration,
        preferences=preferences,
        figures=figures,
        filename_prefix="aerophysics-detached-shock",
        metrics=metrics,
    )
    if shape is not None:
        st.download_button(
            "shock-shape x,y CSVをダウンロード",
            detached_shock_shape_csv(shape, preferences),
            file_name="aerophysics-detached-shock-shape.csv",
            mime="text/csv",
            key="detached_shock_shape_csv",
        )
    with st.expander("モデルの前提・適用範囲"):
        st.write(
            "座標原点はnose curvature center、x正方向は上流です。body vertexは"
            "x=Rn、shock vertexはx=Rn+Δです。平行afterbodyのためBillig双曲線の"
            "漸近角にはMach角を使用します。希薄流、実在気体補正、shock-fittingは"
            "含みません。"
        )
