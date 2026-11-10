"""Single-case marching boundary-layer page."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import streamlit as st

from aerophysics.exceptions import (
    BoundaryLayerConvergenceError,
    BoundaryLayerSeparationError,
)
from aerophysics.gui.body_flow_adapter import body_flow_inputs_match
from aerophysics.gui.components import (
    PLOTLY_CONFIG,
    calculation_button,
    finite_number,
    pop_pending_configuration,
    render_configuration_import,
    render_reset_button,
    render_result_bundle,
)
from aerophysics.gui.config import make_configuration
from aerophysics.gui.csv_inputs import edge_csv_template, parse_edge_csv
from aerophysics.gui.documentation import render_documentation_link
from aerophysics.gui.figures import (
    marching_boundary_layer_figures,
    marching_boundary_layer_profile_figure,
)
from aerophysics.gui.marching_adapter import (
    EdgeDefinition,
    body_case_edge,
    calculate_marching_boundary_layer,
    profile_rows,
    total_conditions,
    uniform_edge,
)
from aerophysics.gui.tables import display_rows, rows_to_csv
from aerophysics.gui.units import QuantityKind, UnitPreferences, from_si, to_si

_TOTAL_PRESSURE, _TOTAL_TEMPERATURE = total_conditions(2.0, 10000.0, 220.0)
_DEFAULTS: dict[str, Any] = {
    "mach": 2.0,
    "pressure": 10000.0,
    "temperature": 220.0,
    "length": 1.0,
    "station_count": 50,
    "cone_half_angle": None,
    "total_pressure": _TOTAL_PRESSURE,
    "total_temperature": _TOTAL_TEMPERATURE,
    "edge_distance": None,
    "edge_pressure": None,
    "edge_radius": None,
    "wall_temperature": None,
    "transition_distance": None,
    "prandtl_number": 0.72,
    "turbulent_prandtl_number": 0.9,
    "edge_sa_ratio": 3.0,
}
_MODELS = {
    "edge_source": "uniform",
    "geometry": "planar",
    "regime": "laminar",
    "wall": "adiabatic",
}
_QUANTITIES: dict[str, QuantityKind] = {
    "pressure": "pressure",
    "temperature": "temperature",
    "length": "length",
    "cone_half_angle": "angle",
    "total_pressure": "pressure",
    "total_temperature": "temperature",
    "wall_temperature": "temperature",
    "transition_distance": "length",
}
_BODY_CASE_KEY = "current_axisymmetric_body_case"


def _preserve_units(preferences: UnitPreferences) -> None:
    stored = st.session_state.get("march_bl_units")
    if stored is not None:
        old = UnitPreferences.from_dict(stored)
        for field, kind in _QUANTITIES.items():
            key = f"march_bl_{field}"
            if key in st.session_state and getattr(old, kind) != getattr(
                preferences, kind
            ):
                value = to_si(st.session_state[key], kind, getattr(old, kind))
                st.session_state[key] = float(
                    from_si(value, kind, getattr(preferences, kind))
                )
    st.session_state["march_bl_units"] = preferences.to_dict()


def _embedded_edge(inputs: dict[str, Any]) -> EdgeDefinition | None:
    if inputs["edge_distance"] is None or inputs["edge_pressure"] is None:
        return None
    radius = inputs["edge_radius"]
    return EdgeDefinition(
        np.asarray(inputs["edge_distance"], dtype=np.float64),
        np.asarray(inputs["edge_pressure"], dtype=np.float64),
        None if radius is None else np.asarray(radius, dtype=np.float64),
        inputs["total_pressure"],
        inputs["total_temperature"],
    )


def _embed(inputs: dict[str, Any], definition: EdgeDefinition) -> None:
    inputs["edge_distance"] = definition.distance.tolist()
    inputs["edge_pressure"] = definition.pressure.tolist()
    inputs["edge_radius"] = (
        None if definition.radius is None else definition.radius.tolist()
    )
    inputs["total_pressure"] = definition.total_pressure
    inputs["total_temperature"] = definition.total_temperature


def render_marching_boundary_layer(preferences: UnitPreferences) -> None:
    """Render the editor and retain results until an explicit new calculation."""
    st.title("境界層マーチング解析")
    render_documentation_link("models/marching_boundary_layer.html")
    with st.expander("モデルの仮定・適用範囲"):
        st.write(
            "鋭い前縁（先端）から壁面に沿って薄層境界層方程式を下流へマーチングします。"
            "気体は AIR（比熱比一定）、粘性は Sutherland、Prandtl 数は一定です。"
            "エッジは全エンタルピー一定の等エントロピー流れで、外部流れ自体は"
            "計算しません。層流、または SA-noft2 乱流モデルを扱います。"
            "指定遷移は SA 生成項を遷移位置まで抑制するだけで、遷移位置や"
            "遷移長さの予測ではありません。"
            "壁面曲率・法線方向圧力勾配・横曲率の補正・排除厚さによる外部流れの"
            "変化・鈍頭のエントロピー層・高温気体効果は扱いません。"
            "剥離点より下流は計算できません。"
            "s は前縁からの壁面に沿う距離、Cf は局所エッジ動圧基準、"
            "qw は壁へ入る向きが正、β>0 は加速流です。"
        )
    imported = pop_pending_configuration("marching_boundary_layer")
    if imported is not None:
        # Clear all local widgets as well when a configuration is queued directly.
        for key in list(st.session_state):
            if isinstance(key, str) and key.startswith("march_bl_"):
                del st.session_state[key]
        st.session_state["march_bl_draft"] = imported
        # A loaded edge stays in use until the user selects the saved case again.
        st.session_state["march_bl_keep_embedded"] = True
    render_configuration_import("marching_boundary_layer", "march_bl")
    render_reset_button("march_bl", "march_bl_payload")
    _preserve_units(preferences)
    draft = st.session_state.get("march_bl_draft", {})
    inputs = {**_DEFAULTS, **draft.get("inputs_si", {})}
    models = {**_MODELS, **draft.get("models", {})}
    previous_source = models["edge_source"]

    def choice(field: str, label: str, labels: dict[str, str]) -> str:
        options = tuple(labels)
        value = st.selectbox(
            label,
            options,
            index=options.index(models[field]),
            format_func=labels.__getitem__,
            key=f"march_bl_{field}",
        )
        models[field] = value
        return str(value)

    def number(
        field: str,
        label: str,
        fallback: float | None = None,
        minimum: float | None = None,
        maximum: float | None = None,
    ) -> float:
        kind = _QUANTITIES.get(field)
        value = inputs[field] if inputs[field] is not None else fallback
        assert value is not None
        unit = str(getattr(preferences, kind)) if kind else ""
        shown = float(from_si(value, kind, unit)) if kind else float(value)
        low = (
            float(from_si(minimum, kind, unit))
            if kind and minimum is not None
            else minimum
        )
        high = (
            float(from_si(maximum, kind, unit))
            if kind and maximum is not None
            else maximum
        )
        entered = finite_number(
            f"{label} [{unit}]" if unit else label,
            shown,
            key=f"march_bl_{field}",
            min_value=low,
            max_value=high,
        )
        result = float(to_si(entered, kind, unit)) if kind else entered
        inputs[field] = result
        return result

    st.subheader("エッジ条件")
    source = choice(
        "edge_source",
        "エッジ条件の入力",
        {
            "uniform": "一様エッジ（平板・円錐）",
            "csv": "圧力分布CSV",
            "body": "保存した軸対称物体流れケース",
        },
    )
    if source != previous_source:
        # Arrays embedded for another source never carry over silently.
        inputs["edge_distance"] = inputs["edge_pressure"] = None
        inputs["edge_radius"] = None
    definition: EdgeDefinition | None = None
    edge_error: ValueError | None = None
    if source == "uniform":
        geometry = choice(
            "geometry", "形状", {"planar": "平板（二次元）", "cone": "鋭い円錐"}
        )
        number("mach", "エッジ Mach M_e", minimum=0.01)
        number("pressure", "エッジ静圧 p_e", minimum=1e-12)
        number("temperature", "エッジ静温 T_e", minimum=1e-12)
        number("length", "壁面長さ L", minimum=1e-12)
        if geometry == "cone":
            number(
                "cone_half_angle",
                "円錐半頂角 θ_c",
                math.radians(10.0),
                minimum=math.pi / 1800,
                maximum=89 * math.pi / 180,
            )
            st.caption(
                "エッジ条件には円錐衝撃波後の円錐表面の値を入力します"
                "（円錐衝撃波ページで求められます）。s は先端からの母線に沿う距離です。"
            )
        else:
            inputs["cone_half_angle"] = None
        inputs["station_count"] = st.number_input(
            "出力点数（等間隔）",
            2,
            1000,
            int(inputs["station_count"]),
            key="march_bl_station_count",
        )
        inputs["edge_distance"] = inputs["edge_pressure"] = None
        inputs["edge_radius"] = None
        try:
            definition = uniform_edge(
                mach=inputs["mach"],
                pressure=inputs["pressure"],
                temperature=inputs["temperature"],
                length=inputs["length"],
                station_count=inputs["station_count"],
                geometry=geometry,
                cone_half_angle=inputs["cone_half_angle"],
            )
            inputs["total_pressure"] = definition.total_pressure
            inputs["total_temperature"] = definition.total_temperature
        except ValueError as error:
            edge_error = error
    elif source == "csv":
        left, right = st.columns(2)
        uploaded = left.file_uploader(
            f"エッジ圧力 CSV: s [{preferences.length}], p [{preferences.pressure}]",
            type="csv",
            key="march_bl_edge_upload",
        )
        right.download_button(
            "エッジ圧力CSVテンプレート",
            edge_csv_template(),
            "aerophysics-boundary-layer-edge.csv",
            "text/csv",
            key="march_bl_template",
        )
        st.caption(
            "列 s,p は必須で、s は 0 から始まり単調増加、p は全圧未満の正の値です。"
            "任意の列 r（壁面半径、s と同じ長さ単位）を加えると軸対称として計算します。"
            "値はアップロード時の表示単位で読み込み、その後に表示単位を変えても"
            "変わりません。圧力変化は十分な点数で滑らかに与えてください。"
        )
        number("total_pressure", "エッジ全圧 p₀", minimum=1e-12)
        number("total_temperature", "全温 T₀", minimum=1e-12)
        try:
            if uploaded is not None:
                # Interpret a file once, in the units active when it was uploaded.
                cached = st.session_state.get("march_bl_edge_upload_si")
                if cached is None or cached[0] != uploaded.file_id:
                    st.session_state.pop("march_bl_edge_upload_si", None)
                    parsed = parse_edge_csv(uploaded.getvalue(), preferences)
                    cached = (
                        uploaded.file_id,
                        parsed.distance.tolist(),
                        parsed.pressure.tolist(),
                        None if parsed.radius is None else parsed.radius.tolist(),
                    )
                    st.session_state["march_bl_edge_upload_si"] = cached
                inputs["edge_distance"] = cached[1]
                inputs["edge_pressure"] = cached[2]
                inputs["edge_radius"] = cached[3]
            elif inputs["edge_distance"] is not None:
                st.caption("設定JSONに保存されたエッジ分布を使用します。")
            definition = _embedded_edge(inputs)
        except ValueError as error:
            edge_error = error
    else:
        case = st.session_state.get(_BODY_CASE_KEY)
        embedded = _embedded_edge(inputs)
        keep = (
            bool(st.session_state.get("march_bl_keep_embedded"))
            and embedded is not None
        )
        if case is not None and keep:
            st.caption(
                "設定JSONに保存されたエッジ分布を使用中です。"
                "現在保存されている物体流れケースとは異なる場合があります。"
            )
            if st.button("現在の物体流れケースに切り替える", key="march_bl_use_case"):
                st.session_state["march_bl_keep_embedded"] = keep = False
        if case is not None and not keep:
            try:
                definition = body_case_edge(case)
                _embed(inputs, definition)
            except ValueError as error:
                edge_error = error
        elif embedded is not None:
            definition = embedded
            if case is None:
                st.caption("設定JSONに保存されたエッジ分布を使用します。")
        else:
            edge_error = ValueError(
                "軸対称物体流れ（特性曲線法）ページで計算し、"
                "「現在の物体流れケースとして保存」を押してください。"
            )
            st.info(str(edge_error))
        st.caption(
            "物体表面流線の圧力・半径と衝撃波後の全圧を使用します。"
            "凸角は角の位置の 0.1% の長さの圧力ランプに置き換わり、"
            "角の近傍（δ99 の数十倍以内）の結果は信頼できません。"
        )
    wall_length = (
        float(definition.distance[-1])
        if definition is not None and definition.distance.size
        else None
    )
    if definition is not None and wall_length is not None:
        p0 = float(from_si(definition.total_pressure, "pressure", preferences.pressure))
        t0 = float(
            from_si(
                definition.total_temperature, "temperature", preferences.temperature
            )
        )
        end = float(from_si(wall_length, "length", preferences.length))
        st.info(
            f"エッジ点数: {definition.distance.size}、"
            f"壁面長さ: {end:.6g} {preferences.length}、"
            f"全圧: {p0:.6g} {preferences.pressure}、"
            f"全温: {t0:.6g} {preferences.temperature}、"
            + ("軸対称" if definition.radius is not None else "二次元")
        )

    st.subheader("境界層の状態と壁面条件")
    regime = choice(
        "regime",
        "境界層の状態",
        {
            "laminar": "層流",
            "turbulent": "乱流（SA-noft2、前縁から）",
            "transitional": "指定位置で遷移（SA-noft2）",
        },
    )
    if regime == "transitional":
        number(
            "transition_distance",
            "遷移位置 s_tr",
            0.5 * (wall_length if wall_length is not None else 1.0),
            minimum=1e-12,
        )
        st.caption(
            "0 < s_tr < 壁面長さ。遷移位置は予測せず、入力値をそのまま使います。"
        )
    else:
        inputs["transition_distance"] = None
    wall = choice("wall", "壁面条件", {"adiabatic": "断熱壁", "isothermal": "等温壁"})
    if wall == "isothermal":
        number("wall_temperature", "壁温 Tw", 300.0, minimum=1e-12)
    else:
        inputs["wall_temperature"] = None
    with st.expander("詳細設定（Prandtl数・SA流入値）"):
        number("prandtl_number", "Prandtl数 Pr", minimum=1e-12)
        number("turbulent_prandtl_number", "乱流Prandtl数 Pr_t", minimum=1e-12)
        number("edge_sa_ratio", "エッジの ν̃/ν_e", minimum=0.0)
        st.caption(
            "Pr_t と ν̃/ν_e は乱流・指定遷移でのみ使われます。"
            "指定遷移では ν̃/ν_e > 0 が必要です。"
        )
    submitted = calculation_button("march_bl", key="march_bl_calculate")
    st.session_state["march_bl_draft"] = {"inputs_si": inputs, "models": models}
    if submitted:
        st.session_state.pop("march_bl_payload", None)
        st.session_state.pop("march_bl_profile_station", None)
        try:
            if edge_error is not None:
                raise edge_error
            if definition is None:
                raise ValueError("エッジ圧力CSVを選択してください。")
            with st.spinner("境界層を計算中…"):
                calculation = calculate_marching_boundary_layer(
                    definition,
                    regime=regime,
                    wall_temperature=inputs["wall_temperature"],
                    transition_distance=inputs["transition_distance"],
                    prandtl_number=inputs["prandtl_number"],
                    turbulent_prandtl_number=inputs["turbulent_prandtl_number"],
                    edge_sa_ratio=inputs["edge_sa_ratio"],
                )
            configuration = make_configuration(
                calculator="marching_boundary_layer",
                mode="single",
                inputs_si=dict(inputs),
                models=dict(models),
                units=preferences,
            )
            st.session_state["march_bl_payload"] = (calculation, configuration)
        except BoundaryLayerSeparationError as error:
            st.error(f"最初の出力点より上流で剥離したため、結果はありません: {error}")
        except (ValueError, BoundaryLayerConvergenceError) as error:
            st.error(str(error))
    payload = st.session_state.get("march_bl_payload")
    if payload is None:
        return
    calculation, configuration = payload
    if not body_flow_inputs_match(inputs, configuration["inputs_si"]) or any(
        models[k] != configuration["models"][k] for k in models
    ):
        st.warning(
            "入力が変更されています。以下は前回計算の結果です。再計算してください。"
        )
    rows = calculation.stations.rows
    separation = calculation.separation_distance
    if separation is not None:
        shown = float(from_si(separation, "length", preferences.length))
        st.error(
            f"逆圧力勾配により s = {shown:.6g} {preferences.length} で剥離しました。"
            "剥離点より上流の結果だけを表示します。"
        )
    adiabatic = configuration["models"]["wall"] == "adiabatic"
    last = rows[-1]

    def shown_value(key: str, kind: QuantityKind) -> float:
        return float(from_si(last[key], kind, str(getattr(preferences, kind))))

    values = (
        ("最終点 Re_s", last["reynolds_number"]),
        (
            f"δ99 [{preferences.length}]",
            shown_value("boundary_layer_thickness", "length"),
        ),
        (f"δ* [{preferences.length}]", shown_value("displacement_thickness", "length")),
        (f"θ [{preferences.length}]", shown_value("momentum_thickness", "length")),
        ("Cf", last["local_skin_friction_coefficient"]),
        (
            f"断熱壁温 [{preferences.temperature}]",
            shown_value("wall_temperature", "temperature"),
        )
        if adiabatic
        else (
            f"qw [{preferences.heat_flux}]",
            shown_value("wall_heat_flux", "heat_flux"),
        ),
    )
    st.caption("以下の指標は計算できた最も下流の点の値です。")
    for column, (label, value) in zip(st.columns(6), values, strict=True):
        column.metric(label, f"{value:.6g}")
    render_result_bundle(
        calculator="marching_boundary_layer",
        result=calculation.stations,
        configuration=configuration,
        preferences=preferences,
        figures=marching_boundary_layer_figures(
            rows,
            preferences,
            adiabatic=adiabatic,
            transition_distance=configuration["inputs_si"]["transition_distance"],
            separation_distance=separation,
        ),
        filename_prefix="marching_boundary_layer_results",
    )
    st.subheader("境界層プロファイル")

    def station_label(index: int) -> str:
        distance = rows[index]["surface_distance"]
        return f"{float(from_si(distance, 'length', preferences.length)):.6g}"

    station = st.selectbox(
        f"表示位置 s [{preferences.length}]",
        range(len(rows)),
        index=len(rows) - 1,
        format_func=station_label,
        key="march_bl_profile_station",
    )
    profile = profile_rows(calculation, int(station))
    st.plotly_chart(
        marching_boundary_layer_profile_figure(
            profile,
            preferences,
            boundary_layer_thickness=rows[station]["boundary_layer_thickness"],
        ),
        width="stretch",
        config=PLOTLY_CONFIG,
    )
    profile_table = display_rows(
        "marching_boundary_layer_profile", profile, preferences
    )
    st.dataframe(profile_table, hide_index=True, width="stretch")
    st.download_button(
        "プロファイルCSVをダウンロード",
        rows_to_csv(profile_table),
        "aerophysics-boundary-layer-profile.csv",
        "text/csv",
        key="march_bl_profile_csv",
    )
    st.info(
        "プロファイルの図は 2δ99 までを表示し、表とCSVは数値計算領域の全点を含みます。"
        "SA 変数 ν̃ の外層分布は強い膨張や小さい ν̃/ν_e では格子収束していません。"
        "検証した条件と限界は「理論・モデルの説明」を参照してください。"
    )
