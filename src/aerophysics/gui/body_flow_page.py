"""Single-point axisymmetric method-of-characteristics page."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import streamlit as st

from aerophysics.exceptions import ShockConvergenceError
from aerophysics.gui.body_flow_adapter import (
    body_flow_inputs_match,
    body_flow_net_rows,
    build_meridian,
    calculate_body_flow,
    resolve_freestream,
)
from aerophysics.gui.components import (
    calculation_button,
    finite_number,
    pop_pending_configuration,
    render_configuration_import,
    render_reset_button,
    render_result_bundle,
)
from aerophysics.gui.config import make_configuration
from aerophysics.gui.csv_inputs import (
    MeridianCSV,
    meridian_csv_template,
    parse_meridian_csv,
)
from aerophysics.gui.documentation import render_documentation_link
from aerophysics.gui.figures import body_flow_figures
from aerophysics.gui.tables import display_rows, rows_to_csv
from aerophysics.gui.units import QuantityKind, UnitPreferences, from_si, to_si

_DEFAULTS: dict[str, Any] = {
    "mach": 3.0,
    "pressure": 101325.0,
    "temperature": 288.15,
    "altitude": None,
    "nose_length": 0.3,
    "base_radius": 0.05,
    "ogive_segments": 32,
    "cylinder_length": 1.0,
    "boattail_length": None,
    "boattail_radius": None,
    "meridian_x": None,
    "meridian_r": None,
    "characteristic_count": 88,
    "fan_angle_step": math.pi / 180,
}
_MODELS = {
    "geometry_source": "parametric",
    "nose": "cone",
    "boattail": False,
    "freestream_basis": "static",
    "altitude_basis": "geometric",
    "pressure_quantity": "cp",
    "net_coloring": "none",
}
_QUANTITIES: dict[str, QuantityKind] = {
    "pressure": "pressure",
    "temperature": "temperature",
    "altitude": "length",
    "nose_length": "length",
    "base_radius": "length",
    "cylinder_length": "length",
    "boattail_length": "length",
    "boattail_radius": "length",
    "fan_angle_step": "angle",
}


def _preserve_units(preferences: UnitPreferences) -> None:
    stored = st.session_state.get("body_flow_units")
    if stored is not None:
        old = UnitPreferences.from_dict(stored)
        for field, kind in _QUANTITIES.items():
            key = f"body_flow_{field}"
            if key in st.session_state and getattr(old, kind) != getattr(
                preferences, kind
            ):
                value = to_si(st.session_state[key], kind, getattr(old, kind))
                st.session_state[key] = float(
                    from_si(value, kind, getattr(preferences, kind))
                )
    st.session_state["body_flow_units"] = preferences.to_dict()


def render_axisymmetric_body(preferences: UnitPreferences) -> None:
    """Render the editor and retain results until an explicit new calculation."""
    st.title("軸対称物体流れ（特性曲線法）")
    render_documentation_link("models/axisymmetric_characteristics.html")
    with st.expander("モデルの仮定・適用範囲"):
        st.write(
            "AIR（比熱比一定）の非粘性・断熱・迎角ゼロの外部流れです。"
            "鋭い先端、凸な折れ線母線、付着した弱い先端衝撃波、超音速状態が必要です。"
            "鈍頭・凹角・剥離・粘性排除・高温化学反応は扱いません。"
            "x は軸方向、r は半径、θ は +x 軸からの流れ角です。"
        )
    imported = pop_pending_configuration("axisymmetric_body")
    if imported is not None:
        # Clear all local widgets as well when a configuration is queued directly.
        for key in list(st.session_state):
            if isinstance(key, str) and key.startswith("body_flow_"):
                del st.session_state[key]
        st.session_state["body_flow_draft"] = imported
    render_configuration_import("axisymmetric_body", "body_flow")
    render_reset_button("body_flow", "body_flow_payload")
    _preserve_units(preferences)
    draft = st.session_state.get("body_flow_draft", {})
    inputs = {**_DEFAULTS, **draft.get("inputs_si", {})}
    models = {**_MODELS, **draft.get("models", {})}

    def choice(field: str, label: str, labels: dict[str, str]) -> str:
        options = tuple(labels)
        value = st.selectbox(
            label,
            options,
            index=options.index(models[field]),
            format_func=labels.__getitem__,
            key=f"body_flow_{field}",
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
            key=f"body_flow_{field}",
            min_value=low,
            max_value=high,
        )
        result = float(to_si(entered, kind, unit)) if kind else entered
        inputs[field] = result
        return result

    source = choice(
        "geometry_source", "形状の入力", {"parametric": "パラメトリック", "csv": "CSV"}
    )
    uploaded = None
    geometry_error: ValueError | None = None
    meridian = None
    if source == "parametric":
        nose = choice(
            "nose", "ノーズ形状", {"cone": "円錐", "tangent_ogive": "接線オジーブ"}
        )
        number("nose_length", "ノーズ長 L_n", minimum=1e-12)
        number("base_radius", "基底半径 R", minimum=1e-12)
        if nose == "tangent_ogive":
            inputs["ogive_segments"] = st.number_input(
                "オジーブの弦数 N",
                2,
                200,
                int(inputs["ogive_segments"]),
                key="body_flow_ogive_segments",
            )
            st.caption(
                "オジーブは円弧角が等間隔の弦で表します。N の細分化で結果が変わります。"
            )
        number("cylinder_length", "円筒長 L_cyl", minimum=0.0)
        models["boattail"] = st.checkbox(
            "ボートテールを付ける", value=models["boattail"], key="body_flow_boattail"
        )
        if models["boattail"]:
            number("boattail_length", "ボートテール長 L_bt", 0.2, minimum=1e-12)
            number("boattail_radius", "後端半径 R_bt", 0.025, minimum=1e-12)
        else:
            inputs["boattail_length"] = inputs["boattail_radius"] = None
        inputs["meridian_x"] = inputs["meridian_r"] = None
        try:
            meridian = build_meridian(
                nose=nose,
                nose_length=inputs["nose_length"],
                base_radius=inputs["base_radius"],
                ogive_segments=inputs["ogive_segments"],
                cylinder_length=inputs["cylinder_length"],
                boattail=models["boattail"],
                boattail_length=inputs["boattail_length"],
                boattail_radius=inputs["boattail_radius"],
            )
        except ValueError as error:
            geometry_error = error
    else:
        left, right = st.columns(2)
        uploaded = left.file_uploader(
            f"母線 CSV: x,r [{preferences.length}]",
            type="csv",
            key="body_flow_meridian_upload",
        )
        right.download_button(
            "母線CSVテンプレート",
            meridian_csv_template(),
            "aerophysics-meridian.csv",
            "text/csv",
            key="body_flow_template",
        )
        st.caption(
            "x,r はアップロード時の表示長さ単位で読み込み、"
            "その後に表示単位を変えても形状は変わりません。"
            "先頭は (0,0)、x は単調増加、以降の r は正。"
            "凸性は計算時に検証します。"
        )
        try:
            if uploaded is not None:
                # Interpret a file once, in the length unit active when it was
                # uploaded; later display-unit changes must not rescale it.
                cached = st.session_state.get("body_flow_meridian_upload_si")
                if cached is not None and cached[0] == uploaded.file_id:
                    meridian = MeridianCSV(np.asarray(cached[1]), np.asarray(cached[2]))
                else:
                    st.session_state.pop("body_flow_meridian_upload_si", None)
                    meridian = parse_meridian_csv(uploaded.getvalue(), preferences)
                    st.session_state["body_flow_meridian_upload_si"] = (
                        uploaded.file_id,
                        meridian.axial_distance.tolist(),
                        meridian.radius.tolist(),
                    )
                inputs["meridian_x"] = meridian.axial_distance.tolist()
                inputs["meridian_r"] = meridian.radius.tolist()
            elif inputs["meridian_x"] is not None and inputs["meridian_r"] is not None:
                meridian = MeridianCSV(
                    np.asarray(inputs["meridian_x"]), np.asarray(inputs["meridian_r"])
                )
                st.caption("設定JSONに保存された母線の頂点を使用します。")
        except ValueError as error:
            geometry_error = error
    if (
        meridian is not None
        and len(meridian.axial_distance) >= 2
        and len(meridian.radius) >= 2
    ):
        tip = math.atan2(meridian.radius[1], meridian.axial_distance[1])
        tip_shown = float(from_si(tip, "angle", preferences.angle))
        text = (
            f"最初の区間の半頂角（有効先端円錐角）: {tip_shown:.6g} {preferences.angle}"
        )
        if source == "parametric" and models["nose"] == "tangent_ogive":
            rho = (inputs["base_radius"] ** 2 + inputs["nose_length"] ** 2) / (
                2 * inputs["base_radius"]
            )
            exact = math.asin(inputs["nose_length"] / rho)
            exact_shown = float(from_si(exact, "angle", preferences.angle))
            text += f"、正確な先端半頂角: {exact_shown:.6g} {preferences.angle}"
        st.info(text)

    basis = choice(
        "freestream_basis",
        "自由流条件",
        {"static": "静圧・静温・Mach", "atmosphere": "標準大気・Mach"},
    )
    number("mach", "自由流 Mach M∞", minimum=1.000000001)
    freestream_error = None
    if basis == "static":
        number("pressure", "静圧 p∞", minimum=1e-12)
        number("temperature", "静温 T∞", minimum=1e-12)
        inputs["altitude"] = None
    else:
        choice(
            "altitude_basis",
            "高度の種類",
            {"geometric": "幾何高度", "geopotential": "ジオポテンシャル高度"},
        )
        number("altitude", "高度", 10000.0)
        st.caption(
            "US Standard Atmosphere 1976、幾何高度 −5～86 km。"
            "入力高度を選択した種類で解釈します。"
        )
        try:
            inputs["pressure"], inputs["temperature"] = resolve_freestream(
                freestream_basis=basis,
                altitude=inputs["altitude"],
                altitude_basis=models["altitude_basis"],
            )
            p_shown = float(
                from_si(inputs["pressure"], "pressure", preferences.pressure)
            )
            t_shown = float(
                from_si(inputs["temperature"], "temperature", preferences.temperature)
            )
            st.info(
                f"静圧: {p_shown:.6g} {preferences.pressure}、"
                f"静温: {t_shown:.6g} {preferences.temperature}"
            )
        except ValueError as error:
            freestream_error = error
    with st.expander("詳細設定（解像度）"):
        inputs["characteristic_count"] = st.number_input(
            "特性線本数",
            4,
            1000,
            int(inputs["characteristic_count"]),
            key="body_flow_characteristic_count",
        )
        number(
            "fan_angle_step",
            "ファン角刻み",
            minimum=math.pi / 3600,
            maximum=math.pi / 36,
        )
        st.caption(
            "収束確認は特性線本数を2倍、角刻みを半分にします。長い物体と粗いファンでは角刻みの縮小が必要なエラーが生じます。"
        )
    submitted = calculation_button("body_flow", key="body_flow_calculate")
    quantity = choice(
        "pressure_quantity",
        "表面圧の表示量",
        {"cp": "Cp", "pressure_rise": "(p−p∞)/p∞"},
    )
    net_coloring = choice(
        "net_coloring",
        "特性曲線網の節点色",
        {"none": "なし", "mach": "Mach", "cp": "Cp"},
    )
    st.session_state["body_flow_draft"] = {"inputs_si": inputs, "models": models}
    if submitted:
        st.session_state.pop("body_flow_payload", None)
        st.session_state.pop("body_flow_net_csv_cache", None)
        try:
            if geometry_error is not None:
                raise geometry_error
            if freestream_error is not None:
                raise freestream_error
            if meridian is None:
                raise ValueError("母線CSVを選択してください。")
            with st.spinner("物体流れを計算中…"):
                calculation = calculate_body_flow(
                    meridian,
                    mach=inputs["mach"],
                    pressure=inputs["pressure"],
                    temperature=inputs["temperature"],
                    characteristic_count=inputs["characteristic_count"],
                    fan_angle_step=inputs["fan_angle_step"],
                )
            configuration = make_configuration(
                calculator="axisymmetric_body",
                mode="single",
                inputs_si=dict(inputs),
                models=dict(models),
                units=preferences,
            )
            st.session_state["body_flow_payload"] = (calculation, configuration)
        except (ValueError, ShockConvergenceError) as error:
            st.error(str(error))
    payload = st.session_state.get("body_flow_payload")
    if payload is None:
        st.session_state.pop("body_flow_net_csv_cache", None)
        return
    calculation, configuration = payload
    if not body_flow_inputs_match(inputs, configuration["inputs_si"]) or any(
        models[k] != configuration["models"][k]
        for k in models
        if k not in {"pressure_quantity", "net_coloring"}
    ):
        st.warning(
            "入力が変更されています。以下は前回計算の結果です。再計算してください。"
        )
    raw = calculation.raw
    rows = calculation.surface.rows
    values = (
        ("先端表面 p/p∞", rows[0]["pressure_ratio"]),
        (
            f"先端衝撃波角 [{preferences.angle}]",
            float(from_si(raw.shock_angle[0], "angle", preferences.angle)),
        ),
        (
            f"衝撃波後の全圧 [{preferences.pressure}]",
            float(from_si(raw.total_pressure, "pressure", preferences.pressure)),
        ),
        (
            f"全温 [{preferences.temperature}]",
            float(
                from_si(raw.total_temperature, "temperature", preferences.temperature)
            ),
        ),
        ("最小 Cp", min(float(row["cp"]) for row in rows)),
        ("物体後端 Cp", rows[-1]["cp"]),
    )
    for column, (label, value) in zip(st.columns(6), values, strict=True):
        column.metric(label, f"{value:.6g}")
    exported = {
        **configuration,
        "models": {
            **configuration["models"],
            "pressure_quantity": quantity,
            "net_coloring": net_coloring,
        },
    }
    figures = body_flow_figures(
        rows,
        calculation.shock_rows,
        preferences,
        pressure_quantity=quantity,
        characteristic_net=raw.characteristic_net,
        net_coloring=net_coloring,
        freestream_pressure=configuration["inputs_si"]["pressure"],
        freestream_mach=configuration["inputs_si"]["mach"],
        heat_capacity_ratio=raw.gas.heat_capacity_ratio,
    )
    render_result_bundle(
        calculator="axisymmetric_body",
        result=calculation.surface,
        configuration=exported,
        preferences=preferences,
        figures=figures,
        filename_prefix="body_flow_results",
    )
    if figures["特性曲線網"].layout.meta["downsampled"]:
        st.caption(
            "特性曲線網は表示用に間引いています（C+ 最大150本、C− 最大80本、"
            "節点色 最大20,000点）。"
            "CSVには全節点を保存します。"
        )
    csv_cache = st.session_state.get("body_flow_net_csv_cache")
    display_units = preferences.to_dict()
    if (
        csv_cache is None
        or csv_cache[0] is not calculation
        or csv_cache[1] != display_units
    ):
        net_rows = body_flow_net_rows(
            raw,
            pressure=configuration["inputs_si"]["pressure"],
            mach=configuration["inputs_si"]["mach"],
        )
        csv_text = rows_to_csv(
            display_rows("axisymmetric_body_net", net_rows, preferences)
        )
        # Retain identity, not a bare id that could be reused after collection.
        csv_cache = (calculation, display_units, csv_text)
        st.session_state["body_flow_net_csv_cache"] = csv_cache
    st.download_button(
        "特性曲線網CSVをダウンロード",
        csv_cache[2],
        "aerophysics-body-net.csv",
        "text/csv",
        key="body_flow_net_csv",
    )
    st.subheader("衝撃波の軌跡")
    shock_rows = display_rows(
        "axisymmetric_body_shock", calculation.shock_rows, preferences
    )
    st.dataframe(shock_rows, hide_index=True, width="stretch")
    st.download_button(
        "衝撃波CSVをダウンロード",
        rows_to_csv(shock_rows),
        "aerophysics-body-shock.csv",
        "text/csv",
        key="body_flow_shock_csv",
    )
    st.info(
        "表面結果は凸な頂点でファンの下流状態だけを保持し、特性曲線網は中間状態も含みます。衝撃波の軌跡は物体後端より手前で終わる場合も、先まで延びる場合もあります。特性曲線網は数値離散化であり、収束済みの流れ場を保証しません。実験に対する精度は確立されていません。"
    )
    if st.button("現在の物体流れケースとして保存", key="body_flow_save_case"):
        st.session_state["current_axisymmetric_body_case"] = raw
        st.success("現在の物体流れケースとして保存しました。")
