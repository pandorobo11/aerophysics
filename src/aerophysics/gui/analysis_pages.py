"""Streamlit pages for advanced viscous-flow and property analyses."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np
import streamlit as st

from aerophysics.boundary_layer_profile import (
    CompressibleVelocityTransformation,
    TemperatureVelocityRelation,
)
from aerophysics.gui.adapters import CalculationResult
from aerophysics.gui.advanced_adapters import (
    BoundaryLayerCase,
    BoundaryProfileCase,
    ProfileCalculation,
    boundary_layer_profiles,
    protrusion_condition,
    protrusion_sweep,
    thermochemistry_condition,
    thermochemistry_sweep,
    viscosity_condition,
    viscosity_sweep,
)
from aerophysics.gui.components import (
    finite_number,
    pop_pending_configuration,
    render_configuration_import,
    render_reset_button,
    render_result_bundle,
)
from aerophysics.gui.config import make_configuration
from aerophysics.gui.csv_inputs import (
    ProfileCSV,
    ShapeCSV,
    parse_profile_csv,
    parse_shape_csv,
    profile_csv_template,
    shape_csv_template,
)
from aerophysics.gui.figures import (
    boundary_layer_profile_figures,
    protrusion_figures,
    protrusion_shape_figure,
    thermochemistry_figures,
    viscosity_figures,
)
from aerophysics.gui.units import UnitPreferences, from_si, to_si


def _defaults(
    configuration: dict[str, object] | None,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    if configuration is None:
        return {}, {}, {}
    return tuple(
        dict(value) if isinstance(value, dict) else {}
        for value in (
            configuration.get("inputs_si"),
            configuration.get("models"),
            configuration.get("sweep_si"),
        )
    )  # type: ignore[return-value]


def _display(value: float, kind: str, unit: str) -> float:
    return float(from_si(value, kind, unit))  # type: ignore[arg-type]


def _si(value: float, kind: str, unit: str) -> float:
    return float(to_si(value, kind, unit))  # type: ignore[arg-type]


def _number(values: Mapping[str, Any], key: str, default: float) -> float:
    value = values.get(key, default)
    return float(value) if isinstance(value, (int, float)) else default


def _array(values: Mapping[str, Any], key: str) -> np.ndarray | None:
    value = values.get(key)
    if not isinstance(value, list):
        return None
    return np.asarray(value, dtype=np.float64)


def _metric(row: Mapping[str, object], label: str, contains: str) -> None:
    heading = next((name for name in row if contains in name), "")
    value = row.get(heading)
    st.metric(label, f"{float(value):.5g}" if isinstance(value, (int, float)) else "—")


def _profile_payload() -> tuple[ProfileCalculation, dict[str, object]] | None:
    value = st.session_state.get("profile_payload")
    if (
        isinstance(value, tuple)
        and len(value) == 2
        and isinstance(value[0], ProfileCalculation)
        and isinstance(value[1], dict)
    ):
        return value
    return None


def _result_payload(
    key: str,
) -> tuple[CalculationResult, dict[str, object]] | None:
    value = st.session_state.get(key)
    if (
        isinstance(value, tuple)
        and len(value) == 2
        and isinstance(value[0], CalculationResult)
        and isinstance(value[1], dict)
    ):
        return value
    return None


_TRANSFORM_LABELS = {
    CompressibleVelocityTransformation.VAN_DRIEST: "Van Driest",
    CompressibleVelocityTransformation.VOLPIANI: "Volpiani",
}


def render_boundary_layer_profile(preferences: UnitPreferences) -> None:
    """Render compressible Spalding-Coles boundary-layer profiles."""
    st.title("圧縮性境界層プロファイル")
    st.caption("滑面・乱流ZPG境界層の平均速度と熱物性分布を予測します。")
    imported = pop_pending_configuration("boundary_layer_profile")
    inputs, models, _ = _defaults(imported)
    render_configuration_import("boundary_layer_profile", "profile")
    render_reset_button("profile", "profile_payload")
    linked = st.session_state.get("current_boundary_layer_case")
    has_linked = isinstance(linked, BoundaryLayerCase)

    with st.form("profile_form"):
        sources = ("manual", "boundary") if has_linked else ("manual",)
        default_source = str(models.get("source", "manual"))
        source = st.radio(
            "入力元",
            sources,
            index=sources.index(default_source) if default_source in sources else 0,
            format_func=lambda value: (
                "手入力" if value == "manual" else "現在の乱流平板境界層ケース"
            ),
            horizontal=True,
            key="profile_source",
        )
        case = linked if source == "boundary" and has_linked else None
        edge_velocity_si = (
            case.edge_velocity
            if isinstance(case, BoundaryLayerCase)
            else _number(inputs, "edge_velocity", 300.0)
        )
        edge_density_si = (
            case.edge_density
            if isinstance(case, BoundaryLayerCase)
            else _number(inputs, "edge_density", 1.0)
        )
        edge_temperature_si = (
            case.edge_temperature
            if isinstance(case, BoundaryLayerCase)
            else _number(inputs, "edge_temperature", 300.0)
        )
        thickness_si = (
            case.boundary_layer_thickness
            if isinstance(case, BoundaryLayerCase)
            else _number(inputs, "boundary_layer_thickness", 0.05)
        )
        shear_si = (
            case.wall_shear_stress
            if isinstance(case, BoundaryLayerCase)
            else _number(inputs, "wall_shear_stress", 85.0)
        )
        disabled = isinstance(case, BoundaryLayerCase)
        edge_velocity = finite_number(
            f"外縁速度 U_e [{preferences.speed}]",
            _display(edge_velocity_si, "speed", preferences.speed),
            key="profile_velocity",
            min_value=1.0e-12,
            disabled=disabled,
        )
        edge_density = finite_number(
            f"外縁密度 ρ_e [{preferences.density}]",
            _display(edge_density_si, "density", preferences.density),
            key="profile_density",
            min_value=1.0e-12,
            disabled=disabled,
        )
        edge_temperature = finite_number(
            f"外縁温度 T_e [{preferences.temperature}]",
            _display(edge_temperature_si, "temperature", preferences.temperature),
            key="profile_temperature",
            disabled=disabled,
        )
        thickness = finite_number(
            f"境界層厚さ δ₉₉ [{preferences.length}]",
            _display(thickness_si, "length", preferences.length),
            key="profile_thickness",
            min_value=1.0e-12,
            disabled=disabled,
        )
        shear = finite_number(
            f"壁面せん断応力 τ_w [{preferences.pressure}]",
            _display(shear_si, "pressure", preferences.pressure),
            key="profile_shear",
            min_value=1.0e-12,
            disabled=disabled,
        )
        comparison = st.selectbox(
            "速度変換",
            ("van_driest", "volpiani", "compare"),
            index=("van_driest", "volpiani", "compare").index(
                str(models.get("transformation", "compare"))
            ),
            format_func=lambda value: {
                "van_driest": "Van Driest",
                "volpiani": "Volpiani",
                "compare": "両モデルを比較",
            }[value],
            key="profile_transformation",
        )
        relation = st.selectbox(
            "温度–速度関係",
            tuple(TemperatureVelocityRelation),
            index=list(TemperatureVelocityRelation).index(
                TemperatureVelocityRelation(
                    str(
                        models.get(
                            "temperature_velocity_relation",
                            TemperatureVelocityRelation.GENERALIZED_REYNOLDS_ANALOGY.value,
                        )
                    )
                )
            ),
            format_func=lambda value: (
                "Generalized Reynolds analogy"
                if value is TemperatureVelocityRelation.GENERALIZED_REYNOLDS_ANALOGY
                else "Walz"
            ),
            key="profile_relation",
        )
        assert relation is not None
        adiabatic = st.checkbox(
            "断熱壁（回復温度を壁温に使用）",
            value=inputs.get("wall_temperature") is None,
            key="profile_adiabatic",
        )
        wall_temperature = None
        if not adiabatic:
            wall_temperature = finite_number(
                f"壁温 T_w [{preferences.temperature}]",
                _display(
                    _number(inputs, "wall_temperature", 250.0),
                    "temperature",
                    preferences.temperature,
                ),
                key="profile_wall_temperature",
            )
        automatic_wake = st.checkbox(
            "wake parameter Πを自動決定",
            value=inputs.get("wake_parameter") is None,
            key="profile_automatic_wake",
        )
        wake_parameter = None
        if not automatic_wake:
            wake_parameter = finite_number(
                "wake parameter Π",
                _number(inputs, "wake_parameter", 0.3),
                key="profile_wake",
                min_value=0.0,
            )
        points = int(
            st.number_input(
                "壁面方向点数",
                51,
                501,
                int(inputs.get("points", 257)),
                1,
                key="profile_points",
            )
        )
        submitted = st.form_submit_button("計算", type="primary")

    if submitted:
        st.session_state.pop("profile_payload", None)
        try:
            transformations = (
                tuple(CompressibleVelocityTransformation)
                if comparison == "compare"
                else (CompressibleVelocityTransformation(str(comparison)),)
            )
            velocity_si = _si(edge_velocity, "speed", preferences.speed)
            density_si = _si(edge_density, "density", preferences.density)
            temperature_si = _si(
                edge_temperature, "temperature", preferences.temperature
            )
            thickness_value_si = _si(thickness, "length", preferences.length)
            shear_value_si = _si(shear, "pressure", preferences.pressure)
            wall_si = (
                None
                if wall_temperature is None
                else _si(wall_temperature, "temperature", preferences.temperature)
            )
            calculation = boundary_layer_profiles(
                edge_velocity=velocity_si,
                edge_density=density_si,
                edge_temperature=temperature_si,
                boundary_layer_thickness=thickness_value_si,
                wall_shear_stress=shear_value_si,
                transformations=transformations,
                temperature_velocity_relation=relation,
                wall_temperature=wall_si,
                wake_parameter=wake_parameter,
                points=points,
            )
            configuration = make_configuration(
                calculator="boundary_layer_profile",
                mode="single",
                inputs_si={
                    "edge_velocity": velocity_si,
                    "edge_density": density_si,
                    "edge_temperature": temperature_si,
                    "boundary_layer_thickness": thickness_value_si,
                    "wall_shear_stress": shear_value_si,
                    "wall_temperature": wall_si,
                    "wake_parameter": wake_parameter,
                    "points": points,
                },
                models={
                    "source": source,
                    "transformation": comparison,
                    "temperature_velocity_relation": relation.value,
                },
                units=preferences,
            )
        except ValueError as error:
            st.error(str(error), icon="🚫")
        else:
            st.session_state["profile_payload"] = (calculation, configuration)

    payload = _profile_payload()
    if payload is None:
        if not has_linked:
            st.info("平板境界層画面で乱流単点ケースを保存すると入力を引き継げます。")
        with st.expander("モデルの前提・適用範囲"):
            st.write("滑面、完全乱流、ゼロ圧力勾配、熱量的完全気体AIRを仮定します。")
        return
    calculation, configuration = payload

    def metrics(row: Mapping[str, object]) -> None:
        columns = st.columns(4)
        for column, item in zip(
            columns,
            (
                ("Re_τ", "Re_τ"),
                ("u_τ", "u_τ"),
                ("δ*", "δ*"),
                ("形状係数 H", "形状係数"),
            ),
            strict=True,
        ):
            with column:
                _metric(row, *item)

    render_result_bundle(
        calculator="boundary_layer_profile",
        result=calculation.result,
        configuration=configuration,
        preferences=preferences,
        figures=boundary_layer_profile_figures(calculation.result.rows, preferences),
        filename_prefix="aerophysics-boundary-layer-profile",
        metrics=metrics,
    )
    profile_options = {
        profile.transformation.value: profile for profile in calculation.profiles
    }
    selected = st.selectbox(
        "突起抗力へ保存するプロファイル",
        tuple(profile_options),
        format_func=lambda value: _TRANSFORM_LABELS[
            CompressibleVelocityTransformation(value)
        ],
        key="profile_transfer_model",
    )
    if st.button("現在の境界層プロファイルとして保存", key="profile_save_case"):
        st.session_state["current_boundary_profile"] = profile_options[selected]
        st.success("SIプロファイルを保存しました。")
    with st.expander("モデルの前提・適用範囲"):
        st.write(
            "Re_τ < 500では尺度分離が弱い旨を警告します。粗面、遷移、圧力勾配、"
            "実在気体効果は含みません。"
        )


_SHAPE_LABELS = {
    "rectangle": "矩形",
    "triangle": "三角形",
    "ellipse": "楕円形",
    "csv": "形状CSV",
}


def render_protrusion_drag(preferences: UnitPreferences) -> None:
    """Render boundary-layer-immersed protrusion drag calculations."""
    st.title("突起抗力")
    st.caption("境界層内の実効動圧を投影面積で積分して直接抗力を推定します。")
    imported = pop_pending_configuration("protrusion_drag")
    inputs, models, sweep = _defaults(imported)
    if imported is not None:
        imported_profile = tuple(
            _array(inputs, key)
            for key in ("profile_height", "profile_velocity", "profile_density")
        )
        if all(value is not None for value in imported_profile):
            st.session_state["protrusion_embedded_profile"] = ProfileCSV(
                imported_profile[0],  # type: ignore[arg-type]
                imported_profile[1],  # type: ignore[arg-type]
                imported_profile[2],  # type: ignore[arg-type]
            )
        imported_shape = tuple(
            _array(inputs, key) for key in ("shape_height", "shape_width")
        )
        if all(value is not None for value in imported_shape):
            st.session_state["protrusion_embedded_shape"] = ShapeCSV(
                imported_shape[0],  # type: ignore[arg-type]
                imported_shape[1],  # type: ignore[arg-type]
            )
    render_configuration_import("protrusion_drag", "protrusion")
    render_reset_button("protrusion", "protrusion_payload")
    st.download_button(
        "プロファイルCSVテンプレート",
        profile_csv_template(),
        file_name="aerophysics-profile-template.csv",
        mime="text/csv",
        key="protrusion_profile_template",
    )
    st.download_button(
        "形状CSVテンプレート",
        shape_csv_template(),
        file_name="aerophysics-shape-template.csv",
        mime="text/csv",
        key="protrusion_shape_template",
    )
    saved = st.session_state.get("current_boundary_profile")
    has_saved = isinstance(saved, BoundaryProfileCase)
    embedded_profile_data = st.session_state.get("protrusion_embedded_profile")
    embedded_profile = isinstance(embedded_profile_data, ProfileCSV)
    embedded_shape_data = st.session_state.get("protrusion_embedded_shape")

    with st.form("protrusion_form"):
        sources = ["power_law"]
        if has_saved:
            sources.append("saved")
        sources.append("csv")
        default_source = str(models.get("profile_source", "power_law"))
        if (
            imported is not None
            and embedded_profile
            and default_source
            in {
                "saved",
                "csv",
            }
        ):
            default_source = "csv"
        source = st.radio(
            "境界層プロファイル",
            tuple(sources),
            index=sources.index(default_source) if default_source in sources else 0,
            format_func=lambda value: {
                "power_law": "1/7乗則",
                "saved": "現在の境界層プロファイル",
                "csv": "プロファイルCSV",
            }[value],
            horizontal=True,
            key="protrusion_source",
        )
        uploaded_profile = (
            st.file_uploader(
                "wall_distance, velocity, density CSV",
                type="csv",
                key="protrusion_profile_upload",
            )
            if source == "csv"
            else None
        )
        linked_profile = saved if source == "saved" and has_saved else None
        edge_velocity_si = (
            linked_profile.edge_velocity
            if isinstance(linked_profile, BoundaryProfileCase)
            else _number(inputs, "edge_velocity", 300.0)
        )
        edge_density_si = (
            linked_profile.edge_density
            if isinstance(linked_profile, BoundaryProfileCase)
            else _number(inputs, "edge_density", 1.0)
        )
        thickness_si = (
            linked_profile.boundary_layer_thickness
            if isinstance(linked_profile, BoundaryProfileCase)
            else _number(inputs, "boundary_layer_thickness", 0.05)
        )
        disabled = isinstance(linked_profile, BoundaryProfileCase)
        edge_velocity = finite_number(
            f"外縁速度 U_e [{preferences.speed}]",
            _display(edge_velocity_si, "speed", preferences.speed),
            key="protrusion_velocity",
            min_value=1.0e-12,
            disabled=disabled,
        )
        edge_density = finite_number(
            f"外縁密度 ρ_e [{preferences.density}]",
            _display(edge_density_si, "density", preferences.density),
            key="protrusion_density",
            min_value=1.0e-12,
            disabled=disabled,
        )
        thickness = finite_number(
            f"境界層厚さ δ [{preferences.length}]",
            _display(thickness_si, "length", preferences.length),
            key="protrusion_thickness",
            min_value=1.0e-12,
            disabled=disabled,
        )
        coefficient = finite_number(
            "自由流抗力係数 C_D",
            _number(inputs, "drag_coefficient", 1.0),
            key="protrusion_cd",
            min_value=0.0,
        )
        height = finite_number(
            f"突起高さ h [{preferences.length}]",
            _display(_number(inputs, "height", 0.01), "length", preferences.length),
            key="protrusion_height",
            min_value=1.0e-12,
        )
        shape = st.selectbox(
            "投影形状",
            tuple(_SHAPE_LABELS),
            index=tuple(_SHAPE_LABELS).index(str(models.get("shape", "rectangle"))),
            format_func=_SHAPE_LABELS.__getitem__,
            key="protrusion_shape",
        )
        width = finite_number(
            f"代表幅 b₀ [{preferences.length}]",
            _display(
                _number(inputs, "base_width", 0.005), "length", preferences.length
            ),
            key="protrusion_width",
            min_value=1.0e-12,
            disabled=shape == "csv",
        )
        uploaded_shape = (
            st.file_uploader(
                "height, width CSV",
                type="csv",
                key="protrusion_shape_upload",
            )
            if shape == "csv"
            else None
        )
        compressible = False
        mach = None
        edge_temperature = None
        wall_temperature = None
        if source == "power_law":
            compressible = st.checkbox(
                "Walz圧縮性近似を適用",
                value=bool(models.get("compressible", False)),
                key="protrusion_compressible",
            )
            if compressible:
                mach = finite_number(
                    "外縁 Mach M_e",
                    _number(inputs, "mach", 2.0),
                    key="protrusion_mach",
                    min_value=0.0,
                )
                edge_temperature = finite_number(
                    f"外縁温度 T_e [{preferences.temperature}]",
                    _display(
                        _number(inputs, "edge_temperature", 250.0),
                        "temperature",
                        preferences.temperature,
                    ),
                    key="protrusion_temperature",
                )
                adiabatic = st.checkbox(
                    "断熱壁",
                    value=inputs.get("wall_temperature") is None,
                    key="protrusion_adiabatic",
                )
                if not adiabatic:
                    wall_temperature = finite_number(
                        f"壁温 T_w [{preferences.temperature}]",
                        _display(
                            _number(inputs, "wall_temperature", 300.0),
                            "temperature",
                            preferences.temperature,
                        ),
                        key="protrusion_wall_temperature",
                    )
        mode = st.radio(
            "計算モード",
            ("single", "sweep"),
            index=0 if (imported or {}).get("mode", "single") == "single" else 1,
            format_func=lambda value: "単点" if value == "single" else "1変数スイープ",
            horizontal=True,
            key="protrusion_mode",
        )
        sweep_field = "height"
        start = stop = 0.0
        points = 101
        if mode == "sweep":
            fields = [
                "height",
                "drag_coefficient",
                "boundary_layer_thickness",
            ]
            if shape != "csv":
                fields.append("base_width")
            if compressible:
                fields.append("mach")
            sweep_field = st.selectbox(
                "スイープ変数",
                tuple(fields),
                index=fields.index(str(sweep.get("field", "height")))
                if str(sweep.get("field", "height")) in fields
                else 0,
                format_func=lambda value: {
                    "height": "突起高さ h",
                    "drag_coefficient": "抗力係数 C_D",
                    "base_width": "代表幅 b₀",
                    "boundary_layer_thickness": "境界層厚さ δ",
                    "mach": "外縁 Mach M_e",
                }[value],
                key="protrusion_sweep_field",
            )
            length_field = sweep_field in {
                "height",
                "base_width",
                "boundary_layer_thickness",
            }
            default_start = _number(sweep, "start", 0.002 if length_field else 0.1)
            default_stop = _number(sweep, "stop", 0.08 if length_field else 2.0)
            if length_field:
                default_start = _display(default_start, "length", preferences.length)
                default_stop = _display(default_stop, "length", preferences.length)
            left, right, count = st.columns(3)
            with left:
                start = finite_number(
                    f"開始 [{preferences.length if length_field else '–'}]",
                    default_start,
                    key="protrusion_sweep_start",
                    min_value=0.0,
                )
            with right:
                stop = finite_number(
                    f"終了 [{preferences.length if length_field else '–'}]",
                    default_stop,
                    key="protrusion_sweep_stop",
                    min_value=0.0,
                )
            with count:
                points = int(
                    st.number_input(
                        "点数",
                        2,
                        501,
                        int(sweep.get("points", 101)),
                        1,
                        key="protrusion_sweep_points",
                    )
                )
        submitted = st.form_submit_button("計算", type="primary")

    if submitted:
        st.session_state.pop("protrusion_payload", None)
        try:
            profile_data: ProfileCSV | None = None
            if isinstance(linked_profile, BoundaryProfileCase):
                profile_data = ProfileCSV(
                    linked_profile.wall_distance,
                    linked_profile.velocity,
                    linked_profile.density,
                )
            elif source == "csv":
                if uploaded_profile is not None:
                    profile_data = parse_profile_csv(
                        uploaded_profile.getvalue(), preferences
                    )
                elif isinstance(embedded_profile_data, ProfileCSV):
                    profile_data = embedded_profile_data
                else:
                    raise ValueError("プロファイルCSVを選択してください")
            shape_data: ShapeCSV | None = None
            if shape == "csv":
                if uploaded_shape is not None:
                    shape_data = parse_shape_csv(uploaded_shape.getvalue(), preferences)
                elif isinstance(embedded_shape_data, ShapeCSV):
                    shape_data = embedded_shape_data
                else:
                    raise ValueError("形状CSVを選択してください")
            velocity_si = _si(edge_velocity, "speed", preferences.speed)
            density_si = _si(edge_density, "density", preferences.density)
            thickness_value_si = _si(thickness, "length", preferences.length)
            height_si = _si(height, "length", preferences.length)
            width_si = _si(width, "length", preferences.length)
            edge_temperature_si = (
                None
                if edge_temperature is None
                else _si(edge_temperature, "temperature", preferences.temperature)
            )
            wall_temperature_si = (
                None
                if wall_temperature is None
                else _si(wall_temperature, "temperature", preferences.temperature)
            )
            profile_height = (
                None if profile_data is None else profile_data.wall_distance
            )
            profile_velocity = None if profile_data is None else profile_data.velocity
            profile_density = None if profile_data is None else profile_data.density
            shape_height = None if shape_data is None else shape_data.height
            shape_width = None if shape_data is None else shape_data.width
            sweep_configuration = None
            if mode == "single":
                result = protrusion_condition(
                    drag_coefficient=coefficient,
                    height=height_si,
                    base_width=width_si,
                    shape=shape,
                    edge_velocity=velocity_si,
                    edge_density=density_si,
                    boundary_layer_thickness=thickness_value_si,
                    profile_height=profile_height,
                    profile_velocity=profile_velocity,
                    profile_density=profile_density,
                    mach=mach,
                    edge_temperature=edge_temperature_si,
                    wall_temperature=wall_temperature_si,
                    shape_height=shape_height,
                    shape_width=shape_width,
                )
            else:
                length_field = sweep_field in {
                    "height",
                    "base_width",
                    "boundary_layer_thickness",
                }
                start_si = (
                    _si(start, "length", preferences.length) if length_field else start
                )
                stop_si = (
                    _si(stop, "length", preferences.length) if length_field else stop
                )
                result = protrusion_sweep(
                    sweep_field=sweep_field,
                    start=start_si,
                    stop=stop_si,
                    points=points,
                    drag_coefficient=coefficient,
                    height=height_si,
                    base_width=width_si,
                    shape=shape,
                    edge_velocity=velocity_si,
                    edge_density=density_si,
                    boundary_layer_thickness=thickness_value_si,
                    profile_height=profile_height,
                    profile_velocity=profile_velocity,
                    profile_density=profile_density,
                    mach=mach,
                    edge_temperature=edge_temperature_si,
                    wall_temperature=wall_temperature_si,
                    shape_height=shape_height,
                    shape_width=shape_width,
                )
                sweep_configuration = {
                    "field": sweep_field,
                    "start": start_si,
                    "stop": stop_si,
                    "points": points,
                }
            configuration = make_configuration(
                calculator="protrusion_drag",
                mode=mode,
                inputs_si={
                    "drag_coefficient": coefficient,
                    "height": height_si,
                    "base_width": width_si,
                    "edge_velocity": velocity_si,
                    "edge_density": density_si,
                    "boundary_layer_thickness": thickness_value_si,
                    "mach": mach,
                    "edge_temperature": edge_temperature_si,
                    "wall_temperature": wall_temperature_si,
                    "profile_height": None
                    if profile_data is None
                    else profile_data.wall_distance.tolist(),
                    "profile_velocity": None
                    if profile_data is None
                    else profile_data.velocity.tolist(),
                    "profile_density": None
                    if profile_data is None
                    else profile_data.density.tolist(),
                    "shape_height": None
                    if shape_data is None
                    else shape_data.height.tolist(),
                    "shape_width": None
                    if shape_data is None
                    else shape_data.width.tolist(),
                },
                models={
                    "profile_source": source,
                    "shape": shape,
                    "compressible": compressible,
                },
                units=preferences,
                sweep_si=sweep_configuration,
            )
        except ValueError as error:
            st.error(str(error), icon="🚫")
        else:
            figure_settings = {
                "height": height_si,
                "base_width": width_si,
                "boundary_layer_thickness": thickness_value_si,
                "shape": shape,
                "shape_height": None if shape_data is None else shape_data.height,
                "shape_width": None if shape_data is None else shape_data.width,
            }
            st.session_state["protrusion_payload"] = (
                result,
                configuration,
                figure_settings,
            )

    raw_payload = st.session_state.get("protrusion_payload")
    if not (
        isinstance(raw_payload, tuple)
        and len(raw_payload) == 3
        and isinstance(raw_payload[0], CalculationResult)
        and isinstance(raw_payload[1], dict)
        and isinstance(raw_payload[2], dict)
    ):
        with st.expander("モデルの前提・適用範囲"):
            st.write("孤立突起の直接抗力のみを評価し、馬蹄渦や壁面干渉は含みません。")
        return
    result, configuration, figure_settings = raw_payload
    config_sweep = configuration.get("sweep_si")
    figure_field = (
        str(config_sweep.get("field", "height"))
        if isinstance(config_sweep, dict)
        else "height"
    )
    figures = (
        protrusion_figures(result.rows, preferences, sweep_field=figure_field)
        if configuration.get("mode") == "sweep"
        else {
            "投影形状": protrusion_shape_figure(
                preferences=preferences,
                **figure_settings,
            )
        }
    )

    def metrics(row: Mapping[str, object]) -> None:
        columns = st.columns(4)
        for column, item in zip(
            columns,
            (
                ("直接抗力 D", "直接抗力"),
                ("実効動圧", "実効動圧"),
                ("遮蔽係数", "遮蔽係数"),
                ("h/δ", "h/δ"),
            ),
            strict=True,
        ):
            with column:
                _metric(row, *item)

    render_result_bundle(
        calculator="protrusion_drag",
        result=result,
        configuration=configuration,
        preferences=preferences,
        figures=figures,
        filename_prefix="aerophysics-protrusion-drag",
        metrics=metrics,
    )
    invalid = sum(row.get("status") != "ok" for row in result.rows)
    if invalid:
        st.warning(f"{invalid}点は入力条件を満たさないため欠損値としました。")
    with st.expander("モデルの前提・適用範囲"):
        st.write(
            "自由流抗力係数を局所実効動圧へ適用する工学推定です。壁面干渉、"
            "馬蹄渦、遷移、複数突起、衝撃波干渉は含みません。"
        )


def render_thermochemistry(preferences: UnitPreferences) -> None:
    """Render frozen-composition dry-air thermochemistry."""
    st.title("熱化学")
    st.caption("NASA7／NASA9多項式による凍結組成乾燥空気の温度依存物性です。")
    imported = pop_pending_configuration("thermochemistry")
    inputs, models, sweep = _defaults(imported)
    render_configuration_import("thermochemistry", "thermo")
    render_reset_button("thermo", "thermo_payload")

    with st.form("thermo_form"):
        mode = st.radio(
            "計算モード",
            ("single", "sweep"),
            index=0 if (imported or {}).get("mode", "sweep") == "single" else 1,
            format_func=lambda value: "単点" if value == "single" else "温度スイープ",
            horizontal=True,
            key="thermo_mode",
        )
        selection = st.selectbox(
            "NASA多項式",
            ("NASA7", "NASA9", "compare"),
            index=("NASA7", "NASA9", "compare").index(
                str(models.get("selection", "compare"))
            ),
            format_func=lambda value: (
                "NASA7とNASA9を比較" if value == "compare" else value
            ),
            key="thermo_selection",
        )
        temperature = finite_number(
            f"温度 T [{preferences.temperature}]",
            _display(
                _number(inputs, "temperature", 300.0),
                "temperature",
                preferences.temperature,
            ),
            key="thermo_temperature",
        )
        pressure = finite_number(
            f"圧力 p [{preferences.pressure}]",
            _display(
                _number(inputs, "pressure", 101_325.0),
                "pressure",
                preferences.pressure,
            ),
            key="thermo_pressure",
            min_value=1.0e-12,
        )
        reference = finite_number(
            f"顕熱基準温度 T_ref [{preferences.temperature}]",
            _display(
                _number(inputs, "reference_temperature", 298.15),
                "temperature",
                preferences.temperature,
            ),
            key="thermo_reference",
        )
        allow_extrapolation = st.checkbox(
            "200–6000 K外への外挿を明示的に許可",
            value=bool(models.get("allow_extrapolation", False)),
            key="thermo_extrapolate",
        )
        start = stop = 0.0
        points = 201
        if mode == "sweep":
            left, right, count = st.columns(3)
            with left:
                start = finite_number(
                    f"開始温度 [{preferences.temperature}]",
                    _display(
                        _number(sweep, "start", 200.0),
                        "temperature",
                        preferences.temperature,
                    ),
                    key="thermo_sweep_start",
                )
            with right:
                stop = finite_number(
                    f"終了温度 [{preferences.temperature}]",
                    _display(
                        _number(sweep, "stop", 6000.0),
                        "temperature",
                        preferences.temperature,
                    ),
                    key="thermo_sweep_stop",
                )
            with count:
                points = int(
                    st.number_input(
                        "点数",
                        2,
                        501,
                        int(sweep.get("points", 201)),
                        1,
                        key="thermo_sweep_points",
                    )
                )
        submitted = st.form_submit_button("計算", type="primary")

    if submitted:
        st.session_state.pop("thermo_payload", None)
        try:
            temperature_si = _si(temperature, "temperature", preferences.temperature)
            pressure_si = _si(pressure, "pressure", preferences.pressure)
            reference_si = _si(reference, "temperature", preferences.temperature)
            selected_models = (
                ("NASA7", "NASA9") if selection == "compare" else (selection,)
            )
            sweep_configuration = None
            if mode == "single":
                result = thermochemistry_condition(
                    temperature=temperature_si,
                    pressure=pressure_si,
                    reference_temperature=reference_si,
                    models=selected_models,
                    allow_extrapolation=allow_extrapolation,
                )
            else:
                start_si = _si(start, "temperature", preferences.temperature)
                stop_si = _si(stop, "temperature", preferences.temperature)
                result = thermochemistry_sweep(
                    start=start_si,
                    stop=stop_si,
                    points=points,
                    pressure=pressure_si,
                    reference_temperature=reference_si,
                    models=selected_models,
                    allow_extrapolation=allow_extrapolation,
                )
                sweep_configuration = {
                    "field": "temperature",
                    "start": start_si,
                    "stop": stop_si,
                    "points": points,
                }
            configuration = make_configuration(
                calculator="thermochemistry",
                mode=mode,
                inputs_si={
                    "temperature": temperature_si,
                    "pressure": pressure_si,
                    "reference_temperature": reference_si,
                },
                models={
                    "selection": selection,
                    "allow_extrapolation": allow_extrapolation,
                },
                units=preferences,
                sweep_si=sweep_configuration,
            )
        except ValueError as error:
            st.error(str(error), icon="🚫")
        else:
            st.session_state["thermo_payload"] = (result, configuration)

    payload = _result_payload("thermo_payload")
    if payload is None:
        with st.expander("モデルの前提・適用範囲"):
            st.write("200–6000 Kの凍結組成理想気体で、解離・反応・電離は含みません。")
        return
    result, configuration = payload

    def metrics(row: Mapping[str, object]) -> None:
        columns = st.columns(4)
        for column, item in zip(
            columns,
            (
                ("c_p", "c_p"),
                ("γ", "比熱比"),
                ("音速 a", "音速"),
                ("顕熱 Δh", "顕熱エンタルピー"),
            ),
            strict=True,
        ):
            with column:
                _metric(row, *item)

    render_result_bundle(
        calculator="thermochemistry",
        result=result,
        configuration=configuration,
        preferences=preferences,
        figures=thermochemistry_figures(result.rows, preferences),
        filename_prefix="aerophysics-thermochemistry",
        metrics=metrics,
    )
    with st.expander("モデルの前提・適用範囲"):
        st.write(
            "N₂/O₂/Ar/CO₂のモル分率を固定した理想混合気体です。NASA標準"
            "エンタルピーと基準温度差の顕熱は異なる量として表示します。"
        )


def render_viscosity(preferences: UnitPreferences) -> None:
    """Render dry-air dynamic-viscosity model calculations."""
    st.title("粘性係数")
    st.caption("Sutherland／Keyes／Blottner-Wilkeによる乾燥空気の動的粘性係数です。")
    imported = pop_pending_configuration("viscosity")
    inputs, models, sweep = _defaults(imported)
    render_configuration_import("viscosity", "viscosity")
    render_reset_button("viscosity", "viscosity_payload")

    selections = ("Sutherland", "Keyes", "Blottner/Wilke", "compare")
    selected_default = str(models.get("selection", "compare"))
    if selected_default not in selections:
        selected_default = "compare"
    scale_default = str(sweep.get("scale", "log"))
    if scale_default not in {"linear", "log"}:
        scale_default = "log"

    with st.form("viscosity_form"):
        mode = st.radio(
            "計算モード",
            ("single", "sweep"),
            index=0 if (imported or {}).get("mode", "sweep") == "single" else 1,
            format_func=lambda value: "単点" if value == "single" else "温度スイープ",
            horizontal=True,
            key="viscosity_mode",
        )
        selection = st.selectbox(
            "粘性モデル",
            selections,
            index=selections.index(selected_default),
            format_func=lambda value: (
                "3モデルを比較" if value == "compare" else value
            ),
            key="viscosity_selection",
        )
        temperature = finite_number(
            f"温度 T [{preferences.temperature}]",
            _display(
                _number(inputs, "temperature", 1000.0),
                "temperature",
                preferences.temperature,
            ),
            key="viscosity_temperature",
        )
        allow_extrapolation = st.checkbox(
            "公称範囲外への外挿を許可",
            value=bool(models.get("allow_extrapolation", False)),
            key="viscosity_extrapolate",
        )
        start = stop = 0.0
        points = 201
        scale = "log"
        if mode == "sweep":
            left, right, count = st.columns(3)
            with left:
                start = finite_number(
                    f"開始温度 [{preferences.temperature}]",
                    _display(
                        _number(sweep, "start", 79.0),
                        "temperature",
                        preferences.temperature,
                    ),
                    key="viscosity_sweep_start",
                )
            with right:
                stop = finite_number(
                    f"終了温度 [{preferences.temperature}]",
                    _display(
                        _number(sweep, "stop", 30_000.0),
                        "temperature",
                        preferences.temperature,
                    ),
                    key="viscosity_sweep_stop",
                )
            with count:
                points = int(
                    st.number_input(
                        "点数",
                        2,
                        501,
                        int(sweep.get("points", 201)),
                        1,
                        key="viscosity_sweep_points",
                    )
                )
            scale = st.radio(
                "温度点の配置",
                ("log", "linear"),
                index=0 if scale_default == "log" else 1,
                format_func=lambda value: "対数" if value == "log" else "線形",
                horizontal=True,
                key="viscosity_sweep_scale",
            )
        submitted = st.form_submit_button("計算", type="primary")

    if submitted:
        st.session_state.pop("viscosity_payload", None)
        try:
            temperature_si = _si(temperature, "temperature", preferences.temperature)
            selected_models = (
                ("Sutherland", "Keyes", "Blottner/Wilke")
                if selection == "compare"
                else (selection,)
            )
            sweep_configuration = None
            if mode == "single":
                result = viscosity_condition(
                    temperature=temperature_si,
                    models=selected_models,
                    allow_extrapolation=allow_extrapolation,
                )
            else:
                start_si = _si(start, "temperature", preferences.temperature)
                stop_si = _si(stop, "temperature", preferences.temperature)
                result = viscosity_sweep(
                    start=start_si,
                    stop=stop_si,
                    points=points,
                    models=selected_models,
                    allow_extrapolation=allow_extrapolation,
                    log_temperature=scale == "log",
                )
                sweep_configuration = {
                    "field": "temperature",
                    "start": start_si,
                    "stop": stop_si,
                    "points": points,
                    "scale": scale,
                }
            configuration = make_configuration(
                calculator="viscosity",
                mode=mode,
                inputs_si={"temperature": temperature_si},
                models={
                    "selection": selection,
                    "allow_extrapolation": allow_extrapolation,
                },
                units=preferences,
                sweep_si=sweep_configuration,
            )
        except ValueError as error:
            st.error(str(error), icon="🚫")
        else:
            st.session_state["viscosity_payload"] = (result, configuration)

    payload = _result_payload("viscosity_payload")
    if payload is None:
        with st.expander("モデルの前提・適用範囲"):
            st.write(
                "Keyesは79–1845 K、Blottner/Wilkeは1000–30000 Kが公称範囲です。"
                "Sutherlandは標準大気と飛行条件の既定モデルです。"
            )
        return
    result, configuration = payload

    def metrics(_: Mapping[str, object]) -> None:
        columns = st.columns(len(result.rows))
        for column, row in zip(columns, result.rows, strict=True):
            viscosity = row.get("dynamic_viscosity")
            difference = row.get("relative_difference")
            value = f"{viscosity:.7g} Pa·s" if isinstance(viscosity, float) else "—"
            delta = (
                f"{difference:+.3f}% vs Sutherland"
                if isinstance(difference, float)
                else None
            )
            with column:
                st.metric(str(row.get("model", "")), value, delta=delta)

    figures: dict[str, Any] = {}
    configured_sweep = configuration.get("sweep_si")
    if configuration["mode"] == "sweep" and isinstance(configured_sweep, dict):
        figures = viscosity_figures(
            result.rows,
            preferences,
            log_temperature=configured_sweep.get("scale") == "log",
        )

    render_result_bundle(
        calculator="viscosity",
        result=result,
        configuration=configuration,
        preferences=preferences,
        figures=figures,
        filename_prefix="aerophysics-viscosity",
        metrics=metrics if configuration["mode"] == "single" else None,
    )
    with st.expander("モデルの前提・適用範囲"):
        st.write(
            "Keyesは79–1845 K、Blottner/Wilkeは1000–30000 Kが公称範囲です。"
            "Sutherlandは標準大気と飛行条件の既定モデルです。Blottner/Wilkeは"
            "N₂/O₂/Ar/CO₂固定組成で、解離・反応・電離を含みません。対数温度軸は"
            "表示温度単位にかかわらず絶対温度Kで表示します。"
        )
