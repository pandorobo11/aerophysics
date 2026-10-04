"""Interactive editor for local flow-state sequences."""

import copy
import csv
import io
from typing import Any, cast

import streamlit as st

from aerophysics.gui.adapters import _SHOCK_GASES
from aerophysics.gui.components import (
    pop_pending_configuration,
    render_configuration_import,
)
from aerophysics.gui.config import (
    dump_configuration,
    make_configuration,
)
from aerophysics.gui.sequence_adapter import (
    BASES,
    INITIAL_FIELDS,
    KINDS,
    calculate_sequence,
    sequence_rows,
)
from aerophysics.gui.units import QuantityKind, UnitPreferences, from_si, to_si

_LABELS = {
    "normal": "垂直衝撃波",
    "oblique": "斜め衝撃波",
    "conical": "円錐衝撃波",
    "expansion": "膨張波",
    "isentropic": "等エントロピー変化",
}
_INITIAL_LABELS = {
    "static_mach": "静圧・静温・Mach",
    "static_velocity": "静圧・静温・速度",
    "total": "全圧・全温・Mach",
    "atmosphere_mach": "標準大気・Mach",
    "atmosphere_velocity": "標準大気・速度",
}
_QUANTITIES: dict[str, QuantityKind] = {
    "pressure": "pressure",
    "temperature": "temperature",
    "altitude": "length",
    "velocity": "speed",
    "angle": "angle",
}
_DEFAULTS = {
    "pressure": 101325.0,
    "temperature": 288.15,
    "mach": 3.0,
    "altitude": 10000.0,
    "velocity": 900.0,
}


def _input(field: str, value: float, key: str, units: UnitPreferences) -> float:
    if field not in _QUANTITIES:
        return float(st.number_input(field, value=float(value), key=key))
    kind = _QUANTITIES[field]
    unit = getattr(units, kind)
    shown = float(from_si(value, kind, unit))
    entered = st.number_input(
        f"{field} [{unit}]", value=shown, key=f"{key}_{unit}", format="%.8g"
    )
    return float(to_si(entered, kind, unit))


def _default_step(kind: str) -> dict[str, Any]:
    step: dict[str, Any] = {"kind": kind, "name": _LABELS[kind]}
    if kind in {"oblique", "conical", "expansion"}:
        step["angle"] = 0.08726646259971647
    if kind == "oblique":
        step["branch"] = "weak"
    if kind == "isentropic":
        step.update(basis="mach", value=2.0, branch=None)
    return step


def _edit_steps(units: UnitPreferences, epoch: int) -> list[dict[str, Any]]:
    steps: list[dict[str, Any]] = st.session_state["sequence_steps"]
    action: tuple[str, int] | None = None
    for index, step in enumerate(steps):
        prefix = f"sequence_{epoch}_{index}"
        with st.container(border=True):
            st.subheader(f"{index + 1}. {_LABELS[step['kind']]}")
            step["name"] = st.text_input(
                "名前", value=step["name"], key=f"{prefix}_name"
            )
            if "angle" in step:
                step["angle"] = _input("angle", step["angle"], f"{prefix}_angle", units)
            if step["kind"] == "oblique":
                step["branch"] = st.selectbox(
                    "衝撃波の枝",
                    ("weak", "strong"),
                    index=("weak", "strong").index(step["branch"]),
                    key=f"{prefix}_branch",
                )
            if step["kind"] == "isentropic":
                bases = ("mach", "pressure_ratio", "area_ratio")
                step["basis"] = st.selectbox(
                    "到達条件",
                    bases,
                    index=bases.index(step["basis"]),
                    key=f"{prefix}_basis",
                )
                st.caption("pressure_ratio = p₂/p₁、area_ratio = A₂/A₁")
                step["value"] = st.number_input(
                    "到達値", value=float(step["value"]), key=f"{prefix}_value"
                )
                step["branch"] = (
                    st.selectbox(
                        "出口の枝",
                        ("supersonic", "subsonic"),
                        index=0 if step["branch"] != "subsonic" else 1,
                        key=f"{prefix}_exit",
                    )
                    if step["basis"] == "area_ratio"
                    else None
                )
            for column, label, operation in zip(
                st.columns(4),
                ("↑", "↓", "複製", "削除"),
                ("up", "down", "copy", "delete"),
                strict=True,
            ):
                if column.button(
                    label,
                    key=f"{prefix}_{operation}",
                    disabled=(operation == "up" and index == 0)
                    or (operation == "down" and index == len(steps) - 1),
                ):
                    action = operation, index
    if action is not None:
        operation, index = action
        if operation == "delete":
            steps.pop(index)
        elif operation == "copy":
            steps.insert(index + 1, copy.deepcopy(steps[index]))
        else:
            other = index + (-1 if operation == "up" else 1)
            steps[index], steps[other] = steps[other], steps[index]
        st.session_state["sequence_epoch"] += 1
        st.rerun()
    return steps


def _display_rows(
    rows: list[dict[str, Any]], units: UnitPreferences
) -> list[dict[str, Any]]:
    kinds: dict[str, QuantityKind] = {
        "pressure": "pressure",
        "total_pressure": "pressure",
        "dynamic_pressure": "pressure",
        "temperature": "temperature",
        "total_temperature": "temperature",
        "density": "density",
        "velocity": "speed",
        "speed_of_sound": "speed",
        "shock_angle": "angle",
    }
    result = []
    for row in rows:
        converted = {}
        for field, value in row.items():
            if field in kinds:
                kind = kinds[field]
                unit = getattr(units, kind)
                converted[f"{field} [{unit}]"] = (
                    None if value is None else float(from_si(value, kind, unit))
                )
            else:
                suffix = {
                    "mass_flux": " [kg/(m² s)]",
                    "total_enthalpy": " [J/kg]",
                    "entropy_change": " [J/(kg K)]",
                }.get(field, "")
                converted[field + suffix] = value
        result.append(converted)
    return result


def render_sequence(preferences: UnitPreferences) -> None:
    """Render inputs without running numerical solvers until Calculate."""
    st.title("連続流れ計算")
    st.caption(
        "各段は局所状態の変換です。円錐段は円錐表面状態を引き継ぎます。流路形状・波の干渉・始動性は計算しません。"
    )
    st.session_state.setdefault("sequence_steps", [])
    st.session_state.setdefault("sequence_epoch", 0)
    imported = pop_pending_configuration("flow_sequence")
    if imported is not None:
        st.session_state["sequence_draft"] = imported
        inputs = cast(dict[str, Any], imported["inputs_si"])
        st.session_state["sequence_steps"] = copy.deepcopy(inputs["steps"])
        st.session_state["sequence_epoch"] += 1
    render_configuration_import("flow_sequence", "sequence")
    previous_units = st.session_state.get("sequence_units")
    if previous_units is not None and previous_units != preferences.to_dict():
        st.session_state["sequence_epoch"] += 1
    st.session_state["sequence_units"] = preferences.to_dict()
    saved = st.session_state.get("sequence_draft", {})
    initial_saved = saved.get("inputs_si", {}).get("initial", {"basis": "static_mach"})
    epoch = st.session_state["sequence_epoch"]
    gas_names = tuple(_SHOCK_GASES)
    gas = st.selectbox(
        "気体モデル",
        gas_names,
        index=gas_names.index(saved.get("models", {}).get("gas_model", "AIR")),
        key=f"sequence_gas_{epoch}",
    )
    basis = st.selectbox(
        "初期状態",
        BASES,
        index=BASES.index(initial_saved["basis"]),
        format_func=lambda item: _INITIAL_LABELS[item],
        key=f"sequence_initial_basis_{epoch}",
    )
    if basis == "total":
        st.caption("pressure と temperature は全圧・全温です。")
    if basis.startswith("atmosphere"):
        st.caption(
            "US Standard Atmosphere 1976、幾何高度 −5～86 km。"
            "p,T を取得し、他の量は選択した気体で算出します。"
        )
    initial: dict[str, Any] = {"basis": basis}
    if basis.startswith("atmosphere"):
        coordinates = ("geometric", "geopotential")
        initial["altitude_basis"] = st.radio(
            "高度の種類",
            coordinates,
            index=coordinates.index(initial_saved.get("altitude_basis", "geometric")),
            format_func=lambda value: (
                "幾何高度" if value == "geometric" else "ジオポテンシャル高度"
            ),
            horizontal=True,
            key=f"sequence_altitude_basis_{epoch}",
        )
        st.caption("altitude の入力値を、選択した高度の種類で解釈します。")
    for field in INITIAL_FIELDS[basis]:
        initial[field] = _input(
            field,
            initial_saved.get(field, _DEFAULTS[field]),
            f"sequence_initial_{epoch}_{basis}_{field}",
            preferences,
        )
    st.session_state["sequence_draft"] = {
        "inputs_si": {"initial": initial},
        "models": {"gas_model": gas},
    }
    kind = st.selectbox("追加する段", KINDS, format_func=lambda item: _LABELS[item])
    if st.button("段を追加"):
        st.session_state["sequence_steps"].append(_default_step(kind))
        st.rerun()
    steps = _edit_steps(preferences, epoch)
    configuration = make_configuration(
        calculator="flow_sequence",
        mode="single",
        inputs_si={"initial": initial, "steps": copy.deepcopy(steps)},
        models={"gas_model": gas},
        units=preferences,
    )
    encoded = dump_configuration(configuration)
    st.download_button(
        "設定JSONを保存", encoded, "aerophysics-flow-sequence.json", "application/json"
    )
    if st.button("計算", type="primary"):
        st.session_state.pop("sequence_result", None)
        try:
            result = calculate_sequence(
                cast(dict[str, Any], configuration["inputs_si"]),
                cast(dict[str, Any], configuration["models"]),
            )
            st.session_state["sequence_result"] = (encoded, result)
        except ValueError as error:
            st.error(str(error))
    payload = st.session_state.get("sequence_result")
    if payload is None:
        return
    previous, result = payload
    if previous != encoded:
        st.warning(
            "入力が変更されています。以下は前回計算の結果です。再計算してください。"
        )
    if not result.complete:
        for stage in result.steps:
            if stage.message:
                st.error(f"{stage.step.name}: {stage.status}: {stage.message}")
    rows = sequence_rows(result)
    for row in rows:
        if row.get("warning"):
            st.warning(f"{row['name']}: {row['warning']}")
    st.dataframe(_display_rows(rows, preferences), hide_index=True)
    metric = st.selectbox(
        "推移グラフ", ("mach", "pressure", "temperature", "total_pressure_recovery")
    )
    st.line_chart(
        [{"stage": row["stage"], metric: row.get(metric)} for row in rows],
        x="stage",
        y=metric,
    )
    st.caption(
        "グラフ・CSVの次元量はSI単位です。entropy_change は初期状態からの累積値です。"
    )
    output = io.StringIO()
    writer = csv.DictWriter(
        output, fieldnames=list(dict.fromkeys(key for row in rows for key in row))
    )
    writer.writeheader()
    writer.writerows(rows)
    st.download_button(
        "結果CSVを保存（SI）",
        output.getvalue(),
        "aerophysics-flow-sequence.csv",
        "text/csv",
    )
