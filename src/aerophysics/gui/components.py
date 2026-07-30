"""Shared Streamlit components."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any

import streamlit as st

from aerophysics.gui.adapters import CalculationResult
from aerophysics.gui.config import (
    ConfigurationError,
    dump_configuration,
    load_configuration,
)
from aerophysics.gui.tables import display_rows, rows_to_csv
from aerophysics.gui.units import UnitPreferences

PLOTLY_CONFIG = {
    "displaylogo": False,
    "scrollZoom": False,
    "toImageButtonOptions": {"format": "png", "scale": 2},
}


def render_unit_sidebar() -> UnitPreferences:
    """Render persistent global unit selectors."""
    st.sidebar.header("表示単位")
    with st.sidebar.expander("単位設定", expanded=False):
        length = st.selectbox("長さ", ("m", "ft"), key="unit_length")
        speed = st.selectbox("速度", ("m/s", "kt"), key="unit_speed")
        pressure = st.selectbox("圧力", ("Pa", "psi"), key="unit_pressure")
        temperature = st.selectbox("温度", ("K", "°F"), key="unit_temperature")
        density = st.selectbox("密度", ("kg/m³", "slug/ft³"), key="unit_density")
        angle = st.selectbox("角度", ("deg", "rad"), key="unit_angle")
    st.sidebar.caption("計算コアへ渡す前にSIへ明示変換します。")
    return UnitPreferences(
        length=length,
        speed=speed,
        pressure=pressure,
        temperature=temperature,
        density=density,
        angle=angle,
    )


def _set_unit_state(preferences: UnitPreferences) -> None:
    for name, value in preferences.to_dict().items():
        st.session_state[f"unit_{name}"] = value


def render_configuration_import(calculator: str, prefix: str) -> None:
    """Validate an uploaded JSON file and queue it for a clean rerun."""
    with st.expander("設定JSONを読み込む", expanded=False):
        uploaded = st.file_uploader(
            "以前に保存した設定JSON",
            type="json",
            key=f"{prefix}_configuration_upload",
        )
        if st.button("設定を適用", key=f"{prefix}_configuration_apply"):
            if uploaded is None:
                st.error("JSONファイルを選択してください。")
                return
            try:
                configuration = load_configuration(uploaded.getvalue())
                if configuration["calculator"] != calculator:
                    raise ConfigurationError(
                        "この画面とは異なる計算設定が指定されています。"
                    )
                units = UnitPreferences.from_dict(configuration["display_units"])
            except (ConfigurationError, ValueError) as error:
                st.error(str(error))
                return
            for key in list(st.session_state):
                if isinstance(key, str) and key.startswith(f"{prefix}_"):
                    del st.session_state[key]
            _set_unit_state(units)
            st.session_state[f"pending_{calculator}_configuration"] = configuration
            st.rerun()


def pop_pending_configuration(calculator: str) -> dict[str, object] | None:
    """Pop a configuration queued by the import component."""
    value = st.session_state.pop(f"pending_{calculator}_configuration", None)
    return value if isinstance(value, dict) else None


def render_reset_button(prefix: str, payload_key: str) -> None:
    """Reset one page to its documented representative defaults."""
    if st.button("代表値に戻す", key=f"{prefix}_reset"):
        for key in list(st.session_state):
            if isinstance(key, str) and key.startswith(f"{prefix}_"):
                del st.session_state[key]
        st.session_state.pop(payload_key, None)
        st.rerun()


def render_warnings(result: CalculationResult) -> None:
    """Render captured applicability warnings."""
    for message in result.warnings:
        st.warning(message, icon="⚠️")


def render_result_bundle(
    *,
    calculator: str,
    result: CalculationResult,
    configuration: dict[str, object],
    preferences: UnitPreferences,
    figures: dict[str, Any],
    filename_prefix: str,
    metrics: Callable[[Mapping[str, object]], None] | None = None,
) -> None:
    """Render metrics, figures, table, and reproducibility downloads."""
    render_warnings(result)
    rows = display_rows(calculator, result.rows, preferences)
    if metrics is not None and rows:
        metrics(rows[0])
    if figures:
        tabs = st.tabs(list(figures))
        for tab, figure in zip(tabs, figures.values(), strict=True):
            with tab:
                st.plotly_chart(
                    figure,
                    width="stretch",
                    config=PLOTLY_CONFIG,
                )
    st.subheader("計算結果")
    st.dataframe(rows, width="stretch", hide_index=True)
    left, right = st.columns(2)
    left.download_button(
        "結果CSVをダウンロード",
        rows_to_csv(rows),
        file_name=f"{filename_prefix}.csv",
        mime="text/csv",
        key=f"{filename_prefix}_csv",
    )
    export_configuration = dict(configuration)
    export_configuration["display_units"] = preferences.to_dict()
    right.download_button(
        "設定JSONをダウンロード",
        dump_configuration(export_configuration),
        file_name=f"{filename_prefix}.json",
        mime="application/json",
        key=f"{filename_prefix}_json",
    )


def finite_number(
    label: str,
    value: float,
    *,
    key: str,
    min_value: float | None = None,
    disabled: bool = False,
    help: str | None = None,
    format: str = "%.6g",
) -> float:
    """Render a consistently configured floating-point input."""
    return float(
        st.number_input(
            label,
            value=float(value),
            min_value=min_value,
            key=key,
            disabled=disabled,
            help=help,
            format=format,
        )
    )
