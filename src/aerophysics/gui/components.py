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
from aerophysics.gui.units import QuantityKind, UnitPreferences, from_si, to_si

PLOTLY_CONFIG = {
    "displaylogo": False,
    "scrollZoom": False,
    "toImageButtonOptions": {"format": "png", "scale": 2},
}


_UNIT_STATE_KEY = "_gui_display_units"


def _current_unit_preferences() -> UnitPreferences:
    """Build preferences from widget state, including not-yet-rendered defaults."""
    defaults = UnitPreferences()
    return UnitPreferences(
        **{
            name: str(st.session_state.get(f"unit_{name}", value))
            for name, value in defaults.to_dict().items()
        }
    )


def _convert_state_value(
    key: str, kind: QuantityKind, old: UnitPreferences, new: UnitPreferences
) -> None:
    value = st.session_state.get(key)
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return
    old_unit = str(getattr(old, kind))
    new_unit = str(getattr(new, kind))
    if old_unit == new_unit:
        return
    si_value = to_si(float(value), kind, old_unit)
    st.session_state[key] = float(from_si(si_value, kind, new_unit))


def _state_matches(key: str, expected: object) -> bool:
    return st.session_state.get(key) == expected


def clear_widget_state(keys: tuple[str, ...]) -> None:
    """Clear dependent widget values after a structural selector changes."""
    for key in keys:
        st.session_state.pop(key, None)


def _convert_display_input_state() -> None:
    """Preserve physical input values when a global display unit changes."""
    new = _current_unit_preferences()
    stored = st.session_state.get(_UNIT_STATE_KEY)
    try:
        old = UnitPreferences.from_dict(stored)
    except ValueError:
        old = UnitPreferences()

    fields: dict[QuantityKind, tuple[str, ...]] = {
        "length": (
            "flight_altitude",
            "flight_length",
            "detached_shock_radius",
            "boundary_distance",
            "profile_thickness",
            "protrusion_thickness",
            "protrusion_height",
            "protrusion_width",
        ),
        "speed": (
            "boundary_velocity",
            "profile_velocity",
            "protrusion_velocity",
        ),
        "pressure": (
            "isentropic_total_pressure",
            "profile_shear",
            "thermo_pressure",
        ),
        "temperature": (
            "isentropic_total_temperature",
            "boundary_temperature",
            "boundary_wall_temperature",
            "profile_temperature",
            "profile_wall_temperature",
            "protrusion_temperature",
            "protrusion_wall_temperature",
            "thermo_temperature",
            "thermo_reference",
            "thermo_sweep_start",
            "thermo_sweep_stop",
            "viscosity_temperature",
            "viscosity_sweep_start",
            "viscosity_sweep_stop",
        ),
        "density": (
            "boundary_density",
            "profile_density",
            "protrusion_density",
        ),
        "angle": (
            "shock_theta",
            "cone_shock_angle",
            "expansion_turn",
        ),
    }
    for kind, keys in fields.items():
        for key in keys:
            _convert_state_value(key, kind, old, new)

    if _state_matches("flight_basis", "velocity"):
        _convert_state_value("flight_motion", "speed", old, new)
        if _state_matches("flight_sweep_field", "motion"):
            for key in ("flight_sweep_start", "flight_sweep_stop"):
                _convert_state_value(key, "speed", old, new)
    if _state_matches("flight_sweep_field", "altitude"):
        for key in ("flight_sweep_start", "flight_sweep_stop"):
            _convert_state_value(key, "length", old, new)
    for prefix in ("boundary",):
        for suffix in ("sweep_start", "sweep_stop"):
            _convert_state_value(f"{prefix}_{suffix}", "length", old, new)
    if _state_matches("shock_sweep_field", "deflection"):
        for key in ("shock_sweep_start", "shock_sweep_stop"):
            _convert_state_value(key, "angle", old, new)
    if _state_matches("cone_shock_sweep_field", "cone_half_angle"):
        for key in ("cone_shock_sweep_start", "cone_shock_sweep_stop"):
            _convert_state_value(key, "angle", old, new)
    if _state_matches("expansion_sweep_field", "turn_angle"):
        for key in ("expansion_sweep_start", "expansion_sweep_stop"):
            _convert_state_value(key, "angle", old, new)
    if st.session_state.get("protrusion_sweep_field") in {
        "height",
        "base_width",
        "boundary_layer_thickness",
    }:
        for key in ("protrusion_sweep_start", "protrusion_sweep_stop"):
            _convert_state_value(key, "length", old, new)

    st.session_state[_UNIT_STATE_KEY] = new.to_dict()


def render_unit_sidebar() -> UnitPreferences:
    """Render persistent global unit selectors."""
    st.sidebar.header("表示単位")
    with st.sidebar.expander("単位設定", expanded=False):
        length = st.selectbox(
            "長さ",
            ("m", "mm", "ft", "in"),
            key="unit_length",
            on_change=_convert_display_input_state,
        )
        area = st.selectbox(
            "面積",
            ("m²", "ft²", "in²"),
            key="unit_area",
            on_change=_convert_display_input_state,
        )
        speed = st.selectbox(
            "速度",
            ("m/s", "kt", "ft/s"),
            key="unit_speed",
            on_change=_convert_display_input_state,
        )
        pressure = st.selectbox(
            "圧力",
            ("Pa", "kPa", "hPa", "psi", "psf"),
            key="unit_pressure",
            on_change=_convert_display_input_state,
        )
        temperature = st.selectbox(
            "温度",
            ("K", "°C", "°F", "°R"),
            key="unit_temperature",
            on_change=_convert_display_input_state,
        )
        density = st.selectbox(
            "密度",
            ("kg/m³", "slug/ft³", "lbm/ft³"),
            key="unit_density",
            on_change=_convert_display_input_state,
        )
        force = st.selectbox(
            "力",
            ("N", "lbf"),
            key="unit_force",
            on_change=_convert_display_input_state,
        )
        inverse_length = st.selectbox(
            "逆長さ",
            ("1/m", "1/ft", "1/in"),
            key="unit_inverse_length",
            on_change=_convert_display_input_state,
        )
        angle = st.selectbox(
            "角度",
            ("deg", "rad"),
            key="unit_angle",
            on_change=_convert_display_input_state,
        )
    st.sidebar.caption("計算コアへ渡す前にSIへ明示変換します。")
    preferences = UnitPreferences(
        length=length,
        area=area,
        speed=speed,
        pressure=pressure,
        temperature=temperature,
        density=density,
        force=force,
        inverse_length=inverse_length,
        angle=angle,
    )
    st.session_state[_UNIT_STATE_KEY] = preferences.to_dict()
    return preferences


def _set_unit_state(preferences: UnitPreferences) -> None:
    for name, value in preferences.to_dict().items():
        st.session_state[f"unit_{name}"] = value
    st.session_state[_UNIT_STATE_KEY] = preferences.to_dict()


def calculation_button(form_key: str) -> bool:
    """Render an explicit calculation trigger without batching input widgets."""
    return st.button(
        "計算",
        type="primary",
        key=f"FormSubmitter:{form_key}-計算",
    )


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
    if key in st.session_state:
        return float(
            st.number_input(
                label,
                min_value=min_value,
                key=key,
                disabled=disabled,
                help=help,
                format=format,
            )
        )
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
