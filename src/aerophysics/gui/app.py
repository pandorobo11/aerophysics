"""Streamlit entry point for the local aerophysics GUI."""

from __future__ import annotations

import streamlit as st

from aerophysics.gui.analysis_pages import (
    render_boundary_layer_profile,
    render_protrusion_drag,
    render_thermochemistry,
)
from aerophysics.gui.components import render_unit_sidebar
from aerophysics.gui.flow_pages import (
    render_expansion,
    render_isentropic,
    render_normal_shock,
)
from aerophysics.gui.pages import (
    render_boundary_layer,
    render_flight,
    render_shock,
)

st.set_page_config(
    page_title="aerophysics GUI",
    page_icon="✈️",
    layout="wide",
)

preferences = render_unit_sidebar()


def flight_page() -> None:
    render_flight(preferences)


def shock_page() -> None:
    render_shock(preferences)


def boundary_layer_page() -> None:
    render_boundary_layer(preferences)


def isentropic_page() -> None:
    render_isentropic(preferences)


def normal_shock_page() -> None:
    render_normal_shock(preferences)


def expansion_page() -> None:
    render_expansion(preferences)


def boundary_layer_profile_page() -> None:
    render_boundary_layer_profile(preferences)


def protrusion_drag_page() -> None:
    render_protrusion_drag(preferences)


def thermochemistry_page() -> None:
    render_thermochemistry(preferences)


pages = {
    "飛行状態": [
        st.Page(
            flight_page,
            title="大気・飛行条件",
            icon="🌤️",
            default=True,
        ),
    ],
    "圧縮性流れ": [
        st.Page(
            isentropic_page,
            title="等エントロピー流れ",
            icon="📊",
        ),
        st.Page(
            normal_shock_page,
            title="垂直衝撃波",
            icon="↕️",
        ),
        st.Page(
            shock_page,
            title="斜め衝撃波",
            icon="〰️",
        ),
        st.Page(
            expansion_page,
            title="Prandtl–Meyer膨張",
            icon="↗️",
        ),
    ],
    "粘性流れ": [
        st.Page(
            boundary_layer_page,
            title="平板境界層",
            icon="📐",
        ),
        st.Page(
            boundary_layer_profile_page,
            title="圧縮性境界層プロファイル",
            icon="📈",
        ),
        st.Page(
            protrusion_drag_page,
            title="突起抗力",
            icon="🧱",
        ),
    ],
    "熱物性": [
        st.Page(
            thermochemistry_page,
            title="熱化学",
            icon="🔥",
        ),
    ],
}
st.navigation(pages).run()
