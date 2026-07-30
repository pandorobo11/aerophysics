"""Streamlit entry point for the local aerophysics GUI."""

from __future__ import annotations

import streamlit as st

from aerophysics.gui.components import render_unit_sidebar
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


pages = {
    "計算": [
        st.Page(
            flight_page,
            title="大気・飛行条件",
            icon="🌤️",
            default=True,
        ),
        st.Page(
            shock_page,
            title="斜め衝撃波",
            icon="〰️",
        ),
        st.Page(
            boundary_layer_page,
            title="平板境界層",
            icon="📐",
        ),
    ]
}
st.navigation(pages).run()
