"""Documentation browser embedded in the Streamlit GUI."""

from __future__ import annotations

import os
from urllib.parse import urljoin, urlsplit

import streamlit as st

DOCS_URL_ENV = "AEROPHYSICS_DOCS_URL"
RELEASES_URL = "https://github.com/pandorobo11/aerophysics/releases/latest"

DOCUMENTATION_TOPICS = {
    "概要": "index.html",
    "クイックスタート": "quickstart.html",
    "気体・標準大気": "gas_and_atmosphere.html",
    "輸送物性": "transport_properties.html",
    "熱化学": "thermochemistry.html",
    "圧縮性流れ": "compressible_flow.html",
    "境界層・突起抗力": "boundary_layers.html",
    "圧縮性速度変換": "compressible_velocity_transformations.html",
    "飛行条件": "flight_conditions.html",
    "単位変換": "unit_conversions.html",
    "検証": "verification.html",
    "APIリファレンス": "api.html",
    "参考文献": "references.html",
}


def documentation_base_url() -> str | None:
    """Return the validated URL of locally served Sphinx documentation."""
    value = os.environ.get(DOCS_URL_ENV, "").strip()
    parsed = urlsplit(value)
    if parsed.scheme not in {"http", "https"} or not parsed.netloc:
        return None
    return value.rstrip("/") + "/"


def documentation_topic_url(filename: str, base_url: str) -> str:
    """Return one documentation page below a trusted base URL."""
    if filename not in DOCUMENTATION_TOPICS.values():
        raise ValueError("unsupported documentation topic")
    return urljoin(base_url, filename)


def render_documentation() -> None:
    """Render local Sphinx documentation or installation guidance."""
    st.title("ドキュメント")
    st.caption("数式、仮定、適用範囲、API、検証資料をGUI内で参照できます。")
    base_url = documentation_base_url()
    if base_url is None:
        st.warning(
            "ローカルHTMLドキュメントが見つかりません。開発環境ではSphinxで"
            "HTMLを生成してから `aerophysics-gui` で起動してください。"
        )
        st.code("uv run sphinx-build -W -b html docs docs/_build/html")
        st.link_button("最新リリースのドキュメントを取得", RELEASES_URL)
        return

    topic = st.selectbox(
        "表示する項目",
        tuple(DOCUMENTATION_TOPICS),
        key="documentation_topic",
    )
    assert topic is not None
    target = documentation_topic_url(DOCUMENTATION_TOPICS[topic], base_url)
    st.link_button("新しいタブで開く", target)
    st.iframe(target, height=900)
