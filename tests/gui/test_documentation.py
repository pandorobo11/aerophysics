"""Tests for the GUI documentation browser and local server."""

from pathlib import Path
from unittest.mock import patch

import pytest
from streamlit.testing.v1 import AppTest

from aerophysics.gui.documentation import (
    DOCS_URL_ENV,
    DOCUMENTATION_TOPICS,
    documentation_base_url,
    documentation_topic_url,
)
from aerophysics.gui.launcher import (
    DOCS_DIR_ENV,
    LOOPBACK_ADDRESS,
    _documentation_directory,
    _start_documentation_server,
    main,
)


def test_documentation_url_validation(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(DOCS_URL_ENV, raising=False)
    assert documentation_base_url() is None
    monkeypatch.setenv(DOCS_URL_ENV, "file:///tmp/docs")
    assert documentation_base_url() is None
    monkeypatch.setenv(DOCS_URL_ENV, "http://127.0.0.1:8123")
    assert documentation_base_url() == "http://127.0.0.1:8123/"
    assert (
        documentation_topic_url(
            "models/shock_waves.html", documentation_base_url() or ""
        )
        == "http://127.0.0.1:8123/models/shock_waves.html"
    )
    with pytest.raises(ValueError, match="unsupported documentation topic"):
        documentation_topic_url("../index.html", "http://127.0.0.1:8123/")


def test_documentation_topics_follow_current_manual_structure() -> None:
    assert DOCUMENTATION_TOPICS == {
        "概要": "index.html",
        "クイックスタート": "getting_started/quickstart.html",
        "気体・標準大気": "models/gas_and_atmosphere.html",
        "輸送物性": "models/transport_properties.html",
        "熱化学": "models/thermochemistry.html",
        "等エントロピー流れ": "models/isentropic_flow.html",
        "衝撃波": "models/shock_waves.html",
        "膨張波": "models/expansion_waves.html",
        "平板境界層": "models/flat_plate_boundary_layer.html",
        "突起抗力": "models/protrusion_drag.html",
        "圧縮性速度変換": "models/compressible_velocity_transformations.html",
        "飛行条件": "models/flight_conditions.html",
        "単位変換": "models/unit_conversions.html",
        "検証": "verification/index.html",
        "APIリファレンス": "api/index.html",
        "参考文献": "references.html",
        "GUIガイド": "guides/gui.html",
    }


def test_documentation_fallback_page(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(DOCS_URL_ENV, raising=False)
    script = """
from aerophysics.gui.documentation import render_documentation
render_documentation()
"""
    app = AppTest.from_string(script, default_timeout=15).run()
    assert not app.exception
    assert app.title[0].value == "ドキュメント"
    assert app.warning
    assert app.code
    assert len(app.get("link_button")) == 1


def test_local_documentation_server(tmp_path: Path) -> None:
    (tmp_path / "index.html").write_text("docs-ready", encoding="utf-8")
    with (
        patch("aerophysics.gui.launcher.ThreadingHTTPServer") as server_class,
        patch("aerophysics.gui.launcher.threading.Thread") as thread_class,
    ):
        server, _thread = _start_documentation_server(tmp_path)
    address, handler = server_class.call_args.args
    assert address == (LOOPBACK_ADDRESS, 0)
    assert callable(handler)
    thread_class.assert_called_once_with(target=server.serve_forever, daemon=True)
    thread_class.return_value.start.assert_called_once_with()


def test_gui_launcher_binds_streamlit_to_loopback(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        "aerophysics.gui.launcher._documentation_directory", lambda: None
    )
    monkeypatch.setattr("aerophysics.gui.launcher.sys.argv", ["aerophysics-gui"])
    with (
        patch(
            "aerophysics.gui.launcher.importlib.util.find_spec", return_value=object()
        ),
        patch("aerophysics.gui.launcher.subprocess.run") as run,
    ):
        run.return_value.returncode = 0
        with pytest.raises(SystemExit, match="0"):
            main()

    command = run.call_args.args[0]
    assert f"--server.address={LOOPBACK_ADDRESS}" in command


def test_documentation_directory_override(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv(DOCS_DIR_ENV, str(tmp_path))
    assert _documentation_directory() != tmp_path
    (tmp_path / "index.html").write_text("docs", encoding="utf-8")
    assert _documentation_directory() == tmp_path.resolve()
