"""Tests for the GUI documentation browser and local server."""

from pathlib import Path
from unittest.mock import patch

import pytest
from streamlit.testing.v1 import AppTest

from aerophysics.gui.documentation import (
    DOCS_URL_ENV,
    documentation_base_url,
    documentation_topic_url,
)
from aerophysics.gui.launcher import (
    _documentation_directory,
    _start_documentation_server,
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
            "compressible_flow.html", documentation_base_url() or ""
        )
        == "http://127.0.0.1:8123/compressible_flow.html"
    )
    with pytest.raises(ValueError, match="unsupported documentation topic"):
        documentation_topic_url("../index.html", "http://127.0.0.1:8123/")


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
    assert address == ("127.0.0.1", 0)
    assert callable(handler)
    thread_class.assert_called_once_with(target=server.serve_forever, daemon=True)
    thread_class.return_value.start.assert_called_once_with()


def test_documentation_directory_override(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("AEROPHYSICS_DOCS_DIR", str(tmp_path))
    assert _documentation_directory() != tmp_path
    (tmp_path / "index.html").write_text("docs", encoding="utf-8")
    assert _documentation_directory() == tmp_path.resolve()
