"""Tests for the GUI documentation browser and local server."""

import re
from pathlib import Path
from unittest.mock import patch
from urllib.request import urlopen

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


def test_page_documentation_link(monkeypatch: pytest.MonkeyPatch) -> None:
    script = """
from aerophysics.gui.documentation import render_documentation_link
render_documentation_link("models/conical_shocks.html")
"""
    monkeypatch.delenv(DOCS_URL_ENV, raising=False)
    app = AppTest.from_string(script, default_timeout=15).run()
    assert not app.exception
    assert len(app.get("link_button")) == 0

    monkeypatch.setenv(DOCS_URL_ENV, "http://127.0.0.1:8123")
    app = AppTest.from_string(script, default_timeout=15).run()
    assert not app.exception
    buttons = app.get("link_button")
    assert len(buttons) == 1
    assert buttons[0].proto.url == "http://127.0.0.1:8123/models/conical_shocks.html"


def test_every_page_documentation_link_is_a_topic() -> None:
    package = Path(__file__).resolve().parents[2] / "src/aerophysics/gui"
    links = {
        match.group(1)
        for path in package.glob("*.py")
        for match in re.finditer(
            r'render_documentation_link\("([^"]+)"\)', path.read_text("utf-8")
        )
    }
    assert len(links) >= 12
    assert links <= set(DOCUMENTATION_TOPICS.values())


def test_local_documentation_server(tmp_path: Path) -> None:
    (tmp_path / "index.html").write_text("docs-ready", encoding="utf-8")
    server, thread = _start_documentation_server(tmp_path)
    try:
        host, port = server.server_address[:2]
        assert host == LOOPBACK_ADDRESS
        with urlopen(f"http://{LOOPBACK_ADDRESS}:{port}/index.html", timeout=5) as page:
            assert page.read() == b"docs-ready"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


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
