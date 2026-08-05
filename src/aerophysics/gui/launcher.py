"""Console launcher for the optional Streamlit application."""

from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

DOCS_URL_ENV = "AEROPHYSICS_DOCS_URL"


class _QuietDocumentationHandler(SimpleHTTPRequestHandler):
    """Serve generated documentation without writing request logs to stderr."""

    def log_message(self, format: str, *args: object) -> None:
        pass


def _documentation_directory() -> Path | None:
    configured = os.environ.get("AEROPHYSICS_DOCS_DIR")
    candidates = [] if configured is None else [Path(configured).expanduser()]
    candidates.append(Path(__file__).resolve().parents[3] / "docs" / "_build" / "html")
    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved.is_dir() and (resolved / "index.html").is_file():
            return resolved
    return None


def _start_documentation_server(
    directory: Path,
) -> tuple[ThreadingHTTPServer, threading.Thread]:
    handler = partial(_QuietDocumentationHandler, directory=str(directory))
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    return server, thread


def main() -> None:
    """Run the bundled Streamlit application."""
    if importlib.util.find_spec("streamlit") is None:
        raise SystemExit(
            "GUI dependencies are not installed. "
            "Install them with: python -m pip install 'aerophysics[gui]'"
        )
    app = Path(__file__).with_name("app.py")
    command = [sys.executable, "-m", "streamlit", "run", str(app), *sys.argv[1:]]
    environment = dict(os.environ)
    docs_server: ThreadingHTTPServer | None = None
    docs_thread: threading.Thread | None = None
    directory = _documentation_directory()
    if directory is not None:
        docs_server, docs_thread = _start_documentation_server(directory)
        port = docs_server.server_address[1]
        environment[DOCS_URL_ENV] = f"http://127.0.0.1:{port}/"
    try:
        try:
            returncode = subprocess.run(
                command,
                check=False,
                env=environment,
            ).returncode
        except KeyboardInterrupt:
            returncode = 130
    finally:
        if docs_server is not None:
            docs_server.shutdown()
            docs_server.server_close()
        if docs_thread is not None:
            docs_thread.join(timeout=2.0)
    raise SystemExit(returncode)


if __name__ == "__main__":
    main()
