"""Console launcher for the optional Streamlit application."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path


def main() -> None:
    """Run the bundled Streamlit application."""
    if importlib.util.find_spec("streamlit") is None:
        raise SystemExit(
            "GUI dependencies are not installed. "
            "Install them with: python -m pip install 'aerophysics[gui]'"
        )
    app = Path(__file__).with_name("app.py")
    command = [sys.executable, "-m", "streamlit", "run", str(app), *sys.argv[1:]]
    raise SystemExit(subprocess.run(command, check=False).returncode)


if __name__ == "__main__":
    main()
