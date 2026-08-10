"""Build hook that bundles the rendered HTML documentation in wheels."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


class CustomBuildHook(BuildHookInterface):
    """Generate Sphinx HTML and add it to the wheel as package resources."""

    def initialize(self, version: str, build_data: dict[str, Any]) -> None:
        """Build the documentation before Hatchling selects wheel files."""
        docs_source = Path(self.root) / "docs"
        temporary_directory = TemporaryDirectory(prefix="aerophysics-docs-")
        docs_output = Path(temporary_directory.name)
        environment = dict(os.environ)
        python_path = [str(Path(self.root) / "src")]
        if existing_python_path := environment.get("PYTHONPATH"):
            python_path.append(existing_python_path)
        environment["PYTHONPATH"] = os.pathsep.join(python_path)

        try:
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "sphinx",
                    "-W",
                    "--keep-going",
                    "-b",
                    "html",
                    str(docs_source),
                    str(docs_output),
                ],
                cwd=self.root,
                env=environment,
                check=True,
            )
        except BaseException:
            temporary_directory.cleanup()
            raise

        self._temporary_directory = temporary_directory
        build_data["force_include"][str(docs_output)] = "aerophysics/_docs"

    def finalize(
        self, version: str, build_data: dict[str, Any], artifact_path: str
    ) -> None:
        """Remove the temporary documentation tree after the wheel is built."""
        temporary_directory = getattr(self, "_temporary_directory", None)
        if temporary_directory is not None:
            temporary_directory.cleanup()
