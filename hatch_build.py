"""Build hook that bundles the rendered HTML documentation in wheels."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


class CustomBuildHook(BuildHookInterface):
    """Generate verification-backed Sphinx HTML for wheel package resources."""

    def initialize(self, version: str, build_data: dict[str, Any]) -> None:
        """Build the documentation before Hatchling selects wheel files."""
        project_root = Path(self.root)
        source_directory = TemporaryDirectory(prefix="aerophysics-docs-source-")
        generated_root = Path(source_directory.name)
        shutil.copytree(
            project_root / "docs",
            generated_root / "docs",
            ignore=shutil.ignore_patterns("_build"),
        )
        shutil.copytree(
            project_root / "src",
            generated_root / "src",
        )
        shutil.copytree(
            project_root / "tests/reference_data",
            generated_root / "tests/reference_data",
        )
        shutil.copy2(project_root / "pyproject.toml", generated_root / "pyproject.toml")

        docs_source = generated_root / "docs"
        output_directory = TemporaryDirectory(prefix="aerophysics-docs-")
        docs_output = Path(output_directory.name)
        environment = dict(os.environ)
        python_path = [str(generated_root / "src")]
        if existing_python_path := environment.get("PYTHONPATH"):
            python_path.append(existing_python_path)
        environment["PYTHONPATH"] = os.pathsep.join(python_path)
        verification_generator = docs_source / "scripts/generate_verification.py"

        try:
            subprocess.run(
                [sys.executable, str(verification_generator)],
                cwd=generated_root,
                env=environment,
                check=True,
            )
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
                cwd=generated_root,
                env=environment,
                check=True,
            )
        except BaseException:
            output_directory.cleanup()
            source_directory.cleanup()
            raise

        self._temporary_directories = (source_directory, output_directory)
        build_data["force_include"][str(docs_output)] = "aerophysics/_docs"

    def finalize(
        self, version: str, build_data: dict[str, Any], artifact_path: str
    ) -> None:
        """Remove the temporary documentation tree after the wheel is built."""
        temporary_directories = getattr(self, "_temporary_directories", ())
        for temporary_directory in temporary_directories:
            temporary_directory.cleanup()
