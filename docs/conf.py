"""Sphinx configuration for aerophysics."""

import tomllib
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as package_version
from pathlib import Path


def _project_version() -> str:
    try:
        return package_version("aerophysics")
    except PackageNotFoundError:
        pyproject = Path(__file__).parents[1] / "pyproject.toml"
        with pyproject.open("rb") as stream:
            project = tomllib.load(stream)["project"]
        return str(project["version"])


project = "aerophysics"
copyright = "2026, aerophysics contributors"
release = _project_version()
version = release

extensions = [
    "numpydoc",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.mathjax",
]
autosummary_generate = True
numpydoc_show_class_members = False
nitpicky = True
nitpick_ignore_regex = [
    ("py:class", r"aerophysics\._array\.Float(Array|Result)"),
    ("py:class", r"aerophysics\.isentropic\.IsentropicGasModel"),
    ("py:class", r"numpy\._typing\..*\.ArrayLike"),
]

html_theme = "pydata_sphinx_theme"
html_title = f"aerophysics {release}"
html_theme_options = {
    "github_url": "https://github.com/pandorobo11/aerophysics",
}

exclude_patterns = ["_build", "_generated/**"]
