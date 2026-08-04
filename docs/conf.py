"""Sphinx configuration for aerophysics."""

from importlib.metadata import version as package_version

project = "aerophysics"
copyright = "2026, aerophysics contributors"
release = package_version("aerophysics")
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
    ("py:class", r"numpy\._typing\..*\.ArrayLike"),
]

html_theme = "pydata_sphinx_theme"
html_title = f"aerophysics {release}"
html_theme_options = {
    "github_url": "https://github.com/pandorobo11/aerophysics",
}

exclude_patterns = ["_build", "_generated/**"]
