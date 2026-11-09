"""Keep the marching boundary-layer solver's dependency boundary explicit.

The solver core may import only the shared foundation modules listed here.
A new dependency is a design decision: extend the allowlist in the same change,
and check that the solver-verification path rules still cover the new input.
Imports under ``if TYPE_CHECKING:`` are ignored because they never run.
"""

import ast
from pathlib import Path

import pytest

PACKAGE = Path(__file__).resolve().parents[1] / "src" / "aerophysics"

SOLVER_CORE = frozenset(
    {
        "marching_boundary_layer",
        "_boundary_layer_march",
        "_boundary_layer_sa",
        "_banded",
    }
)
PRIVATE_CORE = frozenset(name for name in SOLVER_CORE if name.startswith("_"))
FOUNDATION = frozenset({"_array", "exceptions", "gas", "transport"})
# Modules the core may use for selected names only.
NAMED_ONLY = {"boundary_layer": frozenset({"BoundaryLayerRegime"})}


def is_type_checking_block(node: ast.AST) -> bool:
    if not isinstance(node, ast.If):
        return False
    test = node.test
    return (isinstance(test, ast.Name) and test.id == "TYPE_CHECKING") or (
        isinstance(test, ast.Attribute) and test.attr == "TYPE_CHECKING"
    )


def runtime_imports(path: Path) -> list[tuple[str, frozenset[str]]]:
    """Return ``(submodule, names)`` imported at runtime from ``aerophysics``."""
    imports: list[tuple[str, frozenset[str]]] = []
    pending: list[ast.AST] = [ast.parse(path.read_text(encoding="utf-8"))]
    while pending:
        node = pending.pop()
        if isinstance(node, ast.If) and is_type_checking_block(node):
            pending.extend(node.orelse)
            continue
        if isinstance(node, ast.ImportFrom) and node.module:
            names = frozenset(alias.name for alias in node.names)
            if node.module == "aerophysics":
                imports.extend((name, frozenset()) for name in names)
            elif node.module.startswith("aerophysics."):
                imports.append((node.module.removeprefix("aerophysics."), names))
        elif isinstance(node, ast.Import):
            imports.extend(
                (alias.name.removeprefix("aerophysics."), frozenset())
                for alias in node.names
                if alias.name.startswith("aerophysics.")
            )
        pending.extend(ast.iter_child_nodes(node))
    return imports


@pytest.mark.parametrize("module", sorted(SOLVER_CORE))
def test_solver_core_imports_only_its_foundation(module: str) -> None:
    violations = []
    for imported, names in runtime_imports(PACKAGE / f"{module}.py"):
        if imported in SOLVER_CORE or imported in FOUNDATION:
            continue
        allowed = NAMED_ONLY.get(imported)
        if allowed is not None and names and names <= allowed:
            continue
        violations.append(f"{imported}: {sorted(names) or 'module'}")
    assert not violations, f"{module} imports outside the solver boundary: {violations}"


def test_private_solver_modules_stay_behind_the_public_module() -> None:
    importers = sorted(
        f"{path.relative_to(PACKAGE).as_posix()} -> {imported}"
        for path in PACKAGE.rglob("*.py")
        if path.stem not in SOLVER_CORE
        for imported, _ in runtime_imports(path)
        if imported in PRIVATE_CORE
    )
    assert not importers, (
        "use aerophysics.marching_boundary_layer instead of private solver "
        f"modules: {importers}"
    )


def test_type_checking_and_runtime_imports(tmp_path: Path) -> None:
    source = tmp_path / "example.py"
    source.write_text(
        "from typing import TYPE_CHECKING\n"
        "import typing\n"
        "from aerophysics.gas import AIR\n"
        "if TYPE_CHECKING:\n"
        "    from aerophysics.shocks import oblique_shock\n"
        "else:\n"
        "    from aerophysics import _boundary_layer_sa\n"
        "if typing.TYPE_CHECKING:\n"
        "    from aerophysics.heating import x\n"
        "else:\n"
        "    from aerophysics import _boundary_layer_march\n"
        "def lazy():\n"
        "    from aerophysics import _banded\n",
        encoding="utf-8",
    )
    assert sorted(runtime_imports(source)) == [
        ("_banded", frozenset()),
        ("_boundary_layer_march", frozenset()),
        ("_boundary_layer_sa", frozenset()),
        ("gas", frozenset({"AIR"})),
    ]
