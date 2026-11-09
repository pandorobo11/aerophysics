"""Guard the path rules that decide when solver verification runs."""

import ast
import os
import shutil
import subprocess
import sys
from fnmatch import fnmatchcase
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
TESTS = ROOT / "tests"
PACKAGE = ROOT / "src" / "aerophysics"
PATHS_FILE = ROOT / "scripts" / "solver-verification-paths.txt"
MARKER = "solver_verification"


def solver_patterns() -> list[str]:
    lines = PATHS_FILE.read_text(encoding="utf-8").splitlines()
    return [line for line in lines if line and not line.startswith("#")]


def is_covered(path: str) -> bool:
    # Bash `[[ path == pattern ]]` also lets `*` match `/`, like fnmatchcase.
    return any(fnmatchcase(path, pattern) for pattern in solver_patterns())


def package_exports() -> dict[str, str]:
    tree = ast.parse((PACKAGE / "__init__.py").read_text(encoding="utf-8"))
    exports: dict[str, str] = {}
    for node in tree.body:
        if isinstance(node, ast.ImportFrom) and node.module:
            module = node.module.removeprefix("aerophysics.")
            for alias in node.names:
                exports[alias.asname or alias.name] = module
    return exports


def imported_modules(tree: ast.AST, exports: dict[str, str]) -> set[str]:
    """Return aerophysics submodules and local test modules imported by ``tree``."""
    modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            if node.module == "aerophysics":
                for alias in node.names:
                    modules.add(exports.get(alias.name, alias.name))
            elif node.module.startswith("aerophysics."):
                modules.add(node.module.removeprefix("aerophysics."))
            elif node.module.startswith("test_"):
                modules.add(f"tests:{node.module}")
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith("aerophysics."):
                    modules.add(alias.name.removeprefix("aerophysics."))
                elif alias.name.startswith("test_"):
                    modules.add(f"tests:{alias.name}")
    return modules


def module_file(module: str) -> Path:
    if module.startswith("tests:"):
        return TESTS / f"{module.removeprefix('tests:')}.py"
    candidate = PACKAGE / f"{module.replace('.', '/')}.py"
    if candidate.exists():
        return candidate
    return PACKAGE / module.replace(".", "/") / "__init__.py"


def solver_inputs() -> set[Path]:
    """Files that can change the outcome of the solver verification tier."""
    exports = package_exports()
    roots = [
        path
        for path in sorted(TESTS.glob("test_*.py"))
        if f"mark.{MARKER}" in path.read_text(encoding="utf-8")
    ]
    assert roots, f"no tests are marked {MARKER}"
    inputs: set[Path] = {PACKAGE / "__init__.py"}
    pending = list(roots)
    while pending:
        path = pending.pop()
        if path in inputs:
            continue
        inputs.add(path)
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Constant)
                and isinstance(node.value, str)
                and node.value.startswith("reference_data/")
            ):
                inputs.add(TESTS / node.value)
        for module in imported_modules(tree, exports):
            target = module_file(module)
            assert target.exists(), f"{path.name} imports unknown module {module}"
            pending.append(target)
    return inputs


def test_solver_verification_inputs_are_covered_by_path_rules() -> None:
    missing = sorted(
        path.relative_to(ROOT).as_posix()
        for path in solver_inputs()
        if not is_covered(path.relative_to(ROOT).as_posix())
    )
    assert not missing, (
        f"add these solver_verification inputs to {PATHS_FILE.name}: {missing}"
    )


def test_solver_path_rules_have_no_stale_entries() -> None:
    candidates = [path for path in ROOT.iterdir() if path.is_file()]
    for directory in ("src", "tests", "scripts", ".github"):
        candidates.extend(
            path for path in (ROOT / directory).rglob("*") if path.is_file()
        )
    tracked = [path.relative_to(ROOT).as_posix() for path in candidates]
    stale = [
        pattern
        for pattern in solver_patterns()
        if not any(fnmatchcase(path, pattern) for path in tracked)
    ]
    assert not stale, f"patterns in {PATHS_FILE.name} match no file: {stale}"


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX shell script")
@pytest.mark.parametrize(
    ("changed", "expected"),
    [
        ("src/aerophysics/gas.py", "yes"),
        ("src/aerophysics/heating.py", "no"),
        ("docs/scripts/generate_viscous_flow_validation.py", "no"),
    ],
)
def test_needs_verification_uses_path_rules(
    tmp_path: Path, changed: str, expected: str
) -> None:
    git = shutil.which("git")
    bash = shutil.which("bash")
    if git is None or bash is None:
        pytest.skip("git and bash are required")
    scripts = tmp_path / "scripts"
    scripts.mkdir()
    for name in ("needs-verification.sh", PATHS_FILE.name):
        shutil.copy2(ROOT / "scripts" / name, scripts / name)
    target = tmp_path / changed
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("before\n", encoding="utf-8")

    def run(*args: str) -> str:
        result = subprocess.run(
            [git, "-c", "user.name=test", "-c", "user.email=test@example.com", *args],
            cwd=tmp_path,
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout

    run("init", "--quiet", "--initial-branch=main")
    run("add", ".")
    run("commit", "--quiet", "-m", "base")
    target.write_text("after\n", encoding="utf-8")

    result = subprocess.run(
        [bash, "scripts/needs-verification.sh", "solver"],
        cwd=tmp_path,
        env={**os.environ, "VERIFICATION_BASE": "main"},
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.stdout.strip() == expected, result.stderr
