"""Security regression checks for the tag-triggered release workflow."""

import re
import tomllib
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = PROJECT_ROOT / ".github" / "workflows" / "release.yml"


def _workflow_text() -> str:
    return WORKFLOW.read_text(encoding="utf-8")


def _job_sections() -> tuple[str, str, str]:
    workflow = _workflow_text()
    build = workflow.split("\n  build:\n", maxsplit=1)[1]
    build, attest = build.split("\n  attest:\n", maxsplit=1)
    attest, publish = attest.split("\n  publish:\n", maxsplit=1)
    return build, attest, publish


def test_release_actions_are_pinned_to_immutable_commits() -> None:
    workflow = _workflow_text()
    uses_lines = re.findall(r"^\s*(?:-\s+)?uses:\s+.+$", workflow, flags=re.MULTILINE)
    uses = re.findall(
        r"^\s*(?:-\s+)?uses:\s+([^\s]+)(?:\s+#\s+(v[^\s]+))?$",
        workflow,
        flags=re.MULTILINE,
    )
    assert uses
    assert len(uses) == len(uses_lines)
    for action_reference, version_comment in uses:
        action, separator, revision = action_reference.partition("@")
        assert action
        assert separator == "@"
        assert re.fullmatch(r"[0-9a-f]{40}", revision)
        assert version_comment


def test_build_cannot_publish_and_checkout_does_not_persist_credentials() -> None:
    workflow = _workflow_text()
    build, attest, publish = _job_sections()

    assert "permissions:\n  contents: read" in workflow
    assert "contents: write" not in build
    assert "checks: read" in build
    assert "attestations: write" not in build
    assert "id-token: write" not in build
    assert "persist-credentials: false" in build

    assert "contents: write" not in attest
    assert "attestations: write" in attest
    assert "id-token: write" in attest
    assert "actions/checkout" not in attest

    assert workflow.count("contents: write") == 1
    assert "contents: write" in publish
    assert "actions/checkout" not in publish
    assert "GH_REPO: ${{ github.repository }}" in publish
    assert "- build" in publish
    assert "- attest" in publish


def test_tagged_sha_must_have_a_successful_ci_gate() -> None:
    build, _, _ = _job_sections()
    assert "commits/${GITHUB_SHA}/check-runs" in build
    assert "check_name=CI%20gate" in build
    assert 'check["head_sha"] == sha' in build
    assert 'check.get("app", {}).get("slug") == "github-actions"' in build
    assert 'check["conclusion"] == "success"' in build
    assert 'git rev-list -n 1 "refs/tags/${GITHUB_REF_NAME}"' in build
    assert '"${tag_commit}" != "${GITHUB_SHA}"' in build


def test_isolated_build_dependencies_are_hash_locked() -> None:
    configuration = tomllib.loads(
        (PROJECT_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    )
    assert any(
        requirement.startswith("hatchling")
        for requirement in configuration["dependency-groups"]["build"]
    )

    build_script = (PROJECT_ROOT / "scripts" / "build-distributions.sh").read_text(
        encoding="utf-8"
    )
    for required_option in (
        "uv export",
        "--all-groups",
        "--all-extras",
        "--no-emit-project",
        "--locked",
        "env -u PYTHONPATH uv build",
        "uv build",
        "--build-constraints",
        "--require-hashes",
    ):
        assert required_option in build_script

    build, _, _ = _job_sections()
    assert 'UV_NO_SYNC: "1"' in build
    assert "PYTHONPATH: ${{ github.workspace }}/src" in build
    assert "uv sync --all-groups --all-extras --locked --no-install-project" in build
    assert "--no-build" in build
    assert "bash scripts/build-distributions.sh" in build
    ci_workflow = (PROJECT_ROOT / ".github" / "workflows" / "ci.yml").read_text(
        encoding="utf-8"
    )
    ci_package = ci_workflow.split("\n  package:\n", maxsplit=1)[1].split(
        "\n  ci-gate:\n", maxsplit=1
    )[0]
    assert 'UV_NO_SYNC: "1"' in ci_package
    assert "PYTHONPATH: ${{ github.workspace }}/src" in ci_package
    assert (
        "uv sync --all-groups --all-extras --locked --no-install-project" in ci_package
    )
    assert "--no-build" in ci_package
    assert "bash scripts/build-distributions.sh" in ci_package
    assert "bash scripts/build-distributions.sh" in (
        PROJECT_ROOT / "scripts" / "check.sh"
    ).read_text(encoding="utf-8")


def test_release_bundle_is_checksummed_attested_and_transferred() -> None:
    build, attest, publish = _job_sections()

    assert "sha256sum --check SHA256SUMS" in build
    assert "actions/upload-artifact@" in build
    assert "if-no-files-found: error" in build

    assert "actions/download-artifact@" in attest
    assert "sha256sum --check SHA256SUMS" in attest
    assert "actions/attest-build-provenance@" in attest
    assert "subject-path: release-bundle/*" in attest

    assert "actions/download-artifact@" in publish
    assert "sha256sum --check SHA256SUMS" in publish
    assert '"SHA256SUMS"' in publish
    assert "actual != expected" in publish


def test_existing_release_assets_can_never_be_clobbered() -> None:
    _, _, publish = _job_sections()
    assert "gh release view" in publish
    assert "already exists; refusing to replace assets" in publish
    assert publish.count("gh release create") == 1
    assert "gh release upload" not in publish
    assert "--clobber" not in publish
    assert '"release-bundle/SHA256SUMS"' in publish


def test_tag_is_dereferenced_again_immediately_before_publication() -> None:
    _, _, publish = _job_sections()
    check_start = publish.index('"repos/${GH_REPO}/git/ref/tags/${GITHUB_REF_NAME}"')
    create_start = publish.index("gh release create")
    assert check_start < create_start
    assert '"repos/${GH_REPO}/git/tags/${tag_sha}"' in publish
    assert '"${tag_type}" == "commit"' in publish
    assert '"${tag_type}" != "tag"' in publish
    assert '"${tag_sha}" != "${GITHUB_SHA}"' in publish
