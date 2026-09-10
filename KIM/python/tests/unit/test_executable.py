from __future__ import annotations

import subprocess
from pathlib import Path

import pytest
from kim.errors import ExecutableError
from kim.executable import (
    ExecutableSelection,
    discover_kamel_git_metadata,
    executable_sha256,
    resolve_executable,
    select_executable,
)


def make_executable(path: Path, contents: bytes = b"#!/bin/sh\nexit 0\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(contents)
    path.chmod(0o755)
    return path


def make_checkout(path: Path) -> Path:
    (path / "KIM" / "src").mkdir(parents=True)
    (path / "CMakeLists.txt").write_text("project(KAMEL)\n")
    (path / "KIM" / "src" / "CMakeLists.txt").write_text("add_executable(KIM.x)\n")
    return path


def test_explicit_path_has_highest_precedence_and_returns_absolute_path(
    tmp_path: Path,
) -> None:
    explicit = make_executable(tmp_path / "explicit" / "KIM.x")
    environment = make_executable(tmp_path / "environment" / "KIM.x")
    path_entry = make_executable(tmp_path / "path" / "KIM.x")

    resolved = resolve_executable(
        explicit,
        environment={"KIM_EXECUTABLE": str(environment), "PATH": str(path_entry.parent)},
        start_directory=tmp_path,
    )

    assert resolved == explicit.resolve()
    assert resolved.is_absolute()


def test_select_executable_reports_explicit_source(tmp_path: Path) -> None:
    executable = make_executable(tmp_path / "explicit" / "KIM.x")

    selection = select_executable(executable, environment={}, start_directory=tmp_path)

    assert selection == ExecutableSelection(path=executable.resolve(), source="explicit")


def test_select_executable_reports_environment_source(tmp_path: Path) -> None:
    executable = make_executable(tmp_path / "environment" / "KIM.x")

    selection = select_executable(
        environment={"KIM_EXECUTABLE": str(executable), "PATH": ""},
        start_directory=tmp_path,
    )

    assert selection.path == executable.resolve()
    assert selection.source == "environment"


def test_select_executable_reports_checkout_source(tmp_path: Path) -> None:
    checkout = make_checkout(tmp_path / "checkout")
    executable = make_executable(checkout / "build" / "install" / "bin" / "KIM.x")

    selection = select_executable(environment={"PATH": ""}, start_directory=checkout / "KIM")

    assert selection.path == executable.resolve()
    assert selection.source == "checkout"


def test_select_executable_reports_path_source(tmp_path: Path) -> None:
    executable = make_executable(tmp_path / "bin" / "KIM.x")

    selection = select_executable(
        environment={"PATH": str(executable.parent)}, start_directory=tmp_path
    )

    assert selection.path == executable.resolve()
    assert selection.source == "path"


def test_invalid_explicit_path_does_not_fall_back_to_environment(tmp_path: Path) -> None:
    environment = make_executable(tmp_path / "environment" / "KIM.x")
    missing = tmp_path / "missing" / "KIM.x"

    with pytest.raises(ExecutableError, match=r"explicit executable.*does not name a file"):
        select_executable(
            missing,
            environment={"KIM_EXECUTABLE": str(environment), "PATH": str(environment.parent)},
            start_directory=tmp_path,
        )


def test_environment_variable_precedes_checkout_and_path(tmp_path: Path) -> None:
    checkout = make_checkout(tmp_path / "checkout")
    checkout_executable = make_executable(checkout / "build" / "install" / "bin" / "KIM.x")
    environment_executable = make_executable(tmp_path / "environment" / "KIM.x")
    path_executable = make_executable(tmp_path / "path" / "KIM.x")

    resolved = resolve_executable(
        environment={
            "KIM_EXECUTABLE": str(environment_executable),
            "PATH": str(path_executable.parent),
        },
        start_directory=checkout / "KIM" / "src",
    )

    assert resolved == environment_executable.resolve()
    assert resolved != checkout_executable.resolve()


def test_checkout_build_precedes_path_lookup(tmp_path: Path) -> None:
    checkout = make_checkout(tmp_path / "checkout")
    checkout_executable = make_executable(checkout / "build" / "install" / "bin" / "KIM.x")
    path_executable = make_executable(tmp_path / "path" / "KIM.x")

    resolved = resolve_executable(
        environment={"PATH": str(path_executable.parent)},
        start_directory=checkout / "KIM" / "src",
    )

    assert resolved == checkout_executable.resolve()


def test_path_lookup_searches_only_for_uppercase_kim_x(tmp_path: Path) -> None:
    executable = make_executable(tmp_path / "bin" / "KIM.x")

    resolved = resolve_executable(
        environment={"PATH": str(executable.parent)},
        start_directory=tmp_path,
    )

    assert resolved == executable.resolve()


def test_missing_executable_reports_all_supported_resolution_methods(tmp_path: Path) -> None:
    with pytest.raises(
        ExecutableError,
        match=r"KIM_EXECUTABLE.*build/install/bin/KIM\.x.*PATH",
    ):
        resolve_executable(environment={"PATH": ""}, start_directory=tmp_path)


@pytest.mark.parametrize("source", ["explicit", "environment"])
def test_selected_non_executable_file_is_rejected(tmp_path: Path, source: str) -> None:
    candidate = tmp_path / "KIM.x"
    candidate.write_text("not executable\n")
    candidate.chmod(0o644)
    arguments = (
        {"explicit": candidate, "environment": {}}
        if source == "explicit"
        else {"environment": {"KIM_EXECUTABLE": str(candidate)}}
    )

    with pytest.raises(ExecutableError, match=r"KIM\.x.*not executable"):
        resolve_executable(start_directory=tmp_path, **arguments)


@pytest.mark.parametrize("source", ["explicit", "environment"])
def test_python_cli_named_kim_is_refused(tmp_path: Path, source: str) -> None:
    cli = make_executable(tmp_path / "kim")
    arguments = (
        {"explicit": cli, "environment": {}}
        if source == "explicit"
        else {"environment": {"KIM_EXECUTABLE": str(cli)}}
    )

    with pytest.raises(ExecutableError, match=r"Python CLI.*kim.*KIM\.x"):
        resolve_executable(start_directory=tmp_path, **arguments)


def test_lowercase_kim_on_path_is_never_resolved(tmp_path: Path) -> None:
    make_executable(tmp_path / "bin" / "kim")

    with pytest.raises(ExecutableError):
        resolve_executable(
            environment={"PATH": str(tmp_path / "bin")},
            start_directory=tmp_path,
        )


def test_executable_sha256_uses_file_content(tmp_path: Path) -> None:
    executable = make_executable(tmp_path / "KIM.x", b"known KIM executable bytes")

    assert executable_sha256(executable) == (
        "09321b9d6b56c0485604bf4795848bdbbbc7d3704cf818ba821226635115f9ef"
    )


def test_git_metadata_is_discovered_from_a_checkout_descendant(tmp_path: Path) -> None:
    checkout = make_checkout(tmp_path / "checkout")
    subprocess.run(["git", "init", "-q", str(checkout)], check=True)
    subprocess.run(
        ["git", "-C", str(checkout), "config", "user.email", "test@example.invalid"],
        check=True,
    )
    subprocess.run(["git", "-C", str(checkout), "config", "user.name", "KIM Test"], check=True)
    subprocess.run(["git", "-C", str(checkout), "add", "."], check=True)
    subprocess.run(["git", "-C", str(checkout), "commit", "-q", "-m", "initial"], check=True)
    expected_commit = subprocess.run(
        ["git", "-C", str(checkout), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()

    clean = discover_kamel_git_metadata(checkout / "KIM" / "src")
    assert clean is not None
    assert clean.root == checkout.resolve()
    assert clean.commit == expected_commit
    assert clean.dirty is False

    (checkout / "CMakeLists.txt").write_text("project(KAMEL changed)\n")
    dirty = discover_kamel_git_metadata(checkout / "KIM")
    assert dirty is not None
    assert dirty.dirty is True


def test_git_metadata_is_optional_outside_a_checkout(tmp_path: Path) -> None:
    assert discover_kamel_git_metadata(tmp_path) is None
