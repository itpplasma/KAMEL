"""Deterministic discovery and provenance for the KIM scientific executable."""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from kim.errors import ExecutableError


@dataclass(frozen=True)
class KamelGitMetadata:
    """Git identity of a detected KAMEL source checkout."""

    root: Path
    commit: str
    dirty: bool


@dataclass(frozen=True)
class ExecutableSelection:
    """The executable selected and the resolver rule that selected it."""

    path: Path
    source: str


def resolve_executable(
    explicit: Path | str | None = None,
    *,
    environment: Mapping[str, str] | None = None,
    start_directory: Path | str | None = None,
) -> Path:
    """Resolve `KIM.x` according to the documented source precedence.

    This compatibility adapter returns only the selected path. Use
    :func:`select_executable` when the selection source is also needed.
    """

    return select_executable(
        explicit,
        environment=environment,
        start_directory=start_directory,
    ).path


def select_executable(
    explicit: Path | str | None = None,
    *,
    environment: Mapping[str, str] | None = None,
    start_directory: Path | str | None = None,
) -> ExecutableSelection:
    """Select `KIM.x` and report which documented rule supplied it."""

    environ = os.environ if environment is None else environment
    if explicit is not None:
        return ExecutableSelection(
            path=_validate_executable(Path(explicit), source="explicit executable"),
            source="explicit",
        )

    configured = environ.get("KIM_EXECUTABLE", "").strip()
    if configured:
        return ExecutableSelection(
            path=_validate_executable(
                Path(configured), source="KIM_EXECUTABLE environment variable"
            ),
            source="environment",
        )

    checkout = next(
        (
            root
            for start in _search_starts(start_directory)
            if (root := _find_kamel_checkout(start)) is not None
        ),
        None,
    )
    if checkout is not None:
        candidate = checkout / "build" / "install" / "bin" / "KIM.x"
        if candidate.exists():
            return ExecutableSelection(
                path=_validate_executable(candidate, source=f"KAMEL checkout {checkout}"),
                source="checkout",
            )

    located = shutil.which("KIM.x", path=environ.get("PATH", ""))
    if located is not None:
        return ExecutableSelection(
            path=_validate_executable(Path(located), source="PATH"),
            source="path",
        )

    raise ExecutableError(
        "KIM.x was not found: provide an explicit executable, set KIM_EXECUTABLE, "
        "build build/install/bin/KIM.x in a KAMEL checkout, or add KIM.x to PATH"
    )


def executable_sha256(path: Path | str) -> str:
    """Return the SHA-256 digest of an executable's exact bytes."""

    candidate = Path(path)
    if not candidate.is_file():
        raise ExecutableError(f"cannot hash missing executable file: {candidate}")
    digest = hashlib.sha256()
    try:
        with candidate.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as error:
        raise ExecutableError(f"could not read executable {candidate}: {error}") from error
    return digest.hexdigest()


def discover_kamel_git_metadata(
    start_directory: Path | str | None = None,
) -> KamelGitMetadata | None:
    """Return checkout revision metadata, or `None` outside a KAMEL checkout."""

    starts = _search_starts(start_directory)
    checkout = next(
        (root for start in starts if (root := _find_kamel_checkout(start)) is not None),
        None,
    )
    if checkout is None:
        return None
    try:
        commit = _run_git(checkout, "rev-parse", "HEAD")
        status = _run_git(checkout, "status", "--porcelain", "--untracked-files=normal")
    except (FileNotFoundError, subprocess.SubprocessError):
        return None
    if not commit:
        return None
    return KamelGitMetadata(root=checkout.resolve(), commit=commit, dirty=bool(status))


def _search_starts(start_directory: Path | str | None) -> tuple[Path, ...]:
    if start_directory is not None:
        return (Path(start_directory),)
    package_location = Path(__file__).resolve()
    return (Path.cwd(), package_location)


def _find_kamel_checkout(start: Path) -> Path | None:
    candidate = start.resolve()
    if candidate.is_file():
        candidate = candidate.parent
    for directory in (candidate, *candidate.parents):
        if (directory / "CMakeLists.txt").is_file() and (
            directory / "KIM" / "src" / "CMakeLists.txt"
        ).is_file():
            return directory
    return None


def _validate_executable(path: Path, *, source: str) -> Path:
    if path.name.casefold() == "kim":
        raise ExecutableError(
            f"Python CLI name 'kim' is unsafe for {source} {path}; "
            "select the scientific executable KIM.x instead"
        )
    candidate = path.expanduser().resolve()
    if candidate.name.casefold() == "kim":
        raise ExecutableError(
            f"Python CLI name 'kim' is unsafe for {source} {candidate}; select KIM.x instead"
        )
    if not candidate.is_file():
        raise ExecutableError(f"{source} does not name a file: {candidate}")
    if not os.access(candidate, os.X_OK):
        raise ExecutableError(f"{candidate} selected by {source} is not executable")
    return candidate


def _run_git(checkout: Path, *arguments: str) -> str:
    completed = subprocess.run(
        ["git", "-C", str(checkout), *arguments],
        check=True,
        capture_output=True,
        text=True,
        timeout=5,
    )
    return completed.stdout.strip()
