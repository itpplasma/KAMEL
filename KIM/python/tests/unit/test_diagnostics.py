from __future__ import annotations

from pathlib import Path

import pytest
from kim import __version__
from kim.diagnostics import diagnose_environment


def make_executable(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\nexit 91\n", encoding="utf-8")
    path.chmod(0o755)
    return path


def test_diagnose_environment_reports_selection_without_running_executable(
    tmp_path: Path,
) -> None:
    executable = make_executable(tmp_path / "KIM.x")

    report = diagnose_environment(executable)

    assert report.schema_version == 1
    assert report.python_version
    assert report.package_version == __version__
    assert report.selection_status == "selected"
    assert report.executable_path == executable.resolve()
    assert report.executable_source == "explicit"
    assert report.selection_error is None
    assert report.runtime_libraries_checked is False
    assert "does not verify runtime libraries" in report.runtime_check


def test_diagnose_environment_reports_missing_executable_without_running_fallback(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    missing = tmp_path / "missing" / "KIM.x"
    fallback = make_executable(tmp_path / "fallback" / "KIM.x")

    monkeypatch.setenv("KIM_EXECUTABLE", str(fallback))
    report = diagnose_environment(missing)

    assert report.selection_status == "unavailable"
    assert report.executable_path is None
    assert report.executable_source is None
    assert report.selection_error is not None
    assert "explicit executable" in report.selection_error
