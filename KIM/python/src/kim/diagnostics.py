"""Read-only diagnostics for the KIM Python environment."""

from __future__ import annotations

import platform
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from kim.errors import ExecutableError
from kim.executable import ExecutableSelection, select_executable

SelectionStatus = Literal["selected", "unavailable"]
ExecutableSource = Literal["explicit", "environment", "checkout", "path"]


@dataclass(frozen=True)
class EnvironmentDiagnosticReport:
    """Versioned summary of executable selection and environment identity."""

    schema_version: int
    python_version: str
    package_version: str
    selection_status: SelectionStatus
    executable_path: Path | None
    executable_source: ExecutableSource | None
    selection_error: str | None
    runtime_libraries_checked: bool
    runtime_check: str


def diagnose_environment(executable: Path | str | None = None) -> EnvironmentDiagnosticReport:
    """Report executable selection without launching or inspecting the solver."""

    # Keep the package version lookup inside the function so importing this module does not
    # create a cycle while ``kim.__init__`` is establishing its public exports.
    from kim import __version__

    try:
        selection = select_executable(executable)
    except ExecutableError as error:
        return EnvironmentDiagnosticReport(
            schema_version=1,
            python_version=platform.python_version(),
            package_version=__version__,
            selection_status="unavailable",
            executable_path=None,
            executable_source=None,
            selection_error=str(error),
            runtime_libraries_checked=False,
            runtime_check="not performed; executable selection does not verify runtime libraries",
        )
    return _selected_report(selection, __version__)


def _selected_report(
    selection: ExecutableSelection, package_version: str
) -> EnvironmentDiagnosticReport:
    return EnvironmentDiagnosticReport(
        schema_version=1,
        python_version=platform.python_version(),
        package_version=package_version,
        selection_status="selected",
        executable_path=selection.path,
        executable_source=selection.source,  # type: ignore[arg-type]
        selection_error=None,
        runtime_libraries_checked=False,
        runtime_check="not performed; executable selection does not verify runtime libraries",
    )
