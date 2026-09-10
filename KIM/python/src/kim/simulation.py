"""Preparation and subprocess orchestration for one KIM simulation."""

from __future__ import annotations

import os
import subprocess
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import h5py
from kim.config import RunType, SimulationConfig
from kim.errors import KimError
from kim.executable import discover_kamel_git_metadata, executable_sha256, resolve_executable
from kim.namelist import write_namelist
from kim.profiles import ProfileSet
from kim.runs import RunFailure, RunManifest, RunPaths, RunRepository, RunStatus

if TYPE_CHECKING:
    from kim.results import DatasetMetadata, PeriodicResult, Result

_OUTPUT_FILES = {
    RunType.ELECTROSTATIC_PERIODIC: "out_ES_periodic.h5",
    RunType.ELECTROSTATIC: "out_ES.h5",
    RunType.FLR2: "out_FLR2.h5",
}
_REQUIRED_DATASETS = {
    RunType.ELECTROSTATIC_PERIODIC: (
        "fields/Phi",
        "fields/jpar",
        "backs/e/r",
        "setup/periodic_scale/dx_asis",
    ),
    RunType.ELECTROSTATIC: ("fields/Phi_m", "fields/jpar", "backs/e/r"),
    RunType.FLR2: ("fields/Phi", "fields/jpar", "backs/e/r"),
}
_ALLOWED_ENVIRONMENT = frozenset({"OMP_NUM_THREADS"})


@dataclass(frozen=True)
class PreparedSimulation:
    """A fully staged simulation that has not yet started."""

    paths: RunPaths
    manifest: RunManifest
    command: tuple[str, str]
    output_file: Path


@dataclass(frozen=True)
class RunResult:
    """Terminal identity, artifacts, and manifest for one execution attempt."""

    paths: RunPaths
    manifest: RunManifest
    output_file: Path

    @property
    def run_id(self) -> str:
        return self.paths.run_id

    @property
    def status(self) -> RunStatus:
        return self.manifest.status

    @property
    def run_directory(self) -> Path:
        return self.paths.root

    @property
    def namelist(self) -> Path:
        return self.paths.namelist

    @property
    def stdout(self) -> Path:
        return self.paths.stdout

    @property
    def stderr(self) -> Path:
        return self.paths.stderr

    def list_datasets(self) -> tuple[str, ...]:
        """List datasets in this run's primary HDF5 output."""

        return self.output.list_datasets()

    def read_dataset(self, path: str) -> object:
        """Read one dataset from this run's primary HDF5 output."""

        return self.output.read_dataset(path)

    def dataset_metadata(self, path: str) -> DatasetMetadata:
        """Read metadata for one dataset in the primary HDF5 output."""

        return self.output.dataset_metadata(path)

    @property
    def output(self) -> Result:
        """Return generic structured access to the primary HDF5 output."""

        from kim.results import Result

        return Result(self.output_file)

    @property
    def periodic(self) -> PeriodicResult:
        """Return the validated periodic view of the primary output."""

        return self.output.periodic


class Simulation:
    """Validate, stage, and run one immutable KIM request."""

    def __init__(
        self,
        config: SimulationConfig,
        *,
        executable: Path | str | None = None,
        runs_directory: Path | str = "runs",
        label: str | None = None,
        timeout: float | None = None,
        parent_sweep_id: str | None = None,
        requested_parameters: Mapping[str, object] | None = None,
    ) -> None:
        if timeout is not None and timeout <= 0:
            raise ValueError("timeout must be positive or None")
        self.config = config
        self.executable = executable
        self.repository = RunRepository(runs_directory)
        self.label = label
        self.timeout = timeout
        self.parent_sweep_id = parent_sweep_id
        self.requested_parameters = requested_parameters
        self._active_paths: RunPaths | None = None
        self._last_result: RunResult | None = None

    def prepare(
        self,
        *,
        environment: Mapping[str, str] | None = None,
    ) -> PreparedSimulation:
        """Validate all external inputs and create a reproducible prepared run."""

        execution_environment = _validate_environment(environment)
        executable = resolve_executable(self.executable)
        profiles = ProfileSet.from_simulation(self.config)
        validation = profiles.validate_for(self.config)

        paths = self.repository.create(
            self.config,
            label=self.label,
            parent_sweep_id=self.parent_sweep_id,
            requested_parameters=self.requested_parameters,
        )
        self._active_paths = paths
        copies = profiles.copy_to(paths.profiles)
        for copied in copies:
            self.repository.record_input(
                paths.run_id,
                name=copied.name,
                source_path=copied.source,
                staged_path=copied.destination,
            )

        run_type = self.config.run.run_type
        output_name = _OUTPUT_FILES[run_type]
        write_namelist(
            self.config,
            paths.namelist,
            profile_directory=paths.profiles.absolute(),
            output_directory=paths.results.absolute(),
            output_file=output_name,
        )
        self.repository.record_input(
            paths.run_id,
            name="namelist",
            source_path=None,
            staged_path=paths.namelist,
        )

        output_file = (
            paths.results / f"m{self.config.setup.m_mode}_n{self.config.setup.n_mode}" / output_name
        )
        expected = output_file.relative_to(paths.root)
        command = (str(executable), str(paths.namelist.absolute()))
        manifest = self.repository.configure_execution(
            paths.run_id,
            executable_path=executable,
            executable_sha256=executable_sha256(executable),
            arguments=command[1:],
            working_directory=paths.root,
            kamel_git=discover_kamel_git_metadata(executable),
            environment=execution_environment,
            expected_outputs=(expected,),
            validation_details=(
                f"validated {validation.point_count} profile points",
                f"profile radius [{validation.radial_minimum:g}, {validation.radial_maximum:g}] cm",
                f"q resonance at r_eff={validation.resonance_radius:g} cm",
            ),
        )
        return PreparedSimulation(
            paths=paths,
            manifest=manifest,
            command=command,
            output_file=output_file,
        )

    def run(
        self,
        *,
        timeout: float | None = None,
        environment: Mapping[str, str] | None = None,
    ) -> RunResult:
        """Execute KIM and return a persisted terminal result."""

        effective_timeout = self.timeout if timeout is None else timeout
        if effective_timeout is not None and effective_timeout <= 0:
            raise ValueError("timeout must be positive or None")
        self._active_paths = None
        self._last_result = None
        try:
            prepared = self.prepare(environment=environment)
        except KeyboardInterrupt:
            self._record_preparation_interruption()
            raise
        except (KimError, OSError, UnicodeError) as error:
            return self._record_preparation_failure(error)
        self.repository.transition(prepared.paths.run_id, RunStatus.RUNNING)
        process_environment = os.environ.copy()
        process_environment.update(prepared.manifest.environment)

        try:
            with (
                prepared.paths.stdout.open("wb") as stdout,
                prepared.paths.stderr.open("wb") as stderr,
            ):
                completed = subprocess.run(
                    list(prepared.command),
                    cwd=prepared.paths.root,
                    env=process_environment,
                    stdout=stdout,
                    stderr=stderr,
                    timeout=effective_timeout,
                    shell=False,
                    check=False,
                )
        except subprocess.TimeoutExpired:
            manifest = self.repository.transition(
                prepared.paths.run_id,
                RunStatus.TIMED_OUT,
                failure=RunFailure(
                    kind="timeout",
                    message=f"KIM exceeded the {effective_timeout:g} second timeout",
                    details={"timeout_seconds": effective_timeout},
                ),
            )
        except KeyboardInterrupt:
            manifest = self.repository.transition(
                prepared.paths.run_id,
                RunStatus.INTERRUPTED,
                failure=RunFailure(kind="interrupted", message="KIM execution was interrupted"),
            )
            self._last_result = RunResult(
                paths=prepared.paths,
                manifest=manifest,
                output_file=prepared.output_file,
            )
            raise
        except OSError as error:
            manifest = self.repository.transition(
                prepared.paths.run_id,
                RunStatus.FAILED,
                failure=RunFailure(
                    kind="execution_error",
                    message=f"could not execute KIM: {error}",
                    details={"exception_type": type(error).__name__},
                ),
            )
        else:
            manifest = self._complete(prepared, completed.returncode)

        result = RunResult(
            paths=prepared.paths, manifest=manifest, output_file=prepared.output_file
        )
        self._last_result = result
        return result

    def _record_preparation_failure(self, error: Exception) -> RunResult:
        """Persist one failed preparation attempt for sweep-level accounting."""

        paths = self._active_paths
        if paths is None:
            paths = self.repository.create(
                self.config,
                label=self.label,
                parent_sweep_id=self.parent_sweep_id,
                requested_parameters=self.requested_parameters,
            )
            self._active_paths = paths
        output_file = (
            paths.results
            / f"m{self.config.setup.m_mode}_n{self.config.setup.n_mode}"
            / _OUTPUT_FILES[self.config.run.run_type]
        )
        manifest = self.repository.inspect(paths.run_id)
        if manifest.status is RunStatus.PREPARED:
            manifest = self.repository.transition(
                paths.run_id,
                RunStatus.FAILED,
                failure=RunFailure(
                    kind="preparation_error",
                    message=str(error),
                    details={"exception_type": type(error).__name__},
                ),
            )
        result = RunResult(paths=paths, manifest=manifest, output_file=output_file)
        self._last_result = result
        return result

    def _record_preparation_interruption(self) -> RunResult:
        """Persist an interruption that occurred before process execution."""

        paths = self._active_paths
        if paths is None:
            paths = self.repository.create(
                self.config,
                label=self.label,
                parent_sweep_id=self.parent_sweep_id,
                requested_parameters=self.requested_parameters,
            )
            self._active_paths = paths
        output_file = (
            paths.results
            / f"m{self.config.setup.m_mode}_n{self.config.setup.n_mode}"
            / _OUTPUT_FILES[self.config.run.run_type]
        )
        manifest = self.repository.transition(
            paths.run_id,
            RunStatus.INTERRUPTED,
            failure=RunFailure(kind="interrupted", message="KIM preparation was interrupted"),
        )
        result = RunResult(paths=paths, manifest=manifest, output_file=output_file)
        self._last_result = result
        return result

    def _complete(self, prepared: PreparedSimulation, exit_code: int) -> RunManifest:
        discovered = _discovered_outputs(prepared)
        if exit_code != 0:
            return self.repository.transition(
                prepared.paths.run_id,
                RunStatus.FAILED,
                exit_code=exit_code,
                discovered_outputs=discovered,
                failure=RunFailure(
                    kind="process_exit",
                    message=f"KIM exited with status {exit_code}",
                    details={"exit_code": exit_code},
                ),
            )
        failure = _validate_output(prepared.output_file, self.config.run.run_type)
        if failure is not None:
            return self.repository.transition(
                prepared.paths.run_id,
                RunStatus.FAILED,
                exit_code=0,
                discovered_outputs=discovered,
                failure=failure,
            )
        return self.repository.transition(
            prepared.paths.run_id,
            RunStatus.SUCCEEDED,
            exit_code=0,
            discovered_outputs=discovered,
        )


def _validate_environment(environment: Mapping[str, str] | None) -> dict[str, str]:
    values = {name: os.environ[name] for name in _ALLOWED_ENVIRONMENT if name in os.environ}
    values.update({str(name): str(value) for name, value in (environment or {}).items()})
    unsupported = set(values) - _ALLOWED_ENVIRONMENT
    if unsupported:
        raise ValueError(
            "unsupported execution environment field(s): " + ", ".join(sorted(unsupported))
        )
    return values


def _discovered_outputs(prepared: PreparedSimulation) -> tuple[Path, ...]:
    if not prepared.output_file.is_file():
        return ()
    return (prepared.output_file.relative_to(prepared.paths.root),)


def _validate_output(path: Path, run_type: RunType) -> RunFailure | None:
    if not path.is_file():
        return RunFailure(
            kind="missing_output",
            message=f"KIM exited successfully but did not create {path.name}",
            details={"expected_output": str(path)},
        )
    try:
        with h5py.File(path, "r") as handle:
            missing = [name for name in _REQUIRED_DATASETS[run_type] if name not in handle]
            invalid = []
            for name in _REQUIRED_DATASETS[run_type]:
                if name in handle:
                    value = handle[name]
                    if not isinstance(value, h5py.Dataset):
                        invalid.append(name)
                    else:
                        value[()]
    except (KeyError, OSError, TypeError, ValueError) as error:
        return RunFailure(
            kind="invalid_output",
            message=f"KIM output is not a readable HDF5 file: {error}",
            details={"output": str(path)},
        )
    if missing:
        return RunFailure(
            kind="invalid_output",
            message="KIM output is missing required datasets: " + ", ".join(missing),
            details={"output": str(path), "missing_datasets": missing},
        )
    if invalid:
        return RunFailure(
            kind="invalid_output",
            message="KIM output contains invalid required datasets: " + ", ".join(invalid),
            details={"output": str(path), "invalid_datasets": invalid},
        )
    if run_type is RunType.ELECTROSTATIC_PERIODIC:
        from kim.errors import ResultError
        from kim.results import Result

        try:
            Result(path).periodic
        except ResultError as error:
            return RunFailure(
                kind="invalid_output",
                message=f"KIM periodic output violates the result contract: {error}",
                details={"output": str(path)},
            )
    return None
