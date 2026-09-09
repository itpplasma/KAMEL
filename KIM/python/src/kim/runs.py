"""Versioned manifests and persistent run-directory state management."""

from __future__ import annotations

import hashlib
import json
import os
import re
import uuid
from collections.abc import Callable, Mapping, Sequence
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Literal

from kim.config import SimulationConfig
from kim.errors import RunError
from kim.executable import KamelGitMetadata
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    JsonValue,
    ValidationError,
    field_validator,
    model_validator,
)

MANIFEST_SCHEMA_VERSION = 1
PYTHON_PACKAGE_VERSION = "0.1.0"
_RUN_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_TERMINAL_STATUSES = frozenset({"succeeded", "failed", "timed_out", "interrupted"})


class RunModel(BaseModel):
    """Strict immutable base for persisted run metadata."""

    model_config = ConfigDict(
        allow_inf_nan=False,
        extra="forbid",
        frozen=True,
        str_strip_whitespace=True,
    )


class RunStatus(str, Enum):
    """Persistent lifecycle states for one simulation run."""

    PREPARED = "prepared"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    TIMED_OUT = "timed_out"
    INTERRUPTED = "interrupted"

    @property
    def terminal(self) -> bool:
        return self.value in _TERMINAL_STATUSES


class InputArtifact(RunModel):
    """Source and staged identity of one immutable run input."""

    name: str = Field(min_length=1)
    source_path: Path | None = None
    staged_path: Path
    sha256: str = Field(pattern=_SHA256.pattern)

    @field_validator("staged_path")
    @classmethod
    def staged_path_must_be_relative(cls, path: Path) -> Path:
        return _relative_artifact_path(path, "staged input")

    @field_validator("source_path")
    @classmethod
    def source_path_must_be_absolute(cls, path: Path | None) -> Path | None:
        if path is not None and not path.is_absolute():
            raise ValueError("input source path must be absolute")
        return path


class RunFailure(RunModel):
    """Structured reason for an unsuccessful terminal state."""

    kind: str = Field(min_length=1)
    message: str = Field(min_length=1)
    details: dict[str, JsonValue] = Field(default_factory=dict)


class RunManifest(RunModel):
    """Complete versioned record of one KIM execution request."""

    schema_version: Literal[1] = MANIFEST_SCHEMA_VERSION
    python_package_version: str = PYTHON_PACKAGE_VERSION
    run_id: str = Field(min_length=1)
    label: str | None = None
    parent_sweep_id: str | None = None
    status: RunStatus = RunStatus.PREPARED
    created_at: datetime
    started_at: datetime | None = None
    finished_at: datetime | None = None
    duration_seconds: float | None = Field(default=None, ge=0.0)
    requested_parameters: dict[str, JsonValue]
    normalized_parameters: dict[str, JsonValue]
    executable_path: Path | None = None
    executable_sha256: str | None = Field(default=None, pattern=_SHA256.pattern)
    arguments: tuple[str, ...] = ()
    working_directory: Path | None = None
    exit_code: int | None = None
    kamel_git_root: Path | None = None
    kamel_git_commit: str | None = Field(default=None, pattern=r"^[0-9a-f]{40}$")
    kamel_git_dirty: bool | None = None
    environment: dict[str, str] = Field(default_factory=dict)
    inputs: tuple[InputArtifact, ...] = ()
    expected_outputs: tuple[Path, ...] = ()
    discovered_outputs: tuple[Path, ...] = ()
    validation_details: tuple[str, ...] = ()
    failure: RunFailure | None = None

    @field_validator("run_id", "parent_sweep_id")
    @classmethod
    def identifiers_must_be_safe(cls, value: str | None) -> str | None:
        if value is not None and not _RUN_ID.fullmatch(value):
            raise ValueError("identifier must contain only letters, numbers, '.', '_' and '-'")
        return value

    @field_validator("created_at", "started_at", "finished_at")
    @classmethod
    def timestamps_must_include_timezone(cls, value: datetime | None) -> datetime | None:
        if value is not None and value.utcoffset() is None:
            raise ValueError("manifest timestamps must include a timezone")
        return value

    @field_validator("expected_outputs", "discovered_outputs")
    @classmethod
    def outputs_must_be_relative(cls, paths: tuple[Path, ...]) -> tuple[Path, ...]:
        return tuple(_relative_artifact_path(path, "output") for path in paths)

    @field_validator("executable_path", "working_directory", "kamel_git_root")
    @classmethod
    def provenance_paths_must_be_absolute(cls, path: Path | None) -> Path | None:
        if path is not None and not path.is_absolute():
            raise ValueError("execution provenance paths must be absolute")
        return path

    @field_validator("environment")
    @classmethod
    def environment_is_allowlisted(cls, values: dict[str, str]) -> dict[str, str]:
        unsupported = set(values) - {"OMP_NUM_THREADS"}
        if unsupported:
            raise ValueError(
                "unsupported execution environment field(s): " + ", ".join(sorted(unsupported))
            )
        return values

    @model_validator(mode="after")
    def lifecycle_is_consistent(self) -> RunManifest:
        if self.status is RunStatus.PREPARED:
            if self.started_at is not None or self.finished_at is not None:
                raise ValueError("prepared manifest cannot have execution timestamps")
        elif self.status is RunStatus.RUNNING:
            if self.started_at is None or self.finished_at is not None:
                raise ValueError("running manifest requires only a start timestamp")
        else:
            if self.finished_at is None:
                raise ValueError("terminal manifest requires a finish timestamp")
            if self.status is not RunStatus.INTERRUPTED and self.started_at is None:
                raise ValueError("completed execution requires a start timestamp")
        if self.status is RunStatus.SUCCEEDED:
            if self.exit_code != 0 or self.failure is not None:
                raise ValueError("successful manifest requires exit_code=0 and no failure")
        if (
            self.status
            in {
                RunStatus.FAILED,
                RunStatus.TIMED_OUT,
                RunStatus.INTERRUPTED,
            }
            and self.failure is None
        ):
            raise ValueError("unsuccessful manifest requires failure details")
        if not self.status.terminal and (
            self.exit_code is not None
            or self.failure is not None
            or self.duration_seconds is not None
            or self.discovered_outputs
        ):
            raise ValueError("non-terminal manifest cannot contain completion details")
        return self


class RunPaths(RunModel):
    """Stable paths allocated for one run directory."""

    run_id: str
    root: Path
    manifest: Path
    request: Path
    namelist: Path
    profiles: Path
    logs: Path
    stdout: Path
    stderr: Path
    results: Path


class RunRepository:
    """Create, inspect, and atomically update reproducible KIM runs."""

    def __init__(
        self,
        root: Path | str,
        *,
        clock: Callable[[], datetime] | None = None,
    ) -> None:
        self.root = Path(root).expanduser().absolute()
        self._clock = clock or _utc_now

    def create(
        self,
        config: SimulationConfig,
        *,
        label: str | None = None,
        parent_sweep_id: str | None = None,
        requested_parameters: Mapping[str, Any] | None = None,
    ) -> RunPaths:
        """Allocate a run directory and persist its normalized request."""

        if parent_sweep_id is not None:
            self._validate_identifier(parent_sweep_id, "parent sweep ID")
        created_at = _as_utc(self._clock())
        paths = self._allocate_paths(created_at, label or config.run.run_type.value)
        paths.profiles.mkdir(parents=True)
        paths.logs.mkdir(parents=True)
        paths.results.mkdir(parents=True)
        paths.stdout.touch()
        paths.stderr.touch()

        normalized = config.model_dump(mode="json")
        request = dict(requested_parameters) if requested_parameters is not None else normalized
        request_bytes = (json.dumps(normalized, indent=2, sort_keys=True) + "\n").encode("utf-8")
        paths.request.write_bytes(request_bytes)
        request_artifact = InputArtifact(
            name="request",
            staged_path=paths.request.relative_to(paths.root),
            sha256=hashlib.sha256(request_bytes).hexdigest(),
        )
        try:
            manifest = RunManifest(
                run_id=paths.run_id,
                label=label,
                parent_sweep_id=parent_sweep_id,
                created_at=created_at,
                requested_parameters=request,
                normalized_parameters=normalized,
                inputs=(request_artifact,),
            )
            self._write_manifest(paths.manifest, manifest)
        except (OSError, ValidationError, RunError) as error:
            raise RunError(f"could not create run {paths.run_id}: {error}") from error
        return paths

    def paths(self, run_id: str) -> RunPaths:
        """Return the canonical layout for an existing safe run ID."""

        self._validate_run_id(run_id)
        root = self.root / run_id
        if not root.is_dir():
            raise RunError(f"run not found: {run_id}")
        return _paths_for(root, run_id)

    def inspect(self, run_id: str) -> RunManifest:
        """Load and validate one persisted manifest."""

        paths = self.paths(run_id)
        try:
            manifest = RunManifest.model_validate_json(paths.manifest.read_text(encoding="utf-8"))
        except OSError as error:
            raise RunError(f"could not read manifest for run {run_id}: {error}") from error
        except ValidationError as error:
            raise RunError(f"invalid manifest for run {run_id}: {error}") from error
        if manifest.run_id != run_id:
            raise RunError(
                f"manifest run ID {manifest.run_id!r} does not match directory {run_id!r}"
            )
        return manifest

    def list_runs(self, *, status: RunStatus | None = None) -> tuple[RunManifest, ...]:
        """List valid manifests newest first, optionally filtered by state."""

        if not self.root.exists():
            return ()
        manifests = []
        for directory in self.root.iterdir():
            if not directory.is_dir() or not (directory / "manifest.json").is_file():
                continue
            manifest = self.inspect(directory.name)
            if status is None or manifest.status is status:
                manifests.append(manifest)
        manifests.sort(key=lambda item: (item.created_at, item.run_id), reverse=True)
        return tuple(manifests)

    def configure_execution(
        self,
        run_id: str,
        *,
        executable_path: Path | str,
        executable_sha256: str,
        arguments: Sequence[str],
        working_directory: Path | str,
        kamel_git: KamelGitMetadata | None = None,
        environment: Mapping[str, str] | None = None,
        expected_outputs: Sequence[Path | str] = (),
        validation_details: Sequence[str] = (),
    ) -> RunManifest:
        """Attach finalized execution provenance to a prepared run."""

        manifest = self._mutable_manifest(run_id)
        if manifest.status is not RunStatus.PREPARED:
            raise RunError("execution provenance can only be changed while a run is prepared")
        updates = {
            "executable_path": Path(executable_path).expanduser().absolute(),
            "executable_sha256": executable_sha256,
            "arguments": tuple(arguments),
            "working_directory": Path(working_directory).expanduser().absolute(),
            "kamel_git_root": kamel_git.root if kamel_git else None,
            "kamel_git_commit": kamel_git.commit if kamel_git else None,
            "kamel_git_dirty": kamel_git.dirty if kamel_git else None,
            "environment": dict(environment or {}),
            "expected_outputs": tuple(Path(path) for path in expected_outputs),
            "validation_details": tuple(validation_details),
        }
        return self._replace(run_id, manifest, updates)

    def record_input(
        self,
        run_id: str,
        *,
        name: str,
        source_path: Path | str | None,
        staged_path: Path | str,
    ) -> InputArtifact:
        """Record the digest and provenance of a staged input file."""

        manifest = self._mutable_manifest(run_id)
        if manifest.status is not RunStatus.PREPARED:
            raise RunError("inputs can only be changed while a run is prepared")
        paths = self.paths(run_id)
        staged = Path(staged_path).resolve()
        try:
            relative = staged.relative_to(paths.root.resolve())
        except ValueError as error:
            raise RunError(
                f"staged input must be inside run directory {paths.root}: {staged}"
            ) from error
        if not staged.is_file():
            raise RunError(f"staged input file not found: {staged}")
        if any(item.name == name or item.staged_path == relative for item in manifest.inputs):
            raise RunError(f"input artifact is already recorded: {name} ({relative})")
        artifact = InputArtifact(
            name=name,
            source_path=(
                Path(source_path).expanduser().absolute() if source_path is not None else None
            ),
            staged_path=relative,
            sha256=_sha256(staged),
        )
        self._replace(run_id, manifest, {"inputs": manifest.inputs + (artifact,)})
        return artifact

    def transition(
        self,
        run_id: str,
        status: RunStatus,
        *,
        exit_code: int | None = None,
        discovered_outputs: Sequence[Path | str] = (),
        failure: RunFailure | None = None,
    ) -> RunManifest:
        """Atomically apply one valid lifecycle transition."""

        manifest = self._mutable_manifest(run_id)
        allowed = {
            RunStatus.PREPARED: {RunStatus.RUNNING, RunStatus.INTERRUPTED},
            RunStatus.RUNNING: {
                RunStatus.SUCCEEDED,
                RunStatus.FAILED,
                RunStatus.TIMED_OUT,
                RunStatus.INTERRUPTED,
            },
        }
        if status not in allowed.get(manifest.status, set()):
            raise RunError(f"invalid run transition: {manifest.status.value} -> {status.value}")

        now = _as_utc(self._clock())
        updates: dict[str, Any] = {"status": status}
        if status is RunStatus.RUNNING:
            updates["started_at"] = now
        else:
            updates["finished_at"] = now
            updates["exit_code"] = 0 if status is RunStatus.SUCCEEDED else exit_code
            updates["failure"] = failure
            updates["discovered_outputs"] = tuple(Path(path) for path in discovered_outputs)
            if manifest.started_at is not None:
                updates["duration_seconds"] = max(0.0, (now - manifest.started_at).total_seconds())
        return self._replace(run_id, manifest, updates)

    def _mutable_manifest(self, run_id: str) -> RunManifest:
        manifest = self.inspect(run_id)
        if manifest.status.terminal:
            raise RunError(f"completed run {run_id} is immutable ({manifest.status.value})")
        return manifest

    def _replace(
        self, run_id: str, manifest: RunManifest, updates: Mapping[str, Any]
    ) -> RunManifest:
        try:
            values = manifest.model_dump()
            values.update(updates)
            replacement = RunManifest.model_validate(values)
        except ValidationError as error:
            raise RunError(f"invalid manifest update for run {run_id}: {error}") from error
        self._write_manifest(self.paths(run_id).manifest, replacement)
        return replacement

    def _allocate_paths(self, created_at: datetime, label: str) -> RunPaths:
        self.root.mkdir(parents=True, exist_ok=True)
        slug = _slug(label)
        timestamp = created_at.strftime("%Y%m%dT%H%M%SZ")
        for _attempt in range(100):
            run_id = f"{timestamp}-{slug}-{uuid.uuid4().hex[:8]}"
            root = self.root / run_id
            try:
                root.mkdir()
            except FileExistsError:
                continue
            return _paths_for(root, run_id)
        raise RunError("could not allocate a unique run ID after 100 attempts")

    @staticmethod
    def _write_manifest(path: Path, manifest: RunManifest) -> None:
        temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
        try:
            with temporary.open("w", encoding="utf-8") as stream:
                stream.write(manifest.model_dump_json(indent=2))
                stream.write("\n")
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, path)
        except OSError as error:
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                pass
            raise RunError(f"atomic manifest write failed for {path}: {error}") from error

    @staticmethod
    def _validate_run_id(run_id: str) -> None:
        RunRepository._validate_identifier(run_id, "run ID")

    @staticmethod
    def _validate_identifier(value: str, label: str) -> None:
        if not _RUN_ID.fullmatch(value):
            raise RunError(f"invalid {label}: {value!r}")


def _paths_for(root: Path, run_id: str) -> RunPaths:
    inputs = root / "inputs"
    logs = root / "logs"
    return RunPaths(
        run_id=run_id,
        root=root,
        manifest=root / "manifest.json",
        request=inputs / "request.json",
        namelist=inputs / "KIM_config.nml",
        profiles=inputs / "profiles",
        logs=logs,
        stdout=logs / "stdout.log",
        stderr=logs / "stderr.log",
        results=root / "results",
    )


def _relative_artifact_path(path: Path, kind: str) -> Path:
    if path.is_absolute() or ".." in path.parts or path == Path("."):
        raise ValueError(f"{kind} path must be a non-empty path relative to the run directory")
    return path


def _slug(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", value.casefold()).strip("-")
    return slug[:48] or "run"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _as_utc(value: datetime) -> datetime:
    if value.utcoffset() is None:
        raise RunError("repository clock must return a timezone-aware datetime")
    return value.astimezone(timezone.utc)
