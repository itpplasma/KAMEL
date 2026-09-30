"""Validated KIM Condor scan plans and immutable per-point staging."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import time
import fcntl
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

import f90nml
import numpy as np
from kim.config import ElectrostaticPeriodicRun, ProfileConfig, SimulationConfig
from kim.condor_client import (
    CondorError,
    CondorSubmitSpec,
    CondorToolConfig,
    parse_condor_submit_output,
    query_job_ads,
    remove_clusters,
    run_condor,
    verify_shared_filesystem,
)
from kim.errors import ProfileError
from kim.profiles import ProfileSet, _read_profile
from kim.sweep import SweepSpec
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

MANIFEST_NAME = "condor_jobs.json"
STATUS_NAME = "condor_status.json"
WORKER_NAME = "condor_worker.py"
_JOB_STATUS_NAMES = {
    1: "idle",
    2: "running",
    3: "removed",
    4: "completed",
    5: "held",
    6: "transferring_output",
    7: "suspended",
}
_TERMINAL_CODES = {3, 4, 5}


class _CondorModel(BaseModel):
    model_config = ConfigDict(
        allow_inf_nan=False,
        arbitrary_types_allowed=True,
        extra="forbid",
        frozen=True,
        str_strip_whitespace=True,
    )


class JparCurrentMetric(_CondorModel):
    """Explicit interpretation of one current column in a direct solver output."""

    current_column: int = Field(ge=1)
    current_unit: str = Field(min_length=1)
    collision_model: str = Field(min_length=1)

    @field_validator("collision_model")
    @classmethod
    def collision_model_is_filename_token(cls, value: str) -> str:
        if value in {".", ".."} or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
            raise ValueError("collision_model must be a safe filename token")
        return value


class KimCondorPlan(_CondorModel):
    """Explicit backend, filesystem, resource, and metric choices for a KIM scan."""

    backend: Literal["kamel_kim_python", "kim_x_namelist"]
    executable: Path
    python_executable: Path
    kim_source_path: Path | None = None
    condor: CondorToolConfig = Field(default_factory=CondorToolConfig)
    shared_filesystem_prefixes: tuple[Path, ...] = ()
    should_transfer_files: Literal["NEVER", "IF_NEEDED", "ALWAYS"] = "NEVER"
    request_cpus: int = Field(default=1, ge=1)
    omp_threads_per_process: int = Field(default=1, ge=1)
    request_memory_mb: int = Field(default=8192, ge=1)
    job_timeout_s: float = Field(default=3600.0, gt=0.0)
    poll_interval_s: float = Field(default=30.0, gt=0.0)
    max_wall_clock_s: float = Field(default=21600.0, gt=0.0)
    exclude_machines: tuple[str, ...] = ()
    base_namelist: Path | None = None
    namelist_overrides: dict[str, Any] = Field(default_factory=dict)
    jpar_current_metric: JparCurrentMetric | None = None
    api_current_component: Literal["real", "imag", "magnitude"] | None = None
    api_current_unit: str | None = None

    @model_validator(mode="after")
    def validate_plan(self) -> KimCondorPlan:
        for name in ("executable", "python_executable"):
            path = getattr(self, name)
            if not path.is_absolute():
                raise ValueError(f"{name} must be an absolute path")
        if self.kim_source_path is not None and not self.kim_source_path.is_absolute():
            raise ValueError("kim_source_path must be an absolute path")
        if self.base_namelist is not None and not self.base_namelist.is_absolute():
            raise ValueError("base_namelist must be an absolute path")
        if self.omp_threads_per_process > self.request_cpus:
            raise ValueError("omp_threads_per_process cannot exceed request_cpus")
        if not self.shared_filesystem_prefixes:
            raise ValueError("shared_filesystem_prefixes must contain at least one absolute path")
        if any(not prefix.is_absolute() for prefix in self.shared_filesystem_prefixes):
            raise ValueError("shared_filesystem_prefixes must contain only absolute paths")
        if any(not machine.strip() for machine in self.exclude_machines):
            raise ValueError("exclude_machines must contain non-empty machine names")
        if self.backend == "kim_x_namelist":
            if self.base_namelist is None:
                raise ValueError("kim_x_namelist backend requires a reviewed base_namelist")
            if self.api_current_component is not None or self.api_current_unit is not None:
                raise ValueError("API current settings are only valid for kamel_kim_python")
        elif self.base_namelist is not None or self.namelist_overrides:
            raise ValueError(
                "base_namelist and namelist_overrides are only valid for kim_x_namelist"
            )
        if self.jpar_current_metric is not None and self.backend != "kim_x_namelist":
            raise ValueError("jpar_current_metric is only valid for kim_x_namelist")
        if (self.api_current_component is None) != (self.api_current_unit is None):
            raise ValueError("api_current_component and api_current_unit must be supplied together")
        if self.api_current_unit is not None and not self.api_current_unit.strip():
            raise ValueError("api_current_unit must be non-empty")
        return self


@dataclass(frozen=True)
class KimCondorJob:
    """One deterministic staged scan point."""

    job: str
    directory: Path
    profile_scale_factor: float
    scan_order_index: int


@dataclass(frozen=True)
class KimCondorScanResults:
    """Ordered per-scale outcomes and an explicitly gated scalar resonance result."""

    curve: tuple[dict[str, object], ...]
    resonance: dict[str, object]
    notes: tuple[str, ...]


def stage_condor_sweep(
    root: str | Path,
    *,
    spec: SweepSpec,
    plan: KimCondorPlan,
) -> tuple[KimCondorJob, ...]:
    """Validate and stage one isolated job directory for each Er scale factor."""

    if not isinstance(spec, SweepSpec):
        raise CondorError("spec must be a validated SweepSpec")
    run_root = Path(root).expanduser().absolute()
    if run_root.is_symlink():
        raise CondorError(f"refusing to stage into a symlink run directory: {run_root}")
    try:
        verify_shared_filesystem(
            run_root,
            allowed_prefixes=plan.shared_filesystem_prefixes,
            should_transfer_files=plan.should_transfer_files,
        )
    except CondorError:
        raise

    _validate_runtime_paths(plan)
    profile_set = ProfileSet.from_simulation(spec.base)
    try:
        profile_set.validate_for(spec.base)
        profiles = profile_set._read_existing()
    except ProfileError as error:
        raise CondorError(f"KIM source profiles failed validation: {error}") from error
    if "Er" not in profiles:
        raise CondorError("KIM profile scaling requires a source Er.dat profile")
    if not np.array_equal(profiles["Er"].radius, profiles["n"].radius):
        raise CondorError("Er.dat radial grid must match the density profile grid")

    reviewed_groups: dict[str, dict[str, Any]] | None = None
    if plan.backend == "kim_x_namelist":
        reviewed_groups = _read_reviewed_namelist(plan)
        _validate_namelist_overrides(reviewed_groups, plan.namelist_overrides)

    factors = tuple(float(factor) for factor in spec.variation.values)
    if not factors or any(not math.isfinite(factor) for factor in factors):
        raise CondorError("SweepSpec must contain finite profile scale factors")
    if plan.api_current_component is not None and not isinstance(
        spec.base.run, ElectrostaticPeriodicRun
    ):
        raise CondorError(
            "API integrated-current metrics require an electrostatic-periodic SweepSpec"
        )

    jobs = tuple(
        KimCondorJob(
            job=f"scale-{index:04d}",
            directory=run_root / f"scale-{index:04d}",
            profile_scale_factor=factor,
            scan_order_index=index,
        )
        for index, factor in enumerate(factors)
    )
    if run_root.exists() and not run_root.is_dir():
        raise CondorError(f"run root is not a directory: {run_root}")
    for job in jobs:
        if job.directory.exists() or job.directory.is_symlink():
            raise CondorError(f"refusing to overwrite existing job directory: {job.directory}")

    run_root.mkdir(parents=True, exist_ok=True)
    for job in jobs:
        _stage_job(
            job, spec=spec, plan=plan, profile_set=profile_set, reviewed_groups=reviewed_groups
        )
    return jobs


def _validate_runtime_paths(plan: KimCondorPlan) -> None:
    if not plan.executable.is_file():
        raise CondorError(f"KIM executable is missing: {plan.executable}")
    if not plan.python_executable.is_file():
        raise CondorError(f"Python executable is missing: {plan.python_executable}")
    if plan.kim_source_path is not None and not plan.kim_source_path.is_dir():
        raise CondorError(f"KIM source path is not a directory: {plan.kim_source_path}")


def _read_reviewed_namelist(plan: KimCondorPlan) -> dict[str, dict[str, Any]]:
    path = plan.base_namelist
    if path is None or not path.is_file():
        raise CondorError(f"reviewed base namelist is missing: {path}")
    try:
        payload = f90nml.read(path).todict()
    except (OSError, ValueError, StopIteration) as error:
        raise CondorError(f"could not read reviewed base namelist {path}: {error}") from error
    groups: dict[str, dict[str, Any]] = {}
    for name, values in payload.items():
        if not isinstance(values, dict):
            raise CondorError(f"namelist group {name!r} is not a group mapping")
        groups[str(name).lower()] = {str(key).lower(): value for key, value in values.items()}
    for group, key in (
        ("kim_io", "profile_location"),
        ("kim_io", "output_path"),
        ("kim_profiles", "input_profile_dir"),
    ):
        if group not in groups or key not in groups[group]:
            raise CondorError(f"reviewed base namelist is missing required {group}.{key}")
    return groups


def _validate_namelist_overrides(
    groups: dict[str, dict[str, Any]], overrides: dict[str, Any]
) -> None:
    for dotted_key, _value in overrides.items():
        parts = dotted_key.split(".") if isinstance(dotted_key, str) else []
        if len(parts) != 2 or not all(parts):
            raise CondorError(f"namelist override must use a group.key path: {dotted_key!r}")
        group, key = (part.lower() for part in parts)
        if group not in groups or key not in groups[group]:
            raise CondorError(
                f"namelist override {dotted_key!r} is not present in reviewed base namelist"
            )


def _stage_job(
    job: KimCondorJob,
    *,
    spec: SweepSpec,
    plan: KimCondorPlan,
    profile_set: ProfileSet,
    reviewed_groups: dict[str, dict[str, Any]] | None,
) -> None:
    job.directory.mkdir(parents=False, exist_ok=False)
    profile_directory = job.directory / "profiles"
    copies = profile_set.copy_to(profile_directory)
    hashes = {copy.destination.name: copy.sha256 for copy in copies}

    if plan.backend == "kim_x_namelist":
        er_path = profile_directory / spec.base.profiles.radial_electric_field_file
        try:
            er_profile = _read_profile("Er", "statV/cm", er_path)
        except ProfileError as error:
            raise CondorError(f"staged Er profile failed validation: {error}") from error
        scaled = np.column_stack((er_profile.radius, er_profile.values * job.profile_scale_factor))
        np.savetxt(er_path, scaled, fmt="%.16e")
        hashes[er_path.name] = _sha256(er_path)

    staged_config = _config_for_job(spec.base, profile_directory)
    mode = "kim-sweep" if plan.backend == "kamel_kim_python" else "kim-run"
    if plan.backend == "kim_x_namelist":
        assert reviewed_groups is not None
        _write_job_namelist(job, profile_directory, plan, reviewed_groups)
        namelist_path = job.directory / "KIM_config.nml"
        hashes[namelist_path.name] = _sha256(namelist_path)

    payload: dict[str, object] = {
        "schema_version": 1,
        "mode": mode,
        "backend": plan.backend,
        "job_identity": "",
        "job": job.job,
        "scan_order_index": job.scan_order_index,
        "profile_scale_factor": job.profile_scale_factor,
        "base_config": staged_config.model_dump(mode="json"),
        "executable": str(plan.executable),
        "executable_sha256": _sha256(plan.executable),
        "python_executable": str(plan.python_executable),
        "kim_source_path": str(plan.kim_source_path) if plan.kim_source_path else None,
        "omp_num_threads": plan.omp_threads_per_process,
        "request_cpus": plan.request_cpus,
        "request_memory_mb": plan.request_memory_mb,
        "job_timeout_s": plan.job_timeout_s,
        "poll_interval_s": plan.poll_interval_s,
        "max_wall_clock_s": plan.max_wall_clock_s,
        "exclude_machines": list(plan.exclude_machines),
        "shared_filesystem_prefixes": [str(path) for path in plan.shared_filesystem_prefixes],
        "should_transfer_files": plan.should_transfer_files,
        "api_current_component": plan.api_current_component,
        "api_current_unit": plan.api_current_unit,
        "jpar_current_metric": (
            plan.jpar_current_metric.model_dump(mode="json")
            if plan.jpar_current_metric is not None
            else None
        ),
        "namelist_file": "KIM_config.nml" if mode == "kim-run" else None,
        "profile_directory": "profiles",
        "output_directory": ".",
        "input_sha256": dict(sorted(hashes.items())),
    }
    identity_data = dict(payload)
    identity_data.pop("job_identity")
    payload["job_identity"] = hashlib.sha256(
        json.dumps(identity_data, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    _write_exclusive_json(job.directory / "condor_job.json", payload)


def _config_for_job(config: SimulationConfig, profile_directory: Path) -> SimulationConfig:
    profiles = ProfileConfig.model_validate(
        {**config.profiles.model_dump(mode="python"), "directory": profile_directory}
    )
    return SimulationConfig.model_validate(
        {**config.model_dump(mode="python"), "profiles": profiles}
    )


def _write_job_namelist(
    job: KimCondorJob,
    profile_directory: Path,
    plan: KimCondorPlan,
    reviewed_groups: dict[str, dict[str, Any]],
) -> None:
    groups = {name: dict(values) for name, values in reviewed_groups.items()}
    for dotted_key, value in plan.namelist_overrides.items():
        group, key = (part.lower() for part in dotted_key.split("."))
        groups[group][key] = value
    profile_path = profile_directory.as_posix().rstrip("/") + "/"
    groups["kim_io"]["profile_location"] = profile_path
    groups["kim_io"]["output_path"] = job.directory.as_posix()
    groups["kim_profiles"]["input_profile_dir"] = profile_path
    destination = job.directory / "KIM_config.nml"
    try:
        f90nml.write(f90nml.Namelist(groups), destination, force=False)
    except (OSError, ValueError, TypeError) as error:
        raise CondorError(f"could not write staged KIM namelist {destination}: {error}") from error


def _write_exclusive_json(path: Path, payload: dict[str, object]) -> None:
    try:
        with path.open("x", encoding="utf-8") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
    except FileExistsError as error:
        raise CondorError(f"refusing to overwrite staged job input: {path}") from error
    except (OSError, TypeError, ValueError) as error:
        raise CondorError(f"could not write staged job input {path}: {error}") from error


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


@contextmanager
def _queue_lock(root: Path):
    """Serialize submit/reconcile operations for one staged scan directory."""

    lock_path = root / ".kim-condor.lock"
    flags = os.O_CREAT | os.O_RDWR | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(str(lock_path), flags, 0o600)
    except OSError as error:
        raise CondorError(f"could not open KIM Condor lock {lock_path}: {error}") from error
    with os.fdopen(descriptor, "r+b") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX)
        except OSError as error:
            raise CondorError(f"could not acquire KIM Condor lock {lock_path}: {error}") from error
        try:
            yield
        finally:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def _write_json_atomic(path: Path, payload: dict[str, object]) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        directory_fd = os.open(str(path.parent), os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        if temporary.exists():
            temporary.unlink()


def _read_staged_jobs(root: Path) -> list[dict[str, object]]:
    if root.is_symlink() or not root.is_dir():
        raise CondorError(f"KIM Condor run directory is missing or unsafe: {root}")
    directories = sorted(
        (path for path in root.iterdir() if path.is_dir() and path.name.startswith("scale-")),
        key=lambda path: path.name,
    )
    if not directories:
        raise CondorError(f"no staged KIM Condor jobs found under {root}")
    jobs = []
    for directory in directories:
        if not re.fullmatch(r"scale-\d{4,}", directory.name):
            raise CondorError(f"unsafe staged KIM job directory name: {directory.name}")
        input_path = directory / "condor_job.json"
        try:
            payload = json.loads(input_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as error:
            raise CondorError(
                f"could not read staged KIM job input {input_path}: {error}"
            ) from error
        if not isinstance(payload, dict) or payload.get("schema_version") != 1:
            raise CondorError(f"unsupported staged KIM job input: {input_path}")
        if payload.get("job") != directory.name:
            raise CondorError(f"staged job name does not match directory: {directory}")
        _verify_staged_identity(directory, payload)
        jobs.append(payload)
    indexes = [job.get("scan_order_index") for job in jobs]
    if indexes != list(range(len(jobs))):
        raise CondorError("staged KIM jobs must contain contiguous scan indexes in order")
    return jobs


def _verify_staged_identity(directory: Path, payload: dict[str, object]) -> None:
    identity = payload.get("job_identity")
    if not isinstance(identity, str) or not re.fullmatch(r"[0-9a-f]{64}", identity):
        raise CondorError(f"staged KIM job has an invalid identity: {directory}")
    raw_hashes = payload.get("input_sha256")
    if not isinstance(raw_hashes, dict):
        raise CondorError(f"staged KIM job has no input checksums: {directory}")
    for name, expected_hash in raw_hashes.items():
        if (
            not isinstance(name, str)
            or Path(name).name != name
            or not isinstance(expected_hash, str)
            or not re.fullmatch(r"[0-9a-f]{64}", expected_hash)
        ):
            raise CondorError(f"staged KIM job has an invalid input checksum: {directory}")
        source = directory / name if name == "KIM_config.nml" else directory / "profiles" / name
        try:
            actual_hash = _sha256(source)
        except OSError as error:
            raise CondorError(f"staged KIM input is missing: {source}") from error
        if actual_hash != expected_hash:
            raise CondorError(f"staged KIM input checksum mismatch: {source}")
    identity_payload = dict(payload)
    identity_payload.pop("job_identity", None)
    actual_identity = hashlib.sha256(
        json.dumps(identity_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    if actual_identity != identity:
        raise CondorError(f"staged KIM job identity does not match its payload: {directory}")
    executable = payload.get("executable")
    executable_hash = payload.get("executable_sha256")
    if not isinstance(executable, str) or not Path(executable).is_absolute():
        raise CondorError(f"staged KIM job has an invalid executable path: {directory}")
    try:
        actual_executable_hash = _sha256(Path(executable))
    except OSError as error:
        raise CondorError(f"KIM executable is unavailable: {executable}") from error
    if actual_executable_hash != executable_hash:
        raise CondorError(f"KIM executable checksum mismatch: {executable}")


def _new_manifest(root: Path, jobs: list[dict[str, object]]) -> dict[str, object]:
    return {
        "schema_version": 1,
        "run_directory": str(root),
        "created_at_utc": _utc_now(),
        "updated_at_utc": _utc_now(),
        "state": "prepared",
        "jobs": [
            {
                "job": job["job"],
                "directory": job["job"],
                "job_identity": job["job_identity"],
                "profile_scale_factor": job["profile_scale_factor"],
                "scan_order_index": job["scan_order_index"],
                "state": "prepared",
                "cluster": None,
                "error": None,
                "queue_status": None,
            }
            for job in jobs
        ],
    }


def _load_manifest(root: Path, jobs: list[dict[str, object]]) -> dict[str, object]:
    path = root / MANIFEST_NAME
    if path.is_symlink():
        raise CondorError(f"refusing to read KIM Condor manifest symlink: {path}")
    if not path.exists():
        return _new_manifest(root, jobs)
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CondorError(f"could not read KIM Condor manifest {path}: {error}") from error
    if (
        not isinstance(manifest, dict)
        or manifest.get("schema_version") != 1
        or not isinstance(manifest.get("jobs"), list)
    ):
        raise CondorError(f"unsupported KIM Condor manifest: {path}")
    entries = manifest["jobs"]
    if len(entries) != len(jobs):
        raise CondorError("KIM Condor manifest does not match the number of staged jobs")
    for entry, staged in zip(entries, jobs):
        if not isinstance(entry, dict):
            raise CondorError("KIM Condor manifest contains an invalid job entry")
        for field_name in ("job", "job_identity", "scan_order_index", "profile_scale_factor"):
            if entry.get(field_name) != staged.get(field_name):
                raise CondorError(
                    f"KIM Condor manifest {field_name} does not match staged job "
                    f"{staged.get('job')}"
                )
        if entry.get("state") not in {
            "prepared",
            "submitting",
            "submitted",
            "adopted",
            "ambiguous",
            "missing",
            "terminal",
            "removed_by_driver",
        }:
            raise CondorError(f"KIM Condor manifest has an invalid state for {entry['job']}")
    return manifest


def _save_manifest(root: Path, manifest: dict[str, object]) -> None:
    manifest["updated_at_utc"] = _utc_now()
    _write_json_atomic(root / MANIFEST_NAME, manifest)


def _submit_description(job_directory: Path, payload: dict[str, object], plan: KimCondorPlan):
    return CondorSubmitSpec(
        executable=plan.python_executable,
        initialdir=job_directory,
        arguments=(WORKER_NAME,),
        request_cpus=plan.request_cpus,
        request_memory_mb=plan.request_memory_mb,
        getenv=False,
        should_transfer_files=plan.should_transfer_files,
        exclude_machines=plan.exclude_machines,
        accounting={
            "KAMELJobIdentity": str(payload["job_identity"]),
            "KIMBackend": str(payload["backend"]),
            "KIMScanIndex": int(payload["scan_order_index"]),
            "KIMProfileScale": float(payload["profile_scale_factor"]),
        },
    )


def _ensure_file(path: Path, content: bytes) -> None:
    if path.exists() or path.is_symlink():
        if path.is_symlink() or not path.is_file() or path.read_bytes() != content:
            raise CondorError(f"existing KIM Condor artifact does not match staged content: {path}")
        return
    try:
        with path.open("xb") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
    except FileExistsError as error:
        raise CondorError(f"refusing to overwrite KIM Condor artifact: {path}") from error


def _ensure_submission_artifacts(
    directory: Path, payload: dict[str, object], plan: KimCondorPlan
) -> str:
    source = Path(__file__).with_name(WORKER_NAME)
    if not source.is_file():
        raise CondorError(f"KIM Condor worker is missing: {source}")
    worker_content = source.read_bytes()
    _ensure_file(directory / WORKER_NAME, worker_content)
    spec = _submit_description(directory, payload, plan)
    rendered = spec.render()
    submit_path = directory / "condor.submit"
    if submit_path.exists() or submit_path.is_symlink():
        if submit_path.is_symlink() or not submit_path.is_file():
            raise CondorError(f"unsafe existing KIM submit description: {submit_path}")
        if submit_path.read_text(encoding="utf-8") != rendered:
            raise CondorError(
                f"existing KIM submit description differs from the plan: {submit_path}"
            )
    else:
        spec.write(directory)
    return rendered


def _cluster_records(ads: tuple[dict[str, object], ...] | list[dict[str, object]]):
    records = {}
    for ad in ads:
        cluster = ad.get("cluster")
        if isinstance(cluster, bool) or not isinstance(cluster, int):
            continue
        previous = records.get(cluster)
        if previous is None or (
            previous.get("ad_source") != "condor_q" and ad.get("ad_source") == "condor_q"
        ):
            records[cluster] = ad
    return records


def _reconcile_manifest(
    manifest: dict[str, object], *, plan: KimCondorPlan, adopt_existing: bool
) -> list[int]:
    entries = manifest["jobs"]
    ambiguous = []
    for entry in entries:
        if entry["state"] == "submitting":
            entry["state"] = "ambiguous"
            entry["error"] = "previous process ended before submission acknowledgement"
        if entry["state"] == "ambiguous":
            ambiguous.append(entry["job"])
    if ambiguous:
        raise CondorError(f"ambiguous KIM Condor submission for {ambiguous}; refusing resubmission")
    clusters = [entry["cluster"] for entry in entries if entry.get("cluster") is not None]
    adopted = []
    if not clusters or not adopt_existing:
        return adopted
    try:
        ads = query_job_ads(clusters, config=plan.condor, include_history=True)
    except CondorError as error:
        raise CondorError(f"could not reconcile recorded KIM clusters: {error}") from error
    by_cluster = _cluster_records(ads)
    for entry in entries:
        cluster = entry.get("cluster")
        if cluster is None:
            continue
        ad = by_cluster.get(cluster)
        if ad is None:
            entry["state"] = "missing"
            entry["error"] = "recorded cluster is absent from queue and history"
            entry["queue_status"] = None
            continue
        if ad.get("proc") != 0 or ad.get("job_identity") != entry["job_identity"]:
            raise CondorError(
                f"recorded cluster {cluster} does not match staged job identity "
                f"for {entry['job']}"
            )
        entry["queue_status"] = _queue_summary(ad)
        if ad.get("ad_source") == "condor_q":
            entry["state"] = "adopted"
            entry["error"] = None
            adopted.append(cluster)
        else:
            entry["state"] = (
                "terminal" if ad.get("job_status_code") in _TERMINAL_CODES else "submitted"
            )
    return adopted


def submit_condor_sweep(
    root: str | Path,
    *,
    plan: KimCondorPlan,
    dry_run: bool = False,
    adopt_existing: bool = True,
) -> dict[str, object]:
    """Submit each staged scan point once, adopting only identity-matched clusters."""

    if not isinstance(dry_run, bool) or not isinstance(adopt_existing, bool):
        raise CondorError("dry_run and adopt_existing must be booleans")
    run_root = Path(root).expanduser().absolute()
    staged = _read_staged_jobs(run_root)
    manifest = _load_manifest(run_root, staged)
    if dry_run:
        jobs = []
        for entry, payload in zip(manifest["jobs"], staged):
            submit_text = _submit_description(
                run_root / str(entry["directory"]), payload, plan
            ).render()
            jobs.append({**entry, "submit_description": submit_text})
        return {"dry_run": True, "jobs": jobs, "adopted_clusters": []}

    adopted: list[int] = []
    with _queue_lock(run_root):
        # Re-read under the lock so concurrent callers observe acknowledged clusters.
        manifest = _load_manifest(run_root, staged)
        try:
            adopted = _reconcile_manifest(manifest, plan=plan, adopt_existing=adopt_existing)
        except CondorError:
            _save_manifest(run_root, manifest)
            raise
        _save_manifest(run_root, manifest)
        for entry, payload in zip(manifest["jobs"], staged):
            if entry["state"] != "prepared":
                continue
            directory = run_root / str(entry["directory"])
            _ensure_submission_artifacts(directory, payload, plan)
            entry["state"] = "submitting"
            entry["error"] = None
            _save_manifest(run_root, manifest)
            try:
                result = run_condor(
                    [str(directory / "condor.submit")],
                    config=plan.condor,
                    tool="condor_submit",
                )
                if result.returncode != 0:
                    raise CondorError(
                        f"condor_submit failed with exit code {result.returncode}: "
                        f"{result.stderr.strip()}"
                    )
                cluster = parse_condor_submit_output(result.stdout)
            except CondorError as error:
                entry["state"] = "ambiguous"
                entry["error"] = str(error)
                manifest["state"] = "ambiguous"
                _save_manifest(run_root, manifest)
                raise CondorError(
                    f"ambiguous KIM Condor submission for {entry['job']}: {error}"
                ) from error
            entry["cluster"] = cluster
            entry["state"] = "submitted"
            entry["error"] = None
            manifest["state"] = "submitting"
            _save_manifest(run_root, manifest)
        manifest["state"] = "submitted"
        _save_manifest(run_root, manifest)
    return {
        "dry_run": False,
        "jobs": [dict(entry) for entry in manifest["jobs"]],
        "adopted_clusters": sorted(adopted),
    }


def _queue_summary(ad: dict[str, object] | None) -> dict[str, object] | None:
    if ad is None:
        return None
    return {
        key: ad.get(key)
        for key in (
            "job_status",
            "job_status_code",
            "exit_code",
            "remote_host",
            "remote_wall_clock_s",
            "num_holds",
            "hold_reason",
            "ad_source",
        )
    }


def _status_record(entry: dict[str, object], ad: dict[str, object] | None) -> dict[str, object]:
    base: dict[str, object] = {
        "job": entry["job"],
        "job_identity": entry["job_identity"],
        "scan_order_index": entry["scan_order_index"],
        "profile_scale_factor": entry["profile_scale_factor"],
        "cluster": entry["cluster"],
        "job_status": "missing" if ad is None else ad.get("job_status"),
        "job_status_code": None if ad is None else ad.get("job_status_code"),
        "exit_code": None if ad is None else ad.get("exit_code"),
        "remote_host": None if ad is None else ad.get("remote_host"),
        "remote_wall_clock_s": None if ad is None else ad.get("remote_wall_clock_s"),
        "num_holds": None if ad is None else ad.get("num_holds"),
        "hold_reason": None if ad is None else ad.get("hold_reason"),
        "failure_kind": None,
        "failure_reason": None,
        "terminal": False,
        "driver_removed": False,
        "ad_source": None if ad is None else ad.get("ad_source"),
    }
    code = base["job_status_code"]
    if ad is None:
        base["failure_kind"] = "condor_missing"
        base["failure_reason"] = "cluster was absent from queue and history"
        base["terminal"] = True
    elif code == 3:
        base["failure_kind"] = "condor_removed"
        base["failure_reason"] = ad.get("hold_reason") or "job was removed from the pool"
        base["terminal"] = True
    elif code == 5:
        base["failure_kind"] = "condor_held"
        base["failure_reason"] = ad.get("hold_reason") or "job is held"
        base["terminal"] = True
    elif code == 4:
        base["terminal"] = True
        if base["exit_code"] is None:
            base["failure_kind"] = "condor_exit_unknown"
            base["failure_reason"] = "completed job has no scheduler exit code"
        elif base["exit_code"] != 0:
            base["failure_kind"] = "condor_nonzero_exit"
            base["failure_reason"] = f"job completed with exit code {base['exit_code']}"
    return base


def _save_status(root: Path, document: dict[str, object]) -> None:
    document["updated_at_utc"] = _utc_now()
    _write_json_atomic(root / STATUS_NAME, document)


def wait_condor_sweep(
    root: str | Path,
    *,
    plan: KimCondorPlan,
) -> tuple[dict[str, object], ...]:
    """Poll all acknowledged jobs to terminal status, removing them at the deadline."""

    run_root = Path(root).expanduser().absolute()
    staged = _read_staged_jobs(run_root)
    manifest = _load_manifest(run_root, staged)
    if any(entry["state"] in {"ambiguous", "submitting"} for entry in manifest["jobs"]):
        raise CondorError("cannot wait for KIM jobs with ambiguous submission state")
    entries = manifest["jobs"]
    if any(entry.get("cluster") is None for entry in entries):
        raise CondorError("all KIM scan points must have acknowledged clusters before waiting")
    clusters = [int(entry["cluster"]) for entry in entries]
    started = time.monotonic()
    deadline = started + plan.max_wall_clock_s
    status_document: dict[str, object] = {
        "schema_version": 1,
        "state": "polling",
        "started_at_utc": _utc_now(),
        "driver_removed_clusters": [],
        "statuses": [],
    }
    latest: list[dict[str, object]] = []
    while True:
        query_error = None
        try:
            ads = query_job_ads(clusters, config=plan.condor, include_history=True)
        except CondorError as error:
            query_error = str(error)
            status_document["state"] = "query_error"
            status_document["query_error"] = query_error
            if latest:
                for status in latest:
                    status["query_error"] = query_error
            else:
                latest = [
                    {
                        "job": entry["job"],
                        "job_identity": entry["job_identity"],
                        "scan_order_index": entry["scan_order_index"],
                        "profile_scale_factor": entry["profile_scale_factor"],
                        "cluster": entry["cluster"],
                        "job_status": "query_error",
                        "job_status_code": None,
                        "exit_code": None,
                        "failure_kind": "condor_query_error",
                        "failure_reason": query_error,
                        "query_error": query_error,
                        "terminal": False,
                        "driver_removed": False,
                    }
                    for entry in entries
                ]
            status_document["statuses"] = latest
            _save_status(run_root, status_document)
        if query_error is None:
            by_cluster = _cluster_records(ads)
            status_document.pop("query_error", None)
            latest = []
            for entry in entries:
                ad = by_cluster.get(int(entry["cluster"]))
                if ad is not None and (
                    ad.get("proc") != 0 or ad.get("job_identity") != entry["job_identity"]
                ):
                    raise CondorError(
                        f"queue ad for cluster {entry['cluster']} does not match {entry['job']}"
                    )
                status = _status_record(entry, ad)
                latest.append(status)
                entry["queue_status"] = _queue_summary(ad)
                if status["terminal"]:
                    entry["state"] = "terminal"
            _save_manifest(run_root, manifest)
            status_document["statuses"] = latest
            status_document["state"] = (
                "completed"
                if all(status["failure_kind"] is None for status in latest)
                else ("failed" if all(status["terminal"] for status in latest) else "polling")
            )
            _save_status(run_root, status_document)
        pending = [status for status in latest if not status["terminal"]]
        if not pending:
            status_document["finished_at_utc"] = _utc_now()
            _save_status(run_root, status_document)
            return tuple(dict(status) for status in latest)
        now = time.monotonic()
        if now >= deadline:
            break
        time.sleep(min(plan.poll_interval_s, max(0.0, deadline - now)))

    pending_clusters = [int(status["cluster"]) for status in latest if not status["terminal"]]
    removal_error = None
    try:
        remove_clusters(
            pending_clusters,
            config=plan.condor,
            reason=f"KIM Condor driver wall-clock deadline exceeded ({plan.max_wall_clock_s:g} s)",
        )
    except CondorError as error:
        removal_error = str(error)
    removed = [] if removal_error else pending_clusters
    status_document["driver_removed_clusters"] = removed
    if removal_error is not None:
        status_document["removal_error"] = removal_error
    for status in latest:
        if not status["terminal"]:
            status["failure_kind"] = "condor_driver_wall_clock"
            status["failure_reason"] = f"driver deadline exceeded after {plan.max_wall_clock_s:g} s"
            status["driver_removed"] = int(status["cluster"]) in removed
            if removal_error is not None:
                status["removal_error"] = removal_error
            status["terminal"] = True
    for entry, status in zip(entries, latest):
        if status["driver_removed"]:
            entry["state"] = "removed_by_driver"
            entry["error"] = status["failure_reason"]
    manifest["state"] = "failed"
    _save_manifest(run_root, manifest)
    status_document["statuses"] = latest
    status_document["state"] = "failed"
    status_document["finished_at_utc"] = _utc_now()
    _save_status(run_root, status_document)
    return tuple(dict(status) for status in latest)


def _read_json_object(path: Path, *, description: str) -> dict[str, object]:
    if path.is_symlink() or not path.is_file():
        raise CondorError(f"{description} is missing or unsafe: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CondorError(f"could not read {description} {path}: {error}") from error
    if not isinstance(payload, dict):
        raise CondorError(f"{description} must contain a JSON object: {path}")
    return payload


def _finite_record_number(record: dict[str, object], name: str, *, context: str) -> float:
    value = record.get(name)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise CondorError(f"{context} has no valid {name}")
    numeric = float(value)
    if not math.isfinite(numeric):
        raise CondorError(f"{context} {name} is non-finite")
    return numeric


def _check_worker_record(
    job_directory: Path,
    job: dict[str, object],
    manifest_entry: dict[str, object],
) -> dict[str, object]:
    record = _read_json_object(
        job_directory / "condor_run_record.json", description="KIM worker record"
    )
    for field_name in (
        "job_identity",
        "backend",
        "scan_order_index",
        "profile_scale_factor",
        "executable_sha256",
    ):
        expected = job.get(field_name)
        if record.get(field_name) != expected:
            raise CondorError(
                f"KIM worker record {field_name} does not match staged job {job.get('job')}"
            )
    if manifest_entry.get("job_identity") != job.get("job_identity"):
        raise CondorError(f"KIM manifest identity does not match staged job {job.get('job')}")
    if record.get("schema_version") != 1:
        raise CondorError(f"unsupported KIM worker record schema for {job.get('job')}")
    exit_code = record.get("exit_code")
    if isinstance(exit_code, bool) or not isinstance(exit_code, int) or exit_code != 0:
        raise CondorError(f"KIM worker for {job.get('job')} did not exit successfully")
    if record.get("status") != "succeeded":
        raise CondorError(f"KIM worker for {job.get('job')} has status {record.get('status')!r}")
    host = record.get("host")
    if not isinstance(host, str) or not host.strip():
        raise CondorError(f"KIM worker record for {job.get('job')} has no execute host")
    return record


def _validate_collection_status(
    job_name: str,
    status: dict[str, object] | None,
    staged_job: dict[str, object],
    manifest_entry: dict[str, object],
) -> dict[str, object]:
    if status is None:
        raise CondorError(f"KIM Condor status for job {job_name} is missing")
    expected = {
        "job_identity": staged_job["job_identity"],
        "scan_order_index": staged_job["scan_order_index"],
        "profile_scale_factor": staged_job["profile_scale_factor"],
        "cluster": manifest_entry["cluster"],
    }
    for field_name, value in expected.items():
        if status.get(field_name) != value:
            raise CondorError(
                f"KIM Condor status {field_name} does not match staged job {job_name}"
            )
    if status.get("terminal") is not True:
        raise CondorError(f"KIM Condor job {job_name} is still pending or lacks terminal status")
    failure = status.get("failure_kind")
    if failure is not None:
        raise CondorError(
            f"KIM Condor job {job_name} has scheduler failure {failure}: "
            f"{status.get('failure_reason')}"
        )
    if status.get("job_status") != "completed" or status.get("job_status_code") != 4:
        raise CondorError(f"KIM Condor job {job_name} did not complete successfully")
    if status.get("exit_code") != 0:
        raise CondorError(
            f"KIM Condor job {job_name} has scheduler exit code {status.get('exit_code')}"
        )
    return status


def _find_jpar_profile(directory: Path, collision_model: str) -> Path:
    pattern = f"jpar_{collision_model}.dat"
    matches = sorted(path for path in directory.rglob(pattern) if path.is_file())
    if not matches:
        raise CondorError(f"no {pattern} output found under {directory}")
    if len(matches) != 1:
        raise CondorError(f"ambiguous {pattern} outputs under {directory}")
    return matches[0]


def _integrate_jpar_profile(path: Path, metric: JparCurrentMetric) -> float:
    try:
        rows = np.loadtxt(path, ndmin=2)
    except (OSError, ValueError) as error:
        raise CondorError(f"could not read declared jpar profile {path}: {error}") from error
    if rows.ndim != 2 or rows.shape[0] < 2 or rows.shape[1] <= metric.current_column:
        raise CondorError(
            f"jpar profile {path} must have at least two rows and column "
            f"{metric.current_column}"
        )
    if not np.all(np.isfinite(rows)):
        raise CondorError(f"jpar profile contains non-finite values: {path}")
    radius = rows[:, 0]
    if np.any(np.diff(radius) <= 0.0):
        raise CondorError(f"jpar radius grid must be strictly increasing: {path}")
    column = rows[:, metric.current_column]
    integral = float(np.sum(0.5 * (column[:-1] + column[1:]) * np.diff(radius)))
    if not math.isfinite(integral):
        raise CondorError(f"integrated jpar profile is non-finite: {path}")
    return integral


def _resonance_from_curve(curve: list[dict[str, object]]) -> dict[str, object]:
    factors = [float(point["profile_scale_factor"]) for point in curve]
    values = [float(point["integrated_parallel_current"]) for point in curve]
    if len(curve) < 3:
        return {
            "status": "insufficient_points",
            "profile_scale_factor": None,
            "current_value": None,
            "current_unit": curve[0].get("current_unit") if curve else None,
            "candidate_profile_scale_factors": [],
        }
    if len(set(factors)) != len(factors):
        return {
            "status": "invalid_scan",
            "profile_scale_factor": None,
            "current_value": None,
            "current_unit": curve[0].get("current_unit"),
            "candidate_profile_scale_factors": [],
        }
    ordered = sorted(zip(factors, values), key=lambda item: item[0])
    ordered_factors = [item[0] for item in ordered]
    ordered_values = [item[1] for item in ordered]
    minima = [
        index
        for index in range(1, len(ordered_values) - 1)
        if ordered_values[index] < ordered_values[index - 1]
        and ordered_values[index] < ordered_values[index + 1]
    ]
    unit = curve[0].get("current_unit")
    if not minima:
        return {
            "status": "no_interior_minimum",
            "profile_scale_factor": None,
            "current_value": None,
            "current_unit": unit,
            "candidate_profile_scale_factors": [],
        }
    if len(minima) > 1:
        return {
            "status": "multiple_minima",
            "profile_scale_factor": None,
            "current_value": None,
            "current_unit": unit,
            "candidate_profile_scale_factors": [ordered_factors[index] for index in minima],
        }
    index = minima[0]
    return {
        "status": "interior_minimum",
        "profile_scale_factor": ordered_factors[index],
        "current_value": ordered_values[index],
        "current_unit": unit,
        "candidate_profile_scale_factors": [ordered_factors[index]],
    }


def collect_condor_sweep(
    root: str | Path,
    *,
    plan: KimCondorPlan,
) -> KimCondorScanResults:
    """Validate completed job records and produce only explicitly defined scalars."""

    run_root = Path(root).expanduser().absolute()
    staged = _read_staged_jobs(run_root)
    manifest = _load_manifest(run_root, staged)
    status_path = run_root / STATUS_NAME
    status_by_job: dict[str, dict[str, object]] = {}
    if status_path.exists() or status_path.is_symlink():
        status_document = _read_json_object(status_path, description="KIM Condor status")
        if status_document.get("schema_version") != 1:
            raise CondorError(f"unsupported KIM Condor status schema: {status_path}")
        raw_statuses = status_document.get("statuses")
        if not isinstance(raw_statuses, list):
            raise CondorError(f"KIM Condor status has no status list: {status_path}")
        for status in raw_statuses:
            if not isinstance(status, dict) or not isinstance(status.get("job"), str):
                raise CondorError(f"KIM Condor status contains an invalid entry: {status_path}")
            if status["job"] in status_by_job:
                raise CondorError(f"KIM Condor status has duplicate job {status['job']}")
            status_by_job[status["job"]] = status
        expected_jobs = {str(job["job"]) for job in staged}
        if status_by_job.keys() != expected_jobs:
            missing = sorted(expected_jobs - status_by_job.keys())
            unexpected = sorted(status_by_job.keys() - expected_jobs)
            raise CondorError(
                f"KIM Condor status job list does not match staged scan; "
                f"missing={missing}, unexpected={unexpected}"
            )
    else:
        raise CondorError("KIM Condor status is missing; wait for the scan before collecting")

    curve: list[dict[str, object]] = []
    notes: list[str] = []
    for staged_job, entry in zip(staged, manifest["jobs"]):
        job_name = str(staged_job["job"])
        if staged_job.get("backend") != plan.backend:
            raise CondorError(f"plan backend does not match staged KIM job {job_name}")
        directory = run_root / str(entry["directory"])
        if entry.get("cluster") is None:
            raise CondorError(f"KIM Condor job {job_name} has no acknowledged cluster")
        status = _validate_collection_status(
            job_name, status_by_job.get(job_name), staged_job, entry
        )
        worker = _check_worker_record(directory, staged_job, entry)
        point: dict[str, object] = {
            "job": job_name,
            "job_identity": staged_job["job_identity"],
            "backend": plan.backend,
            "scan_order_index": staged_job["scan_order_index"],
            "profile_scale_factor": float(staged_job["profile_scale_factor"]),
            "cluster": entry["cluster"],
            "worker_status": worker["status"],
            "exit_code": worker["exit_code"],
            "worker_host": worker["host"],
            "wall_clock_s": worker.get("wall_clock_s"),
            "scheduler_status": status.get("job_status"),
            "remote_host": status.get("remote_host"),
            "integrated_parallel_current_real": None,
            "integrated_parallel_current_imag": None,
            "integrated_parallel_current": None,
            "current_unit": None,
            "jpar_profile": None,
        }
        if plan.backend == "kamel_kim_python":
            run_type = staged_job["base_config"]["run"]["run_type"]
            if run_type == "electrostatic_periodic":
                real = _finite_record_number(
                    worker, "integrated_parallel_current_real", context=job_name
                )
                imag = _finite_record_number(
                    worker, "integrated_parallel_current_imag", context=job_name
                )
                point["integrated_parallel_current_real"] = real
                point["integrated_parallel_current_imag"] = imag
                component = plan.api_current_component
                if component is not None:
                    if (
                        worker.get("api_current_component") != component
                        or worker.get("api_current_unit") != plan.api_current_unit
                    ):
                        raise CondorError(f"API current metric provenance mismatch for {job_name}")
                    if component == "real":
                        scalar = real
                    elif component == "imag":
                        scalar = imag
                    else:
                        scalar = math.hypot(real, imag)
                    recorded_scalar = _finite_record_number(
                        worker, "api_current_value", context=job_name
                    )
                    if recorded_scalar != scalar:
                        raise CondorError(
                            f"API scalar current does not match raw complex value for {job_name}"
                        )
                    point["integrated_parallel_current"] = scalar
                    point["current_unit"] = plan.api_current_unit
                elif worker.get("api_current_value") is not None:
                    raise CondorError(
                        f"API worker supplied an undeclared scalar current for {job_name}"
                    )
            elif plan.api_current_component is not None:
                raise CondorError(
                    "API integrated-current metrics require an electrostatic-periodic scan"
                )
        elif plan.jpar_current_metric is not None:
            expected_metric = plan.jpar_current_metric.model_dump(mode="json")
            if staged_job.get("jpar_current_metric") != expected_metric:
                raise CondorError(f"declared jpar metric does not match staged job {job_name}")
            profile = _find_jpar_profile(directory, plan.jpar_current_metric.collision_model)
            scalar = _integrate_jpar_profile(profile, plan.jpar_current_metric)
            point["integrated_parallel_current"] = scalar
            point["current_unit"] = plan.jpar_current_metric.current_unit
            point["jpar_profile"] = str(profile)
        curve.append(point)

    if plan.backend == "kamel_kim_python" and plan.api_current_component is None:
        notes.append(
            "raw API integrated-current real and imaginary values are preserved; "
            "no scalar shielding-current metric was selected"
        )
    elif plan.backend == "kim_x_namelist" and plan.jpar_current_metric is None:
        notes.append(
            "no jpar current-column metric was declared; no scalar shielding-current "
            "curve or resonance minimum is selected"
        )

    if any(point["integrated_parallel_current"] is None for point in curve):
        resonance: dict[str, object] = {
            "status": "insufficient_scalars",
            "profile_scale_factor": None,
            "current_value": None,
            "current_unit": None,
            "candidate_profile_scale_factors": [],
        }
    else:
        scalar_values = [float(point["integrated_parallel_current"]) for point in curve]
        if not all(math.isfinite(value) for value in scalar_values):
            raise CondorError("KIM shielding-current curve contains non-finite values")
        resonance = _resonance_from_curve(curve)
    notes.append("curve points retain requested scan order; minimum analysis sorts by scale factor")
    if plan.jpar_current_metric is not None:
        notes.append(
            "namelist current values are trapezoid-integrated over the declared profile's "
            "radius grid; units and column meaning are caller-declared"
        )
    return KimCondorScanResults(curve=tuple(curve), resonance=resonance, notes=tuple(notes))
