"""Validated KIM Condor scan plans and immutable per-point staging."""

from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import f90nml
import numpy as np
from kim.config import ElectrostaticPeriodicRun, ProfileConfig, SimulationConfig
from kim.condor_client import CondorError, CondorToolConfig, verify_shared_filesystem
from kim.errors import ProfileError
from kim.profiles import ProfileSet, _read_profile
from kim.sweep import SweepSpec
from pydantic import BaseModel, ConfigDict, Field, model_validator


class _CondorModel(BaseModel):
    model_config = ConfigDict(
        allow_inf_nan=False,
        arbitrary_types_allowed=True,
        extra="forbid",
        frozen=True,
    )


class JparCurrentMetric(_CondorModel):
    """Explicit interpretation of one current column in a direct solver output."""

    current_column: int = Field(ge=1)
    current_unit: str = Field(min_length=1)
    collision_model: str = Field(min_length=1)


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
