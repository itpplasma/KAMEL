"""Validated sequential scalar and radial-electric-field sweeps."""

from __future__ import annotations

import os
import re
import uuid
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Literal

import numpy as np
from kim.config import KimModel, ProfileConfig, SimulationConfig
from kim.errors import SweepError
from kim.profiles import ProfileSet, _read_profile
from kim.runs import RunStatus
from kim.simulation import RunResult, Simulation
from pydantic import BaseModel, ConfigDict, Field, JsonValue, ValidationError, model_validator

SweepValue = int | float


class SweepModel(BaseModel):
    """Strict immutable base for sweep requests and persisted metadata."""

    model_config = ConfigDict(
        allow_inf_nan=False,
        extra="forbid",
        frozen=True,
        str_strip_whitespace=True,
    )


class LinearRange(SweepModel):
    """An inclusive, evenly spaced scalar value source."""

    start: float
    stop: float
    count: int = Field(ge=2)

    def generate(self) -> tuple[float, ...]:
        return tuple(float(value) for value in np.linspace(self.start, self.stop, self.count))


class ParameterSweep(SweepModel):
    """One validated scalar model field varied over an ordered value sequence."""

    base: SimulationConfig
    parameter: str = Field(min_length=1)
    values: tuple[SweepValue, ...] | LinearRange
    continue_on_failure: bool = True

    @model_validator(mode="after")
    def validate_parameter_and_values(self) -> ParameterSweep:
        _resolve_parameter(self.base, self.parameter)
        if not self.expanded_values:
            raise ValueError("parameter sweep requires at least one value")
        return self

    @property
    def expanded_values(self) -> tuple[SweepValue, ...]:
        if isinstance(self.values, LinearRange):
            return self.values.generate()
        return self.values

    def configurations(self) -> tuple[SimulationConfig, ...]:
        """Return independently validated child configurations in value order."""

        return tuple(
            _replace_parameter(self.base, self.parameter, value) for value in self.expanded_values
        )


class ProfileScale(SweepModel):
    """Scale one supported profile's values without modifying its radial grid."""

    kind: Literal["profile_scale"] = "profile_scale"
    profile: Literal["Er"] = "Er"
    values: tuple[float, ...] = Field(min_length=1)


class SweepSpec(SweepModel):
    """An immutable base request with one typed profile variation."""

    base: SimulationConfig
    variation: ProfileScale
    continue_on_failure: bool = True


class SweepStatus(str, Enum):
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    PARTIAL_FAILURE = "partial_failure"
    FAILED = "failed"
    INTERRUPTED = "interrupted"


class SweepManifest(SweepModel):
    """Persistent ordered identity and progress for a sequential sweep."""

    schema_version: Literal[1] = 1
    sweep_id: str
    label: str | None = None
    status: SweepStatus = SweepStatus.RUNNING
    created_at: datetime
    finished_at: datetime | None = None
    base_parameters: dict[str, JsonValue]
    variation: dict[str, JsonValue]
    values: tuple[JsonValue, ...]
    continue_on_failure: bool
    child_run_ids: tuple[str, ...] = ()
    child_statuses: tuple[str, ...] = ()


class SweepResult(SweepModel):
    """Completed sweep identity and ordered child run results."""

    sweep_id: str
    status: SweepStatus
    manifest: Path
    children: tuple[RunResult, ...]

    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)


def run_sweep(
    spec: ParameterSweep | SweepSpec,
    *,
    executable: Path | str | None = None,
    runs_directory: Path | str = "runs",
    label: str | None = None,
    timeout: float | None = None,
    environment: Mapping[str, str] | None = None,
) -> SweepResult:
    """Prepare and execute sweep children sequentially in deterministic order."""

    runs_root = Path(runs_directory).expanduser().absolute()
    sweep_id = _sweep_id(label)
    sweep_root = runs_root / "sweeps" / sweep_id
    sweep_root.mkdir(parents=True, exist_ok=False)

    if isinstance(spec, ParameterSweep):
        values = spec.expanded_values
        configurations = spec.configurations()
        descriptors = tuple(
            {"kind": "parameter", "parameter": spec.parameter, "value": value} for value in values
        )
        variation: dict[str, JsonValue] = {
            "kind": "parameter",
            "parameter": spec.parameter,
        }
        continue_on_failure = spec.continue_on_failure
    else:
        values = spec.variation.values
        configurations = _profile_scale_configurations(spec, sweep_root)
        descriptors = tuple(
            {
                "kind": "profile_scale",
                "profile": spec.variation.profile,
                "factor": factor,
            }
            for factor in values
        )
        variation = {
            "kind": "profile_scale",
            "profile": spec.variation.profile,
        }
        continue_on_failure = spec.continue_on_failure

    manifest_path = sweep_root / "manifest.json"
    manifest = SweepManifest(
        sweep_id=sweep_id,
        label=label,
        created_at=datetime.now(timezone.utc),
        base_parameters=spec.base.model_dump(mode="json"),
        variation=variation,
        values=values,
        continue_on_failure=continue_on_failure,
    )
    _write_manifest(manifest_path, manifest)

    children: list[RunResult] = []
    for index, (config, descriptor) in enumerate(zip(configurations, descriptors, strict=True)):
        requested: dict[str, Any] = {
            "configuration": config.model_dump(mode="json"),
            "variation": descriptor,
        }
        simulation = Simulation(
            config,
            executable=executable,
            runs_directory=runs_root,
            label=f"{label or 'sweep'}-{index:04d}",
            timeout=timeout,
            parent_sweep_id=sweep_id,
            requested_parameters=requested,
        )
        try:
            child = simulation.run(environment=environment)
        except KeyboardInterrupt:
            if simulation._last_result is not None:
                children.append(simulation._last_result)
            manifest = manifest.model_copy(
                update={
                    "status": SweepStatus.INTERRUPTED,
                    "finished_at": datetime.now(timezone.utc),
                    "child_run_ids": tuple(item.run_id for item in children),
                    "child_statuses": tuple(item.status.value for item in children),
                }
            )
            _write_manifest(manifest_path, manifest)
            raise
        children.append(child)
        manifest = manifest.model_copy(
            update={
                "child_run_ids": tuple(item.run_id for item in children),
                "child_statuses": tuple(item.status.value for item in children),
            }
        )
        _write_manifest(manifest_path, manifest)
        if child.status is not RunStatus.SUCCEEDED and not continue_on_failure:
            break

    status = _final_status(children, requested_count=len(values))
    manifest = manifest.model_copy(
        update={"status": status, "finished_at": datetime.now(timezone.utc)}
    )
    _write_manifest(manifest_path, manifest)
    return SweepResult(
        sweep_id=sweep_id,
        status=status,
        manifest=manifest_path,
        children=tuple(children),
    )


def _resolve_parameter(config: SimulationConfig, public_path: str) -> tuple[str, ...]:
    parts = tuple(public_path.split("."))
    if not parts or any(not part for part in parts):
        raise SweepError(f"invalid sweep parameter path: {public_path!r}")
    internal = ("run", *parts) if parts[0] in {"periodic", "terms"} else parts
    current: Any = config
    for index, part in enumerate(internal):
        if not isinstance(current, KimModel) or part not in type(current).model_fields:
            raise SweepError(f"unknown sweep parameter: {public_path}")
        field = type(current).model_fields[part]
        if index == len(internal) - 1:
            metadata = field.json_schema_extra or {}
            if not isinstance(metadata, Mapping) or metadata.get("sweepable") is not True:
                raise SweepError(f"parameter is not sweepable: {public_path}")
        else:
            current = getattr(current, part)
    return internal


def _replace_parameter(
    base: SimulationConfig, public_path: str, value: SweepValue
) -> SimulationConfig:
    internal = _resolve_parameter(base, public_path)
    data = base.model_dump(mode="python")
    target = data
    for part in internal[:-1]:
        target = target[part]
    target[internal[-1]] = value
    try:
        return SimulationConfig.model_validate(data)
    except ValidationError as error:
        raise SweepError(
            f"sweep parameter {public_path} has invalid value {value!r}: {error}"
        ) from error


def _profile_scale_configurations(
    spec: SweepSpec, sweep_root: Path
) -> tuple[SimulationConfig, ...]:
    source_profiles = ProfileSet.from_simulation(spec.base)
    source_profiles.validate_for(spec.base)
    configurations = []
    for index, factor in enumerate(spec.variation.values):
        destination = sweep_root / "inputs" / f"{index:04d}-profiles"
        source_profiles.copy_to(destination)
        electric_field = destination / spec.base.profiles.radial_electric_field_file
        if not electric_field.is_file():
            raise SweepError(f"Er profile is required for scaling: {electric_field}")
        profile = _read_profile("Er", "statV/cm", electric_field)
        values = np.column_stack((profile.radius, profile.values))
        values[:, 1] *= factor
        np.savetxt(electric_field, values, fmt="%.16e")
        profiles = ProfileConfig.model_validate(
            {**spec.base.profiles.model_dump(mode="python"), "directory": destination}
        )
        configurations.append(
            SimulationConfig.model_validate(
                {**spec.base.model_dump(mode="python"), "profiles": profiles}
            )
        )
    return tuple(configurations)


def _final_status(children: Sequence[RunResult], *, requested_count: int) -> SweepStatus:
    successful = sum(child.status is RunStatus.SUCCEEDED for child in children)
    if successful == requested_count:
        return SweepStatus.SUCCEEDED
    if successful == 0:
        return SweepStatus.FAILED
    return SweepStatus.PARTIAL_FAILURE


def _sweep_id(label: str | None) -> str:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    raw_label = label or "sweep"
    slug = re.sub(r"[^a-z0-9]+", "-", raw_label.casefold()).strip("-")
    return f"{timestamp}-{slug[:48] or 'sweep'}-{uuid.uuid4().hex[:8]}"


def _write_manifest(path: Path, manifest: SweepManifest) -> None:
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("w", encoding="utf-8") as stream:
            stream.write(manifest.model_dump_json(indent=2))
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except OSError as error:
        temporary.unlink(missing_ok=True)
        raise SweepError(f"could not write sweep manifest {path}: {error}") from error
