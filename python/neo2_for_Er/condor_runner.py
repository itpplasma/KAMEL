"""Stage NEO-2 surface jobs for bounded HTCondor execution."""

from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .condor import CondorError, CondorSubmitSpec, CondorToolConfig, verify_shared_filesystem

MANIFEST_NAME = "condor_jobs.json"
STATUS_NAME = "condor_status.json"
JOB_INPUT_NAME = "condor_job.json"
WORKER_NAME = "condor_worker.py"


class Neo2CondorError(RuntimeError):
    """A NEO-2 Condor plan, staged job, or result is invalid."""


def detect_parallel_runtime(executable: Path) -> dict[str, bool]:
    """Report OpenMP/MPI linkage when `ldd` can inspect the solver binary."""

    try:
        result = subprocess.run(
            ["ldd", str(executable)], capture_output=True, text=True, timeout=30.0, check=False
        )
    except (OSError, subprocess.SubprocessError):
        return {"omp": False, "mpi": False}
    text = f"{result.stdout}\n{result.stderr}"
    return {
        "omp": any(name in text for name in ("libgomp", "libiomp5", "libomp")),
        "mpi": any(name in text for name in ("libmpi", "libmpich", "libopen-pal", "libpmi")),
    }


@dataclass(frozen=True, slots=True)
class Neo2CondorPlan:
    """Resources, bounds, and provenance for one staged surface set."""

    executable: Path
    python_executable: Path = field(default_factory=lambda: Path(sys.executable))
    condor: CondorToolConfig = field(default_factory=CondorToolConfig)
    shared_filesystem_prefixes: tuple[Path, ...] = ()
    should_transfer_files: str = "NEVER"
    shot: int | None = None
    requested_time_s: float | None = None
    omp_threads_per_process: int = 1
    request_cpus: int = 1
    request_memory_mb: int = 30_720
    surface_timeout_s: float = 300.0
    max_closure_periods: int = 10_000
    poll_interval_s: float = 30.0
    max_wall_clock_s: float = 6.0 * 3600.0
    exclude_machines: tuple[str, ...] = ("faepop43",)
    parallel_runtime: dict[str, bool] | None = None
    serial_solver_ok: bool = False

    def __post_init__(self) -> None:
        for name in ("executable", "python_executable"):
            path = Path(getattr(self, name))
            if not path.is_absolute() or not path.is_file():
                raise Neo2CondorError(f"{name} must be an existing absolute file path")
            object.__setattr__(self, name, path)
        if not isinstance(self.condor, CondorToolConfig):
            raise Neo2CondorError("condor must be a CondorToolConfig")
        if self.should_transfer_files not in {"NEVER", "IF_NEEDED", "ALWAYS"}:
            raise Neo2CondorError("should_transfer_files must be NEVER, IF_NEEDED, or ALWAYS")
        prefixes = tuple(Path(prefix) for prefix in self.shared_filesystem_prefixes)
        if self.should_transfer_files == "NEVER" and not prefixes:
            raise Neo2CondorError(
                "shared_filesystem_prefixes must be set when should_transfer_files is NEVER"
            )
        for prefix in prefixes:
            if not prefix.is_absolute():
                raise Neo2CondorError(f"shared filesystem prefix must be absolute: {prefix}")
        object.__setattr__(self, "shared_filesystem_prefixes", prefixes)
        object.__setattr__(self, "exclude_machines", tuple(self.exclude_machines))
        if any(not isinstance(name, str) or not name for name in self.exclude_machines):
            raise Neo2CondorError("exclude_machines must contain non-empty machine names")

        for name in ("omp_threads_per_process", "request_cpus", "request_memory_mb"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise Neo2CondorError(f"{name} must be a strictly positive integer")
        if self.request_cpus < self.omp_threads_per_process:
            raise Neo2CondorError("request_cpus cannot be below omp_threads_per_process")
        if not isinstance(self.serial_solver_ok, bool):
            raise Neo2CondorError("serial_solver_ok must be a boolean")

        runtime = self.parallel_runtime
        if runtime is None:
            runtime = detect_parallel_runtime(Path(self.executable))
        if not isinstance(runtime, dict):
            raise Neo2CondorError("parallel_runtime must be a mapping or None")
        unknown = set(runtime) - {"omp", "mpi"}
        if unknown or any(not isinstance(value, bool) for value in runtime.values()):
            raise Neo2CondorError("parallel_runtime must contain only boolean omp/mpi values")
        runtime = {"omp": bool(runtime.get("omp", False)), "mpi": bool(runtime.get("mpi", False))}
        object.__setattr__(self, "parallel_runtime", runtime)
        if runtime["omp"] and self.omp_threads_per_process == 1 and not self.serial_solver_ok:
            raise Neo2CondorError(
                "solver links OpenMP but omp_threads_per_process is 1; raise the thread count "
                "or declare serial_solver_ok=True"
            )

        if (
            isinstance(self.max_closure_periods, bool)
            or not isinstance(self.max_closure_periods, int)
            or self.max_closure_periods <= 0
        ):
            raise Neo2CondorError("max_closure_periods must be a strictly positive integer")
        for name in ("surface_timeout_s", "poll_interval_s", "max_wall_clock_s"):
            value = getattr(self, name)
            if isinstance(value, bool) or not math.isfinite(float(value)) or float(value) <= 0.0:
                raise Neo2CondorError(f"{name} must be strictly positive and finite")
            object.__setattr__(self, name, float(value))
        if self.shot is not None and (
            isinstance(self.shot, bool) or not isinstance(self.shot, int) or self.shot <= 0
        ):
            raise Neo2CondorError("shot must be a strictly positive integer or None")
        if self.requested_time_s is not None:
            value = self.requested_time_s
            if isinstance(value, bool) or not math.isfinite(float(value)) or float(value) < 0.0:
                raise Neo2CondorError("requested_time_s must be a non-negative finite number")
            object.__setattr__(self, "requested_time_s", float(value))


def _surface_directory_name(boozer_s: float) -> str:
    """Match the existing `local_runner.stage_surfaces` directory naming contract."""

    formatted = f"{boozer_s:.9e}".replace("e", "d").replace("d-", "m")
    return f"s{formatted}"


def _staged_surface_jobs(root: Path) -> tuple[tuple[Path, ...], np.ndarray]:
    rows_path = root / "surfaces.dat"
    jobs_list_path = root / "jobs_list.txt"
    try:
        rows = np.loadtxt(rows_path, ndmin=2)
        lines = jobs_list_path.read_text(encoding="utf-8").splitlines()
    except (OSError, ValueError, UnicodeError) as error:
        raise Neo2CondorError(
            f"could not read staged surface inputs under {root}: {error}"
        ) from error
    if rows.ndim != 2 or rows.shape[1] != 7 or not np.all(np.isfinite(rows)):
        raise Neo2CondorError("surfaces.dat must contain finite seven-column surface rows")
    names = [line.strip() for line in lines if line.strip()]
    if not names or any(
        Path(name).is_absolute()
        or Path(name).name != name
        or name in {".", ".."}
        or "/" in name
        or "\\" in name
        for name in names
    ):
        raise Neo2CondorError("jobs_list.txt contains an unsafe or empty surface job name")
    expected_names = [_surface_directory_name(float(row[0])) for row in rows]
    if names != expected_names or len(set(names)) != len(names):
        raise Neo2CondorError("jobs_list.txt does not match surfaces.dat exactly once and in order")
    resolved_root = root.resolve()
    jobs: list[Path] = []
    for name in names:
        job = root / name
        if not job.is_dir() or job.resolve().parent != resolved_root:
            raise Neo2CondorError(f"staged surface job is missing or escapes its work root: {job}")
        if not (job / "neo2.in").is_file():
            raise Neo2CondorError(f"staged surface input is missing: {job / 'neo2.in'}")
        jobs.append(job)
    return tuple(jobs), rows


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_name(f".{path.name}.{path.parent.name}.tmp")
    try:
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def stage_neo2_condor_jobs(
    work_directory: str | Path,
    *,
    plan: Neo2CondorPlan,
) -> tuple[Path, ...]:
    """Add a worker payload and submit description to each already staged surface."""

    root = Path(work_directory).expanduser().absolute()
    try:
        verify_shared_filesystem(
            root,
            allowed_prefixes=plan.shared_filesystem_prefixes,
            should_transfer_files=plan.should_transfer_files,
        )
    except CondorError as error:
        raise Neo2CondorError(str(error)) from error
    jobs, rows = _staged_surface_jobs(root)
    worker_source = Path(__file__).with_name(WORKER_NAME)
    if not worker_source.is_file():
        raise Neo2CondorError(f"Condor worker payload is missing: {worker_source}")
    artifacts = (WORKER_NAME, JOB_INPUT_NAME, "condor.submit")
    for job in jobs:
        for name in artifacts:
            if (job / name).exists():
                raise Neo2CondorError(f"refusing to overwrite existing payload: {job / name}")

    try:
        executable_hash = _sha256(plan.executable)
    except OSError as error:
        raise Neo2CondorError(f"could not fingerprint NEO-2 executable: {error}") from error
    staged: list[Path] = []
    for job, row in zip(jobs, rows, strict=True):
        surface = {
            "boozer_s": float(row[0]),
            "r_eff_cm": float(row[1]),
            "r_beg_cm": float(row[2]),
            "z_beg_cm": float(row[3]),
            "ti_eV": float(row[4]),
            "ne_cm3": float(row[5]),
            "kappa_cm_inv": float(row[6]),
        }
        shutil.copy2(worker_source, job / WORKER_NAME)
        input_record: dict[str, object] = {
            "schema_version": 1,
            "mode": "neo2-surface",
            "executable": str(plan.executable),
            "executable_sha256": executable_hash,
            "omp_num_threads": plan.omp_threads_per_process,
            "solver_runtime": plan.parallel_runtime,
            "surface_timeout_s": plan.surface_timeout_s,
            "max_closure_periods": plan.max_closure_periods,
            "surface": surface,
            "shot": plan.shot,
            "requested_time_s": plan.requested_time_s,
        }
        _write_json(job / JOB_INPUT_NAME, input_record)
        spec = CondorSubmitSpec(
            executable=plan.python_executable,
            initialdir=job,
            arguments=(WORKER_NAME,),
            request_cpus=plan.request_cpus,
            request_memory_mb=plan.request_memory_mb,
            should_transfer_files=plan.should_transfer_files,
            exclude_machines=plan.exclude_machines,
            accounting={
                "KAMELPhase": "neo2",
                "KAMELExecutableSha256": executable_hash[:16],
                "KAMELBoozerS": float(row[0]),
                **({"KAMELShot": plan.shot} if plan.shot is not None else {}),
                **(
                    {"KAMELRequestedTimeS": plan.requested_time_s}
                    if plan.requested_time_s is not None
                    else {}
                ),
            },
        )
        try:
            spec.write(job)
        except CondorError as error:
            raise Neo2CondorError(str(error)) from error
        staged.append(job)
    return tuple(staged)
