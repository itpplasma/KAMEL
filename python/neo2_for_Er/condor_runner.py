"""Stage NEO-2 surface jobs for bounded HTCondor execution."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import h5py
import numpy as np

from .condor import (
    CondorError,
    CondorSubmitSpec,
    CondorToolConfig,
    parse_condor_submit_output,
    query_job_ads,
    remove_clusters,
    run_condor,
    verify_shared_filesystem,
)

MANIFEST_NAME = "condor_jobs.json"
STATUS_NAME = "condor_status.json"
JOB_INPUT_NAME = "condor_job.json"
WORKER_NAME = "condor_worker.py"
_CONDOR_STATUS_NAMES = {
    1: "idle",
    2: "running",
    3: "removed",
    4: "completed",
    5: "held",
    6: "transferring_output",
    7: "suspended",
}


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
        if self.should_transfer_files != "NEVER":
            raise Neo2CondorError(
                "NEO-2 Condor currently requires should_transfer_files='NEVER'; "
                "the worker and solver runtime are not staged as a transfer bundle"
            )
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


@dataclass(frozen=True, slots=True)
class Neo2CondorResults:
    """Validated successful profile points plus an outcome for every surface."""

    profile: np.ndarray | None
    surface_records: tuple[dict[str, object], ...]


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


def _surface_fields(row: np.ndarray) -> dict[str, float]:
    return {
        "boozer_s": float(row[0]),
        "r_eff_cm": float(row[1]),
        "r_beg_cm": float(row[2]),
        "z_beg_cm": float(row[3]),
        "ti_eV": float(row[4]),
        "ne_cm3": float(row[5]),
        "kappa_cm_inv": float(row[6]),
    }


def _job_identity(root: Path, executable_hash: str, surface: dict[str, float]) -> str:
    identity_payload = {
        "work_directory": str(root.resolve()),
        "executable_sha256": executable_hash,
        "surface": surface,
    }
    return hashlib.sha256(
        json.dumps(identity_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _submit_spec(
    plan: Neo2CondorPlan,
    job: Path,
    surface: dict[str, float],
    executable_hash: str,
    job_identity: str,
) -> CondorSubmitSpec:
    return CondorSubmitSpec(
        executable=plan.python_executable,
        initialdir=job,
        arguments=(WORKER_NAME,),
        request_cpus=plan.request_cpus,
        request_memory_mb=plan.request_memory_mb,
        getenv=False,
        should_transfer_files=plan.should_transfer_files,
        exclude_machines=plan.exclude_machines,
        accounting={
            "KAMELPhase": "neo2",
            "KAMELExecutableSha256": executable_hash[:16],
            "KAMELBoozerS": surface["boozer_s"],
            "KAMELJobIdentity": job_identity,
            **({"KAMELShot": plan.shot} if plan.shot is not None else {}),
            **(
                {"KAMELRequestedTimeS": plan.requested_time_s}
                if plan.requested_time_s is not None
                else {}
            ),
        },
    )


def _write_json(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_name(f".{path.name}.{path.parent.name}.{os.getpid()}.tmp")
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
        surface = _surface_fields(row)
        job_identity = _job_identity(root, executable_hash, surface)
        shutil.copy2(worker_source, job / WORKER_NAME)
        input_record: dict[str, object] = {
            "schema_version": 1,
            "mode": "neo2-surface",
            "job_identity": job_identity,
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
        spec = _submit_spec(plan, job, surface, executable_hash, job_identity)
        try:
            spec.write(job)
        except CondorError as error:
            raise Neo2CondorError(str(error)) from error
        staged.append(job)
    return tuple(staged)


def _utc_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _json_hash(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _plan_record(plan: Neo2CondorPlan) -> dict[str, object]:
    """Return a safe, deterministic snapshot of the scheduler-relevant plan."""

    return {
        "executable": str(plan.executable),
        "executable_sha256": _sha256(plan.executable),
        "python_executable": str(plan.python_executable),
        "condor_bin_directory": str(plan.condor.condor_bin_directory),
        "command_timeout_s": plan.condor.command_timeout_s,
        "shared_filesystem_prefixes": [str(path) for path in plan.shared_filesystem_prefixes],
        "should_transfer_files": plan.should_transfer_files,
        "shot": plan.shot,
        "requested_time_s": plan.requested_time_s,
        "omp_threads_per_process": plan.omp_threads_per_process,
        "request_cpus": plan.request_cpus,
        "request_memory_mb": plan.request_memory_mb,
        "surface_timeout_s": plan.surface_timeout_s,
        "max_closure_periods": plan.max_closure_periods,
        "poll_interval_s": plan.poll_interval_s,
        "max_wall_clock_s": plan.max_wall_clock_s,
        "exclude_machines": list(plan.exclude_machines),
        "parallel_runtime": dict(plan.parallel_runtime or {}),
    }


def _staged_job_records(jobs: tuple[Path, ...], plan: Neo2CondorPlan) -> list[dict[str, object]]:
    if not jobs:
        raise Neo2CondorError("there are no staged surface jobs")
    try:
        rows = np.loadtxt(jobs[0].parent / "surfaces.dat", ndmin=2)
        executable_hash = _sha256(plan.executable)
        worker_source = Path(__file__).with_name(WORKER_NAME).read_bytes()
    except (OSError, ValueError) as error:
        raise Neo2CondorError(f"could not validate staged job inputs: {error}") from error
    if rows.shape != (len(jobs), 7) or not np.all(np.isfinite(rows)):
        raise Neo2CondorError("surfaces.dat changed after Condor payload staging")
    records: list[dict[str, object]] = []
    for job, row in zip(jobs, rows, strict=True):
        input_path = job / JOB_INPUT_NAME
        submit_path = job / "condor.submit"
        try:
            payload = json.loads(input_path.read_text(encoding="utf-8"))
            input_hash = _sha256(input_path)
            submit_hash = _sha256(submit_path)
            submit_text = submit_path.read_text(encoding="utf-8")
            worker_bytes = (job / WORKER_NAME).read_bytes()
        except (OSError, UnicodeError, json.JSONDecodeError) as error:
            raise Neo2CondorError(
                f"could not validate staged Condor files in {job}: {error}"
            ) from error
        if not isinstance(payload, dict) or payload.get("schema_version") != 1:
            raise Neo2CondorError(f"invalid staged NEO-2 job input: {input_path}")
        identity = payload.get("job_identity")
        if not isinstance(identity, str) or not re.fullmatch(r"[0-9a-f]{64}", identity):
            raise Neo2CondorError(f"staged job has no valid identity: {input_path}")
        if payload.get("mode") != "neo2-surface":
            raise Neo2CondorError(f"staged job has an unsupported mode: {input_path}")
        surface = _surface_fields(row)
        if payload.get("surface") != surface:
            raise Neo2CondorError(
                f"staged surface metadata does not match surfaces.dat: {input_path}"
            )
        if (
            payload.get("executable") != str(plan.executable)
            or payload.get("executable_sha256") != executable_hash
        ):
            raise Neo2CondorError(
                f"staged solver executable no longer matches the plan: {input_path}"
            )
        for name, expected_value in (
            ("omp_num_threads", plan.omp_threads_per_process),
            ("solver_runtime", plan.parallel_runtime),
            ("surface_timeout_s", plan.surface_timeout_s),
            ("max_closure_periods", plan.max_closure_periods),
            ("shot", plan.shot),
            ("requested_time_s", plan.requested_time_s),
        ):
            if payload.get(name) != expected_value:
                raise Neo2CondorError(f"staged job {name} does not match the requested plan")
        expected_identity = _job_identity(job.parent, executable_hash, surface)
        if identity != expected_identity:
            raise Neo2CondorError(f"staged job identity does not match its surface: {input_path}")
        if worker_bytes != worker_source:
            raise Neo2CondorError(
                f"staged worker differs from the package worker: {job / WORKER_NAME}"
            )
        expected_submit = _submit_spec(
            plan, job, surface, executable_hash, expected_identity
        ).render()
        if submit_text != expected_submit:
            raise Neo2CondorError(
                f"staged submit description does not match the plan: {submit_path}"
            )
        records.append(
            {
                "job_directory": str(job.resolve()),
                "job_identity": identity,
                "input_sha256": input_hash,
                "submit_sha256": submit_hash,
                "surface": payload.get("surface"),
                "cluster": None,
                "state": "prepared",
                "queue_status": None,
            }
        )
    if len({record["job_identity"] for record in records}) != len(records):
        raise Neo2CondorError("staged jobs contain duplicate job identities")
    return records


def _new_manifest(root: Path, plan: Neo2CondorPlan, jobs: tuple[Path, ...]) -> dict[str, object]:
    plan_record = _plan_record(plan)
    return {
        "schema_version": 1,
        "kind": "neo2-condor-manifest",
        "work_directory": str(root.resolve()),
        "created_at_utc": _utc_now(),
        "updated_at_utc": _utc_now(),
        "plan": plan_record,
        "plan_fingerprint": _json_hash(plan_record),
        "submission_complete": False,
        "jobs": _staged_job_records(jobs, plan),
    }


def _load_manifest(root: Path, plan: Neo2CondorPlan, jobs: tuple[Path, ...]) -> dict[str, object]:
    path = root / MANIFEST_NAME
    if not path.exists():
        return _new_manifest(root, plan, jobs)
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise Neo2CondorError(f"could not read Condor manifest {path}: {error}") from error
    if (
        not isinstance(document, dict)
        or document.get("schema_version") != 1
        or document.get("kind") != "neo2-condor-manifest"
        or not isinstance(document.get("jobs"), list)
    ):
        raise Neo2CondorError(f"unsupported or malformed Condor manifest: {path}")
    if document.get("work_directory") != str(root.resolve()):
        raise Neo2CondorError("Condor manifest belongs to a different work directory")
    expected_plan = _plan_record(plan)
    if document.get("plan_fingerprint") != _json_hash(expected_plan):
        raise Neo2CondorError("Condor manifest plan does not match the requested plan")
    expected_jobs = _staged_job_records(jobs, plan)
    actual_jobs = document["jobs"]
    if len(actual_jobs) != len(expected_jobs):
        raise Neo2CondorError("Condor manifest job count does not match staged surfaces")
    for index, (actual, expected) in enumerate(zip(actual_jobs, expected_jobs, strict=True)):
        if not isinstance(actual, dict):
            raise Neo2CondorError(f"Condor manifest job {index} is malformed")
        for field_name in ("job_directory", "job_identity", "input_sha256", "submit_sha256"):
            if actual.get(field_name) != expected[field_name]:
                raise Neo2CondorError(
                    f"Condor manifest job {index} no longer matches staged {field_name}"
                )
        if actual.get("state") not in {
            "prepared",
            "submitting",
            "submitted",
            "adopted",
            "ambiguous",
        }:
            raise Neo2CondorError(f"Condor manifest job {index} has an unknown state")
        cluster = actual.get("cluster")
        if cluster is not None and (
            isinstance(cluster, bool) or not isinstance(cluster, int) or cluster <= 0
        ):
            raise Neo2CondorError(f"Condor manifest job {index} has an invalid cluster id")
        if actual.get("state") in {"submitted", "adopted"} and cluster is None:
            raise Neo2CondorError(f"Condor manifest job {index} is submitted without a cluster id")
        if actual.get("state") in {"prepared", "submitting"} and cluster is not None:
            raise Neo2CondorError(
                f"Condor manifest job {index} has a cluster before acknowledgement"
            )
        # Keep the current validated surface payload; scheduler state comes from the manifest.
        for field_name in ("surface", "queue_status"):
            if field_name not in actual:
                actual[field_name] = expected[field_name]
    return document


def _save_manifest(root: Path, manifest: dict[str, object]) -> None:
    manifest["updated_at_utc"] = _utc_now()
    _write_json(root / MANIFEST_NAME, manifest)


def _ambiguous_jobs(manifest: dict[str, object]) -> list[str]:
    return [
        str(job.get("job_directory"))
        for job in manifest["jobs"]
        if isinstance(job, dict) and job.get("state") in {"ambiguous", "submitting"}
    ]


def _record_by_cluster(
    ads: tuple[dict[str, object], ...] | list[dict[str, object]],
) -> dict[int, dict[str, object]]:
    """Prefer live queue ads over history when both describe one cluster."""

    result: dict[int, dict[str, object]] = {}
    for ad in ads:
        cluster = ad.get("cluster")
        if isinstance(cluster, bool) or not isinstance(cluster, int):
            continue
        current = result.get(cluster)
        if current is None or (
            current.get("ad_source") != "condor_q" and ad.get("ad_source") == "condor_q"
        ):
            result[cluster] = ad
    return result


def _root_and_jobs(
    work_directory: str | Path, *, plan: Neo2CondorPlan
) -> tuple[Path, tuple[Path, ...]]:
    root = Path(work_directory).expanduser().absolute()
    try:
        verify_shared_filesystem(
            root,
            allowed_prefixes=plan.shared_filesystem_prefixes,
            should_transfer_files=plan.should_transfer_files,
        )
    except CondorError as error:
        raise Neo2CondorError(str(error)) from error
    jobs, _rows = _staged_surface_jobs(root)
    for job in jobs:
        if not all(
            (job / name).is_file() for name in (JOB_INPUT_NAME, WORKER_NAME, "condor.submit")
        ):
            raise Neo2CondorError(f"Condor payload is incomplete in {job}")
    return root, jobs


def submit_neo2_condor_jobs(
    work_directory: str | Path,
    *,
    plan: Neo2CondorPlan,
    dry_run: bool = False,
    adopt_existing: bool = True,
) -> dict[str, object]:
    """Submit each staged surface once, durably recording ambiguity and cluster ids."""

    if not isinstance(dry_run, bool) or not isinstance(adopt_existing, bool):
        raise Neo2CondorError("dry_run and adopt_existing must be booleans")
    root, jobs = _root_and_jobs(work_directory, plan=plan)
    manifest = _load_manifest(root, plan, jobs)
    if dry_run:
        preview = json.loads(json.dumps(manifest))
        preview["dry_run"] = True
        preview["jobs"] = [
            {
                **job,
                "cluster": None,
                "state": "dry_run",
            }
            for job in preview["jobs"]
        ]
        return preview

    ambiguous = _ambiguous_jobs(manifest)
    if ambiguous:
        for job in manifest["jobs"]:
            if isinstance(job, dict) and job.get("state") == "submitting":
                job["state"] = "ambiguous"
                job["error"] = (
                    "previous process ended during submission; acknowledgement is unknown"
                )
        _save_manifest(root, manifest)
        raise Neo2CondorError(
            f"ambiguous Condor submission for {ambiguous}; refusing automatic resubmission"
        )

    known_clusters = [
        int(job["cluster"])
        for job in manifest["jobs"]
        if isinstance(job, dict)
        and job.get("state") in {"submitted", "adopted"}
        and job.get("cluster") is not None
    ]
    if adopt_existing and known_clusters:
        try:
            ads = query_job_ads(known_clusters, config=plan.condor, include_history=True)
        except CondorError as error:
            manifest["last_reconciliation_error"] = str(error)
            _save_manifest(root, manifest)
            raise Neo2CondorError(f"could not reconcile known Condor jobs: {error}") from error
        by_cluster = _record_by_cluster(ads)
        for job in manifest["jobs"]:
            if not isinstance(job, dict) or job.get("cluster") is None:
                continue
            ad = by_cluster.get(int(job["cluster"]))
            if ad is None:
                job["queue_status"] = {"job_status": "missing", "ad_source": None}
                continue
            if ad.get("proc") != 0 or ad.get("job_identity") != job.get("job_identity"):
                job["state"] = "ambiguous"
                job["error"] = "scheduler ad identity does not match the staged surface"
                _save_manifest(root, manifest)
                raise Neo2CondorError(
                    f"Condor cluster {job['cluster']} identity does not match staged job; "
                    "refusing to continue"
                )
            job["state"] = "adopted"
            job["queue_status"] = {
                key: ad.get(key)
                for key in (
                    "job_status_code",
                    "job_status",
                    "exit_code",
                    "remote_host",
                    "remote_wall_clock_s",
                    "num_holds",
                    "hold_reason",
                    "ad_source",
                )
            }
        _save_manifest(root, manifest)

    for job in manifest["jobs"]:
        if not isinstance(job, dict):
            raise Neo2CondorError("Condor manifest contains a malformed job record")
        if job.get("cluster") is not None:
            continue
        if job.get("state") != "prepared":
            raise Neo2CondorError(
                f"refusing to submit job {job.get('job_directory')} in state {job.get('state')!r}"
            )
        job["state"] = "submitting"
        job["submit_started_at_utc"] = _utc_now()
        _save_manifest(root, manifest)
        submit_file = Path(str(job["job_directory"])) / "condor.submit"
        try:
            result = run_condor([str(submit_file)], config=plan.condor, tool="condor_submit")
            if result.returncode != 0:
                raise CondorError(
                    f"condor_submit exited {result.returncode}: {result.stderr.strip()}"
                )
            cluster = parse_condor_submit_output(f"{result.stdout}\n{result.stderr}")
        except CondorError as error:
            job["state"] = "ambiguous"
            job["error"] = str(error)
            job["submit_finished_at_utc"] = _utc_now()
            _save_manifest(root, manifest)
            raise Neo2CondorError(
                f"ambiguous Condor submission for {job['job_directory']}: {error}; "
                "refusing automatic retry"
            ) from error
        job["cluster"] = cluster
        job["state"] = "submitted"
        job["submitted_at_utc"] = _utc_now()
        job["submit_acknowledgement"] = (result.stdout + result.stderr).strip()
        _save_manifest(root, manifest)

    manifest["submission_complete"] = all(
        isinstance(job, dict)
        and job.get("cluster") is not None
        and job.get("state") in {"submitted", "adopted"}
        for job in manifest["jobs"]
    )
    _save_manifest(root, manifest)
    result = json.loads(json.dumps(manifest))
    result["dry_run"] = False
    result["adopted_clusters"] = sorted(
        int(job["cluster"])
        for job in manifest["jobs"]
        if isinstance(job, dict)
        and job.get("state") == "adopted"
        and job.get("cluster") is not None
    )
    return result


def _status_record(job: dict[str, object], ad: dict[str, object] | None) -> dict[str, object]:
    record: dict[str, object] = {
        "job_directory": job.get("job_directory"),
        "job_identity": job.get("job_identity"),
        "cluster": job.get("cluster"),
        "surface": job.get("surface"),
        "job_status": "missing",
        "job_status_code": None,
        "exit_code": None,
        "remote_host": None,
        "remote_wall_clock_s": None,
        "num_holds": None,
        "hold_reason": None,
        "ad_source": None,
        "failure_kind": "condor_missing",
        "failure_reason": "no queue or history ad found for the submitted cluster",
        "terminal": True,
    }
    if ad is None:
        return record
    if ad.get("proc") != 0 or ad.get("job_identity") != job.get("job_identity"):
        record.update(
            {
                "job_status": "identity_mismatch",
                "failure_kind": "condor_identity_mismatch",
                "failure_reason": "scheduler ad identity does not match the staged surface",
            }
        )
        return record
    code = ad.get("job_status_code")
    if (
        isinstance(code, bool)
        or not isinstance(code, int)
        or code not in _CONDOR_STATUS_NAMES
        or ad.get("job_status") != _CONDOR_STATUS_NAMES[code]
    ):
        record.update(
            {
                "job_status": "invalid",
                "failure_kind": "condor_invalid_status",
                "failure_reason": f"unsupported scheduler status code: {code!r}",
            }
        )
        return record
    status = str(ad.get("job_status"))
    exit_code = ad.get("exit_code")
    terminal = code in {3, 4, 5}
    record.update(
        {
            "job_status": status,
            "job_status_code": code,
            "exit_code": exit_code,
            "remote_host": ad.get("remote_host"),
            "remote_wall_clock_s": ad.get("remote_wall_clock_s"),
            "num_holds": ad.get("num_holds"),
            "hold_reason": ad.get("hold_reason"),
            "ad_source": ad.get("ad_source"),
            "failure_kind": None,
            "failure_reason": None,
            "terminal": terminal,
        }
    )
    if code == 3:
        record["failure_kind"] = "condor_removed"
        record["failure_reason"] = ad.get("hold_reason") or "job was removed from the pool"
    elif code == 5:
        record["failure_kind"] = "condor_held"
        record["failure_reason"] = ad.get("hold_reason") or "job is held"
    elif code == 4 and exit_code != 0:
        record["failure_kind"] = (
            "condor_nonzero_exit" if exit_code is not None else "condor_exit_unknown"
        )
        record["failure_reason"] = (
            f"job completed with exit code {exit_code}"
            if exit_code is not None
            else "job completed but the scheduler did not report an exit code"
        )
    return record


def _write_status(root: Path, document: dict[str, object]) -> None:
    document["updated_at_utc"] = _utc_now()
    _write_json(root / STATUS_NAME, document)


def wait_neo2_condor_jobs(
    work_directory: str | Path,
    *,
    plan: Neo2CondorPlan,
) -> tuple[dict[str, object], ...]:
    """Poll known clusters to terminal states, removing pending jobs at the deadline."""

    root, jobs = _root_and_jobs(work_directory, plan=plan)
    manifest = _load_manifest(root, plan, jobs)
    if _ambiguous_jobs(manifest):
        raise Neo2CondorError("cannot wait for jobs with ambiguous submission state")
    records = [job for job in manifest["jobs"] if isinstance(job, dict)]
    if len(records) != len(jobs) or any(job.get("cluster") is None for job in records):
        raise Neo2CondorError("all surfaces must have an acknowledged cluster before waiting")
    clusters = [int(job["cluster"]) for job in records]
    started = time.monotonic()
    deadline = started + plan.max_wall_clock_s
    status_document: dict[str, object] = {
        "schema_version": 1,
        "kind": "neo2-condor-status",
        "work_directory": str(root.resolve()),
        "plan_fingerprint": manifest["plan_fingerprint"],
        "started_at_utc": _utc_now(),
        "deadline_s": plan.max_wall_clock_s,
        "clusters": clusters,
        "driver_removed_clusters": [],
        "statuses": [],
        "state": "polling",
    }
    latest: list[dict[str, object]] = []
    while True:
        try:
            ads = query_job_ads(clusters, config=plan.condor, include_history=True)
        except CondorError as error:
            status_document["state"] = "query_error"
            status_document["query_error"] = str(error)
            if latest:
                for status in latest:
                    status["query_error"] = str(error)
            else:
                latest = [
                    {
                        "job_directory": job.get("job_directory"),
                        "job_identity": job.get("job_identity"),
                        "cluster": job.get("cluster"),
                        "surface": job.get("surface"),
                        "job_status": "query_error",
                        "job_status_code": None,
                        "exit_code": None,
                        "failure_kind": "condor_query_error",
                        "failure_reason": str(error),
                        "query_error": str(error),
                        "terminal": False,
                    }
                    for job in records
                ]
            status_document["statuses"] = latest
            _write_status(root, status_document)
            raise Neo2CondorError(f"Condor status query failed: {error}") from error
        by_cluster = _record_by_cluster(ads)
        latest = [_status_record(job, by_cluster.get(int(job["cluster"]))) for job in records]
        for job, status in zip(records, latest, strict=True):
            job["queue_status"] = {
                key: status.get(key)
                for key in (
                    "job_status",
                    "job_status_code",
                    "exit_code",
                    "remote_host",
                    "remote_wall_clock_s",
                    "num_holds",
                    "hold_reason",
                    "failure_kind",
                )
            }
        _save_manifest(root, manifest)
        status_document["statuses"] = latest
        pending = [status for status in latest if not status["terminal"]]
        status_document["state"] = (
            "polling"
            if pending
            else (
                "completed"
                if all(status["failure_kind"] is None for status in latest)
                else "failed"
            )
        )
        _write_status(root, status_document)
        if not pending:
            break
        now = time.monotonic()
        if now >= deadline:
            pending_clusters = sorted(
                {int(status["cluster"]) for status in pending if status.get("cluster") is not None}
            )
            status_document["state"] = "driver_deadline_removal_requested"
            status_document["driver_removal_requested_clusters"] = pending_clusters
            status_document["removal_reason"] = "KAMEL NEO-2 Condor wall-clock deadline exceeded"
            for status in pending:
                status["deadline_exceeded"] = True
            _write_status(root, status_document)
            try:
                remove_clusters(
                    pending_clusters,
                    config=plan.condor,
                    reason="KAMEL NEO-2 Condor wall-clock deadline exceeded",
                )
            except CondorError as error:
                status_document["state"] = "driver_deadline_removal_error"
                status_document["removal_error"] = str(error)
                for status in pending:
                    status["failure_kind"] = "condor_driver_deadline_remove_error"
                    status["failure_reason"] = str(error)
                    status["driver_removal_requested"] = True
                for job, status in zip(records, latest, strict=True):
                    if status.get("failure_kind") == "condor_driver_deadline_remove_error":
                        job["queue_status"] = {
                            "job_status": status["job_status"],
                            "failure_kind": status["failure_kind"],
                        }
                _save_manifest(root, manifest)
                _write_status(root, status_document)
                break
            for status in pending:
                status["job_status"] = "removed"
                status["job_status_code"] = 3
                status["failure_kind"] = "condor_driver_wall_clock"
                status["failure_reason"] = "removed after KAMEL driver wall-clock deadline"
                status["terminal"] = True
                status["driver_removal_requested"] = True
            status_document["driver_removed_clusters"] = pending_clusters
            for job, status in zip(records, latest, strict=True):
                job["queue_status"] = {
                    "job_status": status["job_status"],
                    "job_status_code": status["job_status_code"],
                    "failure_kind": status["failure_kind"],
                }
            _save_manifest(root, manifest)
            status_document["state"] = "driver_deadline_removed"
            _write_status(root, status_document)
            break
        time.sleep(min(plan.poll_interval_s, max(0.0, deadline - now)))

    status_document["finished_at_utc"] = _utc_now()
    status_document["wall_clock_s"] = time.monotonic() - started
    if status_document["state"] == "polling":
        status_document["state"] = (
            "completed" if all(status["failure_kind"] is None for status in latest) else "failed"
        )
    _write_status(root, status_document)
    return tuple(latest)


def _read_json_object(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError(f"could not read {path.name}: {error}") from error
    if not isinstance(payload, dict):
        raise ValueError(f"{path.name} must contain a JSON object")
    return payload


def _surface_failure(
    record: dict[str, object], kind: str, reason: str, *, status: str = "failed"
) -> None:
    record["status"] = status
    record["failure_kind"] = kind
    record["failure_reason"] = reason


def _load_condor_statuses(root: Path, manifest: dict[str, object]) -> dict[str, dict[str, object]]:
    status_path = root / STATUS_NAME
    try:
        document = _read_json_object(status_path)
    except ValueError:
        return {}
    if (
        document.get("schema_version") != 1
        or document.get("kind") != "neo2-condor-status"
        or document.get("work_directory") != str(root.resolve())
        or document.get("plan_fingerprint") != manifest.get("plan_fingerprint")
    ):
        return {}
    statuses = document.get("statuses")
    if not isinstance(statuses, list):
        return {}
    result: dict[str, dict[str, object]] = {}
    for status in statuses:
        if not isinstance(status, dict):
            continue
        identity = status.get("job_identity")
        if isinstance(identity, str):
            if identity in result:
                return {}
            result[identity] = status
    return result


def _single_finite_real(value: object) -> float:
    array = np.asarray(value)
    if array.size != 1 or np.iscomplexobj(array):
        raise ValueError("expected one real scalar")
    result = float(array.reshape(-1)[0])
    if not np.isfinite(result):
        raise ValueError("value is not finite")
    return result


def _validate_hdf5_outputs(job: Path, expected_boozer_s: float) -> float:
    try:
        with h5py.File(job / "neo2_config.h5", "r") as config:
            boozer_s = _single_finite_real(config["settings/boozer_s"][()])
        with h5py.File(job / "fulltransp.h5", "r") as transport:
            k_cof = _single_finite_real(transport["k_cof"][()])
    except (OSError, KeyError, TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"invalid NEO-2 HDF5 output: {error}") from error
    if not np.isclose(boozer_s, expected_boozer_s):
        raise ValueError(
            f"neo2_config.h5 settings/boozer_s={boozer_s} does not match "
            f"the requested surface {expected_boozer_s}"
        )
    return k_cof


def _validate_worker_record(
    record: dict[str, object],
    job_input: dict[str, object],
    plan: Neo2CondorPlan,
    executable_hash: str,
) -> tuple[str | None, str | None]:
    if record.get("schema_version") != 1 or record.get("mode") != "neo2-surface":
        return "worker_provenance_mismatch", "unsupported worker record schema or mode"
    provenance = {
        "job_identity": job_input.get("job_identity"),
        "executable": str(plan.executable),
        "executable_sha256": executable_hash,
        "omp_num_threads": str(plan.omp_threads_per_process),
        "surface": job_input.get("surface"),
    }
    for field_name, expected in provenance.items():
        if record.get(field_name) != expected:
            return (
                "worker_provenance_mismatch",
                f"worker {field_name} does not match the staged job and plan",
            )
    host = record.get("host")
    if not isinstance(host, str) or not host.strip():
        return "worker_provenance_mismatch", "worker record has no execute host"
    closure_period = record.get("closure_period")
    if closure_period is not None and (
        isinstance(closure_period, bool)
        or not isinstance(closure_period, int)
        or closure_period < 0
    ):
        return "worker_provenance_mismatch", "worker closure_period is not a non-negative integer"
    if closure_period is not None and closure_period > plan.max_closure_periods:
        return (
            "closure_period_exceeded",
            f"reported closure period {closure_period} exceeds configured limit "
            f"{plan.max_closure_periods}",
        )
    status = record.get("status")
    exit_code = record.get("exit_code")
    valid_exit_code = not isinstance(exit_code, bool) and isinstance(exit_code, int)
    if status != "succeeded" or not valid_exit_code or exit_code != 0:
        if status == "timeout" or exit_code == 124:
            return "worker_timeout", "NEO-2 worker exceeded its per-surface timeout"
        if status == "launch_failure" or exit_code == 127:
            return "worker_launch_failure", "NEO-2 solver could not be launched"
        if status == "nonzero_exit" and valid_exit_code:
            return "worker_nonzero_exit", f"NEO-2 solver exited with code {exit_code}"
        return (
            "worker_failed",
            f"worker did not report success (status={status!r}, exit={exit_code!r})",
        )
    return None, None


def collect_neo2_condor_results(
    work_directory: str | Path,
    *,
    plan: Neo2CondorPlan,
) -> Neo2CondorResults:
    """Validate each scheduler/worker/HDF5 outcome and retain partial successes."""

    root, jobs = _root_and_jobs(work_directory, plan=plan)
    manifest = _load_manifest(root, plan, jobs)
    statuses = _load_condor_statuses(root, manifest)
    executable_hash = _sha256(plan.executable)
    successful_points: list[tuple[float, float]] = []
    surface_records: list[dict[str, object]] = []

    for job in manifest["jobs"]:
        if not isinstance(job, dict):
            raise Neo2CondorError("Condor manifest contains a malformed surface record")
        job_directory = Path(str(job["job_directory"]))
        identity = str(job["job_identity"])
        surface = job.get("surface")
        base: dict[str, object] = {
            "job_directory": str(job_directory),
            "job_identity": identity,
            "cluster": job.get("cluster"),
            "surface": surface,
            "r_eff_cm": surface.get("r_eff_cm") if isinstance(surface, dict) else None,
            "boozer_s": surface.get("boozer_s") if isinstance(surface, dict) else None,
            "status": "failed",
            "failure_kind": None,
            "failure_reason": None,
        }
        queue_status = statuses.get(identity)
        if queue_status is None:
            _surface_failure(
                base,
                "condor_status_missing",
                "no matching per-surface record in condor_status.json",
            )
            surface_records.append(base)
            continue
        base["condor_status"] = queue_status
        if queue_status.get("cluster") != job.get("cluster"):
            _surface_failure(
                base,
                "condor_status_identity_mismatch",
                "Condor status cluster does not match the submitted surface cluster",
            )
            surface_records.append(base)
            continue
        if queue_status.get("terminal") is not True:
            _surface_failure(
                base,
                "condor_pending",
                f"Condor job is not terminal (state={queue_status.get('job_status')!r})",
                status="pending",
            )
            surface_records.append(base)
            continue
        if queue_status.get("failure_kind") is not None:
            _surface_failure(
                base,
                str(queue_status.get("failure_kind")),
                str(
                    queue_status.get("failure_reason") or "Condor job did not complete successfully"
                ),
            )
            surface_records.append(base)
            continue
        if queue_status.get("job_status") != "completed" or queue_status.get("exit_code") != 0:
            _surface_failure(
                base,
                "condor_not_successful",
                "Condor status is not completed with exit code zero",
            )
            surface_records.append(base)
            continue

        try:
            job_input = _read_json_object(job_directory / JOB_INPUT_NAME)
            worker_record = _read_json_object(job_directory / "condor_run_record.json")
        except ValueError as error:
            kind = (
                "worker_record_missing"
                if not (job_directory / "condor_run_record.json").exists()
                else "worker_record_invalid"
            )
            _surface_failure(base, kind, str(error))
            surface_records.append(base)
            continue
        base["worker_host"] = worker_record.get("host")
        base["closure_period"] = worker_record.get("closure_period")
        worker_failure, worker_reason = _validate_worker_record(
            worker_record, job_input, plan, executable_hash
        )
        if worker_failure is not None:
            _surface_failure(base, worker_failure, str(worker_reason))
            base["worker_record"] = worker_record
            surface_records.append(base)
            continue
        try:
            expected_boozer_s = float(surface["boozer_s"])
            k_cof = _validate_hdf5_outputs(job_directory, expected_boozer_s)
            r_eff_cm = float(surface["r_eff_cm"])
            if not np.isfinite(r_eff_cm):
                raise ValueError("requested r_eff_cm is not finite")
        except (KeyError, TypeError, ValueError, OverflowError) as error:
            _surface_failure(base, "invalid_hdf5_output", str(error))
            surface_records.append(base)
            continue
        base["status"] = "succeeded"
        base["k_cof"] = k_cof
        successful_points.append((r_eff_cm, k_cof))
        surface_records.append(base)

    profile = (
        np.asarray(sorted(successful_points, key=lambda point: point[0]), dtype=float)
        if successful_points
        else None
    )
    return Neo2CondorResults(profile=profile, surface_records=tuple(surface_records))
