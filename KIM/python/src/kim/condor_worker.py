#!/usr/bin/env python3
"""Standard-library Condor worker for one staged KIM scan point."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import re
import signal
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

EXIT_SUCCESS = 0
EXIT_SOLVER_ERROR = 1
EXIT_INVALID_JOB = 2
EXIT_TIMEOUT = 124
EXIT_LAUNCH_FAILURE = 127
_PROCESS_GROUP_TERM_GRACE_S = 5.0
_JOB_INPUT = "condor_job.json"
_RUN_RECORD = "condor_run_record.json"


class InvalidJobError(ValueError):
    """The staged job input is missing, malformed, or unsupported."""


class WorkerLaunchError(RuntimeError):
    """A staged solver or optional API dependency could not be launched."""


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _read_job(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as error:
        raise InvalidJobError(f"staged job input is missing: {path}") from error
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise InvalidJobError(f"could not read staged job input {path}: {error}") from error
    if not isinstance(payload, dict):
        raise InvalidJobError("staged job input must be a JSON object")
    if payload.get("schema_version") != 1:
        raise InvalidJobError("unsupported KIM Condor job schema version")
    backend = payload.get("backend")
    expected_mode = {"kamel_kim_python": "kim-sweep", "kim_x_namelist": "kim-run"}.get(backend)
    if expected_mode is None or payload.get("mode") != expected_mode:
        raise InvalidJobError("job mode does not match an explicit supported backend")
    identity = payload.get("job_identity")
    if not isinstance(identity, str) or not re.fullmatch(r"[0-9a-f]{64}", identity):
        raise InvalidJobError("job_identity must be a lowercase SHA-256 digest")
    executable = payload.get("executable")
    if not isinstance(executable, str) or not Path(executable).is_absolute():
        raise InvalidJobError("executable must be an absolute path")
    executable_hash = payload.get("executable_sha256")
    if not isinstance(executable_hash, str) or not re.fullmatch(r"[0-9a-f]{64}", executable_hash):
        raise InvalidJobError("executable_sha256 must be a lowercase SHA-256 digest")
    factor = payload.get("profile_scale_factor")
    if (
        isinstance(factor, bool)
        or not isinstance(factor, (int, float))
        or not math.isfinite(factor)
    ):
        raise InvalidJobError("profile_scale_factor must be finite")
    index = payload.get("scan_order_index")
    if isinstance(index, bool) or not isinstance(index, int) or index < 0:
        raise InvalidJobError("scan_order_index must be a non-negative integer")
    timeout = payload.get("job_timeout_s")
    if isinstance(timeout, bool) or not isinstance(timeout, (int, float)):
        raise InvalidJobError("job_timeout_s must be a positive finite number")
    if not math.isfinite(timeout) or timeout <= 0.0:
        raise InvalidJobError("job_timeout_s must be a positive finite number")
    omp_threads = payload.get("omp_num_threads")
    if isinstance(omp_threads, bool) or not isinstance(omp_threads, int) or omp_threads <= 0:
        raise InvalidJobError("omp_num_threads must be a positive integer")
    if backend == "kamel_kim_python":
        if not isinstance(payload.get("base_config"), dict):
            raise InvalidJobError("API job input must contain a base_config object")
        component = payload.get("api_current_component")
        unit = payload.get("api_current_unit")
        if (component is None) != (unit is None):
            raise InvalidJobError("API current component and unit must be supplied together")
        if component not in {None, "real", "imag", "magnitude"}:
            raise InvalidJobError("unsupported API current component")
        if unit is not None and (not isinstance(unit, str) or not unit.strip()):
            raise InvalidJobError("API current unit must be a non-empty string")
    else:
        name = payload.get("namelist_file")
        if not isinstance(name, str) or Path(name).name != name or name in {".", ".."}:
            raise InvalidJobError("namelist_file must be a local filename")
        if not (path.parent / name).is_file():
            raise InvalidJobError(f"staged KIM namelist is missing: {name}")
    source_path = payload.get("kim_source_path")
    if source_path is not None and (
        not isinstance(source_path, str) or not Path(source_path).is_absolute()
    ):
        raise InvalidJobError("kim_source_path must be null or an absolute path")
    return payload


def _load_api_runner():
    """Import the optional KIM API only for API-backed jobs."""

    from kim.config import ElectrostaticPeriodicRun, SimulationConfig
    from kim.sweep import ProfileScale, SweepSpec, run_sweep

    return SimulationConfig, ElectrostaticPeriodicRun, ProfileScale, SweepSpec, run_sweep


def _run_api(job: dict[str, Any], directory: Path) -> dict[str, object]:
    source_path = job.get("kim_source_path")
    if source_path is not None:
        sys.path.insert(0, str(source_path))
    (
        SimulationConfig,
        ElectrostaticPeriodicRun,
        ProfileScale,
        SweepSpec,
        run_sweep,
    ) = _load_api_runner()
    config = SimulationConfig.model_validate(job["base_config"])
    _verify_executable(job)
    spec = SweepSpec(
        base=config,
        variation=ProfileScale(values=(float(job["profile_scale_factor"]),)),
        continue_on_failure=False,
    )
    environment = {"OMP_NUM_THREADS": str(job["omp_num_threads"])}
    result = run_sweep(
        spec,
        executable=job["executable"],
        runs_directory=directory / "kim-runs",
        label=str(job.get("job", "condor-job")),
        timeout=float(job["job_timeout_s"]),
        environment=environment,
    )
    children = getattr(result, "children", ())
    if len(children) != 1:
        raise RuntimeError(f"KIM API sweep returned {len(children)} children; expected one")
    child = children[0]
    status = getattr(child, "status", None)
    status_name = getattr(status, "value", status)
    if status_name == "timed_out":
        return {
            "status": "timeout",
            "exit_code": EXIT_TIMEOUT,
            "failure_reason": f"KIM API child run exceeded {job['job_timeout_s']} s",
        }
    if status_name != "succeeded":
        raise RuntimeError(f"KIM API child run ended with status {status_name!r}")
    record: dict[str, object] = {
        "status": "succeeded",
        "exit_code": EXIT_SUCCESS,
        "api_current_component": job.get("api_current_component"),
        "api_current_unit": job.get("api_current_unit"),
    }
    if isinstance(config.run, ElectrostaticPeriodicRun):
        current = child.periodic.integrated_parallel_current(region="as_is")
        current = complex(current)
        if not math.isfinite(current.real) or not math.isfinite(current.imag):
            raise RuntimeError("KIM API integrated parallel current is non-finite")
        record["integrated_parallel_current_real"] = float(current.real)
        record["integrated_parallel_current_imag"] = float(current.imag)
        component = job.get("api_current_component")
        if component == "real":
            scalar = current.real
        elif component == "imag":
            scalar = current.imag
        elif component == "magnitude":
            scalar = abs(current)
        else:
            scalar = None
        record["api_current_value"] = float(scalar) if scalar is not None else None
    return record


def _terminate_process_group(process: subprocess.Popen[str]) -> int:
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return int(process.wait())
    deadline = time.monotonic() + _PROCESS_GROUP_TERM_GRACE_S
    while time.monotonic() < deadline:
        process.poll()
        try:
            os.killpg(process.pid, 0)
        except ProcessLookupError:
            return int(process.wait())
        except PermissionError:
            pass
        time.sleep(min(0.05, max(0.0, deadline - time.monotonic())))
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    return int(process.wait())


def _run_namelist(job: dict[str, Any], directory: Path) -> dict[str, object]:
    _verify_executable(job)
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(job["omp_num_threads"])
    command = [str(job["executable"])]
    try:
        process = subprocess.Popen(
            command,
            cwd=str(directory),
            env=environment,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            start_new_session=True,
        )
    except OSError as error:
        raise WorkerLaunchError(f"could not launch KIM executable: {error}") from error
    try:
        output, _ = process.communicate(timeout=float(job["job_timeout_s"]))
    except subprocess.TimeoutExpired:
        _terminate_process_group(process)
        output, _ = process.communicate()
        _emit_solver_output(output)
        return {
            "status": "timeout",
            "exit_code": EXIT_TIMEOUT,
            "solver_exit_code": process.returncode,
            "failure_reason": f"KIM executable exceeded {job['job_timeout_s']} s",
        }
    _emit_solver_output(output)
    code = int(process.returncode)
    if code != 0:
        return {
            "status": "solver_error",
            "exit_code": code,
            "failure_reason": f"KIM executable exited with code {code}",
        }
    return {"status": "succeeded", "exit_code": EXIT_SUCCESS}


def _emit_solver_output(output: str | None) -> None:
    if output:
        sys.stdout.write(output)
        if not output.endswith("\n"):
            sys.stdout.write("\n")
        sys.stdout.flush()


def _verify_executable(job: dict[str, Any]) -> None:
    executable = Path(job["executable"])
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise WorkerLaunchError(f"KIM executable is missing or not executable: {executable}")
    digest = hashlib.sha256()
    try:
        with executable.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    except OSError as error:
        raise WorkerLaunchError(f"could not read KIM executable {executable}: {error}") from error
    if digest.hexdigest() != job["executable_sha256"]:
        raise InvalidJobError("KIM executable content does not match staged executable_sha256")


def _record_payload(
    job: dict[str, Any] | None,
    *,
    result: dict[str, object],
    started_at: str,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "schema_version": 1,
        "job_identity": job.get("job_identity") if job else None,
        "backend": job.get("backend") if job else None,
        "scan_order_index": job.get("scan_order_index") if job else None,
        "profile_scale_factor": job.get("profile_scale_factor") if job else None,
        "executable": job.get("executable") if job else None,
        "executable_sha256": job.get("executable_sha256") if job else None,
        "host": platform.node() or "unknown",
        "python_version": platform.python_version(),
        "started_at_utc": started_at,
        "finished_at_utc": _utc_now(),
    }
    payload.update(result)
    return payload


def _write_record(path: Path, payload: dict[str, object]) -> None:
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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job-input", default=_JOB_INPUT)
    arguments = parser.parse_args(argv)
    directory = Path.cwd()
    started_at = _utc_now()
    job: dict[str, Any] | None = None
    try:
        job_path = directory / arguments.job_input
        if Path(arguments.job_input).name != arguments.job_input:
            raise InvalidJobError("job input must be a local filename")
        job = _read_job(job_path)
        if job["backend"] == "kamel_kim_python":
            result = _run_api(job, directory)
        else:
            result = _run_namelist(job, directory)
        code = int(result["exit_code"])
    except InvalidJobError as error:
        result = {
            "status": "invalid_job",
            "exit_code": EXIT_INVALID_JOB,
            "failure_reason": str(error),
        }
        code = EXIT_INVALID_JOB
    except ImportError as error:
        result = {
            "status": "launch_failure",
            "exit_code": EXIT_LAUNCH_FAILURE,
            "failure_reason": str(error),
        }
        code = EXIT_LAUNCH_FAILURE
    except WorkerLaunchError as error:
        result = {
            "status": "launch_failure",
            "exit_code": EXIT_LAUNCH_FAILURE,
            "failure_reason": str(error),
        }
        code = EXIT_LAUNCH_FAILURE
    except Exception as error:
        result = {
            "status": "solver_error",
            "exit_code": EXIT_SOLVER_ERROR,
            "failure_reason": f"{type(error).__name__}: {error}",
        }
        code = EXIT_SOLVER_ERROR
    record = _record_payload(job, result=result, started_at=started_at)
    try:
        _write_record(directory / _RUN_RECORD, record)
    except (OSError, TypeError, ValueError) as error:
        print(f"could not persist KIM worker record: {error}", file=sys.stderr)
        return EXIT_SOLVER_ERROR
    return code


if __name__ == "__main__":
    raise SystemExit(main())
