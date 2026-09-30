#!/usr/bin/env python3
"""Standard-library HTCondor payload for one staged NEO-2 surface."""

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
import threading
import time
from pathlib import Path
from typing import TextIO

EXIT_TIMEOUT = 124
EXIT_LAUNCH_FAILURE = 127

_PERIOD_RE = re.compile(r"^\s*period:\s*(?P<period>\d+)\s*$", re.MULTILINE)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_job(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError(f"could not read job input {path}: {error}") from error
    if not isinstance(payload, dict):
        raise ValueError(f"job input must contain a JSON object: {path}")
    if payload.get("schema_version") != 1 or payload.get("mode") != "neo2-surface":
        raise ValueError(f"unsupported NEO-2 job input schema or mode: {path}")
    executable = payload.get("executable")
    if not isinstance(executable, str) or not Path(executable).is_absolute():
        raise ValueError("job input executable must be an absolute path")
    sha256 = payload.get("executable_sha256")
    if not isinstance(sha256, str) or not re.fullmatch(r"[0-9a-f]{64}", sha256):
        raise ValueError("job input executable_sha256 must be a lowercase SHA-256 digest")
    identity = payload.get("job_identity")
    if not isinstance(identity, str) or not re.fullmatch(r"[0-9a-f]{64}", identity):
        raise ValueError("job input job_identity must be a lowercase SHA-256 digest")
    timeout = payload.get("surface_timeout_s")
    if isinstance(timeout, bool) or not isinstance(timeout, (int, float)):
        raise ValueError("job input surface_timeout_s must be a positive finite number")
    if not math.isfinite(float(timeout)) or float(timeout) <= 0.0:
        raise ValueError("job input surface_timeout_s must be a positive finite number")
    threads = payload.get("omp_num_threads")
    if isinstance(threads, bool) or not isinstance(threads, int) or threads <= 0:
        raise ValueError("job input omp_num_threads must be a positive integer")
    surface = payload.get("surface")
    if not isinstance(surface, dict):
        raise ValueError("job input surface must be a JSON object")
    return payload


def _write_record(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _forward_lines(stream: TextIO, destination: TextIO, state: dict[str, object]) -> None:
    byte_count = 0
    for line in iter(stream.readline, ""):
        byte_count += len(line.encode("utf-8"))
        destination.write(line)
        destination.flush()
        if destination is sys.stdout:
            match = _PERIOD_RE.match(line.rstrip("\r\n"))
            if match is not None:
                state["closure_period"] = int(match.group("period"))
    state["bytes"] = byte_count
    stream.close()


def _terminate_process_group(process: subprocess.Popen[str]) -> int:
    """Terminate a timed-out solver and descendants, escalating after five seconds."""

    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        return int(process.wait(timeout=5.0))
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        return int(process.wait())


def run_bounded(
    command: list[str],
    *,
    timeout_s: float,
    omp_threads: int,
    record_path: Path,
    job_input_path: Path,
) -> int:
    """Run one solver process group and leave a provenance record on every outcome."""

    job = _read_job(job_input_path)
    if not command or command[0] != job["executable"]:
        raise ValueError("solver command does not match the staged job executable")
    if isinstance(timeout_s, bool) or not math.isfinite(float(timeout_s)) or timeout_s <= 0.0:
        raise ValueError("timeout_s must be positive and finite")
    if isinstance(omp_threads, bool) or not isinstance(omp_threads, int) or omp_threads <= 0:
        raise ValueError("omp_threads must be a positive integer")
    if omp_threads != job["omp_num_threads"]:
        raise ValueError("OMP thread count does not match the staged job input")
    if float(timeout_s) != float(job["surface_timeout_s"]):
        raise ValueError("solver timeout does not match the staged job input")

    executable = Path(command[0])
    expected_hash = str(job["executable_sha256"])
    record: dict[str, object] = {
        "schema_version": 1,
        "mode": "neo2-surface",
        "job_identity": job["job_identity"],
        "command": command,
        "executable": str(executable),
        "executable_sha256": expected_hash,
        "host": platform.node(),
        "omp_num_threads": str(omp_threads),
        "timeout_s": float(timeout_s),
        "surface": job["surface"],
        "start_time_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    started = time.monotonic()
    try:
        actual_hash = _sha256(executable)
        if actual_hash != expected_hash:
            raise OSError("solver executable SHA-256 differs from staged job input")
        environment = os.environ.copy()
        environment["OMP_NUM_THREADS"] = str(omp_threads)
        process = subprocess.Popen(
            command,
            cwd=Path.cwd(),
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,
        )
    except OSError as error:
        record.update(
            {
                "status": "launch_failure",
                "exit_code": EXIT_LAUNCH_FAILURE,
                "error": str(error),
                "wall_clock_s": time.monotonic() - started,
                "stdout_bytes": 0,
                "stderr_bytes": 0,
                "closure_period": None,
            }
        )
        _write_record(record_path, record)
        sys.stderr.write(f"could not start {executable}: {error}\n")
        return EXIT_LAUNCH_FAILURE

    stdout_state: dict[str, object] = {"bytes": 0, "closure_period": None}
    stderr_state: dict[str, object] = {"bytes": 0}
    readers = (
        threading.Thread(
            target=_forward_lines, args=(process.stdout, sys.stdout, stdout_state), daemon=True
        ),
        threading.Thread(
            target=_forward_lines, args=(process.stderr, sys.stderr, stderr_state), daemon=True
        ),
    )
    for reader in readers:
        reader.start()

    try:
        exit_code = int(process.wait(timeout=float(timeout_s)))
        status = "succeeded" if exit_code == 0 else "nonzero_exit"
    except subprocess.TimeoutExpired:
        _terminate_process_group(process)
        exit_code, status = EXIT_TIMEOUT, "timeout"
    for reader in readers:
        reader.join(timeout=5.0)

    record.update(
        {
            "status": status,
            "exit_code": exit_code,
            "wall_clock_s": time.monotonic() - started,
            "stdout_bytes": stdout_state["bytes"],
            "stderr_bytes": stderr_state["bytes"],
            "closure_period": stdout_state["closure_period"],
            "finish_time_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        }
    )
    _write_record(record_path, record)
    return exit_code


def main(argv: list[str] | None = None) -> int:
    """Run the staged surface named by ``condor_job.json`` in the current directory."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job-input", default="condor_job.json")
    parser.add_argument("--record", default="condor_run_record.json")
    arguments = parser.parse_args(argv)
    try:
        job_input = Path(arguments.job_input)
        job = _read_job(job_input)
        return run_bounded(
            [str(job["executable"])],
            timeout_s=float(job["surface_timeout_s"]),
            omp_threads=int(job["omp_num_threads"]),
            record_path=Path(arguments.record),
            job_input_path=job_input,
        )
    except ValueError as error:
        sys.stderr.write(f"invalid NEO-2 Condor job input: {error}\n")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
