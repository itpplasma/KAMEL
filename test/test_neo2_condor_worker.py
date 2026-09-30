import hashlib
import json
import sys
import time
from pathlib import Path

from neo2_for_Er.condor_worker import main


def _write_solver(path: Path, source: str) -> Path:
    path.write_text(f"#!{sys.executable}\n{source}", encoding="utf-8")
    path.chmod(0o755)
    return path


def _write_job(
    directory: Path,
    executable: Path,
    *,
    timeout_s: float = 3.0,
    omp_threads: int = 2,
) -> Path:
    path = directory / "condor_job.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "mode": "neo2-surface",
                "job_identity": "a" * 64,
                "executable": str(executable),
                "executable_sha256": (
                    hashlib.sha256(executable.read_bytes()).hexdigest()
                    if executable.is_file()
                    else "0" * 64
                ),
                "omp_num_threads": omp_threads,
                "surface_timeout_s": timeout_s,
                "surface": {
                    "boozer_s": 0.25,
                    "r_eff_cm": 42.0,
                    "r_beg_cm": 150.0,
                    "z_beg_cm": 0.0,
                    "ti_eV": 1000.0,
                    "ne_cm3": 1.0e13,
                    "kappa_cm_inv": -0.5,
                },
            }
        ),
        encoding="utf-8",
    )
    return path


def test_worker_records_success_and_requested_openmp_threads(tmp_path: Path, monkeypatch, capsys):
    executable = _write_solver(
        tmp_path / "neo2.x",
        "from pathlib import Path\n"
        "import os\n"
        "Path('omp.txt').write_text(os.environ['OMP_NUM_THREADS'])\n"
        "print('period: 12')\n",
    )
    _write_job(tmp_path, executable, omp_threads=3)
    monkeypatch.chdir(tmp_path)

    assert main([]) == 0

    record = json.loads((tmp_path / "condor_run_record.json").read_text(encoding="utf-8"))
    assert (tmp_path / "omp.txt").read_text(encoding="utf-8") == "3"
    assert "period: 12" in capsys.readouterr().out
    assert record["mode"] == "neo2-surface"
    assert record["job_identity"] == "a" * 64
    assert record["status"] == "succeeded"
    assert record["exit_code"] == 0
    assert record["host"]
    assert record["executable_sha256"] == hashlib.sha256(executable.read_bytes()).hexdigest()
    assert record["omp_num_threads"] == "3"
    assert record["closure_period"] == 12


def test_worker_preserves_nonzero_solver_exit(tmp_path: Path, monkeypatch):
    executable = _write_solver(tmp_path / "neo2.x", "raise SystemExit(9)\n")
    _write_job(tmp_path, executable)
    monkeypatch.chdir(tmp_path)

    assert main([]) == 9

    record = json.loads((tmp_path / "condor_run_record.json").read_text(encoding="utf-8"))
    assert record["status"] == "nonzero_exit"
    assert record["exit_code"] == 9


def test_worker_records_launch_failure(tmp_path: Path, monkeypatch):
    _write_job(tmp_path, tmp_path / "missing-neo2.x")
    monkeypatch.chdir(tmp_path)

    assert main([]) == 127

    record = json.loads((tmp_path / "condor_run_record.json").read_text(encoding="utf-8"))
    assert record["status"] == "launch_failure"
    assert record["exit_code"] == 127


def test_worker_timeout_terminates_solver_process_group(tmp_path: Path, monkeypatch):
    marker = tmp_path / "child-survived.txt"
    executable = _write_solver(
        tmp_path / "neo2.x",
        "import subprocess, sys, time\n"
        f"subprocess.Popen([sys.executable, '-c', \"import time; time.sleep(0.8); "
        f"open({str(marker)!r}, 'w').write('alive')\"])\n"
        "time.sleep(30)\n",
    )
    _write_job(tmp_path, executable, timeout_s=0.15)
    monkeypatch.chdir(tmp_path)

    assert main([]) == 124
    time.sleep(1.0)

    record = json.loads((tmp_path / "condor_run_record.json").read_text(encoding="utf-8"))
    assert record["status"] == "timeout"
    assert record["exit_code"] == 124
    assert not marker.exists()
