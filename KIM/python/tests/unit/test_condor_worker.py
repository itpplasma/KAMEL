from __future__ import annotations

import json
import hashlib
import runpy
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import pytest
from kim import BuiltinPlasma, PlasmaIsotope, SimulationConfig

import kim.condor_worker as worker

FIXTURES = Path(__file__).parents[1] / "fixtures"


def _config(directory: Path) -> SimulationConfig:
    runpy.run_path(str(FIXTURES / "generate_profiles.py"))["generate_profiles"](directory)
    return SimulationConfig.electrostatic_periodic(
        profiles=directory,
        plasma=BuiltinPlasma(isotope=PlasmaIsotope.DEUTERIUM),
        btor=-17_977.413,
        major_radius=165.0,
        m_mode=7,
        n_mode=2,
        frequency=0.0,
        br_boundary_real=1.0,
        br_boundary_imag=0.0,
        radial_minimum=1.0,
        plasma_radius=9.0,
    )


def _executable(path: Path, source: str) -> Path:
    path.write_text(f"#!{sys.executable}\n{source}", encoding="utf-8")
    path.chmod(0o755)
    return path


def _write_job(
    directory: Path,
    *,
    backend: str,
    executable: Path,
    config: SimulationConfig | None = None,
    timeout_s: float = 3.0,
    factor: float = -0.25,
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 1,
        "mode": "kim-sweep" if backend == "kamel_kim_python" else "kim-run",
        "backend": backend,
        "job_identity": "b" * 64,
        "job": "scale-0000",
        "scan_order_index": 0,
        "profile_scale_factor": factor,
        "base_config": config.model_dump(mode="json") if config is not None else None,
        "executable": str(executable),
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "python_executable": sys.executable,
        "kim_source_path": None,
        "omp_num_threads": 3,
        "request_cpus": 4,
        "request_memory_mb": 8192,
        "job_timeout_s": timeout_s,
        "poll_interval_s": 0.01,
        "max_wall_clock_s": 5.0,
        "should_transfer_files": "NEVER",
        "api_current_component": "imag" if backend == "kamel_kim_python" else None,
        "api_current_unit": "A m" if backend == "kamel_kim_python" else None,
        "jpar_current_metric": None,
        "namelist_file": "KIM_config.nml" if backend == "kim_x_namelist" else None,
        "profile_directory": "profiles",
        "output_directory": ".",
        "input_sha256": {},
    }
    (directory / "condor_job.json").write_text(json.dumps(payload), encoding="utf-8")
    if backend == "kim_x_namelist":
        (directory / "KIM_config.nml").write_text("&kim_config /\n", encoding="utf-8")


def _record(directory: Path) -> dict[str, object]:
    return json.loads((directory / "condor_run_record.json").read_text(encoding="utf-8"))


def test_api_worker_runs_one_profile_scale_with_timeout(tmp_path: Path, monkeypatch) -> None:
    profiles = tmp_path / "profiles"
    config = _config(profiles)
    executable = _executable(tmp_path / "KIM.x", "pass\n")
    _write_job(
        tmp_path,
        backend="kamel_kim_python",
        executable=executable,
        config=config,
        timeout_s=47.0,
    )
    monkeypatch.chdir(tmp_path)
    calls = []

    class PeriodicResult:
        def integrated_parallel_current(self, *, region: str) -> complex:
            assert region == "as_is"
            return complex(1.25, -2.5)

    child = SimpleNamespace(periodic=PeriodicResult(), status="succeeded")

    def fake_run_sweep(spec, **kwargs):
        calls.append((spec, kwargs))
        return SimpleNamespace(children=(child,))

    monkeypatch.setattr("kim.sweep.run_sweep", fake_run_sweep)

    assert worker.main([]) == 0

    spec, kwargs = calls[0]
    assert spec.base.setup.m_mode == 7
    assert spec.variation.values == (-0.25,)
    assert kwargs["timeout"] == 47.0
    assert kwargs["environment"] == {"OMP_NUM_THREADS": "3"}
    record = _record(tmp_path)
    assert record["status"] == "succeeded"
    assert record["integrated_parallel_current_real"] == 1.25
    assert record["integrated_parallel_current_imag"] == -2.5
    assert record["api_current_value"] == -2.5
    assert record["api_current_unit"] == "A m"


def test_namelist_worker_runs_staged_config(tmp_path: Path, monkeypatch) -> None:
    executable = _executable(
        tmp_path / "KIM.x",
        "from pathlib import Path\n"
        "Path('solver-ran.txt').write_text(str(Path.cwd()))\n"
        "assert Path('KIM_config.nml').is_file()\n",
    )
    _write_job(tmp_path, backend="kim_x_namelist", executable=executable)
    monkeypatch.chdir(tmp_path)

    assert worker.main([]) == 0

    assert (tmp_path / "solver-ran.txt").read_text() == str(tmp_path)
    record = _record(tmp_path)
    assert record["status"] == "succeeded"
    assert record["backend"] == "kim_x_namelist"


@pytest.mark.parametrize(
    ("backend", "source", "timeout_s", "expected_status", "expected_exit"),
    [
        ("kim_x_namelist", "raise SystemExit(9)\n", 2.0, "solver_error", 9),
        (
            "kim_x_namelist",
            "import time\ntime.sleep(5)\n",
            0.05,
            "timeout",
            worker.EXIT_TIMEOUT,
        ),
        ("kamel_kim_python", None, 0.05, "timeout", worker.EXIT_TIMEOUT),
    ],
)
def test_worker_records_nonzero_exit_and_timeout(
    tmp_path: Path,
    monkeypatch,
    backend: str,
    source: str | None,
    timeout_s: float,
    expected_status: str,
    expected_exit: int,
) -> None:
    executable = _executable(tmp_path / "KIM.x", source or "pass\n")
    _write_job(
        tmp_path,
        backend=backend,
        executable=executable,
        config=_config(tmp_path / "profiles") if backend == "kamel_kim_python" else None,
        timeout_s=timeout_s,
    )
    monkeypatch.chdir(tmp_path)
    if backend == "kamel_kim_python":
        monkeypatch.setattr(
            "kim.sweep.run_sweep",
            lambda *_args, **_kwargs: SimpleNamespace(
                children=(SimpleNamespace(status="timed_out"),)
            ),
        )
    else:
        monkeypatch.setattr(worker, "_PROCESS_GROUP_TERM_GRACE_S", 0.05)

    assert worker.main([]) == expected_exit

    record = _record(tmp_path)
    assert record["status"] == expected_status
    assert record["exit_code"] == expected_exit


def test_worker_records_missing_api_dependency(tmp_path: Path, monkeypatch) -> None:
    executable = _executable(tmp_path / "KIM.x", "pass\n")
    _write_job(
        tmp_path,
        backend="kamel_kim_python",
        executable=executable,
        config=_config(tmp_path / "profiles"),
    )
    monkeypatch.chdir(tmp_path)

    def missing_api():
        raise ImportError("No module named kim")

    monkeypatch.setattr(worker, "_load_api_runner", missing_api)

    assert worker.main([]) == worker.EXIT_LAUNCH_FAILURE

    record = _record(tmp_path)
    assert record["status"] == "launch_failure"
    assert "No module named kim" in record["failure_reason"]


@pytest.mark.parametrize("job_input", [None, {"schema_version": 99, "mode": "kim-run"}])
def test_worker_rejects_missing_or_invalid_job_input(
    tmp_path: Path, monkeypatch, job_input: dict[str, object] | None
) -> None:
    if job_input is not None:
        (tmp_path / "condor_job.json").write_text(json.dumps(job_input), encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    assert worker.main([]) == worker.EXIT_INVALID_JOB

    record = _record(tmp_path)
    assert record["status"] == "invalid_job"
