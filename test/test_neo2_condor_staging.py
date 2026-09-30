import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from neo2_for_Er import Neo2CondorError, Neo2CondorPlan, stage_neo2_condor_jobs
from neo2_for_Er.local_runner import stage_surfaces


def _prepare_surfaces(root: Path, *, count: int = 2) -> tuple[tuple[Path, ...], np.ndarray]:
    rows = np.array(
        [
            [0.25, 42.0, 150.0, 0.0, 1000.0, 1.0e13, -0.5],
            [0.50, 55.0, 160.0, 2.0, 1200.0, 1.2e13, -0.4],
        ][:count]
    )
    np.savetxt(root / "surfaces.dat", rows)
    template = root / "template"
    template.mkdir()
    (template / "neo2.in").write_text(
        "s=<boozer_s> k=<conl_over_mfp> R=<rbeg> Z=<zbeg>\n", encoding="utf-8"
    )
    return stage_surfaces(root, template), rows


def _plan(root: Path, *, runtime: dict[str, bool] | None = None, **overrides) -> Neo2CondorPlan:
    executable = root / "neo_2_par.x"
    executable.write_text(f"#!{sys.executable}\npass\n", encoding="utf-8")
    executable.chmod(0o755)
    values = {
        "executable": executable,
        "python_executable": Path(sys.executable),
        "shared_filesystem_prefixes": (root,),
        "parallel_runtime": runtime or {"omp": False, "mpi": False},
    }
    values.update(overrides)
    return Neo2CondorPlan(**values)


def test_plan_rejects_oversubscribed_openmp(tmp_path: Path):
    executable = tmp_path / "neo_2_par.x"
    executable.write_text("solver\n", encoding="utf-8")

    with pytest.raises(Neo2CondorError, match="links OpenMP"):
        Neo2CondorPlan(
            executable=executable,
            python_executable=Path(sys.executable),
            shared_filesystem_prefixes=(tmp_path,),
            parallel_runtime={"omp": True, "mpi": False},
            omp_threads_per_process=1,
            request_cpus=1,
        )


def test_stage_writes_payload_and_submit_description(tmp_path: Path):
    root = tmp_path / "work"
    root.mkdir()
    jobs, rows = _prepare_surfaces(root)
    original_inputs = {job: (job / "neo2.in").read_bytes() for job in jobs}
    plan = _plan(root, omp_threads_per_process=2, request_cpus=2, request_memory_mb=4096)

    staged = stage_neo2_condor_jobs(root, plan=plan)

    assert staged == jobs
    for job, row in zip(staged, rows, strict=True):
        job_input = json.loads((job / "condor_job.json").read_text(encoding="utf-8"))
        assert job_input["schema_version"] == 1
        assert job_input["mode"] == "neo2-surface"
        assert len(job_input["job_identity"]) == 64
        assert job_input["surface"]["boozer_s"] == pytest.approx(row[0])
        assert job_input["surface"]["r_eff_cm"] == pytest.approx(row[1])
        assert (
            job_input["executable_sha256"]
            == hashlib.sha256(plan.executable.read_bytes()).hexdigest()
        )
        assert (job / "condor_worker.py").read_bytes() == (
            Path(__file__).parents[1] / "python/neo2_for_Er/condor_worker.py"
        ).read_bytes()
        submit = (job / "condor.submit").read_text(encoding="utf-8")
        assert f"Executable = {plan.python_executable}\n" in submit
        assert 'Arguments = "condor_worker.py"\n' in submit
        assert "request_cpus = 2\n" in submit
        assert "request_memory = 4096\n" in submit
        assert "Getenv = false\n" in submit
        assert "should_transfer_files = NEVER\n" in submit
        assert f'+KAMELJobIdentity = "{job_input["job_identity"]}"\n' in submit
        assert (job / "neo2.in").read_bytes() == original_inputs[job]


@pytest.mark.parametrize("jobs_list", ["../escape\n", "s0d250000000e-01\n"])
def test_stage_rejects_unsafe_or_mismatched_jobs_list(tmp_path: Path, jobs_list: str):
    root = tmp_path / "work"
    root.mkdir()
    jobs, _rows = _prepare_surfaces(root)
    (root / "jobs_list.txt").write_text(jobs_list, encoding="utf-8")

    with pytest.raises(Neo2CondorError, match="jobs_list.txt"):
        stage_neo2_condor_jobs(root, plan=_plan(root))

    assert all(not (job / "condor_job.json").exists() for job in jobs)


@pytest.mark.parametrize("existing_file", ["condor_worker.py", "condor_job.json", "condor.submit"])
def test_stage_refuses_existing_payload(tmp_path: Path, existing_file: str):
    root = tmp_path / "work"
    root.mkdir()
    jobs, _rows = _prepare_surfaces(root)
    (jobs[0] / existing_file).write_text("preserve me\n", encoding="utf-8")

    with pytest.raises(Neo2CondorError, match="refusing to overwrite"):
        stage_neo2_condor_jobs(root, plan=_plan(root))

    assert (jobs[0] / existing_file).read_text(encoding="utf-8") == "preserve me\n"
    assert not (jobs[1] / "condor_job.json").exists()


def test_stage_rejects_root_outside_shared_filesystem_before_writing(tmp_path: Path):
    root = tmp_path / "private" / "work"
    root.mkdir(parents=True)
    jobs, _rows = _prepare_surfaces(root)
    plan = _plan(root, shared_filesystem_prefixes=(tmp_path / "shared",))

    with pytest.raises(Neo2CondorError, match="outside the declared shared filesystem prefixes"):
        stage_neo2_condor_jobs(root, plan=plan)

    assert all(not (job / "condor_job.json").exists() for job in jobs)
