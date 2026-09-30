import hashlib
import json
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

import neo2_for_Er.condor_runner as runner
from neo2_for_Er import Neo2CondorPlan, collect_neo2_condor_results, stage_neo2_condor_jobs
from neo2_for_Er.local_runner import stage_surfaces


def _prepare_completed_run(
    root: Path,
    monkeypatch,
    rows: list[list[float]],
    *,
    max_closure_periods: int = 10,
) -> tuple[tuple[Path, ...], Neo2CondorPlan]:
    root.mkdir()
    np.savetxt(root / "surfaces.dat", rows)
    template = root / "template"
    template.mkdir()
    (template / "neo2.in").write_text(
        "s=<boozer_s> k=<conl_over_mfp> R=<rbeg> Z=<zbeg>\n", encoding="utf-8"
    )
    jobs = stage_surfaces(root, template)
    executable = root / "neo_2_par.x"
    executable.write_text(f"#!{sys.executable}\npass\n", encoding="utf-8")
    executable.chmod(0o755)
    plan = Neo2CondorPlan(
        executable=executable,
        python_executable=Path(sys.executable),
        shared_filesystem_prefixes=(root,),
        parallel_runtime={"omp": False, "mpi": False},
        max_closure_periods=max_closure_periods,
    )
    stage_neo2_condor_jobs(root, plan=plan)

    next_cluster = iter(range(501, 501 + len(jobs)))

    def submit(arguments, *, config, tool):
        cluster = next(next_cluster)
        return subprocess.CompletedProcess(
            arguments, 0, f"1 job submitted to cluster {cluster}\n", ""
        )

    monkeypatch.setattr(runner, "run_condor", submit)
    runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=False)
    manifest = json.loads((root / runner.MANIFEST_NAME).read_text(encoding="utf-8"))
    identities = [job["job_identity"] for job in manifest["jobs"]]
    clusters = [job["cluster"] for job in manifest["jobs"]]
    ads = tuple(
        {
            "cluster": cluster,
            "proc": 0,
            "job_status_code": 4,
            "job_status": "completed",
            "exit_code": 0,
            "remote_host": "slot1@execute.example",
            "remote_wall_clock_s": 1.25,
            "num_holds": 0,
            "hold_reason": None,
            "request_cpus": 1,
            "request_memory_mb": 30_720,
            "job_identity": identity,
            "ad_source": "condor_history",
        }
        for cluster, identity in zip(clusters, identities)
    )
    monkeypatch.setattr(runner, "query_job_ads", lambda *_args, **_kwargs: ads)
    runner.wait_neo2_condor_jobs(root, plan=plan)

    executable_hash = hashlib.sha256(executable.read_bytes()).hexdigest()
    for job in jobs:
        job_input = json.loads((job / "condor_job.json").read_text(encoding="utf-8"))
        worker_record = {
            "schema_version": 1,
            "mode": "neo2-surface",
            "job_identity": job_input["job_identity"],
            "status": "succeeded",
            "exit_code": 0,
            "executable": str(executable),
            "executable_sha256": executable_hash,
            "omp_num_threads": str(plan.omp_threads_per_process),
            "host": "execute.example",
            "surface": job_input["surface"],
            "closure_period": 1,
        }
        (job / "condor_run_record.json").write_text(json.dumps(worker_record), encoding="utf-8")
        with h5py.File(job / "neo2_config.h5", "w") as config:
            config.create_dataset("settings/boozer_s", data=job_input["surface"]["boozer_s"])
        with h5py.File(job / "fulltransp.h5", "w") as transport:
            transport.create_dataset("k_cof", data=job_input["surface"]["boozer_s"] / 2.0)
    return jobs, plan


def _record_for(job: Path, results):
    return next(record for record in results.surface_records if record["job_directory"] == str(job))


def test_collects_valid_points_in_radius_order(tmp_path: Path, monkeypatch):
    rows = [
        [0.5, 99.0, 150.0, 0.0, 1000.0, 1.0e13, -0.5],
        [0.25, 11.0, 160.0, 2.0, 1200.0, 1.2e13, -0.4],
    ]
    root = tmp_path / "run"
    _jobs, plan = _prepare_completed_run(root, monkeypatch, rows)

    results = collect_neo2_condor_results(root, plan=plan)

    np.testing.assert_allclose(results.profile, [[11.0, 0.125], [99.0, 0.25]])
    assert [record["status"] for record in results.surface_records] == ["succeeded", "succeeded"]


def test_collection_retains_partial_surface_failure(tmp_path: Path, monkeypatch):
    rows = [
        [0.25, 42.0, 150.0, 0.0, 1000.0, 1.0e13, -0.5],
        [0.5, 55.0, 160.0, 2.0, 1200.0, 1.2e13, -0.4],
    ]
    root = tmp_path / "run"
    jobs, plan = _prepare_completed_run(root, monkeypatch, rows)
    worker_path = jobs[1] / "condor_run_record.json"
    worker_record = json.loads(worker_path.read_text(encoding="utf-8"))
    worker_record.update(status="nonzero_exit", exit_code=9)
    worker_path.write_text(json.dumps(worker_record), encoding="utf-8")

    results = collect_neo2_condor_results(root, plan=plan)

    np.testing.assert_allclose(results.profile, [[42.0, 0.125]])
    failed = _record_for(jobs[1], results)
    assert failed["status"] == "failed"
    assert failed["failure_kind"] == "worker_nonzero_exit"


def test_collection_attributes_scheduler_nonzero_exit_to_worker_record(tmp_path: Path, monkeypatch):
    rows = [[0.25, 42.0, 150.0, 0.0, 1000.0, 1.0e13, -0.5]]
    root = tmp_path / "run"
    jobs, plan = _prepare_completed_run(root, monkeypatch, rows)
    worker_path = jobs[0] / "condor_run_record.json"
    worker_record = json.loads(worker_path.read_text(encoding="utf-8"))
    worker_record.update(status="nonzero_exit", exit_code=9)
    worker_path.write_text(json.dumps(worker_record), encoding="utf-8")
    status_path = root / runner.STATUS_NAME
    status_document = json.loads(status_path.read_text(encoding="utf-8"))
    status_document["statuses"][0].update(
        exit_code=9,
        failure_kind="condor_nonzero_exit",
        failure_reason="job completed with exit code 9",
    )
    status_path.write_text(json.dumps(status_document), encoding="utf-8")

    results = collect_neo2_condor_results(root, plan=plan)

    record = _record_for(jobs[0], results)
    assert record["failure_kind"] == "worker_nonzero_exit"
    assert record["failure_reason"] == "NEO-2 solver exited with code 9"


@pytest.mark.parametrize(
    ("field_name", "wrong_value"),
    [
        ("job_identity", "f" * 64),
        ("executable_sha256", "0" * 64),
        ("omp_num_threads", "999"),
    ],
)
def test_collection_rejects_mismatched_worker_provenance(
    tmp_path: Path, monkeypatch, field_name: str, wrong_value: str
):
    rows = [[0.25, 42.0, 150.0, 0.0, 1000.0, 1.0e13, -0.5]]
    root = tmp_path / "run"
    jobs, plan = _prepare_completed_run(root, monkeypatch, rows)
    worker_path = jobs[0] / "condor_run_record.json"
    worker_record = json.loads(worker_path.read_text(encoding="utf-8"))
    worker_record[field_name] = wrong_value
    worker_path.write_text(json.dumps(worker_record), encoding="utf-8")

    results = collect_neo2_condor_results(root, plan=plan)

    assert results.profile is None
    record = _record_for(jobs[0], results)
    assert record["failure_kind"] == "worker_provenance_mismatch"
    assert field_name in record["failure_reason"]


@pytest.mark.parametrize(
    "invalid_output",
    [
        "missing_k_cof",
        "nonfinite_k_cof",
        "mismatched_boozer_s",
        "nonfinite_boozer_s",
    ],
)
def test_collection_rejects_invalid_hdf5_output(tmp_path: Path, monkeypatch, invalid_output: str):
    rows = [[0.25, 42.0, 150.0, 0.0, 1000.0, 1.0e13, -0.5]]
    root = tmp_path / "run"
    jobs, plan = _prepare_completed_run(root, monkeypatch, rows)
    if invalid_output in {"mismatched_boozer_s", "nonfinite_boozer_s"}:
        with h5py.File(jobs[0] / "neo2_config.h5", "w") as config:
            boozer_s = 0.5 if invalid_output == "mismatched_boozer_s" else np.nan
            config.create_dataset("settings/boozer_s", data=boozer_s)
    else:
        with h5py.File(jobs[0] / "fulltransp.h5", "w") as transport:
            if invalid_output == "nonfinite_k_cof":
                transport.create_dataset("k_cof", data=np.nan)

    results = collect_neo2_condor_results(root, plan=plan)

    assert results.profile is None
    record = _record_for(jobs[0], results)
    assert record["failure_kind"] == "invalid_hdf5_output"
    assert record["failure_reason"]


def test_collection_marks_excess_closure_period(tmp_path: Path, monkeypatch):
    rows = [[0.25, 42.0, 150.0, 0.0, 1000.0, 1.0e13, -0.5]]
    root = tmp_path / "run"
    jobs, plan = _prepare_completed_run(root, monkeypatch, rows, max_closure_periods=2)
    worker_path = jobs[0] / "condor_run_record.json"
    worker_record = json.loads(worker_path.read_text(encoding="utf-8"))
    worker_record["closure_period"] = 3
    worker_path.write_text(json.dumps(worker_record), encoding="utf-8")

    results = collect_neo2_condor_results(root, plan=plan)

    assert results.profile is None
    record = _record_for(jobs[0], results)
    assert record["failure_kind"] == "closure_period_exceeded"
    assert record["closure_period"] == 3


def test_collection_does_not_treat_pending_condor_job_as_completed(tmp_path: Path, monkeypatch):
    rows = [[0.25, 42.0, 150.0, 0.0, 1000.0, 1.0e13, -0.5]]
    root = tmp_path / "run"
    jobs, plan = _prepare_completed_run(root, monkeypatch, rows)
    status_path = root / runner.STATUS_NAME
    status_document = json.loads(status_path.read_text(encoding="utf-8"))
    status_document["statuses"][0].update(
        job_status="running",
        job_status_code=2,
        exit_code=None,
        failure_kind=None,
        failure_reason=None,
        terminal=False,
    )
    status_path.write_text(json.dumps(status_document), encoding="utf-8")

    results = collect_neo2_condor_results(root, plan=plan)

    assert results.profile is None
    record = _record_for(jobs[0], results)
    assert record["status"] == "pending"
    assert record["failure_kind"] == "condor_pending"


def test_collection_rejects_boolean_exit_code(tmp_path: Path, monkeypatch):
    rows = [[0.25, 42.0, 150.0, 0.0, 1000.0, 1.0e13, -0.5]]
    root = tmp_path / "run"
    jobs, plan = _prepare_completed_run(root, monkeypatch, rows)
    worker_path = jobs[0] / "condor_run_record.json"
    worker_record = json.loads(worker_path.read_text(encoding="utf-8"))
    worker_record["exit_code"] = False
    worker_path.write_text(json.dumps(worker_record), encoding="utf-8")

    results = collect_neo2_condor_results(root, plan=plan)

    assert results.profile is None
    assert _record_for(jobs[0], results)["failure_kind"] == "worker_failed"


def test_collection_rejects_status_for_wrong_cluster(tmp_path: Path, monkeypatch):
    rows = [[0.25, 42.0, 150.0, 0.0, 1000.0, 1.0e13, -0.5]]
    root = tmp_path / "run"
    jobs, plan = _prepare_completed_run(root, monkeypatch, rows)
    status_path = root / runner.STATUS_NAME
    status_document = json.loads(status_path.read_text(encoding="utf-8"))
    status_document["statuses"][0]["cluster"] += 100
    status_path.write_text(json.dumps(status_document), encoding="utf-8")

    results = collect_neo2_condor_results(root, plan=plan)

    assert results.profile is None
    assert _record_for(jobs[0], results)["failure_kind"] == "condor_status_identity_mismatch"
