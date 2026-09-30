import json
import os
import stat
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pytest

import neo2_for_Er.condor_runner as runner
from neo2_for_Er import Neo2CondorPlan, stage_neo2_condor_jobs
from neo2_for_Er.local_runner import stage_surfaces


def _staged(root: Path, *, count: int = 1, **plan_values):
    root.mkdir()
    rows = np.array(
        [
            [0.25 * (index + 1), 42.0 + index, 150.0, 0.0, 1000.0, 1.0e13, -0.5]
            for index in range(count)
        ]
    )
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
        **plan_values,
    )
    stage_neo2_condor_jobs(root, plan=plan)
    return jobs, plan


def _ad(
    cluster: int,
    status: int,
    exit_code: int | None = 0,
    *,
    job_identity: str | None = None,
) -> dict[str, object]:
    return {
        "cluster": cluster,
        "proc": 0,
        "job_status_code": status,
        "job_status": {
            1: "idle",
            2: "running",
            3: "removed",
            4: "completed",
            5: "held",
            6: "transferring_output",
            7: "suspended",
        }[status],
        "exit_code": exit_code,
        "remote_host": "slot1@execute.example",
        "remote_wall_clock_s": 2.5,
        "num_holds": 0,
        "hold_reason": None,
        "request_cpus": 1,
        "request_memory_mb": 1024,
        "job_identity": job_identity,
        "ad_source": "condor_q",
    }


def test_submit_records_clusters_and_adopts_live_jobs(tmp_path: Path, monkeypatch):
    root = tmp_path / "run"
    jobs, plan = _staged(root, count=2)
    clusters = iter((101, 102))
    submitted = []

    def submit(arguments, *, config, tool):
        assert tool == "condor_submit"
        submitted.append(arguments[0])
        cluster = next(clusters)
        return subprocess.CompletedProcess(
            arguments, 0, f"1 job submitted to cluster {cluster}\n", ""
        )

    monkeypatch.setattr(runner, "run_condor", submit)
    first = runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=False)

    assert [entry["cluster"] for entry in first["jobs"]] == [101, 102]
    assert all(entry["state"] == "submitted" for entry in first["jobs"])
    assert len(submitted) == len(jobs)
    assert json.loads((root / runner.MANIFEST_NAME).read_text())["submission_complete"] is True

    monkeypatch.setattr(
        runner, "run_condor", lambda *_args, **_kwargs: pytest.fail("duplicate submit")
    )
    identity = first["jobs"][0]["job_identity"]
    monkeypatch.setattr(
        runner,
        "query_job_ads",
        lambda clusters, **_kwargs: (_ad(101, 2, job_identity=identity),),
    )
    restarted = runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=True)

    assert restarted["adopted_clusters"] == [101]
    assert restarted["jobs"][0]["state"] == "adopted"
    assert restarted["jobs"][1]["cluster"] == 102


def test_concurrent_submitters_do_not_duplicate_condor_jobs(tmp_path: Path, monkeypatch):
    root = tmp_path / "run"
    _jobs, plan = _staged(root)
    original_load = runner._load_manifest

    def slow_load(*args, **kwargs):
        manifest = original_load(*args, **kwargs)
        time.sleep(0.1)
        return manifest

    monkeypatch.setattr(runner, "_load_manifest", slow_load)
    submitted = []
    submitted_lock = threading.Lock()

    def submit(arguments, **_kwargs):
        with submitted_lock:
            submitted.append(arguments[0])
            cluster = 451 + len(submitted)
        return subprocess.CompletedProcess(
            arguments, 0, f"1 job submitted to cluster {cluster}\n", ""
        )

    monkeypatch.setattr(runner, "run_condor", submit)
    start = threading.Barrier(3)

    def submitter():
        start.wait(timeout=5)
        return runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=False)

    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(submitter) for _ in range(2)]
        start.wait(timeout=5)
        results = [future.result(timeout=10) for future in futures]

    assert len(submitted) == 1
    assert all(result["jobs"][0]["cluster"] is not None for result in results)


def test_manifest_is_synced_to_disk_before_condor_submit(tmp_path: Path, monkeypatch):
    root = tmp_path / "run"
    _jobs, plan = _staged(root)
    synced_kinds = []
    original_fsync = os.fsync

    def track_fsync(descriptor):
        synced_kinds.append(stat.S_ISDIR(os.fstat(descriptor).st_mode))
        return original_fsync(descriptor)

    monkeypatch.setattr(runner.os, "fsync", track_fsync)

    def submit(arguments, **_kwargs):
        assert False in synced_kinds
        assert True in synced_kinds
        return subprocess.CompletedProcess(arguments, 0, "1 job submitted to cluster 99\n", "")

    monkeypatch.setattr(runner, "run_condor", submit)

    runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=False)

    assert False in synced_kinds
    assert True in synced_kinds


def test_submit_fails_closed_on_ambiguous_acknowledgement(tmp_path: Path, monkeypatch):
    root = tmp_path / "run"
    _jobs, plan = _staged(root)
    calls = []

    def submit(arguments, *, config, tool):
        calls.append(tool)
        return subprocess.CompletedProcess(arguments, 1, "", "acknowledgement lost")

    monkeypatch.setattr(runner, "run_condor", submit)

    with pytest.raises(runner.Neo2CondorError, match="ambiguous"):
        runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=False)

    manifest = json.loads((root / runner.MANIFEST_NAME).read_text(encoding="utf-8"))
    assert manifest["jobs"][0]["state"] == "ambiguous"
    with pytest.raises(runner.Neo2CondorError, match="ambiguous"):
        runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=False)
    assert calls == ["condor_submit"]


def test_submit_reconciles_stale_submitting_job_by_unique_identity(tmp_path: Path, monkeypatch):
    root = tmp_path / "run"
    _jobs, plan = _staged(root)
    monkeypatch.setattr(
        runner,
        "run_condor",
        lambda arguments, **_kwargs: subprocess.CompletedProcess(
            arguments, 1, "", "acknowledgement lost"
        ),
    )
    with pytest.raises(runner.Neo2CondorError, match="ambiguous"):
        runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=False)
    manifest_path = root / runner.MANIFEST_NAME
    stale_manifest = json.loads(manifest_path.read_text())
    stale_manifest["jobs"][0]["state"] = "submitting"
    manifest_path.write_text(json.dumps(stale_manifest), encoding="utf-8")
    identity = stale_manifest["jobs"][0]["job_identity"]
    monkeypatch.setattr(
        runner,
        "query_job_ads_by_identity",
        lambda requested_identity, **_kwargs: (_ad(441, 2, job_identity=requested_identity),),
    )
    monkeypatch.setattr(
        runner, "run_condor", lambda *_args, **_kwargs: pytest.fail("duplicate submit")
    )
    monkeypatch.setattr(
        runner,
        "query_job_ads",
        lambda clusters, **_kwargs: (_ad(441, 2, job_identity=identity),),
    )

    result = runner.submit_neo2_condor_jobs(root, plan=plan)

    assert result["jobs"][0]["cluster"] == 441
    assert result["jobs"][0]["state"] == "adopted"
    assert result["adopted_clusters"] == [441]


def test_submit_keeps_unmatched_ambiguous_job_from_being_resubmitted(tmp_path: Path, monkeypatch):
    root = tmp_path / "run"
    _jobs, plan = _staged(root)
    submits = []

    def failed_submit(arguments, **_kwargs):
        submits.append(arguments[0])
        return subprocess.CompletedProcess(arguments, 1, "", "acknowledgement lost")

    monkeypatch.setattr(runner, "run_condor", failed_submit)
    with pytest.raises(runner.Neo2CondorError, match="ambiguous"):
        runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=False)
    monkeypatch.setattr(runner, "query_job_ads_by_identity", lambda *_args, **_kwargs: ())

    with pytest.raises(runner.Neo2CondorError, match="ambiguous"):
        runner.submit_neo2_condor_jobs(root, plan=plan)

    manifest = json.loads((root / runner.MANIFEST_NAME).read_text(encoding="utf-8"))
    assert manifest["jobs"][0]["state"] == "ambiguous"
    assert "expected exactly one" in manifest["jobs"][0]["error"]
    assert len(submits) == 1


def test_submit_rejects_multiple_process_matches_for_one_identity(tmp_path: Path, monkeypatch):
    root = tmp_path / "run"
    _jobs, plan = _staged(root)
    monkeypatch.setattr(
        runner,
        "run_condor",
        lambda arguments, **_kwargs: subprocess.CompletedProcess(
            arguments, 1, "", "acknowledgement lost"
        ),
    )
    with pytest.raises(runner.Neo2CondorError, match="ambiguous"):
        runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=False)
    identity = json.loads((root / runner.MANIFEST_NAME).read_text())["jobs"][0]["job_identity"]
    proc_zero = _ad(441, 2, job_identity=identity)
    proc_one = {**_ad(441, 2, job_identity=identity), "proc": 1}
    monkeypatch.setattr(
        runner,
        "query_job_ads_by_identity",
        lambda *_args, **_kwargs: (proc_zero, proc_one),
    )

    with pytest.raises(runner.Neo2CondorError, match="ambiguous"):
        runner.submit_neo2_condor_jobs(root, plan=plan)

    manifest = json.loads((root / runner.MANIFEST_NAME).read_text(encoding="utf-8"))
    assert manifest["jobs"][0]["state"] == "ambiguous"
    assert "2 unique scheduler jobs" in manifest["jobs"][0]["error"]


def test_submit_refuses_to_adopt_cluster_with_a_different_identity(tmp_path: Path, monkeypatch):
    root = tmp_path / "run"
    _jobs, plan = _staged(root, count=2)
    clusters = iter((111, 112))
    submitted = []

    def submit(arguments, *, config, tool):
        submitted.append(tool)
        return subprocess.CompletedProcess(
            arguments, 0, f"1 job submitted to cluster {next(clusters)}\n", ""
        )

    monkeypatch.setattr(runner, "run_condor", submit)
    initial = runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=False)
    wrong_identity = "f" * 64
    monkeypatch.setattr(
        runner,
        "query_job_ads",
        lambda *_args, **_kwargs: (_ad(111, 2, job_identity=wrong_identity),),
    )

    with pytest.raises(runner.Neo2CondorError, match="identity does not match"):
        runner.submit_neo2_condor_jobs(root, plan=plan)

    manifest = json.loads((root / runner.MANIFEST_NAME).read_text(encoding="utf-8"))
    assert manifest["jobs"][0]["state"] == "ambiguous"
    assert [job["cluster"] for job in initial["jobs"]] == [111, 112]
    assert submitted == ["condor_submit", "condor_submit"]


def test_submit_fails_closed_after_interrupted_submission(tmp_path: Path, monkeypatch):
    root = tmp_path / "run"
    _jobs, plan = _staged(root)
    calls = []

    def interrupted_submit(arguments, *, config, tool):
        calls.append(tool)
        raise RuntimeError("simulated client interruption")

    monkeypatch.setattr(runner, "run_condor", interrupted_submit)
    with pytest.raises(RuntimeError, match="simulated client interruption"):
        runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=False)

    manifest = json.loads((root / runner.MANIFEST_NAME).read_text(encoding="utf-8"))
    assert manifest["jobs"][0]["state"] == "submitting"
    with pytest.raises(runner.Neo2CondorError, match="ambiguous"):
        runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=False)
    assert manifest["jobs"][0]["cluster"] is None
    assert calls == ["condor_submit"]


def test_submit_rejects_payload_modified_after_staging(tmp_path: Path, monkeypatch):
    root = tmp_path / "run"
    jobs, plan = _staged(root)
    job_input = jobs[0] / "condor_job.json"
    payload = json.loads(job_input.read_text(encoding="utf-8"))
    payload["surface"]["ne_cm3"] *= 2.0
    job_input.write_text(json.dumps(payload), encoding="utf-8")
    monkeypatch.setattr(
        runner, "run_condor", lambda *_args, **_kwargs: pytest.fail("submitted tampered payload")
    )

    with pytest.raises(runner.Neo2CondorError, match="surface metadata does not match"):
        runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=False)

    assert not (root / runner.MANIFEST_NAME).exists()


def test_dry_run_does_not_call_condor_submit(tmp_path: Path, monkeypatch):
    root = tmp_path / "run"
    _jobs, plan = _staged(root)
    monkeypatch.setattr(runner, "run_condor", lambda *_args, **_kwargs: pytest.fail("submitted"))

    result = runner.submit_neo2_condor_jobs(root, plan=plan, dry_run=True, adopt_existing=False)

    assert result["dry_run"] is True
    assert result["jobs"][0]["cluster"] is None
    assert result["jobs"][0]["state"] == "dry_run"
    assert not (root / runner.MANIFEST_NAME).exists()


def test_wait_classifies_terminal_job_states(tmp_path: Path, monkeypatch):
    root = tmp_path / "run"
    _jobs, plan = _staged(root, count=4)
    clusters = iter((201, 202, 203, 204))

    def submit(arguments, *, config, tool):
        return subprocess.CompletedProcess(
            arguments, 0, f"1 job submitted to cluster {next(clusters)}\n", ""
        )

    monkeypatch.setattr(runner, "run_condor", submit)
    runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=False)
    manifest = json.loads((root / runner.MANIFEST_NAME).read_text(encoding="utf-8"))
    identities = [job["job_identity"] for job in manifest["jobs"]]
    monkeypatch.setattr(
        runner,
        "query_job_ads",
        lambda *_args, **_kwargs: (
            _ad(201, 4, 0, job_identity=identities[0]),
            _ad(202, 5, 1, job_identity=identities[1]),
            _ad(203, 3, 1, job_identity=identities[2]),
            _ad(204, 4, 7, job_identity=identities[3]),
        ),
    )

    statuses = runner.wait_neo2_condor_jobs(root, plan=plan)

    assert [status["job_status"] for status in statuses] == [
        "completed",
        "held",
        "removed",
        "completed",
    ]
    assert [status["failure_kind"] for status in statuses] == [
        None,
        "condor_held",
        "condor_removed",
        "condor_nonzero_exit",
    ]


def test_wait_records_missing_cluster_as_distinct_failure(tmp_path: Path, monkeypatch):
    root = tmp_path / "run"
    _jobs, plan = _staged(root)
    monkeypatch.setattr(
        runner,
        "run_condor",
        lambda arguments, **_kwargs: subprocess.CompletedProcess(
            arguments, 0, "1 job submitted to cluster 251\n", ""
        ),
    )
    runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=False)
    monkeypatch.setattr(runner, "query_job_ads", lambda *_args, **_kwargs: ())

    statuses = runner.wait_neo2_condor_jobs(root, plan=plan)

    assert statuses[0]["job_status"] == "missing"
    assert statuses[0]["failure_kind"] == "condor_missing"
    assert statuses[0]["terminal"] is True


def test_wait_retries_query_errors_then_removes_pending_jobs_at_deadline(
    tmp_path: Path, monkeypatch
):
    root = tmp_path / "run"
    _jobs, plan = _staged(root, max_wall_clock_s=0.04, poll_interval_s=0.01)
    monkeypatch.setattr(
        runner,
        "run_condor",
        lambda arguments, **_kwargs: subprocess.CompletedProcess(
            arguments, 0, "1 job submitted to cluster 261\n", ""
        ),
    )
    runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=False)

    query_errors = []

    def query_failure(*_args, **_kwargs):
        query_errors.append("condor_q unavailable")
        raise runner.CondorError("condor_q unavailable")

    monkeypatch.setattr(runner, "query_job_ads", query_failure)
    removed = []
    monkeypatch.setattr(
        runner, "remove_clusters", lambda clusters, **_kwargs: removed.extend(clusters)
    )

    statuses = runner.wait_neo2_condor_jobs(root, plan=plan)

    assert len(query_errors) >= 2
    assert removed == [261]
    assert statuses[0]["failure_kind"] == "condor_driver_wall_clock"
    assert statuses[0]["query_error"] == "condor_q unavailable"
    status_document = json.loads((root / runner.STATUS_NAME).read_text(encoding="utf-8"))
    assert status_document["driver_removed_clusters"] == [261]
    assert status_document["query_error"] == "condor_q unavailable"


def test_wait_deadline_removes_pending_jobs_and_persists_status(tmp_path: Path, monkeypatch):
    root = tmp_path / "run"
    _jobs, plan = _staged(root, max_wall_clock_s=0.01, poll_interval_s=0.05)
    monkeypatch.setattr(
        runner,
        "run_condor",
        lambda arguments, **_kwargs: subprocess.CompletedProcess(
            arguments, 0, "1 job submitted to cluster 301\n", ""
        ),
    )
    runner.submit_neo2_condor_jobs(root, plan=plan, adopt_existing=False)
    manifest = json.loads((root / runner.MANIFEST_NAME).read_text(encoding="utf-8"))
    identity = manifest["jobs"][0]["job_identity"]
    monkeypatch.setattr(
        runner,
        "query_job_ads",
        lambda *_args, **_kwargs: (_ad(301, 2, job_identity=identity),),
    )
    removed = []
    monkeypatch.setattr(
        runner,
        "remove_clusters",
        lambda clusters, **_kwargs: removed.extend(clusters),
    )

    statuses = runner.wait_neo2_condor_jobs(root, plan=plan)

    assert removed == [301]
    assert statuses[0]["job_status"] == "removed"
    assert statuses[0]["failure_kind"] == "condor_driver_wall_clock"
    status_document = json.loads((root / runner.STATUS_NAME).read_text(encoding="utf-8"))
    assert status_document["driver_removed_clusters"] == [301]
    assert status_document["statuses"][0]["failure_kind"] == "condor_driver_wall_clock"
