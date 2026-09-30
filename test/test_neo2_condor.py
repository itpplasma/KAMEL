import subprocess
import sys
from pathlib import Path

import pytest

import neo2_for_Er.condor as condor
from neo2_for_Er.condor import (
    CondorError,
    CondorSubmitSpec,
    CondorToolConfig,
    parse_condor_submit_output,
    query_job_ads,
    query_job_ads_by_identity,
    run_condor,
    verify_shared_filesystem,
)


def _write_tool(path: Path, source: str) -> Path:
    path.write_text(f"#!{sys.executable}\n{source}", encoding="utf-8")
    path.chmod(0o755)
    return path


def test_submit_description_renders_resources_and_escaped_accounting(tmp_path: Path):
    spec = CondorSubmitSpec(
        executable=tmp_path / "neo2",
        initialdir=tmp_path,
        request_cpus=2,
        exclude_machines=("faepop43",),
        accounting={"MSKLabel": 'scan "red"'},
    )

    rendered = spec.render()

    assert "request_cpus = 2\n" in rendered
    assert 'requirements = (TARGET.Machine != "faepop43")\n' in rendered
    assert '+MSKLabel = "scan \\"red\\""\n' in rendered
    assert "queue 1\n" in rendered


def test_submit_description_does_not_forward_environment_by_default(tmp_path: Path):
    spec = CondorSubmitSpec(executable=tmp_path / "python", initialdir=tmp_path)

    assert "Getenv = false\n" in spec.render()


def test_submit_spec_refuses_dangling_symlink_destination(tmp_path: Path):
    outside = tmp_path / "outside.submit"
    destination = tmp_path / "job"
    destination.mkdir()
    (destination / "condor.submit").symlink_to(outside)
    spec = CondorSubmitSpec(executable=tmp_path / "python", initialdir=destination)

    with pytest.raises(CondorError, match="refusing to overwrite"):
        spec.write(destination)

    assert not outside.exists()


def test_submit_output_parses_cluster_id():
    assert parse_condor_submit_output("1 job(s) submitted to cluster 418.\n") == 418


def test_query_job_ads_normalizes_json_running_record(tmp_path: Path):
    tools = tmp_path / "bin"
    tools.mkdir()
    ad = {
        "ClusterId": 418,
        "ProcId": 0,
        "JobStatus": 2,
        "ExitCode": 0,
        "RemoteHost": "slot1@execute.example",
        "RemoteWallClockTime": 12.5,
        "NumHolds": 0,
        "HoldReason": "",
        "RequestCpus": 2,
        "RequestMemory": 8192,
        "ImageSize": 1024,
    }
    _write_tool(tools / "condor_q", f"import json\nprint(json.dumps([{ad!r}]))\n")
    _write_tool(tools / "condor_history", "print('[]')\n")

    records = query_job_ads([418], config=CondorToolConfig(tools))

    assert records == (
        {
            "cluster": 418,
            "proc": 0,
            "job_status_code": 2,
            "job_status": "running",
            "exit_code": 0,
            "remote_host": "slot1@execute.example",
            "remote_wall_clock_s": 12.5,
            "num_holds": 0,
            "hold_reason": "",
            "request_cpus": 2,
            "request_memory_mb": 8192,
            "job_identity": None,
            "ad_source": "condor_q",
        },
    )


def test_query_by_identity_uses_a_quoted_scheduler_constraint(tmp_path: Path, monkeypatch):
    identity = "a" * 64
    calls = []

    def query(arguments, *, config, tool):
        calls.append((tool, arguments))
        return subprocess.CompletedProcess(arguments, 0, "[]", "")

    monkeypatch.setattr(condor, "run_condor", query)

    assert query_job_ads_by_identity(identity, config=CondorToolConfig(tmp_path)) == ()

    expected = ["-constraint", f'KAMELJobIdentity == "{identity}"']
    assert len(calls) == 2
    assert all(arguments[-2:] == expected for _tool, arguments in calls)


def test_query_history_error_is_not_reported_as_missing_jobs(tmp_path: Path, monkeypatch):
    def query(arguments, *, config, tool):
        if tool == "condor_q":
            return subprocess.CompletedProcess(arguments, 0, "[]", "")
        return subprocess.CompletedProcess(arguments, 1, "", "history unavailable")

    monkeypatch.setattr(condor, "run_condor", query)

    with pytest.raises(CondorError, match="condor_history failed.*history unavailable"):
        query_job_ads([418], config=CondorToolConfig(tmp_path))


def test_run_condor_enforces_command_timeout(tmp_path: Path):
    tools = tmp_path / "bin"
    tools.mkdir()
    _write_tool(tools / "condor_status", "import time\ntime.sleep(2)\n")
    config = CondorToolConfig(tools, command_timeout_s=0.05)

    with pytest.raises(CondorError, match="condor_status exceeded"):
        run_condor([], config=config, tool="condor_status")


def test_transfer_never_requires_allowed_shared_prefix(tmp_path: Path):
    shared = tmp_path / "shared"
    run_root = shared / "runs" / "one"
    run_root.mkdir(parents=True)

    evidence = verify_shared_filesystem(
        run_root,
        allowed_prefixes=(shared,),
        should_transfer_files="NEVER",
    )

    assert evidence["status"] == "accepted"
    assert evidence["matched_prefix"] == str(shared)
    with pytest.raises(CondorError, match="outside the declared shared filesystem prefixes"):
        verify_shared_filesystem(
            tmp_path / "private",
            allowed_prefixes=(shared,),
            should_transfer_files="NEVER",
        )
    assert (
        verify_shared_filesystem(
            tmp_path / "private", allowed_prefixes=(), should_transfer_files="ALWAYS"
        )["status"]
        == "not_applicable"
    )


def test_condor_tool_config_rejects_nonpositive_timeout(tmp_path: Path):
    with pytest.raises(CondorError, match="strictly positive"):
        CondorToolConfig(tmp_path, command_timeout_s=0.0)
