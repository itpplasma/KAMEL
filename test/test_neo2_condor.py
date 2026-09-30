import sys
from pathlib import Path

import pytest

from neo2_for_Er.condor import (
    CondorError,
    CondorSubmitSpec,
    CondorToolConfig,
    parse_condor_submit_output,
    query_job_ads,
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
