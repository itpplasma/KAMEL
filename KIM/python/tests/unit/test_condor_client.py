from __future__ import annotations

import sys
from pathlib import Path

import pytest

import kim.condor_client as client
from kim.condor_client import (
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


def test_submit_spec_renders_eight_gib_and_escapes_values(tmp_path: Path) -> None:
    spec = CondorSubmitSpec(
        executable=tmp_path / "KIM solver.x",
        initialdir=tmp_path / "job directory",
        request_cpus=4,
        request_memory_mb=8192,
        accounting={"KIMLabel": 'scan "alpha"'},
    )

    rendered = spec.render()

    assert f'Executable = "{tmp_path / "KIM solver.x"}"\n' in rendered
    assert f'Initialdir = "{tmp_path / "job directory"}"\n' in rendered
    assert "request_cpus = 4\n" in rendered
    assert "request_memory = 8192\n" in rendered
    assert '+KIMLabel = "scan \\"alpha\\""\n' in rendered
    assert "queue 1\n" in rendered


def test_submit_spec_requires_absolute_executable_and_job_directory(tmp_path: Path) -> None:
    with pytest.raises(CondorError, match="absolute"):
        CondorSubmitSpec(executable=Path("KIM.x"), initialdir=tmp_path)

    with pytest.raises(CondorError, match="absolute"):
        CondorSubmitSpec(executable=tmp_path / "KIM.x", initialdir=Path("job"))


def test_submit_spec_rejects_newlines_in_requirements(tmp_path: Path) -> None:
    with pytest.raises(CondorError, match="requirements"):
        CondorSubmitSpec(
            executable=tmp_path / "KIM.x",
            initialdir=tmp_path,
            requirements=("TARGET.Memory > 1024\nqueue 99",),
        )


def test_submit_output_parses_cluster_id() -> None:
    assert parse_condor_submit_output("1 job(s) submitted to cluster 423.\n") == 423


def test_query_job_ads_normalizes_queue_and_history_json(tmp_path: Path) -> None:
    tools = tmp_path / "bin"
    tools.mkdir()
    queue_ad = {
        "ClusterId": 423,
        "ProcId": 0,
        "JobStatus": 2,
        "ExitCode": 0,
        "RemoteHost": "slot1@execute.example",
        "RemoteWallClockTime": 12.5,
        "NumHolds": 0,
        "HoldReason": "",
        "RequestCpus": 4,
        "RequestMemory": 8192,
        "KAMELJobIdentity": "a" * 64,
    }
    history_ad = {**queue_ad, "JobStatus": 4, "ExitCode": 0}
    _write_tool(tools / "condor_q", f"import json\nprint(json.dumps([{queue_ad!r}]))\n")
    _write_tool(tools / "condor_history", f"import json\nprint(json.dumps([{history_ad!r}]))\n")

    records = query_job_ads([423], config=CondorToolConfig(tools))

    assert [(record["job_status"], record["ad_source"]) for record in records] == [
        ("running", "condor_q"),
        ("completed", "condor_history"),
    ]
    assert all(record["job_identity"] == "a" * 64 for record in records)
    assert records[0]["request_memory_mb"] == 8192
    assert records[1]["exit_code"] == 0


def test_run_condor_enforces_command_timeout(tmp_path: Path) -> None:
    tools = tmp_path / "bin"
    tools.mkdir()
    _write_tool(tools / "condor_status", "import time\ntime.sleep(2)\n")

    with pytest.raises(CondorError, match="condor_status exceeded"):
        run_condor([], config=CondorToolConfig(tools, command_timeout_s=0.05), tool="condor_status")


def test_query_nonzero_command_result_raises_condor_error(tmp_path: Path) -> None:
    tools = tmp_path / "bin"
    tools.mkdir()
    _write_tool(tools / "condor_q", "raise SystemExit(7)\n")

    with pytest.raises(CondorError, match="condor_q failed with exit code 7"):
        query_job_ads([423], config=CondorToolConfig(tools), include_history=False)


def test_never_transfer_rejects_run_outside_shared_prefix(tmp_path: Path) -> None:
    shared = tmp_path / "shared"
    root = tmp_path / "private"
    shared.mkdir()
    root.mkdir()

    with pytest.raises(CondorError, match="outside the declared shared filesystem prefixes"):
        verify_shared_filesystem(
            root,
            allowed_prefixes=(shared,),
            should_transfer_files="NEVER",
        )


def test_condor_history_nonzero_result_is_not_treated_as_no_history(
    tmp_path: Path,
) -> None:
    tools = tmp_path / "bin"
    tools.mkdir()
    _write_tool(tools / "condor_q", "print('[]')\n")
    _write_tool(tools / "condor_history", "raise SystemExit(3)\n")

    with pytest.raises(CondorError, match="condor_history failed with exit code 3"):
        query_job_ads([423], config=CondorToolConfig(tools))


def test_condor_module_is_package_local() -> None:
    assert Path(client.__file__).resolve().is_relative_to(Path(__file__).parents[2].resolve())
