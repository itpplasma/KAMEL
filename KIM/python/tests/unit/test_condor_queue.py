from __future__ import annotations

import json
import runpy
import subprocess
import sys
from pathlib import Path

import pytest
from kim import BuiltinPlasma, PlasmaIsotope, ProfileScale, SimulationConfig, SweepSpec
from kim.condor import KimCondorPlan, stage_condor_sweep
from kim.condor_client import CondorError, CondorToolConfig

import kim.condor as condor

FIXTURES = Path(__file__).parents[1] / "fixtures"


def _staged(root: Path, *, factors: tuple[float, ...] = (-1.0, 0.0, 1.0), **plan_values):
    source = root.parent / "source-profiles"
    runpy.run_path(str(FIXTURES / "generate_profiles.py"))["generate_profiles"](source)
    config = SimulationConfig.electrostatic_periodic(
        profiles=source,
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
    executable = root.parent / "KIM.x"
    executable.write_text(f"#!{sys.executable}\npass\n", encoding="utf-8")
    executable.chmod(0o755)
    plan = KimCondorPlan(
        backend="kamel_kim_python",
        executable=executable,
        python_executable=Path(sys.executable),
        condor=CondorToolConfig(root.parent / "condor-bin"),
        shared_filesystem_prefixes=(root.parent,),
        **plan_values,
    )
    stage_condor_sweep(
        root, spec=SweepSpec(base=config, variation=ProfileScale(values=factors)), plan=plan
    )
    return plan


def _ad(
    cluster: int,
    status: int,
    identity: str,
    *,
    exit_code: int | None = 0,
    source: str = "condor_q",
) -> dict[str, object]:
    return {
        "cluster": cluster,
        "proc": 0,
        "job_status_code": status,
        "job_status": {1: "idle", 2: "running", 3: "removed", 4: "completed", 5: "held"}[status],
        "exit_code": exit_code,
        "remote_host": "slot1@execute.example",
        "remote_wall_clock_s": 1.5,
        "num_holds": 1 if status == 5 else 0,
        "hold_reason": "test hold" if status == 5 else None,
        "job_identity": identity,
        "ad_source": source,
    }


def _submit_fake(monkeypatch, clusters: list[int]):
    submitted = []

    def run(arguments, **_kwargs):
        submitted.append(arguments)
        cluster = clusters[len(submitted) - 1]
        return subprocess.CompletedProcess(
            arguments, 0, f"1 job submitted to cluster {cluster}\n", ""
        )

    monkeypatch.setattr(condor, "run_condor", run, raising=False)
    return submitted


def test_submit_tracks_each_scale_factor_cluster(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "run"
    plan = _staged(root)
    _submit_fake(monkeypatch, [401, 402, 403])

    result = condor.submit_condor_sweep(root, plan=plan, adopt_existing=False)

    assert [job["cluster"] for job in result["jobs"]] == [401, 402, 403]
    assert [job["profile_scale_factor"] for job in result["jobs"]] == [-1.0, 0.0, 1.0]
    manifest = json.loads((root / condor.MANIFEST_NAME).read_text(encoding="utf-8"))
    assert [job["state"] for job in manifest["jobs"]] == ["submitted"] * 3
    assert all((root / job["job"] / "condor_worker.py").is_file() for job in result["jobs"])


def test_submit_adopts_known_live_clusters(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "run"
    plan = _staged(root, factors=(1.0,))
    _submit_fake(monkeypatch, [417])
    first = condor.submit_condor_sweep(root, plan=plan, adopt_existing=False)
    identity = first["jobs"][0]["job_identity"]
    monkeypatch.setattr(
        condor,
        "query_job_ads",
        lambda clusters, **_kwargs: (_ad(417, 2, identity),),
        raising=False,
    )
    monkeypatch.setattr(
        condor, "run_condor", lambda *_args, **_kwargs: pytest.fail("duplicate submit")
    )

    restarted = condor.submit_condor_sweep(root, plan=plan)

    assert restarted["jobs"][0]["cluster"] == 417
    assert restarted["jobs"][0]["state"] == "adopted"
    assert restarted["adopted_clusters"] == [417]


def test_submit_fails_closed_after_ambiguous_ack(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "run"
    plan = _staged(root, factors=(1.0,))
    calls = []

    def lost_ack(arguments, **_kwargs):
        calls.append(arguments)
        return subprocess.CompletedProcess(arguments, 1, "", "acknowledgement lost")

    monkeypatch.setattr(condor, "run_condor", lost_ack, raising=False)
    with pytest.raises(CondorError, match="ambiguous"):
        condor.submit_condor_sweep(root, plan=plan, adopt_existing=False)

    with pytest.raises(CondorError, match="ambiguous"):
        condor.submit_condor_sweep(root, plan=plan)

    manifest = json.loads((root / condor.MANIFEST_NAME).read_text(encoding="utf-8"))
    assert manifest["jobs"][0]["state"] == "ambiguous"
    assert len(calls) == 1


def test_dry_run_does_not_submit(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "run"
    plan = _staged(root, factors=(1.0,))
    monkeypatch.setattr(
        condor, "run_condor", lambda *_args, **_kwargs: pytest.fail("submitted"), raising=False
    )

    result = condor.submit_condor_sweep(root, plan=plan, dry_run=True)

    assert result["dry_run"] is True
    assert result["jobs"][0]["state"] == "prepared"
    assert not (root / condor.MANIFEST_NAME).exists()


def test_wait_distinguishes_condor_terminal_states(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "run"
    plan = _staged(root)
    _submit_fake(monkeypatch, [421, 422, 423])
    submitted = condor.submit_condor_sweep(root, plan=plan, adopt_existing=False)
    identities = [job["job_identity"] for job in submitted["jobs"]]
    ads = (
        _ad(421, 4, identities[0], exit_code=0, source="condor_history"),
        _ad(422, 5, identities[1], exit_code=None),
        _ad(423, 3, identities[2], exit_code=None, source="condor_history"),
    )
    monkeypatch.setattr(condor, "query_job_ads", lambda *_args, **_kwargs: ads, raising=False)

    statuses = condor.wait_condor_sweep(root, plan=plan)

    assert [status["failure_kind"] for status in statuses] == [
        None,
        "condor_held",
        "condor_removed",
    ]
    assert all(status["terminal"] for status in statuses)
    persisted = json.loads((root / condor.STATUS_NAME).read_text(encoding="utf-8"))
    assert persisted["state"] == "failed"


def test_deadline_removes_pending_clusters_and_persists_status(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "run"
    plan = _staged(root, factors=(1.0,), max_wall_clock_s=0.04, poll_interval_s=0.01)
    _submit_fake(monkeypatch, [431])
    submitted = condor.submit_condor_sweep(root, plan=plan, adopt_existing=False)
    cluster = submitted["jobs"][0]["cluster"]
    identity = submitted["jobs"][0]["job_identity"]
    monkeypatch.setattr(
        condor,
        "query_job_ads",
        lambda *_args, **_kwargs: (_ad(cluster, 1, identity),),
        raising=False,
    )
    removed = []
    monkeypatch.setattr(
        condor,
        "remove_clusters",
        lambda clusters, **_kwargs: removed.extend(clusters),
        raising=False,
    )

    statuses = condor.wait_condor_sweep(root, plan=plan)

    assert removed == [cluster]
    assert statuses[0]["failure_kind"] == "condor_driver_wall_clock"
    assert statuses[0]["terminal"] is True
    persisted = json.loads((root / condor.STATUS_NAME).read_text(encoding="utf-8"))
    assert persisted["driver_removed_clusters"] == [cluster]
    reread_manifest = json.loads((root / condor.MANIFEST_NAME).read_text(encoding="utf-8"))
    assert reread_manifest["jobs"][0]["cluster"] == cluster
