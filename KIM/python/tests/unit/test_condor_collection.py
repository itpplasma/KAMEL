from __future__ import annotations

import json
import runpy
import subprocess
import sys
from pathlib import Path

import f90nml
import numpy as np
import pytest
from kim import BuiltinPlasma, PlasmaIsotope, ProfileScale, SimulationConfig, SweepSpec
from kim.condor import JparCurrentMetric, KimCondorPlan, stage_condor_sweep
from kim.condor_client import CondorError, CondorToolConfig

import kim.condor as condor

FIXTURES = Path(__file__).parents[1] / "fixtures"


def _stage(
    root: Path,
    monkeypatch,
    *,
    backend: str = "kamel_kim_python",
    current_component: str | None = None,
    metric: JparCurrentMetric | None = None,
    factors: tuple[float, ...] = (-1.0, 0.0, 1.0),
):
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
    reviewed = None
    if backend == "kim_x_namelist":
        from kim.namelist import dumps_namelist

        reviewed = root.parent / "reviewed.nml"
        reviewed.write_text(dumps_namelist(config), encoding="utf-8")
        namelist = f90nml.read(reviewed)
        namelist["kim_config"]["type_of_run"] = "electrostatic"
        namelist["kim_config"]["collision_model"] = (
            metric.collision_model if metric is not None else "FokkerPlanck"
        )
        namelist["kim_io"]["hdf5_output"] = False
        f90nml.write(namelist, reviewed, force=True)
    plan = KimCondorPlan(
        backend=backend,
        executable=executable,
        python_executable=Path(sys.executable),
        condor=CondorToolConfig(root.parent / "condor-bin"),
        shared_filesystem_prefixes=(root.parent,),
        base_namelist=reviewed,
        jpar_current_metric=metric,
        api_current_component=current_component,
        api_current_unit="A m" if current_component is not None else None,
    )
    jobs = stage_condor_sweep(
        root,
        spec=SweepSpec(base=config, variation=ProfileScale(values=factors)),
        plan=plan,
    )
    cluster = iter(range(501, 501 + len(jobs)))
    monkeypatch.setattr(
        condor,
        "run_condor",
        lambda arguments, **_kwargs: subprocess.CompletedProcess(
            arguments, 0, f"1 job submitted to cluster {next(cluster)}\n", ""
        ),
    )
    submission = condor.submit_condor_sweep(root, plan=plan, adopt_existing=False)
    (root / condor.STATUS_NAME).write_text(
        json.dumps(
            {
                "schema_version": 1,
                "state": "completed",
                "statuses": [
                    {
                        "job": entry["job"],
                        "job_identity": entry["job_identity"],
                        "scan_order_index": entry["scan_order_index"],
                        "profile_scale_factor": entry["profile_scale_factor"],
                        "cluster": entry["cluster"],
                        "job_status": "completed",
                        "job_status_code": 4,
                        "exit_code": 0,
                        "failure_kind": None,
                        "failure_reason": None,
                        "terminal": True,
                        "driver_removed": False,
                    }
                    for entry in submission["jobs"]
                ],
            }
        ),
        encoding="utf-8",
    )
    return plan, jobs, submission


def _write_records(jobs, *, currents: tuple[complex, ...] = (3 + 8j, 1 + 2j, 2 + 6j)) -> None:
    for job, current in zip(jobs, currents):
        payload = json.loads((job.directory / "condor_job.json").read_text(encoding="utf-8"))
        component = payload["api_current_component"]
        scalar = (
            None if component is None else current.imag if component == "imag" else current.real
        )
        record = {
            "schema_version": 1,
            "job_identity": payload["job_identity"],
            "backend": payload["backend"],
            "scan_order_index": payload["scan_order_index"],
            "profile_scale_factor": payload["profile_scale_factor"],
            "executable_sha256": payload["executable_sha256"],
            "status": "succeeded",
            "exit_code": 0,
            "host": "execute.example",
            "integrated_parallel_current_real": current.real,
            "integrated_parallel_current_imag": current.imag,
            "api_current_component": component,
            "api_current_unit": payload["api_current_unit"],
            "api_current_value": scalar,
            "started_at_utc": "2026-09-30T00:00:00+00:00",
            "finished_at_utc": "2026-09-30T00:00:01+00:00",
        }
        (job.directory / "condor_run_record.json").write_text(json.dumps(record), encoding="utf-8")


def test_collection_preserves_complex_api_current(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "run"
    plan, jobs, _submission = _stage(root, monkeypatch)
    currents = (3 + 8j, 1 + 2j, 2 + 6j)
    _write_records(jobs, currents=currents)

    results = condor.collect_condor_sweep(root, plan=plan)

    assert [
        (point["integrated_parallel_current_real"], point["integrated_parallel_current_imag"])
        for point in results.curve
    ] == [(value.real, value.imag) for value in currents]


def test_collection_withholds_minimum_without_metric(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "run"
    plan, jobs, _submission = _stage(root, monkeypatch)
    _write_records(jobs)

    results = condor.collect_condor_sweep(root, plan=plan)

    assert results.resonance["status"] == "insufficient_scalars"
    assert results.resonance["profile_scale_factor"] is None


def test_collection_calculates_minimum_only_with_declared_metric(
    tmp_path: Path, monkeypatch
) -> None:
    root = tmp_path / "run"
    plan, jobs, _submission = _stage(root, monkeypatch, current_component="imag")
    _write_records(jobs)

    results = condor.collect_condor_sweep(root, plan=plan)

    assert [point["integrated_parallel_current"] for point in results.curve] == [8.0, 2.0, 6.0]
    assert results.resonance["status"] == "interior_minimum"
    assert results.resonance["profile_scale_factor"] == 0.0
    assert results.resonance["current_unit"] == "A m"


def test_collection_integrates_declared_namelist_current_column(
    tmp_path: Path, monkeypatch
) -> None:
    root = tmp_path / "run"
    plan, jobs, _submission = _stage(
        root,
        monkeypatch,
        backend="kim_x_namelist",
        metric=JparCurrentMetric(
            current_column=1,
            current_unit="statA/cm^2",
            collision_model="FokkerPlanck",
        ),
    )
    _write_records(jobs)
    for job, values in zip(jobs, ([3.0, 3.0, 3.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0])):
        fields = job.directory / "fields"
        fields.mkdir()
        np.savetxt(
            fields / "jpar.dat",
            np.column_stack(([1.0, 2.0, 3.0], values)),
        )

    results = condor.collect_condor_sweep(root, plan=plan)

    assert [point["integrated_parallel_current"] for point in results.curve] == [6.0, 2.0, 4.0]
    assert results.resonance["profile_scale_factor"] == 0.0
    assert all(point["current_unit"] == "statA/cm^2" for point in results.curve)


@pytest.mark.parametrize("invalid_case", ["pending", "missing", "invalid_output", "missing_status"])
def test_collection_rejects_pending_missing_or_invalid_output(
    tmp_path: Path, monkeypatch, invalid_case: str
) -> None:
    root = tmp_path / invalid_case
    metric = (
        JparCurrentMetric(current_column=1, current_unit="A/cm", collision_model="FokkerPlanck")
        if invalid_case == "invalid_output"
        else None
    )
    plan, jobs, submission = _stage(
        root,
        monkeypatch,
        backend="kim_x_namelist" if metric else "kamel_kim_python",
        metric=metric,
        factors=(1.0,),
    )
    if invalid_case != "missing":
        _write_records(jobs, currents=(2 + 1j,))
    if invalid_case == "pending":
        entry = submission["jobs"][0]
        (root / condor.STATUS_NAME).write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "statuses": [
                        {
                            "job": entry["job"],
                            "job_identity": entry["job_identity"],
                            "scan_order_index": entry["scan_order_index"],
                            "profile_scale_factor": entry["profile_scale_factor"],
                            "cluster": entry["cluster"],
                            "terminal": False,
                            "failure_kind": None,
                            "job_status": "running",
                        }
                    ],
                }
            ),
            encoding="utf-8",
        )
    if invalid_case == "missing_status":
        (root / condor.STATUS_NAME).unlink()
    if invalid_case == "invalid_output":
        fields = jobs[0].directory / "fields"
        fields.mkdir()
        np.savetxt(
            fields / "jpar.dat",
            np.column_stack(([1.0, 1.0], [2.0, 3.0])),
        )

    with pytest.raises(CondorError):
        condor.collect_condor_sweep(root, plan=plan)
