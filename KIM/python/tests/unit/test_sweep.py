from __future__ import annotations

import json
import runpy
import sys
from pathlib import Path

import numpy as np
import pytest
from kim import (
    BuiltinPlasma,
    LinearRange,
    ParameterSweep,
    PlasmaIsotope,
    ProfileScale,
    RunStatus,
    SimulationConfig,
    SweepError,
    SweepSpec,
    run_sweep,
)

FIXTURES = Path(__file__).parents[1] / "fixtures"


def configuration(directory: Path) -> SimulationConfig:
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


def fake_executable(tmp_path: Path) -> Path:
    source = (FIXTURES / "fake_kim.py").read_text(encoding="utf-8")
    executable = tmp_path / "KIM.x"
    executable.write_text(f"#!{sys.executable}\n" + source.split("\n", 1)[1], encoding="utf-8")
    executable.chmod(0o755)
    return executable


def test_linear_range_includes_endpoints_and_is_deterministic() -> None:
    values = LinearRange(start=512, stop=2048, count=4)

    assert values.generate() == (512.0, 1024.0, 1536.0, 2048.0)


def test_explicit_scalar_values_resolve_public_periodic_schema_path(tmp_path: Path) -> None:
    spec = ParameterSweep(
        base=configuration(tmp_path / "profiles"),
        parameter="periodic.n_rg",
        values=[8, 12, 16],
    )

    configs = spec.configurations()

    assert [item.run.periodic.n_rg for item in configs] == [8, 12, 16]
    assert spec.values == (8, 12, 16)


def test_scalar_variation_does_not_mutate_base_configuration(tmp_path: Path) -> None:
    base = configuration(tmp_path / "profiles")
    original = base.grid.l_space_dim
    spec = ParameterSweep(base=base, parameter="grid.l_space_dim", values=[50, 60])

    configs = spec.configurations()

    assert base.grid.l_space_dim == original
    assert [item.grid.l_space_dim for item in configs] == [50, 60]


@pytest.mark.parametrize(
    ("parameter", "message"),
    [
        ("profiles.directory", "not sweepable"),
        ("periodic.missing", "unknown sweep parameter"),
        ("missing.value", "unknown sweep parameter"),
    ],
)
def test_unknown_and_unsweepable_paths_are_rejected(
    tmp_path: Path, parameter: str, message: str
) -> None:
    with pytest.raises(SweepError, match=message):
        ParameterSweep(
            base=configuration(tmp_path / parameter.replace(".", "-")),
            parameter=parameter,
            values=[1],
        )


def test_invalid_integer_sweep_value_is_rejected_before_execution(tmp_path: Path) -> None:
    spec = ParameterSweep(
        base=configuration(tmp_path / "profiles"),
        parameter="periodic.n_rg",
        values=[8, 8.5],
    )

    with pytest.raises(SweepError, match=r"periodic\.n_rg.*8\.5"):
        spec.configurations()


def test_scalar_sweep_runs_sequentially_and_records_ordered_children(tmp_path: Path) -> None:
    runs = tmp_path / "runs"
    result = run_sweep(
        ParameterSweep(
            base=configuration(tmp_path / "profiles"),
            parameter="periodic.n_rg",
            values=[8, 12, 16],
        ),
        executable=fake_executable(tmp_path),
        runs_directory=runs,
        label="grid scan",
    )

    assert result.status == "succeeded"
    assert [child.status for child in result.children] == [RunStatus.SUCCEEDED] * 3
    child_values = [
        child.manifest.normalized_parameters["run"]["periodic"]["n_rg"] for child in result.children
    ]
    assert child_values == [8, 12, 16]
    assert [child.manifest.parent_sweep_id for child in result.children] == [result.sweep_id] * 3
    persisted = json.loads(result.manifest.read_text())
    assert persisted["child_run_ids"] == [child.run_id for child in result.children]
    assert persisted["values"] == [8, 12, 16]
    assert persisted["base_parameters"]["run"]["periodic"]["n_rg"] == 96
    assert persisted["status"] == "succeeded"


def test_sweep_continues_after_child_failure_and_reports_partial_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_KIM_FAIL_N_RG", "12")
    result = run_sweep(
        ParameterSweep(
            base=configuration(tmp_path / "profiles"),
            parameter="periodic.n_rg",
            values=[8, 12, 16],
        ),
        executable=fake_executable(tmp_path),
        runs_directory=tmp_path / "runs",
    )

    assert [child.status for child in result.children] == [
        RunStatus.SUCCEEDED,
        RunStatus.FAILED,
        RunStatus.SUCCEEDED,
    ]
    assert result.status == "partial_failure"


def test_stop_on_failure_omits_unstarted_children(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_KIM_FAIL_N_RG", "12")
    result = run_sweep(
        ParameterSweep(
            base=configuration(tmp_path / "profiles"),
            parameter="periodic.n_rg",
            values=[8, 12, 16],
            continue_on_failure=False,
        ),
        executable=fake_executable(tmp_path),
        runs_directory=tmp_path / "runs",
    )

    assert len(result.children) == 2
    assert result.status == "partial_failure"


def test_er_profile_scale_preserves_source_and_stages_each_transformation(
    tmp_path: Path,
) -> None:
    base = configuration(tmp_path / "profiles")
    source = np.loadtxt(base.profiles.directory / "Er.dat")
    spec = SweepSpec(
        base=base,
        variation=ProfileScale(profile="Er", values=[-1.0, 0.0, 2.0]),
    )

    result = run_sweep(
        spec,
        executable=fake_executable(tmp_path),
        runs_directory=tmp_path / "runs",
    )

    np.testing.assert_array_equal(np.loadtxt(base.profiles.directory / "Er.dat"), source)
    assert [child.status for child in result.children] == [RunStatus.SUCCEEDED] * 3
    for child, factor in zip(result.children, (-1.0, 0.0, 2.0), strict=True):
        staged = np.loadtxt(child.run_directory / "inputs/profiles/Er.dat")
        np.testing.assert_allclose(staged[:, 0], source[:, 0], rtol=0.0, atol=0.0)
        np.testing.assert_allclose(staged[:, 1], factor * source[:, 1])
        assert child.manifest.requested_parameters["variation"] == {
            "kind": "profile_scale",
            "profile": "Er",
            "factor": factor,
        }


def test_profile_scale_rejects_unsupported_profiles_and_nonfinite_values(
    tmp_path: Path,
) -> None:
    with pytest.raises(ValueError, match="Er"):
        ProfileScale(profile="n", values=[1.0])
    with pytest.raises(ValueError, match="finite"):
        ProfileScale(profile="Er", values=[float("nan")])
