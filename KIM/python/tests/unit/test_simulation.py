from __future__ import annotations

import json
import runpy
import sys
import time
from pathlib import Path

import f90nml
import h5py
import pytest
from kim import (
    BuiltinPlasma,
    ElectrostaticRun,
    Flr2Run,
    PlasmaIsotope,
    RunStatus,
    Simulation,
    SimulationConfig,
)

FIXTURES = Path(__file__).parents[1] / "fixtures"


def configuration(profile_directory: Path) -> SimulationConfig:
    runpy.run_path(str(FIXTURES / "generate_profiles.py"))["generate_profiles"](profile_directory)
    return SimulationConfig.electrostatic_periodic(
        profiles=profile_directory,
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


def simulation(tmp_path: Path, **kwargs: object) -> Simulation:
    return Simulation(
        configuration(tmp_path / "source-profiles"),
        executable=fake_executable(tmp_path),
        runs_directory=tmp_path / "runs",
        **kwargs,
    )


def test_prepare_stages_validated_inputs_and_execution_contract(tmp_path: Path) -> None:
    prepared = simulation(tmp_path).prepare()

    assert prepared.manifest.status is RunStatus.PREPARED
    assert prepared.command == (
        str(prepared.manifest.executable_path),
        str(prepared.paths.namelist),
    )
    assert prepared.output_file == prepared.paths.results / "m7_n2" / "out_ES_periodic.h5"
    assert {path.name for path in prepared.paths.profiles.iterdir()} == {
        "n.dat",
        "Te.dat",
        "Ti.dat",
        "q.dat",
        "Er.dat",
        "Vz.dat",
    }
    nml = f90nml.read(prepared.paths.namelist)
    assert Path(nml["kim_io"]["profile_location"]) == prepared.paths.profiles
    assert Path(nml["kim_io"]["output_path"]) == prepared.paths.results
    assert nml["kim_io"]["h5_out_file"] == "out_ES_periodic.h5"
    assert {item.name for item in prepared.manifest.inputs} == {
        "request",
        "namelist",
        "n",
        "Te",
        "Ti",
        "q",
        "Er",
        "Vz",
    }


def test_run_uses_exact_command_and_run_directory_and_separate_logs(tmp_path: Path) -> None:
    result = simulation(tmp_path).run(environment={"OMP_NUM_THREADS": "3"})

    assert result.status is RunStatus.SUCCEEDED
    invocation = json.loads((result.run_directory / "invocation.json").read_text())
    assert invocation == {
        "argv": [str(result.manifest.executable_path), str(result.namelist)],
        "cwd": str(result.run_directory),
    }
    assert result.stdout.read_text() == "fake KIM stdout\n"
    assert result.stderr.read_text() == ""
    assert result.output_file.is_file()
    assert result.manifest.exit_code == 0
    assert result.manifest.environment == {"OMP_NUM_THREADS": "3"}
    assert result.manifest.discovered_outputs == (Path("results/m7_n2/out_ES_periodic.h5"),)


def test_stderr_is_preserved_on_success(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("FAKE_KIM_MODE", "stderr")
    result = simulation(tmp_path).run()

    assert result.status is RunStatus.SUCCEEDED
    assert result.stderr.read_text() == "fake KIM diagnostic\n"


def test_nonzero_exit_records_failure_and_retains_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_KIM_MODE", "nonzero")
    result = simulation(tmp_path).run()

    assert result.status is RunStatus.FAILED
    assert result.manifest.exit_code == 7
    assert result.manifest.failure is not None
    assert result.manifest.failure.kind == "process_exit"
    assert "fake KIM failure" in result.stderr.read_text()
    assert result.namelist.is_file()
    assert (result.run_directory / "inputs/profiles/n.dat").is_file()


@pytest.mark.parametrize("mode", ["missing_output", "partial_output"])
def test_zero_exit_requires_complete_periodic_hdf5(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode: str
) -> None:
    monkeypatch.setenv("FAKE_KIM_MODE", mode)
    result = simulation(tmp_path).run()

    assert result.status is RunStatus.FAILED
    assert result.manifest.exit_code == 0
    assert result.manifest.failure is not None
    expected_kind = "missing_output" if mode == "missing_output" else "invalid_output"
    assert result.manifest.failure.kind == expected_kind


def test_explicit_timeout_records_state_and_partial_logs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_KIM_MODE", "timeout")
    started = time.monotonic()
    result = simulation(tmp_path, timeout=0.5).run()

    assert time.monotonic() - started < 3.0
    assert result.status is RunStatus.TIMED_OUT
    assert result.manifest.exit_code is None
    assert result.manifest.failure is not None
    assert result.manifest.failure.kind == "timeout"
    assert result.manifest.failure.details == {"timeout_seconds": 0.5}
    assert "fake KIM waiting" in result.stderr.read_text()
    assert result.namelist.is_file()


def test_timeout_defaults_to_none(tmp_path: Path) -> None:
    assert simulation(tmp_path).timeout is None


@pytest.mark.parametrize(
    ("run", "output_name", "potential_dataset"),
    [
        (ElectrostaticRun(), "out_ES.h5", "fields/Phi_m"),
        (Flr2Run(), "out_FLR2.h5", "fields/Phi"),
    ],
)
def test_stable_run_types_have_explicit_output_contracts(
    tmp_path: Path,
    run: ElectrostaticRun | Flr2Run,
    output_name: str,
    potential_dataset: str,
) -> None:
    config = configuration(tmp_path / "source-profiles").model_copy(update={"run": run})
    result = Simulation(
        config,
        executable=fake_executable(tmp_path),
        runs_directory=tmp_path / "runs",
    ).run()

    assert result.status is RunStatus.SUCCEEDED
    assert result.output_file.name == output_name
    with h5py.File(result.output_file) as handle:
        assert potential_dataset in handle


def test_environment_rejects_arbitrary_process_controls(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="unsupported execution environment"):
        simulation(tmp_path).run(environment={"LD_PRELOAD": "/tmp/library"})
