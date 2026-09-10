from __future__ import annotations

import hashlib
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest
from kim.config import BuiltinPlasma, PlasmaIsotope, SimulationConfig
from kim.errors import RunError
from kim.executable import KamelGitMetadata
from kim.runs import RunFailure, RunRepository, RunStatus


def configuration(profile_directory: Path) -> SimulationConfig:
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
        radial_minimum=3.0,
        plasma_radius=63.0,
    )


class Clock:
    def __init__(self) -> None:
        self.current = datetime(2026, 9, 9, 12, 30, tzinfo=timezone.utc)

    def __call__(self) -> datetime:
        value = self.current
        self.current += timedelta(seconds=2)
        return value


def test_create_allocates_unique_run_ids_and_complete_layout(tmp_path: Path) -> None:
    repository = RunRepository(tmp_path / "runs", clock=Clock())
    config = configuration(tmp_path / "source-profiles")

    first = repository.create(config, label="Periodic Er Scan")
    second = repository.create(config, label="Periodic Er Scan")

    assert first.run_id != second.run_id
    assert first.run_id.startswith("20260909T123000Z-periodic-er-scan-")
    for run in (first, second):
        assert run.root.is_dir()
        assert run.manifest == run.root / "manifest.json"
        assert run.request == run.root / "inputs" / "request.json"
        assert run.profiles == run.root / "inputs" / "profiles"
        assert run.stdout == run.root / "logs" / "stdout.log"
        assert run.stderr == run.root / "logs" / "stderr.log"
        assert run.results == run.root / "results"
        assert run.request.is_file()
        assert run.profiles.is_dir()
        assert run.results.is_dir()


def test_initial_manifest_is_versioned_and_records_normalized_request(tmp_path: Path) -> None:
    repository = RunRepository(tmp_path / "runs", clock=Clock())
    config = configuration(tmp_path / "profiles")
    run = repository.create(config, label="test", parent_sweep_id="sweep-123")

    manifest = repository.inspect(run.run_id)

    assert manifest.schema_version == 1
    assert manifest.python_package_version == "0.1.0"
    assert manifest.status is RunStatus.PREPARED
    assert manifest.label == "test"
    assert manifest.parent_sweep_id == "sweep-123"
    assert manifest.requested_parameters == config.model_dump(mode="json")
    assert manifest.normalized_parameters == config.model_dump(mode="json")
    request = next(item for item in manifest.inputs if item.name == "request")
    assert request.source_path is None
    assert request.staged_path == Path("inputs/request.json")
    assert request.sha256 == hashlib.sha256(run.request.read_bytes()).hexdigest()


def test_execution_provenance_is_recorded_while_prepared(tmp_path: Path) -> None:
    repository = RunRepository(tmp_path / "runs", clock=Clock())
    run = repository.create(configuration(tmp_path / "profiles"))
    executable = tmp_path / "KIM.x"
    executable.write_bytes(b"solver")
    git = KamelGitMetadata(root=tmp_path, commit="a" * 40, dirty=True)

    updated = repository.configure_execution(
        run.run_id,
        executable_path=executable,
        executable_sha256="b" * 64,
        arguments=(str(executable), "inputs/KIM_config.nml"),
        working_directory=run.root,
        kamel_git=git,
        environment={"OMP_NUM_THREADS": "4"},
        expected_outputs=(Path("results/m7_n2/out_ES_periodic.h5"),),
    )

    assert updated.executable_path == executable.absolute()
    assert updated.executable_sha256 == "b" * 64
    assert updated.arguments == (str(executable), "inputs/KIM_config.nml")
    assert updated.working_directory == run.root.absolute()
    assert updated.kamel_git_commit == "a" * 40
    assert updated.kamel_git_dirty is True
    assert updated.environment == {"OMP_NUM_THREADS": "4"}
    assert updated.expected_outputs == (Path("results/m7_n2/out_ES_periodic.h5"),)


def test_input_artifact_digest_and_paths_are_recorded(tmp_path: Path) -> None:
    repository = RunRepository(tmp_path / "runs", clock=Clock())
    run = repository.create(configuration(tmp_path / "profiles"))
    source = tmp_path / "source" / "n.dat"
    source.parent.mkdir()
    source.write_text("0 1e13\n1 9e12\n")
    staged = run.profiles / "n.dat"
    staged.write_bytes(source.read_bytes())

    artifact = repository.record_input(
        run.run_id,
        name="n",
        source_path=source,
        staged_path=staged,
    )

    assert artifact.source_path == source.absolute()
    assert artifact.staged_path == Path("inputs/profiles/n.dat")
    assert artifact.sha256 == hashlib.sha256(staged.read_bytes()).hexdigest()
    assert repository.inspect(run.run_id).inputs[-1] == artifact


def test_state_transitions_record_timing_and_result_details(tmp_path: Path) -> None:
    repository = RunRepository(tmp_path / "runs", clock=Clock())
    run = repository.create(configuration(tmp_path / "profiles"))

    running = repository.transition(run.run_id, RunStatus.RUNNING)
    succeeded = repository.transition(
        run.run_id,
        RunStatus.SUCCEEDED,
        exit_code=0,
        discovered_outputs=(Path("results/result.h5"),),
    )

    assert running.started_at is not None
    assert succeeded.finished_at is not None
    assert succeeded.duration_seconds == pytest.approx(2.0)
    assert succeeded.exit_code == 0
    assert succeeded.discovered_outputs == (Path("results/result.h5"),)


@pytest.mark.parametrize(
    "terminal",
    [RunStatus.FAILED, RunStatus.TIMED_OUT, RunStatus.INTERRUPTED],
)
def test_running_run_can_reach_each_non_success_terminal_state(
    tmp_path: Path, terminal: RunStatus
) -> None:
    repository = RunRepository(tmp_path / terminal.value, clock=Clock())
    run = repository.create(configuration(tmp_path / "profiles"))
    repository.transition(run.run_id, RunStatus.RUNNING)

    manifest = repository.transition(
        run.run_id,
        terminal,
        exit_code=3 if terminal is RunStatus.FAILED else None,
        failure=RunFailure(kind=terminal.value, message="simulated failure"),
    )

    assert manifest.status is terminal
    assert manifest.failure is not None
    assert manifest.failure.message == "simulated failure"


def test_prepared_run_can_be_marked_interrupted_before_launch(tmp_path: Path) -> None:
    repository = RunRepository(tmp_path / "runs", clock=Clock())
    run = repository.create(configuration(tmp_path / "profiles"))

    manifest = repository.transition(
        run.run_id,
        RunStatus.INTERRUPTED,
        failure=RunFailure(kind="interrupted", message="staging interrupted"),
    )

    assert manifest.started_at is None
    assert manifest.finished_at is not None
    assert manifest.duration_seconds is None


@pytest.mark.parametrize(
    ("source", "destination"),
    [
        (RunStatus.PREPARED, RunStatus.SUCCEEDED),
        (RunStatus.PREPARED, RunStatus.TIMED_OUT),
        (RunStatus.RUNNING, RunStatus.PREPARED),
    ],
)
def test_invalid_state_transitions_are_rejected(
    tmp_path: Path, source: RunStatus, destination: RunStatus
) -> None:
    repository = RunRepository(tmp_path / f"{source.value}-{destination.value}", clock=Clock())
    run = repository.create(configuration(tmp_path / "profiles"))
    if source is RunStatus.RUNNING:
        repository.transition(run.run_id, RunStatus.RUNNING)

    with pytest.raises(RunError, match=r"invalid run transition"):
        repository.transition(run.run_id, destination)


def test_completed_runs_reject_all_repository_mutation(tmp_path: Path) -> None:
    repository = RunRepository(tmp_path / "runs", clock=Clock())
    run = repository.create(configuration(tmp_path / "profiles"))
    repository.transition(run.run_id, RunStatus.RUNNING)
    repository.transition(run.run_id, RunStatus.SUCCEEDED, exit_code=0)
    staged = run.root / "extra.input"
    staged.write_text("data")

    with pytest.raises(RunError, match="immutable"):
        repository.transition(run.run_id, RunStatus.FAILED)
    with pytest.raises(RunError, match="immutable"):
        repository.record_input(run.run_id, name="extra", source_path=staged, staged_path=staged)
    with pytest.raises(RunError, match="immutable"):
        repository.configure_execution(
            run.run_id,
            executable_path=tmp_path / "KIM.x",
            executable_sha256="a" * 64,
            arguments=(),
            working_directory=run.root,
        )


def test_listing_is_newest_first_and_supports_status_filter(tmp_path: Path) -> None:
    repository = RunRepository(tmp_path / "runs", clock=Clock())
    first = repository.create(configuration(tmp_path / "profiles"), label="first")
    second = repository.create(configuration(tmp_path / "profiles"), label="second")
    repository.transition(first.run_id, RunStatus.RUNNING)
    repository.transition(
        first.run_id,
        RunStatus.FAILED,
        exit_code=1,
        failure=RunFailure(kind="process_exit", message="exit status 1"),
    )

    all_runs = repository.list_runs()
    prepared = repository.list_runs(status=RunStatus.PREPARED)

    assert [run.run_id for run in all_runs] == [second.run_id, first.run_id]
    assert [run.run_id for run in prepared] == [second.run_id]
    assert repository.inspect(first.run_id).status is RunStatus.FAILED


def test_unknown_or_unsafe_run_id_is_rejected(tmp_path: Path) -> None:
    repository = RunRepository(tmp_path / "runs")

    with pytest.raises(RunError, match="run not found"):
        repository.inspect("missing")
    with pytest.raises(RunError, match="invalid run ID"):
        repository.inspect("../outside")


def test_manifest_run_id_must_match_its_directory(tmp_path: Path) -> None:
    repository = RunRepository(tmp_path / "runs", clock=Clock())
    run = repository.create(configuration(tmp_path / "profiles"))
    contents = run.manifest.read_text().replace(run.run_id, "different-run", 1)
    run.manifest.write_text(contents)

    with pytest.raises(RunError, match="does not match directory"):
        repository.inspect(run.run_id)


def test_interrupted_atomic_write_leaves_last_complete_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    repository = RunRepository(tmp_path / "runs", clock=Clock())
    run = repository.create(configuration(tmp_path / "profiles"))
    before = run.manifest.read_bytes()

    def interrupted_replace(source: Path, destination: Path) -> None:
        raise OSError("simulated interruption")

    monkeypatch.setattr("kim.runs.os.replace", interrupted_replace)
    with pytest.raises(RunError, match="atomic manifest write failed"):
        repository.transition(run.run_id, RunStatus.RUNNING)

    assert run.manifest.read_bytes() == before
    assert repository.inspect(run.run_id).status is RunStatus.PREPARED
    assert list(run.root.glob(".manifest.json.*.tmp")) == []
