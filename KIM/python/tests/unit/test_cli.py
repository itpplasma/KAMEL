from __future__ import annotations

import json
import runpy
import sys
from pathlib import Path

import pytest
from kim import BuiltinPlasma, PlasmaIsotope, SimulationConfig
from kim.cli import app
from typer.testing import CliRunner

FIXTURES = Path(__file__).parents[1] / "fixtures"
runner = CliRunner()


def configuration_file(tmp_path: Path) -> Path:
    profiles = tmp_path / "profiles"
    runpy.run_path(str(FIXTURES / "generate_profiles.py"))["generate_profiles"](profiles)
    config = SimulationConfig.electrostatic_periodic(
        profiles=profiles,
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
    path = tmp_path / "request.json"
    path.write_text(config.model_dump_json(indent=2), encoding="utf-8")
    return path


def fake_executable(tmp_path: Path) -> Path:
    source = (FIXTURES / "fake_kim.py").read_text(encoding="utf-8")
    executable = tmp_path / "KIM.x"
    executable.write_text(f"#!{sys.executable}\n" + source.split("\n", 1)[1], encoding="utf-8")
    executable.chmod(0o755)
    return executable


@pytest.mark.parametrize(
    "command",
    ["parameters", "validate", "run", "sweep", "status", "inspect", "result"],
)
def test_every_command_has_help(command: str) -> None:
    result = runner.invoke(app, [command, "--help"])

    assert result.exit_code == 0
    assert "Usage:" in result.stdout
    assert "Traceback" not in result.output


def test_parameters_supports_table_and_json_schema() -> None:
    table = runner.invoke(app, ["parameters"])
    schema = runner.invoke(app, ["parameters", "--format", "json-schema"])

    assert table.exit_code == 0
    assert "setup.m_mode" in table.stdout
    assert "periodic.n_rg" in table.stdout
    assert schema.exit_code == 0
    parsed = json.loads(schema.stdout)
    assert parsed["title"] == "SimulationConfig"
    assert "$defs" in parsed


def test_validate_accepts_json_and_profile_override(tmp_path: Path) -> None:
    config = configuration_file(tmp_path)

    result = runner.invoke(app, ["validate", str(config), "--format", "json"])

    assert result.exit_code == 0
    payload = json.loads(result.stdout)
    assert payload["valid"] is True
    assert payload["run_type"] == "electrostatic_periodic"
    assert payload["profile_points"] == 11


def test_validate_error_is_concise_and_has_no_traceback(tmp_path: Path) -> None:
    invalid = tmp_path / "invalid.json"
    invalid.write_text('{"unknown": true}', encoding="utf-8")

    result = runner.invoke(app, ["validate", str(invalid)])

    assert result.exit_code == 1
    assert "Error:" in result.stderr
    assert "unknown" in result.stderr
    assert "errors.pydantic.dev" not in result.stderr
    assert "Traceback" not in result.output


def run_once(tmp_path: Path) -> tuple[Path, dict[str, object]]:
    runs = tmp_path / "runs"
    result = runner.invoke(
        app,
        [
            "run",
            str(configuration_file(tmp_path)),
            "--executable",
            str(fake_executable(tmp_path)),
            "--runs-dir",
            str(runs),
            "--format",
            "json",
        ],
    )
    assert result.exit_code == 0, result.output
    return runs, json.loads(result.stdout)


def test_run_returns_machine_readable_identity_and_artifacts(tmp_path: Path) -> None:
    runs, payload = run_once(tmp_path)

    assert payload["status"] == "succeeded"
    assert (runs / str(payload["run_id"]) / "manifest.json").is_file()
    assert Path(str(payload["output_file"])).is_file()


def test_run_failure_uses_nonzero_exit_without_traceback(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_KIM_MODE", "nonzero")
    result = runner.invoke(
        app,
        [
            "run",
            str(configuration_file(tmp_path)),
            "--executable",
            str(fake_executable(tmp_path)),
            "--runs-dir",
            str(tmp_path / "runs"),
            "--format",
            "json",
        ],
    )

    assert result.exit_code == 3
    assert json.loads(result.stdout)["status"] == "failed"
    assert "Traceback" not in result.output


def test_status_and_inspect_read_existing_run_without_execution(tmp_path: Path) -> None:
    runs, execution = run_once(tmp_path)
    run_id = str(execution["run_id"])

    status = runner.invoke(app, ["status", run_id, "--runs-dir", str(runs), "--format", "json"])
    inspection = runner.invoke(
        app, ["inspect", run_id, "--runs-dir", str(runs), "--format", "json"]
    )

    assert json.loads(status.stdout) == {"run_id": run_id, "status": "succeeded"}
    manifest = json.loads(inspection.stdout)
    assert manifest["run_id"] == run_id
    assert manifest["status"] == "succeeded"
    assert manifest["executable_path"].endswith("KIM.x")


def test_result_lists_and_reads_datasets_as_json(tmp_path: Path) -> None:
    runs, execution = run_once(tmp_path)
    base = [str(execution["run_id"]), "--runs-dir", str(runs), "--format", "json"]

    listed = runner.invoke(app, ["result", *base, "--list"])
    read = runner.invoke(app, ["result", *base, "--dataset", "fields/Phi"])

    assert "fields/Phi" in json.loads(listed.stdout)["datasets"]
    payload = json.loads(read.stdout)
    assert payload["path"] == "fields/Phi"
    assert payload["value"] == {"real": [1.0, 1.0, 1.0], "imag": [1.0, 1.0, 1.0]}


def test_scalar_sweep_command_uses_validated_api(tmp_path: Path) -> None:
    runs = tmp_path / "runs"
    result = runner.invoke(
        app,
        [
            "sweep",
            str(configuration_file(tmp_path)),
            "--parameter",
            "periodic.n_rg",
            "--values",
            "8",
            "--values",
            "12",
            "--executable",
            str(fake_executable(tmp_path)),
            "--runs-dir",
            str(runs),
            "--format",
            "json",
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.stdout)
    assert payload["status"] == "succeeded"
    assert len(payload["child_run_ids"]) == 2


def test_er_profile_sweep_command_uses_typed_transform(tmp_path: Path) -> None:
    result = runner.invoke(
        app,
        [
            "sweep",
            str(configuration_file(tmp_path)),
            "--scale-profile",
            "Er",
            "--values",
            "-1",
            "--values",
            "0",
            "--values",
            "1",
            "--executable",
            str(fake_executable(tmp_path)),
            "--runs-dir",
            str(tmp_path / "runs"),
            "--format",
            "json",
        ],
    )

    assert result.exit_code == 0, result.output
    assert len(json.loads(result.stdout)["child_run_ids"]) == 3


def test_sweep_requires_exactly_one_variation_kind(tmp_path: Path) -> None:
    result = runner.invoke(
        app,
        ["sweep", str(configuration_file(tmp_path)), "--values", "1"],
    )

    assert result.exit_code == 1
    assert "exactly one" in result.stderr
    assert "Traceback" not in result.output


def test_missing_run_error_is_concise(tmp_path: Path) -> None:
    result = runner.invoke(app, ["inspect", "missing", "--runs-dir", str(tmp_path / "runs")])

    assert result.exit_code == 1
    assert "run not found" in result.stderr
    assert "Traceback" not in result.output
