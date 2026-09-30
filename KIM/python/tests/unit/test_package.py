from pathlib import Path

import pytest
from kim import SimulationConfig


def test_package_import_and_cli_help() -> None:
    import kim
    from kim.cli import app
    from typer.testing import CliRunner

    assert kim.__version__ == "0.1.0"

    result = CliRunner().invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "Run and inspect KIM plasma simulations." in result.stdout


def test_package_exports_condor_api() -> None:
    import kim

    for name in (
        "KimCondorPlan",
        "JparCurrentMetric",
        "KimCondorJob",
        "KimCondorScanResults",
        "stage_condor_sweep",
        "submit_condor_sweep",
        "wait_condor_sweep",
        "collect_condor_sweep",
    ):
        assert hasattr(kim, name)


@pytest.mark.parametrize(
    ("filename", "run_type"),
    [
        ("request-periodic.json", "electrostatic_periodic"),
        ("request-electrostatic.json", "electrostatic"),
        ("request-flr2.json", "flr2"),
    ],
)
def test_documented_request_examples_are_valid(filename: str, run_type: str) -> None:
    examples = Path(__file__).parents[2] / "examples"

    config = SimulationConfig.model_validate_json((examples / filename).read_text(encoding="utf-8"))

    assert config.run.run_type == run_type
