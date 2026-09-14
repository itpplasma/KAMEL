"""Smoke-test the packaged example and CLI from an installed KIM wheel."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from kim import ProfileSet, SimulationConfig, create_example


def validate_case(request_path: Path) -> None:
    """Validate a copied request while running from its case directory."""

    case_directory = request_path.parent
    previous_directory = Path.cwd()
    try:
        os.chdir(case_directory)
        config = SimulationConfig.model_validate_json(
            Path("request.json").read_text(encoding="utf-8")
        )
        validation = ProfileSet.from_simulation(config).validate_for(config)
    finally:
        os.chdir(previous_directory)
    assert validation.point_count == 65


def validate_with_cli(cli: Path, request_path: Path) -> None:
    """Validate a copied request through the installed console entry point."""

    result = subprocess.run(
        [str(cli), "validate", "request.json", "--format", "json"],
        cwd=request_path.parent,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise AssertionError(f"installed kim validate failed:\n{result.stdout}\n{result.stderr}")
    assert json.loads(result.stdout)["valid"] is True


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="kim-installed-example-") as temporary_directory:
        root = Path(temporary_directory)
        api_case = root / "api-case"
        request_path = create_example(api_case)
        validate_case(request_path)

        cli_case = root / "cli-case"
        cli = Path(sys.executable).with_name("kim")
        result = subprocess.run(
            [str(cli), "init", str(cli_case), "--example", "periodic", "--format", "json"],
            cwd=root,
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise AssertionError(f"installed kim init failed:\n{result.stdout}\n{result.stderr}")
        payload = json.loads(result.stdout)
        assert payload["example"] == "periodic"
        assert Path(payload["request_path"]) == cli_case / "request.json"
        validate_case(cli_case / "request.json")
        validate_with_cli(cli, cli_case / "request.json")


if __name__ == "__main__":
    main()
