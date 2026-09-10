from __future__ import annotations

import json
from importlib import resources
from pathlib import Path

import pytest
from kim import ProfileSet, SimulationConfig

EXAMPLE = resources.files("kim.example_data").joinpath("periodic")
REFERENCE = Path(__file__).parents[1] / "fixtures" / "periodic_reference.json"


def test_packaged_periodic_request_validates_profiles_and_resonance(tmp_path: Path) -> None:
    request = json.loads(EXAMPLE.joinpath("request.json").read_text(encoding="utf-8"))

    profile_directory = tmp_path / "profiles"
    profile_directory.mkdir()
    for resource in EXAMPLE.joinpath("profiles").iterdir():
        (profile_directory / resource.name).write_bytes(resource.read_bytes())
    request["profiles"]["directory"] = str(profile_directory)
    config = SimulationConfig.model_validate(request)
    validation = ProfileSet.from_simulation(config).validate_for(config)

    expected = json.loads(REFERENCE.read_text(encoding="utf-8"))
    assert validation.point_count == 65
    assert validation.resonance_radius == pytest.approx(
        expected["expected"]["profile_linear_q_crossing_cm"],
        abs=expected["tolerances"]["resonance_absolute_cm"],
    )
    assert config.profiles.directory == profile_directory


def test_periodic_guide_json_agrees_with_packaged_explicit_values() -> None:
    guide = Path(__file__).parents[2] / "docs" / "request-json.md"
    snippet = guide.read_text(encoding="utf-8").split("```json", 1)[1].split("```", 1)[0]
    guide_request = json.loads(snippet)
    request = json.loads(EXAMPLE.joinpath("request.json").read_text(encoding="utf-8"))

    def assert_explicit_values(explicit: dict[str, object], packaged: dict[str, object]) -> None:
        for key, value in explicit.items():
            assert key in packaged
            if isinstance(value, dict):
                assert isinstance(packaged[key], dict)
                assert_explicit_values(value, packaged[key])
            else:
                assert packaged[key] == value

    assert_explicit_values(guide_request, request)
