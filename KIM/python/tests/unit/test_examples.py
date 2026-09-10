from __future__ import annotations

import json
from importlib import resources
from pathlib import Path

import pytest
from kim import ConfigurationError, ProfileSet, SimulationConfig, create_example, examples

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


def test_create_example_copies_a_complete_independent_case(tmp_path: Path) -> None:
    destination = tmp_path / "case with spaces"

    request_path = create_example(destination)

    assert request_path == destination / "request.json"
    assert sorted(path.relative_to(destination).as_posix() for path in destination.rglob("*")) == [
        "README.md",
        "profiles",
        "profiles/Er.dat",
        "profiles/Te.dat",
        "profiles/Ti.dat",
        "profiles/Vz.dat",
        "profiles/n.dat",
        "profiles/q.dat",
        "request.json",
    ]
    request = json.loads(request_path.read_text(encoding="utf-8"))
    assert request["profiles"]["directory"] == "./profiles"

    second = tmp_path / "second"
    create_example(second)
    (second / "profiles" / "n.dat").write_text("changed\n", encoding="utf-8")
    assert (destination / "profiles" / "n.dat").read_text(encoding="utf-8") != "changed\n"


@pytest.mark.parametrize("kind", ["directory", "file", "symlink"])
def test_create_example_refuses_existing_destinations(tmp_path: Path, kind: str) -> None:
    destination = tmp_path / "case"
    if kind == "directory":
        destination.mkdir()
    elif kind == "file":
        destination.write_text("existing\n", encoding="utf-8")
    else:
        destination.symlink_to(tmp_path / "missing")

    with pytest.raises(ConfigurationError, match="destination already exists"):
        create_example(destination)


def test_create_example_rejects_unknown_name(tmp_path: Path) -> None:
    with pytest.raises(ConfigurationError, match="unknown example"):
        create_example(tmp_path / "case", name="missing")


def test_create_example_cleans_partial_copy_and_preserves_unrelated_files(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    destination = tmp_path / "case"
    unrelated = tmp_path / "unrelated.txt"
    unrelated.write_text("keep\n", encoding="utf-8")
    original_copy = examples.shutil.copyfileobj
    calls = 0

    def fail_during_copy(source: object, target: object, *args: object, **kwargs: object) -> None:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("simulated copy failure")
        original_copy(source, target, *args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(examples.shutil, "copyfileobj", fail_during_copy)

    with pytest.raises(OSError, match="simulated copy failure"):
        create_example(destination)

    assert not destination.exists()
    assert unrelated.read_text(encoding="utf-8") == "keep\n"


def test_create_example_preserves_entry_added_during_failed_copy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    destination = tmp_path / "case"
    original_copy = examples.shutil.copyfileobj
    calls = 0

    def fail_after_unrelated_entry(
        source: object, target: object, *args: object, **kwargs: object
    ) -> None:
        nonlocal calls
        calls += 1
        if calls == 1:
            (destination / "unrelated.txt").write_text("keep\n", encoding="utf-8")
        if calls == 2:
            raise OSError("simulated copy failure")
        original_copy(source, target, *args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(examples.shutil, "copyfileobj", fail_after_unrelated_entry)

    with pytest.raises(OSError, match="simulated copy failure"):
        create_example(destination)

    assert (destination / "unrelated.txt").read_text(encoding="utf-8") == "keep\n"


def test_create_example_nested_failure_removes_owned_ancestors(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    ancestor = tmp_path / "ancestor"
    ancestor.mkdir()
    unrelated = ancestor / "keep.txt"
    unrelated.write_text("keep\n", encoding="utf-8")
    destination = ancestor / "created" / "nested" / "case"

    def fail_copy(source: object, target: object, *args: object, **kwargs: object) -> None:
        raise OSError("simulated copy failure")

    monkeypatch.setattr(examples.shutil, "copyfileobj", fail_copy)

    with pytest.raises(OSError, match="simulated copy failure"):
        create_example(destination)

    assert not destination.exists()
    assert not (ancestor / "created").exists()
    assert unrelated.read_text(encoding="utf-8") == "keep\n"
