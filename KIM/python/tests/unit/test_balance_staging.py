from __future__ import annotations

import hashlib
import json
import subprocess
from dataclasses import FrozenInstanceError
from pathlib import Path
from types import MappingProxyType

import numpy as np
import pytest
from kim import (
    BalanceInput,
    BalanceMetadata,
    ExperimentalInputError,
    balance_adoption,
    read_balance_profiles,
    read_marsf_profiles,
)
from kim.balance_adoption import stage_balance_marsf_quartet
from kim.importers.experimental import ExperimentalProfile

_PROFILE_ROLES = (
    "density",
    "electron_temperature",
    "ion_temperature",
    "toroidal_rotation",
)

_MARSF_FILENAMES = {
    "density": "PROFDEN.IN",
    "electron_temperature": "PROFTE.IN",
    "ion_temperature": "PROFTI.IN",
    "toroidal_rotation": "PROFROT.IN",
}

_R0_CM = 165.0
_EQUILIBRIUM_PROVENANCE = "synthetic-equilibrium-for-staging-test"


def metadata(**updates: object) -> BalanceMetadata:
    values: dict[str, object] = {
        "source": "synthetic-balance-staging-test-case",
        "coordinate": "rho_pol",
        "coordinate_unit": "1",
        "density_unit": "1/m^3",
        "electron_temperature_unit": "eV",
        "ion_temperature_unit": "eV",
        "toroidal_rotation_unit": "rad/s",
    }
    values.update(updates)
    return BalanceMetadata(**values)


def write_balance_profiles(directory: Path) -> dict[str, Path]:
    directory.mkdir()
    grids = {
        "density": [0.0, 0.5, 1.0],
        "electron_temperature": [0.0, 0.4, 1.0],
        "ion_temperature": [0.0, 0.3, 1.0],
        "toroidal_rotation": [0.0, 0.25, 1.0],
    }
    values = {
        "density": [1.0e19, 2.0e19, 3.0e19],
        "electron_temperature": [100.0, 80.0, 40.0],
        "ion_temperature": [90.0, 70.0, 30.0],
        "toroidal_rotation": [0.0, -2.5e4, -5.0e4],
    }
    paths: dict[str, Path] = {}
    for role in _PROFILE_ROLES:
        path = directory / f"{role}.dat"
        path.write_text(
            "\n".join(
                f"{coordinate:.8e} {value:.8e}"
                for coordinate, value in zip(grids[role], values[role], strict=True)
            )
            + "\n",
            encoding="utf-8",
        )
        paths[role] = path
    return paths


def read_case(paths: dict[str, Path]) -> BalanceInput:
    return read_balance_profiles(
        density=paths["density"],
        electron_temperature=paths["electron_temperature"],
        ion_temperature=paths["ion_temperature"],
        toroidal_rotation=paths["toroidal_rotation"],
        metadata=metadata(),
    )


def _read_rows(path: Path) -> np.ndarray:
    lines = path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 4
    assert lines[0].strip()
    with pytest.raises(ValueError):
        float(lines[0].split()[0])

    rows: list[tuple[float, float]] = []
    for line in lines[1:]:
        columns = line.split()
        assert len(columns) == 2
        coordinate, value = (float(column) for column in columns)
        assert np.isfinite(coordinate)
        assert np.isfinite(value)
        rows.append((coordinate, value))
    return np.asarray(rows, dtype=np.float64)


def _report(result: object) -> tuple[Path, dict[str, object]]:
    report = getattr(result, "report")
    assert isinstance(report, Path)
    assert report.is_file()
    return report, json.loads(report.read_text(encoding="utf-8"))


def test_stages_exact_marsf_quartet_with_explicit_metadata(tmp_path: Path) -> None:
    paths = write_balance_profiles(tmp_path / "balance")
    source = read_case(paths)

    result = stage_balance_marsf_quartet(
        source,
        tmp_path / "marsf",
        major_radius_cm=_R0_CM,
        equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
    )

    assert {path.name for path in result.directory.iterdir()} == {
        *_MARSF_FILENAMES.values(),
        "staging_report.json",
    }
    assert result.metadata.source == source.metadata.source
    assert result.metadata.coordinate == "sqrt_psiN"
    assert result.metadata.coordinate_unit == "1"
    assert result.metadata.density_unit == "1/m^3"
    assert result.metadata.electron_temperature_unit == "eV"
    assert result.metadata.ion_temperature_unit == "eV"
    assert result.metadata.toroidal_velocity_unit == "cm/s"
    assert result.metadata.equilibrium_provenance == _EQUILIBRIUM_PROVENANCE

    for filename in _MARSF_FILENAMES.values():
        rows = _read_rows(result.directory / filename)
        assert rows.shape == (3, 2)


def test_staged_quartet_is_readable_and_preserves_role_grids_and_values(
    tmp_path: Path,
) -> None:
    paths = write_balance_profiles(tmp_path / "balance")
    source = read_case(paths)
    source_bytes = {role: path.read_bytes() for role, path in paths.items()}
    source_arrays = {
        role: (
            profile.coordinate.copy(),
            profile.values.copy(),
        )
        for role, profile in source.profiles.items()
    }

    result = stage_balance_marsf_quartet(
        source,
        tmp_path / "marsf",
        major_radius_cm=_R0_CM,
        equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
    )
    staged = read_marsf_profiles(result.directory, result.metadata)

    for role in _PROFILE_ROLES:
        original = source.profiles[role]
        actual = staged.profiles["toroidal_velocity" if role == "toroidal_rotation" else role]
        np.testing.assert_array_equal(actual.coordinate, original.coordinate)
        expected_values = original.values
        if role == "toroidal_rotation":
            expected_values = expected_values * _R0_CM
        np.testing.assert_array_equal(actual.values, expected_values)

    assert staged.profiles["density"].units == "1/m^3"
    assert staged.profiles["electron_temperature"].units == "eV"
    assert staged.profiles["ion_temperature"].units == "eV"
    assert staged.profiles["toroidal_velocity"].units == "cm/s"
    assert {role: path.read_bytes() for role, path in paths.items()} == source_bytes
    for role, (coordinates, values) in source_arrays.items():
        np.testing.assert_array_equal(source.profiles[role].coordinate, coordinates)
        np.testing.assert_array_equal(source.profiles[role].values, values)


def test_staging_reports_source_and_derived_hashes_and_scientific_operations(
    tmp_path: Path,
) -> None:
    paths = write_balance_profiles(tmp_path / "balance")
    source = read_case(paths)

    result = stage_balance_marsf_quartet(
        source,
        tmp_path / "marsf",
        major_radius_cm=_R0_CM,
        equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
    )
    report_path, report = _report(result)

    assert report_path == result.directory / "staging_report.json"
    assert report["source_basenames"] == {role: path.name for role, path in paths.items()}
    expected_source_hashes = {
        role: hashlib.sha256(path.read_bytes()).hexdigest() for role, path in paths.items()
    }
    assert report["source_hashes"] == expected_source_hashes
    assert report["major_radius_cm"] == _R0_CM
    assert report["coordinate_mapping"] == {
        "source": "rho_pol",
        "target": "sqrt_psiN",
        "operation": "preserve",
    }
    assert report["equilibrium_provenance"] == _EQUILIBRIUM_PROVENANCE
    assert all(operation["quantity"] != "q" for operation in report["operations"])

    rotation_operation = next(
        operation
        for operation in report["operations"]
        if operation["quantity"] == "toroidal_velocity"
    )
    assert rotation_operation == {
        "quantity": "toroidal_velocity",
        "source_unit": "rad/s",
        "target_unit": "cm/s",
        "factor": _R0_CM,
        "operation": "omega_to_v_phi",
        "parameters": {"major_radius_cm": _R0_CM},
    }
    assert report["derived_hashes"] == {
        filename: hashlib.sha256((result.directory / filename).read_bytes()).hexdigest()
        for filename in _MARSF_FILENAMES.values()
    }
    assert str(tmp_path) not in report_path.read_text(encoding="utf-8")


def test_staging_result_is_immutable(tmp_path: Path) -> None:
    source = read_case(write_balance_profiles(tmp_path / "balance"))
    result = stage_balance_marsf_quartet(
        source,
        tmp_path / "marsf",
        major_radius_cm=_R0_CM,
        equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
    )

    with pytest.raises((FrozenInstanceError, AttributeError, TypeError)):
        result.directory = tmp_path / "other"  # type: ignore[misc]


@pytest.mark.parametrize("major_radius_cm", [None, np.nan, 0.0, -1.0])
def test_staging_requires_a_finite_positive_major_radius(
    tmp_path: Path, major_radius_cm: float | None
) -> None:
    source = read_case(write_balance_profiles(tmp_path / "balance"))

    with pytest.raises(ExperimentalInputError, match="major_radius_cm"):
        stage_balance_marsf_quartet(
            source,
            tmp_path / "marsf",
            major_radius_cm=major_radius_cm,
            equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
        )


def test_staging_requires_nonempty_equilibrium_provenance(tmp_path: Path) -> None:
    source = read_case(write_balance_profiles(tmp_path / "balance"))

    with pytest.raises(ExperimentalInputError, match="equilibrium_provenance"):
        stage_balance_marsf_quartet(
            source,
            tmp_path / "marsf",
            major_radius_cm=_R0_CM,
            equilibrium_provenance="",
        )


def test_staging_refuses_existing_destination_without_modifying_it(tmp_path: Path) -> None:
    source = read_case(write_balance_profiles(tmp_path / "balance"))
    destination = tmp_path / "marsf"
    destination.mkdir()
    sentinel = destination / "keep.txt"
    sentinel.write_text("keep", encoding="utf-8")

    with pytest.raises(ExperimentalInputError, match="already exists"):
        stage_balance_marsf_quartet(
            source,
            destination,
            major_radius_cm=_R0_CM,
            equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
        )

    assert sentinel.read_text(encoding="utf-8") == "keep"
    assert {path.name for path in destination.iterdir()} == {"keep.txt"}


def test_staging_refuses_existing_destination_symlink(tmp_path: Path) -> None:
    source = read_case(write_balance_profiles(tmp_path / "balance"))
    destination = tmp_path / "marsf"
    destination.symlink_to(tmp_path / "missing-target")

    with pytest.raises(ExperimentalInputError, match="already exists"):
        stage_balance_marsf_quartet(
            source,
            destination,
            major_radius_cm=_R0_CM,
            equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
        )

    assert destination.is_symlink()


def test_staging_failure_removes_partial_output_and_temporary_artifacts(tmp_path: Path) -> None:
    paths = write_balance_profiles(tmp_path / "balance")
    source = read_case(paths)
    invalid_profile = ExperimentalProfile(
        name="ion_temperature",
        units="eV",
        path=source.source_files["ion_temperature"],
        coordinate=np.asarray([0.0, 0.3, 0.2], dtype=np.float64),
        values=np.asarray([90.0, 70.0, 30.0], dtype=np.float64),
    )
    invalid_profiles = dict(source.profiles)
    invalid_profiles["ion_temperature"] = invalid_profile
    invalid_source = BalanceInput(
        metadata=source.metadata,
        profiles=MappingProxyType(invalid_profiles),
        source_files=source.source_files,
        source_hashes=source.source_hashes,
    )
    destination = tmp_path / "marsf"
    before = set(tmp_path.iterdir())

    with pytest.raises(ExperimentalInputError, match="strictly increasing"):
        stage_balance_marsf_quartet(
            invalid_source,
            destination,
            major_radius_cm=_R0_CM,
            equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
        )

    assert not destination.exists()
    assert set(tmp_path.iterdir()) == before


def test_staging_write_failure_removes_partial_output_and_temporary_artifacts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = read_case(write_balance_profiles(tmp_path / "balance"))
    destination = tmp_path / "marsf"
    before = set(tmp_path.iterdir())
    # This private seam mirrors preparation._write_profile and gives every quartet
    # write a deterministic failure-injection point without patching filesystem primitives.
    original_write_profile = balance_adoption._write_profile
    writes = 0

    def fail_after_first_write(*args: object, **kwargs: object) -> object:
        nonlocal writes
        writes += 1
        if writes == 2:
            raise OSError("simulated failure")
        return original_write_profile(*args, **kwargs)

    monkeypatch.setattr(balance_adoption, "_write_profile", fail_after_first_write)
    with pytest.raises((OSError, ExperimentalInputError)):
        stage_balance_marsf_quartet(
            source,
            destination,
            major_radius_cm=_R0_CM,
            equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
        )

    assert writes == 2
    assert not destination.exists()
    assert not list(tmp_path.glob(f".{destination.name}-*"))
    assert set(tmp_path.iterdir()) == before


def test_staging_does_not_launch_external_processes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = read_case(write_balance_profiles(tmp_path / "balance"))

    def fail_process_launch(*_args: object, **_kwargs: object) -> None:
        pytest.fail("BALANCE quartet staging must not launch an external process")

    for launcher in ("run", "Popen", "call", "check_call", "check_output"):
        monkeypatch.setattr(subprocess, launcher, fail_process_launch)
    stage_balance_marsf_quartet(
        source,
        tmp_path / "marsf",
        major_radius_cm=_R0_CM,
        equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
    )
