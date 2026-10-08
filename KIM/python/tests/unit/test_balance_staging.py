from __future__ import annotations

import ctypes
import errno
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


def test_staging_uses_r_big_from_equilibrium_output(tmp_path: Path) -> None:
    source = read_case(write_balance_profiles(tmp_path / "balance"))
    parameters = tmp_path / "btor_rbig.dat"
    parameters.write_text("-17573.19212 169.9117661\n", encoding="utf-8")

    result = stage_balance_marsf_quartet(
        source,
        tmp_path / "marsf-from-equilibrium",
        equilibrium_parameters_file=parameters,
        equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
    )

    rotation = np.loadtxt(result.directory / "PROFROT.IN", skiprows=1)
    np.testing.assert_allclose(
        rotation[:, 1], source.profiles["toroidal_rotation"].values * 169.9117661
    )
    report = json.loads(result.report.read_text(encoding="utf-8"))
    assert report["schema_version"] == 2
    assert report["equilibrium_parameters"] == {
        "source_basename": "btor_rbig.dat",
        "sha256": hashlib.sha256(parameters.read_bytes()).hexdigest(),
        "btor_gauss": -17573.19212,
        "r_big_cm": 169.9117661,
    }


@pytest.mark.parametrize("major_radius_cm", [np.nan, True, "169.9117661"])
def test_equilibrium_radius_consistency_check_rejects_invalid_values(
    tmp_path: Path, major_radius_cm: object
) -> None:
    source = read_case(write_balance_profiles(tmp_path / "balance"))
    parameters = tmp_path / "btor_rbig.dat"
    parameters.write_text("-17573.19212 169.9117661\n", encoding="utf-8")

    with pytest.raises(ExperimentalInputError, match="major_radius_cm must be finite and positive"):
        stage_balance_marsf_quartet(
            source,
            tmp_path / "marsf-invalid-radius",
            equilibrium_parameters_file=parameters,
            major_radius_cm=major_radius_cm,
            equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
        )


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
        "source_unit": "1",
        "target": "sqrt_psiN",
        "target_unit": "1",
        "operation": "rho_pol = sqrt(psi_pol_norm)",
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


def test_staging_normalizes_equilibrium_provenance_before_metadata_and_report(
    tmp_path: Path,
) -> None:
    source = read_case(write_balance_profiles(tmp_path / "balance"))

    result = stage_balance_marsf_quartet(
        source,
        tmp_path / "marsf",
        major_radius_cm=_R0_CM,
        equilibrium_provenance="  synthetic-equilibrium  ",
    )

    assert result.metadata.equilibrium_provenance == "synthetic-equilibrium"
    report = json.loads(result.report.read_text(encoding="utf-8"))
    assert report["equilibrium_provenance"] == "synthetic-equilibrium"


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


def test_staging_rejects_source_file_changed_after_read(tmp_path: Path) -> None:
    paths = write_balance_profiles(tmp_path / "balance")
    source = read_case(paths)
    paths["density"].write_text(
        "0.0 1.0e19\n0.5 2.0e19\n1.0 4.0e19\n",
        encoding="utf-8",
    )
    destination = tmp_path / "marsf"

    with pytest.raises(ExperimentalInputError, match="(?i)(hash|changed)"):
        stage_balance_marsf_quartet(
            source,
            destination,
            major_radius_cm=_R0_CM,
            equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
        )

    assert not destination.exists()
    assert not list(tmp_path.glob(f".{destination.name}-*"))


def test_staging_refuses_concurrent_destination_reservation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = read_case(write_balance_profiles(tmp_path / "balance"))
    destination = tmp_path / "marsf"
    original_rename = balance_adoption._rename_directory_noreplace

    def create_concurrent_destination(staging: Path, path: Path) -> None:
        path.mkdir()
        original_rename(staging, path)

    monkeypatch.setattr(
        balance_adoption, "_rename_directory_noreplace", create_concurrent_destination
    )

    with pytest.raises(ExperimentalInputError, match="already exists"):
        stage_balance_marsf_quartet(
            source,
            destination,
            major_radius_cm=_R0_CM,
            equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
        )

    assert destination.is_dir()
    assert not list(destination.iterdir())
    assert not list(tmp_path.glob(f".{destination.name}-*"))


def _write_publication_tree(path: Path, marker: str) -> None:
    path.mkdir()
    (path / "marker.txt").write_text(marker, encoding="utf-8")


def test_exclusive_publication_publishes_one_complete_directory(tmp_path: Path) -> None:
    staging = tmp_path / ".staging"
    destination = tmp_path / "published"
    _write_publication_tree(staging, "winner")

    balance_adoption._publish_staging_directory(staging, destination)

    assert not staging.exists()
    assert (destination / "marker.txt").read_text(encoding="utf-8") == "winner"


def test_exclusive_publication_loser_preserves_complete_winner(
    tmp_path: Path,
) -> None:
    winner_staging = tmp_path / ".winner"
    loser_staging = tmp_path / ".loser"
    destination = tmp_path / "published"
    _write_publication_tree(winner_staging, "winner")
    _write_publication_tree(loser_staging, "loser")

    balance_adoption._publish_staging_directory(winner_staging, destination)
    with pytest.raises(ExperimentalInputError, match="already exists"):
        balance_adoption._publish_staging_directory(loser_staging, destination)

    assert (destination / "marker.txt").read_text(encoding="utf-8") == "winner"
    assert (loser_staging / "marker.txt").read_text(encoding="utf-8") == "loser"


def test_exclusive_publication_preserves_foreign_destination(
    tmp_path: Path,
) -> None:
    staging = tmp_path / ".staging"
    destination = tmp_path / "published"
    _write_publication_tree(staging, "staged")
    destination.mkdir()
    sentinel = destination / "foreign.txt"
    sentinel.write_text("foreign", encoding="utf-8")

    with pytest.raises(ExperimentalInputError, match="already exists"):
        balance_adoption._publish_staging_directory(staging, destination)

    assert sentinel.read_text(encoding="utf-8") == "foreign"
    assert (staging / "marker.txt").read_text(encoding="utf-8") == "staged"


def test_publication_error_cleans_only_owned_staging(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = read_case(write_balance_profiles(tmp_path / "balance"))
    destination = tmp_path / "marsf"
    foreign = tmp_path / "foreign"
    foreign.mkdir()
    sentinel = foreign / "keep.txt"
    sentinel.write_text("keep", encoding="utf-8")

    def fail_publication(*_args: object, **_kwargs: object) -> None:
        raise OSError("simulated publication failure")

    monkeypatch.setattr(balance_adoption, "_rename_directory_noreplace", fail_publication)
    with pytest.raises(ExperimentalInputError, match="publication"):
        stage_balance_marsf_quartet(
            source,
            destination,
            major_radius_cm=_R0_CM,
            equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
        )

    assert not destination.exists()
    assert sentinel.read_text(encoding="utf-8") == "keep"
    assert not list(tmp_path.glob(f".{destination.name}-*"))


def test_linux_renameat2_uses_linux_at_fdcwd_and_ctypes_signature(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[tuple[object, ...]] = []

    class FakeRenameAt2:
        argtypes: object = None
        restype: object = None

        def __call__(self, *args: object) -> int:
            calls.append(args)
            return 0

    class FakeLibc:
        renameat2 = FakeRenameAt2()

    monkeypatch.setattr(balance_adoption.sys, "platform", "linux")
    monkeypatch.setattr(balance_adoption.ctypes, "CDLL", lambda *_args, **_kwargs: FakeLibc())

    balance_adoption._rename_directory_noreplace(tmp_path / ".staging", tmp_path / "published")

    renameat2 = FakeLibc.renameat2
    assert renameat2.argtypes == [
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    ]
    assert renameat2.restype is ctypes.c_int
    assert calls[0][0] == -100
    assert calls[0][2] == -100
    assert calls[0][4] == 1


@pytest.mark.parametrize("failure", ["missing-symbol", "enosys"])
def test_linux_unsupported_renameat2_is_explicit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    class FakeRenameAt2:
        argtypes: object = None
        restype: object = None

        def __call__(self, *_args: object) -> int:
            return -1

    class FakeLibc:
        if failure == "enosys":
            renameat2 = FakeRenameAt2()

    monkeypatch.setattr(balance_adoption.sys, "platform", "linux")
    monkeypatch.setattr(balance_adoption.ctypes, "CDLL", lambda *_args, **_kwargs: FakeLibc())
    if failure == "enosys":
        monkeypatch.setattr(balance_adoption.ctypes, "get_errno", lambda: errno.ENOSYS)

    staging = tmp_path / ".staging"
    _write_publication_tree(staging, "staged")
    with pytest.raises(ExperimentalInputError, match="unsupported"):
        balance_adoption._publish_staging_directory(staging, tmp_path / "published")

    assert staging.exists()


def test_staging_rejects_balance_arrays_altered_after_read(tmp_path: Path) -> None:
    source = read_case(write_balance_profiles(tmp_path / "balance"))
    original = source.profiles["density"]
    altered_values = original.values.copy()
    altered_values[0] += 1.0e12
    altered_profiles = dict(source.profiles)
    altered_profiles["density"] = ExperimentalProfile(
        name=original.name,
        units=original.units,
        path=original.path,
        coordinate=original.coordinate,
        values=altered_values,
    )
    altered_source = BalanceInput(
        metadata=source.metadata,
        profiles=MappingProxyType(altered_profiles),
        source_files=source.source_files,
        source_hashes=source.source_hashes,
    )

    with pytest.raises(ExperimentalInputError, match="(?i)(snapshot|array)"):
        stage_balance_marsf_quartet(
            altered_source,
            tmp_path / "marsf",
            major_radius_cm=_R0_CM,
            equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
        )


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
