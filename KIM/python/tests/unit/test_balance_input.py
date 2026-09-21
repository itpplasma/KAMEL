from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from kim import BalanceMetadata, ExperimentalInputError, read_balance_profiles
from pydantic import ValidationError

_PROFILE_ROLES = (
    "density",
    "electron_temperature",
    "ion_temperature",
    "toroidal_rotation",
)

_PROFILE_FILENAMES = {
    "density": "source-a.raw",
    "electron_temperature": "source-b.profile",
    "ion_temperature": "source-c.dat",
    "toroidal_rotation": "source-d.input",
}


def metadata(**updates: object) -> BalanceMetadata:
    values: dict[str, object] = {
        "source": "synthetic-balance-test-case",
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
        path = directory / _PROFILE_FILENAMES[role]
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


def read_case(paths: dict[str, Path], *, declared_metadata: BalanceMetadata | None = None):
    if declared_metadata is None:
        declared_metadata = metadata()
    return read_balance_profiles(
        density=paths["density"],
        electron_temperature=paths["electron_temperature"],
        ion_temperature=paths["ion_temperature"],
        toroidal_rotation=paths["toroidal_rotation"],
        metadata=declared_metadata,
    )


def test_reads_four_explicit_balance_profile_roles_and_preserves_grids(
    tmp_path: Path,
) -> None:
    paths = write_balance_profiles(tmp_path / "balance")
    declared_metadata = metadata()

    result = read_case(paths, declared_metadata=declared_metadata)

    assert set(result.profiles) == set(_PROFILE_ROLES)
    assert result.metadata == declared_metadata
    assert {role: result.source_files[role] for role in _PROFILE_ROLES} == {
        role: path.absolute() for role, path in paths.items()
    }
    np.testing.assert_allclose(result.profiles["density"].coordinate, [0.0, 0.5, 1.0])
    np.testing.assert_allclose(result.profiles["density"].values, [1.0e19, 2.0e19, 3.0e19])
    np.testing.assert_allclose(result.profiles["electron_temperature"].coordinate, [0.0, 0.4, 1.0])
    np.testing.assert_allclose(result.profiles["electron_temperature"].values, [100.0, 80.0, 40.0])
    np.testing.assert_allclose(result.profiles["ion_temperature"].coordinate, [0.0, 0.3, 1.0])
    np.testing.assert_allclose(result.profiles["ion_temperature"].values, [90.0, 70.0, 30.0])
    np.testing.assert_allclose(result.profiles["toroidal_rotation"].coordinate, [0.0, 0.25, 1.0])
    np.testing.assert_allclose(result.profiles["toroidal_rotation"].values, [0.0, -2.5e4, -5.0e4])
    assert result.profiles["density"].units == "1/m^3"
    assert result.profiles["electron_temperature"].units == "eV"
    assert result.profiles["ion_temperature"].units == "eV"
    assert result.profiles["toroidal_rotation"].units == "rad/s"


def test_balance_profiles_have_at_least_two_rows_and_finite_two_column_data(
    tmp_path: Path,
) -> None:
    paths = write_balance_profiles(tmp_path / "balance")

    result = read_case(paths)

    for profile in result.profiles.values():
        assert len(profile.coordinate) >= 2
        assert profile.coordinate.ndim == 1
        assert profile.values.ndim == 1
        assert np.isfinite(profile.coordinate).all()
        assert np.isfinite(profile.values).all()


@pytest.mark.parametrize(
    ("role", "contents", "message"),
    [
        ("density", "0.0 1.0e19\n", "at least two data rows"),
        ("density", "0.0\n0.5 2.0e19\n", "exactly two numeric columns"),
        ("electron_temperature", "0.0 1.0\n0.5 2.0 extra\n", "exactly two numeric columns"),
        ("ion_temperature", "0.0 1.0\n0.5 not-a-number\n", "exactly two numeric columns"),
        ("toroidal_rotation", "0.0 1.0\n0.5 nan\n", "finite"),
        ("density", "0.0 1.0\n0.5 inf\n", "finite"),
        ("density", "0.0 1.0e19\nnan 2.0e19\n", "finite"),
    ],
)
def test_rejects_missing_rows_non_numeric_or_non_finite_data(
    tmp_path: Path, role: str, contents: str, message: str
) -> None:
    paths = write_balance_profiles(tmp_path / "balance")
    paths[role].write_text(contents, encoding="utf-8")

    with pytest.raises(ExperimentalInputError, match=message):
        read_case(paths)


@pytest.mark.parametrize("role", _PROFILE_ROLES)
@pytest.mark.parametrize(
    "coordinates",
    [
        pytest.param(("0.0", "0.5", "0.5"), id="duplicate-coordinate"),
        pytest.param(("0.0", "0.75", "0.5"), id="non-monotonic-coordinate"),
    ],
)
def test_rejects_duplicate_or_non_monotonic_rho_pol(
    tmp_path: Path, role: str, coordinates: tuple[str, str, str]
) -> None:
    paths = write_balance_profiles(tmp_path / "balance")
    paths[role].write_text(
        "\n".join(f"{coordinate} {value}" for coordinate, value in zip(coordinates, [1, 2, 3]))
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ExperimentalInputError, match="strictly increasing"):
        read_case(paths)


def test_requires_each_profile_path_to_exist(tmp_path: Path) -> None:
    paths = write_balance_profiles(tmp_path / "balance")
    paths["ion_temperature"].unlink()

    with pytest.raises(ExperimentalInputError, match=r"source-c\.dat.*missing"):
        read_case(paths)


def test_requires_explicit_metadata(tmp_path: Path) -> None:
    paths = write_balance_profiles(tmp_path / "balance")

    with pytest.raises(TypeError, match="metadata"):
        read_balance_profiles(
            density=paths["density"],
            electron_temperature=paths["electron_temperature"],
            ion_temperature=paths["ion_temperature"],
            toroidal_rotation=paths["toroidal_rotation"],
        )


def test_reader_does_not_mutate_source_file_bytes(tmp_path: Path) -> None:
    paths = write_balance_profiles(tmp_path / "balance")
    original = {role: path.read_bytes() for role, path in paths.items()}

    read_case(paths)

    assert {role: path.read_bytes() for role, path in paths.items()} == original


def test_balance_metadata_is_explicit_frozen_and_does_not_infer_units() -> None:
    declared = metadata()
    assert declared.coordinate == "rho_pol"
    assert declared.coordinate_unit == "1"
    assert declared.density_unit == "1/m^3"
    assert declared.electron_temperature_unit == "eV"
    assert declared.ion_temperature_unit == "eV"
    assert declared.toroidal_rotation_unit == "rad/s"

    with pytest.raises(ValidationError):
        declared.density_unit = "1/cm^3"

    values = declared.model_dump()
    values.pop("density_unit")
    with pytest.raises(ValidationError, match="density_unit"):
        BalanceMetadata(**values)

    with pytest.raises(ValidationError, match="extra"):
        BalanceMetadata(**declared.model_dump(), guess_units=True)


def test_accepts_alternate_density_unit_spelling_and_preserves_declaration(
    tmp_path: Path,
) -> None:
    paths = write_balance_profiles(tmp_path / "balance")
    declared = metadata(density_unit="m^-3")

    result = read_case(paths, declared_metadata=declared)

    assert result.metadata.density_unit == "m^-3"
    assert result.profiles["density"].units == "m^-3"
