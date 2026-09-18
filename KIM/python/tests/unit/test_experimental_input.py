from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from kim import (
    ExperimentalInputError,
    MarsFMetadata,
    read_marsf_profiles,
)
from pydantic import ValidationError

_FILES = {
    "density": "PROFDEN.IN",
    "electron_temperature": "PROFTE.IN",
    "ion_temperature": "PROFTI.IN",
    "toroidal_velocity": "PROFROT.IN",
}
FIXTURE = Path(__file__).parents[1] / "fixtures" / "experimental"


def metadata(**updates: object) -> MarsFMetadata:
    values: dict[str, object] = {
        "source": "mars-f-test-case",
        "coordinate": "sqrt_psiN",
        "coordinate_unit": "1",
        "density_unit": "1/m^3",
        "electron_temperature_unit": "eV",
        "ion_temperature_unit": "eV",
        "toroidal_velocity_unit": "m/s",
        "equilibrium_provenance": "approved-test-equilibrium",
    }
    values.update(updates)
    return MarsFMetadata(**values)


def write_case(directory: Path) -> dict[str, bytes]:
    directory.mkdir()
    coordinates = {
        "density": [0.0, 0.5, 1.0],
        "electron_temperature": [0.0, 0.4, 1.0],
        "ion_temperature": [0.0, 0.4, 1.0],
        "toroidal_velocity": [0.0, 0.25, 1.0],
    }
    values = {
        "density": [1.0e19, 2.0e19, 3.0e19],
        "electron_temperature": [100.0, 80.0, 40.0],
        "ion_temperature": [90.0, 70.0, 30.0],
        "toroidal_velocity": [0.0, -2.5e3, -5.0e3],
    }
    for name, filename in _FILES.items():
        rows = ["MARS-F profile: original source values"]
        rows.extend(
            f"{radius:.8e} {value:.8e}"
            for radius, value in zip(coordinates[name], values[name], strict=True)
        )
        (directory / filename).write_text("\n".join(rows) + "\n")
    return {filename: (directory / filename).read_bytes() for filename in _FILES.values()}


def test_reads_all_marsf_profiles_without_interpolation_or_conversion(tmp_path: Path) -> None:
    original = write_case(tmp_path / "marsf")

    result = read_marsf_profiles(tmp_path / "marsf", metadata())

    assert result.metadata.source == "mars-f-test-case"
    assert set(result.profiles) == set(_FILES)
    np.testing.assert_allclose(result.profiles["density"].coordinate, [0.0, 0.5, 1.0])
    np.testing.assert_allclose(result.profiles["density"].values, [1.0e19, 2.0e19, 3.0e19])
    np.testing.assert_allclose(result.profiles["toroidal_velocity"].coordinate, [0.0, 0.25, 1.0])
    assert result.profiles["toroidal_velocity"].units == "m/s"
    assert {path.name for path in result.source_files.values()} == set(_FILES.values())
    assert original == {
        filename: (tmp_path / "marsf" / filename).read_bytes() for filename in original
    }


def test_reads_the_checked_in_marsf_fixture() -> None:
    result = read_marsf_profiles(FIXTURE, metadata())

    assert result.profiles["density"].values[1] == pytest.approx(2.0e19)
    assert result.profiles["electron_temperature"].values[1] == pytest.approx(80.0)
    assert result.profiles["ion_temperature"].values[1] == pytest.approx(70.0)
    assert result.profiles["toroidal_velocity"].values[-1] == pytest.approx(-5.0e3)
    np.testing.assert_allclose(result.profiles["electron_temperature"].coordinate, [0.0, 0.4, 1.0])
    np.testing.assert_allclose(result.profiles["toroidal_velocity"].coordinate, [0.0, 0.25, 1.0])


@pytest.mark.parametrize(
    ("filename", "contents", "message"),
    [
        (
            "PROFTE.IN",
            "header\n0 1 2\n1 2\n",
            "exactly two numeric columns",
        ),
        (
            "PROFTI.IN",
            "header\n0 nan\n1 2\n",
            "row 2.*finite",
        ),
        (
            "PROFROT.IN",
            "header\n0 1\n0 2\n",
            "row 3.*strictly increasing",
        ),
        (
            "PROFDEN.IN",
            "header\n0 1\n",
            "at least two data rows",
        ),
        (
            "PROFDEN.IN",
            "0 1\n1 2\n",
            "header line",
        ),
    ],
)
def test_rejects_malformed_marsf_profile(
    tmp_path: Path, filename: str, contents: str, message: str
) -> None:
    write_case(tmp_path / "marsf")
    (tmp_path / "marsf" / filename).write_text(contents)

    with pytest.raises(ExperimentalInputError, match=message):
        read_marsf_profiles(tmp_path / "marsf", metadata())


def test_required_marsf_profile_must_exist(tmp_path: Path) -> None:
    write_case(tmp_path / "marsf")
    (tmp_path / "marsf" / "PROFTI.IN").unlink()

    with pytest.raises(ExperimentalInputError, match=r"PROFTI\.IN.*required"):
        read_marsf_profiles(tmp_path / "marsf", metadata())


def test_marsf_metadata_is_explicit_and_strict() -> None:
    with pytest.raises(ValidationError, match="density_unit"):
        metadata(density_unit="")
    with pytest.raises(ValidationError, match="coordinate_unit"):
        metadata(coordinate="r_eff", coordinate_unit="1")
    with pytest.raises(ValidationError, match="equilibrium_provenance"):
        metadata(equilibrium_provenance="")
    with pytest.raises(ValidationError, match="extra"):
        MarsFMetadata.model_validate(metadata().model_dump() | {"guess_units": True})


def test_marsf_reader_requires_a_directory(tmp_path: Path) -> None:
    with pytest.raises(ExperimentalInputError, match="directory"):
        read_marsf_profiles(tmp_path / "missing", metadata())
