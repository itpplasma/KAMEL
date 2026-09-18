"""Strict, read-only readers for approved experimental input formats."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Literal, Mapping

import numpy as np
from kim.errors import ExperimentalInputError
from numpy.typing import NDArray
from pydantic import BaseModel, ConfigDict, Field, model_validator


class MarsFMetadata(BaseModel):
    """Conventions declared by one MARS-F profile source.

    The reader records these declarations but does not convert them.  A later
    preparation step is responsible for applying approved transformations.
    """

    model_config = ConfigDict(
        extra="forbid",
        frozen=True,
        str_strip_whitespace=True,
        validate_default=True,
    )

    schema_version: Literal[1] = 1
    source: str = Field(min_length=1)
    coordinate: Literal["sqrt_psiN", "r_eff"]
    coordinate_unit: Literal["1", "cm", "m"]
    density_unit: Literal["1/m^3", "m^-3", "1/cm^3", "cm^-3"]
    electron_temperature_unit: Literal["eV", "keV", "K"]
    ion_temperature_unit: Literal["eV", "keV", "K"]
    toroidal_velocity_unit: Literal["m/s", "cm/s"]
    equilibrium_provenance: str = Field(min_length=1)

    @model_validator(mode="after")
    def coordinate_unit_is_compatible(self) -> MarsFMetadata:
        if self.coordinate == "sqrt_psiN" and self.coordinate_unit != "1":
            raise ValueError("coordinate_unit must be 1 for a sqrt_psiN coordinate")
        if self.coordinate == "r_eff" and self.coordinate_unit not in {"cm", "m"}:
            raise ValueError("coordinate_unit must be cm or m for an r_eff coordinate")
        return self


@dataclass(frozen=True)
class ExperimentalProfile:
    """One untouched two-column profile from an experimental source."""

    name: str
    units: str
    path: Path
    coordinate: NDArray[np.float64]
    values: NDArray[np.float64]


@dataclass(frozen=True)
class MarsFInput:
    """The four source profiles and their explicit source declarations."""

    directory: Path
    metadata: MarsFMetadata
    profiles: Mapping[str, ExperimentalProfile]
    source_files: Mapping[str, Path]


_PROFILE_SPECS: tuple[tuple[str, str, str], ...] = (
    ("density", "PROFDEN.IN", "density_unit"),
    ("electron_temperature", "PROFTE.IN", "electron_temperature_unit"),
    ("ion_temperature", "PROFTI.IN", "ion_temperature_unit"),
    ("toroidal_velocity", "PROFROT.IN", "toroidal_velocity_unit"),
)


def read_marsf_profiles(directory: Path | str, metadata: MarsFMetadata) -> MarsFInput:
    """Read the approved MARS-F profile quartet without scientific reduction.

    Each source file must contain one header line followed by at least two
    rows with exactly two finite numeric values.  Every profile keeps its own
    coordinate grid because interpolation and coordinate reduction belong to
    the preparation stage, not this reader.
    """

    if not isinstance(metadata, MarsFMetadata):
        raise ExperimentalInputError("metadata must be a MarsFMetadata instance")

    root = Path(directory)
    if not root.is_dir():
        raise ExperimentalInputError(f"MARS-F input directory does not exist: {root}")

    profiles: dict[str, ExperimentalProfile] = {}
    source_files: dict[str, Path] = {}
    for name, filename, unit_field in _PROFILE_SPECS:
        path = root / filename
        if not path.is_file():
            raise ExperimentalInputError(f"{path}: required MARS-F profile is missing")
        profile = _read_profile(name, getattr(metadata, unit_field), path)
        profiles[name] = profile
        source_files[name] = path.absolute()

    return MarsFInput(
        directory=root.absolute(),
        metadata=metadata,
        profiles=MappingProxyType(profiles),
        source_files=MappingProxyType(source_files),
    )


def _read_profile(name: str, units: str, path: Path) -> ExperimentalProfile:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise ExperimentalInputError(f"{path}: unable to read MARS-F profile") from error

    if not lines or not lines[0].strip():
        raise ExperimentalInputError(f"{path}: first line must be a MARS-F header line")
    # The legacy reader skipped this line without interpreting its grammar;
    # A06 does not approve a stronger header contract.
    if _is_numeric_pair(lines[0]):
        raise ExperimentalInputError(f"{path}: first line must be a MARS-F header line")

    rows: list[tuple[float, float]] = []
    for line_number, line in enumerate(lines[1:], start=2):
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "!")):
            continue
        columns = stripped.split()
        if len(columns) != 2:
            raise ExperimentalInputError(
                f"{path}: row {line_number} must contain exactly two numeric columns"
            )
        try:
            coordinate, value = (_parse_float(column) for column in columns)
        except ValueError as error:
            raise ExperimentalInputError(
                f"{path}: row {line_number} must contain exactly two numeric columns"
            ) from error
        if not np.isfinite(coordinate) or not np.isfinite(value):
            raise ExperimentalInputError(f"{path}: row {line_number} {name} values must be finite")
        if rows and coordinate <= rows[-1][0]:
            raise ExperimentalInputError(
                f"{path}: row {line_number} coordinate must be strictly increasing"
            )
        rows.append((coordinate, value))

    if len(rows) < 2:
        raise ExperimentalInputError(f"{path}: profile must contain at least two data rows")

    data = np.asarray(rows, dtype=np.float64)
    return ExperimentalProfile(
        name=name,
        units=units,
        path=path.absolute(),
        coordinate=data[:, 0],
        values=data[:, 1],
    )


def _is_numeric_pair(line: str) -> bool:
    columns = line.split()
    if len(columns) != 2:
        return False
    try:
        _parse_float(columns[0])
        _parse_float(columns[1])
    except ValueError:
        return False
    return True


def _parse_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


__all__ = [
    "ExperimentalProfile",
    "MarsFInput",
    "MarsFMetadata",
    "read_marsf_profiles",
]
