"""Strict, read-only readers for BALANCE profile sources."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Literal, Mapping

import numpy as np
from kim.errors import ExperimentalInputError
from kim.importers.experimental import ExperimentalProfile
from numpy.typing import NDArray
from pydantic import BaseModel, ConfigDict, Field


class BalanceMetadata(BaseModel):
    """Explicit conventions declared by one BALANCE profile source.

    The reader records these declarations but deliberately does not infer or
    convert any scientific quantity.  ``rho_pol`` is the BALANCE spelling of
    the dimensionless ``sqrt_psiN`` coordinate used by KIM preparation.
    """

    model_config = ConfigDict(
        extra="forbid",
        frozen=True,
        str_strip_whitespace=True,
        validate_default=True,
    )

    schema_version: Literal[1] = 1
    source: str = Field(min_length=1)
    coordinate: Literal["rho_pol", "sqrt_psiN"]
    coordinate_unit: Literal["1"]
    density_unit: Literal["1/m^3", "m^-3", "1/cm^3", "cm^-3"]
    electron_temperature_unit: Literal["eV"]
    ion_temperature_unit: Literal["eV"]
    toroidal_rotation_unit: Literal["rad/s"]


@dataclass(frozen=True)
class BalanceInput:
    """The four source profiles and their explicit source declarations."""

    metadata: BalanceMetadata
    profiles: Mapping[str, ExperimentalProfile]
    source_files: Mapping[str, Path]
    source_hashes: Mapping[str, str]


_PROFILE_SPECS: tuple[tuple[str, str], ...] = (
    ("density", "density_unit"),
    ("electron_temperature", "electron_temperature_unit"),
    ("ion_temperature", "ion_temperature_unit"),
    ("toroidal_rotation", "toroidal_rotation_unit"),
)


def read_balance_profiles(
    *,
    density: Path | str,
    electron_temperature: Path | str,
    ion_temperature: Path | str,
    toroidal_rotation: Path | str,
    metadata: BalanceMetadata,
) -> BalanceInput:
    """Read four explicitly supplied, headerless BALANCE profile files.

    Each file must contain at least two rows of exactly two finite numeric
    columns.  Profiles retain independent source grids; interpolation,
    conversion, and coordinate mapping belong to later preparation steps.
    """

    if not isinstance(metadata, BalanceMetadata):
        raise ExperimentalInputError("metadata must be a BalanceMetadata instance")

    candidates = {
        "density": Path(density),
        "electron_temperature": Path(electron_temperature),
        "ion_temperature": Path(ion_temperature),
        "toroidal_rotation": Path(toroidal_rotation),
    }
    absolute_paths = {name: _resolve_source_path(path) for name, path in candidates.items()}
    identities = {name: _file_identity(path) for name, path in absolute_paths.items()}
    duplicate_identities = {
        identity
        for identity in identities.values()
        if list(identities.values()).count(identity) > 1
    }
    if duplicate_identities:
        duplicate_names = ", ".join(
            str(absolute_paths[name])
            for name, identity in identities.items()
            if identity in duplicate_identities
        )
        raise ExperimentalInputError(f"BALANCE profile paths must be distinct: {duplicate_names}")

    profiles: dict[str, ExperimentalProfile] = {}
    source_hashes: dict[str, str] = {}
    for name, unit_field in _PROFILE_SPECS:
        path = absolute_paths[name]
        units = getattr(metadata, unit_field)
        profile, source_hash = _read_profile(name, units, path)
        profiles[name] = profile
        source_hashes[name] = source_hash

    return BalanceInput(
        metadata=metadata,
        profiles=MappingProxyType(profiles),
        source_files=MappingProxyType(absolute_paths),
        source_hashes=MappingProxyType(source_hashes),
    )


def _resolve_source_path(candidate: Path) -> Path:
    try:
        path = candidate.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise ExperimentalInputError(
            f"{candidate}: BALANCE profile source is missing or unavailable"
        ) from error
    try:
        is_file = path.is_file()
    except OSError as error:
        raise ExperimentalInputError(f"{path}: unable to inspect BALANCE profile") from error
    if not is_file:
        raise ExperimentalInputError(f"{path}: BALANCE profile source is missing")
    return path


def _file_identity(path: Path) -> tuple[int, int]:
    try:
        file_stat = path.stat()
    except OSError as error:
        raise ExperimentalInputError(f"{path}: unable to inspect BALANCE profile") from error
    return file_stat.st_dev, file_stat.st_ino


def _read_profile(name: str, units: str, path: Path) -> tuple[ExperimentalProfile, str]:
    return _read_profile_snapshot(name, units, path)


def _read_profile_snapshot(name: str, units: str, path: Path) -> tuple[ExperimentalProfile, str]:
    """Read, hash, and parse one profile from the same immutable byte snapshot."""

    try:
        source_bytes = path.read_bytes()
    except (OSError, UnicodeError) as error:
        raise ExperimentalInputError(f"{path}: unable to read BALANCE profile") from error
    try:
        source_hash = _sha256(source_bytes)
    except OSError as error:
        raise ExperimentalInputError(f"{path}: unable to hash BALANCE profile") from error
    try:
        lines = source_bytes.decode("utf-8").splitlines()
    except UnicodeError as error:
        raise ExperimentalInputError(f"{path}: unable to read BALANCE profile") from error
    return _parse_profile_lines(name, units, path, lines), source_hash


def _parse_profile_lines(
    name: str, units: str, path: Path, lines: list[str]
) -> ExperimentalProfile:
    rows: list[tuple[float, float]] = []
    for line_number, line in enumerate(lines, start=1):
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
    data.setflags(write=False)
    coordinate = data[:, 0]
    values = data[:, 1]
    return ExperimentalProfile(
        name=name,
        units=units,
        path=path,
        coordinate=coordinate,
        values=values,
    )


def _parse_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def _sha256(source_bytes: bytes) -> str:
    return hashlib.sha256(source_bytes).hexdigest()


__all__ = ["BalanceInput", "BalanceMetadata", "read_balance_profiles"]
