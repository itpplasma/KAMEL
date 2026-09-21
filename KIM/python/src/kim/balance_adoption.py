"""Pure, explicit transformations for BALANCE profile adoption.

This module is intentionally narrower than :mod:`kim.conventions`.  It accepts
only the conventions used by the approved BALANCE source files and returns
detached arrays together with the operation needed to explain each conversion.
It does not read or write files, map radial coordinates, or infer scientific
metadata from values.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
from kim.conventions import ConventionError, ConversionOperation
from kim.errors import ExperimentalInputError
from kim.importers.balance import BalanceInput
from kim.importers.experimental import (
    ExperimentalProfile,
    MarsFMetadata,
    read_marsf_profiles,
)
from numpy.typing import ArrayLike, NDArray

_BALANCE_COORDINATES = {"rho_pol", "sqrt_psiN"}
_BALANCE_COORDINATE_UNIT = "1"
_BALANCE_PROFILE_ROLES = (
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


@dataclass(frozen=True)
class BalanceTransformation:
    """Detached BALANCE arrays and their explicit conversion operation."""

    coordinate: NDArray[np.float64]
    values: NDArray[np.float64]
    operation: ConversionOperation


@dataclass(frozen=True)
class StagedMarsFQuartet:
    """An atomically staged, reader-compatible MARS-F profile quartet."""

    directory: Path
    metadata: MarsFMetadata
    report: Path


def stage_balance_marsf_quartet(
    source: BalanceInput,
    destination: Path | str,
    *,
    major_radius_cm: float,
    equilibrium_provenance: str,
) -> StagedMarsFQuartet:
    """Stage BALANCE profiles as a derived, reader-compatible MARS-F quartet.

    The source arrays stay in their declared units except for the approved
    angular-rotation-to-linear-velocity conversion.  Staging happens in a
    temporary sibling directory and is committed with one directory rename;
    the existing preparation path remains responsible for later CGS density
    conversion and equilibrium mapping.
    """

    staging_directory: Path | None = None
    try:
        source_hashes = _validate_staging_arguments(source, major_radius_cm, equilibrium_provenance)
        final_directory = Path(destination).absolute()
        if final_directory.exists() or final_directory.is_symlink():
            raise ExperimentalInputError(
                f"staged quartet destination already exists: {final_directory}"
            )
        final_directory.parent.mkdir(parents=True, exist_ok=True)

        staging_directory = Path(
            tempfile.mkdtemp(prefix=f".{final_directory.name}-", dir=final_directory.parent)
        )
        metadata = MarsFMetadata(
            source=source.metadata.source,
            coordinate="sqrt_psiN",
            coordinate_unit="1",
            density_unit=source.metadata.density_unit,
            electron_temperature_unit=source.metadata.electron_temperature_unit,
            ion_temperature_unit=source.metadata.ion_temperature_unit,
            toroidal_velocity_unit="cm/s",
            equilibrium_provenance=equilibrium_provenance,
        )

        operations: list[dict[str, object]] = [
            _preservation_operation("density", source.metadata.density_unit),
            _preservation_operation("electron_temperature", "eV"),
            _preservation_operation("ion_temperature", "eV"),
        ]
        rotation = source.profiles["toroidal_rotation"]
        rotation_result = convert_balance_rotation(
            rotation.coordinate,
            rotation.values,
            coordinate=source.metadata.coordinate,
            coordinate_unit=source.metadata.coordinate_unit,
            source_unit=source.metadata.toroidal_rotation_unit,
            r0_cm=major_radius_cm,
        )
        operations.append(rotation_result.operation.model_dump(mode="json"))

        profiles = {
            "density": source.profiles["density"],
            "electron_temperature": source.profiles["electron_temperature"],
            "ion_temperature": source.profiles["ion_temperature"],
            "toroidal_rotation": ExperimentalProfile(
                name="toroidal_velocity",
                units="cm/s",
                path=rotation.path,
                coordinate=rotation_result.coordinate,
                values=rotation_result.values,
            ),
        }
        for role in _BALANCE_PROFILE_ROLES:
            profile = profiles[role]
            _write_profile(
                staging_directory / _MARSF_FILENAMES[role],
                profile.coordinate,
                profile.values,
            )

        # Keep the adapter coupled to the existing strict reader contract.
        read_marsf_profiles(staging_directory, metadata)
        derived_hashes = {
            filename: _sha256(staging_directory / filename)
            for filename in _MARSF_FILENAMES.values()
        }
        report_payload = {
            "schema_version": 1,
            "source_basenames": {
                role: Path(source.source_files[role]).name for role in _BALANCE_PROFILE_ROLES
            },
            "source_hashes": {role: source_hashes[role] for role in _BALANCE_PROFILE_ROLES},
            "derived_hashes": derived_hashes,
            "major_radius_cm": float(major_radius_cm),
            "coordinate_mapping": {
                "source": source.metadata.coordinate,
                "source_unit": source.metadata.coordinate_unit,
                "target": "sqrt_psiN",
                "target_unit": "1",
                "operation": (
                    "rho_pol = sqrt(psi_pol_norm)"
                    if source.metadata.coordinate == "rho_pol"
                    else "preserve"
                ),
            },
            "equilibrium_provenance": equilibrium_provenance,
            "operations": operations,
        }
        report_name = "staging_report.json"
        (staging_directory / report_name).write_text(
            json.dumps(report_payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )

        _commit_staging_directory(staging_directory, final_directory)
        return StagedMarsFQuartet(
            directory=final_directory,
            metadata=metadata,
            report=final_directory / report_name,
        )
    except (MemoryError, KeyboardInterrupt):
        raise
    except ExperimentalInputError:
        raise
    except Exception as error:
        raise ExperimentalInputError("unable to stage BALANCE MARS-F quartet") from error
    finally:
        if staging_directory is not None and staging_directory.exists():
            shutil.rmtree(staging_directory, ignore_errors=True)


def _validate_staging_arguments(
    source: BalanceInput, major_radius_cm: float, equilibrium_provenance: str
) -> dict[str, str]:
    if not isinstance(source, BalanceInput):
        raise ExperimentalInputError("source must be a BalanceInput instance")
    try:
        radius = float(major_radius_cm)
    except (OverflowError, TypeError, ValueError) as error:
        raise ExperimentalInputError("major_radius_cm must be finite and positive") from error
    if not np.isfinite(radius) or radius <= 0.0:
        raise ExperimentalInputError("major_radius_cm must be finite and positive")
    if not isinstance(equilibrium_provenance, str) or not equilibrium_provenance.strip():
        raise ExperimentalInputError("equilibrium_provenance must be nonempty")
    if source.metadata.coordinate_unit != _BALANCE_COORDINATE_UNIT:
        raise ExperimentalInputError("BALANCE coordinate_unit must be 1")
    source_hashes: dict[str, str] = {}
    for role in _BALANCE_PROFILE_ROLES:
        if role not in source.profiles:
            raise ExperimentalInputError(f"BALANCE source is missing {role} profile")
        if role not in source.source_files:
            raise ExperimentalInputError(f"BALANCE source is missing {role} source file")
        _validate_source_profile(source.profiles[role], role)
        path = Path(source.source_files[role])
        if not path.is_file():
            raise ExperimentalInputError(f"{path}: BALANCE profile source is missing")
        try:
            current_hash = _sha256(path)
        except OSError as error:
            raise ExperimentalInputError(f"{path}: unable to verify BALANCE source hash") from error
        if source.source_hashes.get(role) != current_hash:
            raise ExperimentalInputError(f"{path}: BALANCE source hash changed since it was read")
        source_hashes[role] = current_hash
    return source_hashes


def _reserve_destination(destination: Path) -> None:
    """Reserve a previously absent destination without replacing a race winner."""

    try:
        destination.mkdir()
    except FileExistsError as error:
        raise ExperimentalInputError(
            f"staged quartet destination already exists: {destination}"
        ) from error


def _commit_staging_directory(staging: Path, destination: Path) -> None:
    """Reserve the destination, move staged children, and clean up on failure."""

    _reserve_destination(destination)
    try:
        for child in tuple(staging.iterdir()):
            target = destination / child.name
            if target.exists() or target.is_symlink():
                raise ExperimentalInputError(f"staged quartet target already exists: {target}")
            child.rename(target)
        staging.rmdir()
    except Exception:
        shutil.rmtree(destination, ignore_errors=True)
        raise


def _validate_source_profile(profile: ExperimentalProfile, role: str) -> None:
    coordinates = np.asarray(profile.coordinate)
    values = np.asarray(profile.values)
    if coordinates.ndim != 1 or values.ndim != 1 or coordinates.size != values.size:
        raise ExperimentalInputError(
            f"{role} profile coordinate and values must be one-dimensional"
        )
    if coordinates.size < 2:
        raise ExperimentalInputError(f"{role} profile must contain at least two rows")
    if not np.issubdtype(coordinates.dtype, np.number) or not np.issubdtype(
        values.dtype, np.number
    ):
        raise ExperimentalInputError(f"{role} profile values must be numeric")
    if not np.all(np.isfinite(coordinates)) or not np.all(np.isfinite(values)):
        raise ExperimentalInputError(f"{role} profile values must be finite")
    if np.any(np.diff(coordinates) <= 0.0):
        raise ExperimentalInputError(f"{role} profile coordinate must be strictly increasing")


def _preservation_operation(quantity: str, unit: str) -> dict[str, object]:
    return {
        "quantity": quantity,
        "source_unit": unit,
        "target_unit": unit,
        "factor": 1.0,
        "operation": "preserve",
    }


def _write_profile(
    path: Path, coordinate: NDArray[np.float64], values: NDArray[np.float64]
) -> None:
    """Write one reader-compatible profile; kept as the atomicity test seam."""

    data = np.column_stack((coordinate, values))
    np.savetxt(
        path,
        data,
        fmt="%.17e",
        header="BALANCE-derived MARS-F profile",
        comments="",
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def convert_balance_density(
    coordinate_values: ArrayLike,
    values: ArrayLike,
    *,
    coordinate: str,
    coordinate_unit: str,
    source_unit: str,
) -> BalanceTransformation:
    """Convert an explicitly declared BALANCE density to ``1/cm^3``."""

    coordinate_array, value_array = _validate_profile(
        coordinate_values,
        values,
        coordinate=coordinate,
        coordinate_unit=coordinate_unit,
    )
    if source_unit in {"1/m^3", "m^-3"}:
        factor = 1.0e-6
        operation = "density_si_to_cgs"
    elif source_unit in {"1/cm^3", "cm^-3"}:
        factor = 1.0
        operation = "preserve"
    else:
        raise ConventionError(
            f"density unit {source_unit!r} is unsupported; declare SI or CGS density explicitly"
        )

    converted = _scaled_values(value_array, factor, "density")
    return BalanceTransformation(
        coordinate=coordinate_array,
        values=converted,
        operation=ConversionOperation(
            quantity="density",
            source_unit=source_unit,
            target_unit="1/cm^3",
            factor=factor,
            operation=operation,
        ),
    )


def convert_balance_temperature(
    coordinate_values: ArrayLike,
    values: ArrayLike,
    *,
    species: Literal["electron", "ion"] | str,
    coordinate: str,
    coordinate_unit: str,
    source_unit: str,
) -> BalanceTransformation:
    """Validate an explicitly declared BALANCE temperature in eV."""

    coordinate_array, value_array = _validate_profile(
        coordinate_values,
        values,
        coordinate=coordinate,
        coordinate_unit=coordinate_unit,
    )
    if species not in {"electron", "ion"}:
        raise ConventionError("temperature species must be explicitly declared as electron or ion")
    if source_unit != "eV":
        raise ConventionError("temperature unit must be declared in eV; conversion is not inferred")

    return BalanceTransformation(
        coordinate=coordinate_array,
        values=_readonly_copy(value_array),
        operation=ConversionOperation(
            quantity=f"{species}_temperature",
            source_unit="eV",
            target_unit="eV",
            factor=1.0,
            operation="preserve",
        ),
    )


def convert_balance_rotation(
    coordinate_values: ArrayLike,
    values: ArrayLike,
    *,
    coordinate: str,
    coordinate_unit: str,
    source_unit: str,
    r0_cm: float | None,
) -> BalanceTransformation:
    """Convert signed angular rotation to signed toroidal velocity in cm/s."""

    coordinate_array, value_array = _validate_profile(
        coordinate_values,
        values,
        coordinate=coordinate,
        coordinate_unit=coordinate_unit,
    )
    if source_unit != "rad/s":
        raise ConventionError(
            "rotation unit must be declared in rad/s; velocity conversion is not inferred"
        )
    try:
        radius = float(r0_cm)
    except (OverflowError, TypeError, ValueError) as error:
        raise ConventionError("r0_cm must be a finite positive radius") from error
    if not np.isfinite(radius) or radius <= 0.0:
        raise ConventionError("r0_cm must be a finite positive radius")

    converted = _scaled_values(value_array, radius, "rotation")
    return BalanceTransformation(
        coordinate=coordinate_array,
        values=converted,
        operation=ConversionOperation(
            quantity="toroidal_velocity",
            source_unit="rad/s",
            target_unit="cm/s",
            factor=radius,
            operation="omega_to_v_phi",
            parameters={"major_radius_cm": radius},
        ),
    )


def convert_balance_q(
    coordinate_values: ArrayLike,
    values: ArrayLike,
    *,
    coordinate: str,
    coordinate_unit: str,
    source_unit: str,
    operation: Literal["preserve", "negate"] | str,
) -> BalanceTransformation:
    """Apply an explicitly selected signed-q operation."""

    coordinate_array, value_array = _validate_profile(
        coordinate_values,
        values,
        coordinate=coordinate,
        coordinate_unit=coordinate_unit,
    )
    if source_unit != "1":
        raise ConventionError("q unit must be declared as dimensionless unit 1")
    if operation not in {"preserve", "negate"}:
        raise ConventionError("q operation must be explicitly selected as preserve or negate")

    factor = 1.0 if operation == "preserve" else -1.0
    converted = _scaled_values(value_array, factor, "q")
    return BalanceTransformation(
        coordinate=coordinate_array,
        values=converted,
        operation=ConversionOperation(
            quantity="q",
            source_unit="1",
            target_unit="1",
            factor=factor,
            operation=operation,
        ),
    )


def _validate_profile(
    coordinate_values: ArrayLike,
    values: ArrayLike,
    *,
    coordinate: str,
    coordinate_unit: str,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    if coordinate not in _BALANCE_COORDINATES:
        raise ConventionError("coordinate must be explicitly declared as rho_pol or sqrt_psiN")
    if coordinate_unit != _BALANCE_COORDINATE_UNIT:
        raise ConventionError("coordinate_unit must be 1 for rho_pol or sqrt_psiN")

    # Strict BALANCE readers enforce strictly increasing coordinates upstream.
    coordinate_array = _as_finite_vector(coordinate_values, "coordinate")
    value_array = _as_finite_vector(values, "profile")
    if coordinate_array.size != value_array.size:
        raise ConventionError("coordinate and values must have equal lengths")
    return coordinate_array, value_array


def _as_finite_vector(values: ArrayLike, name: str) -> NDArray[np.float64]:
    try:
        raw = np.asarray(values)
    except (TypeError, ValueError) as error:
        raise ConventionError(f"{name} values must be numeric and finite") from error
    if np.iscomplexobj(raw) or (
        raw.dtype.kind == "O" and any(np.iscomplexobj(item) for item in raw.flat)
    ):
        raise ConventionError(f"{name} values must be real; complex values are unsupported")
    try:
        array = np.asarray(raw, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise ConventionError(f"{name} values must be numeric and finite") from error
    if array.ndim != 1:
        raise ConventionError(f"{name} values must be one-dimensional")
    if not np.all(np.isfinite(array)):
        raise ConventionError(f"{name} values must be finite")
    result = np.array(array, dtype=np.float64, copy=True)
    result.setflags(write=False)
    return result


def _scaled_values(values: NDArray[np.float64], factor: float, name: str) -> NDArray[np.float64]:
    converted = values * factor
    if not np.all(np.isfinite(converted)):
        raise ConventionError(f"converted {name} values must be finite")
    result = np.array(converted, dtype=np.float64, copy=True)
    result.setflags(write=False)
    return result


def _readonly_copy(values: NDArray[np.float64]) -> NDArray[np.float64]:
    result = np.array(values, dtype=np.float64, copy=True)
    result.setflags(write=False)
    return result


__all__ = [
    "BalanceTransformation",
    "StagedMarsFQuartet",
    "convert_balance_density",
    "convert_balance_q",
    "convert_balance_rotation",
    "convert_balance_temperature",
    "stage_balance_marsf_quartet",
]
