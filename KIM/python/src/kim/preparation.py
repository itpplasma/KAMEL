"""Preparation of explicit experimental profiles for a KIM case."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal, Mapping, Sequence

import numpy as np
from kim import executable as executable_module
from kim.config import ProfileConfig, SimulationConfig
from kim.conventions import ConversionOperation, SourceMetadata, convert_quantity
from kim.errors import ExperimentalInputError
from kim.importers import experimental as experimental_module
from kim.importers.balance import BalanceMetadata
from kim.importers.experimental import ExperimentalProfile, MarsFInput, MarsFProfileSnapshot
from numpy.typing import NDArray

QOperation = Literal["preserve", "negate"]

_BALANCE_PROFILE_ROLES = (
    "density",
    "electron_temperature",
    "ion_temperature",
    "toroidal_rotation",
)
_MARSF_PROFILE_FILENAMES = {
    "density": "PROFDEN.IN",
    "electron_temperature": "PROFTE.IN",
    "ion_temperature": "PROFTI.IN",
    "toroidal_rotation": "PROFROT.IN",
}
_MARSF_PROFILE_UNIT_FIELDS = {
    "density": "density_unit",
    "electron_temperature": "electron_temperature_unit",
    "ion_temperature": "ion_temperature_unit",
    "toroidal_velocity": "toroidal_velocity_unit",
}
_STAGING_REPORT_V1_KEYS = frozenset(
    {
        "schema_version",
        "source_basenames",
        "source_hashes",
        "source_metadata",
        "derived_hashes",
        "major_radius_cm",
        "coordinate_mapping",
        "equilibrium_provenance",
        "operations",
    }
)
_STAGING_REPORT_V2_KEYS = _STAGING_REPORT_V1_KEYS | {"equilibrium_parameters"}
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"JSON constant {value!r} is not allowed")


@dataclass(frozen=True)
class PreparedExperimentalCase:
    """Files and validated request produced by experimental preparation."""

    directory: Path
    profiles: Path
    equilibrium: Path
    request: Path
    report: Path
    config: SimulationConfig
    equilibrium_parameters: Path | None = None


@dataclass(frozen=True)
class EquilibriumParameters:
    """Magnetic-axis scalars emitted by the equilibrium preprocessor."""

    btor_gauss: float
    r_big_cm: float
    sha256: str


def read_equilibrium_parameters(path: Path | str) -> EquilibriumParameters:
    """Read the exact ``btor_rbig.dat`` output record from an equilibrium calculation."""

    candidate = Path(path)
    if candidate.is_symlink() or not candidate.is_file():
        raise ExperimentalInputError(
            f"{candidate}: btor_rbig.dat must be a regular equilibrium output file"
        )
    raw_bytes, digest = _read_bytes_snapshot(candidate, "btor_rbig.dat equilibrium output")
    try:
        columns = raw_bytes.decode("utf-8").split()
    except UnicodeError as error:
        raise ExperimentalInputError(f"{candidate}: btor_rbig.dat is not valid text") from error
    if len(columns) != 2:
        raise ExperimentalInputError(
            f"{candidate}: btor_rbig.dat must contain exactly btor and r_big"
        )
    try:
        btor_gauss, r_big_cm = (
            float(column.replace("D", "E").replace("d", "e")) for column in columns
        )
    except (OverflowError, ValueError) as error:
        raise ExperimentalInputError(f"{candidate}: btor_rbig.dat values are invalid") from error
    if not np.isfinite(btor_gauss) or btor_gauss == 0.0:
        raise ExperimentalInputError(f"{candidate}: btor_rbig.dat btor must be finite and nonzero")
    if not np.isfinite(r_big_cm) or r_big_cm <= 0.0:
        raise ExperimentalInputError(
            f"{candidate}: btor_rbig.dat r_big must be finite and positive"
        )
    return EquilibriumParameters(btor_gauss, r_big_cm, digest)


@dataclass(frozen=True)
class EquilibriumCalculation:
    """Paired profile and magnetic-axis outputs from one equilibrium run."""

    equilibrium_file: Path
    parameters_file: Path
    parameters: EquilibriumParameters
    equilibrium_sha256: str
    generator_provenance: Mapping[str, object]


def simulation_config_for_equilibrium(
    config: SimulationConfig, parameters: EquilibriumParameters
) -> SimulationConfig:
    """Return a validated KIM config using equilibrium-produced ``btor`` and ``r_big``."""

    if not isinstance(config, SimulationConfig):
        raise ExperimentalInputError("config must be a SimulationConfig instance")
    if not isinstance(parameters, EquilibriumParameters):
        raise ExperimentalInputError("parameters must be an EquilibriumParameters instance")
    payload = config.model_dump(mode="python")
    payload["setup"]["btor"] = parameters.btor_gauss
    payload["setup"]["major_radius"] = parameters.r_big_cm
    try:
        return SimulationConfig.model_validate(payload)
    except Exception as error:
        raise ExperimentalInputError(
            "equilibrium parameters do not form a valid KIM config"
        ) from error


def run_equilibrium_calculation(
    equilibrium_executable: Path | str,
    equilibrium_input_files: Sequence[Path | str],
    destination: Path | str,
    *,
    timeout_seconds: float = 3600.0,
) -> EquilibriumCalculation:
    """Run the supplied equilibrium preprocessor once and retain its paired outputs."""

    if not np.isfinite(timeout_seconds) or timeout_seconds <= 0.0:
        raise ExperimentalInputError("equilibrium timeout must be positive and finite")
    command = _resolve_generator_command(equilibrium_executable)
    generator_path = _resolved_executable_path(command[0])
    if generator_path is None or not os.access(generator_path, os.X_OK):
        raise ExperimentalInputError(
            f"equilibrium_executable does not name an executable file: {command[0]}"
        )
    _validate_equilibrium_input_names(equilibrium_input_files, generator_path)

    work_directory = Path(destination).absolute()
    if work_directory.exists() or work_directory.is_symlink():
        raise ExperimentalInputError(
            f"equilibrium calculation destination already exists: {work_directory}"
        )
    if not work_directory.parent.is_dir():
        raise ExperimentalInputError(
            f"equilibrium calculation parent does not exist: {work_directory.parent}"
        )
    work_directory.mkdir(mode=0o700)

    input_hashes: dict[str, str] = {}
    for input_file in equilibrium_input_files:
        input_path = Path(input_file)
        raw_bytes, digest = _read_bytes_snapshot(input_path, "equilibrium input file")
        (work_directory / input_path.name).write_bytes(raw_bytes)
        input_hashes[input_path.name] = digest

    generator_bytes, generator_hash = _read_bytes_snapshot(
        generator_path, "equilibrium generator executable"
    )
    verified_generator = work_directory / (
        ".verified-equilibrium-generator" + generator_path.suffix
    )
    try:
        verified_generator.write_bytes(generator_bytes)
        verified_generator.chmod(0o700)
    except OSError as error:
        raise ExperimentalInputError("unable to stage equilibrium generator snapshot") from error
    # The working directory can be anchored through /proc/self/fd on Linux.
    # That descriptor closes at exec, so the interpreter must reopen a staged
    # script relative to its inherited cwd rather than through the parent fd.
    _run_equilibrium_generator([f"./{verified_generator.name}"], work_directory, timeout_seconds)

    equilibrium_file = work_directory / "equil_r_q_psi.dat"
    parameters_file = work_directory / "btor_rbig.dat"
    for output, label in (
        (equilibrium_file, "equil_r_q_psi.dat"),
        (parameters_file, "btor_rbig.dat"),
    ):
        if output.is_symlink() or not output.is_file():
            raise ExperimentalInputError(
                f"equilibrium executable must write a regular {label} output file"
            )

    equilibrium_bytes, equilibrium_hash = _read_bytes_snapshot(
        equilibrium_file, "generated equilibrium file"
    )
    _read_equilibrium_bytes(equilibrium_file, equilibrium_bytes)
    parameters = read_equilibrium_parameters(parameters_file)
    generator_provenance = {
        "source_command": command,
        "source_executable_sha256": generator_hash,
        "executed_executable": str(verified_generator),
        "executed_sha256": generator_hash,
        "input_sha256": input_hashes,
        "timeout_seconds": float(timeout_seconds),
    }
    return EquilibriumCalculation(
        equilibrium_file=equilibrium_file,
        parameters_file=parameters_file,
        parameters=parameters,
        equilibrium_sha256=equilibrium_hash,
        generator_provenance=generator_provenance,
    )


def prepare_marsf_case(
    source: MarsFInput,
    config: SimulationConfig,
    destination: Path | str,
    *,
    equilibrium_file: Path | str | None = None,
    equilibrium_parameters_file: Path | str | None = None,
    equilibrium_executable: Path | str | None = None,
    equilibrium_calculation: EquilibriumCalculation | None = None,
    equilibrium_input_files: Sequence[Path | str] = (),
    equilibrium_timeout_seconds: float = 3600.0,
    q_operation: QOperation = "preserve",
    upstream_staging_report: Path | str | None = None,
) -> PreparedExperimentalCase:
    """Stage a MARS-F source as a runnable, provenance-preserving KIM case.

    The equilibrium table must either be supplied directly or generated by an
    explicitly supplied executable.  The executable is run without a shell in
    a fresh staging directory; input files are copied there before launch.
    """

    if not isinstance(source, MarsFInput):
        raise ExperimentalInputError("source must be a MarsFInput instance")
    if not isinstance(config, SimulationConfig):
        raise ExperimentalInputError("config must be a SimulationConfig instance")
    if equilibrium_calculation is not None:
        if not isinstance(equilibrium_calculation, EquilibriumCalculation):
            raise ExperimentalInputError(
                "equilibrium_calculation must be an EquilibriumCalculation"
            )
        if equilibrium_executable is not None:
            raise ExperimentalInputError(
                "equilibrium_calculation and equilibrium_executable cannot both be supplied"
            )
        if equilibrium_file is not None and Path(equilibrium_file).absolute() != (
            equilibrium_calculation.equilibrium_file
        ):
            raise ExperimentalInputError("equilibrium_file does not match equilibrium_calculation")
        if (
            equilibrium_parameters_file is not None
            and Path(equilibrium_parameters_file).absolute()
            != equilibrium_calculation.parameters_file
        ):
            raise ExperimentalInputError(
                "equilibrium_parameters_file does not match equilibrium_calculation"
            )
        equilibrium_file = equilibrium_calculation.equilibrium_file
        equilibrium_parameters_file = equilibrium_calculation.parameters_file
    if not isinstance(q_operation, str) or q_operation not in {"preserve", "negate"}:
        raise ExperimentalInputError(
            "q operation must be explicitly selected as preserve or negate"
        )
    if equilibrium_file is not None and equilibrium_executable is not None:
        raise ExperimentalInputError("provide equilibrium_file or equilibrium_executable, not both")
    if equilibrium_file is None and equilibrium_executable is None:
        raise ExperimentalInputError(
            "an explicit equilibrium_file or equilibrium_executable is required"
        )
    if equilibrium_parameters_file is not None and equilibrium_file is None:
        raise ExperimentalInputError(
            "btor_rbig.dat must be paired with its precomputed equilibrium_file"
        )
    if not np.isfinite(equilibrium_timeout_seconds) or equilibrium_timeout_seconds <= 0.0:
        raise ExperimentalInputError("equilibrium_timeout_seconds must be positive and finite")

    equilibrium_parameters = (
        read_equilibrium_parameters(equilibrium_parameters_file)
        if equilibrium_parameters_file is not None
        else None
    )
    if equilibrium_parameters is not None:
        config = simulation_config_for_equilibrium(config, equilibrium_parameters)
    if equilibrium_calculation is not None:
        _, equilibrium_hash = _read_bytes_snapshot(
            equilibrium_calculation.equilibrium_file, "calculated equilibrium table"
        )
        if equilibrium_hash != equilibrium_calculation.equilibrium_sha256:
            raise ExperimentalInputError("calculated equilibrium table changed after generation")
        if equilibrium_parameters.sha256 != equilibrium_calculation.parameters.sha256:
            raise ExperimentalInputError("calculated btor_rbig.dat changed after generation")

    verified_snapshots = _snapshot_marsf_source(source)
    validated_upstream = _validate_upstream_staging_report(
        source,
        config,
        upstream_staging_report,
        verified_snapshots,
        equilibrium_parameters_file=equilibrium_parameters_file,
        equilibrium_parameters=equilibrium_parameters,
    )
    generator_command: list[str] | None = None
    generator_path: Path | None = None
    if equilibrium_file is None:
        assert equilibrium_executable is not None
        generator_command = _resolve_generator_command(equilibrium_executable)
        generator_path = _resolved_executable_path(generator_command[0])
        if generator_path is None or not os.access(generator_path, os.X_OK):
            raise ExperimentalInputError(
                f"equilibrium_executable does not name an executable file: "
                f"{generator_command[0]}"
            )
        _validate_equilibrium_input_names(equilibrium_input_files, generator_path)

    final_directory = Path(destination).absolute()
    if final_directory.exists() or final_directory.is_symlink():
        raise ExperimentalInputError(f"prepared case already exists: {final_directory}")
    final_directory.parent.mkdir(parents=True, exist_ok=True)

    staging_directory = Path(
        tempfile.mkdtemp(prefix=f".{final_directory.name}-", dir=final_directory.parent)
    )
    try:
        prepared = _prepare_in_directory(
            source,
            config,
            staging_directory,
            equilibrium_file=equilibrium_file,
            equilibrium_parameters_file=equilibrium_parameters_file,
            equilibrium_parameters=equilibrium_parameters,
            equilibrium_executable=equilibrium_executable,
            equilibrium_calculation=equilibrium_calculation,
            equilibrium_input_files=equilibrium_input_files,
            equilibrium_timeout_seconds=equilibrium_timeout_seconds,
            q_operation=q_operation,
            upstream_staging=validated_upstream,
            source_snapshots=verified_snapshots,
            final_directory=final_directory,
            generator_command=generator_command,
            generator_path=generator_path,
        )
        _commit_staging_directory(staging_directory, final_directory)
    except Exception:
        shutil.rmtree(staging_directory, ignore_errors=True)
        raise

    return PreparedExperimentalCase(
        directory=final_directory,
        profiles=final_directory / "profiles",
        equilibrium=final_directory / "equilibrium" / prepared.equilibrium.name,
        request=final_directory / "request.json",
        report=final_directory / "conversion_report.json",
        config=prepared.config,
        equilibrium_parameters=(
            final_directory / "equilibrium" / "btor_rbig.dat"
            if equilibrium_parameters is not None
            else None
        ),
    )


def _snapshot_marsf_source(source: MarsFInput) -> dict[str, MarsFProfileSnapshot]:
    """Read and verify every source profile from one exact byte snapshot."""

    snapshots: dict[str, MarsFProfileSnapshot] = {}
    for role, unit_field in _MARSF_PROFILE_UNIT_FIELDS.items():
        try:
            path = source.source_files[role]
            supplied = source.profiles[role]
        except KeyError as error:
            raise ExperimentalInputError(f"MARS-F source is missing {role} profile") from error
        snapshot = experimental_module.read_marsf_profile_snapshot(
            role, getattr(source.metadata, unit_field), path
        )
        if supplied.name != snapshot.profile.name or supplied.units != snapshot.profile.units:
            raise ExperimentalInputError(f"MARS-F {role} profile metadata changed after reading")
        try:
            coordinates_match = np.array_equal(supplied.coordinate, snapshot.profile.coordinate)
            values_match = np.array_equal(supplied.values, snapshot.profile.values)
        except (TypeError, ValueError) as error:
            raise ExperimentalInputError(f"MARS-F {role} profile arrays are invalid") from error
        if not coordinates_match or not values_match:
            raise ExperimentalInputError(f"MARS-F {role} profile arrays changed after reading")
        snapshots[role] = snapshot
    return snapshots


def _validate_upstream_staging_report(
    source: MarsFInput,
    config: SimulationConfig,
    report_argument: Path | str | None,
    snapshots: Mapping[str, MarsFProfileSnapshot],
    *,
    equilibrium_parameters_file: Path | str | None,
    equilibrium_parameters: EquilibriumParameters | None,
) -> dict[str, object] | None:
    """Read and validate one exact-byte BALANCE staging report.

    The report is linked by the digest embedded in the MARS-F metadata.  No
    report discovery or caller-supplied mapping is accepted here: the path and
    its bytes are the complete provenance boundary.
    """

    linked_digest = source.metadata.upstream_staging_sha256
    if linked_digest is not None and _SHA256_PATTERN.fullmatch(linked_digest) is None:
        raise ExperimentalInputError("MARS-F metadata contains an invalid staging report digest")
    if report_argument is None:
        # Direct preparation remains a supported path.  A staging link is only
        # consumed when its report is explicitly supplied; the output records
        # the absent upstream report as null rather than auto-discovering it.
        return None
    if not isinstance(report_argument, (Path, str)):
        raise ExperimentalInputError("upstream_staging_report must be a path")
    if linked_digest is None:
        raise ExperimentalInputError(
            "upstream staging report was supplied without a linked metadata digest"
        )
    report_path = Path(report_argument)
    try:
        report_bytes = report_path.read_bytes()
    except (OSError, UnicodeError) as error:
        raise ExperimentalInputError(
            f"unable to read upstream staging report: {report_path}"
        ) from error
    report_digest = hashlib.sha256(report_bytes).hexdigest()
    if report_digest != linked_digest:
        raise ExperimentalInputError(
            "upstream staging report digest does not match MARS-F metadata"
        )
    try:
        payload = json.loads(report_bytes, parse_constant=_reject_json_constant)
    except (
        TypeError,
        UnicodeError,
        ValueError,
        OverflowError,
        RecursionError,
        json.JSONDecodeError,
    ) as error:
        raise ExperimentalInputError("upstream staging report is not valid JSON") from error
    if not isinstance(payload, dict):
        raise ExperimentalInputError("upstream staging report must be a JSON object")
    schema_version = payload.get("schema_version")
    if type(schema_version) is not int or (
        (schema_version == 1 and set(payload) != _STAGING_REPORT_V1_KEYS)
        or (schema_version == 2 and set(payload) != _STAGING_REPORT_V2_KEYS)
        or schema_version not in {1, 2}
    ):
        raise ExperimentalInputError("upstream staging report schema is invalid")

    source_metadata_payload = payload["source_metadata"]
    if not isinstance(source_metadata_payload, dict):
        raise ExperimentalInputError("upstream staging report source_metadata must be an object")
    source_metadata_keys = {
        "schema_version",
        "source",
        "coordinate",
        "coordinate_unit",
        "density_unit",
        "electron_temperature_unit",
        "ion_temperature_unit",
        "toroidal_rotation_unit",
    }
    if set(source_metadata_payload) != source_metadata_keys:
        raise ExperimentalInputError("upstream staging report source_metadata is not complete")
    if type(source_metadata_payload["schema_version"]) is not int:
        raise ExperimentalInputError(
            "upstream staging report source_metadata schema version is invalid"
        )
    if any(
        not isinstance(source_metadata_payload[field], str)
        for field in source_metadata_keys - {"schema_version"}
    ):
        raise ExperimentalInputError("upstream staging report source_metadata types are invalid")
    try:
        balance_metadata = BalanceMetadata.model_validate(source_metadata_payload)
    except Exception as error:
        raise ExperimentalInputError(
            "upstream staging report source_metadata is invalid"
        ) from error
    if source_metadata_payload != balance_metadata.model_dump(mode="json"):
        raise ExperimentalInputError("upstream staging report source_metadata is not canonical")

    source_basenames = _strict_mapping(payload, "source_basenames", _BALANCE_PROFILE_ROLES)
    if any(
        not isinstance(value, str) or not value or Path(value).name != value
        for value in source_basenames.values()
    ):
        raise ExperimentalInputError("upstream staging report source_basenames are invalid")
    source_hashes = _strict_hash_mapping(payload, "source_hashes", _BALANCE_PROFILE_ROLES)
    derived_hashes = _strict_hash_mapping(
        payload, "derived_hashes", tuple(_MARSF_PROFILE_FILENAMES.values())
    )
    reported_parameters = payload.get("equilibrium_parameters")
    if (schema_version == 2) != (reported_parameters is not None):
        raise ExperimentalInputError(
            "upstream staging report schema does not match its equilibrium parameters"
        )
    if equilibrium_parameters_file is None:
        if reported_parameters is not None:
            raise ExperimentalInputError(
                "upstream staging equilibrium parameters were not supplied to preparation"
            )
    else:
        if equilibrium_parameters is None or not isinstance(reported_parameters, dict):
            raise ExperimentalInputError("upstream staging equilibrium parameters are invalid")
        expected_parameters = {
            "source_basename": Path(equilibrium_parameters_file).name,
            "sha256": equilibrium_parameters.sha256,
            "btor_gauss": equilibrium_parameters.btor_gauss,
            "r_big_cm": equilibrium_parameters.r_big_cm,
        }
        if reported_parameters != expected_parameters:
            raise ExperimentalInputError(
                "upstream staging equilibrium parameters do not match the calculation output"
            )
        if float(config.setup.btor) != equilibrium_parameters.btor_gauss:
            raise ExperimentalInputError("KIM btor does not match the equilibrium calculation")
        if float(config.setup.major_radius) != equilibrium_parameters.r_big_cm:
            raise ExperimentalInputError(
                "KIM major_radius does not match the equilibrium calculation"
            )

    major_radius = payload["major_radius_cm"]
    try:
        major_radius_value = float(major_radius)
    except (OverflowError, TypeError, ValueError) as error:
        raise ExperimentalInputError(
            "upstream staging report major_radius_cm is invalid"
        ) from error
    if (
        type(major_radius) not in {int, float}
        or not np.isfinite(major_radius_value)
        or major_radius_value <= 0.0
    ):
        raise ExperimentalInputError("upstream staging report major_radius_cm is invalid")
    major_radius = major_radius_value
    if major_radius != float(config.setup.major_radius):
        raise ExperimentalInputError(
            "upstream staging major radius does not match the KIM configuration"
        )

    coordinate_mapping = payload["coordinate_mapping"]
    if not isinstance(coordinate_mapping, dict):
        raise ExperimentalInputError("upstream staging report coordinate_mapping is invalid")
    expected_mapping = {
        "source": balance_metadata.coordinate,
        "source_unit": balance_metadata.coordinate_unit,
        "target": "sqrt_psiN",
        "target_unit": "1",
        "operation": (
            "rho_pol = sqrt(psi_pol_norm)"
            if balance_metadata.coordinate == "rho_pol"
            else "preserve"
        ),
    }
    if coordinate_mapping != expected_mapping:
        raise ExperimentalInputError("upstream staging report coordinate mapping is inconsistent")

    equilibrium_provenance = payload["equilibrium_provenance"]
    if (
        not isinstance(equilibrium_provenance, str)
        or not equilibrium_provenance.strip()
        or equilibrium_provenance != source.metadata.equilibrium_provenance
    ):
        raise ExperimentalInputError(
            "upstream staging report equilibrium provenance is inconsistent"
        )

    operations = payload["operations"]
    expected_operations = [
        {
            "quantity": "density",
            "source_unit": balance_metadata.density_unit,
            "target_unit": balance_metadata.density_unit,
            "factor": 1.0,
            "operation": "preserve",
        },
        {
            "quantity": "electron_temperature",
            "source_unit": "eV",
            "target_unit": "eV",
            "factor": 1.0,
            "operation": "preserve",
        },
        {
            "quantity": "ion_temperature",
            "source_unit": "eV",
            "target_unit": "eV",
            "factor": 1.0,
            "operation": "preserve",
        },
        {
            "quantity": "toroidal_velocity",
            "source_unit": "rad/s",
            "target_unit": "cm/s",
            "factor": major_radius,
            "operation": "omega_to_v_phi",
            "parameters": {"major_radius_cm": major_radius},
        },
    ]
    if not _operations_have_strict_semantics(operations, expected_operations):
        raise ExperimentalInputError("upstream staging report operations are inconsistent")

    expected_source_metadata = {
        "source": balance_metadata.source,
        "coordinate": "sqrt_psiN",
        "coordinate_unit": "1",
        "density_unit": balance_metadata.density_unit,
        "electron_temperature_unit": balance_metadata.electron_temperature_unit,
        "ion_temperature_unit": balance_metadata.ion_temperature_unit,
        "toroidal_velocity_unit": "cm/s",
        "equilibrium_provenance": equilibrium_provenance,
    }
    for field, expected in expected_source_metadata.items():
        if getattr(source.metadata, field) != expected:
            raise ExperimentalInputError(f"MARS-F metadata does not match upstream staging {field}")

    current_derived_hashes = {
        filename: snapshots["toroidal_velocity" if role == "toroidal_rotation" else role].sha256
        for role, filename in _MARSF_PROFILE_FILENAMES.items()
    }
    if derived_hashes != current_derived_hashes:
        raise ExperimentalInputError(
            "upstream staging derived hashes do not match the current MARS-F quartet"
        )

    # Reconstruct only fields that passed validation; arbitrary JSON members
    # never cross into the conversion report.
    return {
        "schema_version": 1,
        "source_basenames": dict(source_basenames),
        "source_hashes": dict(source_hashes),
        "source_metadata": balance_metadata.model_dump(mode="json"),
        "derived_hashes": dict(derived_hashes),
        "major_radius_cm": major_radius,
        "equilibrium_parameters": reported_parameters,
        "coordinate_mapping": dict(expected_mapping),
        "equilibrium_provenance": equilibrium_provenance,
        "operations": expected_operations,
        "provenance_sha256": report_digest,
    }


def _strict_mapping(
    payload: Mapping[str, object], key: str, expected_keys: Sequence[str]
) -> dict[str, str]:
    value = payload[key]
    if not isinstance(value, dict) or set(value) != set(expected_keys):
        raise ExperimentalInputError(f"upstream staging report {key} is invalid")
    if any(not isinstance(item, str) for item in value.values()):
        raise ExperimentalInputError(f"upstream staging report {key} is invalid")
    return dict(value)


def _strict_hash_mapping(
    payload: Mapping[str, object], key: str, expected_keys: Sequence[str]
) -> dict[str, str]:
    value = _strict_mapping(payload, key, expected_keys)
    if any(_SHA256_PATTERN.fullmatch(item) is None for item in value.values()):
        raise ExperimentalInputError(f"upstream staging report {key} contains invalid hashes")
    return value


def _operations_have_strict_semantics(
    operations: object, expected: list[dict[str, object]]
) -> bool:
    if not isinstance(operations, list) or len(operations) != len(expected):
        return False
    for actual, reference in zip(operations, expected, strict=True):
        if not isinstance(actual, dict) or set(actual) != set(reference):
            return False
        for key, expected_value in reference.items():
            actual_value = actual[key]
            if isinstance(expected_value, str):
                if not isinstance(actual_value, str) or actual_value != expected_value:
                    return False
            if isinstance(expected_value, float):
                if type(actual_value) not in {int, float}:
                    return False
                try:
                    numeric_value = float(actual_value)
                except (OverflowError, TypeError, ValueError):
                    return False
                if not np.isfinite(numeric_value) or numeric_value != expected_value:
                    return False
            if isinstance(expected_value, dict):
                if not isinstance(actual_value, dict) or actual_value != expected_value:
                    return False
                for nested_value in actual_value.values():
                    if type(nested_value) not in {int, float}:
                        return False
                    try:
                        numeric_value = float(nested_value)
                    except (OverflowError, TypeError, ValueError):
                        return False
                    if not np.isfinite(numeric_value):
                        return False
    return True


def _commit_staging_directory(staging: Path, destination: Path) -> None:
    try:
        destination.mkdir()
    except FileExistsError as error:
        raise ExperimentalInputError(f"prepared case already exists: {destination}") from error
    try:
        for child in staging.iterdir():
            target = destination / child.name
            if target.exists() or target.is_symlink():
                raise ExperimentalInputError(f"prepared case target already exists: {target}")
            child.rename(target)
        staging.rmdir()
    except Exception:
        shutil.rmtree(destination, ignore_errors=True)
        raise


def _prepare_in_directory(
    source: MarsFInput,
    config: SimulationConfig,
    directory: Path,
    *,
    equilibrium_file: Path | str | None,
    equilibrium_parameters_file: Path | str | None,
    equilibrium_parameters: EquilibriumParameters | None,
    equilibrium_executable: Path | str | None,
    equilibrium_calculation: EquilibriumCalculation | None,
    equilibrium_input_files: Sequence[Path | str],
    equilibrium_timeout_seconds: float,
    q_operation: QOperation,
    upstream_staging: dict[str, object] | None,
    source_snapshots: Mapping[str, MarsFProfileSnapshot],
    final_directory: Path,
    generator_command: Sequence[str] | None,
    generator_path: Path | None,
) -> PreparedExperimentalCase:
    source_directory = directory / "source"
    equilibrium_directory = directory / "equilibrium"
    profiles_directory = directory / "profiles"
    source_directory.mkdir()
    equilibrium_directory.mkdir()
    profiles_directory.mkdir()

    source_hashes: dict[str, str] = {}
    for _name in _MARSF_PROFILE_UNIT_FIELDS:
        try:
            snapshot = source_snapshots[_name]
        except KeyError as error:
            raise ExperimentalInputError(f"MARS-F source is missing {_name} snapshot") from error
        source_path = snapshot.profile.path
        target = source_directory / source_path.name
        target.write_bytes(snapshot.raw_bytes)
        source_hashes[f"source/{source_path.name}"] = snapshot.sha256
    (source_directory / "metadata.json").write_text(
        json.dumps(source.metadata.model_dump(mode="json"), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    generated = equilibrium_file is None
    generator_provenance: dict[str, object] | None = (
        dict(equilibrium_calculation.generator_provenance)
        if equilibrium_calculation is not None
        else None
    )
    if generated:
        assert equilibrium_executable is not None
        assert generator_command is not None
        assert generator_path is not None
        input_names: set[str] = set()
        for input_file in equilibrium_input_files:
            input_path = Path(input_file)
            if not input_path.is_file():
                raise ExperimentalInputError(f"equilibrium input file does not exist: {input_path}")
            input_name = input_path.name.casefold()
            if input_name in input_names:
                raise ExperimentalInputError(
                    f"duplicate equilibrium input basename: {input_path.name}"
                )
            input_names.add(input_name)
            input_bytes, input_hash = _read_bytes_snapshot(input_path, "equilibrium input file")
            (equilibrium_directory / input_path.name).write_bytes(input_bytes)
            source_hashes[f"equilibrium_input/{input_path.name}"] = input_hash
        generator_bytes, generator_hash = _read_bytes_snapshot(
            generator_path, "equilibrium generator executable"
        )
        verified_generator = equilibrium_directory / (
            ".verified-equilibrium-generator" + generator_path.suffix
        )
        try:
            verified_generator.write_bytes(generator_bytes)
            verified_generator.chmod(0o700)
        except OSError as error:
            raise ExperimentalInputError(
                "unable to create the verified equilibrium generator copy"
            ) from error
        verified_command = [f"./{verified_generator.name}"]
        reported_verified_path = f"equilibrium/{verified_generator.name}"
        generator_provenance = {
            "command": [reported_verified_path],
            "source_command": generator_command,
            "executable": str(generator_path),
            "sha256": generator_hash,
            "command_file_hashes": {str(generator_path): generator_hash},
            "executed_command": [reported_verified_path],
            "executed_executable": reported_verified_path,
            "executed_sha256": generator_hash,
            "hash_basis": "exact bytes copied to executed_executable",
            "timeout_seconds": equilibrium_timeout_seconds,
        }
        _run_equilibrium_generator(
            verified_command, equilibrium_directory, equilibrium_timeout_seconds
        )
        equilibrium_source = equilibrium_directory / "equil_r_q_psi.dat"
        if equilibrium_source.is_symlink() or not equilibrium_source.is_file():
            raise ExperimentalInputError(
                "equilibrium executable must write a regular equil_r_q_psi.dat file"
            )
        equilibrium_bytes, equilibrium_source_hash = _read_bytes_snapshot(
            equilibrium_source, "generated equilibrium file"
        )
        equilibrium_operation = "generated by supplied equilibrium executable"
    else:
        equilibrium_path = Path(equilibrium_file).absolute()
        equilibrium_bytes, equilibrium_source_hash = _read_bytes_snapshot(
            equilibrium_path, "equilibrium file"
        )
        staged_equilibrium = equilibrium_directory / equilibrium_path.name
        staged_equilibrium.write_bytes(equilibrium_bytes)
        equilibrium_source = staged_equilibrium
        equilibrium_operation = (
            "generated by supplied equilibrium executable"
            if equilibrium_calculation is not None
            else "copied from explicit equilibrium_file"
        )

    equilibrium_hash_name = equilibrium_source.name
    source_hashes[f"equilibrium/{equilibrium_hash_name}"] = equilibrium_source_hash
    staged_parameters: Path | None = None
    if equilibrium_parameters_file is not None:
        assert equilibrium_parameters is not None
        parameter_bytes, parameter_hash = _read_bytes_snapshot(
            Path(equilibrium_parameters_file), "btor_rbig.dat equilibrium output"
        )
        if parameter_hash != equilibrium_parameters.sha256:
            raise ExperimentalInputError("btor_rbig.dat changed after its values were read")
        staged_parameters = equilibrium_directory / "btor_rbig.dat"
        staged_parameters.write_bytes(parameter_bytes)
        source_hashes["equilibrium/btor_rbig.dat"] = parameter_hash
    radius, q, psi_n = _read_equilibrium_bytes(equilibrium_source, equilibrium_bytes)
    verified_profiles = {role: snapshot.profile for role, snapshot in source_snapshots.items()}
    output_grid, profile_values, coordinate_operation = _prepare_profiles(
        source, verified_profiles, radius, q, psi_n, profiles_directory, config.profiles
    )
    _validate_radial_coverage(config, output_grid)
    q_values = (
        q if source.metadata.coordinate == "sqrt_psiN" else _interpolate(q, radius, output_grid)
    )
    q_factor = 1.0 if q_operation == "preserve" else -1.0
    q_values = q_values * q_factor
    _write_profile(profiles_directory / config.profiles.safety_factor_file, output_grid, q_values)
    prepared_output_hashes = {
        f"profiles/{filename}": _sha256(profiles_directory / filename)
        for filename in (
            config.profiles.density_file,
            config.profiles.electron_temperature_file,
            config.profiles.ion_temperature_file,
            config.profiles.toroidal_velocity_file,
            config.profiles.safety_factor_file,
        )
    }

    prepared_config = config.model_copy(
        update={
            "profiles": config.profiles.model_copy(
                update={"directory": final_directory / "profiles"}
            )
        }
    )
    request_payload = prepared_config.model_dump(mode="json")
    request_payload["profiles"]["directory"] = "./profiles"
    request_path = directory / "request.json"
    request_path.write_text(
        json.dumps(request_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    q_conversion = ConversionOperation(
        quantity="q",
        source_unit="1",
        target_unit="1",
        factor=q_factor,
        operation=q_operation,
    )
    report = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source": source.metadata.model_dump(mode="json"),
        "target": "KIM-CGS-r_eff",
        "equilibrium": {
            "source": str(Path(equilibrium_file).absolute()) if equilibrium_file else None,
            "staged_file": f"equilibrium/{equilibrium_source.name}",
            "operation": equilibrium_operation,
            "source_hash": equilibrium_source_hash,
        },
        "equilibrium_parameters": (
            {
                "staged_file": "equilibrium/btor_rbig.dat",
                "source_hash": equilibrium_parameters.sha256,
                "btor_gauss": equilibrium_parameters.btor_gauss,
                "r_big_cm": equilibrium_parameters.r_big_cm,
            }
            if equilibrium_parameters is not None
            else None
        ),
        "generator": generator_provenance,
        "coordinate_operation": (
            "natural cubic interpolation from sqrt_psiN to equilibrium r_eff"
            if source.metadata.coordinate == "sqrt_psiN"
            else coordinate_operation["method"]
        ),
        "coordinate_mapping": coordinate_operation,
        "source_hashes": source_hashes,
        "operations": [*_unit_operations(source), q_conversion.model_dump(mode="json")],
        "output_grid_points": int(output_grid.size),
        "prepared_output_hashes": prepared_output_hashes,
        "software": _software_identity(),
        "git": _git_identity(),
        "comparison": {"domain": None},
        "upstream_staging": upstream_staging,
    }
    report_path = directory / "conversion_report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return PreparedExperimentalCase(
        directory=directory,
        profiles=profiles_directory,
        equilibrium=equilibrium_source,
        request=request_path,
        report=report_path,
        config=prepared_config,
        equilibrium_parameters=staged_parameters,
    )


def _run_equilibrium_generator(
    command: Sequence[str], working_directory: Path, timeout_seconds: float
) -> None:
    if not command:
        raise ExperimentalInputError("equilibrium_executable must not be empty")
    try:
        completed = subprocess.run(
            command,
            cwd=working_directory,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as error:
        raise ExperimentalInputError(
            f"equilibrium generator timed out after {timeout_seconds:g} seconds"
        ) from error
    except OSError as error:
        raise ExperimentalInputError(
            f"unable to execute equilibrium generator: {command[0]}"
        ) from error
    if completed.returncode != 0:
        details = completed.stderr.strip() or completed.stdout.strip()
        raise ExperimentalInputError(
            f"equilibrium generator failed with exit code {completed.returncode}: {details}"
        )
    (working_directory / "generator.stdout").write_text(completed.stdout, encoding="utf-8")
    (working_directory / "generator.stderr").write_text(completed.stderr, encoding="utf-8")


def _resolve_generator_command(executable: Path | str) -> list[str]:
    candidate = Path(executable).expanduser()
    if candidate.is_absolute():
        return [str(candidate)]
    if (
        (Path.cwd() / candidate).is_file()
        or os.sep in str(candidate)
        or (os.altsep and os.altsep in str(candidate))
    ):
        return [str((Path.cwd() / candidate).resolve())]
    return [str(executable)]


def _validate_equilibrium_input_names(
    input_files: Sequence[Path | str], generator_path: Path
) -> None:
    """Reject input names owned by generated equilibrium preparation artifacts."""

    reserved_names = {
        ".verified-equilibrium-generator" + generator_path.suffix,
        "btor_rbig.dat",
        "equil_r_q_psi.dat",
        "generator.stdout",
        "generator.stderr",
    }
    reserved_names = {name.casefold() for name in reserved_names}
    input_names: set[str] = set()
    for input_file in input_files:
        input_path = Path(input_file)
        if not input_path.is_file():
            raise ExperimentalInputError(f"equilibrium input file does not exist: {input_path}")
        input_name = input_path.name.casefold()
        if input_name in input_names:
            raise ExperimentalInputError(f"duplicate equilibrium input basename: {input_path.name}")
        if input_name in reserved_names:
            raise ExperimentalInputError(
                f"equilibrium input basename is reserved by preparation: {input_path.name}"
            )
        input_names.add(input_name)


def _resolved_executable_path(command: str) -> Path | None:
    candidate = Path(command)
    if candidate.is_file():
        return candidate.absolute()
    located = shutil.which(command)
    return Path(located).absolute() if located else None


def _software_identity() -> dict[str, str] | None:
    try:
        version = importlib.metadata.version("kamel-kim")
    except Exception:
        return None
    return {"name": "kamel-kim", "version": version}


def _git_identity() -> dict[str, object] | None:
    try:
        metadata = executable_module.discover_kamel_git_metadata()
    except Exception:
        return None
    if metadata is None:
        return None
    return {"commit": metadata.commit, "dirty": metadata.dirty}


def _read_bytes_snapshot(path: Path, label: str) -> tuple[bytes, str]:
    try:
        raw_bytes = path.read_bytes()
    except (OSError, UnicodeError) as error:
        raise ExperimentalInputError(f"{path}: unable to read {label}") from error
    return raw_bytes, hashlib.sha256(raw_bytes).hexdigest()


def _read_equilibrium(
    path: Path,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Compatibility wrapper that parses one equilibrium byte snapshot."""

    raw_bytes, _digest = _read_bytes_snapshot(path, "equilibrium file")
    return _read_equilibrium_bytes(path, raw_bytes)


def _read_equilibrium_bytes(
    path: Path, raw_bytes: bytes
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    try:
        lines = raw_bytes.decode("utf-8").splitlines()
    except UnicodeError as error:
        raise ExperimentalInputError(f"{path}: unable to read equilibrium file") from error
    rows: list[tuple[float, float, float]] = []
    for line_number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        columns = stripped.split()
        if len(columns) < 3:
            raise ExperimentalInputError(f"{path}: row {line_number} needs r_eff, q, and psi")
        try:
            row = tuple(float(value.replace("D", "E").replace("d", "e")) for value in columns[:3])
        except (OverflowError, ValueError) as error:
            raise ExperimentalInputError(f"{path}: row {line_number} is not numeric") from error
        if not np.all(np.isfinite(row)):
            raise ExperimentalInputError(f"{path}: row {line_number} contains non-finite values")
        rows.append(row)
    if len(rows) < 2:
        raise ExperimentalInputError(f"{path}: equilibrium table needs at least two rows")
    data = np.asarray(rows, dtype=np.float64)
    if np.any(np.diff(data[:, 0]) <= 0):
        raise ExperimentalInputError(f"{path}: radius column must be strictly increasing")
    if data[-1, 2] == 0:
        raise ExperimentalInputError(f"{path}: final psi value must not be zero")
    psi_n = data[:, 2] / data[-1, 2]
    if np.any(np.diff(psi_n) <= 0):
        raise ExperimentalInputError(f"{path}: normalized psi column must be strictly increasing")
    return data[:, 0], data[:, 1], psi_n


def _prepare_profiles(
    source: MarsFInput,
    profiles: Mapping[str, ExperimentalProfile],
    radius: NDArray[np.float64],
    q: NDArray[np.float64],
    psi_n: NDArray[np.float64],
    profiles_directory: Path,
    profile_config: ProfileConfig,
) -> tuple[NDArray[np.float64], dict[str, NDArray[np.float64]], dict[str, str]]:
    source_metadata = _conversion_metadata(source)
    if source.metadata.coordinate == "sqrt_psiN":
        output_grid = radius
        coordinate = psi_n
        coordinate_operation = {
            "source_coordinate": "sqrt_psiN",
            "target_coordinate": "r_eff",
            "method": "natural cubic interpolation",
        }
        squared_coordinates: dict[str, NDArray[np.float64]] = {}
        with np.errstate(over="ignore", invalid="ignore"):
            for name, profile in profiles.items():
                if np.any(profile.coordinate < 0.0):
                    raise ExperimentalInputError("sqrt_psiN coordinate must be nonnegative")
                squared_coordinate = np.square(profile.coordinate)
                if not np.all(np.isfinite(squared_coordinate)):
                    raise ExperimentalInputError(
                        f"{name} squared sqrt_psiN coordinate is not finite"
                    )
                squared_coordinates[name] = squared_coordinate
        profile_values = {
            name: _convert_profile(
                _interpolate(profile.values, squared_coordinates[name], coordinate),
                name,
                source_metadata,
            )
            for name, profile in profiles.items()
        }
        for name in profiles:
            _require_covered(squared_coordinates[name], coordinate, "sqrt_psiN profile")
    else:
        scale = 100.0 if source.metadata.coordinate_unit == "m" else 1.0
        output_grid = profiles["density"].coordinate * scale
        coordinate_operation = {
            "source_coordinate": "r_eff",
            "target_coordinate": "r_eff",
            "method": (
                "converted explicit r_eff grid from m to cm"
                if scale != 1.0
                else "preserved explicit r_eff grid"
            ),
        }
        for profile in profiles.values():
            if not np.array_equal(profile.coordinate, profiles["density"].coordinate):
                raise ExperimentalInputError("r_eff MARS-F profiles must share one coordinate grid")
        _require_covered(radius, output_grid, "r_eff profile")
        profile_values = {
            name: _convert_profile(profile.values, name, source_metadata)
            for name, profile in profiles.items()
        }
    _validate_profile_filenames(profile_config)
    for name, values in profile_values.items():
        _write_profile(
            profiles_directory / _profile_filename(name, profile_config), output_grid, values
        )
    return output_grid, profile_values, coordinate_operation


def _validate_radial_coverage(config: SimulationConfig, output_grid: NDArray[np.float64]) -> None:
    if output_grid[0] > config.grid.radial_minimum or output_grid[-1] < config.grid.plasma_radius:
        raise ExperimentalInputError(
            "prepared profile grid does not cover the configured radial minimum and plasma radius"
        )


def _conversion_metadata(source: MarsFInput) -> SourceMetadata:
    metadata = source.metadata
    if metadata.electron_temperature_unit != "eV" or metadata.ion_temperature_unit != "eV":
        raise ExperimentalInputError(
            "only eV electron and ion temperatures are approved for KIM preparation"
        )
    return SourceMetadata(
        source=metadata.source,
        coordinate="r_eff",
        radius_unit="cm",
        density_unit=metadata.density_unit,
        temperature_unit="eV",
        magnetic_field_unit="G",
        magnetic_field_sign="preserve",
        electric_field_unit="statV/cm",
        velocity_unit=metadata.toroidal_velocity_unit,
        frequency_unit="rad/s",
        frequency_convention="signed_omega",
        time_phase="exp(-i omega t)",
        perturbation_phase="preserve",
        equilibrium_provenance=metadata.equilibrium_provenance,
        q_unit="1",
        mode_convention="signed",
        fourier_phase="exp(+i k r)",
        complex_phase="preserve",
        resonance_convention="q=-m/n",
    )


def _convert_profile(
    values: NDArray[np.float64], name: str, metadata: SourceMetadata
) -> NDArray[np.float64]:
    quantity = {
        "density": "density",
        "electron_temperature": "temperature",
        "ion_temperature": "temperature",
        "toroidal_velocity": "toroidal_velocity",
    }[name]
    return convert_quantity(values, quantity, metadata)


def _unit_operations(source: MarsFInput) -> list[dict[str, object]]:
    metadata = source.metadata
    operations: list[dict[str, object]] = []
    if metadata.coordinate == "r_eff":
        operations.append(
            {
                "quantity": "radius",
                "source_unit": metadata.coordinate_unit,
                "target_unit": "cm",
                "factor": 100.0 if metadata.coordinate_unit == "m" else 1.0,
            }
        )
    if metadata.density_unit in {"1/m^3", "m^-3"}:
        operations.append(
            {
                "quantity": "density",
                "source_unit": metadata.density_unit,
                "target_unit": "1/cm^3",
                "factor": 1.0e-6,
            }
        )
    else:
        operations.append(
            {
                "quantity": "density",
                "source_unit": metadata.density_unit,
                "target_unit": "1/cm^3",
                "factor": 1.0,
            }
        )
    velocity_factor = 100.0 if metadata.toroidal_velocity_unit == "m/s" else 1.0
    operations.append(
        {
            "quantity": "toroidal_velocity",
            "source_unit": metadata.toroidal_velocity_unit,
            "target_unit": "cm/s",
            "factor": velocity_factor,
        }
    )
    for name, unit in (
        ("electron_temperature", metadata.electron_temperature_unit),
        ("ion_temperature", metadata.ion_temperature_unit),
    ):
        operations.append(
            {"quantity": name, "source_unit": unit, "target_unit": "eV", "factor": 1.0}
        )
    return operations


def _require_covered(
    source_grid: NDArray[np.float64], target_grid: NDArray[np.float64], label: str
) -> None:
    tolerance = 1.0e-12
    if target_grid[0] < source_grid[0] - tolerance or target_grid[-1] > source_grid[-1] + tolerance:
        raise ExperimentalInputError(f"{label} range is outside equilibrium/profile range")


def _interpolate(
    values: NDArray[np.float64],
    source_grid: NDArray[np.float64],
    target_grid: NDArray[np.float64],
) -> NDArray[np.float64]:
    _require_covered(source_grid, target_grid, "profile")
    if source_grid.size == 2:
        return np.interp(target_grid, source_grid, values)
    h = np.diff(source_grid)
    alpha = 3.0 * ((values[2:] - values[1:-1]) / h[1:] - (values[1:-1] - values[:-2]) / h[:-1])
    diagonal = np.ones(source_grid.size)
    upper = np.zeros(source_grid.size - 1)
    lower = np.zeros(source_grid.size - 1)
    rhs = np.zeros(source_grid.size)
    for index in range(1, source_grid.size - 1):
        lower[index - 1] = h[index - 1]
        diagonal[index] = 2.0 * (h[index - 1] + h[index])
        upper[index] = h[index]
        rhs[index] = alpha[index - 1]
    for index in range(1, source_grid.size):
        factor = lower[index - 1] / diagonal[index - 1]
        diagonal[index] -= factor * upper[index - 1]
        rhs[index] -= factor * rhs[index - 1]
    curvature = np.zeros(source_grid.size)
    for index in range(source_grid.size - 2, -1, -1):
        curvature[index] = (rhs[index] - upper[index] * curvature[index + 1]) / diagonal[index]
    evaluation_grid = np.clip(target_grid, source_grid[0], source_grid[-1])
    interval = np.searchsorted(source_grid, evaluation_grid, side="right") - 1
    interval = np.clip(interval, 0, source_grid.size - 2)
    distance = evaluation_grid - source_grid[interval]
    slope = (values[interval + 1] - values[interval]) / h[interval] - h[interval] * (
        curvature[interval + 1] + 2.0 * curvature[interval]
    ) / 3.0
    return (
        values[interval]
        + slope * distance
        + curvature[interval] * distance**2
        + (curvature[interval + 1] - curvature[interval]) * distance**3 / (3.0 * h[interval])
    )


def _write_profile(
    path: Path, coordinate: NDArray[np.float64], values: NDArray[np.float64]
) -> None:
    np.savetxt(path, np.column_stack((coordinate, values)), fmt="%.15e")


def _profile_filename(name: str, config: ProfileConfig) -> str:
    return {
        "density": config.density_file,
        "electron_temperature": config.electron_temperature_file,
        "ion_temperature": config.ion_temperature_file,
        "toroidal_velocity": config.toroidal_velocity_file,
    }[name]


def _validate_profile_filenames(config: ProfileConfig) -> None:
    filenames = (
        config.density_file,
        config.electron_temperature_file,
        config.ion_temperature_file,
        config.toroidal_velocity_file,
        config.radial_electric_field_file,
        config.safety_factor_file,
    )
    if len(set(filenames)) != len(filenames):
        raise ExperimentalInputError("profile filenames must be unique")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


__all__ = ["PreparedExperimentalCase", "prepare_marsf_case"]
