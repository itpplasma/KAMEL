"""Opt-in, solver-free characterization of an explicit BALANCE case.

This module is deliberately a small orchestration boundary.  It stages the
four source profiles, prepares a KIM input case, and measures that case
against a read-only QL-Balance oracle.  It never launches either scientific
solver; the only process that may be started is the explicitly supplied
equilibrium preprocessor.  Failed work may retain an incomplete private
staging container so cleanup never follows a swapped path.
"""

from __future__ import annotations

import hashlib
import json
import os
import secrets
import stat
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Mapping, Sequence

import numpy as np
from kim.acceptance_comparison import ComparisonResult, compare_profiles
from kim.balance_adoption import (
    _publish_staging_directory_at,
    stage_balance_marsf_quartet,
)
from kim.config import SimulationConfig
from kim.errors import ExperimentalInputError
from kim.importers.balance import BalanceMetadata, read_balance_profiles
from kim.importers.experimental import read_marsf_profiles
from kim.namelist import load_namelist_bytes
from kim.preparation import (
    prepare_marsf_case,
    read_equilibrium_parameters,
    run_equilibrium_calculation,
    simulation_config_for_equilibrium,
)
from kim.ql_balance_oracle import read_ql_balance_oracle

BalanceStatus = Literal["MEASURED", "PASS", "FAIL", "PARTIALLY_EVALUATED", "UNAVAILABLE"]

_BALANCE_ROLES = (
    "density",
    "electron_temperature",
    "ion_temperature",
    "toroidal_rotation",
)
_PROFILE_ROLES = ("n", "Te", "Ti", "Vz", "q")
_ROLE_ALIASES = {
    "density": "n",
    "n": "n",
    "electron_temperature": "Te",
    "Te": "Te",
    "ion_temperature": "Ti",
    "Ti": "Ti",
    "toroidal_rotation": "Vz",
    "toroidal_velocity": "Vz",
    "Vz": "Vz",
    "safety_factor": "q",
    "q": "q",
}
_PROFILE_CONFIG_FIELDS = {
    "n": "density_file",
    "Te": "electron_temperature_file",
    "Ti": "ion_temperature_file",
    "Vz": "toroidal_velocity_file",
    "q": "safety_factor_file",
}
_PROFILE_UNITS = {"n": "1/cm^3", "Te": "eV", "Ti": "eV", "Vz": "cm/s", "q": "1"}
_ORACLE_FIELDS = {"n": "n", "Te": "Te", "Ti": "Ti", "Vz": "Vz", "q": "q"}
_SOURCE_ROLES = {
    "n": "density",
    "Te": "electron_temperature",
    "Ti": "ion_temperature",
    "Vz": "toroidal_rotation",
    "q": "equilibrium_q",
}
_ORACLE_DATASETS = {
    "n": "preprocprof/n",
    "Te": "preprocprof/Te",
    "Ti": "preprocprof/Ti",
    "Vz": "preprocprof/Vz",
    "q": "preprocprof/q",
}
_METRICS = ("absolute_rms", "absolute_max", "relative_rms", "relative_max")
_LIMITATIONS = {
    "input_adoption_only": True,
    "no_physical_response_validation": True,
    "solvers_not_run": True,
    "retained_staging_possible": True,
    "retained_staging_success_state": "empty_container",
    "retained_staging_failure_state": "incomplete_container_possible",
}

_PAYLOAD_NAME = "payload"
_RETAINED_STAGING_REASON = "safe ownership policy / portable rename limitation"


@dataclass(frozen=True)
class BalanceCharacterizationRequest:
    """Explicit inputs for :func:`characterize_balance`.

    ``equilibrium_file`` and ``equilibrium_parameters_file`` are the canonical,
    colocated outputs of one precomputed calculation for synthetic cases.  The
    production route supplies ``original_equilibrium`` and
    ``equilibrium_executable``; both outputs are generated once from those
    inputs before profile staging.

    ``destination.parent`` must already exist; this workflow never creates a
    destination parent.  A successful report exposes the empty private
    staging container; remove it only when no characterization is running.
    Failed work may retain an incomplete container.  ``expected_sha256`` is
    optional.  When supplied, it must cover every canonical identity in this
    request; hashless characterization remains a valid measured outcome and
    is reported as such.

    The legacy positional argument order is preserved. The BALANCE route keeps
    q exactly as emitted by ``fouriermodes.x``; ``domains`` and both
    interpolation settings remain runtime-required. ``relative_floors`` is
    optional; omitting it reports absolute metrics only.
    """

    density: Path | str
    electron_temperature: Path | str
    ion_temperature: Path | str
    toroidal_rotation: Path | str
    metadata: BalanceMetadata | Mapping[str, object] | Path | str
    equilibrium_provenance: str
    config: SimulationConfig | Path | str
    oracle: Path | str
    destination: Path | str
    # Deprecated compatibility input; the equilibrium calculation remains authoritative.
    major_radius_cm: float | None = None
    q_operation: Literal["preserve"] = "preserve"
    domains: Mapping[str, Sequence[float]] | None = None
    relative_floors: Mapping[str, float] | None = None
    interpolation_direction: str | None = None
    interpolation_method: str | None = None
    equilibrium_file: Path | str | None = None
    original_equilibrium: Path | str | None = None
    equilibrium_executable: Path | str | None = None
    equilibrium_inputs: Sequence[Path | str] = ()
    equilibrium_timeout_seconds: float = 3600.0
    tolerances: Mapping[str, object] | None = None
    expected_sha256: Mapping[str, str] | None = None
    metadata_path: Path | str | None = None
    config_path: Path | str | None = None
    metadata_snapshot: bytes | None = None
    config_snapshot: bytes | None = None
    metadata_sha256: str | None = None
    config_sha256: str | None = None
    equilibrium_parameters_file: Path | str | None = None


@dataclass(frozen=True)
class BalanceCharacterization:
    """Machine-readable outcome of one characterization request."""

    status: BalanceStatus
    report: Mapping[str, object]
    report_path: Path | None = None

    @property
    def overall_pass(self) -> bool | None:
        value = self.report.get("overall_pass")
        return value if isinstance(value, bool) else None


def characterize_balance(
    request: BalanceCharacterizationRequest | None = None, **kwargs: object
) -> BalanceCharacterization:
    """Run the explicit BALANCE characterization workflow.

    Errors are represented as ``UNAVAILABLE`` rather than converted into a
    successful-looking empty result.  The CLI turns that status into a
    non-zero exit code while keeping its JSON output parseable.
    """

    if request is None:
        try:
            request = BalanceCharacterizationRequest(**kwargs)  # type: ignore[arg-type]
        except Exception as error:
            return _unavailable("invalid characterization request", error)
    elif kwargs:
        return _unavailable("request and keyword arguments cannot be combined")

    output_root: Path | None = None
    parent_fd: int | None = None
    container_name: str | None = None
    container_fd: int | None = None
    payload_fd: int | None = None
    payload_inode: tuple[int, int] | None = None
    published = False
    final_destination: Path | None = None
    owns_output = False
    try:
        normalized = _normalise_request(request)
        final_destination = normalized["destination"]
        actual_hashes = _preflight(normalized)
        parent_fd = _open_directory_anchor(final_destination.parent)
        _validate_parent_anchor(final_destination.parent, parent_fd)
        _ensure_new_destination(
            final_destination,
            normalized["protected_roots"],
            normalized["source_paths"],
            parent_fd,
        )
        (
            container_name,
            container_fd,
            payload_fd,
            output_root,
            payload_inode,
        ) = _create_private_output(parent_fd, final_destination.name)
        output_root = _validate_private_entry(container_fd, _PAYLOAD_NAME, payload_fd)
        report_private_root = output_root.absolute()
        owns_output = True
        _validate_parent_anchor(final_destination.parent, parent_fd)
        output_root = _validate_private_entry(container_fd, _PAYLOAD_NAME, payload_fd)
        snapshots = _snapshot_nonprofile_inputs(normalized, output_root, actual_hashes)

        calculation = None
        if snapshots["equilibrium_file"] is None:
            calculation = run_equilibrium_calculation(
                snapshots["equilibrium_executable"],
                snapshots["equilibrium_inputs"],
                output_root / "equilibrium-calculation",
                timeout_seconds=normalized["equilibrium_timeout_seconds"],
            )
            equilibrium_file = calculation.equilibrium_file
            equilibrium_parameters_file = calculation.parameters_file
            equilibrium_sha256 = calculation.equilibrium_sha256
        else:
            equilibrium_file = snapshots["equilibrium_file"]
            equilibrium_parameters_file = snapshots["equilibrium_parameters_file"]
            equilibrium_sha256 = _hash_file(equilibrium_file, "equilibrium table")
        equilibrium_parameters = read_equilibrium_parameters(equilibrium_parameters_file)
        legacy_major_radius = normalized["major_radius_cm"]
        if (
            legacy_major_radius is not None
            and legacy_major_radius != equilibrium_parameters.r_big_cm
        ):
            raise ExperimentalInputError(
                "deprecated major_radius_cm does not match the equilibrium calculation"
            )
        normalized["config"] = simulation_config_for_equilibrium(
            normalized["config"], equilibrium_parameters
        )
        equilibrium_calculation = {
            "equilibrium_file": equilibrium_file,
            "equilibrium_sha256": equilibrium_sha256,
            "parameters_file": equilibrium_parameters_file,
            "parameters_sha256": equilibrium_parameters.sha256,
            "btor_gauss": equilibrium_parameters.btor_gauss,
            "r_big_cm": equilibrium_parameters.r_big_cm,
            "generator": calculation.generator_provenance if calculation is not None else None,
        }

        source = read_balance_profiles(
            density=normalized["density"],
            electron_temperature=normalized["electron_temperature"],
            ion_temperature=normalized["ion_temperature"],
            toroidal_rotation=normalized["toroidal_rotation"],
            metadata=normalized["metadata"],
        )
        _validate_parent_anchor(final_destination.parent, parent_fd)
        output_root = _validate_private_entry(container_fd, _PAYLOAD_NAME, payload_fd)
        staged = stage_balance_marsf_quartet(
            source,
            output_root / "staged-marsf",
            equilibrium_parameters_file=equilibrium_parameters_file,
            equilibrium_provenance=normalized["equilibrium_provenance"],
        )
        _validate_parent_anchor(final_destination.parent, parent_fd)
        output_root = _validate_private_entry(container_fd, _PAYLOAD_NAME, payload_fd)
        marsf = read_marsf_profiles(staged.directory, staged.metadata)
        actual_hashes = _verify_unchanged_inputs(normalized, actual_hashes)

        _validate_parent_anchor(final_destination.parent, parent_fd)
        output_root = _validate_private_entry(container_fd, _PAYLOAD_NAME, payload_fd)
        prepared = prepare_marsf_case(
            marsf,
            normalized["config"],
            output_root / "prepared",
            equilibrium_file=equilibrium_file,
            equilibrium_parameters_file=equilibrium_parameters_file,
            equilibrium_calculation=calculation,
            q_operation=normalized["q_operation"],
            upstream_staging_report=staged.report,
        )
        _validate_parent_anchor(final_destination.parent, parent_fd)
        output_root = _validate_private_entry(container_fd, _PAYLOAD_NAME, payload_fd)

        # The oracle is intentionally opened only after preparation.  Its
        # arrays never enter staging or KIM profile generation.
        oracle = read_ql_balance_oracle(snapshots["oracle"])
        prepared_profile_files = _prepared_profile_files(prepared)
        comparisons = _compare_all(normalized, prepared, oracle, prepared_profile_files)
        status, threshold_decisions, overall_pass = _status_for(comparisons)
        prepared_profile_hashes = {
            role: _hash_file(prepared.profiles / filename, "prepared profile")
            for role, filename in prepared_profile_files.items()
        }

        actual_hashes = _verify_unchanged_inputs(normalized, actual_hashes)
        _validate_parent_anchor(final_destination.parent, parent_fd)
        output_root = _validate_private_entry(container_fd, _PAYLOAD_NAME, payload_fd)
        report_private_root = output_root.absolute()
        _rewrite_report_paths(prepared.report, report_private_root, final_destination)
        report = _success_report(
            normalized,
            actual_hashes,
            staged,
            prepared,
            comparisons,
            status,
            threshold_decisions,
            overall_pass,
            prepared_profile_files,
            prepared_profile_hashes,
            snapshots=snapshots,
            equilibrium_calculation=equilibrium_calculation,
            private_root=report_private_root,
            final_root=final_destination,
            retained_staging=_retained_staging_info(container_fd, "empty_container"),
        )
        report_path = output_root / "characterization_report.json"
        report_path.write_text(_json_dump(report), encoding="utf-8")
        _fsync_file_and_directory(report_path, output_root)
        actual_hashes = _verify_unchanged_inputs(normalized, actual_hashes)
        _validate_parent_anchor(final_destination.parent, parent_fd)
        output_root = _validate_private_entry(container_fd, _PAYLOAD_NAME, payload_fd)
        assert container_fd is not None and parent_fd is not None and payload_inode is not None
        _publish_staging_directory_at(
            container_fd,
            _PAYLOAD_NAME,
            parent_fd,
            final_destination.name,
            expected_source_inode=payload_inode,
        )
        published = True
        owns_output = False
        return BalanceCharacterization(
            status,
            report,
            final_destination / "characterization_report.json",
        )
    except (MemoryError, KeyboardInterrupt):
        if not published and owns_output and container_fd is not None and payload_fd is not None:
            _retain_private_container(container_fd, payload_fd)
        raise
    except Exception as error:
        retained_staging = (
            _retained_staging_info(container_fd, "incomplete_container")
            if container_fd is not None
            else None
        )
        if not published and owns_output and container_fd is not None and payload_fd is not None:
            _retain_private_container(container_fd, payload_fd)
        return _unavailable(
            "characterization unavailable",
            error,
            retained_staging=retained_staging,
        )
    finally:
        if payload_fd is not None:
            try:
                os.close(payload_fd)
            except OSError:
                pass
        if container_fd is not None:
            try:
                os.close(container_fd)
            except OSError:
                pass
        if parent_fd is not None:
            os.close(parent_fd)


def _normalise_request(request: BalanceCharacterizationRequest) -> dict[str, Any]:
    if not isinstance(request, BalanceCharacterizationRequest):
        raise ExperimentalInputError("request must be a BalanceCharacterizationRequest")
    if request.q_operation != "preserve":
        raise ExperimentalInputError(
            "BALANCE characterization preserves Fouriers q without a sign change"
        )
    if not isinstance(request.interpolation_direction, str) or not request.interpolation_direction:
        raise ExperimentalInputError("interpolation_direction must be a non-empty string")
    if not isinstance(request.interpolation_method, str) or not request.interpolation_method:
        raise ExperimentalInputError("interpolation_method must be a non-empty string")
    _reject_conflicting_reference(
        request.metadata_path,
        request.metadata,
        "metadata",
    )
    _reject_conflicting_reference(request.config_path, request.config, "config")
    metadata_reference = (
        request.metadata_path
        if request.metadata_path is not None
        else request.metadata if isinstance(request.metadata, (Path, str)) else None
    )
    config_reference = (
        request.config_path
        if request.config_path is not None
        else request.config if isinstance(request.config, (Path, str)) else None
    )
    metadata_alias = _absolute_alias(metadata_reference) if metadata_reference is not None else None
    config_alias = _absolute_alias(config_reference) if config_reference is not None else None
    metadata_source = (
        _resolve_source_file(metadata_alias, "BALANCE metadata")
        if metadata_reference is not None
        else None
    )
    config_source = (
        _resolve_source_file(config_alias, "KIM configuration")
        if config_reference is not None
        else None
    )
    metadata_bytes, metadata_digest = _source_snapshot(
        metadata_source,
        request.metadata_snapshot,
        request.metadata_sha256,
        "BALANCE metadata",
    )
    config_bytes, config_digest = _source_snapshot(
        config_source,
        request.config_snapshot,
        request.config_sha256,
        "KIM configuration",
    )
    metadata = _load_metadata(request.metadata, metadata_bytes, metadata_source)
    config = _load_config(request.config, config_bytes, config_source)
    destination = _resolve_destination(request.destination)
    if not request.equilibrium_provenance.strip():
        raise ExperimentalInputError("equilibrium_provenance must be nonempty")
    if request.equilibrium_file is not None and (
        request.original_equilibrium is not None or request.equilibrium_executable is not None
    ):
        raise ExperimentalInputError(
            "equilibrium_file is mutually exclusive with "
            "original_equilibrium/equilibrium_executable"
        )
    if request.equilibrium_file is not None and request.equilibrium_parameters_file is None:
        raise ExperimentalInputError("equilibrium_file requires its paired btor_rbig.dat output")
    if request.equilibrium_file is None and (
        request.original_equilibrium is None or request.equilibrium_executable is None
    ):
        raise ExperimentalInputError(
            "provide equilibrium_file, or both original_equilibrium and equilibrium_executable"
        )
    if request.equilibrium_file is None and request.equilibrium_parameters_file is not None:
        raise ExperimentalInputError(
            "btor_rbig.dat must be generated with the equilibrium calculation"
        )
    if request.equilibrium_file is not None and request.equilibrium_inputs:
        raise ExperimentalInputError("equilibrium_inputs require equilibrium_executable")
    timeout = _strict_real(request.equilibrium_timeout_seconds, "equilibrium_timeout_seconds")
    if not np.isfinite(timeout) or timeout <= 0.0:
        raise ExperimentalInputError("equilibrium_timeout_seconds must be finite and positive")
    legacy_major_radius = None
    if request.major_radius_cm is not None:
        legacy_major_radius = _strict_real(request.major_radius_cm, "major_radius_cm")
        if not np.isfinite(legacy_major_radius) or legacy_major_radius <= 0.0:
            raise ExperimentalInputError("major_radius_cm must be finite and positive")
    domains = _normalise_domains(request.domains)
    floors = _normalise_floors(request.relative_floors)
    tolerances = _normalise_tolerances(request.tolerances, domains)
    if (
        floors is None
        and tolerances is not None
        and any(
            metric.startswith("relative_")
            for metric_mapping in tolerances.values()
            for metric in metric_mapping
        )
    ):
        raise ExperimentalInputError("relative tolerances require relative_floors")
    balance_aliases = {
        "density": _absolute_alias(request.density),
        "electron_temperature": _absolute_alias(request.electron_temperature),
        "ion_temperature": _absolute_alias(request.ion_temperature),
        "toroidal_rotation": _absolute_alias(request.toroidal_rotation),
    }
    balance_paths = {
        "density": _resolve_source_file(balance_aliases["density"], "BALANCE density"),
        "electron_temperature": _resolve_source_file(
            balance_aliases["electron_temperature"], "BALANCE electron temperature"
        ),
        "ion_temperature": _resolve_source_file(
            balance_aliases["ion_temperature"], "BALANCE ion temperature"
        ),
        "toroidal_rotation": _resolve_source_file(
            balance_aliases["toroidal_rotation"], "BALANCE toroidal rotation"
        ),
    }
    original_alias = (
        _absolute_alias(request.original_equilibrium)
        if request.original_equilibrium is not None
        else None
    )
    equilibrium_alias = (
        _absolute_alias(request.equilibrium_file) if request.equilibrium_file is not None else None
    )
    equilibrium_parameters_alias = (
        _absolute_alias(request.equilibrium_parameters_file)
        if request.equilibrium_parameters_file is not None
        else None
    )
    original_equilibrium = (
        _resolve_source_file(original_alias, "original equilibrium")
        if original_alias is not None
        else None
    )
    equilibrium_file = (
        _resolve_source_file(equilibrium_alias, "reduced equilibrium")
        if equilibrium_alias is not None
        else None
    )
    equilibrium_parameters_file = (
        _resolve_source_file(equilibrium_parameters_alias, "equilibrium btor_rbig.dat")
        if equilibrium_parameters_alias is not None
        else None
    )
    if (
        equilibrium_parameters_file is not None
        and equilibrium_parameters_file.name.casefold() != "btor_rbig.dat"
    ):
        raise ExperimentalInputError("equilibrium parameters file must be named btor_rbig.dat")
    if equilibrium_file is not None:
        if equilibrium_file.name.casefold() != "equil_r_q_psi.dat":
            raise ExperimentalInputError("equilibrium table file must be named equil_r_q_psi.dat")
        if equilibrium_parameters_file is None or (
            equilibrium_file.parent != equilibrium_parameters_file.parent
        ):
            raise ExperimentalInputError(
                "precomputed equilibrium outputs must share the same "
                "equilibrium calculation directory"
            )
    equilibrium_input_aliases = tuple(_absolute_alias(item) for item in request.equilibrium_inputs)
    equilibrium_input_paths = tuple(
        _resolve_source_file(item, "equilibrium input") for item in equilibrium_input_aliases
    )
    all_input_paths = tuple(
        ([original_equilibrium] if original_equilibrium is not None else [])
        + list(equilibrium_input_paths)
    )
    input_basenames = [item.name.casefold() for item in all_input_paths]
    if len(input_basenames) != len(set(input_basenames)):
        raise ExperimentalInputError(
            "equilibrium input basenames must be unique case-insensitively"
        )
    equilibrium_inputs = all_input_paths
    equilibrium_executable = (
        _resolve_source_file(
            _absolute_alias(request.equilibrium_executable), "equilibrium executable"
        )
        if request.equilibrium_executable is not None
        else None
    )
    equilibrium_executable_alias = (
        _absolute_alias(request.equilibrium_executable)
        if request.equilibrium_executable is not None
        else None
    )
    oracle_alias = _absolute_alias(request.oracle)
    oracle = _resolve_source_file(oracle_alias, "QL-Balance oracle")
    source_paths = tuple(
        [
            *balance_paths.values(),
            *(
                path
                for path in (
                    original_equilibrium,
                    equilibrium_file,
                    equilibrium_parameters_file,
                )
                if path
            ),
        ]
        + list(equilibrium_input_paths)
        + ([equilibrium_executable] if equilibrium_executable is not None else [])
        + [oracle]
        + ([metadata_source] if metadata_source is not None else [])
        + ([config_source] if config_source is not None else [])
    )
    # Every explicit source file protects its resolved parent directory.  In
    # particular, this prevents an output tree from becoming an alias inside
    # an equilibrium/oracle/preprocessor directory as well as a BALANCE tree.
    protected_roots = tuple(sorted({path.parent for path in source_paths}, key=str))
    snapshot_hashes = {
        key: digest
        for key, digest in (("metadata", metadata_digest), ("config", config_digest))
        if digest is not None
    }
    return {
        "metadata": metadata,
        "metadata_source": metadata_source,
        "metadata_alias": metadata_alias,
        "config": config,
        "config_source": config_source,
        "config_alias": config_alias,
        "snapshot_hashes": snapshot_hashes,
        "destination": destination,
        "balance_paths": balance_paths,
        "balance_aliases": balance_aliases,
        **balance_paths,
        "original_equilibrium": original_equilibrium,
        "original_equilibrium_alias": original_alias,
        "equilibrium_file": equilibrium_file,
        "equilibrium_file_alias": equilibrium_alias,
        "equilibrium_parameters_file": equilibrium_parameters_file,
        "equilibrium_parameters_alias": equilibrium_parameters_alias,
        "equilibrium_executable": equilibrium_executable,
        "equilibrium_executable_alias": equilibrium_executable_alias,
        "equilibrium_input_aliases": equilibrium_input_aliases,
        "equilibrium_input_paths": equilibrium_input_paths,
        "equilibrium_inputs": equilibrium_inputs,
        "oracle": oracle,
        "oracle_alias": oracle_alias,
        "equilibrium_provenance": request.equilibrium_provenance.strip(),
        "major_radius_cm": legacy_major_radius,
        "q_operation": request.q_operation,
        "domains": domains,
        "relative_floors": floors,
        "interpolation_direction": request.interpolation_direction,
        "interpolation_method": request.interpolation_method,
        "tolerances": tolerances,
        "tolerances_input": _json_value(request.tolerances),
        "expected_sha256": _normalise_expected(
            request.expected_sha256,
            equilibrium_executable=equilibrium_executable,
            equilibrium_input_paths=equilibrium_input_paths,
            equilibrium_parameters_file=equilibrium_parameters_file,
            metadata_source=metadata_source,
            config_source=config_source,
        ),
        "protected_roots": protected_roots,
        "source_paths": source_paths,
        "equilibrium_timeout_seconds": timeout,
    }


def _reject_conflicting_reference(
    explicit_path: Path | str | None,
    value: object,
    label: str,
) -> None:
    if explicit_path is None or not isinstance(value, (Path, str)):
        return
    explicit = _resolve_source_file(_absolute_alias(explicit_path), label)
    value_path = _resolve_source_file(_absolute_alias(value), label)
    if explicit != value_path:
        raise ExperimentalInputError(
            f"{label} path and {label} reference resolve to different files"
        )


def _source_snapshot(
    source: Path | None,
    snapshot: bytes | None,
    declared_digest: str | None,
    label: str,
) -> tuple[bytes | None, str | None]:
    if source is None:
        if snapshot is not None or declared_digest is not None:
            raise ExperimentalInputError(f"{label} snapshot requires a file path")
        return None, None
    if snapshot is None:
        try:
            snapshot = source.read_bytes()
        except (OSError, ExperimentalInputError) as error:
            raise ExperimentalInputError(f"unable to read {label}: {source}") from error
    elif not isinstance(snapshot, bytes):
        raise ExperimentalInputError(f"{label} snapshot must be bytes")
    digest = hashlib.sha256(snapshot).hexdigest()
    if declared_digest is not None:
        _validate_digest(declared_digest, f"{label} snapshot")
        if declared_digest != digest:
            raise ExperimentalInputError(f"{label} snapshot digest does not match its bytes")
    return snapshot, digest


def _validate_digest(value: str, label: str) -> None:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(char not in "0123456789abcdef" for char in value)
    ):
        raise ExperimentalInputError(f"{label} is not a SHA-256 digest")


def _strict_json_payload(payload: bytes, label: str) -> object:
    try:
        return json.loads(
            payload,
            object_pairs_hook=_reject_duplicate_json_keys,
            parse_constant=lambda token: (_ for _ in ()).throw(
                ValueError(f"JSON constant {token!r} is not allowed")
            ),
        )
    except (TypeError, ValueError, json.JSONDecodeError) as error:
        raise ExperimentalInputError(f"{label} must be valid strict JSON") from error


def _reject_duplicate_json_keys(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON object key: {key!r}")
        result[key] = value
    return result


def _load_metadata(
    value: BalanceMetadata | Mapping[str, object] | Path | str,
    snapshot: bytes | None = None,
    source: Path | None = None,
) -> BalanceMetadata:
    try:
        if snapshot is not None:
            parsed = BalanceMetadata.model_validate(_strict_json_payload(snapshot, "metadata"))
            if isinstance(value, BalanceMetadata) and parsed != value:
                raise ExperimentalInputError("BALANCE metadata object differs from file snapshot")
            if isinstance(value, Mapping) and parsed != BalanceMetadata.model_validate(value):
                raise ExperimentalInputError("BALANCE metadata mapping differs from file snapshot")
            return parsed
        if isinstance(value, BalanceMetadata):
            return value
        if isinstance(value, (Path, str)):
            path = source or Path(value)
            return BalanceMetadata.model_validate(
                _strict_json_payload(path.read_bytes(), "metadata")
            )
        return BalanceMetadata.model_validate(value)
    except Exception as error:
        if isinstance(error, ExperimentalInputError):
            raise
        raise ExperimentalInputError("unable to load BALANCE metadata") from error


def _load_config(
    value: SimulationConfig | Path | str,
    snapshot: bytes | None = None,
    source: Path | None = None,
) -> SimulationConfig:
    try:
        path = source or Path(value) if isinstance(value, (Path, str)) else source
        if snapshot is not None:
            if path is None:
                raise ExperimentalInputError("configuration snapshot requires a file path")
            if path.suffix.casefold() == ".json":
                parsed = SimulationConfig.model_validate(
                    _strict_json_payload(snapshot, "configuration")
                )
            else:
                parsed = load_namelist_bytes(snapshot, str(path))
            if isinstance(value, SimulationConfig) and parsed != value:
                raise ExperimentalInputError("KIM configuration object differs from file snapshot")
            return parsed
        if isinstance(value, SimulationConfig):
            return value
        assert path is not None
        if path.suffix.casefold() == ".json":
            return SimulationConfig.model_validate(
                _strict_json_payload(path.read_bytes(), "configuration")
            )
        return SimulationConfig.from_namelist(path)
    except Exception as error:
        if isinstance(error, ExperimentalInputError):
            raise
        raise ExperimentalInputError(f"unable to load KIM configuration: {path}") from error


def _resolve_source_file(value: Path | str, label: str) -> Path:
    candidate = Path(value).expanduser()
    try:
        resolved = candidate.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise ExperimentalInputError(f"{candidate}: required {label} is missing") from error
    if not resolved.is_file():
        raise ExperimentalInputError(f"{candidate}: required {label} is not a regular file")
    return resolved


def _absolute_alias(value: Path | str) -> Path:
    return Path(value).expanduser().absolute()


def _resolve_destination(value: Path | str) -> Path:
    candidate = Path(value).expanduser().absolute()
    _reject_dangling_symlink_components(candidate)
    if candidate.is_symlink():
        raise ExperimentalInputError(
            f"characterization destination must not be a symlink: {candidate}"
        )
    try:
        # ``resolve(strict=False)`` resolves every existing parent (including
        # symlinks) and appends the lexical remainder, while also normalizing
        # ``..`` through that parent.  It does not create the destination.
        resolved = candidate.resolve(strict=False)
    except (OSError, RuntimeError) as error:
        raise ExperimentalInputError(
            f"unable to resolve characterization destination: {candidate}"
        ) from error
    if not resolved.parent.exists():
        raise ExperimentalInputError(
            f"characterization destination parent must already exist: {resolved.parent}"
        )
    if not resolved.parent.is_dir():
        raise ExperimentalInputError(
            f"characterization destination parent is not a directory: {resolved.parent}"
        )
    return resolved


def _reject_dangling_symlink_components(candidate: Path) -> None:
    current = Path(candidate.anchor)
    for component in candidate.parts[1:]:
        current /= component
        if current.is_symlink() and not current.exists():
            raise ExperimentalInputError(
                f"characterization destination contains a dangling symlink: {current}"
            )


def _open_directory_anchor(path: Path) -> int:
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_CLOEXEC", 0)
    nofollow = getattr(os, "O_NOFOLLOW", 0)
    descriptor: int | None = None
    try:
        descriptor = os.open(path.anchor, flags | nofollow)
        for component in path.parts[1:]:
            next_descriptor = os.open(component, flags | nofollow, dir_fd=descriptor)
            old_descriptor = descriptor
            descriptor = None
            try:
                os.close(old_descriptor)
            except OSError:
                os.close(next_descriptor)
                raise
            descriptor = next_descriptor
        assert descriptor is not None
        return descriptor
    except OSError as error:
        if descriptor is not None:
            try:
                os.close(descriptor)
            except OSError:
                pass
        raise ExperimentalInputError(f"unable to anchor destination parent: {path}") from error


def _validate_parent_anchor(path: Path, descriptor: int) -> None:
    try:
        path_stat = os.stat(path, follow_symlinks=False)
        descriptor_stat = os.fstat(descriptor)
    except OSError as error:
        raise ExperimentalInputError(f"destination parent changed: {path}") from error
    if (path_stat.st_dev, path_stat.st_ino) != (
        descriptor_stat.st_dev,
        descriptor_stat.st_ino,
    ):
        raise ExperimentalInputError(f"destination parent changed: {path}")


def _descriptor_path(descriptor: int) -> Path:
    """Return a path rooted at an open directory descriptor.

    Linux exposes descriptors as traversable ``/proc/self/fd`` paths.  macOS
    exposes ``/dev/fd`` for descriptor operations, but does not make that path
    traversable in all supported Python/runtime combinations; in that case
    ``F_GETPATH`` is the platform-native descriptor anchor fallback.  The
    descriptor remains open and callers validate its inode before publication.
    """

    if sys.platform.startswith("linux"):
        candidate = Path(f"/proc/self/fd/{descriptor}")
    elif sys.platform == "darwin":
        candidate = Path(f"/dev/fd/{descriptor}")
    else:
        raise ExperimentalInputError(
            "private characterization paths require a supported descriptor namespace"
        )
    try:
        candidate_stat = os.stat(candidate, follow_symlinks=True)
        descriptor_stat = os.fstat(descriptor)
        if not stat.S_ISDIR(candidate_stat.st_mode) or (
            candidate_stat.st_dev,
            candidate_stat.st_ino,
        ) != (descriptor_stat.st_dev, descriptor_stat.st_ino):
            raise OSError("descriptor namespace does not identify the open directory")
        # Validate that Path operations can actually traverse this namespace.
        os.listdir(candidate)
        return candidate
    except OSError:
        if sys.platform != "darwin":
            raise ExperimentalInputError(
                f"descriptor path is unavailable for directory descriptor {descriptor}"
            )
        try:
            import fcntl

            raw = fcntl.fcntl(descriptor, fcntl.F_GETPATH, b"\0" * 1024)
            fallback = Path(os.fsdecode(raw).split("\0", 1)[0])
            fallback_stat = os.stat(fallback, follow_symlinks=False)
            descriptor_stat = os.fstat(descriptor)
            if (fallback_stat.st_dev, fallback_stat.st_ino) != (
                descriptor_stat.st_dev,
                descriptor_stat.st_ino,
            ):
                raise OSError("descriptor path changed")
            return fallback
        except (ImportError, OSError, ValueError) as error:
            raise ExperimentalInputError(
                f"descriptor path is unavailable for directory descriptor {descriptor}"
            ) from error


def _canonical_descriptor_path(descriptor: int) -> Path:
    try:
        return _descriptor_path(descriptor).resolve(strict=True)
    except (OSError, RuntimeError, ExperimentalInputError) as error:
        raise ExperimentalInputError(
            f"unable to resolve retained staging container for descriptor {descriptor}"
        ) from error


def _retained_staging_info(container_fd: int | None, state: str) -> dict[str, object]:
    if container_fd is None:
        return {
            "path": None,
            "reason": "private staging path unavailable: no retained descriptor",
            "state": state,
            "cleanup_safe_when_no_characterization_is_running": True,
        }
    try:
        path = _canonical_descriptor_path(container_fd)
    except ExperimentalInputError as error:
        return {
            "path": None,
            "reason": str(error),
            "state": state,
            "cleanup_safe_when_no_characterization_is_running": True,
        }
    return {
        "path": str(path),
        "reason": _RETAINED_STAGING_REASON,
        "state": state,
        "cleanup_safe_when_no_characterization_is_running": True,
    }


def _validate_private_entry(parent_fd: int, name: str, private_fd: int) -> Path:
    """Validate a fixed private child and return its descriptor-derived path."""

    try:
        entry_stat = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        private_stat = os.fstat(private_fd)
    except OSError as error:
        raise ExperimentalInputError(f"private staging entry changed: {name}") from error
    if not stat.S_ISDIR(entry_stat.st_mode) or (entry_stat.st_dev, entry_stat.st_ino) != (
        private_stat.st_dev,
        private_stat.st_ino,
    ):
        raise ExperimentalInputError(f"private staging entry changed: {name}")
    return _descriptor_path(private_fd)


def _create_private_output(
    parent_fd: int, destination_name: str
) -> tuple[str, int, int, Path, tuple[int, int]]:
    for _ in range(128):
        container_name = f".{destination_name}-{secrets.token_hex(8)}"
        try:
            os.mkdir(container_name, mode=0o700, dir_fd=parent_fd)
        except FileExistsError:
            continue
        container_fd: int | None = None
        payload_fd: int | None = None
        try:
            flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_CLOEXEC", 0)
            flags |= getattr(os, "O_NOFOLLOW", 0)
            container_fd = os.open(container_name, flags, dir_fd=parent_fd)
            os.mkdir(_PAYLOAD_NAME, mode=0o700, dir_fd=container_fd)
            payload_fd = os.open(_PAYLOAD_NAME, flags, dir_fd=container_fd)
            payload_root = _validate_private_entry(container_fd, _PAYLOAD_NAME, payload_fd)
            payload_stat = os.fstat(payload_fd)
        except (OSError, ExperimentalInputError) as error:
            if payload_fd is not None:
                try:
                    os.close(payload_fd)
                except OSError:
                    pass
            if container_fd is not None:
                try:
                    os.close(container_fd)
                except OSError:
                    pass
            raise ExperimentalInputError(
                "unable to anchor private characterization container"
            ) from error
        assert container_fd is not None and payload_fd is not None
        return (
            container_name,
            container_fd,
            payload_fd,
            payload_root,
            (payload_stat.st_dev, payload_stat.st_ino),
        )
    raise ExperimentalInputError("unable to reserve a private characterization container")


def _retain_private_container(container_fd: int, payload_fd: int) -> None:
    """Intentionally retain incomplete private staging; descriptors close in ``finally``."""

    # The container is deliberately not recursively deleted.  Its private
    # mode-0700 name remains the ownership boundary for any partial files.
    del container_fd, payload_fd


def _fsync_file_and_directory(report_path: Path, directory: Path) -> None:
    try:
        with report_path.open("rb") as report_file:
            os.fsync(report_file.fileno())
        directory_fd = os.open(
            directory,
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_CLOEXEC", 0),
        )
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except OSError:
        # Directory fsync is not available on every supported filesystem.  The
        # report write itself is complete before publication; continue when the
        # platform cannot durably fsync the directory.
        return


def _strict_real(value: object, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, (int, float, np.integer, np.floating)
    ):
        raise ExperimentalInputError(f"{name} must be a real int or float")
    try:
        return float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ExperimentalInputError(f"{name} must be a real int or float") from error


def _normalise_domains(domains: Mapping[str, Sequence[float]]) -> dict[str, tuple[float, float]]:
    if not isinstance(domains, Mapping) or not domains:
        raise ExperimentalInputError("domains must be a non-empty named mapping")
    result: dict[str, tuple[float, float]] = {}
    for name, endpoints in domains.items():
        if not isinstance(name, str) or not name.strip():
            raise ExperimentalInputError("domain names must be non-empty strings")
        try:
            if len(endpoints) != 2:
                raise ValueError
            lower = _strict_real(endpoints[0], f"domain {name!r} endpoint")
            upper = _strict_real(endpoints[1], f"domain {name!r} endpoint")
        except ExperimentalInputError:
            raise
        except (TypeError, ValueError, OverflowError, IndexError) as error:
            raise ExperimentalInputError(
                f"domain {name!r} must contain two finite endpoints"
            ) from error
        if not np.isfinite(lower) or not np.isfinite(upper):
            raise ExperimentalInputError(f"domain {name!r} endpoints must be finite")
        if lower >= upper:
            raise ExperimentalInputError(f"domain {name!r} must be non-empty with lower < upper")
        if name in result:
            raise ExperimentalInputError(f"duplicate domain name: {name}")
        result[name] = (lower, upper)
    return result


def _normalise_floors(floors: Mapping[str, float] | None) -> dict[str, float] | None:
    if floors is None:
        return None
    if not isinstance(floors, Mapping):
        raise ExperimentalInputError("relative_floors must be a mapping")
    result: dict[str, float] = {}
    for key, value in floors.items():
        if key not in _ROLE_ALIASES:
            raise ExperimentalInputError(f"relative_floors contains unknown profile: {key!r}")
        role = _ROLE_ALIASES[key]
        if role in result:
            raise ExperimentalInputError(f"relative_floors specifies {role} more than once")
        try:
            floor = _strict_real(value, f"relative floor for {role}")
        except ExperimentalInputError as error:
            raise ExperimentalInputError(f"relative floor for {role} is invalid") from error
        if not np.isfinite(floor) or floor <= 0.0:
            raise ExperimentalInputError(f"relative floor for {role} must be finite and positive")
        result[role] = floor
    missing = tuple(role for role in _PROFILE_ROLES if role not in result)
    if missing:
        raise ExperimentalInputError(f"relative_floors is missing profiles: {', '.join(missing)}")
    return result


def _normalise_tolerances(
    tolerances: Mapping[str, object] | None, domains: Mapping[str, tuple[float, float]]
) -> dict[str, dict[str, float]] | None:
    if tolerances is None:
        return None
    if not isinstance(tolerances, Mapping) or not tolerances:
        raise ExperimentalInputError("tolerances must be a non-empty mapping when supplied")
    # Accepted JSON forms are flat metrics, profile -> metrics, or domain ->
    # (metrics/profile -> metrics).  Every numeric value remains caller data.
    if all(key in _METRICS for key in tolerances):
        flat = _metrics(tolerances)
        return {f"{domain}/{role}": flat for domain in domains for role in _PROFILE_ROLES}
    result: dict[str, dict[str, float]] = {}
    keys = set(tolerances)
    if keys.issubset(set(_PROFILE_ROLES) | set(_ROLE_ALIASES)):
        by_role: dict[str, dict[str, float]] = {}
        for key, value in tolerances.items():
            role = _ROLE_ALIASES[key]
            if role in by_role:
                raise ExperimentalInputError(f"tolerances specifies {role} more than once")
            by_role[role] = _metrics(value)
        for domain in domains:
            for role in _PROFILE_ROLES:
                if role in by_role:
                    result[f"{domain}/{role}"] = by_role[role]
        return result
    if not keys.issubset(domains):
        raise ExperimentalInputError("tolerances must use metric, profile, or named-domain keys")
    for domain, value in tolerances.items():
        if not isinstance(value, Mapping) or not value:
            raise ExperimentalInputError(f"tolerances for domain {domain!r} must be a mapping")
        value_keys = set(value)
        if value_keys.issubset(_METRICS):
            metrics = _metrics(value)
            for role in _PROFILE_ROLES:
                result[f"{domain}/{role}"] = metrics
        elif value_keys.issubset(set(_PROFILE_ROLES) | set(_ROLE_ALIASES)):
            for role_key, role_value in value.items():
                role = _ROLE_ALIASES[role_key]
                target = f"{domain}/{role}"
                if target in result:
                    raise ExperimentalInputError(f"tolerances specifies {role} more than once")
                result[target] = _metrics(role_value)
        else:
            raise ExperimentalInputError(f"tolerances for domain {domain!r} are invalid")
    return result


def _metrics(value: object) -> dict[str, float]:
    if not isinstance(value, Mapping) or not value:
        raise ExperimentalInputError("each tolerance value must be a non-empty metric mapping")
    if any(key not in _METRICS for key in value):
        raise ExperimentalInputError("tolerances contain an unknown metric")
    result: dict[str, float] = {}
    for key, item in value.items():
        try:
            limit = _strict_real(item, f"tolerance for {key}")
        except ExperimentalInputError as error:
            raise ExperimentalInputError(f"tolerance for {key} is invalid") from error
        if not np.isfinite(limit) or limit < 0.0:
            raise ExperimentalInputError(f"tolerance for {key} must be finite and non-negative")
        result[str(key)] = limit
    return result


def _normalise_expected(
    expected: Mapping[str, str] | None,
    *,
    equilibrium_executable: Path | None,
    equilibrium_input_paths: Sequence[Path],
    equilibrium_parameters_file: Path | None,
    metadata_source: Path | None,
    config_source: Path | None,
) -> dict[str, str] | None:
    if expected is None:
        return None
    if not isinstance(expected, Mapping):
        raise ExperimentalInputError("expected_sha256 must be a mapping")
    aliases = {
        "density": "density",
        "electron_temperature": "electron_temperature",
        "ion_temperature": "ion_temperature",
        "toroidal_rotation": "toroidal_rotation",
        "equilibrium": "equilibrium",
        "equilibrium_parameters": "equilibrium_parameters",
        "equilibrium_executable": "equilibrium_executable",
        "oracle": "oracle",
        "metadata": "metadata",
        "config": "config",
    }
    input_keys = {
        f"equilibrium_input:{path.name.casefold()}": f"equilibrium_input:{path.name.casefold()}"
        for path in equilibrium_input_paths
    }
    aliases.update(input_keys)
    if equilibrium_executable is None and "equilibrium_executable" in expected:
        raise ExperimentalInputError(
            "expected_sha256 includes an unavailable equilibrium executable"
        )
    if equilibrium_parameters_file is None and "equilibrium_parameters" in expected:
        raise ExperimentalInputError("expected_sha256 includes unavailable equilibrium parameters")
    if metadata_source is None and "metadata" in expected:
        raise ExperimentalInputError("expected_sha256 includes unavailable metadata")
    if config_source is None and "config" in expected:
        raise ExperimentalInputError("expected_sha256 includes unavailable config")
    result: dict[str, str] = {}
    for key, value in expected.items():
        if key not in aliases:
            raise ExperimentalInputError(f"expected_sha256 contains unknown key: {key!r}")
        canonical = aliases[key]
        if canonical in result:
            raise ExperimentalInputError(f"expected_sha256 specifies {canonical} more than once")
        if (
            not isinstance(value, str)
            or len(value) != 64
            or any(char not in "0123456789abcdef" for char in value)
        ):
            raise ExperimentalInputError(f"expected_sha256 for {canonical} is not a SHA-256 digest")
        result[canonical] = value
    required = {
        "density",
        "electron_temperature",
        "ion_temperature",
        "toroidal_rotation",
        "equilibrium",
        "oracle",
    }
    if equilibrium_parameters_file is not None:
        required.add("equilibrium_parameters")
    if equilibrium_executable is not None:
        required.add("equilibrium_executable")
    if metadata_source is not None:
        required.add("metadata")
    if config_source is not None:
        required.add("config")
    required.update(input_keys.values())
    if set(result) != required:
        raise ExperimentalInputError(
            "expected_sha256 must contain all four BALANCE sources, equilibrium, and oracle"
        )
    return result


def _preflight(data: Mapping[str, Any]) -> dict[str, str]:
    paths = data["balance_paths"]
    actual: dict[str, str] = {
        role: _hash_alias(data["balance_aliases"][role], path, "BALANCE profile")
        for role, path in paths.items()
    }
    equilibrium_path = data["equilibrium_file"] or data["original_equilibrium"]
    assert equilibrium_path is not None
    equilibrium_alias = data["equilibrium_file_alias"] or data["original_equilibrium_alias"]
    assert equilibrium_alias is not None
    actual["equilibrium"] = _hash_alias(equilibrium_alias, equilibrium_path, "equilibrium input")
    if data["equilibrium_parameters_file"] is not None:
        actual["equilibrium_parameters"] = _hash_alias(
            data["equilibrium_parameters_alias"],
            data["equilibrium_parameters_file"],
            "equilibrium btor_rbig.dat",
        )
    actual["oracle"] = _hash_alias(data["oracle_alias"], data["oracle"], "QL-Balance oracle")
    if data["metadata_source"] is not None:
        actual["metadata"] = _hash_alias(
            data["metadata_alias"],
            data["metadata_source"],
            "BALANCE metadata",
        )
    if data["config_source"] is not None:
        actual["config"] = _hash_alias(
            data["config_alias"],
            data["config_source"],
            "KIM configuration",
        )
    if data["equilibrium_executable"] is not None:
        actual["equilibrium_executable"] = _hash_alias(
            data["equilibrium_executable_alias"],
            data["equilibrium_executable"],
            "equilibrium executable",
        )
    for alias, path in zip(
        data["equilibrium_input_aliases"], data["equilibrium_input_paths"], strict=True
    ):
        actual[f"equilibrium_input:{path.name.casefold()}"] = _hash_alias(
            alias, path, "equilibrium input"
        )
    expected = data["expected_sha256"]
    snapshot_mismatches = {
        key: {"snapshot": digest, "actual": actual[key]}
        for key, digest in data["snapshot_hashes"].items()
        if actual.get(key) != digest
    }
    if snapshot_mismatches:
        raise ExperimentalInputError(
            f"source file changed after its snapshot was read: {snapshot_mismatches}"
        )
    if expected is not None:
        mismatches = {
            key: {"expected": expected[key], "actual": actual[key]}
            for key in actual
            if expected[key] != actual[key]
        }
        if mismatches:
            raise ExperimentalInputError(f"expected SHA-256 identity mismatch: {mismatches}")
    return actual


def _snapshot_nonprofile_inputs(
    data: Mapping[str, Any], private_root: Path, verified_hashes: Mapping[str, str]
) -> dict[str, Any]:
    snapshot_root = private_root / "input-snapshots"
    snapshot_root.mkdir()
    entries: dict[str, dict[str, str]] = {}

    def copy_input(
        key: str,
        alias: Path,
        source: Path,
        relative: Path,
        label: str,
        *,
        executable: bool = False,
    ) -> Path:
        current = _resolve_source_file(alias, label)
        if current != source:
            raise ExperimentalInputError(f"{alias}: {label} changed its resolved path")
        try:
            payload = source.read_bytes()
        except OSError as error:
            raise ExperimentalInputError(f"unable to snapshot {label}: {source}") from error
        digest = hashlib.sha256(payload).hexdigest()
        if digest != verified_hashes[key]:
            raise ExperimentalInputError(f"{source}: {label} changed before snapshotting")
        destination = snapshot_root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(payload)
        if executable:
            destination.chmod(stat.S_IMODE(source.stat().st_mode))
        entries[key] = {
            "original_path": str(source),
            "snapshot_path": str(destination.absolute()),
            "sha256": digest,
        }
        return destination

    if data["equilibrium_file"] is not None:
        equilibrium_file = copy_input(
            "equilibrium",
            data["equilibrium_file_alias"],
            data["equilibrium_file"],
            Path("equilibrium") / data["equilibrium_file"].name,
            "reduced equilibrium",
        )
        equilibrium_parameters_file = copy_input(
            "equilibrium_parameters",
            data["equilibrium_parameters_alias"],
            data["equilibrium_parameters_file"],
            Path("equilibrium") / data["equilibrium_parameters_file"].name,
            "equilibrium btor_rbig.dat",
        )
        equilibrium_inputs: tuple[Path, ...] = ()
    else:
        assert data["original_equilibrium"] is not None
        original = copy_input(
            "equilibrium",
            data["original_equilibrium_alias"],
            data["original_equilibrium"],
            Path("equilibrium-inputs") / data["original_equilibrium"].name,
            "original equilibrium",
        )
        explicit_inputs = tuple(
            copy_input(
                f"equilibrium_input:{source.name.casefold()}",
                alias,
                source,
                Path("equilibrium-inputs") / source.name,
                "equilibrium input",
            )
            for alias, source in zip(
                data["equilibrium_input_aliases"], data["equilibrium_input_paths"], strict=True
            )
        )
        equilibrium_file = None
        equilibrium_parameters_file = None
        equilibrium_inputs = (original, *explicit_inputs)

    equilibrium_executable = None
    if data["equilibrium_executable"] is not None:
        equilibrium_executable = copy_input(
            "equilibrium_executable",
            data["equilibrium_executable_alias"],
            data["equilibrium_executable"],
            Path("equilibrium-executable") / data["equilibrium_executable"].name,
            "equilibrium executable",
            executable=True,
        )
    oracle = copy_input(
        "oracle",
        data["oracle_alias"],
        data["oracle"],
        Path("oracle") / data["oracle"].name,
        "QL-Balance oracle",
    )
    return {
        "directory": snapshot_root,
        "entries": entries,
        "equilibrium_file": equilibrium_file,
        "equilibrium_parameters_file": equilibrium_parameters_file,
        "equilibrium_executable": equilibrium_executable,
        "equilibrium_inputs": equilibrium_inputs,
        "oracle": oracle,
    }


def _rewrite_report_paths(report_path: Path, private_root: Path, final_root: Path) -> None:
    try:
        payload = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise ExperimentalInputError(f"unable to read prepared report: {report_path}") from error

    def rewrite(value: object) -> object:
        if isinstance(value, str):
            try:
                return str(final_root / Path(value).relative_to(private_root))
            except ValueError:
                return value
        if isinstance(value, Mapping):
            return {str(key): rewrite(item) for key, item in value.items()}
        if isinstance(value, list):
            return [rewrite(item) for item in value]
        return value

    report_path.write_text(
        json.dumps(rewrite(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _hash_file(path: Path, label: str) -> str:
    try:
        if not path.is_file() or path.is_symlink():
            raise OSError("not a regular file")
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise ExperimentalInputError(
            f"{path}: required {label} is missing or unavailable"
        ) from error


def _hash_alias(alias: Path, resolved: Path, label: str) -> str:
    current = _resolve_source_file(alias, label)
    if current != resolved:
        raise ExperimentalInputError(f"{alias}: {label} changed its resolved path")
    return _hash_file(resolved, label)


def _ensure_new_destination(
    destination: Path,
    protected_roots: Sequence[Path],
    source_paths: Sequence[Path],
    parent_fd: int,
) -> None:
    try:
        os.stat(destination.name, dir_fd=parent_fd, follow_symlinks=False)
    except FileNotFoundError:
        pass
    except OSError as error:
        raise ExperimentalInputError(
            f"unable to inspect characterization destination: {destination}"
        ) from error
    else:
        raise ExperimentalInputError(f"characterization destination already exists: {destination}")
    for root in protected_roots:
        try:
            destination.relative_to(root)
            raise ExperimentalInputError(
                "characterization destination must not be a source directory or its descendant"
            )
        except ValueError:
            pass
        try:
            root.relative_to(destination)
            raise ExperimentalInputError(
                "characterization destination must not contain a source directory"
            )
        except ValueError:
            pass
    for source_path in source_paths:
        try:
            source_path.relative_to(destination)
        except ValueError:
            continue
        raise ExperimentalInputError("characterization destination must not contain a source file")


def _verify_unchanged_inputs(data: Mapping[str, Any], before: Mapping[str, str]) -> dict[str, str]:
    after = _preflight(data)
    if dict(after) != dict(before):
        raise ExperimentalInputError("an input changed during characterization")
    return after


def _prepared_profile_files(prepared: Any) -> dict[str, str]:
    try:
        return {
            role: str(getattr(prepared.config.profiles, field))
            for role, field in _PROFILE_CONFIG_FIELDS.items()
        }
    except (AttributeError, TypeError) as error:
        raise ExperimentalInputError(
            "prepared configuration does not expose profile filenames"
        ) from error


def _compare_all(
    data: Mapping[str, Any],
    prepared: Any,
    oracle: Any,
    prepared_profile_files: Mapping[str, str],
) -> dict[str, dict[str, dict[str, Any]]]:
    prepared_values = {
        role: _read_prepared_profile(prepared.profiles / prepared_profile_files[role])
        for role in _PROFILE_ROLES
    }
    result: dict[str, dict[str, dict[str, Any]]] = {}
    for domain_name, domain in data["domains"].items():
        result[domain_name] = {}
        for role in _PROFILE_ROLES:
            oracle_values = getattr(oracle, _ORACLE_FIELDS[role])
            tolerance = (data["tolerances"] or {}).get(f"{domain_name}/{role}")
            comparison = compare_profiles(
                oracle_radius_cm=oracle.r_out,
                oracle_values=oracle_values,
                prepared_radius_cm=prepared_values[role][0],
                prepared_values=prepared_values[role][1],
                domain_cm=domain,
                interpolation_direction=data["interpolation_direction"],
                method=data["interpolation_method"],
                relative_floor=(
                    data["relative_floors"][role] if data["relative_floors"] is not None else None
                ),
                tolerances=tolerance,
                resonance=(
                    (data["config"].setup.m_mode, data["config"].setup.n_mode)
                    if role == "q"
                    else None
                ),
            )
            result[domain_name][role] = _comparison_payload(comparison)
    return result


def _read_prepared_profile(path: Path) -> tuple[np.ndarray, np.ndarray]:
    try:
        data = np.loadtxt(path, comments="#", ndmin=2)
    except Exception as error:
        raise ExperimentalInputError(f"unable to read prepared profile: {path}") from error
    if data.ndim != 2 or data.shape[1] != 2 or data.shape[0] < 2:
        raise ExperimentalInputError(f"prepared profile must have two columns: {path}")
    if not np.all(np.isfinite(data)) or not np.all(np.diff(data[:, 0]) > 0.0):
        raise ExperimentalInputError(f"prepared profile is invalid: {path}")
    return data[:, 0], data[:, 1]


def _comparison_payload(comparison: ComparisonResult) -> dict[str, Any]:
    return {
        "measurements": {
            key: (
                float(value)
                if (value := getattr(comparison.measurements, key)) is not None
                else None
            )
            for key in _METRICS
        },
        "exclusions": _json_value(comparison.exclusions),
        "warnings": list(comparison.warnings),
        "threshold_decisions": (
            dict(comparison.threshold_decisions)
            if comparison.threshold_decisions is not None
            else None
        ),
        "overall_pass": comparison.overall_pass,
        "resonance": _json_value(comparison.resonance),
        "comparison_radius_cm": comparison.comparison_radius_cm.tolist(),
    }


def _status_for(
    comparisons: Mapping[str, Mapping[str, Mapping[str, Any]]],
) -> tuple[BalanceStatus, dict[str, bool] | None, bool | None]:
    decisions: dict[str, bool] = {}
    unthresholded = False
    for domain, profiles in comparisons.items():
        for role, payload in profiles.items():
            per_metric = payload["threshold_decisions"]
            if per_metric is None:
                unthresholded = True
                continue
            decisions.update(
                {f"{domain}/{role}/{key}": bool(value) for key, value in per_metric.items()}
            )
    if not decisions:
        return "MEASURED", None, None
    if unthresholded:
        return "PARTIALLY_EVALUATED", decisions, None
    overall = all(decisions.values())
    return ("PASS" if overall else "FAIL"), decisions, overall


def _success_report(
    data: Mapping[str, Any],
    actual_hashes: Mapping[str, str],
    staged: Any,
    prepared: Any,
    comparisons: Mapping[str, Mapping[str, Mapping[str, Any]]],
    status: BalanceStatus,
    threshold_decisions: Mapping[str, bool] | None,
    overall_pass: bool | None,
    prepared_profile_files: Mapping[str, str],
    prepared_profile_hashes: Mapping[str, str],
    snapshots: Mapping[str, Any],
    *,
    equilibrium_calculation: Mapping[str, Any],
    private_root: Path,
    final_root: Path,
    retained_staging: Mapping[str, object],
) -> dict[str, Any]:
    source_hashes = {key: actual_hashes[key] for key in _BALANCE_ROLES}
    report_comparisons: dict[str, dict[str, dict[str, Any]]] = {}
    for domain, profiles in comparisons.items():
        report_comparisons[domain] = {}
        for role, payload in profiles.items():
            report_comparisons[domain][role] = {
                "role": {
                    "source": _SOURCE_ROLES[role],
                    "prepared": prepared_profile_files[role],
                    "oracle": _ORACLE_DATASETS[role],
                },
                "units": {"prepared": _PROFILE_UNITS[role], "oracle": _PROFILE_UNITS[role]},
                **payload,
            }
    generator = equilibrium_calculation["generator"]
    if isinstance(generator, Mapping):
        generator = dict(generator)
        for key in ("executed_executable",):
            if isinstance(generator.get(key), str):
                path = Path(generator[key])
                try:
                    generator[key] = str(final_root / path.relative_to(private_root))
                except ValueError:
                    pass
        source_command = generator.get("source_command")
        if isinstance(source_command, list):
            generator["source_command"] = [
                (
                    str(final_root / Path(item).relative_to(private_root))
                    if isinstance(item, str)
                    and Path(item).is_absolute()
                    and Path(item).is_relative_to(private_root)
                    else item
                )
                for item in source_command
            ]
    equilibrium_report = {
        "equilibrium_file": str(
            _published_path(equilibrium_calculation["equilibrium_file"], private_root, final_root)
        ),
        "equilibrium_sha256": equilibrium_calculation["equilibrium_sha256"],
        "parameters_file": str(
            _published_path(equilibrium_calculation["parameters_file"], private_root, final_root)
        ),
        "parameters_sha256": equilibrium_calculation["parameters_sha256"],
        "btor_gauss": equilibrium_calculation["btor_gauss"],
        "r_big_cm": equilibrium_calculation["r_big_cm"],
        "generator": generator,
        "used_for": [
            "KIM setup btor and major_radius",
            "BALANCE angular rotation to toroidal velocity conversion",
        ],
    }
    return {
        "schema_version": 1,
        "status": status,
        "report_path": str(final_root / "characterization_report.json"),
        "retained_staging": dict(retained_staging),
        "overall_pass": overall_pass,
        "threshold_decisions": threshold_decisions,
        "source_hashes": source_hashes,
        "expected_sha256": data["expected_sha256"],
        "expected_hashes_requested": data["expected_sha256"] is not None,
        "actual_sha256": dict(actual_hashes),
        "request_provenance": {
            "metadata": _provenance_entry(data, "metadata"),
            "config": _provenance_entry(data, "config"),
        },
        "input_provenance": {
            key: {
                **entry,
                "snapshot_path": str(
                    _published_path(Path(entry["snapshot_path"]), private_root, final_root)
                ),
            }
            for key, entry in snapshots["entries"].items()
        },
        "staged": {
            "directory": str(_published_path(staged.directory, private_root, final_root)),
            "report": str(_published_path(staged.report, private_root, final_root)),
        },
        "prepared": {
            "directory": str(_published_path(prepared.directory, private_root, final_root)),
            "profiles": str(_published_path(prepared.profiles, private_root, final_root)),
            "equilibrium": str(_published_path(prepared.equilibrium, private_root, final_root)),
            "request": str(_published_path(prepared.request, private_root, final_root)),
            "report": str(_published_path(prepared.report, private_root, final_root)),
            "profile_hashes": {
                prepared_profile_files[role]: digest
                for role, digest in prepared_profile_hashes.items()
            },
            "profile_hashes_by_role": dict(prepared_profile_hashes),
            "equilibrium_parameters": (
                {
                    "file": str(
                        _published_path(prepared.equilibrium_parameters, private_root, final_root)
                    ),
                    "sha256": equilibrium_calculation["parameters_sha256"],
                    "btor_gauss": prepared.config.setup.btor,
                    "r_big_cm": prepared.config.setup.major_radius,
                }
                if prepared.equilibrium_parameters is not None
                else None
            ),
        },
        "equilibrium_calculation": equilibrium_report,
        "comparison_configuration": {
            "domains": data["domains"],
            "interpolation_direction": data["interpolation_direction"],
            "interpolation_method": data["interpolation_method"],
            "relative_floors": data["relative_floors"],
            "q_operation": data["q_operation"],
            "btor_gauss": data["config"].setup.btor,
            "r_big_cm": data["config"].setup.major_radius,
            "resonance_mode": {
                "m": data["config"].setup.m_mode,
                "n": data["config"].setup.n_mode,
                "convention": "q=-m/n",
            },
            "tolerances": data["tolerances_input"],
        },
        "comparisons": report_comparisons,
        "limitations": dict(_LIMITATIONS),
    }


def _published_path(path: Path, private_root: Path, final_root: Path) -> Path:
    try:
        relative = path.relative_to(private_root)
    except ValueError as error:
        raise ExperimentalInputError(f"generated path escaped private output: {path}") from error
    return final_root / relative


def _provenance_entry(data: Mapping[str, Any], key: str) -> dict[str, object]:
    source = data[f"{key}_source"]
    if source is None:
        return {"mode": "in_memory"}
    return {
        "mode": "file_snapshot",
        "path": str(source),
        "sha256": data["snapshot_hashes"][key],
    }


def _unavailable(
    message: str,
    error: Exception | None = None,
    *,
    retained_staging: Mapping[str, object] | None = None,
) -> BalanceCharacterization:
    details = str(error) if error is not None and str(error) else message
    report = {
        "schema_version": 1,
        "status": "UNAVAILABLE",
        "overall_pass": None,
        "threshold_decisions": None,
        "error": {"message": details},
        "retained_staging": dict(
            retained_staging
            if retained_staging is not None
            else _retained_staging_info(None, "not_created")
        ),
        "limitations": dict(_LIMITATIONS),
    }
    return BalanceCharacterization("UNAVAILABLE", report)


def _json_dump(value: object) -> str:
    return json.dumps(_json_value(value), indent=2, sort_keys=True) + "\n"


def _json_value(value: object) -> object:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    if hasattr(value, "__dataclass_fields__"):
        return {key: _json_value(getattr(value, key)) for key in value.__dataclass_fields__}
    return str(value)


__all__ = [
    "BalanceCharacterization",
    "BalanceCharacterizationRequest",
    "BalanceStatus",
    "characterize_balance",
]
