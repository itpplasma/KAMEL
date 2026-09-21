from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import stat
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path
from types import MappingProxyType

import numpy as np
import pytest
from kim import (
    BalanceMetadata,
    BuiltinPlasma,
    ExperimentalInputError,
    ExperimentalProfile,
    MarsFInput,
    PlasmaIsotope,
    SimulationConfig,
)
from kim import executable as executable_module
from kim import preparation as preparation_module
from kim import (
    prepare_marsf_case,
    read_balance_profiles,
    read_marsf_profiles,
)
from kim.balance_adoption import stage_balance_marsf_quartet

_BALANCE_ROLES = (
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
_SOURCE_METADATA = {
    "source": "synthetic-aug-balance",
    "coordinate": "rho_pol",
    "coordinate_unit": "1",
    "density_unit": "1/m^3",
    "electron_temperature_unit": "eV",
    "ion_temperature_unit": "eV",
    "toroidal_rotation_unit": "rad/s",
}


def _balance_metadata(**updates: object) -> BalanceMetadata:
    values = dict(_SOURCE_METADATA)
    values.update(updates)
    return BalanceMetadata(**values)


def _write_balance_profiles(directory: Path, *, density_offset: float = 0.0) -> dict[str, Path]:
    directory.mkdir(parents=True)
    grids = {
        "density": [0.0, 0.5, 1.0],
        "electron_temperature": [0.0, 0.4, 1.0],
        "ion_temperature": [0.0, 0.3, 1.0],
        # Align angular rotation with the normalized equilibrium nodes so the
        # R0 conversion is checked without an extra cubic interpolation.
        "toroidal_rotation": [0.0, 0.5, 1.0],
    }
    values = {
        "density": [1.0e19 + density_offset, 2.0e19, 3.0e19],
        "electron_temperature": [100.0, 80.0, 40.0],
        "ion_temperature": [90.0, 70.0, 30.0],
        "toroidal_rotation": [0.0, -2.5e4, -5.0e4],
    }
    paths: dict[str, Path] = {}
    for role in _BALANCE_ROLES:
        path = directory / f"aug-{role}.profile"
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


def _read_balance_case(paths: dict[str, Path]):
    return read_balance_profiles(
        density=paths["density"],
        electron_temperature=paths["electron_temperature"],
        ion_temperature=paths["ion_temperature"],
        toroidal_rotation=paths["toroidal_rotation"],
        metadata=_balance_metadata(),
    )


def _write_equilibrium(path: Path) -> None:
    path.write_text(
        "# radius q psi\n" "0.0 1.0 0.0\n" "10.0 1.5 0.25\n" "20.0 2.0 1.0\n",
        encoding="utf-8",
    )


def _config() -> SimulationConfig:
    return SimulationConfig.electrostatic_periodic(
        profiles=Path("unused"),
        plasma=BuiltinPlasma(isotope=PlasmaIsotope.DEUTERIUM),
        btor=-17_977.413,
        major_radius=_R0_CM,
        m_mode=7,
        n_mode=2,
        frequency=0.0,
        br_boundary_real=1.0,
        br_boundary_imag=0.0,
        radial_minimum=3.0,
        plasma_radius=19.0,
    )


def _stage_balance_case(tmp_path: Path, *, density_offset: float = 0.0):
    paths = _write_balance_profiles(
        tmp_path / f"balance-{density_offset:g}", density_offset=density_offset
    )
    source = _read_balance_case(paths)
    staged = stage_balance_marsf_quartet(
        source,
        tmp_path / f"staged-{density_offset:g}",
        major_radius_cm=_R0_CM,
        equilibrium_provenance="original-aug-equilibrium",
    )
    marsf = read_marsf_profiles(staged.directory, staged.metadata)
    return paths, source, staged, marsf


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _report(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_staging_report_variant(
    source_report: Path, destination: Path, mutate: Callable[[dict[str, object]], None]
) -> Path:
    payload = _report(source_report)
    mutate(payload)
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return destination


def _relink_marsf_source(source, report: Path):
    metadata = source.metadata.model_copy(update={"upstream_staging_sha256": _sha256(report)})
    return replace(source, metadata=metadata)


def _patch_identity_seams(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    *,
    commit: str | None,
    software_available: bool = True,
) -> None:
    git_metadata = (
        None
        if commit is None
        else executable_module.KamelGitMetadata(
            root=tmp_path / "synthetic-checkout", commit=commit, dirty=False
        )
    )
    monkeypatch.setattr(
        executable_module,
        "discover_kamel_git_metadata",
        lambda *args, **kwargs: git_metadata,
    )
    if software_available:
        monkeypatch.setattr(importlib.metadata, "version", lambda package: "test-version")
    else:

        def missing_version(package: str) -> str:
            raise importlib.metadata.PackageNotFoundError(package)

        monkeypatch.setattr(importlib.metadata, "version", missing_version)


def test_composes_balance_staging_provenance_into_prepared_report(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    balance_paths, balance_source, staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)

    _patch_identity_seams(
        monkeypatch,
        tmp_path,
        commit="0123456789abcdef0123456789abcdef01234567",
    )

    prepared = prepare_marsf_case(
        marsf,
        _config(),
        tmp_path / "prepared",
        equilibrium_file=equilibrium,
        q_operation="negate",
        upstream_staging_report=staged.report,
    )
    report = _report(prepared.report)

    assert report["schema_version"] == 1
    upstream = report["upstream_staging"]
    assert upstream["schema_version"] == 1
    assert upstream["source_basenames"] == {role: path.name for role, path in balance_paths.items()}
    assert upstream["source_hashes"] == {
        role: _sha256(path) for role, path in balance_paths.items()
    }
    assert upstream["source_metadata"] == balance_source.metadata.model_dump(mode="json")
    assert upstream["derived_hashes"] == {
        filename: _sha256(staged.directory / filename) for filename in _MARSF_FILENAMES.values()
    }

    density_operation = next(
        operation for operation in report["operations"] if operation["quantity"] == "density"
    )
    assert density_operation == {
        "quantity": "density",
        "source_unit": "1/m^3",
        "target_unit": "1/cm^3",
        "factor": 1.0e-6,
    }
    rotation_operation = next(
        operation
        for operation in upstream["operations"]
        if operation["quantity"] == "toroidal_velocity"
    )
    assert rotation_operation["source_unit"] == "rad/s"
    assert rotation_operation["target_unit"] == "cm/s"
    assert rotation_operation["factor"] == _R0_CM
    assert rotation_operation["operation"] == "omega_to_v_phi"
    assert rotation_operation["parameters"] == {"major_radius_cm": _R0_CM}

    assert report["equilibrium"]["source_hash"] == _sha256(equilibrium)
    assert report["equilibrium"]["operation"] == "copied from explicit equilibrium_file"
    assert report["coordinate_operation"] == {
        "source_coordinate": "sqrt_psiN",
        "target_coordinate": "r_eff",
        "method": "natural cubic interpolation",
    }
    q_operation = next(
        operation for operation in report["operations"] if operation["quantity"] == "q"
    )
    assert q_operation["operation"] == "negate"
    assert q_operation["factor"] == -1.0

    assert report["prepared_output_hashes"] == {
        f"profiles/{filename}": _sha256(prepared.profiles / filename)
        for filename in ("n.dat", "Te.dat", "Ti.dat", "Vz.dat", "q.dat")
    }
    assert report["software"] == {"name": "kamel-kim", "version": "test-version"}
    assert report["git"] == {
        "commit": "0123456789abcdef0123456789abcdef01234567",
        "dirty": False,
    }
    assert report["comparison"]["domain"] is None
    upstream_digest = _sha256(staged.report)
    assert upstream["provenance_sha256"] == upstream_digest
    assert upstream_digest == upstream_digest.lower()

    np.testing.assert_array_equal(
        np.loadtxt(prepared.profiles / "n.dat")[:, 1], [1.0e13, 2.0e13, 3.0e13]
    )
    np.testing.assert_array_equal(
        np.loadtxt(prepared.profiles / "Vz.dat")[:, 1], [0.0, -4.125e6, -8.25e6]
    )


def test_upstream_staging_radius_must_match_kim_configuration(
    tmp_path: Path,
) -> None:
    _balance_paths, _balance_source, staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)
    mismatched_config = _config().model_copy(
        update={"setup": _config().setup.model_copy(update={"major_radius": 180.0})}
    )
    destination = tmp_path / "prepared"

    with pytest.raises(ExperimentalInputError):
        prepare_marsf_case(
            marsf,
            mismatched_config,
            destination,
            equilibrium_file=equilibrium,
            upstream_staging_report=staged.report,
        )

    assert not destination.exists()


def test_mutated_marsf_arrays_are_rejected_before_preparation(
    tmp_path: Path,
) -> None:
    _balance_paths, _balance_source, staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)
    density = marsf.profiles["density"]
    mutable_values = density.values.copy()
    mutable_values[0] += 1.0e12
    mutable_profiles = dict(marsf.profiles)
    mutable_profiles["density"] = ExperimentalProfile(
        name=density.name,
        units=density.units,
        path=density.path,
        coordinate=density.coordinate,
        values=mutable_values,
    )
    mutable_source = MarsFInput(
        directory=marsf.directory,
        metadata=marsf.metadata,
        profiles=MappingProxyType(mutable_profiles),
        source_files=marsf.source_files,
    )
    destination = tmp_path / "prepared"

    with pytest.raises(ExperimentalInputError):
        prepare_marsf_case(
            mutable_source,
            _config(),
            destination,
            equilibrium_file=equilibrium,
            upstream_staging_report=staged.report,
        )

    assert not destination.exists()


def test_replaced_marsf_file_is_rejected_before_direct_preparation(
    tmp_path: Path,
) -> None:
    _balance_paths, _balance_source, staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)
    (staged.directory / "PROFDEN.IN").write_text(
        "MARS-F profile\n0.0 1.1e13\n0.5 2.0e13\n1.0 3.0e13\n",
        encoding="utf-8",
    )
    destination = tmp_path / "prepared"

    with pytest.raises(ExperimentalInputError):
        prepare_marsf_case(marsf, _config(), destination, equilibrium_file=equilibrium)

    assert not destination.exists()


@pytest.mark.parametrize(
    "report_bytes",
    [b"\xff\xfe"],
)
def test_malformed_upstream_report_bytes_raise_experimental_input_error(
    tmp_path: Path, report_bytes: bytes
) -> None:
    _balance_paths, _balance_source, staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)
    malformed = tmp_path / "malformed-staging-report.json"
    malformed.write_bytes(report_bytes)
    relinked_source = _relink_marsf_source(marsf, malformed)

    with pytest.raises(ExperimentalInputError):
        prepare_marsf_case(
            relinked_source,
            _config(),
            tmp_path / "prepared",
            equilibrium_file=equilibrium,
            upstream_staging_report=malformed,
        )


def test_nested_staging_schema_types_are_strict(
    tmp_path: Path,
) -> None:
    _balance_paths, _balance_source, staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)
    malformed = tmp_path / "bool-schema-version.json"

    def mutate(payload: dict[str, object]) -> None:
        source_metadata = payload["source_metadata"]
        assert isinstance(source_metadata, dict)
        source_metadata["schema_version"] = True

    _write_staging_report_variant(staged.report, malformed, mutate)
    relinked_source = _relink_marsf_source(marsf, malformed)

    with pytest.raises(ExperimentalInputError):
        prepare_marsf_case(
            relinked_source,
            _config(),
            tmp_path / "prepared",
            equilibrium_file=equilibrium,
            upstream_staging_report=malformed,
        )


def test_huge_staging_radius_is_rejected_as_experimental_input_error(
    tmp_path: Path,
) -> None:
    _balance_paths, _balance_source, staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)
    malformed = tmp_path / "huge-radius.json"

    def mutate(payload: dict[str, object]) -> None:
        payload["major_radius_cm"] = 10**1000

    _write_staging_report_variant(staged.report, malformed, mutate)
    relinked_source = _relink_marsf_source(marsf, malformed)

    with pytest.raises(ExperimentalInputError):
        prepare_marsf_case(
            relinked_source,
            _config(),
            tmp_path / "prepared",
            equilibrium_file=equilibrium,
            upstream_staging_report=malformed,
        )


def test_huge_staging_operation_factor_is_rejected_as_experimental_input_error(
    tmp_path: Path,
) -> None:
    _balance_paths, _balance_source, staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)
    malformed = tmp_path / "huge-operation-factor.json"

    def mutate(payload: dict[str, object]) -> None:
        operations = payload["operations"]
        assert isinstance(operations, list)
        rotation = operations[-1]
        assert isinstance(rotation, dict)
        rotation["factor"] = 10**1000

    _write_staging_report_variant(staged.report, malformed, mutate)
    relinked_source = _relink_marsf_source(marsf, malformed)
    destination = tmp_path / "prepared"

    with pytest.raises(ExperimentalInputError):
        prepare_marsf_case(
            relinked_source,
            _config(),
            destination,
            equilibrium_file=equilibrium,
            upstream_staging_report=malformed,
        )

    assert not destination.exists()


def test_direct_marsf_preparation_records_absent_upstream_provenance_honestly(
    tmp_path: Path,
) -> None:
    _balance_paths, _balance_source, staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)

    prepared = prepare_marsf_case(
        marsf, _config(), tmp_path / "prepared", equilibrium_file=equilibrium
    )

    report = _report(prepared.report)
    assert report["upstream_staging"] is None
    assert report["comparison"]["domain"] is None


def test_missing_software_and_git_identity_are_recorded_as_absent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _balance_paths, _balance_source, _staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)
    _patch_identity_seams(monkeypatch, tmp_path, commit=None, software_available=False)

    prepared = prepare_marsf_case(
        marsf, _config(), tmp_path / "prepared", equilibrium_file=equilibrium
    )

    report = _report(prepared.report)
    assert report["software"] is None
    assert report["git"] is None


def test_missing_upstream_staging_report_is_rejected(tmp_path: Path) -> None:
    _balance_paths, _balance_source, _staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)

    with pytest.raises(ExperimentalInputError):
        prepare_marsf_case(
            marsf,
            _config(),
            tmp_path / "prepared",
            equilibrium_file=equilibrium,
            upstream_staging_report=tmp_path / "missing-staging-report.json",
        )


def test_malformed_upstream_staging_report_is_rejected(tmp_path: Path) -> None:
    _balance_paths, _balance_source, _staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)
    malformed = tmp_path / "malformed-staging-report.json"
    malformed.write_text("{not-json", encoding="utf-8")

    with pytest.raises(ExperimentalInputError):
        prepare_marsf_case(
            marsf,
            _config(),
            tmp_path / "prepared",
            equilibrium_file=equilibrium,
            upstream_staging_report=malformed,
        )


@pytest.mark.parametrize(
    "variant",
    ["missing-required-field", "unsupported-schema-version"],
)
def test_rejects_valid_json_upstream_reports_with_invalid_schema(
    tmp_path: Path, variant: str
) -> None:
    _balance_paths, _balance_source, staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)
    malformed = tmp_path / f"{variant}.json"

    def mutate(payload: dict[str, object]) -> None:
        if variant == "missing-required-field":
            payload.pop("source_basenames")
        else:
            payload["schema_version"] = 999

    _write_staging_report_variant(staged.report, malformed, mutate)

    with pytest.raises(ExperimentalInputError):
        prepare_marsf_case(
            marsf,
            _config(),
            tmp_path / "prepared",
            equilibrium_file=equilibrium,
            upstream_staging_report=malformed,
        )


def test_arbitrary_upstream_json_is_rejected(tmp_path: Path) -> None:
    _balance_paths, _balance_source, staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)
    augmented = tmp_path / "augmented-staging-report.json"

    def mutate(payload: dict[str, object]) -> None:
        payload["arbitrary_extra"] = {"must_not": "cross_boundary"}

    _write_staging_report_variant(staged.report, augmented, mutate)

    with pytest.raises(ExperimentalInputError):
        prepare_marsf_case(
            marsf,
            _config(),
            tmp_path / "prepared",
            equilibrium_file=equilibrium,
            upstream_staging_report=augmented,
        )


@pytest.mark.parametrize("tampered_field", ["source_metadata", "source_basenames"])
def test_tampered_upstream_identity_is_rejected_even_when_derived_hashes_match(
    tmp_path: Path, tampered_field: str
) -> None:
    _balance_paths, balance_source, staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)
    tampered = tmp_path / f"tampered-{tampered_field}.json"

    def mutate(payload: dict[str, object]) -> None:
        if tampered_field == "source_metadata":
            source_metadata = payload.setdefault(
                "source_metadata", balance_source.metadata.model_dump(mode="json")
            )
            assert isinstance(source_metadata, dict)
            source_metadata["density_unit"] = "1/cm^3"
        else:
            source_basenames = payload["source_basenames"]
            assert isinstance(source_basenames, dict)
            source_basenames["density"] = "tampered-density.profile"

    _write_staging_report_variant(staged.report, tampered, mutate)
    original_derived_hashes = _report(staged.report)["derived_hashes"]
    assert _report(tampered)["derived_hashes"] == original_derived_hashes
    assert _sha256(tampered) != _sha256(staged.report)

    with pytest.raises(ExperimentalInputError):
        prepare_marsf_case(
            marsf,
            _config(),
            tmp_path / "prepared",
            equilibrium_file=equilibrium,
            upstream_staging_report=tampered,
        )


def test_mismatched_upstream_staging_report_is_rejected(tmp_path: Path) -> None:
    _balance_paths_a, _balance_source_a, staged_a, _marsf_a = _stage_balance_case(
        tmp_path / "case-a"
    )
    _balance_paths_b, _balance_source_b, _staged_b, marsf_b = _stage_balance_case(
        tmp_path / "case-b", density_offset=1.0e18
    )
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)

    with pytest.raises(ExperimentalInputError):
        prepare_marsf_case(
            marsf_b,
            _config(),
            tmp_path / "prepared",
            equilibrium_file=equilibrium,
            upstream_staging_report=staged_a.report,
        )


@pytest.mark.parametrize(
    ("operation_index", "field"),
    [
        (0, "quantity"),
        (0, "source_unit"),
        (0, "target_unit"),
        (0, "operation"),
        (3, "quantity"),
        (3, "source_unit"),
        (3, "target_unit"),
        (3, "operation"),
    ],
)
def test_tampered_staging_operation_semantics_are_rejected_before_preparation(
    tmp_path: Path, operation_index: int, field: str
) -> None:
    _balance_paths, _balance_source, staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)
    tampered = tmp_path / f"tampered-operation-{operation_index}-{field}.json"

    def mutate(payload: dict[str, object]) -> None:
        operations = payload["operations"]
        assert isinstance(operations, list)
        operation = operations[operation_index]
        assert isinstance(operation, dict)
        operation[field] = "tampered"

    _write_staging_report_variant(staged.report, tampered, mutate)
    relinked_source = _relink_marsf_source(marsf, tampered)
    destination = tmp_path / "prepared"

    with pytest.raises(ExperimentalInputError):
        prepare_marsf_case(
            relinked_source,
            _config(),
            destination,
            equilibrium_file=equilibrium,
            upstream_staging_report=tampered,
        )

    assert not destination.exists()


@pytest.mark.parametrize("mutation", ["duplicate", "missing", "extra"])
def test_staging_operation_role_set_must_be_exact(tmp_path: Path, mutation: str) -> None:
    _balance_paths, _balance_source, staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)
    tampered = tmp_path / f"tampered-operation-{mutation}.json"

    def mutate(payload: dict[str, object]) -> None:
        operations = payload["operations"]
        assert isinstance(operations, list)
        if mutation == "duplicate":
            operations[0] = dict(operations[1])
        elif mutation == "missing":
            operations.pop()
        else:
            operations.append(dict(operations[-1]))

    _write_staging_report_variant(staged.report, tampered, mutate)
    relinked_source = _relink_marsf_source(marsf, tampered)

    with pytest.raises(ExperimentalInputError):
        prepare_marsf_case(
            relinked_source,
            _config(),
            tmp_path / "prepared",
            equilibrium_file=equilibrium,
            upstream_staging_report=tampered,
        )


def test_generated_equilibrium_provenance_records_inputs_method_and_identity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _balance_paths, _balance_source, staged, marsf = _stage_balance_case(tmp_path)
    generator = tmp_path / "synthetic-equilibrium-generator.py"
    generator.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "if Path('equilibrium-control.in').read_text(encoding='utf-8') != 'synthetic control\\n':\n"
        "    raise SystemExit('unexpected copied equilibrium control')\n"
        "Path('equil_r_q_psi.dat').write_text(\n"
        "    '# radius q psi\\n0.0 1.0 0.0\\n10.0 1.5 0.25\\n20.0 2.0 1.0\\n',\n"
        "    encoding='utf-8',\n"
        ")\n",
        encoding="utf-8",
    )
    generator.chmod(generator.stat().st_mode | stat.S_IXUSR)
    equilibrium_input = tmp_path / "equilibrium-control.in"
    equilibrium_input.write_text("synthetic control\n", encoding="utf-8")
    generator_bytes = generator.read_bytes()
    equilibrium_input_bytes = equilibrium_input.read_bytes()
    generator_digest = _sha256(generator)
    equilibrium_input_digest = _sha256(equilibrium_input)
    _patch_identity_seams(
        monkeypatch,
        tmp_path,
        commit="fedcba9876543210fedcba9876543210fedcba98",
    )

    prepared = prepare_marsf_case(
        marsf,
        _config(),
        tmp_path / "prepared",
        equilibrium_executable=generator,
        equilibrium_input_files=(equilibrium_input,),
        q_operation="negate",
        upstream_staging_report=staged.report,
    )
    report = _report(prepared.report)
    generated_equilibrium = prepared.equilibrium

    assert report["equilibrium"]["operation"] == "generated by supplied equilibrium executable"
    assert report["equilibrium"]["source_hash"] == _sha256(generated_equilibrium)
    assert report["source_hashes"][f"equilibrium_input/{equilibrium_input.name}"] == _sha256(
        equilibrium_input
    )
    assert report["generator"]["executable"] == str(generator.absolute())
    assert report["generator"]["sha256"] == generator_digest
    assert report["generator"]["command_file_hashes"][str(generator.absolute())] == generator_digest
    assert report["source_hashes"][f"equilibrium_input/{equilibrium_input.name}"] == (
        equilibrium_input_digest
    )
    assert report["software"] == {"name": "kamel-kim", "version": "test-version"}
    assert report["git"] == {
        "commit": "fedcba9876543210fedcba9876543210fedcba98",
        "dirty": False,
    }
    assert generator.read_bytes() == generator_bytes
    assert equilibrium_input.read_bytes() == equilibrium_input_bytes


def test_generator_executes_verified_snapshot_when_original_changes_after_hash(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _balance_paths, _balance_source, _staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "unused-equilibrium.dat"
    generator = tmp_path / "switching-generator.py"
    generator_a = (
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "Path('equil_r_q_psi.dat').write_text("
        "'# radius q psi\\n0.0 1.0 0.0\\n10.0 1.5 0.25\\n20.0 2.0 1.0\\n', "
        "encoding='utf-8')\n"
    ).encode("utf-8")
    generator_b = generator_a.replace(b"0.0 1.0 0.0", b"0.0 9.0 0.0")
    generator.write_bytes(generator_a)
    generator.chmod(generator.stat().st_mode | stat.S_IXUSR)
    generator_digest = _sha256(generator)

    real_run = preparation_module._run_equilibrium_generator

    def replace_original_then_run(
        command: list[str], working_directory: Path, timeout_seconds: float
    ) -> None:
        generator.write_bytes(generator_b)
        return real_run(command, working_directory, timeout_seconds)

    monkeypatch.setattr(preparation_module, "_run_equilibrium_generator", replace_original_then_run)
    prepared = prepare_marsf_case(
        marsf,
        _config(),
        tmp_path / "prepared",
        equilibrium_executable=generator,
    )

    np.testing.assert_array_equal(np.loadtxt(prepared.equilibrium)[:, 1], [1.0, 1.5, 2.0])
    report = _report(prepared.report)
    assert report["generator"]["executable"] == str(generator.absolute())
    assert report["generator"]["sha256"] == generator_digest
    assert report["generator"]["executed_sha256"] == generator_digest
    assert (prepared.directory / report["generator"]["executed_executable"]).is_file()
    assert generator.read_bytes() == generator_b
    assert equilibrium.exists() is False


def test_generator_bare_path_command_is_resolved_and_verified(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _balance_paths, _balance_source, _staged, marsf = _stage_balance_case(tmp_path)
    executable_directory = tmp_path / "bin"
    executable_directory.mkdir()
    generator = executable_directory / "marsf-generator"
    generator.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "Path('equil_r_q_psi.dat').write_text("
        "'# radius q psi\\n0.0 1.0 0.0\\n10.0 1.5 0.25\\n20.0 2.0 1.0\\n', "
        "encoding='utf-8')\n",
        encoding="utf-8",
    )
    generator.chmod(generator.stat().st_mode | stat.S_IXUSR)
    work_directory = tmp_path / "work"
    work_directory.mkdir()
    monkeypatch.chdir(work_directory)
    monkeypatch.setenv("PATH", f"{executable_directory}{os.pathsep}{os.environ['PATH']}")

    prepared = prepare_marsf_case(
        marsf,
        _config(),
        tmp_path / "prepared",
        equilibrium_executable=generator.name,
    )

    report = _report(prepared.report)
    assert report["generator"]["executable"] == str(generator.absolute())
    assert report["generator"]["sha256"] == _sha256(generator)
    assert (prepared.directory / report["generator"]["executed_executable"]).is_file()


def test_optional_identity_discovery_failures_record_null(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _balance_paths, _balance_source, _staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)

    def unavailable_version(_package: str) -> str:
        raise RuntimeError("metadata backend unavailable")

    def unavailable_git(*_args: object, **_kwargs: object) -> None:
        raise PermissionError("git metadata unavailable")

    monkeypatch.setattr(importlib.metadata, "version", unavailable_version)
    monkeypatch.setattr(executable_module, "discover_kamel_git_metadata", unavailable_git)
    prepared = prepare_marsf_case(
        marsf, _config(), tmp_path / "prepared", equilibrium_file=equilibrium
    )

    report = _report(prepared.report)
    assert report["software"] is None
    assert report["git"] is None
