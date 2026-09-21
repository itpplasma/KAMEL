from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from kim import (
    BalanceMetadata,
    BuiltinPlasma,
    ExperimentalInputError,
    PlasmaIsotope,
    SimulationConfig,
)
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
        "toroidal_rotation": [0.0, 0.25, 1.0],
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


def test_composes_balance_staging_provenance_into_prepared_report(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    balance_paths, balance_source, staged, marsf = _stage_balance_case(tmp_path)
    equilibrium = tmp_path / "original-equilibrium.dat"
    _write_equilibrium(equilibrium)

    monkeypatch.setattr(
        preparation_module,
        "_software_identity",
        lambda: {"name": "kamel-kim", "version": "test-version"},
        raising=False,
    )
    monkeypatch.setattr(
        preparation_module,
        "_git_identity",
        lambda: {"commit": "0123456789abcdef" * 5, "dirty": False},
        raising=False,
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

    upstream = report["upstream_staging"]
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
    assert report["git"] == {"commit": "0123456789abcdef" * 5, "dirty": False}
    assert report["comparison"]["domain"] is None

    np.testing.assert_allclose(
        np.loadtxt(prepared.profiles / "n.dat")[:, 1], [1.0e13, 2.0e13, 3.0e13]
    )
    np.testing.assert_allclose(
        np.loadtxt(prepared.profiles / "Vz.dat")[:, 1], [0.0, -4.125e6, -8.25e6]
    )


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
    monkeypatch.setattr(preparation_module, "_software_identity", lambda: None, raising=False)
    monkeypatch.setattr(preparation_module, "_git_identity", lambda: None, raising=False)

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

    with pytest.raises(ExperimentalInputError, match="upstream.*staging.*report"):
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

    with pytest.raises(ExperimentalInputError, match="upstream.*staging.*report"):
        prepare_marsf_case(
            marsf,
            _config(),
            tmp_path / "prepared",
            equilibrium_file=equilibrium,
            upstream_staging_report=malformed,
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

    with pytest.raises(ExperimentalInputError, match="upstream.*staging.*(mismatch|hash)"):
        prepare_marsf_case(
            marsf_b,
            _config(),
            tmp_path / "prepared",
            equilibrium_file=equilibrium,
            upstream_staging_report=staged_a.report,
        )
