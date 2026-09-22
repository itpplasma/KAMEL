from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pytest
from kim import (
    BalanceMetadata,
    BuiltinPlasma,
    ExperimentalInputError,
    PlasmaIsotope,
    SimulationConfig,
    compare_profiles,
)
from kim import preparation as preparation_module
from kim import (
    prepare_marsf_case,
    read_balance_profiles,
    read_marsf_profiles,
    read_ql_balance_oracle,
)
from kim.balance_adoption import stage_balance_marsf_quartet

_R0_CM = 200.0
_EQUILIBRIUM_PROVENANCE = "synthetic-aug-equilibrium"
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


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _balance_metadata() -> BalanceMetadata:
    return BalanceMetadata(
        source="synthetic-aug-balance",
        coordinate="rho_pol",
        coordinate_unit="1",
        density_unit="1/m^3",
        electron_temperature_unit="eV",
        ion_temperature_unit="eV",
        toroidal_rotation_unit="rad/s",
    )


def _write_balance_profiles(directory: Path) -> dict[str, Path]:
    directory.mkdir()
    # Values are linear in x = rho_pol**2.  Each role deliberately uses a
    # different source grid so the reader never gets a shared-grid shortcut.
    grids = {
        "density": [0.0, 0.5, 1.0],
        "electron_temperature": [0.0, 0.5, 0.75, 1.0],
        "ion_temperature": [0.0, 0.25, 0.5, 1.0],
        "toroidal_rotation": [0.0, 0.4, 0.5, 1.0],
    }
    values = {
        "density": lambda x: 1.0e18 + 2.0e18 * x,
        "electron_temperature": lambda x: 100.0 - 60.0 * x,
        "ion_temperature": lambda x: 8.0 - 4.0 * x,
        "toroidal_rotation": lambda x: 10.0 - 100.0 * x,
    }
    paths: dict[str, Path] = {}
    for role in _BALANCE_ROLES:
        path = directory / f"{role}.dat"
        path.write_text(
            "\n".join(f"{rho:.17g} {values[role](rho * rho):.17g}" for rho in grids[role]) + "\n",
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
    # r_eff = 1 + 4*psi_pol_norm and q is positive in the supplied table.
    path.write_text(
        "# r_eff q psi_pol\n" "1.0 2.0 0.0\n" "3.0 3.5 0.5\n" "5.0 5.0 1.0\n",
        encoding="utf-8",
    )


def _config() -> SimulationConfig:
    return SimulationConfig.electrostatic_periodic(
        profiles=Path("unused"),
        plasma=BuiltinPlasma(isotope=PlasmaIsotope.DEUTERIUM),
        btor=-18_000.0,
        major_radius=_R0_CM,
        m_mode=7,
        n_mode=2,
        frequency=0.0,
        br_boundary_real=1.0,
        br_boundary_imag=0.0,
        radial_minimum=1.0,
        plasma_radius=5.0,
    )


def _write_ql_balance_oracle(path: Path) -> None:
    # Independent oracle data: these arrays are not read from the prepared
    # case and the HDF file is created only after staging and preparation.
    r_out = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    psi = (r_out - 1.0) / 4.0
    with h5py.File(path, "w") as handle:
        datasets = {
            "preprocprof/r_out": r_out,
            "preprocprof/n": 1.0e12 + 2.0e12 * psi,
            "preprocprof/Te": 100.0 - 60.0 * psi,
            "preprocprof/Ti": 8.0 - 4.0 * psi,
            "preprocprof/Vz": 2000.0 - 20000.0 * psi,
            "preprocprof/q": -2.0 - 3.0 * psi,
            "preprocprof/equil/r": np.array([1.0, 3.0, 5.0]),
            "preprocprof/equil/psi_pol_norm": np.array([0.0, 0.5, 1.0]),
            "preprocprof/equil/q": np.array([2.0, 3.5, 5.0]),
        }
        for dataset_path, values in datasets.items():
            handle.create_dataset(dataset_path, data=values)


def test_synthetic_balance_to_ql_balance_acceptance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    launched: list[object] = []

    def unexpected_process(*args: object, **kwargs: object) -> None:
        launched.append((args, kwargs))
        raise AssertionError("synthetic acceptance must not launch a process")

    monkeypatch.setattr(preparation_module.subprocess, "run", unexpected_process)
    monkeypatch.setattr(
        preparation_module.executable_module,
        "discover_kamel_git_metadata",
        lambda: None,
    )

    balance_paths = _write_balance_profiles(tmp_path / "balance")
    source_bytes = {role: path.read_bytes() for role, path in balance_paths.items()}
    source = _read_balance_case(balance_paths)
    assert len({tuple(profile.coordinate) for profile in source.profiles.values()}) == 4

    staged = stage_balance_marsf_quartet(
        source,
        tmp_path / "staged-marsf",
        major_radius_cm=_R0_CM,
        equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
    )
    marsf = read_marsf_profiles(staged.directory, staged.metadata)
    assert {path.name for path in staged.directory.iterdir()} == {
        *_MARSF_FILENAMES.values(),
        "staging_report.json",
    }
    equilibrium = tmp_path / "equilibrium.dat"
    _write_equilibrium(equilibrium)
    equilibrium_bytes = equilibrium.read_bytes()

    prepared = prepare_marsf_case(
        marsf,
        _config(),
        tmp_path / "prepared",
        equilibrium_file=equilibrium,
        q_operation="negate",
        upstream_staging_report=staged.report,
    )

    assert launched == []
    assert {role: path.read_bytes() for role, path in balance_paths.items()} == source_bytes
    assert equilibrium.read_bytes() == equilibrium_bytes

    staged_report = json.loads(staged.report.read_text(encoding="utf-8"))
    assert staged_report["source_hashes"] == {
        role: hashlib.sha256(path.read_bytes()).hexdigest() for role, path in balance_paths.items()
    }
    assert staged_report["derived_hashes"] == {
        filename: _sha256(staged.directory / filename) for filename in _MARSF_FILENAMES.values()
    }
    assert staged_report["coordinate_mapping"] == {
        "source": "rho_pol",
        "source_unit": "1",
        "target": "sqrt_psiN",
        "target_unit": "1",
        "operation": "rho_pol = sqrt(psi_pol_norm)",
    }
    assert staged_report["operations"][-1] == {
        "quantity": "toroidal_velocity",
        "source_unit": "rad/s",
        "target_unit": "cm/s",
        "factor": _R0_CM,
        "operation": "omega_to_v_phi",
        "parameters": {"major_radius_cm": _R0_CM},
    }

    np.testing.assert_allclose(
        marsf.profiles["toroidal_velocity"].values,
        np.array([2000.0, -1200.0, -3000.0, -18000.0]),
        rtol=0.0,
        atol=1.0e-9,
    )

    prepared_values = {
        name: np.loadtxt(prepared.profiles / filename)[:, 1]
        for name, filename in {
            "n": "n.dat",
            "Te": "Te.dat",
            "Ti": "Ti.dat",
            "Vz": "Vz.dat",
            "q": "q.dat",
        }.items()
    }
    prepared_radius = np.loadtxt(prepared.profiles / "n.dat")[:, 0]
    np.testing.assert_array_equal(prepared_radius, [1.0, 3.0, 5.0])
    np.testing.assert_array_equal(prepared_values["n"], [1.0e12, 2.0e12, 3.0e12])
    np.testing.assert_array_equal(prepared_values["Te"], [100.0, 70.0, 40.0])
    np.testing.assert_array_equal(prepared_values["Ti"], [8.0, 6.0, 4.0])
    np.testing.assert_array_equal(prepared_values["Vz"], [2000.0, -8000.0, -18000.0])
    np.testing.assert_array_equal(prepared_values["q"], [-2.0, -3.5, -5.0])

    prepared_report = json.loads(prepared.report.read_text(encoding="utf-8"))
    assert prepared_report["coordinate_operation"] == {
        "source_coordinate": "sqrt_psiN",
        "target_coordinate": "r_eff",
        "method": "natural cubic interpolation",
    }
    assert prepared_report["equilibrium"]["source_hash"] == _sha256(equilibrium)
    assert prepared_report["upstream_staging"]["provenance_sha256"] == _sha256(staged.report)
    assert prepared_report["upstream_staging"]["source_hashes"] == staged_report["source_hashes"]
    assert prepared_report["upstream_staging"]["derived_hashes"] == staged_report["derived_hashes"]
    assert prepared_report["source_hashes"] == {
        **{
            f"source/{filename}": staged_report["derived_hashes"][filename]
            for filename in _MARSF_FILENAMES.values()
        },
        f"equilibrium/{equilibrium.name}": _sha256(equilibrium),
    }
    assert prepared_report["prepared_output_hashes"] == {
        f"profiles/{filename}": _sha256(prepared.profiles / filename)
        for filename in ("n.dat", "Te.dat", "Ti.dat", "Vz.dat", "q.dat")
    }
    density_operation = next(
        operation
        for operation in prepared_report["operations"]
        if operation["quantity"] == "density"
    )
    assert density_operation == {
        "quantity": "density",
        "source_unit": "1/m^3",
        "target_unit": "1/cm^3",
        "factor": 1.0e-6,
    }
    q_operation = next(
        operation for operation in prepared_report["operations"] if operation["quantity"] == "q"
    )
    assert q_operation["operation"] == "negate"
    assert q_operation["factor"] == -1.0

    with pytest.raises(ExperimentalInputError, match="already exists"):
        stage_balance_marsf_quartet(
            source,
            staged.directory,
            major_radius_cm=_R0_CM,
            equilibrium_provenance=_EQUILIBRIUM_PROVENANCE,
        )
    with pytest.raises(ExperimentalInputError, match="already exists"):
        prepare_marsf_case(
            marsf,
            _config(),
            prepared.directory,
            equilibrium_file=equilibrium,
            q_operation="negate",
            upstream_staging_report=staged.report,
        )
    assert staged.report.is_file()
    assert prepared.report.is_file()

    oracle_path = tmp_path / "ql-balance-oracle.h5"
    _write_ql_balance_oracle(oracle_path)
    oracle = read_ql_balance_oracle(oracle_path)
    assert oracle.r_out.flags.writeable is False

    oracle_profiles = {
        "n": oracle.n,
        "Te": oracle.Te,
        "Ti": oracle.Ti,
        "Vz": oracle.Vz,
        "q": oracle.q,
    }
    floors = {"n": 1.0e12, "Te": 1.0, "Ti": 1.0, "Vz": 1.0, "q": 0.1}
    for name, oracle_values in oracle_profiles.items():
        comparison = compare_profiles(
            oracle_radius_cm=oracle.r_out,
            oracle_values=oracle_values,
            prepared_radius_cm=prepared_radius,
            prepared_values=prepared_values[name],
            domain_cm=(2.0, 4.0),
            interpolation_direction="prepared_to_oracle",
            method="linear",
            relative_floor=floors[name],
            tolerances=None,
            resonance=(7, 2) if name == "q" else None,
        )
        assert comparison.threshold_decisions is None
        assert comparison.overall_pass is None
        assert comparison.comparison_radius_cm.tolist() == [2.0, 3.0, 4.0]
        assert comparison.measurements.absolute_max == 0.0
        assert comparison.measurements.absolute_rms == 0.0
        assert comparison.measurements.relative_max == 0.0
        assert comparison.measurements.relative_rms == 0.0
        assert comparison.exclusions.oracle.outside_domain_points_cm == (1.0, 5.0)
        assert comparison.exclusions.prepared.outside_domain_points_cm == (1.0, 5.0)
        intervals = comparison.exclusions.oracle.outside_domain_intervals_cm
        assert tuple(
            (item.lower_cm, item.upper_cm, item.lower_inclusive, item.upper_inclusive)
            for item in intervals
        ) == ((1.0, 2.0, True, False), (4.0, 5.0, False, True))
        assert comparison.exclusions.prepared.outside_domain_intervals_cm == intervals
        assert comparison.exclusions.oracle.outside_shared_overlap_points_cm == ()
        assert comparison.exclusions.prepared.outside_shared_overlap_points_cm == ()
        assert comparison.exclusions.oracle.outside_shared_overlap_intervals_cm == ()
        assert comparison.exclusions.prepared.outside_shared_overlap_intervals_cm == ()
        assert comparison.warnings == ("domain-exclusions",)
        if name == "q":
            assert comparison.resonance is not None
            assert comparison.resonance.target_q == -3.5
            assert comparison.resonance.reference_crossing_radii_cm == (3.0,)
            assert comparison.resonance.candidate_crossing_radii_cm == (3.0,)
            assert comparison.resonance.reference_crossing_covered is True
            assert comparison.resonance.candidate_crossing_covered is True
        else:
            assert comparison.resonance is None
