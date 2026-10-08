from __future__ import annotations

import hashlib
import json
import os
import shutil
import stat
from dataclasses import replace
from pathlib import Path

import h5py
import numpy as np
import pytest
from kim import BuiltinPlasma, PlasmaIsotope, SimulationConfig
from kim.balance_characterization import BalanceCharacterizationRequest, characterize_balance
from kim.cli import app
from kim.errors import ExperimentalInputError
from typer.testing import CliRunner

runner = CliRunner()


def test_destination_rejects_dangling_final_and_parent_symlinks(tmp_path: Path) -> None:
    from kim.balance_characterization import _resolve_destination

    dangling_target = tmp_path / "missing-target"
    dangling_parent = tmp_path / "dangling-parent"
    dangling_parent.symlink_to(dangling_target, target_is_directory=True)
    with pytest.raises(ExperimentalInputError, match="dangling symlink"):
        _resolve_destination(dangling_parent / "output")

    dangling_destination = tmp_path / "dangling-output"
    dangling_destination.symlink_to(dangling_target, target_is_directory=True)
    with pytest.raises(ExperimentalInputError, match="dangling symlink"):
        _resolve_destination(dangling_destination)

    valid_parent = tmp_path / "valid-parent"
    valid_parent_target = tmp_path / "valid-parent-target"
    valid_parent_target.mkdir()
    valid_parent.symlink_to(valid_parent_target, target_is_directory=True)
    assert _resolve_destination(valid_parent / "output") == valid_parent_target / "output"


def test_missing_destination_parent_is_unavailable_and_anchor_failures_close_fds(
    tmp_path: Path,
) -> None:
    from kim.balance_characterization import _open_directory_anchor, _resolve_destination

    missing_parent = tmp_path / "not-created" / "output"
    with pytest.raises(ExperimentalInputError, match="must already exist"):
        _resolve_destination(missing_parent)

    missing_anchor = tmp_path / "missing-anchor" / "child"
    before = len(os.listdir("/dev/fd"))
    for _ in range(32):
        with pytest.raises(ExperimentalInputError, match="unable to anchor"):
            _open_directory_anchor(missing_anchor)
    after = len(os.listdir("/dev/fd"))
    assert after <= before + 1


def test_atomic_publish_rejects_replaced_private_entry(tmp_path: Path) -> None:
    from kim.balance_adoption import _publish_staging_directory_at

    staging = tmp_path / "staging"
    staging.mkdir()
    original_inode = staging.stat()
    foreign = tmp_path / "foreign"
    foreign.mkdir()
    sentinel = foreign / "foreign-sentinel"
    sentinel.write_text("keep")
    backup = tmp_path / "staging-backup"
    parent_fd = os.open(tmp_path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        staging.rename(backup)
        staging.symlink_to(foreign, target_is_directory=True)
        with pytest.raises(ExperimentalInputError, match="private staging entry changed"):
            _publish_staging_directory_at(
                parent_fd,
                staging.name,
                parent_fd,
                "destination",
                expected_source_inode=(original_inode.st_dev, original_inode.st_ino),
            )
    finally:
        os.close(parent_fd)
    assert not (tmp_path / "destination").exists()
    assert sentinel.read_text() == "keep"
    assert list(foreign.iterdir()) == [sentinel]
    staging.unlink()
    backup.rmdir()


def test_private_container_anchor_failure_retains_only_owned_container(
    tmp_path: Path, monkeypatch
) -> None:
    import kim.balance_characterization as characterization_module

    parent_fd = os.open(tmp_path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    original_validate = characterization_module._validate_private_entry

    def fail_anchor(parent: int, name: str, descriptor: int) -> Path:
        raise ExperimentalInputError("synthetic tuple-assignment anchor failure")

    monkeypatch.setattr(characterization_module, "_validate_private_entry", fail_anchor)
    try:
        with pytest.raises(ExperimentalInputError, match="private characterization container"):
            characterization_module._create_private_output(parent_fd, "tuple-anchor")
    finally:
        os.close(parent_fd)
        monkeypatch.setattr(characterization_module, "_validate_private_entry", original_validate)

    retained = list(tmp_path.glob(".tuple-anchor-*"))
    assert len(retained) == 1
    assert (retained[0] / "payload").is_dir()


@pytest.mark.parametrize(
    ("value", "message"),
    [
        (True, "real"),
        ("2", "real"),
        (float("nan"), "finite"),
        (4.0, "non-empty"),
    ],
)
def test_domains_require_strict_finite_real_nonempty_endpoints(value: object, message: str) -> None:
    from kim.balance_characterization import _normalise_domains

    with pytest.raises(ExperimentalInputError, match=message):
        _normalise_domains({"core": (value, 4.0)})


@pytest.mark.parametrize("value", [True, "1.0"])
def test_floors_and_tolerances_reject_bool_and_string_numbers(value: object) -> None:
    from kim.balance_characterization import _metrics, _normalise_floors

    with pytest.raises(ExperimentalInputError, match="floor"):
        _normalise_floors({"n": value, "Te": 1.0, "Ti": 1.0, "Vz": 1.0, "q": 1.0})
    with pytest.raises(ExperimentalInputError, match="tolerance"):
        _metrics({"absolute_max": value})


def test_tolerance_profile_aliases_cannot_overwrite_one_semantic_profile() -> None:
    from kim.balance_characterization import _normalise_tolerances

    with pytest.raises(ExperimentalInputError, match="more than once"):
        _normalise_tolerances(
            {"n": {"absolute_max": 1.0}, "density": {"absolute_max": 2.0}},
            {"core": (1.0, 4.0)},
        )


@pytest.mark.parametrize("value", [True, "165"])
def test_major_radius_and_timeout_use_strict_real_validation(value: object) -> None:
    from kim.balance_characterization import _strict_real

    with pytest.raises(ExperimentalInputError):
        _strict_real(value, "major_radius_cm")
    with pytest.raises(ExperimentalInputError):
        _strict_real(value, "equilibrium_timeout_seconds")


def test_missing_required_input_is_unavailable_without_creating_output(tmp_path: Path) -> None:
    request = BalanceCharacterizationRequest(
        density=tmp_path / "missing-n.dat",
        electron_temperature=tmp_path / "missing-te.dat",
        ion_temperature=tmp_path / "missing-ti.dat",
        toroidal_rotation=tmp_path / "missing-vt.dat",
        metadata={
            "source": "synthetic",
            "coordinate": "rho_pol",
            "coordinate_unit": "1",
            "density_unit": "1/m^3",
            "electron_temperature_unit": "eV",
            "ion_temperature_unit": "eV",
            "toroidal_rotation_unit": "rad/s",
        },
        equilibrium_provenance="synthetic-equilibrium",
        equilibrium_file=tmp_path / "missing-equilibrium.dat",
        equilibrium_parameters_file=tmp_path / "btor_rbig.dat",
        config=tmp_path / "missing-request.json",
        oracle=tmp_path / "missing-oracle.h5",
        destination=tmp_path / "characterization",
        domains={"all": (1.0, 2.0)},
        relative_floors={"n": 1.0, "Te": 1.0, "Ti": 1.0, "Vz": 1.0, "q": 1.0},
        interpolation_direction="prepared_to_oracle",
        interpolation_method="linear",
        q_operation="preserve",
    )

    result = characterize_balance(request)

    assert result.status == "UNAVAILABLE"
    assert result.report["status"] == "UNAVAILABLE"
    assert result.report["limitations"] == {
        "input_adoption_only": True,
        "no_physical_response_validation": True,
        "solvers_not_run": True,
        "retained_staging_possible": True,
        "retained_staging_success_state": "empty_container",
        "retained_staging_failure_state": "incomplete_container_possible",
    }
    assert result.report["retained_staging"] == {
        "path": None,
        "reason": "private staging path unavailable: no retained descriptor",
        "state": "not_created",
        "cleanup_safe_when_no_characterization_is_running": True,
    }
    assert not request.destination.exists()


def test_synthetic_characterization_measures_without_running_a_solver(
    tmp_path: Path, monkeypatch
) -> None:
    sources = tmp_path / "sources"
    balance = sources / "balance"
    balance.mkdir(parents=True)
    grids = {
        "density": [0.0, 0.5, 1.0],
        "electron_temperature": [0.0, 0.5, 1.0],
        "ion_temperature": [0.0, 0.5, 1.0],
        "toroidal_rotation": [0.0, 0.5, 1.0],
    }
    values = {
        "density": lambda rho: 1.0e18 + 2.0e18 * rho * rho,
        "electron_temperature": lambda rho: 100.0 - 60.0 * rho * rho,
        "ion_temperature": lambda rho: 8.0 - 4.0 * rho * rho,
        "toroidal_rotation": lambda rho: 10.0 - 100.0 * rho * rho,
    }
    paths = {}
    for role, grid in grids.items():
        path = balance / f"{role}.dat"
        path.write_text("\n".join(f"{rho} {values[role](rho)}" for rho in grid) + "\n")
        paths[role] = path
    equilibrium = sources / "equilibrium" / "equil_r_q_psi.dat"
    equilibrium.parent.mkdir()
    equilibrium.write_text("# r q psi\n1 2 0\n3 3.5 .5\n5 5 1\n")
    equilibrium_parameters = equilibrium.parent / "btor_rbig.dat"
    equilibrium_parameters.write_text("-17573.19212 200.0\n", encoding="utf-8")
    oracle_path = sources / "oracle" / "oracle.h5"
    oracle_path.parent.mkdir()
    radius = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    psi = (radius - 1.0) / 4.0
    with h5py.File(oracle_path, "w") as handle:
        for dataset, array in {
            "preprocprof/r_out": radius,
            "preprocprof/n": 1.0e12 + 2.0e12 * psi,
            "preprocprof/Te": 100.0 - 60.0 * psi,
            "preprocprof/Ti": 8.0 - 4.0 * psi,
            "preprocprof/Vz": 2000.0 - 20000.0 * psi,
            "preprocprof/q": 2.0 + 3.0 * psi,
            "preprocprof/equil/r": np.array([1.0, 3.0, 5.0]),
            "preprocprof/equil/psi_pol_norm": np.array([0.0, 0.5, 1.0]),
            "preprocprof/equil/q": np.array([2.0, 3.5, 5.0]),
        }.items():
            handle.create_dataset(dataset, data=array)
    config = SimulationConfig.electrostatic_periodic(
        profiles=Path("unused"),
        plasma=BuiltinPlasma(isotope=PlasmaIsotope.DEUTERIUM),
        btor=-18000.0,
        major_radius=165.0,
        m_mode=7,
        n_mode=2,
        frequency=0.0,
        br_boundary_real=1.0,
        br_boundary_imag=0.0,
        radial_minimum=1.0,
        plasma_radius=5.0,
    )
    config = config.model_copy(
        update={
            "profiles": config.profiles.model_copy(
                update={
                    "density_file": "custom-density.dat",
                    "electron_temperature_file": "custom-te.dat",
                    "ion_temperature_file": "custom-ti.dat",
                    "toroidal_velocity_file": "custom-vz.dat",
                    "safety_factor_file": "custom-q.dat",
                }
            )
        }
    )
    source_bytes = {key: path.read_bytes() for key, path in paths.items()}
    monkeypatch.setattr(
        "kim.preparation.subprocess.run",
        lambda *a, **k: (_ for _ in ()).throw(AssertionError("solver launched")),
    )
    request = BalanceCharacterizationRequest(
        density=paths["density"],
        electron_temperature=paths["electron_temperature"],
        ion_temperature=paths["ion_temperature"],
        toroidal_rotation=paths["toroidal_rotation"],
        metadata={
            "source": "synthetic",
            "coordinate": "rho_pol",
            "coordinate_unit": "1",
            "density_unit": "1/m^3",
            "electron_temperature_unit": "eV",
            "ion_temperature_unit": "eV",
            "toroidal_rotation_unit": "rad/s",
        },
        equilibrium_provenance="synthetic-equilibrium",
        equilibrium_file=equilibrium,
        equilibrium_parameters_file=equilibrium_parameters,
        config=config,
        oracle=oracle_path,
        destination=tmp_path / "characterization",
        domains={"core": (2.0, 4.0)},
        relative_floors={"n": 1.0e12, "Te": 1.0, "Ti": 1.0, "Vz": 1.0, "q": 0.1},
        interpolation_direction="prepared_to_oracle",
        interpolation_method="linear",
        major_radius_cm=200.0,
    )
    from kim.balance_characterization import _normalise_request

    with pytest.raises(ExperimentalInputError, match="preserves Fouriers q without a sign change"):
        _normalise_request(replace(request, q_operation="negate"))
    assert request.q_operation == "preserve"
    assert request.q_operation == "preserve"

    positional_legacy_request = BalanceCharacterizationRequest(
        request.density,
        request.electron_temperature,
        request.ion_temperature,
        request.toroidal_rotation,
        request.metadata,
        request.equilibrium_provenance,
        request.config,
        request.oracle,
        request.destination,
        200.0,
        "preserve",
        request.domains,
        request.relative_floors,
        request.interpolation_direction,
        request.interpolation_method,
        equilibrium_file=request.equilibrium_file,
        equilibrium_parameters_file=request.equilibrium_parameters_file,
    )
    assert _normalise_request(positional_legacy_request)["major_radius_cm"] == 200.0

    unrelated_equilibrium_directory = sources / "unrelated-equilibrium"
    unrelated_equilibrium_directory.mkdir()
    unrelated_parameters = unrelated_equilibrium_directory / "btor_rbig.dat"
    unrelated_parameters.write_bytes(equilibrium_parameters.read_bytes())
    with pytest.raises(ExperimentalInputError, match="same equilibrium calculation directory"):
        _normalise_request(replace(request, equilibrium_parameters_file=unrelated_parameters))

    renamed_equilibrium = equilibrium.with_name("other-equilibrium-table.dat")
    renamed_equilibrium.write_bytes(equilibrium.read_bytes())
    with pytest.raises(ExperimentalInputError, match="equil_r_q_psi.dat"):
        _normalise_request(replace(request, equilibrium_file=renamed_equilibrium))

    cwd_before = Path.cwd()
    fd_before = len(os.listdir("/dev/fd"))
    result = characterize_balance(request)
    fd_after = len(os.listdir("/dev/fd"))

    assert result.status == "MEASURED"
    assert Path.cwd() == cwd_before
    assert fd_after <= fd_before + 1
    assert result.report_path == request.destination / "characterization_report.json"
    assert result.report["threshold_decisions"] is None
    assert result.report["overall_pass"] is None
    floor_free = characterize_balance(
        replace(
            request,
            destination=tmp_path / "floor-free-characterization",
            relative_floors=None,
        )
    )
    assert floor_free.status == "MEASURED"
    assert floor_free.report["comparison_configuration"]["relative_floors"] is None
    assert floor_free.report["threshold_decisions"] is None
    assert floor_free.report["overall_pass"] is None
    for profiles in floor_free.report["comparisons"].values():
        for comparison in profiles.values():
            measurements = comparison["measurements"]
            assert np.isfinite(measurements["absolute_rms"])
            assert np.isfinite(measurements["absolute_max"])
            assert measurements["relative_rms"] is None
            assert measurements["relative_max"] is None
            assert any(
                "relative metrics unavailable" in warning for warning in comparison["warnings"]
            )
            assert comparison["overall_pass"] is None
    relative_tolerance_without_floors = characterize_balance(
        replace(
            request,
            destination=tmp_path / "relative-tolerance-without-floors",
            relative_floors=None,
            tolerances={"core": {"n": {"relative_rms": 0.1}}},
        )
    )
    assert relative_tolerance_without_floors.status == "UNAVAILABLE"
    assert (
        "relative tolerances require relative_floors"
        in relative_tolerance_without_floors.report["error"]["message"]
    )
    assert not (tmp_path / "relative-tolerance-without-floors").exists()
    assert result.report["expected_hashes_requested"] is False
    assert json.loads(Path(result.report["staged"]["report"]).read_text())["schema_version"] == 2
    assert result.report["equilibrium_calculation"]["btor_gauss"] == pytest.approx(-17573.19212)
    assert result.report["equilibrium_calculation"]["r_big_cm"] == pytest.approx(200.0)
    equilibrium_report = result.report["equilibrium_calculation"]
    assert (
        equilibrium_report["equilibrium_sha256"]
        == hashlib.sha256(Path(equilibrium_report["equilibrium_file"]).read_bytes()).hexdigest()
    )
    assert (
        equilibrium_report["parameters_sha256"]
        == hashlib.sha256(Path(equilibrium_report["parameters_file"]).read_bytes()).hexdigest()
    )
    prepared_request = json.loads(Path(result.report["prepared"]["request"]).read_text())
    assert prepared_request["setup"]["btor"] == pytest.approx(-17573.19212)
    assert prepared_request["setup"]["major_radius"] == pytest.approx(200.0)
    staged_rotation = np.loadtxt(
        Path(result.report["staged"]["directory"]) / "PROFROT.IN", skiprows=1
    )
    np.testing.assert_allclose(
        staged_rotation[:, 1],
        np.loadtxt(paths["toroidal_rotation"])[:, 1] * 200.0,
    )
    assert result.report["comparisons"]["core"]["q"]["resonance"]["target_q"] == -3.5
    assert result.report["comparisons"]["core"]["n"]["role"]["prepared"] == "custom-density.dat"
    assert "custom-density.dat" in result.report["prepared"]["profile_hashes"]
    assert result.report_path.is_file()
    assert {key: path.read_bytes() for key, path in paths.items()} == source_bytes
    retained = result.report["retained_staging"]
    retained_path = Path(retained["path"])
    assert retained["reason"] == "safe ownership policy / portable rename limitation"
    assert retained["state"] == "empty_container"
    assert retained["cleanup_safe_when_no_characterization_is_running"] is True
    assert retained_path.is_dir()
    assert stat.S_IMODE(retained_path.stat().st_mode) == 0o700
    assert list(retained_path.iterdir()) == []
    assert retained_path.parent == request.destination.parent
    assert not retained_path.is_relative_to(request.destination)
    for source_root in (balance, equilibrium.parent, oracle_path.parent):
        assert not retained_path.is_relative_to(source_root)
    persisted = json.loads(result.report_path.read_text(encoding="utf-8"))
    assert persisted["retained_staging"] == retained

    mismatched_legacy_radius = characterize_balance(
        replace(
            request,
            destination=tmp_path / "mismatched-legacy-radius",
            major_radius_cm=201.0,
        )
    )
    assert mismatched_legacy_radius.status == "UNAVAILABLE"
    assert (
        "deprecated major_radius_cm does not match"
        in mismatched_legacy_radius.report["error"]["message"]
    )
    assert not (tmp_path / "mismatched-legacy-radius").exists()

    from kim import balance_characterization as characterization_module

    anchor_parent = tmp_path / "anchor-parent"
    anchor_parent.mkdir()
    redirected_parent = tmp_path / "redirected-parent"
    redirected_parent.mkdir()
    redirected_sentinel = redirected_parent / "foreign-sentinel"
    redirected_sentinel.write_text("keep")
    original_validate_anchor = characterization_module._validate_parent_anchor
    anchor_checks = 0
    anchor_backup = tmp_path / "anchor-parent-backup"

    def swap_parent_after_private_reservation(path: Path, descriptor: int) -> None:
        nonlocal anchor_checks
        original_validate_anchor(path, descriptor)
        anchor_checks += 1
        if anchor_checks == 2:
            anchor_parent.rename(anchor_backup)
            anchor_parent.symlink_to(redirected_parent, target_is_directory=True)

    monkeypatch.setattr(
        characterization_module,
        "_validate_parent_anchor",
        swap_parent_after_private_reservation,
    )
    swapped_parent_destination = anchor_parent / "characterization"
    swapped_parent = characterize_balance(replace(request, destination=swapped_parent_destination))
    assert swapped_parent.status == "UNAVAILABLE"
    assert anchor_checks == 2
    assert redirected_sentinel.read_text() == "keep"
    assert list(redirected_parent.iterdir()) == [redirected_sentinel]
    assert not (redirected_parent / "characterization").exists()
    retained_containers = list(anchor_backup.glob(".characterization-*"))
    assert len(retained_containers) == 1
    assert (retained_containers[0] / "payload").is_dir()
    anchor_parent.unlink()
    anchor_backup.rename(anchor_parent)
    monkeypatch.setattr(
        characterization_module,
        "_validate_parent_anchor",
        original_validate_anchor,
    )

    outer_swap_foreign = tmp_path / "outer-container-foreign"
    outer_swap_foreign.mkdir()
    outer_swap_sentinel = outer_swap_foreign / "foreign-sentinel"
    outer_swap_sentinel.write_text("keep")
    outer_swap_backup = tmp_path / "outer-container-backup"
    outer_anchor_checks = 0

    def replace_outer_container(path: Path, descriptor: int) -> None:
        nonlocal outer_anchor_checks
        original_validate_anchor(path, descriptor)
        outer_anchor_checks += 1
        if outer_anchor_checks == 2:
            outer_container = next(tmp_path.glob(".outer-container-characterization-*"))
            outer_container.rename(outer_swap_backup)
            outer_container.symlink_to(outer_swap_foreign, target_is_directory=True)

    monkeypatch.setattr(
        characterization_module,
        "_validate_parent_anchor",
        replace_outer_container,
    )
    outer_swapped_destination = tmp_path / "outer-container-characterization"
    outer_swapped = characterize_balance(replace(request, destination=outer_swapped_destination))
    assert outer_swapped.status == "MEASURED"
    assert outer_swapped_destination.is_dir()
    assert outer_swap_sentinel.read_text() == "keep"
    assert list(outer_swap_foreign.iterdir()) == [outer_swap_sentinel]
    assert not (outer_swap_backup / "payload").exists()
    outer_container_link = next(tmp_path.glob(".outer-container-characterization-*"))
    assert outer_container_link.is_symlink()
    outer_container_link.unlink()
    outer_swap_backup.rmdir()
    monkeypatch.setattr(
        characterization_module,
        "_validate_parent_anchor",
        original_validate_anchor,
    )

    request_files = sources / "request"
    request_files.mkdir()
    metadata_file = request_files / "metadata.json"
    metadata_file.write_text(
        '{"source":"synthetic","coordinate":"rho_pol","coordinate_unit":"1",'
        '"density_unit":"1/m^3","electron_temperature_unit":"eV",'
        '"ion_temperature_unit":"eV","toroidal_rotation_unit":"rad/s"}\n'
    )
    config_file = request_files / "config.json"
    config_file.write_text(config.model_dump_json() + "\n")
    path_request = replace(request, metadata=metadata_file, config=config_file)
    conflicting_metadata = request_files / "other-metadata.json"
    conflicting_metadata.write_bytes(metadata_file.read_bytes())
    conflicting_config = request_files / "other-config.json"
    conflicting_config.write_bytes(config_file.read_bytes())
    for label, conflicting_request in (
        (
            "metadata",
            replace(
                path_request,
                metadata_path=conflicting_metadata,
                destination=tmp_path / "conflicting-metadata-characterization",
            ),
        ),
        (
            "config",
            replace(
                path_request,
                config_path=conflicting_config,
                destination=tmp_path / "conflicting-config-characterization",
            ),
        ),
    ):
        conflicting_result = characterize_balance(conflicting_request)
        assert conflicting_result.status == "UNAVAILABLE", label
        assert not conflicting_request.destination.exists()
    request_source_destination = request_files / "nested-output"
    request_source_result = characterize_balance(
        replace(path_request, destination=request_source_destination)
    )
    assert request_source_result.status == "UNAVAILABLE"
    assert not request_source_destination.exists()
    path_result = characterize_balance(
        replace(path_request, destination=tmp_path / "path-characterization")
    )
    assert path_result.status == "MEASURED"
    assert (
        path_result.report["actual_sha256"]["metadata"]
        == hashlib.sha256(metadata_file.read_bytes()).hexdigest()
    )
    assert path_result.report["request_provenance"]["metadata"]["mode"] == "file_snapshot"

    original_validate_anchor = characterization_module._validate_parent_anchor
    snapshot_swap_foreign = tmp_path / "private-entry-snapshot-foreign"
    snapshot_swap_foreign.mkdir()
    snapshot_swap_sentinel = snapshot_swap_foreign / "foreign-sentinel"
    snapshot_swap_sentinel.write_text("keep")
    snapshot_swap_backup_name = "payload-backup"
    snapshot_anchor_checks = 0

    def replace_private_before_snapshot(path: Path, descriptor: int) -> None:
        nonlocal snapshot_anchor_checks
        original_validate_anchor(path, descriptor)
        snapshot_anchor_checks += 1
        if snapshot_anchor_checks == 2:
            private_container = next(tmp_path.glob(".private-entry-snapshot-characterization-*"))
            payload = private_container / "payload"
            payload.rename(private_container / snapshot_swap_backup_name)
            payload.symlink_to(snapshot_swap_foreign, target_is_directory=True)

    monkeypatch.setattr(
        characterization_module,
        "_validate_parent_anchor",
        replace_private_before_snapshot,
    )
    snapshot_swapped = characterize_balance(
        replace(request, destination=tmp_path / "private-entry-snapshot-characterization")
    )
    assert snapshot_swapped.status == "UNAVAILABLE"
    assert snapshot_anchor_checks == 2
    assert snapshot_swap_sentinel.read_text() == "keep"
    assert list(snapshot_swap_foreign.iterdir()) == [snapshot_swap_sentinel]
    snapshot_container = next(tmp_path.glob(".private-entry-snapshot-characterization-*"))
    snapshot_payload = snapshot_container / "payload"
    assert snapshot_payload.is_symlink()
    snapshot_payload.unlink()
    (snapshot_container / snapshot_swap_backup_name).rmdir()
    monkeypatch.setattr(
        characterization_module,
        "_validate_parent_anchor",
        original_validate_anchor,
    )

    original_publish = characterization_module._publish_staging_directory_at
    concurrent_destination = tmp_path / "concurrent-characterization"

    def inject_concurrent_winner(
        source_fd: int,
        source_name: str,
        destination_fd: int,
        destination: str,
        **kwargs: object,
    ) -> None:
        destination_path = characterization_module._descriptor_path(destination_fd) / destination
        destination_path.mkdir()
        (destination_path / "foreign-sentinel").write_text("keep")
        original_publish(source_fd, source_name, destination_fd, destination, **kwargs)

    monkeypatch.setattr(
        characterization_module,
        "_publish_staging_directory_at",
        inject_concurrent_winner,
    )
    concurrent = characterize_balance(replace(request, destination=concurrent_destination))
    assert concurrent.status == "UNAVAILABLE"
    assert (concurrent_destination / "foreign-sentinel").read_text() == "keep"
    concurrent_containers = list(tmp_path.glob(f".{concurrent_destination.name}-*"))
    assert len(concurrent_containers) == 1
    assert (concurrent_containers[0] / "payload").is_dir()
    assert Path(concurrent.report["retained_staging"]["path"]) == concurrent_containers[0]
    assert concurrent.report["retained_staging"]["state"] == "incomplete_container"

    def fail_publish(
        source_fd: int,
        source_name: str,
        destination_fd: int,
        destination: str,
        **kwargs: object,
    ) -> None:
        raise OSError("synthetic publication failure")

    monkeypatch.setattr(characterization_module, "_publish_staging_directory_at", fail_publish)
    failed_publication_destination = tmp_path / "failed-publication-characterization"
    failed_publication = characterize_balance(
        replace(request, destination=failed_publication_destination)
    )
    assert failed_publication.status == "UNAVAILABLE"
    assert not failed_publication_destination.exists()
    failed_containers = list(tmp_path.glob(f".{failed_publication_destination.name}-*"))
    assert len(failed_containers) == 1
    assert (failed_containers[0] / "payload").is_dir()
    assert Path(failed_publication.report["retained_staging"]["path"]) == failed_containers[0]
    assert failed_publication.report["retained_staging"]["state"] == "incomplete_container"

    original_verify_inputs = characterization_module._verify_unchanged_inputs
    publish_swap_foreign = tmp_path / "private-entry-publish-foreign"
    publish_swap_foreign.mkdir()
    publish_swap_sentinel = publish_swap_foreign / "foreign-sentinel"
    publish_swap_sentinel.write_text("keep")
    publish_swap_backup_name = "payload-backup"
    verify_calls = 0

    def replace_private_before_publish(data: object, before: object) -> object:
        nonlocal verify_calls
        result = original_verify_inputs(data, before)
        verify_calls += 1
        if verify_calls == 3:
            private_container = next(tmp_path.glob(".private-entry-publish-characterization-*"))
            payload = private_container / "payload"
            payload.rename(private_container / publish_swap_backup_name)
            payload.symlink_to(publish_swap_foreign, target_is_directory=True)
        return result

    monkeypatch.setattr(
        characterization_module,
        "_verify_unchanged_inputs",
        replace_private_before_publish,
    )
    publish_swapped_destination = tmp_path / "private-entry-publish-characterization"
    publish_swapped = characterize_balance(
        replace(request, destination=publish_swapped_destination)
    )
    assert publish_swapped.status == "UNAVAILABLE"
    assert verify_calls == 3
    assert not publish_swapped_destination.exists()
    assert publish_swap_sentinel.read_text() == "keep"
    assert list(publish_swap_foreign.iterdir()) == [publish_swap_sentinel]
    publish_container = next(tmp_path.glob(".private-entry-publish-characterization-*"))
    publish_payload = publish_container / "payload"
    assert publish_payload.is_symlink()
    publish_payload.unlink()
    shutil.rmtree(publish_container / publish_swap_backup_name)
    monkeypatch.setattr(
        characterization_module,
        "_verify_unchanged_inputs",
        original_verify_inputs,
    )

    post_anchor_parent = tmp_path / "post-anchor-parent"
    post_anchor_parent.mkdir()
    post_redirected_parent = tmp_path / "post-redirected-parent"
    post_redirected_parent.mkdir()
    post_sentinel = post_redirected_parent / "foreign-sentinel"
    post_sentinel.write_text("keep")
    post_anchor_backup = tmp_path / "post-anchor-parent-backup"

    def publish_then_swap_parent(
        source_fd: int,
        source_name: str,
        destination_fd: int,
        destination: str,
        **kwargs: object,
    ) -> None:
        original_publish(source_fd, source_name, destination_fd, destination, **kwargs)
        post_anchor_parent.rename(post_anchor_backup)
        post_anchor_parent.symlink_to(post_redirected_parent, target_is_directory=True)

    monkeypatch.setattr(
        characterization_module,
        "_publish_staging_directory_at",
        publish_then_swap_parent,
    )
    post_swapped_destination = post_anchor_parent / "characterization"
    post_swapped = characterize_balance(replace(request, destination=post_swapped_destination))
    assert post_swapped.status == "MEASURED"
    assert post_sentinel.read_text() == "keep"
    assert list(post_redirected_parent.iterdir()) == [post_sentinel]
    assert not post_swapped_destination.exists()
    assert (post_anchor_backup / "characterization" / "characterization_report.json").is_file()
    post_anchor_parent.unlink()
    post_anchor_backup.rename(post_anchor_parent)
    monkeypatch.setattr(
        characterization_module,
        "_publish_staging_directory_at",
        original_publish,
    )

    original_stage = characterization_module.stage_balance_marsf_quartet
    for label, source_file, replacement in (
        ("metadata", metadata_file, b"{}\n"),
        ("config", config_file, b"{}\n"),
    ):
        original_bytes = source_file.read_bytes()

        def stage_then_mutate_snapshot(*args: object, **kwargs: object):
            staged_result = original_stage(*args, **kwargs)
            source_file.write_bytes(replacement)
            return staged_result

        monkeypatch.setattr(
            characterization_module,
            "stage_balance_marsf_quartet",
            stage_then_mutate_snapshot,
        )
        mutated_snapshot = characterize_balance(
            replace(
                path_request,
                destination=tmp_path / f"mutated-{label}-characterization",
            )
        )
        assert mutated_snapshot.status == "UNAVAILABLE"
        assert not (tmp_path / f"mutated-{label}-characterization").exists()
        source_file.write_bytes(original_bytes)
        monkeypatch.setattr(
            characterization_module,
            "stage_balance_marsf_quartet",
            original_stage,
        )

    expected = {
        "density": hashlib.sha256(paths["density"].read_bytes()).hexdigest(),
        "electron_temperature": hashlib.sha256(
            paths["electron_temperature"].read_bytes()
        ).hexdigest(),
        "ion_temperature": hashlib.sha256(paths["ion_temperature"].read_bytes()).hexdigest(),
        "toroidal_rotation": hashlib.sha256(paths["toroidal_rotation"].read_bytes()).hexdigest(),
        "equilibrium": hashlib.sha256(equilibrium.read_bytes()).hexdigest(),
        "equilibrium_parameters": hashlib.sha256(equilibrium_parameters.read_bytes()).hexdigest(),
        "oracle": hashlib.sha256(oracle_path.read_bytes()).hexdigest(),
    }
    path_expected = {
        **expected,
        "metadata": hashlib.sha256(metadata_file.read_bytes()).hexdigest(),
        "config": hashlib.sha256(config_file.read_bytes()).hexdigest(),
    }
    path_hashed = characterize_balance(
        replace(
            path_request,
            destination=tmp_path / "path-hashed-characterization",
            expected_sha256=path_expected,
        )
    )
    assert path_hashed.status == "MEASURED"
    assert path_hashed.report["actual_sha256"] == path_expected
    assert path_hashed.report["expected_hashes_requested"] is True

    cli_arguments = [
        "characterize-balance",
        str(paths["density"]),
        str(paths["electron_temperature"]),
        str(paths["ion_temperature"]),
        str(paths["toroidal_rotation"]),
        str(metadata_file),
        str(config_file),
        str(oracle_path),
        str(tmp_path / "cli-characterization"),
        "--equilibrium-file",
        str(equilibrium),
        "--equilibrium-parameters-file",
        str(equilibrium_parameters),
        "--equilibrium-provenance",
        "synthetic-equilibrium",
        "--domains",
        '{"core":[2,4]}',
        "--relative-floors",
        '{"n":1e12,"Te":1,"Ti":1,"Vz":1,"q":0.1}',
        "--interpolation-direction",
        "prepared_to_oracle",
        "--interpolation-method",
        "linear",
        "--format",
        "json",
    ]
    cli_measured = runner.invoke(app, cli_arguments)
    assert cli_measured.exit_code == 0, cli_measured.stdout
    cli_payload = json.loads(cli_measured.stdout)
    assert cli_payload["status"] == "MEASURED"
    assert cli_payload["overall_pass"] is None
    assert cli_payload["staged"]["directory"] == str(
        tmp_path / "cli-characterization" / "staged-marsf"
    )
    assert cli_payload["report_path"] == str(
        tmp_path / "cli-characterization" / "characterization_report.json"
    )
    cli_retained = Path(cli_payload["retained_staging"]["path"])
    assert cli_payload["retained_staging"]["state"] == "empty_container"
    assert cli_retained.is_dir()
    assert list(cli_retained.iterdir()) == []
    assert cli_retained.parent == tmp_path
    assert cli_payload["request_provenance"]["metadata"]["sha256"] == path_expected["metadata"]
    assert Path(cli_payload["prepared"]["report"]).is_file()

    floor_free_arguments = cli_arguments.copy()
    floor_free_floor_index = floor_free_arguments.index("--relative-floors")
    del floor_free_arguments[floor_free_floor_index : floor_free_floor_index + 2]
    floor_free_arguments[8] = str(tmp_path / "cli-floor-free-characterization")
    cli_floor_free = runner.invoke(app, floor_free_arguments)
    assert cli_floor_free.exit_code == 0, cli_floor_free.stdout
    floor_free_payload = json.loads(cli_floor_free.stdout)
    assert floor_free_payload["status"] == "MEASURED"
    assert floor_free_payload["comparisons"]["core"]["Vz"]["measurements"]["relative_rms"] is None

    cli_partial_arguments = [*cli_arguments[:-2]]
    cli_partial_arguments[8] = str(tmp_path / "cli-partial-characterization")
    cli_partial_arguments.extend(
        [
            "--tolerances",
            '{"core":{"n":{"absolute_max":1e100}}}',
            "--format",
            "json",
        ]
    )
    cli_partial = runner.invoke(app, cli_partial_arguments)
    assert cli_partial.exit_code == 4
    cli_partial_payload = json.loads(cli_partial.stdout)
    assert cli_partial_payload["status"] == "PARTIALLY_EVALUATED"
    assert cli_partial_payload["overall_pass"] is None
    assert cli_partial_payload["threshold_decisions"] == {"core/n/absolute_max": True}

    hashed = characterize_balance(
        replace(request, destination=tmp_path / "hashed-characterization", expected_sha256=expected)
    )
    assert hashed.status == "MEASURED"
    assert hashed.report["expected_sha256"] == expected
    assert hashed.report["actual_sha256"] == expected

    unknown_hash_key = characterize_balance(
        replace(
            request,
            destination=tmp_path / "unknown-hash-key-characterization",
            expected_sha256={**expected, "unexpected": "0" * 64},
        )
    )
    assert unknown_hash_key.status == "UNAVAILABLE"
    assert not (tmp_path / "unknown-hash-key-characterization").exists()

    mismatched = characterize_balance(
        replace(
            request,
            destination=tmp_path / "mismatch-characterization",
            expected_sha256={**expected, "oracle": "0" * 64},
        )
    )
    assert mismatched.status == "UNAVAILABLE"
    assert not (tmp_path / "mismatch-characterization").exists()

    existing = tmp_path / "existing-characterization"
    existing.mkdir()
    sentinel = existing / "sentinel.txt"
    sentinel.write_text("keep")
    refused = characterize_balance(replace(request, destination=existing))
    assert refused.status == "UNAVAILABLE"
    assert sentinel.read_text() == "keep"

    balance_alias = tmp_path / "balance-alias"
    balance_alias.symlink_to(balance, target_is_directory=True)
    aliased_destination = balance_alias / "nested-output"
    aliased = characterize_balance(replace(request, destination=aliased_destination))
    assert aliased.status == "UNAVAILABLE"
    assert not aliased_destination.exists()

    dotdot_destination = tmp_path / "not-created" / ".." / "dotdot-characterization"
    dotdot = characterize_balance(replace(request, destination=dotdot_destination))
    assert dotdot.status == "MEASURED"
    assert (tmp_path / "dotdot-characterization" / "characterization_report.json").is_file()

    for name, invalid_request in (
        (
            "missing-equilibrium-parameters-characterization",
            replace(
                request,
                destination=tmp_path / "missing-equilibrium-parameters-characterization",
                equilibrium_parameters_file=None,
            ),
        ),
        (
            "string-timeout-characterization",
            replace(
                request,
                destination=tmp_path / "string-timeout-characterization",
                equilibrium_timeout_seconds="2.0",
            ),
        ),
    ):
        invalid_result = characterize_balance(invalid_request)
        assert invalid_result.status == "UNAVAILABLE"
        assert not (tmp_path / name).exists()

    decisions = characterize_balance(
        replace(
            request,
            destination=tmp_path / "tolerated-characterization",
            tolerances={
                "absolute_rms": 0.0,
                "absolute_max": 0.0,
                "relative_rms": 0.0,
                "relative_max": 0.0,
            },
        )
    )
    assert decisions.status == "PASS"
    assert decisions.overall_pass is True
    assert decisions.report["threshold_decisions"]

    partial = characterize_balance(
        replace(
            request,
            destination=tmp_path / "partially-evaluated-characterization",
            tolerances={"core": {"n": {"absolute_max": 1.0e100}}},
        )
    )
    assert partial.status == "PARTIALLY_EVALUATED"
    assert partial.overall_pass is None
    assert partial.report["threshold_decisions"] == {"core/n/absolute_max": True}
    assert partial.report["comparisons"]["core"]["n"]["threshold_decisions"] == {
        "absolute_max": True
    }

    monkeypatch.undo()
    equilibrium_input = sources / "preprocessor" / "FourierModes.INP"
    equilibrium_input.parent.mkdir()
    equilibrium_input.write_text("synthetic-input\n")
    generator = sources / "preprocessor" / "equilibrium-preprocessor.py"
    run_counter = tmp_path / "equilibrium-run-count.txt"
    generator.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "Path('equil_r_q_psi.dat').write_text('# r q psi\\n1 2 0\\n3 3.5 .5\\n5 5 1\\n')\n"
        "Path('btor_rbig.dat').write_text('-16000 211\\n')\n"
        f"Path({str(run_counter)!r}).open('a').write('x')\n"
    )
    generator.chmod(0o755)
    generator_hashes = {
        "equilibrium_executable": hashlib.sha256(generator.read_bytes()).hexdigest(),
        "equilibrium_input:fouriermodes.inp": hashlib.sha256(
            equilibrium_input.read_bytes()
        ).hexdigest(),
    }
    original_equilibrium = sources / "original-equilibrium.gfile"
    original_equilibrium.write_bytes(equilibrium.read_bytes())
    generated_request = replace(
        request,
        destination=tmp_path / "generated-characterization",
        equilibrium_file=None,
        equilibrium_parameters_file=None,
        original_equilibrium=original_equilibrium,
        equilibrium_executable=generator,
        equilibrium_inputs=(equilibrium_input,),
        major_radius_cm=None,
        expected_sha256={
            **{key: value for key, value in expected.items() if key != "equilibrium_parameters"},
            **generator_hashes,
        },
    )
    generated = characterize_balance(generated_request)
    assert generated.status == "MEASURED"
    generated_expected = {
        **{key: value for key, value in expected.items() if key != "equilibrium_parameters"},
        **generator_hashes,
    }
    assert generated.report["actual_sha256"] == generated_expected
    assert generated.report["equilibrium_calculation"]["btor_gauss"] == -16000.0
    assert generated.report["equilibrium_calculation"]["r_big_cm"] == 211.0
    generated_equilibrium = generated.report["equilibrium_calculation"]
    assert (
        generated_equilibrium["equilibrium_sha256"]
        == hashlib.sha256(Path(generated_equilibrium["equilibrium_file"]).read_bytes()).hexdigest()
    )
    assert (
        generated_equilibrium["parameters_sha256"]
        == hashlib.sha256(Path(generated_equilibrium["parameters_file"]).read_bytes()).hexdigest()
    )
    assert generated.report["equilibrium_calculation"]["generator"] is not None
    assert run_counter.read_text() == "x"
    generated_request_json = json.loads(Path(generated.report["prepared"]["request"]).read_text())
    assert generated_request_json["setup"]["btor"] == -16000.0
    assert generated_request_json["setup"]["major_radius"] == 211.0
    prepared_report = json.loads(Path(generated.report["prepared"]["report"]).read_text())
    assert prepared_report["equilibrium"]["operation"] == (
        "generated by supplied equilibrium executable"
    )
    assert prepared_report["generator"] is not None
    for key in (
        "equilibrium",
        "equilibrium_input:fouriermodes.inp",
        "equilibrium_executable",
        "oracle",
    ):
        entry = generated.report["input_provenance"][key]
        assert Path(entry["snapshot_path"]).is_file()
        assert str(tmp_path / "generated-characterization") in entry["snapshot_path"]
    assert str(tmp_path / "generated-characterization") in generated.report["prepared"]["report"]

    # Every explicit source-file parent is protected, not just the BALANCE
    # profile directory.  These destinations must be rejected before staging
    # or preprocessor execution.
    for label, source_parent in (
        ("equilibrium", original_equilibrium.parent),
        ("oracle", oracle_path.parent),
        ("equilibrium-input", equilibrium_input.parent),
        ("equilibrium-executable", generator.parent),
    ):
        source_parent_destination = source_parent / f"{label}-nested-output"
        rejected = characterize_balance(
            replace(generated_request, destination=source_parent_destination)
        )
        assert rejected.status == "UNAVAILABLE"
        assert not source_parent_destination.exists()

    original_stage = characterization_module.stage_balance_marsf_quartet
    stage_calls: list[object] = []

    def unexpected_stage(*args, **kwargs):
        stage_calls.append(args)
        raise AssertionError("staging must not run after hash mismatch")

    monkeypatch.setattr(characterization_module, "stage_balance_marsf_quartet", unexpected_stage)
    bad_generator_hash = characterize_balance(
        replace(
            request,
            destination=tmp_path / "bad-generator-hash-characterization",
            equilibrium_file=None,
            equilibrium_parameters_file=None,
            original_equilibrium=original_equilibrium,
            equilibrium_executable=generator,
            equilibrium_inputs=(equilibrium_input,),
            expected_sha256={
                **{
                    key: value for key, value in expected.items() if key != "equilibrium_parameters"
                },
                **generator_hashes,
                "equilibrium_executable": "0" * 64,
            },
        )
    )
    assert bad_generator_hash.status == "UNAVAILABLE"
    assert stage_calls == []

    for label, source_file in (
        ("equilibrium", original_equilibrium),
        ("equilibrium-input", equilibrium_input),
        ("equilibrium-executable", generator),
        ("oracle", oracle_path),
    ):
        original_bytes = source_file.read_bytes()

        def stage_then_mutate(*args, _source_file=source_file, **kwargs):
            staged_result = original_stage(*args, **kwargs)
            _source_file.write_bytes(b"mutated-after-preflight\n")
            return staged_result

        monkeypatch.setattr(
            characterization_module,
            "stage_balance_marsf_quartet",
            stage_then_mutate,
        )
        mutated = characterize_balance(
            replace(
                generated_request,
                destination=tmp_path / f"mutated-{label}-characterization",
                expected_sha256=None,
            )
        )
        assert mutated.status == "UNAVAILABLE"
        assert not (tmp_path / f"mutated-{label}-characterization").exists()
        source_file.write_bytes(original_bytes)
        monkeypatch.setattr(
            characterization_module,
            "stage_balance_marsf_quartet",
            original_stage,
        )

    def interrupt_after_claiming_output(*args, **kwargs):
        destination = Path(args[1])
        destination.parent.mkdir(parents=True, exist_ok=True)
        raise KeyboardInterrupt

    interrupted_destination = tmp_path / "interrupted-characterization"
    monkeypatch.setattr(
        characterization_module,
        "stage_balance_marsf_quartet",
        interrupt_after_claiming_output,
    )
    with pytest.raises(KeyboardInterrupt):
        characterize_balance(replace(request, destination=interrupted_destination))
    assert not interrupted_destination.exists()


def test_characterize_balance_cli_emits_parseable_unavailable_json(tmp_path: Path) -> None:
    result = runner.invoke(
        app,
        [
            "characterize-balance",
            str(tmp_path / "density"),
            str(tmp_path / "te"),
            str(tmp_path / "ti"),
            str(tmp_path / "vz"),
            str(tmp_path / "metadata.json"),
            str(tmp_path / "config.json"),
            str(tmp_path / "oracle.h5"),
            str(tmp_path / "output"),
            "--equilibrium-file",
            str(tmp_path / "equilibrium.dat"),
            "--equilibrium-parameters-file",
            str(tmp_path / "btor_rbig.dat"),
            "--equilibrium-provenance",
            "synthetic-equilibrium",
            "--major-radius-cm",
            "200",
            "--domains",
            '{"all":[1,2]}',
            "--relative-floors",
            '{"n":1,"Te":1,"Ti":1,"Vz":1,"q":1}',
            "--interpolation-direction",
            "prepared_to_oracle",
            "--interpolation-method",
            "linear",
            "--format",
            "json",
        ],
    )

    assert result.exit_code != 0
    payload = __import__("json").loads(result.stdout)
    assert payload["status"] == "UNAVAILABLE"
    assert payload["retained_staging"]["path"] is None
    assert payload["retained_staging"]["state"] == "not_created"
