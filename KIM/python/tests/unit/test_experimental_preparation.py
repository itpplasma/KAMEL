from __future__ import annotations

import json
import stat
from pathlib import Path

import numpy as np
import pytest
from kim import (
    BuiltinPlasma,
    ExperimentalInputError,
    MarsFMetadata,
    PlasmaIsotope,
    ProfileConfig,
    SimulationConfig,
    prepare_marsf_case,
    read_marsf_profiles,
)


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


def write_marsf_case(directory: Path) -> None:
    directory.mkdir()
    profiles = {
        "PROFDEN.IN": ([0.0, 0.5, 1.0], [1.0e19, 2.0e19, 3.0e19]),
        "PROFTE.IN": ([0.0, 0.5, 1.0], [100.0, 80.0, 40.0]),
        "PROFTI.IN": ([0.0, 0.5, 1.0], [90.0, 70.0, 30.0]),
        "PROFROT.IN": ([0.0, 0.5, 1.0], [0.0, -2.5e3, -5.0e3]),
    }
    for filename, (coordinates, values) in profiles.items():
        rows = ["MARS-F profile: original source values"]
        rows.extend(
            f"{coordinate:.8e} {value:.8e}" for coordinate, value in zip(coordinates, values)
        )
        (directory / filename).write_text("\n".join(rows) + "\n", encoding="utf-8")


def write_equilibrium(path: Path) -> None:
    path.write_text(
        "# radius q psi\n" "0.0 1.0 0.0 0 0\n" "10.0 1.5 0.25 0 0\n" "20.0 2.0 1.0 0 0\n",
        encoding="utf-8",
    )


def config() -> SimulationConfig:
    return SimulationConfig.electrostatic_periodic(
        profiles=Path("unused"),
        plasma=BuiltinPlasma(isotope=PlasmaIsotope.DEUTERIUM),
        btor=-17_977.413,
        major_radius=165.0,
        m_mode=7,
        n_mode=2,
        frequency=0.0,
        br_boundary_real=1.0,
        br_boundary_imag=0.0,
        radial_minimum=3.0,
        plasma_radius=19.0,
    )


def test_prepares_marsf_profiles_with_explicit_mapping_and_report(tmp_path: Path) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    write_marsf_case(source_directory)
    write_equilibrium(equilibrium)
    source = read_marsf_profiles(source_directory, metadata())

    prepared = prepare_marsf_case(
        source, config(), tmp_path / "prepared", equilibrium_file=equilibrium
    )

    assert prepared.config.profiles.directory == prepared.profiles
    assert prepared.equilibrium == prepared.directory / "equilibrium" / equilibrium.name
    np.testing.assert_allclose(
        np.loadtxt(prepared.profiles / "n.dat"),
        [[0.0, 1.0e13], [10.0, 2.0e13], [20.0, 3.0e13]],
    )
    np.testing.assert_allclose(
        np.loadtxt(prepared.profiles / "Vz.dat")[:, 1], [0.0, -2.5e5, -5.0e5]
    )
    np.testing.assert_allclose(np.loadtxt(prepared.profiles / "q.dat")[:, 1], [1.0, 1.5, 2.0])

    request = json.loads(prepared.request.read_text(encoding="utf-8"))
    assert request["profiles"]["directory"] == "./profiles"
    report = json.loads(prepared.report.read_text(encoding="utf-8"))
    assert report["target"] == "KIM-CGS-r_eff"
    assert report["equilibrium"]["source"] == str(equilibrium.absolute())
    assert set(report["source_hashes"]) == {
        "source/PROFDEN.IN",
        "source/PROFTE.IN",
        "source/PROFTI.IN",
        "source/PROFROT.IN",
        "equilibrium/equil_r_q_psi.dat",
    }


def test_rejects_profile_extrapolation(tmp_path: Path) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    write_marsf_case(source_directory)
    equilibrium.write_text(
        "# radius q psi\n10.0 1.0 0.25\n20.0 2.0 1.0\n",
        encoding="utf-8",
    )
    (source_directory / "PROFDEN.IN").write_text(
        "header\n0.0 1.0e19\n0.5 2.0e19\n0.8 3.0e19\n",
        encoding="utf-8",
    )
    source = read_marsf_profiles(source_directory, metadata())

    with pytest.raises(ExperimentalInputError, match="outside equilibrium"):
        prepare_marsf_case(source, config(), tmp_path / "prepared", equilibrium_file=equilibrium)

    assert not (tmp_path / "prepared").exists()


def test_generates_equilibrium_with_supplied_executable(tmp_path: Path) -> None:
    source_directory = tmp_path / "marsf"
    write_marsf_case(source_directory)
    generator = tmp_path / "generator.py"
    generator.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "Path('equil_r_q_psi.dat').write_text("
        "'# radius q psi\\n0 1 0\\n10 1.5 .25\\n20 2 1\\n')\n",
        encoding="utf-8",
    )
    generator.chmod(generator.stat().st_mode | stat.S_IXUSR)
    source = read_marsf_profiles(source_directory, metadata())

    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.chdir(tmp_path)
    try:
        prepared = prepare_marsf_case(
            source,
            config(),
            tmp_path / "prepared",
            equilibrium_executable="generator.py",
        )
    finally:
        monkeypatch.undo()

    assert prepared.equilibrium.is_file()
    assert "generated by supplied equilibrium executable" in prepared.report.read_text(
        encoding="utf-8"
    )
    report = json.loads(prepared.report.read_text(encoding="utf-8"))
    assert report["generator"]["executable"] == str(generator.absolute())
    assert len(report["generator"]["sha256"]) == 64


def test_refuses_to_overwrite_prepared_case(tmp_path: Path) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    write_marsf_case(source_directory)
    write_equilibrium(equilibrium)
    destination = tmp_path / "prepared"
    destination.mkdir()
    source = read_marsf_profiles(source_directory, metadata())

    with pytest.raises(ExperimentalInputError, match="already exists"):
        prepare_marsf_case(source, config(), destination, equilibrium_file=equilibrium)


def test_rejects_nonexistent_destination_symlink(tmp_path: Path) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    write_marsf_case(source_directory)
    write_equilibrium(equilibrium)
    destination = tmp_path / "prepared"
    destination.symlink_to(tmp_path / "missing-target")
    source = read_marsf_profiles(source_directory, metadata())

    with pytest.raises(ExperimentalInputError, match="already exists"):
        prepare_marsf_case(source, config(), destination, equilibrium_file=equilibrium)


def test_matches_fortran_spline_interpolation_between_grid_points(tmp_path: Path) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    write_marsf_case(source_directory)
    equilibrium.write_text(
        "# radius q psi\n0.0 1.0 0.0\n10.0 1.5 0.36\n20.0 2.0 1.0\n",
        encoding="utf-8",
    )
    source = read_marsf_profiles(source_directory, metadata())

    prepared = prepare_marsf_case(
        source, config(), tmp_path / "prepared", equilibrium_file=equilibrium
    )

    assert np.loadtxt(prepared.profiles / "n.dat")[1, 1] == pytest.approx(2.3206328889e13)


def test_accepts_negative_polarity_equilibrium(tmp_path: Path) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    write_marsf_case(source_directory)
    equilibrium.write_text(
        "# radius q psi\n0.0 1.0 0.0\n10.0 1.5 -0.25\n20.0 2.0 -1.0\n",
        encoding="utf-8",
    )
    source = read_marsf_profiles(source_directory, metadata())

    prepared = prepare_marsf_case(
        source, config(), tmp_path / "prepared", equilibrium_file=equilibrium
    )

    np.testing.assert_allclose(np.loadtxt(prepared.profiles / "q.dat")[:, 1], [1.0, 1.5, 2.0])


def test_uses_profile_filenames_from_the_request(tmp_path: Path) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    write_marsf_case(source_directory)
    write_equilibrium(equilibrium)
    source = read_marsf_profiles(source_directory, metadata())
    base_config = config()
    custom_profiles = ProfileConfig(
        directory=Path("unused"),
        density_file="density.profile",
        electron_temperature_file="electron-temperature.profile",
        ion_temperature_file="ion-temperature.profile",
        toroidal_velocity_file="rotation.profile",
        safety_factor_file="safety-factor.profile",
    )
    custom_config = base_config.model_copy(update={"profiles": custom_profiles})

    prepared = prepare_marsf_case(
        source, custom_config, tmp_path / "prepared", equilibrium_file=equilibrium
    )

    assert (prepared.profiles / "density.profile").is_file()
    assert (prepared.profiles / "electron-temperature.profile").is_file()
    assert (prepared.profiles / "ion-temperature.profile").is_file()
    assert (prepared.profiles / "rotation.profile").is_file()
    assert (prepared.profiles / "safety-factor.profile").is_file()
    assert not (prepared.profiles / "n.dat").exists()


def test_rejects_colliding_profile_filenames(tmp_path: Path) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    write_marsf_case(source_directory)
    write_equilibrium(equilibrium)
    source = read_marsf_profiles(source_directory, metadata())
    base_config = config()
    custom_profiles = ProfileConfig(
        directory=Path("unused"),
        density_file="same.dat",
        safety_factor_file="same.dat",
    )
    custom_config = base_config.model_copy(update={"profiles": custom_profiles})

    with pytest.raises(ExperimentalInputError, match="profile filenames must be unique"):
        prepare_marsf_case(
            source, custom_config, tmp_path / "prepared", equilibrium_file=equilibrium
        )


def test_rejects_sqrt_psiN_coordinates_outside_the_declared_domain(tmp_path: Path) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    write_marsf_case(source_directory)
    write_equilibrium(equilibrium)
    (source_directory / "PROFDEN.IN").write_text(
        "header\n-0.1 1.0e19\n0.5 2.0e19\n1.0 3.0e19\n",
        encoding="utf-8",
    )
    source = read_marsf_profiles(source_directory, metadata())

    with pytest.raises(ExperimentalInputError, match="sqrt_psiN coordinate"):
        prepare_marsf_case(source, config(), tmp_path / "prepared", equilibrium_file=equilibrium)
