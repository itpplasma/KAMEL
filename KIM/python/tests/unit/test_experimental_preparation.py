from __future__ import annotations

import hashlib
import json
import os
import stat
import sys
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
from kim.preparation import (
    read_equilibrium_parameters,
    run_equilibrium_calculation,
    simulation_config_for_equilibrium,
)
from kim.profiles import ProfileSet


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


def test_reads_btor_and_r_big_from_equilibrium_output(tmp_path: Path) -> None:
    parameters = tmp_path / "btor_rbig.dat"
    parameters.write_text("-1.757319212D+04 1.699117661D+02\n", encoding="utf-8")

    result = read_equilibrium_parameters(parameters)

    assert result.btor_gauss == pytest.approx(-17573.19212)
    assert result.r_big_cm == pytest.approx(169.9117661)
    assert result.sha256 == hashlib.sha256(parameters.read_bytes()).hexdigest()


def test_config_uses_equilibrium_calculated_btor_and_r_big(tmp_path: Path) -> None:
    parameters_path = tmp_path / "btor_rbig.dat"
    parameters_path.write_text("-17573.19212 169.9117661\n", encoding="utf-8")
    parameters = read_equilibrium_parameters(parameters_path)

    result = simulation_config_for_equilibrium(config(), parameters)

    assert result.setup.btor == pytest.approx(-17573.19212)
    assert result.setup.major_radius == pytest.approx(169.9117661)
    assert config().setup.btor == -17_977.413
    assert config().setup.major_radius == 165.0


def test_runs_equilibrium_once_and_requires_paired_outputs(tmp_path: Path) -> None:
    input_file = tmp_path / "equilibrium.inp"
    input_file.write_text("source input\n", encoding="utf-8")
    generator = tmp_path / "equilibrium-generator.py"
    generator.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "Path('equil_r_q_psi.dat').write_text('# r q psi\\n1 2 0\\n3 3.5 .5\\n5 5 1\\n')\n"
        "Path('btor_rbig.dat').write_text('-17573.19212 169.9117661\\n')\n"
    )
    generator.chmod(0o755)

    result = run_equilibrium_calculation(
        generator,
        (input_file,),
        tmp_path / "calculation",
    )

    assert result.equilibrium_file.is_file()
    assert result.parameters_file.is_file()
    assert result.parameters.btor_gauss == pytest.approx(-17573.19212)
    assert result.parameters.r_big_cm == pytest.approx(169.9117661)
    assert (
        result.equilibrium_sha256
        == hashlib.sha256(result.equilibrium_file.read_bytes()).hexdigest()
    )
    assert result.generator_provenance["input_sha256"] == {
        input_file.name: hashlib.sha256(input_file.read_bytes()).hexdigest()
    }


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="requires Linux descriptor paths")
def test_runs_equilibrium_script_in_descriptor_anchored_directory(tmp_path: Path) -> None:
    generator = tmp_path / "generator.py"
    generator.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "Path('equil_r_q_psi.dat').write_text('# r q psi\\n1 2 0\\n3 3.5 .5\\n5 5 1\\n')\n"
        "Path('btor_rbig.dat').write_text('-18000 165\\n')\n"
    )
    generator.chmod(0o755)
    descriptor = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        result = run_equilibrium_calculation(
            generator, (), Path(f"/proc/self/fd/{descriptor}") / "calculation"
        )
        assert result.equilibrium_file.is_file()
        assert result.parameters.r_big_cm == 165.0
        assert (
            result.generator_provenance["executed_sha256"]
            == hashlib.sha256(generator.read_bytes()).hexdigest()
        )
    finally:
        os.close(descriptor)


@pytest.mark.parametrize(
    "contents",
    ["-18000 0\n", "-18000 -165\n", "-18000 nan\n", "-18000 165 extra\n"],
)
def test_rejects_invalid_equilibrium_parameter_output(tmp_path: Path, contents: str) -> None:
    parameters = tmp_path / "btor_rbig.dat"
    parameters.write_text(contents, encoding="utf-8")

    with pytest.raises(ExperimentalInputError, match="btor_rbig.dat"):
        read_equilibrium_parameters(parameters)


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


def q_operation_record(report_path: Path) -> dict[str, object]:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    q_operations = [operation for operation in report["operations"] if operation["quantity"] == "q"]
    assert len(q_operations) == 1
    return q_operations[0]


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
    np.testing.assert_array_equal(np.loadtxt(prepared.profiles / "q.dat")[:, 1], [1.0, 1.5, 2.0])

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


def test_preparation_uses_equilibrium_btor_and_r_big_outputs(tmp_path: Path) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    parameters = tmp_path / "btor_rbig.dat"
    write_marsf_case(source_directory)
    write_equilibrium(equilibrium)
    parameters.write_text("-17573.19212 169.9117661\n", encoding="utf-8")
    source = read_marsf_profiles(source_directory, metadata())

    prepared = prepare_marsf_case(
        source,
        config(),
        tmp_path / "prepared-from-equilibrium",
        equilibrium_file=equilibrium,
        equilibrium_parameters_file=parameters,
    )

    assert prepared.config.setup.btor == pytest.approx(-17573.19212)
    assert prepared.config.setup.major_radius == pytest.approx(169.9117661)
    staged_parameters = prepared.directory / "equilibrium" / "btor_rbig.dat"
    assert staged_parameters.read_bytes() == parameters.read_bytes()
    report = json.loads(prepared.report.read_text(encoding="utf-8"))
    assert report["equilibrium_parameters"] == {
        "staged_file": "equilibrium/btor_rbig.dat",
        "source_hash": hashlib.sha256(parameters.read_bytes()).hexdigest(),
        "btor_gauss": -17573.19212,
        "r_big_cm": 169.9117661,
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

    np.testing.assert_array_equal(np.loadtxt(prepared.profiles / "q.dat")[:, 1], [1.0, 1.5, 2.0])


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


def test_accepts_sqrt_psiN_profiles_extended_beyond_lcfs(tmp_path: Path) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    write_marsf_case(source_directory)
    write_equilibrium(equilibrium)
    for filename in ("PROFDEN.IN", "PROFTE.IN", "PROFTI.IN", "PROFROT.IN"):
        (source_directory / filename).write_text(
            "header\n0.0 1.0\n0.5 2.0\n1.02 3.0\n",
            encoding="utf-8",
        )
    source = read_marsf_profiles(source_directory, metadata())

    prepared = prepare_marsf_case(
        source, config(), tmp_path / "prepared", equilibrium_file=equilibrium
    )

    np.testing.assert_array_equal(np.loadtxt(prepared.profiles / "n.dat")[:, 0], [0.0, 10.0, 20.0])
    report = json.loads(prepared.report.read_text(encoding="utf-8"))
    assert report["coordinate_operation"] == (
        "natural cubic interpolation from sqrt_psiN to equilibrium r_eff"
    )
    assert report["coordinate_mapping"] == {
        "source_coordinate": "sqrt_psiN",
        "target_coordinate": "r_eff",
        "method": "natural cubic interpolation",
    }


@pytest.mark.parametrize("unit", ["cm", "m"])
def test_r_eff_report_preserves_version_one_coordinate_operation(tmp_path: Path, unit: str) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    write_marsf_case(source_directory)
    write_equilibrium(equilibrium)
    radii = [0.0, 10.0, 20.0] if unit == "cm" else [0.0, 0.1, 0.2]
    for path in source_directory.iterdir():
        rows = np.loadtxt(path, skiprows=1)
        rows[:, 0] = radii
        np.savetxt(path, rows, header="MARS-F profile")
    source = read_marsf_profiles(
        source_directory, metadata(coordinate="r_eff", coordinate_unit=unit)
    )

    prepared = prepare_marsf_case(
        source, config(), tmp_path / "prepared", equilibrium_file=equilibrium
    )

    report = json.loads(prepared.report.read_text(encoding="utf-8"))
    expected = (
        "preserved explicit r_eff grid"
        if unit == "cm"
        else "converted explicit r_eff grid from m to cm"
    )
    assert report["schema_version"] == 1
    assert report["coordinate_operation"] == expected
    assert report["coordinate_mapping"]["method"] == expected


def test_rejects_nonfinite_squared_sqrt_psiN_coordinates(tmp_path: Path) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    write_marsf_case(source_directory)
    write_equilibrium(equilibrium)
    (source_directory / "PROFDEN.IN").write_text(
        "header\n0.0 1.0e19\n0.5 2.0e19\n1.0e200 3.0e19\n",
        encoding="utf-8",
    )
    source = read_marsf_profiles(source_directory, metadata())

    with pytest.raises(ExperimentalInputError, match="squared sqrt_psiN coordinate is not finite"):
        prepare_marsf_case(source, config(), tmp_path / "prepared", equilibrium_file=equilibrium)

    assert not (tmp_path / "prepared").exists()


def test_rejects_negative_sqrt_psiN_coordinates(tmp_path: Path) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    write_marsf_case(source_directory)
    write_equilibrium(equilibrium)
    (source_directory / "PROFDEN.IN").write_text(
        "header\n-0.1 1.0e19\n0.5 2.0e19\n1.0 3.0e19\n",
        encoding="utf-8",
    )
    source = read_marsf_profiles(source_directory, metadata())

    with pytest.raises(ExperimentalInputError, match="sqrt_psiN coordinate must be nonnegative"):
        prepare_marsf_case(source, config(), tmp_path / "prepared", equilibrium_file=equilibrium)


def test_default_q_operation_preserves_supplied_equilibrium_q(tmp_path: Path) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    write_marsf_case(source_directory)
    write_equilibrium(equilibrium)
    source = read_marsf_profiles(source_directory, metadata())

    prepared = prepare_marsf_case(
        source,
        config(),
        tmp_path / "prepared",
        equilibrium_file=equilibrium,
    )

    np.testing.assert_array_equal(np.loadtxt(prepared.profiles / "q.dat")[:, 1], [1.0, 1.5, 2.0])
    q_record = q_operation_record(prepared.report)
    assert q_record["factor"] == 1.0
    assert q_record["operation"] == "preserve"
    assert q_record["quantity"] == "q"
    assert q_record["source_unit"] == "1"
    assert q_record["target_unit"] == "1"


@pytest.mark.parametrize(
    ("q_operation", "expected_q"),
    [
        ("preserve", [1.0, 1.5, 2.0]),
        ("negate", [-1.0, -1.5, -2.0]),
    ],
)
def test_explicit_q_operation_controls_written_equilibrium_q(
    tmp_path: Path, q_operation: str, expected_q: list[float]
) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    write_marsf_case(source_directory)
    write_equilibrium(equilibrium)
    source = read_marsf_profiles(source_directory, metadata())

    prepared = prepare_marsf_case(
        source,
        config(),
        tmp_path / "prepared",
        equilibrium_file=equilibrium,
        q_operation=q_operation,
    )

    np.testing.assert_array_equal(np.loadtxt(prepared.profiles / "q.dat")[:, 1], expected_q)


@pytest.mark.parametrize(
    ("q_operation", "factor"),
    [("preserve", 1.0), ("negate", -1.0)],
)
def test_q_operation_is_recorded_in_conversion_report(
    tmp_path: Path, q_operation: str, factor: float
) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    write_marsf_case(source_directory)
    write_equilibrium(equilibrium)
    source = read_marsf_profiles(source_directory, metadata())

    prepared = prepare_marsf_case(
        source,
        config(),
        tmp_path / "prepared",
        equilibrium_file=equilibrium,
        q_operation=q_operation,
    )

    q_record = q_operation_record(prepared.report)
    assert q_record["factor"] == factor
    assert q_record["operation"] == q_operation
    assert q_record["quantity"] == "q"
    assert q_record["source_unit"] == "1"
    assert q_record["target_unit"] == "1"


def test_invalid_q_operation_is_rejected_before_preparation(tmp_path: Path) -> None:
    source_directory = tmp_path / "marsf"
    write_marsf_case(source_directory)
    generator = tmp_path / "generator.py"
    counter = tmp_path / "generator-called"
    generator.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        f"Path({str(counter)!r}).write_text('called\\n')\n"
        "Path('equil_r_q_psi.dat').write_text('# radius q psi\\n0 1 0\\n10 1.5 .25\\n20 2 1\\n')\n",
        encoding="utf-8",
    )
    generator.chmod(generator.stat().st_mode | stat.S_IXUSR)
    equilibrium_input = tmp_path / "equilibrium-input.dat"
    equilibrium_input.write_bytes(b"equilibrium input remains untouched\n")
    source_bytes = {path.name: path.read_bytes() for path in source_directory.iterdir()}
    generator_bytes = generator.read_bytes()
    equilibrium_input_bytes = equilibrium_input.read_bytes()
    source = read_marsf_profiles(source_directory, metadata())

    with pytest.raises(ExperimentalInputError, match="q operation"):
        prepare_marsf_case(
            source,
            config(),
            tmp_path / "prepared",
            equilibrium_executable=generator,
            equilibrium_input_files=(equilibrium_input,),
            q_operation="infer",
        )

    assert {path.name: path.read_bytes() for path in source_directory.iterdir()} == source_bytes
    assert generator.read_bytes() == generator_bytes
    assert equilibrium_input.read_bytes() == equilibrium_input_bytes
    assert not counter.exists()
    assert not (tmp_path / "prepared").exists()


def test_negated_analytic_q_profile_has_signed_kim_resonance(tmp_path: Path) -> None:
    source_directory = tmp_path / "marsf"
    equilibrium = tmp_path / "equil_r_q_psi.dat"
    write_marsf_case(source_directory)
    equilibrium.write_text(
        "# radius q psi\n" "0.000 2.000 0.000\n" "58.496 3.500 0.500\n" "100.000 4.500 1.000\n",
        encoding="utf-8",
    )
    source = read_marsf_profiles(source_directory, metadata())
    base_config = config()
    resonant_config = base_config.model_copy(
        update={
            "grid": base_config.grid.model_copy(update={"plasma_radius": 80.0}),
        }
    )

    prepared = prepare_marsf_case(
        source,
        resonant_config,
        tmp_path / "prepared",
        equilibrium_file=equilibrium,
        q_operation="negate",
    )

    np.testing.assert_array_equal(np.loadtxt(prepared.profiles / "q.dat")[:, 1], [-2.0, -3.5, -4.5])
    validation = ProfileSet.from_simulation(prepared.config).validate_for(prepared.config)
    assert validation.resonance_radius == pytest.approx(58.496)
