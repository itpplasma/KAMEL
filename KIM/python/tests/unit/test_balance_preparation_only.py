"""BALANCE preparation works without an external reference dataset."""

import hashlib
import json
import subprocess
from pathlib import Path

import pytest
from kim import BuiltinPlasma, PlasmaIsotope, SimulationConfig
from kim.balance_characterization import BalanceCharacterizationRequest, characterize_balance
from kim.cli import app
from typer.testing import CliRunner


@pytest.fixture
def preparation_inputs(tmp_path, monkeypatch):
    sources = tmp_path / "sources"
    sources.mkdir()
    profiles = {}
    for role, value in {
        "density": 1e19,
        "electron_temperature": 100,
        "ion_temperature": 80,
        "toroidal_rotation": 10,
    }.items():
        profiles[role] = sources / f"{role}.dat"
        profiles[role].write_text(f"0 {value}\n0.5 {value}\n1 {value}\n")
    equilibrium = sources / "equil_r_q_psi.dat"
    equilibrium.write_text("# r q psi\n1 2 0\n3 3 .5\n5 4 1\n")
    parameters = sources / "btor_rbig.dat"
    parameters.write_text("-18000 165\n")
    metadata = {
        "source": "synthetic",
        "coordinate": "rho_pol",
        "coordinate_unit": "1",
        "density_unit": "1/m^3",
        "electron_temperature_unit": "eV",
        "ion_temperature_unit": "eV",
        "toroidal_rotation_unit": "rad/s",
    }
    config = SimulationConfig.electrostatic_periodic(
        profiles=Path("unused"),
        plasma=BuiltinPlasma(isotope=PlasmaIsotope.DEUTERIUM),
        btor=-18000,
        major_radius=165,
        m_mode=7,
        n_mode=2,
        frequency=0,
        br_boundary_real=1,
        br_boundary_imag=0,
        radial_minimum=1,
        plasma_radius=5,
    )

    def unexpected(*args, **kwargs):
        pytest.fail("preparation must not read a reference, compare profiles, or launch a solver")

    monkeypatch.setattr("kim.balance_characterization.read_ql_balance_oracle", unexpected)
    monkeypatch.setattr("kim.balance_characterization._compare_all", unexpected)
    run = subprocess.run

    def git_only(command, *args, **kwargs):
        if command[0] not in ("git", "./.verified-equilibrium-generator"):
            unexpected()
        return run(command, *args, **kwargs)

    monkeypatch.setattr("kim.preparation.subprocess.run", git_only)
    return dict(
        **profiles,
        metadata=metadata,
        config=config,
        equilibrium_provenance="synthetic-equilibrium",
        equilibrium_file=equilibrium,
        equilibrium_parameters_file=parameters,
        destination=tmp_path / "prepared-output",
    )


@pytest.mark.parametrize("verify_hashes", [False, True])
@pytest.mark.parametrize("route", ["precomputed", "generated"])
def test_api_prepares_without_reference_or_comparison_settings(
    preparation_inputs, verify_hashes, route
):
    inputs = preparation_inputs
    if route == "generated":
        sources = inputs["equilibrium_file"].parent
        generator = sources / "generator"
        generator.write_text(
            "#!/bin/sh\n"
            "printf '# r q psi\\n1 2 0\\n3 3 .5\\n5 4 1\\n' > equil_r_q_psi.dat\n"
            "printf '%s\\n' '-18000 165' > btor_rbig.dat\n"
        )
        generator.chmod(0o755)
        original = sources / "original-equilibrium"
        original.write_text("synthetic preprocessor input\n")
        inputs.pop("equilibrium_file")
        inputs.pop("equilibrium_parameters_file")
        inputs.update(original_equilibrium=original, equilibrium_executable=generator)
    if verify_hashes:
        identities = {
            key: inputs[key]
            for key in ("density", "electron_temperature", "ion_temperature", "toroidal_rotation")
        }
        if route == "precomputed":
            identities.update(
                equilibrium=inputs["equilibrium_file"],
                equilibrium_parameters=inputs["equilibrium_parameters_file"],
            )
        else:
            identities.update(
                equilibrium=inputs["original_equilibrium"],
                equilibrium_executable=inputs["equilibrium_executable"],
            )
        inputs["expected_sha256"] = {
            key: hashlib.sha256(path.read_bytes()).hexdigest() for key, path in identities.items()
        }
    request = BalanceCharacterizationRequest(**inputs)
    result = characterize_balance(request)
    assert result.status == "PREPARED", result.report
    assert result.overall_pass is None
    assert result.report["threshold_decisions"] is None
    assert result.report["comparison_performed"] is False
    assert result.report["comparison_configuration"] is None
    assert result.report["comparisons"] == {}
    assert "oracle" not in result.report["actual_sha256"]
    assert "oracle" not in result.report["input_provenance"]
    prepared = Path(result.report["prepared"]["directory"])
    assert (prepared / "request.json").is_file()
    for name in ("n.dat", "Te.dat", "Ti.dat", "Vz.dat", "q.dat"):
        assert (prepared / "profiles" / name).is_file()
    assert json.loads(result.report_path.read_text())["status"] == "PREPARED"

    duplicate = characterize_balance(request)
    assert duplicate.status == "UNAVAILABLE"
    assert (prepared / "request.json").is_file()


def test_cli_prepares_without_reference_or_comparison_options(preparation_inputs):
    inputs = preparation_inputs
    sources = inputs["equilibrium_file"].parent
    metadata = sources / "metadata.json"
    metadata.write_text(json.dumps(inputs["metadata"]))
    config = sources / "request.json"
    config.write_text(inputs["config"].model_dump_json())
    arguments = [
        "characterize-balance",
        *(
            str(inputs[key])
            for key in ("density", "electron_temperature", "ion_temperature", "toroidal_rotation")
        ),
        str(metadata),
        str(config),
        str(inputs["destination"]),
        "--equilibrium-file",
        str(inputs["equilibrium_file"]),
        "--equilibrium-parameters-file",
        str(inputs["equilibrium_parameters_file"]),
        "--equilibrium-provenance",
        inputs["equilibrium_provenance"],
    ]
    result = CliRunner().invoke(app, arguments)
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["status"] == "PREPARED"


@pytest.mark.parametrize(
    "option",
    [
        {"tolerances": {"n": {"absolute_max": 1}}},
        {"domains": {"all": (1, 5)}},
        {"relative_floors": {"n": 1}},
        {"interpolation_direction": "prepared_to_oracle"},
        {"interpolation_method": "linear"},
    ],
)
def test_comparison_options_require_a_reference(preparation_inputs, option):
    result = characterize_balance(**preparation_inputs, **option)
    assert result.status == "UNAVAILABLE"
    assert "require a reference" in result.report["error"]["message"]
    assert not preparation_inputs["destination"].exists()


def test_expected_reference_hash_requires_a_reference(preparation_inputs):
    result = characterize_balance(**preparation_inputs, expected_sha256={"oracle": "0" * 64})
    assert result.status == "UNAVAILABLE"
    assert "unavailable oracle" in result.report["error"]["message"]
    assert not preparation_inputs["destination"].exists()


def test_explicit_missing_reference_is_not_silently_skipped(preparation_inputs):
    result = characterize_balance(
        **preparation_inputs,
        oracle=preparation_inputs["destination"].parent / "missing-reference.h5",
        domains={"all": (1, 5)},
        interpolation_direction="prepared_to_oracle",
        interpolation_method="linear",
    )
    assert result.status == "UNAVAILABLE"
    assert not preparation_inputs["destination"].exists()
