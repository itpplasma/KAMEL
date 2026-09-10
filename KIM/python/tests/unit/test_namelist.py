from __future__ import annotations

from pathlib import Path

import f90nml
import pytest
from kim.config import (
    BuiltinPlasma,
    ExplicitPlasma,
    GridConfig,
    IonSpecies,
    PeriodicConfig,
    PhysicsConfig,
    PlasmaIsotope,
    SimulationConfig,
)
from kim.errors import ConfigurationError
from kim.namelist import dumps_namelist, load_namelist

FIXTURES = Path(__file__).parents[1] / "fixtures" / "namelists"


def physical_case() -> dict[str, object]:
    return {
        "profiles": Path("/source/profiles"),
        "plasma": BuiltinPlasma(isotope=PlasmaIsotope.DEUTERIUM),
        "btor": -17_977.413,
        "major_radius": 165.0,
        "m_mode": 7,
        "n_mode": 2,
        "frequency": 0.0,
        "br_boundary_real": 1.0,
        "br_boundary_imag": 0.0,
        "radial_minimum": 3.0,
        "plasma_radius": 63.0,
    }


def configuration(run_type: str) -> SimulationConfig:
    factory = getattr(SimulationConfig, run_type)
    options: dict[str, object] = {}
    if run_type == "electrostatic_periodic":
        options = {
            "grid": GridConfig(
                radial_minimum=3.0,
                plasma_radius=63.0,
                resonance_width=2.5,
                resonance_amplification=70.0,
                maximum_step_scale=2.5,
                l_space_dim=225,
                rg_space_dim=750,
                rg_spacing="non-equidistant",
                l_spacing="non-equidistant",
                larmor_skip_factor=100.0,
                gauss_nodes_x=11,
                gauss_nodes_x_prime=10,
                gauss_nodes_theta=7,
            ),
            "periodic": PeriodicConfig(
                as_is_width_scale=2.25,
                transition_width_scale=13.53,
                wavenumber_cutoff_scale=44.6,
                n_rg=2048,
            ),
        }
    return factory(**physical_case(), **options)


@pytest.mark.parametrize(
    ("run_type", "factory"),
    [
        ("electrostatic_periodic", SimulationConfig.electrostatic_periodic),
        ("electrostatic", SimulationConfig.electrostatic),
        ("flr2", SimulationConfig.flr2),
    ],
)
def test_canonical_namelist_matches_golden_fixture(run_type: str, factory: object) -> None:
    rendered = dumps_namelist(configuration(run_type))

    assert rendered == (FIXTURES / f"{run_type}.nml").read_text()


def test_groups_follow_fortran_read_order_and_fields_are_translated() -> None:
    rendered = dumps_namelist(configuration("electrostatic_periodic"))
    groups = [line[1:].strip() for line in rendered.splitlines() if line.startswith("&")]
    parsed = f90nml.reads(rendered)

    assert groups == [
        "KIM_CONFIG",
        "WKB_DISPERSION",
        "KIM_IO",
        "KIM_SETUP",
        "KIM_GRID",
        "KIM_PROFILES",
        "KIM_PERIODIC",
    ]
    assert parsed["kim_setup"]["r0"] == 165.0
    assert parsed["kim_grid"]["r_min"] == 3.0
    assert parsed["kim_grid"]["r_plas"] == 63.0
    assert parsed["kim_config"]["type_of_run"] == "electrostatic_periodic"
    assert parsed["kim_io"]["hdf5_input"] is False
    assert parsed["kim_io"]["hdf5_output"] is True
    assert parsed["kim_io"]["profile_location"] == "/source/profiles/"
    assert parsed["kim_profiles"]["input_profile_dir"] == "/source/profiles/"


def test_current_physics_controls_are_serialized_and_round_trip(tmp_path: Path) -> None:
    base = configuration("electrostatic_periodic")
    config = base.model_copy(
        update={
            "physics": PhysicsConfig(
                electron_ifunc_conservation_model=1,
                ion_ifunc_conservation_model=3,
                ion_temperature_gradient_model="zero_Tprime",
            ),
            "run": base.run.model_copy(
                update={
                    "periodic": PeriodicConfig(
                        bparallel_ratio_real=0.25,
                        bparallel_ratio_imag=-0.5,
                    )
                }
            ),
        }
    )
    rendered = dumps_namelist(config)
    parsed = f90nml.reads(rendered)

    assert parsed["kim_config"]["electron_ifunc_conservation_model"] == 1
    assert parsed["kim_config"]["ion_ifunc_conservation_model"] == 3
    assert parsed["kim_config"]["ion_temperature_gradient_model"] == "zero_Tprime"
    assert parsed["kim_periodic"]["periodic_bparallel_ratio"] == "(0.25, -0.5)"

    path = tmp_path / "KIM_config.nml"
    path.write_text(rendered)
    assert load_namelist(path) == config


@pytest.mark.parametrize(
    "factory",
    [
        SimulationConfig.electrostatic_periodic,
        SimulationConfig.electrostatic,
        SimulationConfig.flr2,
    ],
)
def test_namelist_round_trip_preserves_supported_configuration(
    tmp_path: Path, factory: object
) -> None:
    expected = configuration(factory.__name__)
    path = tmp_path / "KIM_config.nml"
    path.write_text(dumps_namelist(expected))

    assert load_namelist(path) == expected
    assert SimulationConfig.from_namelist(path) == expected


def test_explicit_species_round_trip() -> None:
    config = SimulationConfig.electrostatic(
        **(
            physical_case()
            | {
                "plasma": ExplicitPlasma(
                    ions=(
                        IonSpecies(mass_number=2, charge_number=1),
                        IonSpecies(mass_number=4, charge_number=2),
                    )
                )
            }
        )
    )
    rendered = dumps_namelist(config)

    assert "&KIM_SPECIES" in rendered
    assert load_namelist_text(rendered) == config


def load_namelist_text(text: str) -> SimulationConfig:
    """Exercise the file-only public importer without retaining test files."""

    from tempfile import NamedTemporaryFile

    with NamedTemporaryFile(mode="w", suffix=".nml") as stream:
        stream.write(text)
        stream.flush()
        return load_namelist(Path(stream.name))


@pytest.mark.parametrize(
    ("group", "entry"),
    [
        ("UNKNOWN_GROUP", "value = 1"),
        ("KIM_GRID", "obsolete_grid_key = 1"),
        ("WKB_DISPERSION", "WKB_dispersion_mode = 'KIM'"),
        ("KIM_GRID", "theta_integration_method = 'RKF45'"),
        ("KIM_PROFILES", "equil_file = 'legacy.dat'"),
        ("KIM_PROFILES", "n_input_file = 'n_of_psiN.dat'"),
    ],
)
def test_import_rejects_unknown_and_obsolete_entries(group: str, entry: str) -> None:
    rendered = dumps_namelist(SimulationConfig.electrostatic(**physical_case()))
    if group == "UNKNOWN_GROUP":
        rendered += f"&{group}\n    {entry}\n/\n"
    else:
        marker = f"&{group}\n"
        rendered = rendered.replace(marker, marker + f"    {entry}\n")

    with pytest.raises(ConfigurationError, match="unsupported"):
        load_namelist_text(rendered)


def test_import_rejects_run_specific_groups_for_a_different_run_type() -> None:
    rendered = dumps_namelist(SimulationConfig.electrostatic(**physical_case()))
    rendered += "&KIM_PERIODIC\n    periodic_n_rg = 96\n/\n"

    with pytest.raises(ConfigurationError, match="not valid for run type"):
        load_namelist_text(rendered)


def test_import_rejects_inactive_species_group() -> None:
    rendered = dumps_namelist(SimulationConfig.electrostatic(**physical_case()))
    rendered += "&KIM_SPECIES\n    ai = 2\n    zi = 1\n/\n"

    with pytest.raises(ConfigurationError, match="not valid for built-in plasma"):
        load_namelist_text(rendered)


def test_import_rejects_hdf5_input() -> None:
    rendered = dumps_namelist(SimulationConfig.electrostatic(**physical_case()))
    rendered = rendered.replace("hdf5_input = .false.", "hdf5_input = .true.")

    with pytest.raises(ConfigurationError, match="hdf5_input"):
        load_namelist_text(rendered)


def test_import_rejects_inconsistent_profile_directories() -> None:
    rendered = dumps_namelist(SimulationConfig.electrostatic(**physical_case()))
    rendered = rendered.replace(
        "input_profile_dir = '/source/profiles/'",
        "input_profile_dir = '/different/profiles/'",
    )

    with pytest.raises(ConfigurationError, match="profile directories must match"):
        load_namelist_text(rendered)


def test_import_rejects_missing_required_group() -> None:
    rendered = dumps_namelist(SimulationConfig.electrostatic(**physical_case()))
    start = rendered.index("&KIM_SETUP")
    end = rendered.index("/\n", start) + 2

    with pytest.raises(ConfigurationError, match="KIM_SETUP"):
        load_namelist_text(rendered[:start] + rendered[end:])
