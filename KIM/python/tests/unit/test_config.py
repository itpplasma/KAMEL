from pathlib import Path

import pytest
from kim.config import (
    BuiltinPlasma,
    CollisionModel,
    ElectrostaticPeriodicRun,
    ElectrostaticRun,
    ExplicitPlasma,
    Flr2Run,
    GridConfig,
    IOConfig,
    IonCollisionModel,
    IonSpecies,
    PeriodicConfig,
    PhysicsConfig,
    PlasmaIsotope,
    ProfileConfig,
    RunType,
    SetupConfig,
    SimulationConfig,
    ThetaIntegration,
)
from kim.errors import ConfigurationError, KimError
from pydantic import ValidationError


def physical_case() -> dict[str, object]:
    return {
        "profiles": Path("/scientific-input/profiles"),
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


def test_physical_case_parameters_are_required() -> None:
    with pytest.raises(ValidationError) as error:
        SetupConfig()

    missing = {item["loc"][0] for item in error.value.errors()}

    assert missing == {
        "btor",
        "major_radius",
        "m_mode",
        "n_mode",
        "frequency",
        "br_boundary_real",
        "br_boundary_imag",
    }

    with pytest.raises(ValidationError, match="Field required"):
        GridConfig()

    with pytest.raises(ValidationError, match="Field required"):
        ProfileConfig()

    with pytest.raises(ValidationError, match="Field required"):
        BuiltinPlasma()


@pytest.mark.parametrize(
    ("factory", "expected_type", "expected_run_class"),
    [
        (
            SimulationConfig.electrostatic_periodic,
            RunType.ELECTROSTATIC_PERIODIC,
            ElectrostaticPeriodicRun,
        ),
        (
            SimulationConfig.electrostatic,
            RunType.ELECTROSTATIC,
            ElectrostaticRun,
        ),
        (SimulationConfig.flr2, RunType.FLR2, Flr2Run),
    ],
)
def test_stable_run_factories_build_complete_configs(
    factory: object,
    expected_type: RunType,
    expected_run_class: type[object],
) -> None:
    config = factory(**physical_case())

    assert config.run.run_type is expected_type
    assert isinstance(config.run, expected_run_class)
    assert config.profiles.coordinate_type == "r_eff"
    assert config.io.hdf5_output is True
    assert config.io.hdf5_input is False


def test_unknown_run_type_is_rejected() -> None:
    valid = SimulationConfig.electrostatic(**physical_case()).model_dump(mode="json")
    valid["run"]["run_type"] = "electromagnetic"

    with pytest.raises(ValidationError, match="union_tag_invalid"):
        SimulationConfig.model_validate(valid)


def test_unknown_fields_are_rejected_at_every_level() -> None:
    valid = SimulationConfig.electrostatic(**physical_case()).model_dump(mode="json")
    valid["unknown_parameter"] = 1

    with pytest.raises(ValidationError, match="extra_forbidden"):
        SimulationConfig.model_validate(valid)

    with pytest.raises(ValidationError, match="extra_forbidden"):
        GridConfig(radial_minimum=3.0, plasma_radius=63.0, unknown_parameter=1)


@pytest.mark.parametrize("field", ["hdf5_input", "hdf5_output"])
def test_api_managed_runs_force_supported_hdf5_modes(field: str) -> None:
    value = field == "hdf5_input"

    with pytest.raises(ValidationError):
        IOConfig(**{field: value})


@pytest.mark.parametrize(
    "changes",
    [
        {"l_space_dim": 0},
        {"rg_space_dim": -1},
        {"gauss_nodes_x": 0},
        {"gauss_nodes_x_prime": -1},
        {"gauss_nodes_theta": 0},
    ],
)
def test_grid_dimensions_must_be_positive(changes: dict[str, int]) -> None:
    with pytest.raises(ValidationError, match="greater than 0"):
        GridConfig(radial_minimum=3.0, plasma_radius=63.0, **changes)


def test_radial_minimum_must_be_inside_plasma_radius() -> None:
    with pytest.raises(ValidationError, match="radial_minimum must be smaller"):
        GridConfig(radial_minimum=63.0, plasma_radius=63.0)


@pytest.mark.parametrize(("field", "value"), [("m_mode", 0), ("n_mode", 0)])
def test_mode_numbers_must_be_nonzero(field: str, value: int) -> None:
    values = {
        "btor": -17_977.413,
        "major_radius": 165.0,
        "m_mode": 7,
        "n_mode": 2,
        "frequency": 0.0,
        "br_boundary_real": 1.0,
        "br_boundary_imag": 0.0,
    }
    values[field] = value

    with pytest.raises(ValidationError, match="must be nonzero"):
        SetupConfig(**values)


def test_nonfinite_numeric_parameters_are_rejected() -> None:
    values = {
        "btor": float("nan"),
        "major_radius": 165.0,
        "m_mode": 7,
        "n_mode": 2,
        "frequency": 0.0,
        "br_boundary_real": 1.0,
        "br_boundary_imag": 0.0,
    }

    with pytest.raises(ValidationError, match="finite number"):
        SetupConfig(**values)


def test_invalid_enumerated_values_are_rejected() -> None:
    with pytest.raises(ValidationError):
        PhysicsConfig(collision_model="BGK")

    with pytest.raises(ValidationError):
        PhysicsConfig(ion_collision_model="Krook")

    with pytest.raises(ValidationError):
        GridConfig(
            radial_minimum=3.0,
            plasma_radius=63.0,
            theta_integration="trapezoidal",
        )


def test_explicit_ion_species_are_validated() -> None:
    with pytest.raises(ValidationError, match="charge_number cannot exceed mass_number"):
        IonSpecies(mass_number=2, charge_number=3)

    with pytest.raises(ValidationError, match="at least 1 item"):
        ExplicitPlasma(ions=())


def test_collisionless_ions_require_supported_physics() -> None:
    with pytest.raises(ValidationError, match="collisionless_kpar_epsilon"):
        PhysicsConfig(
            ion_collision_model=IonCollisionModel.COLLISIONLESS,
        )

    with pytest.raises(ValidationError, match="requires collision_model=FokkerPlanck"):
        PhysicsConfig(
            collision_model=CollisionModel.KROOK,
            ion_collision_model=IonCollisionModel.COLLISIONLESS,
            collisionless_kpar_epsilon=1.0e-5,
        )


def test_collisionless_ions_require_gauss_legendre_integration() -> None:
    config = SimulationConfig.electrostatic_periodic(
        **physical_case(),
        physics=PhysicsConfig(
            ion_collision_model=IonCollisionModel.COLLISIONLESS,
            collisionless_kpar_epsilon=1.0e-5,
        ),
        grid=GridConfig(
            radial_minimum=3.0,
            plasma_radius=63.0,
            theta_integration=ThetaIntegration.GAUSS_LEGENDRE,
        ),
    )
    invalid = config.model_dump(mode="json")
    invalid["grid"]["theta_integration"] = "RKF45"

    with pytest.raises(ValidationError, match="requires GaussLegendre"):
        SimulationConfig.model_validate(invalid)


def test_fokker_planck_cannot_be_combined_with_collisions_off() -> None:
    valid = SimulationConfig.electrostatic(**physical_case()).model_dump(mode="json")
    valid["setup"]["collisions_off"] = True

    with pytest.raises(ValidationError, match="collisions_off"):
        SimulationConfig.model_validate(valid)


def test_both_species_cannot_be_disabled() -> None:
    with pytest.raises(ValidationError, match="cannot both be disabled"):
        PhysicsConfig(
            turn_off_ions=True,
            turn_off_electrons=True,
        )


def test_flr2_enforces_current_fortran_constraints() -> None:
    base = SimulationConfig.flr2(**physical_case()).model_dump(mode="json")

    nonzero_frequency = base | {"setup": base["setup"] | {"frequency": 1.0}}
    with pytest.raises(ValidationError, match="frequency = 0"):
        SimulationConfig.model_validate(nonzero_frequency)

    two_ions = base | {
        "plasma": {
            "source": "explicit",
            "ions": [
                {"mass_number": 2, "charge_number": 1},
                {"mass_number": 3, "charge_number": 1},
            ],
        }
    }
    with pytest.raises(ValidationError, match="exactly one ion species"):
        SimulationConfig.model_validate(two_ions)


def test_periodic_dimensions_are_positive() -> None:
    with pytest.raises(ValidationError, match="greater than 0"):
        PeriodicConfig(n_rg=0)


def test_json_schema_exposes_units_descriptions_and_sweepability() -> None:
    setup_schema = SetupConfig.model_json_schema()["properties"]
    grid_schema = GridConfig.model_json_schema()["properties"]
    periodic_schema = PeriodicConfig.model_json_schema()["properties"]

    assert setup_schema["btor"]["units"] == "G"
    assert setup_schema["btor"]["description"]
    assert grid_schema["radial_minimum"]["units"] == "cm"
    assert grid_schema["l_space_dim"]["sweepable"] is True
    assert periodic_schema["n_rg"]["sweepable"] is True


def test_configuration_errors_share_a_public_base_class() -> None:
    assert issubclass(ConfigurationError, KimError)
