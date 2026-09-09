"""Validated configuration models for supported KIM simulations."""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Annotated, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    PositiveFloat,
    PositiveInt,
    ValidationInfo,
    field_validator,
    model_validator,
)


class KimModel(BaseModel):
    """Common validation policy for public KIM models."""

    model_config = ConfigDict(
        allow_inf_nan=False,
        extra="forbid",
        frozen=True,
        str_strip_whitespace=True,
        validate_default=True,
    )


class RunType(str, Enum):
    """Run types supported by the stable Python API."""

    ELECTROSTATIC_PERIODIC = "electrostatic_periodic"
    ELECTROSTATIC = "electrostatic"
    FLR2 = "flr2"


class CollisionModel(str, Enum):
    """Electron collision models implemented by KIM."""

    FOKKER_PLANCK = "FokkerPlanck"
    KROOK = "Krook"


class IonCollisionModel(str, Enum):
    """Ion collision models implemented by stable KIM run types."""

    FOKKER_PLANCK = "FokkerPlanck"
    COLLISIONLESS = "collisionless"


class PlasmaIsotope(str, Enum):
    """Built-in single-ion plasma definitions provided by KIM."""

    HYDROGEN = "H"
    DEUTERIUM = "D"


class GridSpacing(str, Enum):
    """Available radial grid-spacing algorithms."""

    EQUIDISTANT = "equidistant"
    NON_EQUIDISTANT = "non-equidistant"
    ADAPTIVE = "adaptive"


class ThetaIntegration(str, Enum):
    """Available theta-integration algorithms."""

    GAUSS_LEGENDRE = "GaussLegendre"
    RKF45 = "RKF45"
    QUADPACK = "QUADPACK"


class QuadpackAlgorithm(str, Enum):
    """Finite-interval QUADPACK algorithms supported by KIM."""

    QAG = "QAG"
    QAGS = "QAGS"


class IonSpecies(KimModel):
    """One explicitly configured ion species."""

    mass_number: PositiveInt = Field(
        description="Ion mass number in proton-mass units.",
        json_schema_extra={"units": "1", "sweepable": False},
    )
    charge_number: PositiveInt = Field(
        description="Positive ion charge number.",
        json_schema_extra={"units": "e", "sweepable": False},
    )

    @model_validator(mode="after")
    def charge_does_not_exceed_mass(self) -> IonSpecies:
        if self.charge_number > self.mass_number:
            raise ValueError("charge_number cannot exceed mass_number")
        return self


class BuiltinPlasma(KimModel):
    """A built-in hydrogen or deuterium plasma."""

    source: Literal["builtin"] = "builtin"
    isotope: PlasmaIsotope = Field(
        description="Built-in single-ion plasma isotope.",
        json_schema_extra={"sweepable": False},
    )


class ExplicitPlasma(KimModel):
    """An explicit list of one or more positive ion species."""

    source: Literal["explicit"] = "explicit"
    ions: tuple[IonSpecies, ...] = Field(
        min_length=1,
        max_length=16,
        description="Ion species written to the KIM_SPECIES namelist.",
        json_schema_extra={"sweepable": False},
    )


PlasmaConfig = Annotated[
    BuiltinPlasma | ExplicitPlasma,
    Field(discriminator="source"),
]


class PhysicsConfig(KimModel):
    """Collision and species-participation controls shared by stable runs."""

    collision_model: CollisionModel = Field(
        default=CollisionModel.FOKKER_PLANCK,
        description="Electron collision model.",
        json_schema_extra={"sweepable": False},
    )
    ion_collision_model: IonCollisionModel = Field(
        default=IonCollisionModel.FOKKER_PLANCK,
        description="Ion collision model.",
        json_schema_extra={"sweepable": False},
    )
    collisionless_kpar_epsilon: PositiveFloat | None = Field(
        default=None,
        description="Causal pole displacement required for collisionless ions.",
        json_schema_extra={"units": "1/cm", "sweepable": True},
    )
    ion_fp_collision_scale: PositiveFloat = Field(
        default=1.0,
        description="Multiplier for computed Fokker-Planck ion collision frequencies.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    collision_frequency_scale: PositiveFloat = Field(
        default=1.0,
        description="Multiplier for computed collision frequencies.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    density_scale: PositiveFloat | None = Field(
        default=None,
        description="Optional multiplier applied to number-density profiles.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    ion_larmor_radius_scale: float = Field(
        default=1.0,
        ge=0.0,
        description="Multiplier applied to the ion Larmor radius.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    artificial_debye_case: Literal[0, 1, 2] = Field(
        default=0,
        description="Artificial Debye-response diagnostic selector.",
        json_schema_extra={"sweepable": False},
    )
    conserve_energy: bool = Field(
        default=True,
        description="Apply the energy-conserving correction to Fokker-Planck I-functions.",
        json_schema_extra={"sweepable": False},
    )
    turn_off_ions: bool = Field(
        default=False,
        description="Exclude ions from the response calculation.",
        json_schema_extra={"sweepable": False},
    )
    turn_off_electrons: bool = Field(
        default=False,
        description="Exclude electrons from the response calculation.",
        json_schema_extra={"sweepable": False},
    )

    @model_validator(mode="after")
    def validate_collision_options(self) -> PhysicsConfig:
        if self.turn_off_ions and self.turn_off_electrons:
            raise ValueError("ions and electrons cannot both be disabled")
        if self.ion_collision_model is IonCollisionModel.COLLISIONLESS:
            if self.collision_model is not CollisionModel.FOKKER_PLANCK:
                raise ValueError(
                    "ion_collision_model=collisionless requires collision_model=FokkerPlanck"
                )
            if self.collisionless_kpar_epsilon is None:
                raise ValueError("collisionless_kpar_epsilon is required for collisionless ions")
        return self


class IOConfig(KimModel):
    """User-facing output controls for API-managed runs."""

    hdf5_input: Literal[False] = Field(
        default=False,
        description="HDF5 input is unsupported by the current Fortran core.",
        json_schema_extra={"sweepable": False},
    )
    hdf5_output: Literal[True] = Field(
        default=True,
        description="API-managed runs require structured HDF5 output.",
        json_schema_extra={"sweepable": False},
    )
    log_level: int = Field(
        default=3,
        ge=-1,
        le=5,
        description="Fortran logger level from silent (-1) through trace (5).",
        json_schema_extra={"sweepable": False},
    )
    data_verbosity: int = Field(
        default=1,
        ge=0,
        description="Amount of optional diagnostic data emitted by KIM.",
        json_schema_extra={"sweepable": False},
    )
    calculate_asymptotics: bool = Field(
        default=False,
        description="Calculate additional asymptotic diagnostics.",
        json_schema_extra={"sweepable": False},
    )
    write_diagnostics_dat: bool = Field(
        default=False,
        description="Write the legacy flat diagnostics file in addition to HDF5.",
        json_schema_extra={"sweepable": False},
    )


class SetupConfig(KimModel):
    """Required physical geometry, mode, frequency, and perturbation."""

    btor: float = Field(
        description="Signed toroidal magnetic field on the magnetic axis.",
        json_schema_extra={"units": "G", "sweepable": True},
    )
    major_radius: PositiveFloat = Field(
        description="Major radius of the magnetic axis.",
        json_schema_extra={"units": "cm", "sweepable": True},
    )
    m_mode: int = Field(
        description="Signed poloidal mode number.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    n_mode: int = Field(
        description="Signed toroidal mode number.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    frequency: float = Field(
        description="Angular frequency of the applied perturbation.",
        json_schema_extra={"units": "1/s", "sweepable": True},
    )
    br_boundary_real: float = Field(
        description="Real part of the constant radial magnetic perturbation.",
        json_schema_extra={"units": "G", "sweepable": True},
    )
    br_boundary_imag: float = Field(
        description="Imaginary part of the constant radial magnetic perturbation.",
        json_schema_extra={"units": "G", "sweepable": True},
    )
    spline_base: Literal[1] = Field(
        default=1,
        description="Spline basis selector; the stable API supports hat functions.",
        json_schema_extra={"sweepable": False},
    )
    magnetic_perturbation_type: Literal[12] = Field(
        default=12,
        description="Radial magnetic-field selector; the stable API supports a constant field.",
        json_schema_extra={"sweepable": False},
    )
    collisions_off: bool = Field(
        default=False,
        description="Disable collisional terms outside the collision-model selection.",
        json_schema_extra={"sweepable": False},
    )
    constant_profile_mode: Literal[0, 1, 2] = Field(
        default=0,
        description="Optional constant-profile diagnostic mode.",
        json_schema_extra={"sweepable": False},
    )
    boundary_condition: Literal[0, 1, 2, 3] = Field(
        default=3,
        description="Fortran boundary-condition selector.",
        json_schema_extra={"sweepable": False},
    )
    cyclotron_harmonics: int = Field(
        default=0,
        ge=0,
        description="Maximum cyclotron harmonic included on either side of zero.",
        json_schema_extra={"units": "1", "sweepable": True},
    )

    @field_validator("btor")
    @classmethod
    def magnetic_field_must_be_nonzero(cls, value: float) -> float:
        if value == 0.0:
            raise ValueError("btor must be nonzero")
        return value

    @field_validator("m_mode", "n_mode")
    @classmethod
    def mode_must_be_nonzero(cls, value: int, info: ValidationInfo) -> int:
        field_name = info.field_name
        if value == 0:
            raise ValueError(f"{field_name} must be nonzero")
        return value


class GridConfig(KimModel):
    """Radial-domain and numerical grid controls."""

    radial_minimum: PositiveFloat = Field(
        description="Minimum effective radius included in the calculation.",
        json_schema_extra={"units": "cm", "sweepable": True},
    )
    plasma_radius: PositiveFloat = Field(
        description="Maximum effective radius included in the calculation.",
        json_schema_extra={"units": "cm", "sweepable": True},
    )
    resonance_width: PositiveFloat = Field(
        default=0.5,
        description="Width parameter for resonance-focused grid refinement.",
        json_schema_extra={"units": "cm", "sweepable": True},
    )
    resonance_amplification: float = Field(
        default=15.0,
        ge=0.0,
        description="Strength of resonance-focused grid refinement.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    maximum_step_scale: PositiveFloat = Field(
        default=1.0,
        description="Scale factor for the largest generated radial step.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    l_space_dim: PositiveInt = Field(
        default=512,
        description="Requested spline-grid dimension.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    rg_space_dim: PositiveInt = Field(
        default=512,
        description="Requested background radial-grid dimension.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    rg_spacing: GridSpacing = Field(
        default=GridSpacing.EQUIDISTANT,
        description="Background radial-grid spacing algorithm.",
        json_schema_extra={"sweepable": False},
    )
    l_spacing: GridSpacing = Field(
        default=GridSpacing.EQUIDISTANT,
        description="Spline radial-grid spacing algorithm.",
        json_schema_extra={"sweepable": False},
    )
    theta_integration: ThetaIntegration = Field(
        default=ThetaIntegration.GAUSS_LEGENDRE,
        description="Theta-integration algorithm.",
        json_schema_extra={"sweepable": False},
    )
    larmor_skip_factor: PositiveFloat = Field(
        default=5.0,
        description="Distance cutoff scale used to omit negligible kernel entries.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    gauss_nodes_x: PositiveInt = Field(
        default=31,
        description="Gauss-Legendre nodes for the first radial coordinate.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    gauss_nodes_x_prime: PositiveInt = Field(
        default=30,
        description="Gauss-Legendre nodes for the second radial coordinate.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    gauss_nodes_theta: PositiveInt = Field(
        default=17,
        description="Gauss-Legendre nodes for theta integration.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    rkf45_absolute_tolerance: PositiveFloat = Field(
        default=1.0e-9,
        description="Absolute tolerance for adaptive RKF45 integration.",
        json_schema_extra={"sweepable": True},
    )
    rkf45_relative_tolerance: PositiveFloat = Field(
        default=1.0e-6,
        description="Relative tolerance for adaptive RKF45 integration.",
        json_schema_extra={"sweepable": True},
    )
    kernel_taper_skip_threshold: float = Field(
        default=1.0e-6,
        ge=0.0,
        le=1.0,
        description="Skip a kernel element below this taper weight.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    quadpack_algorithm: QuadpackAlgorithm = Field(
        default=QuadpackAlgorithm.QAG,
        description="Finite-interval QUADPACK algorithm.",
        json_schema_extra={"sweepable": False},
    )
    quadpack_key: int = Field(
        default=6,
        ge=1,
        le=6,
        description="Gauss-Kronrod rule selector for QAG.",
        json_schema_extra={"sweepable": False},
    )
    quadpack_limit: PositiveInt = Field(
        default=500,
        description="Maximum QUADPACK subdivisions.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    quadpack_absolute_tolerance: PositiveFloat = Field(
        default=1.0e-10,
        description="Absolute QUADPACK integration tolerance.",
        json_schema_extra={"sweepable": True},
    )
    quadpack_relative_tolerance: PositiveFloat = Field(
        default=1.0e-10,
        description="Relative QUADPACK integration tolerance.",
        json_schema_extra={"sweepable": True},
    )
    quadpack_use_u_substitution: bool = Field(
        default=True,
        description="Apply the u=sin(theta/2) substitution for QUADPACK.",
        json_schema_extra={"sweepable": False},
    )

    @model_validator(mode="after")
    def validate_grid(self) -> GridConfig:
        if self.radial_minimum >= self.plasma_radius:
            raise ValueError("radial_minimum must be smaller than plasma_radius")
        if self.gauss_nodes_x == self.gauss_nodes_x_prime:
            raise ValueError("gauss_nodes_x and gauss_nodes_x_prime must differ")
        return self


class ProfileConfig(KimModel):
    """File names and source directory for `r_eff` profiles."""

    directory: Path = Field(
        description="Directory containing the source profile files.",
        json_schema_extra={"sweepable": False},
    )
    coordinate_type: Literal["r_eff"] = Field(
        default="r_eff",
        description="Profile radial coordinate supported by the MVP.",
        json_schema_extra={"units": "cm", "sweepable": False},
    )
    density_file: str = Field(
        default="n.dat",
        min_length=1,
        description="Electron-density profile filename.",
        json_schema_extra={"units": "1/cm^3", "sweepable": False},
    )
    electron_temperature_file: str = Field(
        default="Te.dat",
        min_length=1,
        description="Electron-temperature profile filename.",
        json_schema_extra={"units": "eV", "sweepable": False},
    )
    ion_temperature_file: str = Field(
        default="Ti.dat",
        min_length=1,
        description="Ion-temperature profile filename.",
        json_schema_extra={"units": "eV", "sweepable": False},
    )
    toroidal_velocity_file: str = Field(
        default="Vz.dat",
        min_length=1,
        description="Toroidal-velocity profile filename.",
        json_schema_extra={"units": "cm/s", "sweepable": False},
    )
    radial_electric_field_file: str = Field(
        default="Er.dat",
        min_length=1,
        description="Radial-electric-field profile filename.",
        json_schema_extra={"units": "statV/cm", "sweepable": False},
    )
    safety_factor_file: str = Field(
        default="q.dat",
        min_length=1,
        description="Safety-factor profile filename.",
        json_schema_extra={"units": "1", "sweepable": False},
    )


class PeriodicConfig(KimModel):
    """Numerical controls for the forced-periodicity solver."""

    as_is_width_scale: PositiveFloat = Field(
        default=5.0,
        description="Unmodified half-width in resonant-species Larmor radii.",
        json_schema_extra={"units": "rho_ref", "sweepable": True},
    )
    transition_width_scale: PositiveFloat = Field(
        default=10.0,
        description="Periodization transition width in Larmor radii per side.",
        json_schema_extra={"units": "rho_ref", "sweepable": True},
    )
    wavenumber_cutoff_scale: PositiveFloat = Field(
        default=5.0,
        description="Fourier cutoff multiplied by the reference Larmor radius.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    n_rg: PositiveInt = Field(
        default=96,
        description="Number of periodic-window radial boundary points.",
        json_schema_extra={"units": "1", "sweepable": True},
    )
    match_global_kernel_approximations: bool = Field(
        default=False,
        description="Use approximations intended only for global-kernel comparisons.",
        json_schema_extra={"sweepable": False},
    )


class Flr2Config(KimModel):
    """Term switches for the standalone FLR2 response."""

    electron_flr: bool = Field(
        default=True,
        description="Include electron finite-Larmor-radius terms.",
        json_schema_extra={"sweepable": False},
    )
    ion_flr: bool = Field(
        default=True,
        description="Include ion finite-Larmor-radius terms.",
        json_schema_extra={"sweepable": False},
    )
    electron_potential: bool = Field(
        default=True,
        description="Include the electron response in the potential equation.",
        json_schema_extra={"sweepable": False},
    )
    ion_potential: bool = Field(
        default=True,
        description="Include the ion response in the potential equation.",
        json_schema_extra={"sweepable": False},
    )
    electron_current: bool = Field(
        default=True,
        description="Include the electron contribution to parallel current.",
        json_schema_extra={"sweepable": False},
    )
    ion_current: bool = Field(
        default=True,
        description="Include the ion contribution to parallel current.",
        json_schema_extra={"sweepable": False},
    )
    include_potential_in_current: bool = Field(
        default=True,
        description="Include potential-response terms in parallel current.",
        json_schema_extra={"sweepable": False},
    )


class ElectrostaticRun(KimModel):
    """Global electrostatic run selection."""

    run_type: Literal[RunType.ELECTROSTATIC] = Field(
        default=RunType.ELECTROSTATIC,
        description="Select the global electrostatic solver.",
    )


class ElectrostaticPeriodicRun(KimModel):
    """Forced-periodicity electrostatic run selection and controls."""

    run_type: Literal[RunType.ELECTROSTATIC_PERIODIC] = Field(
        default=RunType.ELECTROSTATIC_PERIODIC,
        description="Select the forced-periodicity electrostatic solver.",
    )
    periodic: PeriodicConfig = Field(default_factory=PeriodicConfig)


class Flr2Run(KimModel):
    """Standalone FLR2 run selection and controls."""

    run_type: Literal[RunType.FLR2] = Field(
        default=RunType.FLR2,
        description="Select the standalone FLR2 solver.",
    )
    terms: Flr2Config = Field(default_factory=Flr2Config)


RunConfig = Annotated[
    ElectrostaticRun | ElectrostaticPeriodicRun | Flr2Run,
    Field(discriminator="run_type"),
]


class SimulationConfig(KimModel):
    """Complete validated request for one supported KIM simulation."""

    run: RunConfig
    plasma: PlasmaConfig
    setup: SetupConfig
    grid: GridConfig
    profiles: ProfileConfig
    physics: PhysicsConfig = Field(default_factory=PhysicsConfig)
    io: IOConfig = Field(default_factory=IOConfig)

    @model_validator(mode="after")
    def validate_cross_group_constraints(self) -> SimulationConfig:
        if (
            self.setup.collisions_off
            and self.physics.collision_model is CollisionModel.FOKKER_PLANCK
        ):
            raise ValueError("collisions_off cannot be true when collision_model is FokkerPlanck")
        if (
            self.physics.ion_collision_model is IonCollisionModel.COLLISIONLESS
            and self.grid.theta_integration is not ThetaIntegration.GAUSS_LEGENDRE
        ):
            raise ValueError("collisionless ion model requires GaussLegendre theta integration")
        if (
            self.physics.ion_collision_model is IonCollisionModel.COLLISIONLESS
            and self.physics.artificial_debye_case != 0
        ):
            raise ValueError("collisionless ions require artificial_debye_case=0")
        if isinstance(self.run, Flr2Run):
            self._validate_flr2()
        return self

    def _validate_flr2(self) -> None:
        if self.physics.collision_model is not CollisionModel.FOKKER_PLANCK:
            raise ValueError("FLR2 requires collision_model=FokkerPlanck")
        if self.setup.frequency != 0.0:
            raise ValueError("FLR2 currently requires frequency = 0")
        if self.number_of_ion_species != 1:
            raise ValueError("FLR2 currently requires exactly one ion species")

        electron_potential = (
            self.run.terms.electron_potential and not self.physics.turn_off_electrons
        )
        ion_potential = self.run.terms.ion_potential and not self.physics.turn_off_ions
        if not electron_potential and not ion_potential:
            raise ValueError("FLR2 requires at least one species in the potential equation")

    @property
    def number_of_ion_species(self) -> int:
        if isinstance(self.plasma, BuiltinPlasma):
            return 1
        return len(self.plasma.ions)

    @classmethod
    def electrostatic_periodic(
        cls,
        *,
        profiles: Path,
        plasma: BuiltinPlasma | ExplicitPlasma,
        btor: float,
        major_radius: float,
        m_mode: int,
        n_mode: int,
        frequency: float,
        br_boundary_real: float,
        br_boundary_imag: float,
        radial_minimum: float,
        plasma_radius: float,
        physics: PhysicsConfig | None = None,
        io: IOConfig | None = None,
        grid: GridConfig | None = None,
        periodic: PeriodicConfig | None = None,
    ) -> SimulationConfig:
        """Build a complete forced-periodicity electrostatic request."""

        return cls._from_run(
            ElectrostaticPeriodicRun(periodic=periodic or PeriodicConfig()),
            profiles=profiles,
            plasma=plasma,
            btor=btor,
            major_radius=major_radius,
            m_mode=m_mode,
            n_mode=n_mode,
            frequency=frequency,
            br_boundary_real=br_boundary_real,
            br_boundary_imag=br_boundary_imag,
            radial_minimum=radial_minimum,
            plasma_radius=plasma_radius,
            physics=physics,
            io=io,
            grid=grid,
        )

    @classmethod
    def electrostatic(
        cls,
        *,
        profiles: Path,
        plasma: BuiltinPlasma | ExplicitPlasma,
        btor: float,
        major_radius: float,
        m_mode: int,
        n_mode: int,
        frequency: float,
        br_boundary_real: float,
        br_boundary_imag: float,
        radial_minimum: float,
        plasma_radius: float,
        physics: PhysicsConfig | None = None,
        io: IOConfig | None = None,
        grid: GridConfig | None = None,
    ) -> SimulationConfig:
        """Build a complete global electrostatic request."""

        return cls._from_run(
            ElectrostaticRun(),
            profiles=profiles,
            plasma=plasma,
            btor=btor,
            major_radius=major_radius,
            m_mode=m_mode,
            n_mode=n_mode,
            frequency=frequency,
            br_boundary_real=br_boundary_real,
            br_boundary_imag=br_boundary_imag,
            radial_minimum=radial_minimum,
            plasma_radius=plasma_radius,
            physics=physics,
            io=io,
            grid=grid,
        )

    @classmethod
    def flr2(
        cls,
        *,
        profiles: Path,
        plasma: BuiltinPlasma | ExplicitPlasma,
        btor: float,
        major_radius: float,
        m_mode: int,
        n_mode: int,
        frequency: float,
        br_boundary_real: float,
        br_boundary_imag: float,
        radial_minimum: float,
        plasma_radius: float,
        physics: PhysicsConfig | None = None,
        io: IOConfig | None = None,
        grid: GridConfig | None = None,
        terms: Flr2Config | None = None,
    ) -> SimulationConfig:
        """Build a complete standalone FLR2 request."""

        return cls._from_run(
            Flr2Run(terms=terms or Flr2Config()),
            profiles=profiles,
            plasma=plasma,
            btor=btor,
            major_radius=major_radius,
            m_mode=m_mode,
            n_mode=n_mode,
            frequency=frequency,
            br_boundary_real=br_boundary_real,
            br_boundary_imag=br_boundary_imag,
            radial_minimum=radial_minimum,
            plasma_radius=plasma_radius,
            physics=physics,
            io=io,
            grid=grid,
        )

    @classmethod
    def _from_run(
        cls,
        run: ElectrostaticRun | ElectrostaticPeriodicRun | Flr2Run,
        *,
        profiles: Path,
        plasma: BuiltinPlasma | ExplicitPlasma,
        btor: float,
        major_radius: float,
        m_mode: int,
        n_mode: int,
        frequency: float,
        br_boundary_real: float,
        br_boundary_imag: float,
        radial_minimum: float,
        plasma_radius: float,
        physics: PhysicsConfig | None,
        io: IOConfig | None,
        grid: GridConfig | None,
    ) -> SimulationConfig:
        if grid is None:
            grid = GridConfig(
                radial_minimum=radial_minimum,
                plasma_radius=plasma_radius,
            )
        elif grid.radial_minimum != radial_minimum or grid.plasma_radius != plasma_radius:
            raise ValueError("factory radial bounds must match the supplied GridConfig")

        return cls(
            run=run,
            plasma=plasma,
            setup=SetupConfig(
                btor=btor,
                major_radius=major_radius,
                m_mode=m_mode,
                n_mode=n_mode,
                frequency=frequency,
                br_boundary_real=br_boundary_real,
                br_boundary_imag=br_boundary_imag,
            ),
            grid=grid,
            profiles=ProfileConfig(directory=profiles),
            physics=physics or PhysicsConfig(),
            io=io or IOConfig(),
        )
