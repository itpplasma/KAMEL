# Experimental Input Preparation

The experimental importer is read-only. To turn a MARS-F profile quartet into a runnable KIM
case, declare the source conventions and call `prepare_marsf_case`:

```python
from pathlib import Path

from kim import (
    BuiltinPlasma,
    MarsFMetadata,
    PlasmaIsotope,
    SimulationConfig,
    prepare_marsf_case,
    read_marsf_profiles,
)

metadata = MarsFMetadata(
    source="mars-f-shot-12345",
    coordinate="sqrt_psiN",
    coordinate_unit="1",
    density_unit="1/m^3",
    electron_temperature_unit="eV",
    ion_temperature_unit="eV",
    toroidal_velocity_unit="m/s",
    equilibrium_provenance="equilibrium-2026-09-18",
)
source = read_marsf_profiles("./marsf-input", metadata)
config = SimulationConfig.electrostatic_periodic(
    profiles=Path("unused"),
    plasma=BuiltinPlasma(isotope=PlasmaIsotope.DEUTERIUM),
    btor=-17977.413,
    major_radius=165.0,
    m_mode=7,
    n_mode=2,
    frequency=0.0,
    br_boundary_real=1.0,
    br_boundary_imag=0.0,
    radial_minimum=3.0,
    plasma_radius=63.0,
)
case = prepare_marsf_case(
    source,
    config,
    "./prepared-case",
    equilibrium_file="./equil_r_q_psi.dat",
)
```

Preparation refuses to replace an existing destination. It copies the untouched source files,
writes `n.dat`, `Te.dat`, `Ti.dat`, `Vz.dat`, and `q.dat` in KIM CGS units, writes the validated
`request.json`, and records source hashes plus every coordinate/unit operation in
`conversion_report.json`. `Er.dat` is intentionally omitted so KIM can calculate radial force
balance from the prepared profiles.

For an equilibrium that has not already been reduced to `equil_r_q_psi.dat`, provide the existing
KAMEL equilibrium executable and its input files. The executable runs without a shell in the
staging directory and must write `equil_r_q_psi.dat` there:

```python
case = prepare_marsf_case(
    source,
    config,
    "./prepared-case",
    equilibrium_executable="./build/install/bin/fouriermodes.x",
    equilibrium_input_files=("./field_divB0.inp", "./fouriermodes.inp"),
)
```

Coordinate mapping uses the explicit equilibrium table and a natural cubic spline. Extrapolation,
implicit unit inference, unsupported temperature units, and ambiguous `r_eff` grids are rejected.
