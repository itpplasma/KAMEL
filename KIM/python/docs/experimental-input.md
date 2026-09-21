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

## CLI Walkthrough

The same preparation is available without writing Python. Store the source declarations in a
version-controlled metadata file such as `marsf-metadata.json`:

```json
{
  "source": "mars-f-shot-12345",
  "coordinate": "sqrt_psiN",
  "coordinate_unit": "1",
  "density_unit": "1/m^3",
  "electron_temperature_unit": "eV",
  "ion_temperature_unit": "eV",
  "toroidal_velocity_unit": "m/s",
  "equilibrium_provenance": "equilibrium-2026-09-18"
}
```

Prepare a new case from an existing equilibrium table:

```bash
kim prepare-marsf ./marsf-input request.json ./prepared-case \
  --metadata ./marsf-metadata.json \
  --equilibrium-file ./equil_r_q_psi.dat \
  --format json
cd ./prepared-case
kim validate request.json
```

The command does not launch `KIM.x`. It prints the staged paths and writes a request whose profile
directory is `./profiles`, so validation and later run commands work from inside the prepared case.
For an equilibrium that must be traced by KAMEL, replace `--equilibrium-file` with the existing
Fortran preprocessor and its control files:

```bash
kim prepare-marsf ./marsf-input request.json ./prepared-case \
  --metadata ./marsf-metadata.json \
  --equilibrium-executable ./build/fouriermodes.x \
  --equilibrium-input ./field_divB0.inp \
  --equilibrium-input ./fouriermodes.inp \
  --equilibrium-timeout 3600 \
  --format json
```

The generated `conversion_report.json` is the preparation provenance record:

- `schema_version` identifies the report schema; `created_at` records its UTC creation time.
- `source` repeats every explicit metadata declaration supplied by the caller.
- `target` names the KIM representation, currently `KIM-CGS-r_eff`.
- `equilibrium` records the original table path when one was supplied, the staged relative path,
  and whether the table was copied or generated.
- `coordinate_operation` describes preservation or explicit `sqrt_psiN`/`r_eff` mapping.
- `source_hashes` contains lowercase SHA-256 values. `source/` identifies copied MARS-F files,
  `equilibrium/` identifies the equilibrium table, and `equilibrium_input/` identifies generator
  control files.
- `operations` lists each quantity's source unit, target unit, and scalar factor.
- `output_grid_points` records the number of rows written to each prepared profile.
- `generator` is `null` for a supplied table. Otherwise it records the command, executable path and
  hash, hashes for command files, and the configured timeout.

The prepared case stores copied MARS-F files and `metadata.json` under `source/`, equilibrium
artifacts and generator control files under `equilibrium/`, converted profiles under `profiles/`,
and the request/report at the case root. Restricted or private inputs are supplied by path at
invocation time and are copied only when the caller is authorized to stage them. The CLI downloads
nothing and performs no format, unit, coordinate, or equilibrium inference.
