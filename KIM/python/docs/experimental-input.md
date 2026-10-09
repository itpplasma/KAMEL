# Experimental Input Preparation

For a ready-to-edit Python script, use [`prepare_marsf.py`](../examples/prepare_marsf.py).
It prepares inputs only and needs no reference HDF5 file. From the repository root:

```bash
python -m pip install -e './KIM/python'
cp KIM/python/examples/request-periodic.json ./request.json
# Edit request.json for your plasma, fields, mode, radial domain, and resolution.
# Edit the four paths and source metadata at the top of prepare_marsf.py.
python KIM/python/examples/prepare_marsf.py
cd prepared-case
kim validate request.json
```

The four paths identify the MARS-F profile directory, a precomputed `equil_r_q_psi.dat`,
your JSON request, and a new output directory. Relative paths use your working directory.
The equilibrium table's first three columns are `r_eff [cm]`, q, and poloidal flux psi;
header/comment lines start with `#`. Use the existing KAMEL equilibrium calculation to
produce it. Preparation normalizes psi by its final tabulated value.
The script preserves its q sign. The request supplies magnetic field and major radius;
these must be consistent with your equilibrium. The source metadata explicitly declares
coordinates, units, and provenance. Confirm those declarations before using your own data;
the provided declarations are examples, not an inference about your files. No solver runs
during preparation. Existing output directories are never replaced.

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

## BALANCE Preparation and Optional Comparison

`characterize_balance` prepares a KIM case from four explicit BALANCE profile paths,
metadata, a KIM configuration, and an equilibrium calculation. An existing QL-Balance HDF5
file is **optional**. In Python, omit `oracle` (or pass `None`) for preparation only.
No comparison domains, interpolation settings, relative floors, or tolerances are needed.

The CLI accepts the destination immediately after the configuration:

```bash
kim characterize-balance density.dat Te.dat Ti.dat vt.dat metadata.json request.json \
  ./prepared-balance \
  --equilibrium-file ./equilibrium/equil_r_q_psi.dat \
  --equilibrium-parameters-file ./equilibrium/btor_rbig.dat \
  --equilibrium-provenance selected-equilibrium
```

Without a reference, success is reported as `PREPARED`, with `comparison_performed: false`,
empty `comparisons`, and `null` comparison configuration, threshold decisions, and overall pass.
The prepared profiles, request, conversion reports, and input provenance are still written.
Preparation does not require or access a reference HDF5 file.

To compare as well, pass `oracle=Path(...)` in Python or add these CLI options:

```bash
--reference-hdf5 ./reference.hdf5 --domains '{"core":[3,60]}' \
  --interpolation-direction prepared_to_oracle --interpolation-method linear
```

The illustrated domain is an example, not an approved AUG comparison domain. Choose domains
explicitly for your case. Comparison settings require a reference; a missing or invalid supplied
reference is an error, not a request to skip comparison. The old CLI form with
`REFERENCE_HDF5 DESTINATION` after the configuration remains supported, as does the Python
positional argument order. Expected hashes, when supplied, cover every supplied input;
the `oracle` hash is required only when a reference is supplied.

### Equilibrium Values

The `characterize-balance` API and CLI take `btor` and `r_big` from the equilibrium calculation.
On the production route, the selected equilibrium executable runs once and must produce both
`equil_r_q_psi.dat` and `btor_rbig.dat`. On the precomputed route, provide both files with
`equilibrium_file` and `equilibrium_parameters_file` (CLI options `--equilibrium-file` and
`--equilibrium-parameters-file`). They must use those canonical filenames and reside in the same
calculation output directory, which the caller identifies with `equilibrium_provenance`.

The two values in `btor_rbig.dat` replace `config.setup.btor` and `config.setup.major_radius` for
the staged KIM request. Its `r_big` is also used for the BALANCE conversion
`v_phi [cm/s] = r_big [cm] * omega [rad/s]`. The characterization report records the values and
hashes of both equilibrium outputs, plus the generator provenance when the calculation is run by
the pipeline. The former `major_radius_cm` API argument and `--major-radius-cm` CLI option remain
available as deprecated consistency checks; if supplied, the value must exactly equal the
calculation's `r_big` and is never used as an input.

The BALANCE `profiles/Vz.dat` interface expects toroidal velocity in `cm/s`. Its reader loads those
values directly; BALANCE converts the linear velocity to its internal rotation frequency by
dividing by `rtor`, and converts back to `cm/s` when updating KIM. For the AUG `vt` source files,
the user confirmed angular rotation in `rad/s`; the pipeline converts it with the selected
calculation's `r_big` before staging the BALANCE velocity profile. The two-column source file has
no embedded unit label, so the confirmed unit is recorded in case metadata.

BALANCE characterization preserves q exactly as `fouriermodes.x` calculates it from the selected
EQDSK. EQDSK coordinate conventions can change q's sign, so this route does not impose a fixed sign
change. The selected AUG MICDU file is read in EFIT format, and its Fouriers output is used without
modification.

When a reference is supplied, characterization reports absolute RMS and maximum errors.
Relative RMS and maximum require
caller-supplied per-profile denominator floors; if `relative_floors` or CLI `--relative-floors` is
omitted, those fields are `null` and the report warns that relative metrics are unavailable. A
relative tolerance is rejected without floors. No pass/fail result is produced unless tolerances
are explicitly supplied.

The oracle's required datasets must be stored in its HDF5 file. External links, virtual datasets,
and external raw storage are rejected because the file hash cannot identify those dependencies.
Same-file aliases are supported. RMS metrics describe evaluated grid samples, without radial
quadrature weighting. Resonance coverage flags describe the requested domain; use the reported
exclusions and comparison radii to assess shared support at a crossing.

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
- `coordinate_mapping` adds structured source/target coordinates and method while the version-1
  `coordinate_operation` string remains compatible with existing readers.
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
