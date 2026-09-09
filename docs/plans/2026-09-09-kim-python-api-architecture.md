# KIM Python API and CLI Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a stable, validated Python automation layer and thin CLI around the authoritative
KIM Fortran executable, with reproducible runs and interfaces suitable for later MCP tools.

**Architecture:** A separately installable Python distribution will own the `kim` import namespace
and `kim` command while invoking `KIM.x` as a subprocess. Typed configuration models will generate
canonical KIM namelists, stage copied `r_eff` profiles, manage immutable run directories, and expose
structured HDF5 results without reproducing scientific calculations in Python.

**Tech Stack:** Python 3.10+, Pydantic 2, Typer, f90nml, NumPy, h5py, pytest, the existing CMake/Ninja
Fortran build, and the existing KIM HDF5 output.

---

## Status and agreed decisions

This document records the architecture agreed after inspecting the KAMEL repository, the current
Fortran implementation, the existing Python code, and the publication Er-scan workflow in
`~/proj/RMP_Er_scan_shielding_KIM`.

The following decisions are settled:

- The scientific implementation remains in Fortran and is executed through `KIM.x`.
- Python will not bind individual Fortran procedures with f2py, ctypes, or similar mechanisms.
- The Python CLI owns the command name `kim`; the scientific executable remains `KIM.x`.
- KIM receives a separate Python distribution rather than being added to the broad KAMELpy
  distribution.
- The distribution supports Python 3.10 and newer.
- Complete typed Python models are the primary configuration interface.
- Importing an existing namelist is a secondary migration and interoperability path.
- The first stable profile coordinate is `r_eff`.
- API-managed runs use HDF5 output by default and retain all raw solver artifacts.
- The MVP supports `electrostatic_periodic`, `electrostatic`, and `flr2`.
- `electrostatic_periodic` is the priority because it is fast and was the main workflow behind the
  recent publication.
- Obsolete Python orchestration and GUI code may be removed after useful behavior has replacement
  coverage.

## What the repository does today

`make KIM` builds `build/install/bin/KIM.x`. The executable accepts zero or one positional namelist
path. It reads the namelist groups in a fixed sequence and resolves relative input and output paths
against its working directory. It appends `m<m>_n<n>/` beneath the configured output directory.

The current run contract has several implications for Python:

- A complete input deck must be staged before launching the process.
- The subprocess working directory is part of the scientific input and must be controlled.
- Profiles must be copied into the run because KIM may create `Er.dat` and `Er_no_Vpol.dat` when the
  electric-field profile is absent.
- Exit status alone is not a sufficient success condition because Fortran failure paths are mixed,
  and an HDF5 file may be created before a periodic solve fails.
- `hdf5_input = .true.` is currently unimplemented and must not be exposed as a supported option.
- The primary density, temperature, and safety-factor readers expect matching radial grids. Only
  the Er reader currently performs interpolation.

For `electrostatic_periodic`, the publication campaign establishes the useful output contract:

- `fields/Phi`: complex electrostatic potential on the periodic window grid.
- `fields/jpar`: complex total parallel-current density on that grid.
- `fields/jpar_e`, `fields/jpar_i`, and per-ion datasets when emitted by the executable.
- `backs/e/r`: the periodic window grid in cm.
- `backs/E0r`: Er on the periodic window grid in statV/cm.
- `setup/periodic_scale/M` and `setup/periodic_scale/dx_asis`: derived periodic geometry.
- `runtime`: solver runtime.

The scientifically meaningful campaign integral is

```text
I_parallel = 2 pi integral(r * j_parallel(r) dr)
```

over `abs(r - r_res) <= dx_asis`. The periodization transition zones contain boundary artifacts, so
the full-window integral must remain a separate diagnostic rather than silently replacing the
as-is-region result.

## Target package layout

The independent distribution should live with the KIM component:

```text
KIM/python/
├── pyproject.toml
├── README.md
├── src/
│   └── kim/
│       ├── __init__.py
│       ├── cli.py
│       ├── config.py
│       ├── errors.py
│       ├── executable.py
│       ├── namelist.py
│       ├── profiles.py
│       ├── results.py
│       ├── runs.py
│       ├── simulation.py
│       └── sweep.py
└── tests/
    ├── fixtures/
    ├── integration/
    └── unit/
```

The proposed distribution name is `kamel-kim`. Its import package is `kim`, and its console script
is also `kim`. Keeping the distribution name distinct from the import and command reduces the risk
of a package-index name collision while preserving the desired user interface.

The root KAMEL build continues to own `KIM.x`. Installing the Python distribution does not compile
or bundle the executable in the MVP.

## Public Python API

The API is organized around serializable request and result objects rather than process-control
objects with hidden mutable state.

```python
from pathlib import Path

from kim import Simulation, SimulationConfig

config = SimulationConfig.electrostatic_periodic(
    profiles=Path("profiles"),
    m_mode=7,
    n_mode=2,
    btor=-17_977.413,
    major_radius=165.0,
    radial_minimum=3.0,
    plasma_radius=63.0,
    periodic_n_rg=1024,
)

result = Simulation(config).run()
print(result.run_id)
print(result.status)
print(result.output_file)
```

The main public types are:

- `SimulationConfig`: complete validated scientific configuration.
- `ProfileSet`: validated paths and metadata for a set of `r_eff` profiles.
- `Simulation`: prepares and executes one immutable request.
- `RunResult`: run identity, status, logs, manifest, and structured output access.
- `SweepSpec`: immutable base configuration, one validated variation, and failure policy.
- `SweepResult`: ordered child-run results and sweep-level metadata.
- `RunRepository`: list and inspect existing runs without executing anything.

The package root should expose only the stable types and convenience functions. Fortran namelist
group names, subprocess details, and h5py handles remain internal implementation concerns.

## Configuration model

Pydantic models provide runtime validation and JSON Schema generation for the future MCP layer.
Models should reject unknown fields. Each field should carry, where applicable:

- Python and JSON type;
- default or required status;
- units;
- human-readable description;
- numerical range or enumerated values;
- whether it is safe to sweep;
- stability status if a field is advanced or experimental.

Internally, the model maps explicitly to the current namelist groups:

```text
KIM_CONFIG
WKB_DISPERSION
KIM_IO
KIM_SETUP
KIM_GRID
KIM_SPECIES
KIM_PROFILES
KIM_PERIODIC
KIM_FLR2
```

Public field names should use clear domain language. The serializer owns translations such as
`major_radius -> r0` and `radial_minimum -> r_min`. The serialized file must use the exact group
order expected by `read_config.f90`.

Configuration creation has two paths:

```python
config = SimulationConfig.electrostatic_periodic(...)
```

and

```python
config = SimulationConfig.from_namelist("KIM_config.nml")
```

`from_namelist` must parse, normalize, and validate the file. It must reject unsupported keys and
must not preserve arbitrary opaque entries that would bypass the MCP schema. A separately named
low-level inspection function may report unknown legacy entries without making them executable.

Operational values are owned by the runner when a run is staged:

- `profile_location` and `input_profile_dir` point to the copied run inputs;
- `output_path` points to the run result directory;
- `hdf5_output` is true;
- `hdf5_input` is false;
- `h5_out_file` is selected from the run type unless explicitly supported by the model.

## Profile model

The MVP accepts `r_eff` profile files in CGS units:

| Profile | Units | MVP requirement |
|---|---|---|
| `n.dat` | 1/cm^3 | required |
| `Te.dat` | eV | required |
| `Ti.dat` | eV | required |
| `q.dat` | dimensionless | required |
| `Er.dat` | statV/cm | optional to the Fortran core, normally present for automation |
| `Vz.dat` | cm/s | optional |

Validation should check two-column numeric data, finite values, strictly increasing coordinates,
matching grids where the Fortran reader requires them, radial coverage, density unit plausibility,
and the requested resonant-surface crossing. Validation errors should identify the file, row, field,
units, and violated constraint.

Profiles are copied by default. A run never symlinks mutable scientific inputs. The manifest records
the source path and SHA-256 digest for each copied file.

Er scaling is a first-class profile transformation:

```python
from kim import ProfileScale, SweepSpec

sweep = SweepSpec(
    base=config,
    variation=ProfileScale(profile="Er", values=[-1.0, 0.0, 1.0]),
)
```

This is distinct from changing a scalar namelist field. Future transforms can use the same explicit
pattern without allowing arbitrary Python code or shell commands.

## Namelist generation

The Python model is authoritative for API-created requests. f90nml may be used to parse and render
Fortran syntax, but the package must control group ordering and supported keys.

Generated namelists should be deterministic so that equivalent requests produce equivalent input
files and hashes. The generated file should include a short comment declaring the Python package
version and manifest schema version; complete provenance belongs in `manifest.json`.

The renderer should have golden text tests for all three stable run types. Those tests validate the
Python-to-Fortran boundary and should fail clearly when the Fortran configuration contract changes.

## Executable discovery and execution

Executable discovery is deterministic, in this order:

1. An explicit `Path` passed to `Simulation` or the CLI.
2. The `KIM_EXECUTABLE` environment variable.
3. `build/install/bin/KIM.x` relative to a detected KAMEL checkout.
4. `KIM.x` on `PATH`.

The resolver must not search for an executable named `kim`, because that name belongs to the Python
CLI and would create recursion.

Execution uses `subprocess` without a shell. The command is a two-element argument vector containing
the resolved executable and staged namelist. Standard output and standard error go to separate log
files. The manifest is updated atomically before launch and after completion.

## Run repository and reproducibility

The default layout is:

```text
runs/
└── 20260909T142530Z-periodic-er-scan-a1b2c3d4/
    ├── manifest.json
    ├── inputs/
    │   ├── request.json
    │   ├── KIM_config.nml
    │   └── profiles/
    ├── logs/
    │   ├── stdout.log
    │   └── stderr.log
    └── results/
        └── m7_n2/
            └── out_ES_periodic.h5
```

Run states are `prepared`, `running`, `succeeded`, `failed`, `timed_out`, and `interrupted`.
Completed run directories are treated as immutable.

The manifest records:

- manifest schema version and Python package version;
- run ID, optional label, parent sweep ID, and state;
- requested and normalized parameters;
- start and finish timestamps and duration;
- executable path and SHA-256 digest;
- argument vector and working directory;
- process exit status;
- KAMEL git commit and dirty-state flag when discoverable;
- environment values that affect execution, initially `OMP_NUM_THREADS`;
- input source paths, staged paths, and digests;
- expected and discovered output artifacts;
- validation or failure details.

A periodic run succeeds only when the process exits with zero and the readable HDF5 output contains
at least `fields/Phi`, `fields/jpar`, `backs/e/r`, and `setup/periodic_scale/dx_asis`. Each other run
type gets its own explicit artifact predicate.

## Result API

`results.py` owns HDF5 lifetime and conversion. It must not return live h5py datasets after closing
the file. Complex compound datasets with `real` and `imag` members become NumPy complex arrays.

The generic layer supports safe discovery:

```python
result.list_datasets()
result.read_dataset("fields/Phi")
result.dataset_metadata("fields/Phi")
```

Typed periodic conveniences sit above it:

```python
periodic = result.periodic
periodic.radius
periodic.potential
periodic.parallel_current_density
periodic.integrated_parallel_current(region="as_is")
periodic.integrated_parallel_current(region="full_window")
```

The typed API validates compatible shapes and gives explicit errors for partial or older files.
Species-resolved current fields remain optional capabilities discoverable from the result.

## Sweep API

The MVP executes sweeps sequentially. A sweep contains an immutable base configuration and exactly
one typed variation:

```python
from kim import LinearRange, ParameterSweep, run_sweep

spec = ParameterSweep(
    base=config,
    parameter="periodic.n_rg",
    values=LinearRange(start=512, stop=2048, count=4),
)
results = run_sweep(spec)
```

Supported value sources are explicit lists and validated linear ranges. Logarithmic ranges may be
included if their integer and positivity semantics are unambiguous. Sweep parameters must resolve
through the model schema and be marked sweepable.

An Er scan uses `ProfileScale` rather than a scalar parameter path. Every child run stores the
resolved scale and digest of its transformed `Er.dat`. The sweep manifest records child order and
continues after a failed child by default while reporting partial failure at sweep level.

Parallel execution, schedulers, resumption, adaptive sampling, and multidimensional sweeps are later
features.

## CLI

Typer provides the console frontend. Commands call the public Python API and contain no namelist,
filesystem, subprocess, or HDF5 business logic.

```text
kim parameters [--format table|json-schema]
kim validate CONFIG
kim run CONFIG [--profiles DIR] [--executable PATH] [--runs-dir DIR]
kim sweep CONFIG --parameter PATH --values VALUE...
kim sweep CONFIG --scale-profile Er --values VALUE...
kim status RUN_ID
kim inspect RUN_ID [--format table|json]
kim result RUN_ID [--list | --dataset PATH]
```

The CLI should support explicit model fields directly over time, but the MVP may require a complete
JSON request or namelist for large scientific configurations. Any `--set` syntax must resolve only
known schema paths and pass through normal validation.

CLI errors go to stderr, use stable nonzero exit codes, and do not expose Python tracebacks unless
`--debug` is supplied.

## Future MCP boundary

The API should later map directly to a small MCP surface:

```text
get_parameters()
run_simulation(request)
run_sweep(request)
list_runs(filters)
inspect_run(run_id)
read_result(run_id, dataset_or_quantity)
```

Every request and response must be JSON-serializable or refer to a run artifact by stable ID. The MCP
server will call the same API used by the CLI and will not receive a shell-execution primitive, raw
namelist mutation function, arbitrary Python callback, or unrestricted filesystem path operation.

## Existing Python code disposition

| Existing component | Decision | Reason |
|---|---|---|
| `python/KIMpy/KIMpy.py` | Replace, then remove | Uses process-global working-directory changes, shell-created symlinks, weak error handling, and no durable run metadata. |
| `python/KIMpy/KIMData.py` | Replace, then remove | Dataset names have drifted, missing-data behavior is inconsistent, and HDF5 ownership is unsafe. |
| `python/KIMgui/` | Remove after CLI replacement | GUI owns stale configuration and execution logic instead of calling a reusable API. |
| `python/KIMpy/kim_gui.py` | Remove | It is an unfinished duplicate GUI prototype. |
| Python WKB/dispersion helpers | Keep outside the stable package pending a separate audit | They reimplement scientific calculations and are outside the MVP run-orchestration boundary. |
| `KIMPoissonSolver` and `KIMElectromagneticSolver` | Keep outside the stable package | They are validation/re-solving tools rather than the executable automation layer. |
| `python/utility/create_parabolic_profiles.py` | Replace with a deterministic test-fixture generator | Its API and sign-dependent floors are unsuitable for a stable public contract, but its purpose is useful. |
| Other KAMELpy packages | Leave unchanged | They belong to KiLCA and QL-Balance workflows and are outside this focused change. |

Deletion should occur only after tests demonstrate that the new package covers the useful KIM run
and result behavior. No compatibility shim should preserve unsafe global `chdir` or shell behavior.

## MVP acceptance criteria

The MVP is complete when all of the following are true:

- `pip install -e KIM/python` installs `import kim` and the `kim` command.
- The three stable run types have complete validated configuration models.
- A canonical namelist can be generated and parsed back without losing supported values.
- A valid `r_eff` profile set can be validated and copied into a run.
- `Simulation.run()` resolves and launches `KIM.x` without a shell.
- Each run has an atomic manifest, preserved inputs, separate logs, and an explicit final state.
- Partial HDF5 output is reported as failure.
- Periodic HDF5 output can be read through both generic and typed result APIs.
- Er scaling and scalar one-dimensional sweeps execute sequentially with child manifests.
- CLI help, validation, run, sweep, status, inspect, and result commands use the same API.
- Unit tests do not require a compiled solver.
- One opt-in integration test runs a small periodic parabolic-profile case against a built `KIM.x`.
- CI runs the Python unit suite and leaves the existing CTest and golden-record jobs intact.

## MVP implementation steps

### Task 1: Scaffold the independent distribution

**Files:**

- Create: `KIM/python/pyproject.toml`
- Create: `KIM/python/README.md`
- Create: `KIM/python/src/kim/__init__.py`
- Create: `KIM/python/src/kim/cli.py`
- Create: `KIM/python/tests/unit/test_package.py`

**Steps:**

1. Write a test importing `kim` and invoking `kim --help` with Typer's `CliRunner`.
2. Run `python -m pytest KIM/python/tests/unit/test_package.py -v` and verify collection fails.
3. Add the src-layout package, Python 3.10 requirement, dependencies, and `kim = "kim.cli:app"`.
4. Export an explicit package version and create an empty Typer application with descriptive help.
5. Install with `python -m pip install -e KIM/python` in the development environment.
6. Re-run the focused test and verify it passes.
7. Commit with `feat(KIM): scaffold standalone Python package`.

### Task 2: Define validated configuration models

**Files:**

- Create: `KIM/python/src/kim/config.py`
- Create: `KIM/python/src/kim/errors.py`
- Create: `KIM/python/tests/unit/test_config.py`

**Steps:**

1. Write failing tests for required fields, enum values, positive dimensions, radial bounds, mode
   constraints, unsupported run types, unknown fields, and forced HDF5 settings.
2. Add tests that JSON Schema includes units and descriptions for representative fields.
3. Run the focused tests and verify they fail because the models do not exist.
4. Implement strict grouped Pydantic models and the top-level `SimulationConfig` discriminated by
   run type.
5. Add constructors for periodic electrostatic, global electrostatic, and FLR2 configurations.
6. Run the focused tests and verify they pass.
7. Commit with `feat(KIM): add validated simulation configuration`.

### Task 3: Generate and import canonical namelists

**Files:**

- Create: `KIM/python/src/kim/namelist.py`
- Create: `KIM/python/tests/unit/test_namelist.py`
- Create: `KIM/python/tests/fixtures/namelists/electrostatic_periodic.nml`
- Create: `KIM/python/tests/fixtures/namelists/electrostatic.nml`
- Create: `KIM/python/tests/fixtures/namelists/flr2.nml`

**Steps:**

1. Write failing tests for exact group order and representative Python-to-Fortran field mappings.
2. Write failing round-trip tests for each stable run type.
3. Write failing tests for unknown keys, `hdf5_input = .true.`, and obsolete template keys.
4. Implement deterministic model-to-namelist serialization.
5. Implement `SimulationConfig.from_namelist` through f90nml and normal model validation.
6. Compare the periodic fixture with the validated campaign template while excluding operational
   paths that the runner owns.
7. Run the focused tests and verify they pass.
8. Commit with `feat(KIM): add canonical namelist boundary`.

### Task 4: Add deterministic parabolic profiles and profile validation

**Files:**

- Create: `KIM/python/src/kim/profiles.py`
- Create: `KIM/python/tests/fixtures/generate_profiles.py`
- Create: `KIM/python/tests/unit/test_profiles.py`

**Steps:**

1. Write failing tests for two-column data, finite values, monotonic radius, compatible grids,
   profile coverage, CGS density plausibility, and a requested q crossing.
2. Write a failing test generating profiles with known values at the resonant surface.
3. Add signed Er and Vz cases so floor logic cannot reverse or clip their sign incorrectly.
4. Implement the test-fixture generator with explicit units and deterministic formatting.
5. Implement `ProfileSet.validate()` and actionable error messages.
6. Implement copying and SHA-256 recording without symlinks.
7. Run the focused tests and verify they pass.
8. Commit with `feat(KIM): validate and stage radial profiles`.

### Task 5: Resolve the scientific executable

**Files:**

- Create: `KIM/python/src/kim/executable.py`
- Create: `KIM/python/tests/unit/test_executable.py`

**Steps:**

1. Write failing tests for explicit path, environment variable, checkout build, PATH lookup, missing
   executable, non-executable file, and refusal to resolve the CLI command `kim`.
2. Implement the ordered resolver returning an absolute `Path`.
3. Add executable SHA-256 calculation and optional KAMEL git metadata discovery.
4. Run the focused tests and verify they pass.
5. Commit with `feat(KIM): add deterministic executable discovery`.

### Task 6: Build the run repository and manifest state machine

**Files:**

- Create: `KIM/python/src/kim/runs.py`
- Create: `KIM/python/tests/unit/test_runs.py`

**Steps:**

1. Write failing tests for unique run IDs, directory layout, atomic manifest writes, state
   transitions, input digests, listing, and inspection.
2. Implement the versioned manifest models and `RunRepository`.
3. Reject invalid transitions and mutation of completed runs.
4. Simulate interruption between writes and verify the last complete manifest remains readable.
5. Run the focused tests and verify they pass.
6. Commit with `feat(KIM): add reproducible run repository`.

### Task 7: Implement simulation orchestration with a fake executable

**Files:**

- Create: `KIM/python/src/kim/simulation.py`
- Create: `KIM/python/tests/fixtures/fake_kim.py`
- Create: `KIM/python/tests/unit/test_simulation.py`

**Steps:**

1. Write a fake executable that records its argument vector and working directory and can emulate
   success, nonzero exit, timeout, missing output, partial HDF5 output, and stderr output.
2. Write failing tests for staging, command construction, working directory, separate logs, exit
   status, timeout, and final manifest state.
3. Implement `Simulation.prepare()` and `Simulation.run()` using `subprocess` without a shell.
4. Implement run-type-specific output predicates, starting with the strict periodic predicate.
5. Verify failed and timed-out runs retain inputs and diagnostics.
6. Run the focused tests and verify they pass.
7. Commit with `feat(KIM): orchestrate executable runs`.

### Task 8: Read generic and periodic HDF5 results

**Files:**

- Create: `KIM/python/src/kim/results.py`
- Create: `KIM/python/tests/unit/test_results.py`

**Steps:**

1. Generate minimal HDF5 fixtures in tests, including real arrays, compound complex arrays,
   attributes, optional species fields, incompatible shapes, and partial files.
2. Write failing tests for dataset discovery, reads after file closure, metadata, and complex decode.
3. Write failing tests for the periodic typed view and both integration regions.
4. Implement generic result access with bounded, explicit paths.
5. Implement `PeriodicResult` shape checks and the campaign-compatible trapezoidal integral.
6. Run the focused tests and verify they pass.
7. Commit with `feat(KIM): add structured HDF5 results`.

### Task 9: Add sequential scalar and Er-profile sweeps

**Files:**

- Create: `KIM/python/src/kim/sweep.py`
- Create: `KIM/python/tests/unit/test_sweep.py`

**Steps:**

1. Write failing tests for explicit values, linear ranges, schema path resolution, unsweepable
   fields, invalid integer values, deterministic child order, and partial failure.
2. Write failing tests for Er scaling at negative, zero, and positive factors without changing the
   source profile.
3. Implement immutable scalar variations and `ProfileScale`.
4. Implement sequential child execution and the sweep manifest.
5. Run the focused tests and verify they pass.
6. Commit with `feat(KIM): add validated parameter sweeps`.

### Task 10: Implement the thin CLI

**Files:**

- Modify: `KIM/python/src/kim/cli.py`
- Create: `KIM/python/tests/unit/test_cli.py`

**Steps:**

1. Write failing command tests for `parameters`, `validate`, `run`, `sweep`, `status`, `inspect`, and
   `result` using temporary run repositories and the fake executable.
2. Test help output and concise validation errors for every command.
3. Implement commands as adapters over the public Python API.
4. Add JSON output modes suitable for scripts and future tool adapters.
5. Ensure expected user errors do not print tracebacks without `--debug`.
6. Run the focused tests and verify they pass.
7. Commit with `feat(KIM): add Python command-line interface`.

### Task 11: Add the periodic integration reference

**Files:**

- Create: `KIM/python/tests/integration/test_periodic_run.py`
- Create: `KIM/python/tests/fixtures/periodic_reference.json`
- Modify: `KIM/python/README.md`

**Steps:**

1. Generate a compact deterministic `r_eff` profile set in a temporary directory.
2. Choose a small periodic configuration that completes quickly while retaining a known resonant
   surface and nontrivial fields.
3. Run it manually against `build/install/bin/KIM.x` and inspect the full HDF5 schema.
4. Record only robust reference quantities such as shapes, finite-value requirements, derived M,
   q-crossing position, and selected integrated scalars with justified tolerances.
5. Write the opt-in integration test, skipped when `KIM.x` is unavailable.
6. Verify it fails for missing mandatory datasets and passes against the current executable.
7. Document how and why the reference may be updated.
8. Commit with `test(KIM): add periodic Python integration reference`.

### Task 12: Integrate Python tests into CI

**Files:**

- Modify: `.github/workflows/ci.yml`
- Modify: `KIM/python/README.md`

**Steps:**

1. Add a Python job or focused step that installs `KIM/python` with test dependencies.
2. Run the unit suite independently of the Fortran build.
3. Run the periodic integration test only where the KIM build and dependencies are already
   available and its runtime is acceptable.
4. Preserve the existing CTest and golden-record behavior.
5. Run the local unit suite and relevant CTest subset.
6. Commit with `ci(KIM): test standalone Python interface`.

### Task 13: Remove replaced obsolete KIM Python code

**Files:**

- Remove after replacement verification: `python/KIMpy/KIMpy.py`
- Remove after replacement verification: `python/KIMpy/KIMData.py`
- Remove: `python/KIMpy/kim_gui.py`
- Remove after CLI replacement: `python/KIMgui/`
- Modify: `python/pyproject.toml`
- Modify any documentation or imports found by `rg`.

**Steps:**

1. Search the repository again for imports and entry points before deletion.
2. Document any remaining behavior that has no replacement and decide whether it is scientifically
   useful or obsolete.
3. Add missing replacement coverage where justified.
4. Remove the obsolete modules and narrow broad package discovery if necessary.
5. Run both KIM Python tests and affected KAMELpy tests.
6. Commit with `refactor(KIM): remove obsolete Python interfaces`.

## Steps after the MVP

### Stabilize all scientific result contracts

- Define typed electrostatic and FLR2 result views after examining representative current outputs.
- Add dataset schema versions or explicit capability metadata to KIM HDF5 files.
- Add complete git revision, timestamp, executable identity, and normalized configuration metadata
  to the HDF5 output where it benefits non-Python consumers.
- Standardize units and descriptions on datasets and attributes.
- Add numerical compatibility tests across supported compiler and HDF5 versions.

### Expand profile support

- Design explicit `sqrt_psiN -> r_eff` conversion inputs without hidden `$CODE`, AUG wall, or local
  environment dependencies.
- Represent GEQDSK and equilibrium conversion provenance in manifests.
- Add controlled interpolation onto a shared grid with a documented numerical method.
- Support generated and transformed profiles as reusable, versioned request objects.

### Expand sweep execution

- Add logarithmic ranges where constraints are well-defined.
- Add multidimensional Cartesian sweeps with safeguards against accidental run explosions.
- Add resume, retry, cancellation, and bounded local parallelism.
- Add scheduler backends only after the local execution contract is stable.
- Add derived sweep tables without discarding individual HDF5 artifacts.

### Add the MCP server

- Create a separate optional distribution or extra depending on the stable `kim` package.
- Generate tool schemas from Pydantic models rather than maintaining duplicate JSON schemas.
- Expose only validated operations and repository-scoped run identifiers.
- Add response-size limits and explicit selection for large arrays.
- Test the server with deterministic tool calls before connecting a local LLM.

### Add local-LLM integration

- Evaluate small local models on a fixed corpus of KIM requests and expected tool calls.
- Measure parameter-selection accuracy, unit handling, refusal of unsupported fields, and recovery
  from validation errors.
- Keep model output behind the same validation layer used by Python and the CLI.
- Do not grant the model shell access or direct source-editing operations.

### Improve developer and user experience

- Add shell completion and richer parameter documentation generated from the schema.
- Add notebook examples that consume stored run results without rerunning simulations.
- Add migration documentation for publication scripts and legacy KIMpy users.
- Consider a higher-level campaign package only after repeated workflows establish stable concepts.

## Open questions

These questions do not block the first scaffolding tasks, but they must be resolved before the
affected interfaces are declared stable:

1. **Distribution name:** Is `kamel-kim` acceptable for packaging, or is another distribution name
   required for publication or an internal package index? The import and command remain `kim`.
2. **Configuration defaults:** Which scientific values may be true API defaults, and which must
   always be explicit? The publication configuration is a reference case, not automatically a
   universal default.
3. **Reference tolerances:** Which derived quantities from the parabolic periodic test should be
   numerically frozen, and what compiler/platform tolerances are scientifically meaningful?
4. **Electrostatic and FLR2 success predicates:** Which exact HDF5 datasets establish a scientifically
   complete run for these two modes?
5. **Er absence:** Should the public API allow omission of `Er.dat` and expose Fortran's force-balance
   generation, or require Er for all MVP automation requests?
6. **Run repository location:** Should the default be `./runs`, a user configuration directory, or
   an explicit required path? MCP deployments will likely require an explicitly scoped repository.
7. **Timeout defaults:** Should simulations have no default timeout, or should each stable run type
   define a conservative default that users can override?
8. **Legacy external users:** Are any KIMpy or Tk GUI entry points used outside this repository? This
   affects release notes and migration documentation, though obsolete implementations are approved
   for removal.
9. **Publication migration:** Should the Er-scan project's collection and plotting code eventually
   consume `RunRepository` and `PeriodicResult`, or remain an independent publication archive?
10. **HDF5 schema ownership:** Should schema-version metadata be added to the Fortran-produced file
    during the MVP, or should the Python manifest carry the first schema version until the result
    contract matures?

## Implementation constraints

- Preserve unrelated working-tree changes.
- Keep the Fortran core authoritative for all physics.
- Avoid arbitrary shell execution and unrestricted path access.
- Do not duplicate configuration, execution, sweep, or result logic in the CLI.
- Keep unit tests independent of a compiled solver through a fake executable.
- Keep the real periodic test compact and separate from the publication dataset.
- Use the publication Er scan as behavioral evidence, not as a checked-in test-data dependency.
- Update this document when a resolved open question changes a public contract.
