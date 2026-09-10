# KIM request JSON guide

A request JSON file describes one KIM simulation. It contains scientific inputs only: the CLI
selects the executable, run repository, timeout, and OpenMP thread count separately. KIM remains the
scientific solver, while the Python layer validates this request and translates it to
`KIM_config.nml`.

## Start with a periodic request

To create a complete, standalone copy of the packaged periodic case, use the public CLI:

```bash
kim init ./my-periodic-case --example periodic
cd my-periodic-case
kim validate request.json
```

The destination must not already exist. The copied request refers to `./profiles`, so validate and
run it after changing into the case directory. The case uses compact demonstration settings for
the reference workflow; check radial, periodic-window, and quadrature convergence before treating
it as a production setup.

The forced-periodicity electrostatic solver is the recommended starting point. Copy
[`request-periodic.json`](../examples/request-periodic.json), then change the profile directory and
physical parameters for your case:

```json
{
  "run": {
    "run_type": "electrostatic_periodic",
    "periodic": {
      "as_is_width_scale": 2.0,
      "transition_width_scale": 4.0,
      "wavenumber_cutoff_scale": 2.0,
      "n_rg": 64,
      "bparallel_ratio_real": 0.0,
      "bparallel_ratio_imag": 0.0
    }
  },
  "plasma": {
    "source": "builtin",
    "isotope": "D"
  },
  "setup": {
    "btor": -17977.413,
    "major_radius": 165.0,
    "m_mode": 7,
    "n_mode": 2,
    "frequency": 0.0,
    "br_boundary_real": 1.0,
    "br_boundary_imag": 0.0
  },
  "grid": {
    "radial_minimum": 3.0,
    "plasma_radius": 63.0,
    "l_space_dim": 64,
    "rg_space_dim": 64,
    "larmor_skip_factor": 20.0,
    "gauss_nodes_x": 11,
    "gauss_nodes_x_prime": 10,
    "gauss_nodes_theta": 7
  },
  "profiles": {
    "directory": "./profiles",
    "coordinate_type": "r_eff"
  }
}
```

Paths are interpreted relative to the directory where `kim` is invoked. Use `--profiles` to replace
the directory stored in a request without editing the JSON file:

```bash
kim validate KIM/python/examples/request-periodic.json --profiles /path/to/profiles
kim run KIM/python/examples/request-periodic.json \
  --profiles /path/to/profiles \
  --executable /path/to/KIM.x \
  --runs-dir runs
```

`--executable` may be omitted when `KIM_EXECUTABLE` is set, a checkout-local
`build/install/bin/KIM.x` exists, or `KIM.x` is on `PATH`.

## Request sections

| Section | Purpose | Required |
| --- | --- | --- |
| `run` | Selects `electrostatic_periodic`, `electrostatic`, or `flr2` and its controls. | Yes |
| `plasma` | Selects built-in hydrogen/deuterium or explicit ion species. | Yes |
| `setup` | Defines geometry, mode numbers, frequency, and boundary perturbation. | Yes |
| `grid` | Defines the radial domain and numerical resolution. | Yes |
| `profiles` | Locates the named `r_eff` profile files. | Yes |
| `physics` | Controls collision models, species response, and scaling. | No; validated defaults apply. |
| `io` | Controls HDF5 output, logging, and diagnostics. | No; validated defaults apply. |

Unknown fields are rejected. Numeric values must be finite, and constraints spanning several
sections are checked together. For example, `radial_minimum` must be smaller than `plasma_radius`,
both radii must lie inside the profile domain, and the safety-factor profile must contain the
signed `q = -m_mode / n_mode` resonance used by KIM.

List every public parameter, its type, units, description, and whether it can be swept with:

```bash
kim parameters
kim parameters --format json-schema > request.schema.json
```

The generated JSON Schema is the authoritative machine-readable contract.

## Profile directory

Version 1 supports profiles in effective radius `r_eff` only. Every file is whitespace-delimited
with two numeric columns: radius followed by the profile value. Radii must be finite and strictly
increasing. Profile values use CGS units:

| File | Quantity | Units | Required |
| --- | --- | --- | --- |
| `n.dat` | Electron density | `1/cm^3` | Yes |
| `Te.dat` | Electron temperature | `eV` | Yes |
| `Ti.dat` | Ion temperature | `eV` | Yes |
| `q.dat` | Safety factor | dimensionless | Yes |
| `Er.dat` | Radial electric field | `statV/cm` | No |
| `Vz.dat` | Toroidal velocity | `cm/s` | No |

The `n.dat`, `Te.dat`, `Ti.dat`, and `q.dat` grids must match. `Er.dat` may use a different grid;
KIM interpolates it. `Vz.dat` is optional, but when present its grid must match the main grid.
File names can be changed with `density_file`, `electron_temperature_file`,
`ion_temperature_file`, `safety_factor_file`, `radial_electric_field_file`, and
`toroidal_velocity_file` inside `profiles`. These values are plain filenames within `directory`;
absolute paths and directory components are rejected so staging cannot read or write outside the
copied profile directory.

The `physics` section exposes the current I-function controls as
`electron_ifunc_conservation_model` and `ion_ifunc_conservation_model` (`-1`, `0`, `1`, `2`, or
`3`), and the ion temperature-gradient selector as `ion_temperature_gradient_model` (`full`,
`zero_A2`, or `zero_Tprime`). A value of `-1` inherits the legacy `conserve_energy` setting. For a
periodic run, `bparallel_ratio_real` and `bparallel_ratio_imag` prescribe the complex
`B_parallel/Br` drive ratio; both default to zero.

## Other run types

For the global electrostatic solver, use:

```json
"run": {"run_type": "electrostatic"}
```

For the standalone FLR2 solver, use:

```json
"run": {"run_type": "flr2", "terms": {}}
```

FLR2 currently requires `frequency` to be `0.0`, the Fokker-Planck collision model, exactly one ion
species, and at least one species enabled in the potential equation. Complete requests are provided
in [`request-electrostatic.json`](../examples/request-electrostatic.json) and
[`request-flr2.json`](../examples/request-flr2.json).

The examples use a built-in deuterium plasma. An explicit single-ion plasma has this form:

```json
"plasma": {
  "source": "explicit",
  "ions": [{"mass_number": 2, "charge_number": 1}]
}
```

## Defaults and normalized requests

Fields omitted from the examples receive validated defaults. To inspect the complete normalized
request, including all defaults, use Python:

```python
from pathlib import Path
from kim import SimulationConfig

request = SimulationConfig.model_validate_json(Path("request.json").read_text())
print(request.model_dump_json(indent=2))
```

Each prepared run stores this normalized form in `inputs/request.json`, along with the generated
namelist and immutable copies of the profile files.

## Sweeps

Scalar sweeps address public fields by dotted path:

```bash
kim sweep request.json \
  --parameter periodic.n_rg \
  --values 64 --values 96 --values 128
```

The `periodic.` prefix is the public shorthand for `run.periodic.`. `kim parameters` identifies
which fields are sweepable. Scale the entire radial electric-field profile with:

```bash
kim sweep request.json \
  --scale-profile Er \
  --values -10 --values -5 --values 0 --values 5 --values 10
```

Every child run receives its own copied inputs, logs, result file, and manifest.
