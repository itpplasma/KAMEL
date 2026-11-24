# Tutorial 2: Run KIM with the Python interface

The `kim` Python package validates a KIM request, writes the namelist for
you, copies the inputs into a new run directory, runs `KIM.x` and records a
manifest with checksums of every input. This tutorial uses both its
command-line tool and its Python API.

The physics case is the one from [Tutorial 1](../01-kim-namelist/README.md),
solved with the forced-periodicity electrostatic solver
(`electrostatic_periodic`). Each run takes less than a second. The case is a
demonstration, not a converged result.

Before you start, build KAMEL and install the `kim` package as described in
[the tutorial overview](../README.md).

## 1. Prepare a working directory

A request is a JSON file. Copy the example request and the profiles:

```sh
mkdir -p "$HOME/kamel-tutorials/kim-python"
cd "$HOME/kamel-tutorials/kim-python"
cp "$KAMEL/tutorials/02-kim-python/request.json" .
cp -r "$KAMEL/tutorials/data/kim-profiles" profiles
cat request.json
```

The request has five sections: `run` selects the solver, `plasma` the ion
species, `setup` the magnetic field, major radius and mode numbers, `grid` the
radial domain, and `profiles` the profile directory. Paths are relative to the
directory where you call `kim`. Settings you leave out get validated defaults.
The [request guide](../../KIM/python/docs/request-json.md) lists every field.

## 2. Validate the request

`kim validate` checks the request and the profiles without running KIM. It
reports the radial range of the profiles and the radius of the resonant
surface:

```sh
cd "$HOME/kamel-tutorials/kim-python"
kim validate request.json
```

Errors name the field and the problem, for example a grid that extends
beyond the profiles.

## 3. Run KIM

```sh
cd "$HOME/kamel-tutorials/kim-python"
kim run request.json --runs-dir runs
```

`kim` looks for `KIM.x` in the `KIM_EXECUTABLE` environment variable, then in
`build/install/bin/` of the KAMEL checkout the package was installed from,
then on `PATH`. Pass `--executable /path/to/KIM.x` to choose one explicitly.

Each run gets its own directory below `runs/`, named by time, run type and a
short hash:

```text
runs/<run-id>/
├── inputs/       request.json, generated KIM_config.nml, copied profiles
├── logs/         stdout.log, stderr.log
├── results/      KIM output (HDF5)
└── manifest.json status, timing, executable, git commit, input checksums
```

## 4. Inspect the run

```sh
cd "$HOME/kamel-tutorials/kim-python"
RUN_ID=$(ls runs | tail -n 1)
kim status "$RUN_ID" --runs-dir runs
kim result "$RUN_ID" --runs-dir runs --list | grep '^fields/'
```

`kim inspect "$RUN_ID" --runs-dir runs` prints the complete manifest, and
`kim result "$RUN_ID" --runs-dir runs --dataset fields/Phi` prints one
dataset.

## 5. Run and read results from Python

The same workflow is available as a Python API. This script loads the
request, runs it, and plots the electrostatic potential around the resonant
surface:

```python
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from kim import Simulation, SimulationConfig

config = SimulationConfig.model_validate_json(Path("request.json").read_text())
run = Simulation(config, runs_directory="runs").run()
print("status:", run.status.value)
print("output:", run.output_file)

periodic = run.output.periodic
print(f"resonance radius: {periodic.resonance_radius:.2f} cm")
print("integrated parallel current:", periodic.integrated_parallel_current("as_is"))

fig, ax = plt.subplots()
ax.plot(periodic.radius, np.abs(periodic.potential))
ax.axvline(periodic.resonance_radius, linestyle="--", color="gray", label="resonance")
ax.set_xlabel("r [cm]")
ax.set_ylabel("|Phi| [statV]")
ax.legend()
fig.savefig("kim_periodic_phi.png", dpi=120)
print("saved kim_periodic_phi.png")
```

The periodic solver works in a window around the resonant surface, so the
radius covers only part of the grid domain. Use `kim.Result(path)` to open
any KIM output file, including ones from earlier runs.

## 6. Check the resolution with a sweep

A sweep runs the same request for several values of one parameter. Run the
periodic solver with three different numbers of radial samples:

```sh
cd "$HOME/kamel-tutorials/kim-python"
kim sweep request.json --parameter periodic.n_rg \
    --values 128 --values 256 --values 512 --runs-dir sweeps
```

Then compare the potentials:

```python
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from kim import Result

fig, ax = plt.subplots()
for output in sorted(Path("sweeps").glob("*/results/m-6_n2/*.h5")):
    periodic = Result(output).periodic
    ax.plot(periodic.radius, np.abs(periodic.potential), label=f"n_rg = {len(periodic.radius)}")
ax.set_xlabel("r [cm]")
ax.set_ylabel("|Phi| [statV]")
ax.legend()
fig.savefig("kim_n_rg_sweep.png", dpi=120)
print("saved kim_n_rg_sweep.png")
```

`kim parameters` lists every parameter that can be swept. `kim sweep
request.json --scale-profile Er --values -1 --values 0 --values 1` scales the
radial electric field profile instead.

## Next steps

- Change `run.run_type` to `electrostatic` to use the global solver of
  Tutorial 1. The request guide also covers the `flr2` run type.
- Replace the profiles with your own. Version 1 of the interface accepts
  profiles in effective radius `r_eff` and CGS units only.
