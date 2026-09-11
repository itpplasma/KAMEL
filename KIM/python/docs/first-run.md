# KIM first result

This walkthrough takes a packaged periodic example from installation to a stored HDF5 result and
one figure. It uses the public `kim` commands and the small plotting script in
[`../examples/plot_periodic.py`](../examples/plot_periodic.py). The plot step reads an existing
result; it does not start another solver run.

## Build and install

Build the scientific executable from the repository root, then install the Python package with the
optional test and plotting dependencies:

```sh
make KIM
python -m pip install -e './KIM/python[test,plot]'
kim --help
```

Use `kim doctor` to check package identity and executable selection. Pass the built executable
explicitly while learning the workflow:

```sh
kim doctor --executable build/install/bin/KIM.x
```

`doctor` checks that the selected path is an executable file. It does not launch `KIM.x` or check
its runtime libraries.

## Create and validate a case

Create the packaged periodic case in a new directory and change into it. `init` never replaces an
existing destination:

```sh
kim init ./my-periodic-case --example periodic
cd my-periodic-case
kim validate request.json
```

The request keeps `profiles.directory` as `./profiles`, so running `validate` and `run` from the
case directory is required unless `--profiles` supplies another directory. The validator checks
the profile files, radial coverage, and the signed resonance `q = -m_mode / n_mode`.

The copied example has 65 profile samples and compact numerical settings (`n_rg = 64`) for a quick
demonstration. Check radial, periodic-window, and quadrature convergence before using it for a
scientific result.

## Run and find the HDF5 output

Run the case with one OpenMP thread and request machine-readable output:

```sh
kim run request.json \
  --executable ../build/install/bin/KIM.x \
  --runs-dir runs \
  --omp-threads 1 \
  --format json | tee run.json
```

The report's `run_id`, `run_directory`, and `output_file` fields identify the persisted run. The
HDF5 result is normally below
`runs/<run-id>/results/m<mode>_n<toroidal-mode>/out_ES_periodic.h5`; use the report's `output_file`
value as the source of truth. The run directory also contains staged immutable inputs, the
normalized request, generated namelist, manifest, and `logs/stdout.log` and `logs/stderr.log`.

Inspect the lifecycle record and HDF5 structure without rerunning anything:

```sh
RUN_ID=$(python -c 'import json; print(json.load(open("run.json"))["run_id"])')
kim status "$RUN_ID" --runs-dir runs
kim inspect "$RUN_ID" --runs-dir runs --format json
kim result "$RUN_ID" --runs-dir runs --list
```

For a failed preparation or solver process, inspect `failure` in the `kim run --format json`
report or manifest, then read both log files. A failed run remains a useful diagnostic record and
can be corrected and rerun under a new run ID.

## Make the first figure

Pass the exact HDF5 path from `output_file` and an explicit destination to the plotting example:

```sh
RESULT=$(python -c 'import json; print(json.load(open("run.json"))["output_file"])')
python ../KIM/python/examples/plot_periodic.py "$RESULT" periodic.png
```

The two panels show the real and imaginary parts of `fields/Phi` and `fields/jpar` against the
effective radius. The vertical marker is the resonant radius and the shaded interval is the
periodic as-is region. KIM writes `Phi` in statV and `jpar` in statA/cm²; the radius and periodic
widths are in cm. The script uses `Result.periodic`, closes the HDF5 access after each read, closes
the figure after saving, and never invokes the solver.

## Profiles and physical conventions

The current Python workflow accepts profiles in effective radius `r_eff` with CGS values:

| File | Quantity | Units | Required |
| --- | --- | --- | --- |
| `n.dat` | Electron density | `1/cm^3` | yes |
| `Te.dat` | Electron temperature | `eV` | yes |
| `Ti.dat` | Ion temperature | `eV` | yes |
| `q.dat` | Safety factor | dimensionless | yes |
| `Er.dat` | Radial electric field | `statV/cm` | optional |
| `Vz.dat` | Toroidal velocity | `cm/s` | optional |

The required profile grids must match and have strictly increasing radii. `Er.dat` can use another
grid and is interpolated. If `Er.dat` is absent, KIM calculates the radial electric field from
force balance; `Vz.dat` is also optional. The validator reports missing files or unit and coverage
problems before execution.

The packaged case is a resolution example, not a convergence claim. Its profile domain and window
settings are intentionally small. A production study should vary radial, periodic-window, and
quadrature settings and retain those requests with the resulting HDF5 files.

## Ordered electric-field sweep

To run an ordered sweep of the full radial electric-field profile, keep repeating `--values` in the
desired order:

```sh
kim sweep request.json \
  --scale-profile Er \
  --values -1 --values 0 --values 1 \
  --executable ../build/install/bin/KIM.x \
  --runs-dir sweeps \
  --omp-threads 1 \
  --format json | tee sweep.json
```

Each child receives its own copied request, profiles, logs, HDF5 output, and manifest. The sweep
report lists child run IDs and statuses. Use `kim inspect` or `kim result` with a child ID to find
its stored result.

## Local smoke record

A local smoke trial of these commands used an Apple M4 Pro on macOS arm64, Python 3.12.14,
`kamel-kim` 0.1.0, and a compatible prebuilt `KIM.x`; the single run completed in about 0.31 s with
one OpenMP thread. This timing is descriptive only. The executable build, machine, and numerical
resolution can change the result and runtime. An external user should also run the synthetic plot
case and report any undocumented steps.
