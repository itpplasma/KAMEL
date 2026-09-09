# KIM Python API

`kamel-kim` provides the supported Python API and command-line interface for the KIM
plasma-response solver. The scientific implementation remains in the `KIM.x` Fortran executable;
this package validates inputs, stages reproducible runs, invokes the executable, reads structured
results, and runs one-dimensional parameter scans.

## Requirements

- Python 3.10 or newer
- A separately built `KIM.x` executable for simulation runs

## Development installation

From the KAMEL repository root:

```bash
python -m pip install -e './KIM/python[test]'
```

Check the installed command:

```bash
kim --help
```

## CLI overview

The CLI accepts a complete JSON request or a supported KIM namelist:

```bash
kim validate request.json
kim run request.json --executable /path/to/KIM.x --runs-dir runs
kim sweep request.json --parameter periodic.n_rg --values 512 --values 1024
kim sweep request.json --scale-profile Er --values -1 --values 0 --values 1
kim status RUN_ID --runs-dir runs
kim inspect RUN_ID --runs-dir runs --format json
kim result RUN_ID --runs-dir runs --list
```

Use `kim parameters --format json-schema` to obtain the structured configuration schema used by
the API and future automation tools.

## Tests

The unit suite is self-contained and uses a fake executable. It does not require a Fortran build:

```bash
python -m pytest tests/unit
```

CI runs this suite in a separate Python 3.10 job so the package's minimum supported Python version
and its orchestration behavior are checked independently of the Fortran toolchain.

The periodic integration reference runs a compact parabolic-profile case against a real `KIM.x`.
It is opt-in so ordinary unit tests do not depend on a compiled solver:

```bash
KIM_RUN_INTEGRATION=1 \
KIM_EXECUTABLE=/absolute/path/to/KIM.x \
python -m pytest tests/integration/test_periodic_run.py -v
```

CI enables this test in the existing Fortran build job, using the executable produced at
`build/install/bin/KIM.x`. Running `python -m pytest` locally executes the unit suite and skips the
real-executable case unless both integration environment variables are set.

The reference checks mandatory HDF5 datasets, field shapes, finite nonzero fields, the derived
Fourier mode count, the q-crossing radius, and both campaign-compatible parallel-current
integrals. Integral tolerances allow small compiler and linear-algebra differences while still
detecting scientific drift.

Update `tests/fixtures/periodic_reference.json` only after intentionally changing KIM physics or
numerics. Run the case with one OpenMP thread, inspect the complete HDF5 schema and logs, record the
new values, justify the change in review, and keep tolerances no wider than cross-platform evidence
requires.

## Legacy interface migration

This package replaces the former `KIMpy` runner, `KIMData` reader, and Tk-based KIM GUIs. Use
`Simulation` or `kim run` instead of creating executable symlinks and changing the process working
directory. Supply profiles explicitly; the API copies them into every run directory rather than
linking mutable external directories. Use `Result`, `kim inspect`, and `kim result` for production
HDF5 output instead of converting legacy text field files.

The experimental dispersion-relation, Poisson, and Poisson-Ampere analysis modules remain in the
top-level KAMELpy distribution because they perform separate scientific analysis and are not
replaced by this orchestration API.
