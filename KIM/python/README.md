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

The normal Python suite is self-contained and uses a fake executable:

```bash
python -m pytest
```

The periodic integration reference runs a compact parabolic-profile case against a real `KIM.x`.
It is opt-in so ordinary unit tests do not depend on a compiled solver:

```bash
KIM_RUN_INTEGRATION=1 \
KIM_EXECUTABLE=/absolute/path/to/KIM.x \
python -m pytest tests/integration/test_periodic_run.py -v
```

The reference checks mandatory HDF5 datasets, field shapes, finite nonzero fields, the derived
Fourier mode count, the q-crossing radius, and both campaign-compatible parallel-current
integrals. Integral tolerances allow small compiler and linear-algebra differences while still
detecting scientific drift.

Update `tests/fixtures/periodic_reference.json` only after intentionally changing KIM physics or
numerics. Run the case with one OpenMP thread, inspect the complete HDF5 schema and logs, record the
new values, justify the change in review, and keep tolerances no wider than cross-platform evidence
requires.
