# KIM Python API

`kamel-kim` provides the supported Python API and command-line interface for the KIM
plasma-response solver. The scientific implementation remains in the `KIM.x` Fortran executable;
this package will validate inputs, stage reproducible runs, invoke the executable, and read its
results.

The package currently contains only the initial API and CLI scaffold. Simulation configuration and
execution will be added in the next implementation tasks.

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
