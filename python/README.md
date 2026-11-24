# KAMELpy

KAMELpy is the Python interface for preparing, running and reading
QL-Balance and KiLCA runs, plus analysis modules from past studies. KIM has
its own package, `kim` (in [`KIM/python`](../KIM/python/README.md)), which is
installed alongside KAMELpy.

## Installation

Python 3.10 or newer is required. From the KAMEL repository root, create a
virtual environment and install both packages in editable mode:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ./KIM/python -e ./python
```

Editable installs find the executables in `build/install/bin/` of this
checkout, so build the codes first (see the [main README](../README.md)).

The experimental dispersion modules in `KIMpy` and `susc_functions` need
`plasmapy` and `cxroots` in addition:

```sh
python -m pip install -e './python[analysis]'
```

## Packages

| Package | Purpose |
| --- | --- |
| `balance_interface` | `QL_Balance_interface`: write QL-Balance inputs, run `ql-balance.x`, read results |
| [`KiLCA_interface`](KiLCA_interface/README.md) | Prepare, run and post-process KiLCA runs |
| `utility` | Profile helpers, e.g. parabolic test profiles around a resonant surface |
| `device_config` | Device parameters (ASDEX Upgrade, MAST-U) |
| [`KIMpy`](KIMpy/README.md) | Experimental KIM dispersion and field-solver analysis |
| [`pyABC`](pyABC/README.md) | Inputs for the analytical bifurcation criterion (unfinished) |

The remaining packages hold analysis scripts from individual studies.

[Tutorial 4](../tutorials/04-ql-balance/README.md) shows a complete
QL-Balance run set up with KAMELpy.
