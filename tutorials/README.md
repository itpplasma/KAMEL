# KAMEL tutorials

Four short walkthroughs that run each code on a small example case. Every
case runs in seconds on a laptop. The cases show how to use the codes; they
are not converged physics results.

| Tutorial | Code | Python needed |
| --- | --- | --- |
| [1. Run KIM from a namelist](01-kim-namelist/README.md) | KIM | Only for plotting |
| [2. Run KIM with the Python interface](02-kim-python/README.md) | KIM | Yes (`kim` package) |
| [3. Run KiLCA](03-kilca/README.md) | KiLCA | Only for plotting |
| [4. Run QL-Balance with KiLCA](04-ql-balance/README.md) | QL-Balance, KiLCA | Yes (KAMELpy) for preparation |

Tutorials 1 and 2 solve the same physics case, once with the executable
directly and once through Python.

## Setup

### Build the codes

From the KAMEL repository root, build all codes as described in the
[main README](../README.md):

```sh
make all
ls build/install/bin
```

### Set up Python

Create a virtual environment and install the two Python packages, the KIM
interface `kim` (from `KIM/python`) and KAMELpy (from `python/`), plus
matplotlib for the plots:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ./KIM/python -e ./python matplotlib
```

Install them in editable mode (`-e`) as shown. The packages then find the
executables in `build/install/bin/` of this checkout without further
configuration.

### Set `KAMEL`

The tutorials refer to the repository root as `$KAMEL` and write their runs
to `$HOME/kamel-tutorials/`. From the repository root:

```sh
export KAMEL=$PWD
```

Then follow the tutorial of your choice. Run the shell commands in a terminal
and the Python scripts with `python` from the directory the tutorial is in at
that point, or paste them into an interactive Python session.

## Testing

`make tutorials` runs every shell command and Python script of the tutorials
in order, in a temporary home directory, and fails if any of them fails. CI
runs it on every pull request, so the tutorials stay in sync with the code.
Keep a block runnable when you edit a tutorial: code blocks marked `sh` or
`python` are executed, blocks marked `text` are not.
