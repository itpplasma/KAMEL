# KAMEL — Kinetic plAsma response ModEL

[![CI](https://github.com/itpplasma/KAMEL/actions/workflows/ci.yml/badge.svg)](https://github.com/itpplasma/KAMEL/actions/workflows/ci.yml)
[![Golden record](https://github.com/itpplasma/KAMEL/actions/workflows/golden.yml/badge.svg)](https://github.com/itpplasma/KAMEL/actions/workflows/golden.yml)

KAMEL is a collection of research codes for kinetic plasma response and
quasilinear transport in fusion plasmas. The repository provides one CMake
build for three coupled solvers, their preprocessing tools, and Python
interfaces.

## Components

| Component | Purpose | Documentation |
| --- | --- | --- |
| KiLCA | Cylindrical kinetic response solver with finite-Larmor-radius effects | [KiLCA/README](KiLCA/README) |
| KIM | Non-local integral plasma-response model, including global and forced-periodic electrostatic solvers | [KIM/README.md](KIM/README.md) |
| QL-Balance | One-dimensional quasilinear transport solver using KIM or KiLCA response fields | [QL-Balance/README.md](QL-Balance/README.md) |

Supporting directories include:

- [`PreProc/`](PreProc/) — Fourier-mode and NEO-2 preprocessing tools.
- [`python/`](python/) — KAMELpy interfaces and analysis utilities.
- [`common/`](common/) — shared equilibrium, HDF5, logging, and numerical code.
- [`test/golden/`](test/golden/) — physics-output regression tests for all three solvers.

KIM run types include electrostatic, forced-periodic electrostatic,
electromagnetic, WKB dispersion and second-order local FLR (`flr2`) runs; see
the [KIM namelist reference](KIM/nmls/README.md) for options and restrictions.
Data flow between the codes, file formats and units are described in
[Data.md](Data.md).

## Prerequisites

The complete build requires:

- CMake 3.24 or newer, Ninja, Git, and a POSIX build environment.
- C, C++, and Fortran compilers with OpenMP support.
- BLAS and LAPACK.
- SuperLU (used by KIM).

The first configuration needs network access. CMake fetches pinned or selected
versions of fortnum, fortio, SUNDIALS, SuiteSparse, and libneo. HDF5 and NetCDF
files are read and written through fortio, so the native HDF5 and NetCDF
libraries are not linked. The complex error function library libcerf is
bundled in KIM. Python is optional for the compiled solvers but required for
KAMELpy and several analysis workflows; its package requirements are listed in
[`python/requirements.txt`](python/requirements.txt).

CI builds on Ubuntu with GNU compilers. The CMake configuration also contains
macOS/Homebrew support for LLVM OpenMP and HDF5.

## Build

From the repository root:

```sh
git clone https://github.com/itpplasma/KAMEL.git
cd KAMEL
make all
```

The Makefile is a thin wrapper around the unified CMake/Ninja build. Common
targets are:

```sh
make KIM          # build KIM.x
make KiLCA        # build the KiLCA executables
make QL-Balance   # build ql-balance.x
make PreProc      # build the Fourier preprocessing program
make clean
```

Executables and libraries are written below `build/install/`; the principal
executables are:

```text
build/install/bin/KIM.x
build/install/bin/KiLCA_Normal_...
build/install/bin/ql-balance.x
```

For direct CMake use, the equivalent release build is:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

To select a libneo branch, tag, or commit, pass `LIBNEO_REF=<ref>` to `make`
or `-DLIBNEO_REF=<ref>` to CMake. To use a local libneo checkout, use
`LIBNEO_PATH=<dir>` or `-DLIBNEO_PATH=<dir>`.

## Run

Each solver expects a prepared run directory; the example namelists describe
the configuration format but are not self-contained plasma datasets.

### KIM

KIM accepts the namelist path as its first argument:

```sh
build/install/bin/KIM.x /path/to/KIM_config.nml
```

For the supported Python workflow (electrostatic, forced-periodic electrostatic, and FLR2), install
and use the separate package:

```sh
python -m pip install './KIM/python'
kim --help
kim validate /path/to/request.json
kim run /path/to/request.json --executable /absolute/path/to/KIM.x
```

See the [Python API guide](KIM/python/README.md) and
[request JSON guide](KIM/python/docs/request-json.md). The Python API initially supports `r_eff`
profiles in CGS units and requires HDF5 output. Additional Fortran run modes are documented in the
[namelist reference](KIM/nmls/README.md).

### KiLCA

Run the versioned `KiLCA_Normal_...` executable from a directory containing a
complete KiLCA input set. The
[`python/KiLCA_interface/`](python/KiLCA_interface/) package is the recommended
way to prepare and manage those runs.

### QL-Balance

QL-Balance reads `balance_conf.nml` from its working directory:

```sh
cd /path/to/run-directory
/path/to/KAMEL/build/install/bin/ql-balance.x
```

The template is [`QL-Balance/namelists/balance_conf.nml`](QL-Balance/namelists/balance_conf.nml).

## Test

Build and run the CTest suite with:

```sh
make test
```

The golden-record suite compares the current physics output with the frozen
`golden-baseline` tag. It is intentionally separate because it performs two
complete builds:

```sh
make golden
```

See [`test/golden/README.md`](test/golden/README.md) before adding a case or
re-blessing the baseline.

## Python interface

Install KAMELpy in editable mode from a virtual environment:

```sh
cd python
make init
make install
```

See [`python/README.md`](python/README.md) and the component-specific Python
READMEs for usage.

## Installation

The build stages its executables and libraries below `build/install`. The
`kim` command belongs to the `kamel-kim` Python package; the solver executable
is `KIM.x`. If an older setup left a `kim` alias or symlink pointing to `KIM.x`,
remove it so that the Python command is found (`command -v kim`).

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for the development workflow and
formatting requirements. KAMEL is academic research software; numerical
changes should include focused tests and must pass the physics golden record
when they intentionally affect solver output.

## License

KAMEL is licensed under the MIT License. See [`LICENSE`](LICENSE).
