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
| KiLCA | Cylindrical kinetic response solver with finite-Larmor-radius effects | [KiLCA/README.md](KiLCA/README.md) |
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

## Publications

The physics models and numerical methods are described in the following
papers. If you use KAMEL in a publication, cite the papers of the codes you
used.

KiLCA:

- M. F. Heyn, I. B. Ivanov, S. V. Kasilov and W. Kernbichler, *Kinetic
  modelling of the interaction of rotating magnetic fields with a radially
  inhomogeneous plasma*, Nucl. Fusion **46**, S159 (2006).
  [doi:10.1088/0029-5515/46/4/S07](https://doi.org/10.1088/0029-5515/46/4/S07)
- M. F. Heyn, I. B. Ivanov, S. V. Kasilov, W. Kernbichler, I. Joseph,
  R. A. Moyer and A. M. Runov, *Kinetic estimate of the shielding of resonant
  magnetic field perturbations by the plasma in DIII-D*, Nucl. Fusion **48**,
  024005 (2008).
  [doi:10.1088/0029-5515/48/2/024005](https://doi.org/10.1088/0029-5515/48/2/024005)
- I. B. Ivanov, M. F. Heyn, S. V. Kasilov and W. Kernbichler, *Kinetic linear
  model of the interaction of helical magnetic perturbations with cylindrical
  plasmas*, Phys. Plasmas **18**, 022501 (2011).
  [doi:10.1063/1.3551740](https://doi.org/10.1063/1.3551740)

QL-Balance:

- M. F. Heyn et al., *Quasilinear modelling of RMP interaction with a tokamak
  plasma: application to ASDEX Upgrade ELM mitigation experiments*, Nucl.
  Fusion **54**, 064005 (2014).
  [doi:10.1088/0029-5515/54/6/064005](https://doi.org/10.1088/0029-5515/54/6/064005)
- M. Markl et al., *Kinetic study of the bifurcation of resonant magnetic
  perturbations for edge localized mode suppression in ASDEX Upgrade*, Nucl.
  Fusion **63**, 126007 (2023).
  [doi:10.1088/1741-4326/acf20c](https://doi.org/10.1088/1741-4326/acf20c)

KIM:

- M. Markl, S. V. Kasilov and C. G. Albert, *Integral ion Larmor radius
  treatment in the linear plasma response to resonant magnetic
  perturbations*, J. Plasma Phys. (2026, accepted).

## Prerequisites

The complete build requires:

- CMake 3.24 or newer, Ninja, Git, and a POSIX build environment.
- C, C++, and Fortran compilers with OpenMP support.
- BLAS and LAPACK.
- SuperLU (used by KIM).

On Ubuntu 24.04 or newer (Ubuntu 22.04 ships CMake 3.22, which is too old):

```sh
sudo apt-get install git make cmake ninja-build gcc g++ gfortran \
    libopenblas-dev liblapack-dev libsuperlu-dev python3-venv
```

On macOS, install the Xcode command line tools (`xcode-select --install`) and
then, with [Homebrew](https://brew.sh):

```sh
brew install cmake ninja gcc llvm libomp openblas lapack superlu
```

The build picks up Homebrew's LLVM `clang` and GNU Fortran automatically. CI
installs exactly these packages, so the lists stay tested.

The first configuration needs network access. CMake fetches pinned or selected
versions of fortnum, fortio, SUNDIALS, SuiteSparse, and libneo. HDF5 and NetCDF
files are read and written through fortio, so the native HDF5 and NetCDF
libraries are not linked. The complex error function library libcerf is
bundled in KIM. Python is optional for the compiled solvers but required for
the Python interfaces (see [Python interface](#python-interface)).

CI builds and tests on Ubuntu (GNU compilers) and on macOS with Homebrew
(GNU Fortran, LLVM OpenMP).

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

Dependency commits are pinned in `cmake/DependencyPins.cmake` and added
through [ItpDeps](https://github.com/itpplasma/cmake-deps). To build another
branch, tag, or commit of a dependency, pass `-DITP_<NAME>_REF=<ref>` to CMake,
e.g. `-DITP_LIBNEO_REF=<ref>`. To use a local checkout, pass
`-DFETCHCONTENT_SOURCE_DIR_<NAME>=<dir>`. With `make`, put these flags in
`CMAKE_ARGS`, e.g. `make all CMAKE_ARGS=-DITP_LIBNEO_REF=<ref>`. They only take
effect on a fresh configure, so run `make clean` first if `build/` exists.

## Run

The [tutorials](tutorials/README.md) walk through complete, tested example runs
of KIM (with and without Python), KiLCA and QL-Balance. Each runs in seconds.

Each solver expects a prepared run directory; the example namelists in
`KIM/nmls/` and `QL-Balance/namelists/` describe the configuration format but
are not self-contained plasma datasets.

### KIM

KIM accepts the namelist path as its first argument:

```sh
build/install/bin/KIM.x /path/to/KIM_config.nml
```

For the supported Python workflow (electrostatic, forced-periodic electrostatic, and FLR2), use
the `kim` package (installation: [Python interface](#python-interface)):

```sh
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

`make pytest` runs the Python checks for the golden-record harness, the
periodic workflow configuration, and the consistency of the KIM namelist
reference with the source.

## Python interface

KAMEL has two Python packages: `kim` for KIM (in `KIM/python`) and KAMELpy
for QL-Balance and KiLCA (in `python/`). Both need Python 3.10 or newer.
Install them in editable mode into a virtual environment from the repository
root:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ./KIM/python -e ./python
```

Editable installs find the executables in `build/install/bin/` of this
checkout. See [`python/README.md`](python/README.md) and the
[KIM Python API guide](KIM/python/README.md) for usage.

## Installation

The build stages its executables and libraries below `build/install`. The
`kim` command belongs to the `kamel-kim` Python package; the solver executable
is `KIM.x`. If an older setup left a `kim` alias or symlink pointing to `KIM.x`,
remove it so that the Python command is found (`command -v kim`).

## Developer documentation

Design notes, implementation plans, benchmark records and verification
reports are collected in [`docs/`](docs/README.md).

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for the development workflow and
formatting requirements. KAMEL is academic research software; numerical
changes should include focused tests and must pass the physics golden record
when they intentionally affect solver output.

## License

KAMEL is licensed under the MIT License. See [`LICENSE`](LICENSE).
