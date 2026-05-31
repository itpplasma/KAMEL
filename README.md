# KAMEL - Kinetic plAsma response ModEL

[![CI](https://github.com/itpplasma/KAMEL/actions/workflows/ci.yml/badge.svg)](https://github.com/itpplasma/KAMEL/actions/workflows/ci.yml)

KAMEL is a scientific computing framework for modeling plasma response to external magnetic perturbations in fusion plasmas.

## Overview

The framework consists of three main codes:

1. **KiLCA** - Kinetic Linear Cylindrical Approximation, plasma response solver using finite Larmor radius formalism
2. **KIM** - KiLCA Integral Model using integral formalism for non-local plasma response
3. **QL-Balance** - Quasilinear 1D radial transport code including anomalous and electromagnetic diffusion for time-evolution studies

## Requirements

### System Dependencies
- **CMake 3.24+**
- **Fortran compiler** (gfortran 10+)
- **C/C++ compiler** with OpenMP support (LLVM clang or gcc/g++)
- **Python 3.8+** with pip
- **NetCDF** with Fortran bindings (`nf-config` must be on PATH)

### Optional System Dependencies
- **Ninja** build system (recommended for faster builds; falls back to Make)
- **OpenBLAS** or **LAPACK/BLAS** (built from source if not found)
- **GSL** (built from source if not found)

### Auto-Fetched Libraries
The following are automatically downloaded and built if not found on the system:
- **HDF5** - Hierarchical Data Format (with Fortran and HL bindings)
- **GSL** - GNU Scientific Library
- **SuiteSparse** - Sparse matrix operations (UMFPACK)
- **SUNDIALS** - Numerical differential equation solvers
- **Zeal** - Mathematical special functions (zero finding)
- **LAPACK/BLAS** - Fallback if not found on system

### Bundled Libraries (built from source)
- **slatec** - Special functions and ODE solvers (`common/math/`)
- **libcerf** - Complex error function (`KIM/src/math/`)

### Python Dependencies
- numpy, scipy, h5py, f90nml, matplotlib

## Building

### Quick Start

```bash
# 1. Clone repository
git clone https://github.com/itpplasma/KAMEL.git
cd KAMEL

# 2. Build all codes (Release mode by default)
make all

# 3. Install Python interface
cd python && make init && make install

# 4. Run tests
make test
```

### Build Commands

```bash
make all            # Build all codes (Release)
CONFIG=Debug make   # Build in Debug mode
make KiLCA          # Build individual components
make KIM
make QL-Balance
make clean          # Clean build
```

**Note:** External dependencies are automatically downloaded and built during the first compilation if not found on the system.

### macOS

On macOS, the build system automatically detects compilers and dependencies:

- **C/C++**: Prefers Homebrew LLVM clang (for OpenMP support). Falls back to any clang on `$PATH`. Apple Clang is not supported (no OpenMP).
- **Fortran**: Uses `gfortran` from `$PATH`.
- **OpenMP**: Located via `find_package(OpenMP)`. Falls back to Homebrew `libomp` if the linker can't find it.
- **BLAS/LAPACK**: Prefers OpenBLAS. On macOS with Homebrew, install with `brew install openblas`.

All compiler and dependency paths can be overridden via CMake variables:
```bash
cmake -S . -B build \
  -DCMAKE_C_COMPILER=/path/to/clang \
  -DCMAKE_CXX_COMPILER=/path/to/clang++ \
  -DCMAKE_Fortran_COMPILER=/path/to/gfortran \
  -DLIBOMP_PREFIX=/path/to/libomp \
  -DHDF5_DIR=/path/to/hdf5/lib/cmake
```

### Nix

A `flake.nix` is provided for building with Nix. Enter the development shell:

```bash
nix develop   # or: nix develop --extra-experimental-features 'nix-command flakes'
cmake -S . -B build -G Ninja -DHDF5_DIR=$HDF5_DIR -DPython_EXECUTABLE=$(which python3)
cmake --build build
```

The Nix shell provides all compilers and dependencies. The build system detects the Nix environment automatically and uses compilers from `$PATH` instead of searching Homebrew paths.

> **Known issue (macOS):** C++ compilation may fail due to a Nix `stdenv` header ordering conflict with `libcxx`. See [#127](https://github.com/itpplasma/KAMEL/issues/127).

### Linux

On Linux, the build uses `gcc`/`g++`/`gfortran` by default. Install dependencies via your package manager:

```bash
# Debian/Ubuntu
sudo apt install cmake gfortran gcc g++ libopenblas-dev libnetcdf-dev libnetcdff-dev python3-numpy

# Fedora
sudo dnf install cmake gcc-gfortran gcc gcc-c++ openblas-devel netcdf-fortran-devel python3-numpy
```

### System-wide Installation

To install KIM so it can be run from anywhere as `kim`:

```bash
make install-kim
# Follow the instructions displayed:
sudo ln -sf /path/to/KAMEL/build/install/bin/KIM.x /usr/local/bin/kim
```

### Tested Configurations
- **Apple Silicon (macOS)**: Homebrew LLVM clang 21 + gfortran 15
- **Debian/Ubuntu**: GNU compiler 12.2.0

## Project Structure

- `/KiLCA/` - Finite Larmor radius plasma response solver
- `/KIM/` - Integral formalism plasma response solver
- `/QL-Balance/` - Quasilinear transport code (supports KiLCA and KIM wave codes)
- `/PreProc/` - Preprocessing utilities (fouriermodes, neo-2 templates)
- `/python/` - Python interface (KAMELpy) for all codes
- `/common/` - Shared utilities: equilibrium handling, math libraries, logger

## References

P. Kravanja, M. Van Barel, O. Ragos, M.N. Vrahatis, F.A. Zafiropoulos,
*ZEAL: A mathematical software package for computing zeros of analytic functions*,
Computer Physics Communications **124** (2000) 212-232.
[doi:10.1016/S0010-4655(99)00429-4](https://doi.org/10.1016/S0010-4655(99)00429-4)

## Contributing

Contributions are welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
