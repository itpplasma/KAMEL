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
- **CMake 3.24+** and **Ninja** build system
- **Fortran compiler** (gfortran 10+ or ifort)
- **C/C++ compiler** (gcc/g++ 10+ or clang/clang++)
- **MPI** (MPICH or OpenMPI)
- **HDF5** with Fortran bindings
- **LAPACK/BLAS** - Linear algebra
- **Python 3.8+** with pip

### Optional System Dependencies
- **SuperLU** (sparse matrix solver, used by KIM if found)
- **Doxygen** (for documentation generation)

### Auto-Fetched Libraries
The following are automatically downloaded and built if not found on the system:
- **fortnum** - Numerical core (special functions, quadrature, ODE, root finding)
- **SuiteSparse** - Sparse matrix operations (UMFPACK)
- **SUNDIALS** - Numerical differential equation solvers
- **NetCDF** - Network Common Data Form (with Fortran bindings)
- **LAPACK/BLAS** - Fallback if not found on system

### Bundled Libraries (built from source)
- **fortnum_amos_compat** - AMOS complex-Bessel ABI backed by fortnum (`common/math/`)
- **libcerf** - Complex error function (`KIM/src/math/`)

### Python Dependencies
- numpy, scipy, h5py, f90nml, matplotlib

## Compilation

```bash
# Build all three codes (Release mode by default)
make all

# Build in Debug mode
CONFIG=Debug make all

# Build individual components
make KiLCA
make KIM
make QL-Balance

# Clean build
make clean
```

**Note:** External dependencies (LAPACK, SuiteSparse, fortnum, SUNDIALS) are automatically downloaded and built during the first compilation if not found on the system.

To pin libneo to a specific branch, tag, or commit, pass `-DLIBNEO_REF=<ref>` to cmake or `LIBNEO_REF=<ref>` to make. To use a local checkout instead of fetching, pass `-DLIBNEO_PATH=<dir>` / `LIBNEO_PATH=<dir>`.

### System-wide Installation

To install KIM so it can be run from anywhere as `kim`:

```bash
# Step 1: Build and prepare KIM
make install-kim

# Step 2: Follow the instructions displayed, which will show:
sudo ln -sf /path/to/KAMEL/build/install/bin/KIM.x /usr/local/bin/kim
```

Alternatively, you can add an alias to your shell configuration file (`.bashrc`, `.zshrc`, etc.):
```bash
alias kim='/path/to/KAMEL/build/install/bin/KIM.x'
```

### Tested Configurations
- **Apple Silicon**: clang 16.0 + gfortran 14.2
- **Debian/Ubuntu**: GNU compiler 12.2.0

## Quick Start

```bash
# 1. Clone repository
git clone https://github.com/itpplasma/KAMEL.git
cd KAMEL

# 2. Build all codes
make all

# 3. Install Python interface
cd python && make init && make install

# 4. Run tests
make test
```

## Project Structure

- `/KiLCA/` - Finite Larmor radius plasma response solver
- `/KIM/` - Integral formalism plasma response solver
- `/QL-Balance/` - Quasilinear transport code (supports KiLCA and KIM wave codes)
- `/PreProc/` - Preprocessing utilities (fouriermodes, neo-2 templates)
- `/python/` - Python interface (KAMELpy) for all codes
- `/common/` - Shared utilities: equilibrium handling, math libraries, logger

## Contributing

Contributions are welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
