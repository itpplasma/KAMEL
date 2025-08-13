# KAMEL - Kinetic plAsma response ModEL

KAMEL is a scientific computing framework for modeling plasma response to external magnetic perturbations in fusion plasmas.

## Overview

The framework consists of three main codes:

1. **KiLCA** - Kinetic Linear Cylindrical Approximation, plasma response solver using finite Larmor radius formalism
2. **KIM** - KiLCA Integral Model using integral formalism for non-local plasma response
3. **QL-Balance** - Quasilinear 1D radial transport code including anomalous and electromagnetic diffusion for time-evolution studies

## Requirements

### System Dependencies
- **MPI** (MPICH or OpenMPI)
- **HDF5** with parallel support
- **Python 3.8+** with pip
- **CMake 3.16+** and **Ninja** build system
- **Fortran compiler** (gfortran 10+ or ifort)
- **C/C++ compiler** (gcc/g++ 10+ or clang/clang++)

### External Libraries
The following are automatically fetched during compilation:
- **LAPACK/BLAS** - Linear algebra operations
- **SuiteSparse** - Sparse matrix operations
- **GSL** - GNU Scientific Library
- **SUNDIALS** - Numerical differential equation solvers
- **FFTW3** - Fast Fourier transforms

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

**Note:** External dependencies (LAPACK, SuiteSparse, GSL, SUNDIALS) are automatically downloaded and built during the first compilation if not found on the system.

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
- `/QL-Balance/` - Quasilinear transport code (requires KiLCA)
- `/PreProc/` - Preprocessing utilities (fouriermodes, neo-2 templates)
- `/python/` - Python interface (KAMELpy) for all codes
- `/matlab/` - MATLAB interfaces and workflow management
- `/template_scripts/` - Standard workflow templates
- `/external/` - External dependencies (auto-fetched)
- `/Documentation/` - Mathematical background and user guides

## License

See LICENSE file for details.
