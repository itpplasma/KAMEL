# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

KAMEL (Kinetic plAsma response ModEL) is a scientific computing framework for modeling plasma response to external magnetic perturbations in fusion plasmas. It consists of three main codes:

1. **KiLCA** - Cylindrical linear plasma response solver using finite Larmor radius formalism
2. **KIM** - KiLCA Integral Model using integral formalism  
3. **QL-Balance** - Quasilinear 1D radial transport code

## Build Commands

### Building the Project
```bash
# Build all three codes (Release mode by default)
make all

# Build in Debug mode
CONFIG=Debug make all

# Build individual components
make KiLCA
make KIM
make QL-Balance

# Clean build (removes build/)
make clean
```

### Running Tests
```bash
# Run all tests via CTest from root
make test  # invokes ctest --test-dir build
```

### Python Package
```bash
# Install Python interface
cd python && make init && make install
```

## Architecture Overview

### Core Components
- **Fortran/C++ cores** - High-performance computation engines in `/KiLCA/`, `/KIM/`, `/QL-Balance/` (executables use `.x` suffix)
- **Python interfaces** - Modern object-oriented wrappers in `/python/` (KAMELpy)
- **MATLAB interfaces** - Research workflow management in `/matlab/`
- **Template scripts** - Standardized workflows in `/template_scripts/`
- **Common utilities** - Shared math/utils in `/common/`

### Key Directories
- `/external/` - External dependencies (SuiteSparse, GSL, LAPACK, etc.)
- `/PreProc/` - Preprocessing utilities (fouriermodes, neo-2 templates)
- `/Documentation/` - LaTeX documentation and mathematical background
- `/utility_scripts/` - Helper scripts in MATLAB and Python
- `/build/` - Build artifacts and compiled binaries (gitignored)

### Data Flow & Workflow
All data exchange uses **HDF5 format** for standardization. The typical workflow:
1. **Prerun**: Generate HDF5 input files using template scripts
2. **Main run**: Execute solver (KiLCA/KIM/QL-Balance) with appropriate template:
   - `linearrun/` - Quasilinear diffusion coefficients
   - `timeevol/` - Dynamic transport evolution
   - `parameterscan/` - Parameter space exploration
3. **Post-processing**: Python/MATLAB analysis and visualization
4. **Metadata**: HDF5 outputs include git version and timestamps for reproducibility

## Python Interface (KAMELpy)

### Core Classes
- **`KIMpy`** - KIM calculations with dispersion relations and collision models
- **`KiLCA_interface`** - Comprehensive KiLCA workflow management with modular components
- **`QL_Balance_interface`** - Complete transport calculations with automatic preprocessing

### Common Pattern
```python
# Object-oriented configuration
interface = KiLCA_interface(shot, time, path, run_type, machine)
interface.set_modes(m_modes, n_modes)
interface.prepare_balance_input(input_file)
interface.run_balance()
```

## MATLAB Interface

### Core Framework
- **`Balance` class** - Comprehensive workflow management with modular preprocessing
- **Blueprint system** - Template-based configuration in `/matlab/balance/KiLCA_interface/blueprints/`
- **Device-specific configs** - AUG, MAST-U support

### Common Pattern
```matlab
bal = Balance(runpath, shot, time, name, hdf5file);
bal.setModes(m, n);
bal.setCoil(cfile, pfile);
bal.setEqui(gfile, fluxdatapath);
bal.setProfiles(neprof, Teprof, Tiprof, vtprof);
bal.setKiLCA(ion_mass);
bal.run();
```

## Configuration Management

### Namelist Files
- **KIM**: `/KIM/nmls/KIM_config.nml`
- **QL-Balance**: `balance_conf.nml` (in run directory)

### Build Configuration
- **Default build type**: Release (configure with `CONFIG=Debug` for debug builds)
- **Build generator**: Ninja (via CMake, ensure Ninja is available on PATH)
- **Build files**: `build/build.ninja`
- **Compilers**: MPI Fortran (`mpif90`), C/C++ with clang-format support
- **Platforms tested**: Apple Silicon (clang 16.0 + gfortran 14.2), Debian (GNU 12.2.0)

## Template-Based Development

### Standard Templates
- **`script_prerun.m`** - Creates input HDF5 files for all run types
- **`linearrun/`** - Quasilinear diffusion coefficient calculations
- **`timeevol/`** - Dynamic transport evolution
- **`parameterscan/`** - Systematic parameter space exploration

### Development Workflow
1. Use `create_proj_dir.py` for directory structure
2. Execute template prerun script for data preparation
3. Use appropriate template for main calculations
4. Post-process with Python/MATLAB visualization

## Key Dependencies

### External Libraries
- **MPI** - Parallel computing
- **LAPACK/BLAS** - Linear algebra
- **SuiteSparse** - Sparse matrix operations
- **HDF5** - Data storage
- **GSL** - GNU Scientific Library
- **SUNDIALS** - Numerical solvers

### Python Dependencies
- **Core**: numpy, scipy, h5py, f90nml
- **Visualization**: matplotlib

## Development Guidelines

### Coding Style
- **Indentation**: 4 spaces (no tabs)
- **Line width**: 100 columns maximum
- **Fortran formatting**: Use `fprettify` (see `.fprettify` config)
  - Module names use `_m` suffix
  - Test files named `test_*.f90`
  - Executables use `.x` suffix
- **C/C++ formatting**: Use `clang-format` (LLVM base style, see `.clang-format`)

### Testing
- **Framework**: CTest enabled at top level
- **Execution**: Tests run from `build/` via `ctest`
- **New tests**: Add via CMake in relevant subproject
  - Name sources `test_*.f90`
  - Register with `add_test` in `CMakeLists.txt`
  - Keep self-contained, avoid interactive I/O
  - Print diagnostics for debugging

### Version Control
- **Commit style**: Conventional Commits (`feat:`, `fix:`, `refactor:`, `chore:`)
  - Optional scope: `feat(KIM): add new solver`
  - Single-topic, buildable commits
- **Pull Requests**:
  - Include purpose, key changes, usage notes
  - Link related issues
  - Provide test output (`ctest` summary)
  - Require green CI before merge
  - Avoid force-push after review (except for rebase)

### Build Tips
- **Reconfigure**: When in doubt, `make clean` then rebuild
- **Dependencies**: QL-Balance requires KiLCA built first
- **HDF5 outputs**: Include git version and timestamps for reproducibility
- **Modular physics**: Extensible for new collision models or zone types