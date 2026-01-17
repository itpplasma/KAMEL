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
- **Common utilities** - Shared math/utils in `/common/`

### Key Directories
- `/PreProc/` - Preprocessing utilities (fouriermodes, neo-2 templates)
- `/build/` - Build artifacts and compiled binaries (gitignored)

### Data Flow & Workflow
All data exchange uses **HDF5 format** for standardization. The typical workflow:
1. **Profile preparation**: Prepare input profiles in CGS units
2. **Main run**: Execute solver (KiLCA/KIM/QL-Balance)
3. **Post-processing**: Python analysis and visualization
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

## KIM Profile Input System

### Unit Requirements (CGS)
All KIM input profiles must be in CGS units:
- **Density**: 1/cm³ (typically 10^12 to 10^15). **NOT** SI units (10^19 1/m³)
- **Temperature**: eV (typically 10 to 20000 eV)
- **Electric field (Er)**: statV/cm (typically ±0.5 statV/cm)
- **Magnetic field**: Gauss (typically ~20000 G for tokamaks)
- **Radial coordinate**: cm (effective radius r_eff)

### Profile Coordinate Types (`KIM_PROFILES` namelist)
```fortran
&KIM_PROFILES
    coord_type = 'auto'           ! 'auto', 'sqrt_psiN', or 'r_eff'
    input_profile_dir = './profiles/'
    geqdsk_file = './gfile.eqdsk'  ! For equilibrium computation
/
```

- `'auto'` - Auto-detect from max radius (>2 cm → r_eff, else sqrt_psiN)
- `'sqrt_psiN'` - Profiles in √(ψ/ψ_max) coordinates, transformed to r_eff using equilibrium
- `'r_eff'` - Profiles already in effective radius [cm], used directly

### Required Profile Files
Located in `profile_location` directory:
- `n.dat` - Electron density (r_eff [cm], n [1/cm³])
- `Te.dat` - Electron temperature (r_eff [cm], Te [eV])
- `Ti.dat` - Ion temperature (r_eff [cm], Ti [eV])
- `q.dat` - Safety factor (r_eff [cm], q)
- `Er.dat` - Radial electric field (r_eff [cm], Er [statV/cm]) - optional, calculated if missing
- `Vz.dat` - Toroidal rotation (r_eff [cm], Vz [cm/s]) - optional

### Automatic Validation Checks
KIM performs these checks on startup:
1. **Density units** - Error if density >10^17 (likely SI instead of CGS)
2. **q vs m_mode sign** - Warning if q>0 with m>0 (no resonance expected)
3. **Radial range** - Error if r_min or r_plas outside profile range
4. **Er interpolation** - Automatic interpolation if Er.dat grid differs from other profiles

### Er Calculation
If `Er.dat` is not provided, KIM calculates it from radial force balance:
```
Er = (Ti/e·n)·dn/dr + (1/e)·dTi/dr + (r·B0·Vz)/(c·q·R0)
```
Output written to `Er_no_Vpol.dat` (without poloidal rotation contribution).

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