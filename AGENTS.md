# AGENTS.md

This file provides guidance to AI coding agents (Claude Code, Codex, etc.) working in this repository.

## Project Overview

KAMEL (Kinetic plAsma response ModEL) is a scientific computing framework for modeling plasma response to external magnetic perturbations in fusion plasmas. It consists of three main codes:

1. **KiLCA** - Cylindrical linear plasma response solver using finite Larmor radius formalism
2. **KIM** - KiLCA Integral Model using integral formalism
3. **QL-Balance** - Quasilinear 1D radial transport code

## Project Structure

- **Top-level build config**: `CMakeLists.txt`, `Makefile`; build outputs in `build/`
- **Solver codes**: `KiLCA/`, `KIM/`, `QL-Balance/` (executables use `.x` suffix, e.g. `KIM.x`)
- **Shared math/utils**: `common/`
- **Python interfaces**: `python/` (KAMELpy) and `KIM/python/` (supported KIM API/CLI)
- **Preprocessing tools**: `PreProc/` (fouriermodes, neo-2 templates)
- **KIM config namelist**: `KIM/nmls/KIM_config.nml`
- **QL-Balance config**: `balance_conf.nml` (in run directory)

## Build, Test, and Development Commands

```bash
# Build all three codes (Release mode by default)
make all

# Build in Debug mode
CONFIG=Debug make all

# Build individual components
make KiLCA
make KIM
make QL-Balance

# Run all tests via CTest from root
make test  # invokes ctest --test-dir build

# Install Python interface
(cd python && make init && make install)

# Install and test the standalone KIM API/CLI
python -m pip install -e './KIM/python[test]'
python -m pytest KIM/python/tests/unit

# Clean build (removes build/)
make clean
```

**Build notes:**
- Build system is CMake + Ninja (`build/build.ninja`); ensure Ninja is on PATH.
- `QL-Balance` depends on `KiLCA` artifacts — build KiLCA first.
- Compilers: MPI Fortran (`mpif90`), C/C++ with clang-format support.
- Platforms tested: Apple Silicon (clang 16.0 + gfortran 14.2), Debian (GNU 12.2.0).
- Reconfigure tip: when in doubt, `make clean` then rebuild.
- With multiple worktrees, confirm that `python -c "import kim; print(kim.__file__)"` resolves to
  the intended worktree before trusting API test results.

## Key Dependencies

### External Libraries

- **MPI** - Parallel computing
- **LAPACK/BLAS** - Linear algebra
- **SuiteSparse** - Sparse matrix operations
- **HDF5** - Data storage
- **fortnum** - Numerical core (special functions, quadrature, ODE, root finding)
- **SUNDIALS** - Numerical solvers

### Python Dependencies

- **Core**: numpy, scipy, h5py, f90nml
- **Visualization**: matplotlib

## Data Flow & Workflow

Solver configuration uses Fortran namelists, and profile inputs are normally two-column text files.
The supported KIM Python interface accepts validated JSON requests, stages immutable input copies,
and records logs and manifests. Scientific solver results use HDF5 and include version/provenance
metadata where supported.

Typical workflow:

1. **Profile preparation**: Prepare input profiles in CGS units
2. **Main run**: Execute solver (KiLCA/KIM/QL-Balance)
3. **Post-processing**: Python analysis and visualization via KAMELpy

## Python Interfaces

The repository contains two distinct Python distributions:

- **`kamel-kim` / `kim` (`KIM/python`)** - Supported KIM configuration, validation, execution,
  sweeps, and HDF5 result API/CLI. Its main entry points are `SimulationConfig`, `Simulation`,
  `Result`, and `run_sweep`.
- **KAMELpy (`python/`)** - KiLCA and QL-Balance workflow interfaces plus experimental KIM
  dispersion and field-solver analysis modules. Do not use the experimental `KIMpy` analysis
  modules as substitutes for the supported orchestration API.

### KAMELpy Pattern

```python
interface = KiLCA_interface(shot, time, path, run_type, machine)
interface.set_modes(m_modes, n_modes)
interface.prepare_balance_input(input_file)
interface.run_balance()
```

For supported KIM automation, start with `KIM/python/README.md` and
`KIM/python/docs/request-json.md`. Keep CLI commands as thin adapters over the public Python API.

## KIM Profile Input System

The Fortran solver supports the coordinate modes below. Version 1 of the supported `kim` Python
API accepts only explicit `r_eff` profiles; it must not silently rely on Fortran auto-detection or
perform an undeclared coordinate conversion.

### Unit Requirements (CGS)

Prepared `r_eff` profiles consumed by KIM use the following units. The Fortran preprocessing path
may accept a normalized `sqrt_psiN` radial coordinate, but transforms it to `r_eff` before solving.

- **Density**: 1/cm^3 (typically 10^12 to 10^15). **NOT** SI units (10^19 1/m^3)
- **Temperature**: eV (typically 10 to 20000 eV)
- **Electric field (Er)**: statV/cm (typically +/-0.5 statV/cm)
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

- `'auto'` - Auto-detect from max radius (>2 cm -> r_eff, else sqrt_psiN)
- `'sqrt_psiN'` - Profiles in sqrt(psi/psi_max) coordinates, transformed to r_eff using equilibrium
- `'r_eff'` - Profiles already in effective radius [cm], used directly

### Required Profile Files

Located in `profile_location` directory:

- `n.dat` - Electron density (r_eff [cm], n [1/cm^3])
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
Er = Ti/(e*n)*dn/dr + (1/e)*dTi/dr + (r*B0*Vz)/(c*q*R0)
```

Output written to `Er_no_Vpol.dat` (without poloidal rotation contribution).

## Coding Style

- **Indentation**: 4 spaces (no tabs)
- **Line width**: 100 columns maximum
- **Fortran formatting**: Use `fprettify` (see `.fprettify` config)
  - Module names use `_m` suffix
  - Test files named `test_*.f90`
  - Executables use `.x` suffix
- **C/C++ formatting**: Use `clang-format` (LLVM base style, see `.clang-format`)

## Testing

- **Framework**: CTest enabled at top level; tests run from `build/` via `ctest`.
- **KIM Python unit tests**: install `./KIM/python[test]`, then run
  `python -m pytest KIM/python/tests/unit`. These tests use a fake executable.
- **KIM Python integration test**: opt in with `KIM_RUN_INTEGRATION=1` and an absolute
  `KIM_EXECUTABLE` path. A skipped integration test is not evidence that the real solver passed.
- **Example**: `QL-Balance/src/test/test_sparse.f90` registered in `QL-Balance/src/test/CMakeLists.txt` with `add_test`.
- **New tests**: Add via CMake in relevant subproject.
  - Name sources `test_*.f90`; register with `add_test` in `CMakeLists.txt`.
  - Keep self-contained, avoid interactive I/O.
  - Print diagnostics for debugging; prefer explicit diagnostics in failing cases.
- **Golden-record regression** (KIM, QL-Balance, KiLCA): an A/B two-build
  comparison lives in `test/golden/` and runs in the dedicated `golden-record`
  GitHub Actions job (`.github/workflows/golden.yml`). It is intentionally
  excluded from local `ctest` / `make test` so normal builds stay fast and
  offline. It compares a frozen `golden-baseline` tag against the PR head; see
  `test/golden/README.md` for how to add a case or re-bless the baseline.

## Commit & Pull Request Guidelines

- **Commit style**: Conventional Commits (`feat:`, `fix:`, `refactor:`, `chore:`)
  - Optional scope: `feat(KIM): add new solver`
  - Single-topic, buildable commits
- **Pull Requests**:
  - Include purpose, key changes, usage notes
  - Link related issues
  - Provide test output (`ctest` summary)
  - Require green CI before merge
  - Avoid force-push after review (except for rebase)

## Agent Guardrails

- Do not revert unrelated local edits.
- Prefer focused, minimal changes over broad refactors.
- When adding features, add or update tests where practical.
- Reconfigure with `make clean` if CMake cache/dependency state looks inconsistent.
- Modular physics: code is extensible for new collision models or zone types.
- For KIM data adoption, preserve original source arrays/files and record units, coordinates,
  sign/phase conventions, transformations, hashes, and provenance in staged artifacts.
- Never guess a scientific unit, coordinate mapping, sign convention, normalization, missing
  experimental input, or comparison tolerance. Require an explicit reviewed contract and
  independent reference values before implementing a scientific conversion.
- Preparation/import operations must not launch the solver or overwrite an existing destination.

## Known Issues

### Stale `forces_nl` in `rhs_balance` Jacobian probing (QL-Balance)

In `rhs_balance_m.f90`, the `rhs_balance` subroutine computes `forces_nl` in a pre-loop over all boundary points but only retains the value from the last point (`ipoi = npoib`). This stale value is then reused for all boundary points inside the Jacobian probing loop (lines ~280-283), producing incorrect nonlinear QL fluxes for the torque computation at interior boundary points. The impact is limited to the Jacobian accuracy for the implicit solver and may cause slower convergence or subtle inaccuracies in the nonlinear torque terms during probing.
