# CODEX.md

This file captures the minimum operating context for coding agents working in this repository.

## Purpose
- KAMEL is a fusion-plasma modeling framework with three main solvers:
  - `KiLCA/` (finite Larmor radius response)
  - `KIM/` (integral formalism)
  - `QL-Balance/` (quasilinear radial transport)

## Build And Test
- Default build:
  - `make all`
- Debug build:
  - `CONFIG=Debug make all`
- Build one component:
  - `make KiLCA`
  - `make KIM`
  - `make QL-Balance`
- Run tests:
  - `make test` (runs `ctest --test-dir build`)
- Clean:
  - `make clean`

Notes:
- Build system is CMake + Ninja (`build/build.ninja`).
- `QL-Balance` depends on `KiLCA` artifacts.

## Key Paths
- Top-level build config: `CMakeLists.txt`, `Makefile`
- Shared code: `common/`
- Preprocessing tools: `PreProc/`
- Python interface (KAMELpy): `python/`
- KIM config namelist: `KIM/nmls/KIM_config.nml`

## Code Style
- Indentation: 4 spaces; max line width: 100.
- Fortran formatting: `fprettify` (config in `.fprettify`).
- C/C++ formatting: `clang-format` (config in `.clang-format`).
- Fortran tests should follow `test_*.f90`.

## Testing Conventions
- Add tests through CMake in the relevant subproject.
- Keep tests self-contained and non-interactive.
- Prefer explicit diagnostics in failing cases.

## Runtime/Data Conventions
- HDF5 is the standard output/input format across components.
- Typical workflow:
  1. Prepare profiles
  2. Run solver(s)
  3. Post-process with Python interfaces

## Dependency Expectations
- MPI Fortran (`mpif90`), HDF5, BLAS/LAPACK, SuiteSparse, GSL, SUNDIALS.
- CMake and Ninja must be available.

## Agent Guardrails
- Do not revert unrelated local edits.
- Prefer focused, minimal changes over broad refactors.
- When adding features, add or update tests where practical.
- Reconfigure with `make clean` if CMake cache/dependency state looks inconsistent.
