# Repository Guidelines

## Project Structure & Module Organization
- Root CMake/Make build: `CMakeLists.txt`, `Makefile`, build outputs in `build/`.
- Codes: `KiLCA/`, `KIM/`, `QL-Balance/` (executables like `KIM.x`).
- Shared math/utils: `common/`.
- Interfaces and tools: `python/` (KAMELpy), `PreProc/`, `utility_scripts/`, `matlab/`.
- Docs and templates: `Documentation/`, `template_scripts/`.

## Build, Test, and Development Commands
- Configure+build (Release default): `make all` or `CONFIG=Debug make all`.
- Component builds: `make KiLCA`, `make KIM`, `make QL-Balance`.
- Generator: Ninja; build files in `build/build.ninja`.
- Run tests: `make test` (invokes `ctest --test-dir build`).
- Python interface: `cd python && make init && make install`.
- Clean: `make clean` (removes `build/`).

## Coding Style & Naming Conventions
- Indent 4 spaces; wrap at 100 cols.
- Fortran: format with `fprettify` (see `.fprettify`); prefer module names with `_m` suffix; test files `test_*.f90`.
- C/C++: format with `clang-format` (LLVM base, `.clang-format`).
- Executables use `.x` suffix (e.g., `test_sparse.x`).

## Testing Guidelines
- CTest is enabled at the top level; tests run from `build/` via `ctest`.
- Example test target: `QL-Balance/src/test/test_sparse.f90` registered in `QL-Balance/src/test/CMakeLists.txt` with `add_test`.
- Add new tests via CMake in the relevant subproject; name sources `test_*.f90` and binaries `*.x`.
- Keep tests self-contained; print diagnostics and avoid interactive I/O.

## Commit & Pull Request Guidelines
- Use Conventional Commits: `feat:`, `refactor:`, `chore:`, optional scope (e.g., `feat(KIM): ...`).
- Commits should be single-topic and buildable.
- PRs: include purpose, key changes, usage/build notes, and link issues. Provide test output snippets (e.g., `ctest` summary) and screenshots only when visual.
- Require green CI before merge; avoid force-push after review except to rebase.

## Configuration & Data
- Compilers/libs: MPI Fortran (`mpif90`), HDF5, LAPACK/BLAS, SuiteSparse, GSL, SUNDIALS.
- Ninja-based builds; ensure Ninja is available on PATH.
- I/O uses HDF5 across codes; templates in `template_scripts/` prepare inputs and post-process outputs.
- Namelists: KIM `KIM/nmls/KIM_config.nml`; QL-Balance `balance_conf.nml` (in run directory).
- Dependency: QL-Balance requires KiLCA to be built first.
- Reconfigure tip: when in doubt, `make clean` then rebuild; HDF5 outputs record git version/timestamps for reproducibility.

## Workflows & Interfaces
- Typical flow: prerun (generate HDF5) → run solver (KiLCA/KIM/QL-Balance) → Python/MATLAB post-processing.
- Python (KAMELpy): install via `python/Makefile`; use OO interfaces in `python/` for analysis and orchestration.
