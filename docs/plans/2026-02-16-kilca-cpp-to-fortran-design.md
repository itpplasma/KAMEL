# KiLCA C/C++ to Fortran Translation — Design Document

## Motivation

Translate KiLCA's ~28K lines of C/C++ (117 files) into idiomatic Fortran 2003+ for
maintainability. The team is more fluent in Fortran; a single-language codebase is easier
to maintain and extend.

## Scope

### What gets translated
- 56 `.cpp` source files, 48 `.h` headers, 13 `.hpp` headers
- The `extern "C"` bridge layer (36+ functions) gradually disappears

### What stays as-is
- All existing Fortran modules (~90K lines): conductivity tensors, Maxwell equations,
  antenna equations, Bessel/AMOS library, SLATEC math, TOMS hypergeometric
- LAPACK/BLAS calls (already Fortran-native)
- External C libraries (SUNDIALS, GSL, zersol) — accessed via `ISO_C_BINDING` wrappers

### Starting point
Fresh branch from `main`. The old `35-translate-c-to-fortran` branch (86 commits, 242
commits behind main) is too diverged to merge but can be used as a reference for
translation patterns.

### End state
A single-language Fortran codebase (plus thin `ISO_C_BINDING` wrappers for
SUNDIALS/GSL/zersol) that produces bit-identical results to the current build on all
existing tests.

## Architecture

### Module Structure

```
KiLCA/src/                           (new Fortran modules replacing C++)
  kilca_settings_m.f90               <- core/settings.cpp/.h
  kilca_spline_m.f90                 <- spline/spline.cpp/.h
  kilca_interp_m.f90                 <- interp/interp.cpp/.h
  kilca_background_m.f90             <- background/background.cpp/.h + calc_back.cpp
  kilca_antenna_m.f90                <- antenna/antenna.cpp/.h (thin wrapper)
  kilca_solver_m.f90                 <- solver/solver.cpp, solver_a.cpp
  kilca_zone_m.f90                   <- mode/zone.cpp/.h (base type)
  kilca_flre_zone_m.f90              <- flre/flre_zone.cpp/.h
  kilca_imhd_zone_m.f90              <- imhd/imhd_zone.cpp/.h + compressible/incompressible
  kilca_hmedium_zone_m.f90           <- hom_medium/hmedium_zone.cpp/.h
  kilca_mode_m.f90                   <- mode/mode.cpp/.h + calc_mode.cpp
  kilca_eigmode_m.f90                <- eigmode/calc_eigmode.cpp/.h
  kilca_core_m.f90                   <- core/core.cpp/.h
  kilca_inout_m.f90                  <- io/inout.cpp
  kilca_post_proc_m.f90              <- post_processing/*.cpp
  kilca_wave_interface_m.f90         <- interface/wave_code_interface.cpp
  kilca_gsl_binding_m.f90            <- ISO_C_BINDING for GSL
  kilca_sundials_binding_m.f90       <- ISO_C_BINDING for SUNDIALS CVODe
  kilca_zersol_binding_m.f90         <- ISO_C_BINDING for zersol
```

### OOP Approach — Fortran 2003 Type Extension

```fortran
type :: zone_t                         ! base type
type, extends(zone_t) :: flre_zone_t
type, extends(zone_t) :: imhd_zone_t
type, extends(zone_t) :: hmedium_zone_t
```

Type-bound procedures for polymorphic dispatch (mirrors C++ virtual methods).

## Migration Strategy — Incremental, Bottom-Up

Translate leaf modules first, replacing one C++ file at a time while keeping the build
green. At each phase boundary, the full test suite must pass.

### Phase 1 — C Library Bindings (standalone, no KiLCA deps)

1. `kilca_gsl_binding_m.f90` — wrap `gsl_odeiv`, `gsl_deriv_central`
2. `kilca_sundials_binding_m.f90` — wrap CVODe create/init/solve, N_Vector ops
3. `kilca_zersol_binding_m.f90` — wrap root-finding entry points

### Phase 2 — Utilities (used by everything, no KiLCA deps)

4. `kilca_spline_m.f90` — cubic spline with derivatives (~500 LOC)
5. `kilca_interp_m.f90` — grid adaptation and interpolation (~230 LOC)
6. `kilca_settings_m.f90` — namelist reading, all config structs (~100 LOC)

### Phase 3 — Infrastructure

7. `kilca_background_m.f90` — profiles, equilibrium ODE via GSL binding (~1.5K LOC)
8. `kilca_antenna_m.f90` — thin wrapper, most logic already Fortran
9. `kilca_solver_m.f90` — basis vector integration via SUNDIALS binding (~1K LOC)

### Phase 4 — Zone Hierarchy

10. `kilca_zone_m.f90` — base type with grid, boundary conditions (~600 LOC)
11. `kilca_hmedium_zone_m.f90` — simplest zone type (~100 LOC)
12. `kilca_imhd_zone_m.f90` — incompressible/compressible MHD (~1.1K LOC)
13. `kilca_flre_zone_m.f90` — FLRE physics (~1.2K LOC)

### Phase 5 — Mode & Orchestration

14. `kilca_mode_m.f90` — mode structure, calc_mode, stitching (~2K LOC)
15. `kilca_eigmode_m.f90` — eigenmode search (~480 LOC)
16. `kilca_core_m.f90` — top-level data container (~200 LOC)

### Phase 6 — I/O & Interface

17. `kilca_inout_m.f90` — text file output (~735 LOC)
18. `kilca_wave_interface_m.f90` — QL-Balance coupling (~900 LOC)
19. `kilca_post_proc_m.f90` — post-processing (~1.2K LOC)

### Phase 7 — Entry Points

20. `main_linear.f90` — replaces `main_linear.cpp`
21. `main_eig_param.f90` — replaces `main_eig_param.cpp`

## Integration Pattern (C++/Fortran Coexistence)

During the transition, newly translated Fortran modules coexist with remaining C++ code.

When a C++ file is replaced (e.g., `spline.cpp` -> `kilca_spline_m.f90`), any remaining
C++ callers need a temporary `ISO_C_BINDING` shim:

```fortran
! Temporary bridge for remaining C++ callers
subroutine spline_evaluate_c(ptr, x, val) bind(C, name="spline_evaluate")
! TODO: remove when caller is translated
```

These shims monotonically decrease as translation progresses:
- Early phases (2-3) need the most shims
- Late phases (5-7) need almost none
- The existing 36+ `extern "C"` bridge functions become direct `use module` calls as
  their C++ callers are translated

At each step:
1. The module compiles and links into the existing build
2. `ctest` passes
3. The old C++ file is removed from `CMakeLists.txt`

## Testing & Validation

### Level 1 — Reference outputs (before any translation)

- Use `create_parabolic_profiles_from_res_surf()` from `python/utility/` to generate
  standard input profiles
- Run the current build on 2-3 representative cases (single-mode FLRE, multi-mode,
  eigenmode search)
- Store all text output files (background-data/, linear-data/, dispersion-data/) as
  reference

### Level 2 — CTest at every module translation

The existing `ctest` suite must pass after each module is replaced.

### Level 3 — Text file diff after each phase

Compare output `.dat` files against reference. KiLCA writes 16-digit scientific notation
(`%.16le`). A Python comparison script checks all values with relative error < 1e-12.

### Level 4 — Progress metric

Track remaining `extern "C"` bridge function count: 36+ -> 0.

## Risks & Mitigations

### SUNDIALS/GSL callback complexity
Both libraries use C function pointer callbacks. In Fortran these require `C_FUNPTR` and
`C_FUNLOC`. Mitigated by translating solvers early (Phase 3) and following SUNDIALS'
own Fortran examples.

### Pointer management during transition
Shims pass Fortran derived type pointers to remaining C++ via `C_LOC`/`C_F_POINTER`.
Mitigated by bottom-up order (shim count peaks early, decreases from there).

### Subtle numerical differences
Different compiler FP optimizations on Fortran vs C++ could cause tiny differences.
Mitigated by text-file diff with 1e-12 tolerance. Can compile with `-ffp-contract=off`
during validation if needed.

### Scope creep
Temptation to refactor while translating. Strict rule: **translate first, refactor later**.
The Fortran code should read like the C++ with Fortran syntax. Improvements go in
separate follow-up PRs.
