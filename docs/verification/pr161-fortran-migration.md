# PR #161 migration repair validation

The repaired branch integrates `main` at `d4d315ef`, including its Lorentz-force
correction. The obsolete C++ torque test is translated to Fortran. Golden baselines
and comparison tolerances are unchanged.

## Repairs and regression coverage

| Repair | Production regression |
| --- | --- |
| Correct zero-based current interpolation offsets | `kilca_current_interpolation` covers every channel, endpoints, and intermediate radii |
| Use complex LAPACK QR query and workspace with an explicit interface | `test_solver` checks analytic independent solutions with and without QR normalization |
| Preserve the C complex ABI on macOS | `test_zersol_bridge` finds a complex root through Fortran callbacks |
| Read native directory entries through a small C adapter | `test_directory` checks enumeration and actual eigenmode cleanup, including retained folders and quoted paths |
| Restore real and imaginary output columns and density filenames | `kilca_complex_quantity_output` checks the production writers |
| Trim padded matrix filename prefixes | `test_inout` checks filenames and complex matrix contents |
| Accept tab-separated settings tokens | `test_back_sett` and `test_eigmode_sett` exercise tabs |
| Release the system-matrix spline on destruction | `test_sysmat_profiles` requires exactly one release of the correct handle |
| Preserve current main's corrected Lorentz force and cylindrical integration | `kilca_lorentz_torque` checks components, phase invariance, and signed integrals |

The interpolation, output, directory/cleanup, matrix filename, settings, and spline
regressions were reproduced before their repairs. A guarded LAPACK query reproduced
the original eight-byte overwrite; the explicit complex interface now prevents
passing the original real query buffer. The original C++ bridge failed compilation
on macOS; its repaired callback test compiles and runs.

## Additional integration failure

Fresh end-to-end validation exposed a gfortran 16.1.0 failure on macOS ARM64:
passing an element of one allocatable component together with a section of another
from a polymorphic object generated a bounds check using an uninitialized register.
This was reduced to an independent small Fortran program and inspected with LLDB.
It was not a numerical solver failure.

The affected zone and stitching calls now pass the first element of each contiguous
column or plane to their existing explicit-shape or assumed-size dummies. Standard
sequence association preserves the original layout and avoids the faulty section
lowering. Independent review checked the extents and interfaces.

`kilca_full_output` now runs the real FLRE/vacuum executable against the local
`flre_m6n2` input deck with all quantity outputs enabled. It checks all 48 complex
current files, three density files, 18 torque files, finite values, radial grids,
imaginary components, and the combined electromagnetic fields. It needs neither
network access nor a second build; the separate golden-record suite is unchanged.

## Fresh local results

On macOS ARM64 with Homebrew gfortran 16.1.0:

```text
cmake --build build --parallel 8       PASS
ctest --test-dir build --output-on-failure
                                      78/78 passed
make pytest                           12/12 passed
```

The final local build uses the normal Release flags, including bounds checking;
temporary unoptimized objects used during diagnosis were rebuilt.

A separate clean current-main C++ build used the same pinned dependency sources.
Both executables completed the full-output deck. Across all 114 linear-data files,
radial grids and shapes matched exactly and every value was finite. The largest
physical-output error, normalized by the peak magnitude of its individual output
column, was `9.146e-10`. The energy-balance residual satisfied the existing absolute
self-consistency threshold in both builds.

This expanded output comparison is not bitwise identity: the unchanged generic
pointwise comparator flags 13 of 168 files (including input/background files in
that count) at `rtol=1e-7`. Its worst relative difference is about `4.12e-5` in an
integrated dissipated-power value near a cancellation. Peak-normalized error remains
below `9.146e-10` for every physical output column, below the deck's `1e-8` solver
tolerance. These differences are recorded rather than hidden by changing golden
tolerances.

Fresh hosted Ubuntu/macOS builds, tests, and strict golden comparisons are attached
to [PR #161](https://github.com/itpplasma/KAMEL/pull/161/checks). The pre-fortnum job
is diagnostic and permits failure; its green job status alone is not evidence of
numerical equivalence.

## Native Fortran interface cleanup (2026-09-07)

The follow-up cleanup replaces internal C linkage with module imports and native
Fortran interfaces. Production KiLCA/QL-Balance `bind(C)` declarations decreased
from 443 to 39, counting procedure definitions, interface declarations, and types
and excluding tests. The remaining declarations describe libc/POSIX, SUNDIALS,
fortnum, and ZerSol interfaces or their callbacks and context types.

Settings readers, constructors, path getters, and file writers accept native
character strings. Adaptive-grid and solver callbacks use native procedure
arguments. The solver and FLRE zone share the actual parameter types. Opaque
handles retain their existing ownership. Legacy external Fortran procedures have
canonical interfaces in `legacy_interfaces_m.f90`, including their true complex
array types and compiler-managed character lengths. Eigenmode orchestration uses
a submodule to avoid a core/eigenmode dependency cycle; legacy FLRE settings data
is compiled separately from its higher-level adapters.

Explicit interfaces also required correcting previously unchecked calls: the
background dimension provider is a subroutine; Maxwell system-index retrieval
takes one argument; conductivity fills contiguous background-array sections;
wave fields and currents return complex arrays; and per-mode dissipation selects
the corresponding scalar core handle. Unused field outputs now have separate
storage instead of aliasing the same output argument repeatedly.

The updated regressions cover native settings strings (including embedded spaces
and short output buffers), callback replacement between integrations, typed
profile test providers, and short paths in all three zone factories. Restoring
the missed C-string scan in an isolated homogeneous-medium object makes the new
zone-path regression fail; the corrected production implementation passes.

Fresh local Release validation uses gfortran 16.1.0 with runtime checking. The
full project builds, all 79 CTests pass, and all 12 Python tests pass. All 114
numeric linear-data files in the full-output FLRE/vacuum case are exactly equal
to the saved pre-cleanup Fortran build at `1f655463`, including their radial grids
and every output column. No golden inputs, baselines, or tolerances changed.

Independent specification and code-quality reviews found no remaining blockers.
Hosted validation must be checked on the final cleanup commit before merging.
