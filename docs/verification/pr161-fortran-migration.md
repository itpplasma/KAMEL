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
