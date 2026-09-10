# PR #161 migration repairs implementation plan

**Goal:** Fix the eight confirmed migration regressions, preserve current main's physics,
and validate the integrated branch with targeted and full regression checks.

**Architecture:** Keep the existing solver interfaces and numerical algorithms. Use a small
portable directory adapter for OS-owned structures and correctly typed LAPACK workspace.
Exercise production routines in CTest; preserve the existing complex output schema.

**Tech stack:** Fortran, C/C++ interoperability, CMake/CTest, Python golden harness.

1. Merge current main. Resolve build registrations and translate the newer Lorentz-torque
   implementation and test into the Fortran branch before deleting obsolete C++ sources.
2. Add current interpolation and complex output tests in `KiLCA/tests/`; demonstrate the
   shifted samples and missing output columns, then fix `flre_quants_m.f90`.
3. Reproduce the default macOS bridge build failure, make Zersol complex conversions portable,
   and register a complex callback round-trip regression test.
4. Test portable directory enumeration and cleanup retention, then replace native `dirent`
   layouts in `mode_data_m.f90` and `main_eig_param.f90` with shared portable routines.
5. Reproduce the LAPACK query overwrite with a guard, enforce a complex workspace interface,
   and add an analytic basis-integration regression covering both QR paths.
6. Test real matrix output and whitespace-separated settings, then fix `inout_m.f90`,
   `background_m.f90`, and `eigmode_sett_m.f90` with minimal changes.
7. Extend the sysmat lifecycle regression to require releasing its spline handle, then fix
   `sysmat_profs_m.f90` destruction.
8. Run targeted tests after each fix, then fresh macOS build and complete CTest, Python
   harness tests, strict golden comparison, formatting checks, and independent review.
9. Commit and push to the PR branch; inspect fresh Linux/macOS CI and resolve failures.

Use failing regression tests before each production fix. Keep independent work on disjoint
files, coordinate CMake registrations centrally, and record actual commands/results.
