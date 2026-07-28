# Collisionless Finite-FLR Remainder Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Restore the analytically summed nonzero-cyclotron-harmonic remainder in the periodic collisionless ion `rho-Phi` kernel and prevent zero-FLR tests from masking regressions.

**Architecture:** Keep the existing zeroth-harmonic collisionless moments. Add the closed-form remainder from Markl Eq. (13.127) only to `rho_phi`, using the existing scaled Bessel factor and the canonical radial Gaussian. Validate complete observable kernel values rather than internal helper arithmetic.

**Tech Stack:** Fortran 2008, CMake/Ninja, CTest.

---

### Task 1: Add finite-FLR static regression tests

**Files:**
- Modify: `KIM/tests/test_collisionless_fourier_kernel.f90`

1. Add a homogeneous finite-FLR diagonal test with nonzero `ks`, `rho_L`, and `kr`.
2. Add a homogeneous finite-FLR off-diagonal test with `kr != krp`.
3. Retain the zero-FLR limit as a separate test.
4. Build and run `test_collisionless_fourier_kernel`; verify the two new tests fail because the current core returns the zeroth-harmonic Bessel factor.

### Task 2: Implement the analytic remainder

**Files:**
- Modify: `KIM/src/asymptotics/collisionless_fourier_kernel.f90`

1. Compute `exp(-rho_L^2*(kr-krp)^2/2)`.
2. Add `(sI0 - gaussian)/lambda_D^2/(8*pi^2)` to `rho_phi`.
3. Explain the decomposition and cite Eqs. (13.125)--(13.127) in comments.
4. Rebuild and verify the new static tests pass.

### Task 3: Strengthen the manufactured oracle

**Files:**
- Modify: `KIM/tests/test_collisionless_fourier_kernel.f90`

1. Update the independent finite-FLR oracle to include the analytic remainder.
2. Check `rho_phi`, `rho_B`, `j_phi`, and `j_B` independently.
3. Run the focused kernel and assembly suites.
4. Temporarily remove or negate the remainder and verify the finite-FLR tests fail, then restore it.

### Task 4: Full verification

1. Format changed Fortran with the repository formatter if available.
2. Build KIM and all tests.
3. Run full local CTest and record any acknowledged baseline failure separately.
4. Review the final diff for scope and derivation consistency.
