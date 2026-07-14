# Collisionless Taper Consistency Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make collisionless-ion kernel assembly use the same threshold-derived Larmor band as the Fokker--Planck kernel assembly.

**Architecture:** Extract the existing FP band-limit formula into one kernel-module helper that reads the configured Larmor factor and taper threshold and scans the available species Larmor radii. Use that helper in both FP and collisionless fillers. The collisionless filler will limit `(l,lp)` work using the shared global distance and will no longer apply its independent per-cell hard cutoff.

**Tech Stack:** Fortran 2008, CMake/CTest, KIM Gauss-Legendre kernel integration.

---

### Task 1: Add a real-path regression test

**Files:**
- Modify: `KIM/tests/test_collisionless_ion_assembly.f90`

1. Assemble the collisionless ion kernel with `kernel_taper_skip_threshold = 1e-6` and save the ion matrix.
2. Assemble it again with `kernel_taper_skip_threshold = exp(-1)`, which changes the shared band from about `18.6 rho_L,max` to `5 rho_L,max`.
3. Assert that the wider threshold produces resolved entries outside the narrow matrix support.
4. Build and run `test_collisionless_ion_assembly`; expect failure because the collisionless filler currently ignores the threshold.

### Task 2: Share and use the FP band limit

**Files:**
- Modify: `KIM/src/kernels/kernel.f90`

1. Extract the existing formula `alpha * rhoT_max * sqrt(log(1/tau))` into a module procedure.
2. Replace the duplicated FP calculation with the helper.
3. Compute the same limit in `Krook_fill_kernel_phi_ions_collisionless`.
4. Restrict the collisionless `(l,lp)` loop to the shared band and remove the per-cell `distance > alpha*rhoT` cutoff.
5. Rebuild and rerun the regression test; expect it to pass.

### Task 3: Verify the affected system

**Files:**
- Verify: `KIM/src/kernels/kernel.f90`
- Verify: `KIM/tests/test_collisionless_ion_assembly.f90`

1. Run the focused collisionless kernel and assembly tests.
2. Run the broader focused KIM CTest selection used by the benchmark.
3. Inspect `git diff --check` and the final diff for unrelated changes.
4. Rerun the collisionless benchmark cases and regenerate the comparison outputs because the ion-matrix support changes materially.
