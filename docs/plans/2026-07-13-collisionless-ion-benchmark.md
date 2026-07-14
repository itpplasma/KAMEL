# Collisionless-Ion Benchmark Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a global electrostatic KIM path that compares collisional Fokker--Planck ions with collisionless analytical Krook/Hamiltonian ions while retaining collisional, zero-FLR Fokker--Planck electrons and finite ion FLR.

**Architecture:** Keep `collision_model = 'FokkerPlanck'` as the electron and global-solver path, and add a default-preserving `ion_collision_model` selector. For the collisionless choice, calculate the existing FP matrices, replace only the ion charge-response matrices with explicitly collisionless, ion-only Krook matrices, assemble `K_total = K_e,FP + sum(K_i,C0)`, and perform one Poisson solve. Do not silently report the old FP ion current matrices as collisionless currents.

**Tech Stack:** Fortran 2008, CMake/Ninja, CTest, OpenMP, KIM's existing Gauss-Legendre kernel integration and UMFPACK Poisson solve.

---

### Task 1: Specify the Collisionless Plasma-Dispersion Argument

**Files:**
- Create: `KIM/tests/test_collisionless_ion_kernel.f90`
- Modify: `KIM/tests/CMakeLists.txt`
- Modify: `KIM/src/kernels/Krook_kernel_plasma_prefacs.f90`

**Step 1: Write the failing test**

Add a small test program that imports a wished-for pure function:

```fortran
use Krook_kernel_plasma_prefacs_m, only: Krook_collisionless_z0
z0 = Krook_collisionless_z0(om_E=2.0_dp, omega=0.5_dp, kpar=-3.0_dp, vT=4.0_dp)
expected = cmplx(-1.5_dp / (3.0_dp * sqrt(2.0_dp) * 4.0_dp), 0.0_dp, dp)
if (abs(z0 - expected) > 1.0e-14_dp) error stop
```

Also register `test_collisionless_ion_kernel` in `KIM/tests/CMakeLists.txt` with the standard `KIM_lib kilca_lib lapack cerf ddeabm slatec` link set.

**Step 2: Run test to verify it fails**

Run:

```bash
cmake --build build-ninja --target test_collisionless_ion_kernel -j 8
```

Expected: compilation fails because `Krook_collisionless_z0` is not defined.

**Step 3: Write minimal implementation**

Add a public pure function whose body is:

```fortran
val = cmplx(-(om_E - omega) / (abs(kpar) * sqrt(2.0_dp) * vT), 0.0_dp, dp)
```

Reject zero `kpar` or `vT` with a clear `error stop` in the non-pure caller rather than hiding a division by zero in the pure formula.

**Step 4: Run test to verify it passes**

Run the target and then:

```bash
ctest --test-dir build-ninja --output-on-failure -R '^test_collisionless_ion_kernel$'
```

Expected: PASS.

### Task 2: Make Krook Prefactors Accept an Explicit Collisionless Argument

**Files:**
- Modify: `KIM/tests/test_collisionless_ion_kernel.f90`
- Modify: `KIM/src/kernels/Krook_kernel_plasma_prefacs.f90`

**Step 1: Write the failing characterization test**

Extend the test to initialize a small in-memory electrostatic plasma and check both forms at one interior cell:

```fortran
z0_collisional = 0.5_dp * (plasma%spec(1)%z0(j) + plasma%spec(1)%z0(j+1))
z0_collisionless = Krook_z0_cc(j, plasma%spec(1), collisionless=.true.)
if (abs(aimag(z0_collisional)) <= 0.0_dp) error stop
if (abs(aimag(z0_collisionless)) > 1.0e-14_dp) error stop
if (minval(plasma%spec(0)%nu) <= 0.0_dp) error stop
```

**Step 2: Run test to verify it fails**

Expected: failure because `Krook_z0_cc` and the explicit collisionless selection do not exist.

**Step 3: Write minimal implementation**

Add `Krook_z0_cc(j, spec, collisionless)` and optional `collisionless` arguments to the four z0-dependent prefactors:

- `Krook_G1_rho_phi`
- `Krook_G2_rho_phi`
- `Krook_G1_rho_B`
- `Krook_G2_rho_B`

Default the optional argument to `.false.` so existing Krook behavior is bitwise unchanged. When true, form each endpoint's collisionless z0 from the local `om_E`, `kp`, `vT`, and global `omega`, then average the endpoints. Do not set `spec%nu = 0` and do not mutate `spec%z0`.

**Step 4: Run targeted and baseline tests**

Run:

```bash
ctest --test-dir build-ninja --output-on-failure -R '^(test_collisionless_ion_kernel|test_flr2_fourier_kernel)$'
```

Expected: both PASS.

### Task 3: Assemble Ion-Only Collisionless Krook Matrices

**Files:**
- Modify: `KIM/tests/test_collisionless_ion_kernel.f90`
- Modify: `KIM/src/kernels/kernel.f90`

**Step 1: Write the failing matrix-composition test**

On the same small in-memory plasma, initialize two `kernel_spl_t` charge kernels and call the wished-for API:

```fortran
call Kphi%init_kernel(xl_grid%npts_b, xl_grid%npts_b)
call KB%init_kernel(xl_grid%npts_b, xl_grid%npts_b)
call Krook_fill_kernel_phi_ions_collisionless(Kphi, KB)
```

Assert:

- `Kllp_e` remains exactly zero;
- every matrix entry is finite;
- at least one ion entry is nonzero;
- `Kllp == sum(Kllp_i, dim=3)` to `1e-12` relative tolerance;
- matrices are symmetric under `(l,lp)` exchange;
- at least one supported off-diagonal ion entry is nonzero with `ion_flr_scale_factor = 1`.

**Step 2: Run test to verify it fails**

Expected: compilation fails because `Krook_fill_kernel_phi_ions_collisionless` is missing.

**Step 3: Write minimal implementation**

Refactor `Krook_calc_kernel_rho_term_by_term` to accept an explicit inclusive species range and a collisionless flag. Preserve the existing `Krook_fill_kernel_phi` defaults. Add `Krook_fill_kernel_phi_ions_collisionless` that computes each ion separately into `Kllp_i(:,:,sp)`, mirrors the lower triangle, and sums only those ion matrices into `Kllp`.

Pass `collisionless=.true.` only to the z0-dependent Krook prefactors. Keep the finite `spec%rho_L` integration support unchanged.

**Step 4: Run targeted and baseline tests**

Run:

```bash
ctest --test-dir build-ninja --output-on-failure -R '^(test_collisionless_ion_kernel|test_kim_diagnostics|test_flr2_fourier_kernel)$'
```

Expected: all PASS.

### Task 4: Add the Default-Preserving Ion Model Selector

**Files:**
- Modify: `KIM/tests/test_collisionless_ion_kernel.f90`
- Modify: `KIM/src/setup/config_mod.f90`
- Modify: `KIM/src/general/read_config.f90`
- Modify: `KIM/src/general/config_display.f90`
- Modify: `KIM/src/util/IO_collection.f90`

**Step 1: Write failing config tests**

Write two small namelists through the test helper. Verify omission yields `ion_collision_model == 'FokkerPlanck'`, and explicit `ion_collision_model = 'collisionless'` is read exactly. Add a subprocess/negative test only if the repository has a stable pattern for validating `error stop`; otherwise keep validation in a focused pure predicate.

**Step 2: Run test to verify it fails**

Expected: failure because the new namelist field does not exist.

**Step 3: Write minimal implementation**

Declare:

```fortran
character(20) :: ion_collision_model = 'FokkerPlanck'
```

Add it to `/KIM_CONFIG/`, configuration display, and HDF5 metadata. Accept only `FokkerPlanck` and `collisionless`. Reject `collisionless` with adaptive theta integration for this first benchmark implementation, and reject it with an electron model other than `collision_model = 'FokkerPlanck'`.

**Step 4: Run targeted tests**

Expected: selector tests and baseline tests PASS.

### Task 5: Wire the Hybrid Charge Response into One Poisson Solve

**Files:**
- Modify: `KIM/tests/test_collisionless_ion_kernel.f90`
- Modify: `KIM/src/electrostatic_poisson/poisson.f90`
- Modify: `KIM/src/kernels/kernel.f90`

**Step 1: Write the failing end-to-end test**

Run two small global in-memory electrostatic solves with identical profiles and drive. Capture the collisional and hybrid `EBdat%Phi` plus exposed kernel-composition diagnostics. Assert for the hybrid case:

- electron collision frequencies remain positive;
- electron FP matrices match the collisional run;
- total charge matrices equal FP electron plus collisionless ion matrices;
- `Phi_m` is finite and nonzero;
- the selector changes the ion matrix and resulting `Phi_m` for the chosen non-degenerate fixture.

**Step 2: Run test to verify it fails**

Expected: the explicit selector still follows the all-FP assembly, so the expected ion-matrix difference assertion fails.

**Step 3: Write minimal implementation**

After `FP_fill_kernels` returns in `run_FP`, call the ion-only collisionless Krook fill when selected and rebuild only:

```fortran
K_rho_phi%Kllp = K_rho_phi%Kllp_e + sum(K_rho_phi%Kllp_i, dim=3)
K_rho_B%Kllp   = K_rho_B%Kllp_e   + sum(K_rho_B%Kllp_i, dim=3)
```

Perform the existing single total Poisson solve and retain `/fields/Phi_m`. For the hybrid case, write only the physically available collisional electron current diagnostic; do not label stale FP ion-current matrices as collisionless current.

**Step 4: Run the focused suite**

Run:

```bash
ctest --test-dir build-ninja --output-on-failure -R '^(test_collisionless_ion_kernel|test_kim_diagnostics|test_flr2_fourier_kernel|test_periodic_vs_global|test_sparse)$'
```

Expected: all PASS.

### Task 6: Create Reproducible AUG 33353 Benchmark Cases

**Files:**
- Create: `benchmark/inputs/33353_2900/manifest.sha256`
- Create: `benchmark/configs/m6_fp.nml`
- Create: `benchmark/configs/m6_collisionless.nml`
- Create: `benchmark/configs/m7_fp.nml`
- Create: `benchmark/configs/m7_collisionless.nml`
- Create: `scripts/prepare_inputs.py`
- Create: `scripts/run_kim.py`
- Create: `benchmark/analyze_phi.py`
- Create: `Makefile`
- Create: `README.md`

**Step 1: Add a failing harness smoke test**

Add a `--dry-run`/`--check-inputs` mode that must verify profile hashes, list exactly four core case IDs, locate resonances by interpolating `q(r) = -m/n`, and reject mismatched physics settings.

**Step 2: Run it and observe the missing-config failure**

Expected: FAIL until all four configs and the manifest exist.

**Step 3: Add immutable inputs and moderate configs**

Use AUG 33353 at 2.9 s from `/Users/markusmarkl/data/AUG/PROFILES/33353/2900/reff`, modes `(6,2)` and `(7,2)`, `ion_flr_scale_factor = 1`, zero electron FLR through the existing FP electron path, `l_space_dim = 128`, `rg_space_dim = 256`, and low-order Gauss quadrature. Keep every field except `ion_collision_model`, mode, and the mode-dependent radial window identical between paired runs. Python scripts must stage/verify inputs, generate configs, and run KIM; the Makefile is the public reproducibility interface with aggregate and per-case targets.

Record the KAMEL base commit and SHA-256 hashes. Copy inputs rather than relying on mutable external paths.

**Step 4: Verify the harness**

Run the dry check and one deliberately reduced smoke case. Expected: hashes pass, resonance is inside the domain, and finite nonzero `Phi_m` is produced.

### Task 7: Converge, Compare, and Report

**Files:**
- Create: `benchmark/results/manifest.json`
- Create: `benchmark/results/summary.csv`
- Create: `benchmark/results/phi_m6.png`
- Create: `benchmark/results/phi_m7.png`
- Create: `benchmark/REPORT.md`

**Step 1: Run the four moderate cases**

Capture configs, logs, timings, peak memory if available, and complex `Phi_m` files under immutable case directories.

**Step 2: Run convergence checks**

Double both radial grids, widen the radial domain, and increase quadrature once if needed. Compare only on the central common window. Require numerical changes to be smaller than the physical FP/FP versus FP/C0 difference, or report the difference as unresolved.

**Step 3: Compute metrics**

For each mode report complex relative L2 difference, normalized maximum difference, peak-magnitude ratio/shift, phase difference above a signal floor, complex value at resonance, and cost ratio.

**Step 4: Run final verification**

Build `KIM_exe`, run the focused CTest set, validate all hashes, rerun analysis from stored raw results, and confirm the generated CSV/plots/report are reproducible.

**Step 5: Record the conclusion**

Treat `m=7` as the primary away-from-zero result and `m=6` as the zero-sensitive result. Report metrics and numerical uncertainty without inventing a post-hoc threshold for “not large.”
