# KIM Radial-Current Fourier Kernels Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement all three harmonic-aware radial-current Fourier kernels and quantify their cyclotron-harmonic convergence.

**Architecture:** Add a focused `radial_current_fourier_kernel_m` module that consumes KIM's equilibrium and susceptibility data without changing the existing charge and parallel-current kernels. Retain the missing `I03` moment, expose per-harmonic and cutoff-summed APIs, then integrate the two presently driven columns into the periodic result as `jrad`; keep the tested `B_parallel` column callable for future use.

**Tech Stack:** Fortran 2008, CMake/Ninja, CTest, fortnum modified-Bessel functions, KIM periodic Fourier assembly.

---

### Task 1: Retain the I03 susceptibility moment

**Files:**
- Modify: `KIM/src/background_equilibrium/species_mod.f90`
- Modify: `KIM/src/background_equilibrium/periodic_background.f90`
- Test: `KIM/tests/test_radial_current_fourier_kernel.f90`
- Modify: `KIM/tests/CMakeLists.txt`

**Step 1: Write the failing test**

Add a focused executable that initializes the small in-memory equilibrium with
`mphi_max=2` and checks allocation and equality of `I03(:,ell)` and
`symbI(0,3)` recomputed from `x1,x2` for `ell=-2:2`.

**Step 2: Run test to verify it fails**

Run: `cmake --build build --target test_radial_current_fourier_kernel && ctest --test-dir build -R '^test_radial_current_fourier_kernel$' --output-on-failure`

Expected: compilation fails because `species_t%I03` does not exist.

**Step 3: Write minimal implementation**

Add `I03` and `I03_cc` to `species_t` and update allocation, population,
diagnostic writing, periodization/copying, and deallocation paths alongside
`I01` and `I13`.

**Step 4: Run test to verify it passes**

Run the focused CTest command above. Expected: PASS.

**Step 5: Commit**

Commit: `feat(KIM): retain I03 susceptibility moment`

### Task 2: Implement scaled arbitrary-order Bessel values

**Files:**
- Create: `KIM/src/asymptotics/radial_current_fourier_kernel.f90`
- Modify: `KIM/src/CMakeLists.txt`
- Test: `KIM/tests/test_radial_current_fourier_kernel.f90`

**Step 1: Write the failing test**

Test scaled `exp(-bplus) I_ell(bcross)` for `ell=-4:4` against direct
`bessel_in` values for moderate arguments, integer-order symmetry, zero
arguments, and finiteness at `bcross=1000` with `bplus>=bcross`.

**Step 2: Verify RED**

Run the focused test. Expected: missing module/procedure failure.

**Step 3: Write minimal implementation**

Implement the arbitrary-order scaled Bessel helper with direct evaluation for
moderate arguments and an overflow-safe large-argument branch.

**Step 4: Verify GREEN**

Run the focused test. Expected: PASS.

**Step 5: Commit**

Commit: `feat(KIM): add scaled harmonic Bessel helper`

### Task 3: Implement radial FLR coefficients

**Files:**
- Modify: `KIM/src/asymptotics/radial_current_fourier_kernel.f90`
- Test: `KIM/tests/test_radial_current_fourier_kernel.f90`

**Step 1: Write the failing tests**

For several non-singular `(ell,ks,kr,ksp,krp)` cases, compare `O0,O2` with an
independent centred response-side derivative and `W0,W2` with an independent
mixed derivative of the ordinary coefficient including phase terms. Add the
documented zero-wavenumber anchors for one and two radial insertions.

**Step 2: Verify RED**

Run the focused test. Expected: missing radial coefficient procedure.

**Step 3: Write minimal implementation**

Compute `Gamma`, `Lambda`, and their analytic source, response, and mixed
derivatives using neighboring scaled Bessel orders. Form the phase-stripped
`O` and `W` coefficients, with a regular adjacent-order zero-wavenumber path.

**Step 4: Verify GREEN**

Run the focused test. Expected: PASS.

**Step 5: Commit**

Commit: `feat(KIM): implement radial FLR coefficients`

### Task 4: Implement all three harmonic kernels

**Files:**
- Modify: `KIM/src/asymptotics/radial_current_fourier_kernel.f90`
- Test: `KIM/tests/test_radial_current_fourier_kernel.f90`

**Step 1: Write the failing tests**

Add independent formula checks for one species/harmonic for `jrad-Phi`,
`jrad-Br`, and `jrad-Bparallel`. Check source moments `(I00,I02)`,
`(I01,I03)`, `(I00,I02)`, normalization, complex prefactors, signed `omega_c`,
and the Fourier/harmonic phase.

**Step 2: Verify RED**

Run the focused test. Expected: missing kernel procedures.

**Step 3: Write minimal implementation**

Add per-species/per-harmonic cores and public species-summed functions. Accept
an optional cutoff, validate it, and sum symmetrically over harmonics.

**Step 4: Verify GREEN**

Run the focused test. Expected: PASS.

**Step 5: Commit**

Commit: `feat(KIM): add radial-current Fourier kernels`

### Task 5: Assemble and expose periodic radial current

**Files:**
- Modify: `KIM/src/electrostatic_poisson/periodic_assembly.f90`
- Modify: `KIM/src/electrostatic_poisson/periodic_solve.f90`
- Modify: `KIM/src/electrostatic_poisson/poisson_periodic.f90`
- Modify: `KIM/src/electrostatic_poisson/fields_mod.f90`
- Modify: `KIM/src/general/kim_solver_m.f90`
- Test: `KIM/tests/test_periodic_assembly.f90`
- Test: `KIM/tests/test_periodic_solve.f90`

**Step 1: Write the failing tests**

Require periodic assembly to produce `KjrPhi` and `KjrBr`, reconstruction to
return `jrad`, and solver results to expose a finite `jrad` vector. Confirm the
assembled matrices equal direct quadrature of the new kernels.

**Step 2: Verify RED**

Run: `cmake --build build --target test_periodic_assembly test_periodic_solve && ctest --test-dir build -R 'test_periodic_(assembly|solve)' --output-on-failure`

Expected: compile failures for missing radial-current arguments/results.

**Step 3: Write minimal implementation**

Thread the two active radial-current matrices through periodic assembly,
reconstruct `jrad=KjrPhi*Phi+KjrBr*Br`, store it in `EBdat`, write it with the
other fields, and copy it into `kim_results_t`. Do not invent a `Bparallel`
field; its kernel remains separately callable.

**Step 4: Verify GREEN**

Run the focused periodic tests. Expected: PASS.

**Step 5: Commit**

Commit: `feat(KIM): expose periodic radial current`

### Task 6: Measure cyclotron-harmonic convergence

**Files:**
- Create: `KIM/tests/test_radial_current_harmonics.f90`
- Modify: `KIM/tests/CMakeLists.txt`
- Optionally create: `docs/analysis/2026-08-13-radial-current-harmonics.md`

**Step 1: Write the convergence harness**

Initialize the representative static deuterium case once with at least
`mphi_max=4`. Evaluate cutoff-summed versions of all three kernels on a fixed
Fourier/grid sample and assemble the active periodic matrices. Print each
shell `K_N-K_(N-1)` relative to the converged norm and the resulting change in
`jrad`.

**Step 2: Run the harness**

Run: `cmake --build build --target test_radial_current_harmonics && ctest --test-dir build -R '^test_radial_current_harmonics$' -V`

Expected: finite results, explicit nonzero `+/-1` contribution, and decreasing
tail or a clearly reported lack of convergence.

**Step 3: Record findings**

Document the case parameters, norm definitions, harmonic shells, and whether
each kernel and reconstructed `jrad` is adequately converged.

**Step 4: Commit**

Commit: `test(KIM): quantify radial-current harmonics`

### Task 7: Final verification

**Files:** all changed files.

**Step 1: Format**

Run `fprettify` on new and modified Fortran sources using the repository
configuration, then run `git diff --check`.

**Step 2: Run focused tests**

Run all radial-current, Fourier-kernel, periodic assembly, periodic solve, and
solver tests with `--output-on-failure`.

**Step 3: Run full tests**

Run: `ctest --test-dir build --output-on-failure`

Expected: 100% pass.

**Step 4: Review branch**

Inspect `git status`, `git diff main...HEAD`, and confirm the original worktree
still has only its pre-existing changes.
