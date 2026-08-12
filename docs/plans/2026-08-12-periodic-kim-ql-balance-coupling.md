# Periodic KIM to QL-Balance Coupling Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Run the configured QL-Balance modes through KIM's forced-periodicity solver, embed the local physical response on the global balance grid with a compact transition, retain drift-kinetic electrons, and consume KIM's integral ion tensor with diagnostic provenance.

**Architecture:** QL-Balance selects `electrostatic_periodic` explicitly through its KIM adapter. KIM returns a bounded local response; the adapter interpolates already-derived physical fields and current only inside the local support, while a shared compact weight applies to fields and the square of that weight applies to transport tensors. The adapter computes the local ion integral tensor through KIM's public harmonic seam and keeps the legacy electron drift-kinetic calculation separate.

**Tech Stack:** Fortran 2018, CMake/CTest, KIM solver API, QL-Balance modules, existing Fortran regression tests.

---

### Task 1: Specify the coupling contract with a failing test

**Files:**
- Modify: `QL-Balance/src/base/periodic_embedding_m.f90`
- Create: `QL-Balance/src/test/test_periodic_response_contract.f90`
- Modify: `QL-Balance/src/test/CMakeLists.txt`

**Step 1: Write the failing test**

Test a local complex field, a constant local tensor, the trusted core, transition support, zero outside support, and the invariant that tensor output is `w*w*D_local`.

**Step 2: Run it to verify it fails**

Run: `cmake --build build-ql290 --target test_periodic_response_contract`

Expected: compile failure because the embedding API is not yet present.

**Step 3: Implement the minimal embedding API**

Add `embed_complex_profile` and `embed_tensor_profile` with bounded linear interpolation of the local response, common field weighting, and `w^2` tensor weighting. Reject invalid grids and non-finite local values.

**Step 4: Run the focused test**

Run: `ctest --test-dir build-ql290 -R test_periodic_response_contract --output-on-failure`

Expected: PASS.

**Step 5: Commit**

```bash
git add QL-Balance/src/base/periodic_embedding_m.f90 QL-Balance/src/test/test_periodic_response_contract.f90 QL-Balance/src/test/CMakeLists.txt
git commit -m "test(QL-Balance): specify periodic response embedding contract"
```

### Task 2: Add explicit periodic KIM selection and local result storage

**Files:**
- Modify: `QL-Balance/src/base/control_mod.f90`
- Modify: `QL-Balance/src/base/read_config.f90`
- Modify: `QL-Balance/src/base/kim_wave_code_adapter.f90`
- Modify: `QL-Balance/src/base/wave_code_data_mod.f90`
- Modify: `QL-Balance/namelists/balance_conf.nml`
- Modify: `QL-Balance/src/test/test_kim_adapter.f90`

**Step 1: Write the failing test**

Add a contract assertion that the coupled KIM path requests `electrostatic_periodic`, not the legacy global electromagnetic run, and that periodic responses have a bounded support flag.

**Step 2: Run it to verify it fails**

Run: `ctest --test-dir build-ql290 -R test_kim_adapter --output-on-failure`

Expected: FAIL because the adapter hard-codes `electromagnetic` and exposes no local-support metadata.

**Step 3: Implement the minimal selection/storage path**

Add `kim_run_type` with a production default of `electrostatic_periodic` for `wave_code='KIM'`, pass it through `kim_solver_t%init`, skip global vacuum continuation for periodic results, and store per-mode local response metadata and local ion tensors in adapter-owned arrays.

**Step 4: Run focused adapter tests**

Run: `ctest --test-dir build-ql290 -R 'test_kim_adapter|test_periodic_embedding' --output-on-failure`

Expected: PASS.

**Step 5: Commit**

```bash
git add QL-Balance/src/base/control_mod.f90 QL-Balance/src/base/read_config.f90 QL-Balance/src/base/kim_wave_code_adapter.f90 QL-Balance/src/base/wave_code_data_mod.f90 QL-Balance/namelists/balance_conf.nml QL-Balance/src/test/test_kim_adapter.f90
git commit -m "feat(QL-Balance): select and retain periodic KIM responses"
```

### Task 3: Add KIM periodic ion-tensor production seam

**Files:**
- Modify: `KIM/src/general/kim_solver_m.f90`
- Modify: `KIM/src/electrostatic_poisson/poisson_periodic.f90`
- Modify: `KIM/src/diagnostics/kim_qldiff_mod.f90`
- Modify: `KIM/src/CMakeSources.in`
- Modify: `QL-Balance/src/base/kim_wave_code_adapter.f90`
- Create: `QL-Balance/src/test/test_periodic_ion_tensor.f90`

**Step 1: Write the failing test**

Verify that a periodic solve exposes finite `(2,2,n)` ion coefficients, that the Bparallel-only input is nonzero when configured, and that zero Bparallel reduces exactly to the Phi/Br tensor.

**Step 2: Run it to verify it fails**

Run: `ctest --test-dir build-ql290 -R test_periodic_ion_tensor --output-on-failure`

Expected: compile or contract failure because `kim_results_t` has no tensor member and periodic KIM currently sets Bparallel to zero.

**Step 3: Implement the minimal production seam**

Add optional periodic Bparallel shape/amplitude inputs with a zero-default benchmark mode, compute the local diagonal radial-wave-number contraction by summing the configured cyclotron harmonics through `calc_dqli_integral_harmonic`, and copy the tensor and Bparallel into `kim_results_t`. Preserve the signed cyclotron frequency, perturbation-frequency detuning, and KIM susceptibility evaluator.

**Step 4: Run the focused tensor test**

Run: `ctest --test-dir build-ql290 -R test_periodic_ion_tensor --output-on-failure`

Expected: PASS.

**Step 5: Commit**

```bash
git add KIM/src/general/kim_solver_m.f90 KIM/src/electrostatic_poisson/poisson_periodic.f90 KIM/src/diagnostics/kim_qldiff_mod.f90 KIM/src/CMakeSources.in QL-Balance/src/base/kim_wave_code_adapter.f90 QL-Balance/src/test/test_periodic_ion_tensor.f90
git commit -m "feat(KIM): expose periodic integral ion tensor"
```

### Task 4: Dispatch species models and embed mode responses

**Files:**
- Modify: `QL-Balance/src/base/get_dql.f90`
- Modify: `QL-Balance/src/base/diff_coeffs.f90`
- Modify: `QL-Balance/src/base/kim_wave_code_adapter.f90`
- Modify: `QL-Balance/src/base/wave_code_data_mod.f90`
- Create: `QL-Balance/src/test/test_periodic_transport_dispatch.f90`
- Modify: `QL-Balance/src/test/CMakeLists.txt`

**Step 1: Write the failing test**

Exercise one stationary mode and assert: electrons use the existing drift-kinetic routine, ions use the KIM integral tensor, physical fields are tapered before global accumulation, the tensor is weighted by `w^2`, and the misalignment components share the same complex weight.

**Step 2: Run it to verify it fails**

Run: `ctest --test-dir build-ql290 -R test_periodic_transport_dispatch --output-on-failure`

Expected: FAIL because `get_dql` currently calls the drift-kinetic routine for both species and accumulates unembedded global arrays.

**Step 3: Implement the minimal dispatch**

For `wave_code='KIM'` and `kim_run_type='electrostatic_periodic'`, use the adapter's embedded fields for the established electron routine and copy the adapter's embedded ion tensor into `dqli11:dqli22`. Keep KiLCA and legacy global KIM behavior unchanged. Add explicit model validation and diagnostic model names.

**Step 4: Run focused and existing transport tests**

Run: `ctest --test-dir build-ql290 -R 'test_periodic_transport_dispatch|test_rhs_balance|test_kim_adapter' --output-on-failure`

Expected: PASS.

**Step 5: Commit**

```bash
git add QL-Balance/src/base/get_dql.f90 QL-Balance/src/base/diff_coeffs.f90 QL-Balance/src/base/kim_wave_code_adapter.f90 QL-Balance/src/base/wave_code_data_mod.f90 QL-Balance/src/test/test_periodic_transport_dispatch.f90 QL-Balance/src/test/CMakeLists.txt
git commit -m "feat(QL-Balance): dispatch periodic electron and ion transport models"
```

### Task 5: Add stationary one-mode integration and diagnostics

**Files:**
- Modify: `QL-Balance/src/base/writeData.f90`
- Modify: `QL-Balance/src/base/diagnostic_mod.f90`
- Create: `QL-Balance/src/test/test_periodic_stationary_mode.f90`
- Modify: `QL-Balance/src/test/CMakeLists.txt`
- Modify: `QL-Balance/README.md`

**Step 1: Write the failing test**

Run a synthetic stationary mode through the embedding path and assert finite localized electron/ion coefficients, zero transport outside support, and diagnostic availability for Br, Bparallel, E, species current, transition metadata, and separate tensors.

**Step 2: Run it to verify it fails**

Run: `ctest --test-dir build-ql290 -R test_periodic_stationary_mode --output-on-failure`

Expected: FAIL because diagnostics do not expose the periodic response contract.

**Step 3: Implement diagnostics and stationary harness**

Add a compact diagnostic record and output fields with model names and transition metadata. Keep HDF5 writes optional and make the test use synthetic arrays when external profile fixtures are unavailable.

**Step 4: Run the complete validation suite**

Run: `cmake --build build-ql290 --parallel 4 && ctest --test-dir build-ql290 --output-on-failure`

Expected: all existing tests plus the new periodic coupling tests pass.

**Step 5: Commit**

```bash
git add QL-Balance/src/base/writeData.f90 QL-Balance/src/base/diagnostic_mod.f90 QL-Balance/src/test/test_periodic_stationary_mode.f90 QL-Balance/src/test/CMakeLists.txt QL-Balance/README.md docs/plans/2026-08-12-periodic-kim-ql-balance-coupling.md
git commit -m "feat(QL-Balance): validate periodic KIM stationary transport"
```
