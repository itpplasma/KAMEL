# Periodic KIM and Integral QL-Balance Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement the approved periodic-KIM local-response contract, the gyrokinetic integral ion tensor including (B_\parallel), and the target-current QL-Balance coupling chain represented by GitHub issues #285--#295.

**Architecture:** Periodic KIM remains a local response service. It returns a typed local window, physical fields, species currents, and the ion integral tensor. QL-Balance owns smooth local-to-global embedding, drift-kinetic electrons, target-current amplitude state, profile evolution, and restart behavior. KiLCA supplies the initial vacuum magnetic-drive shapes.

**Tech Stack:** Fortran 2008/2018, CMake/CTest, HDF5, existing KIM/QL-Balance interfaces, Mathematica symbolic oracle copied from the KIM documentation repository, Python/f90nml configuration, and existing golden-record harness.

---

## Constraints and invariants

- Preserve the dirty main checkout and its unrelated edits; all implementation occurs in the isolated `feature/periodic-kim-ql-balance-integral` worktree.
- Follow the signed mode convention tracked by #231. Do not introduce new absolute-value resonance logic.
- Extract/reconcile the finite-harmonic work tracked by #258 before relying on nonzero cyclotron harmonics.
- Target-current normalization is settled: use the existing `I_par_toroidal` interpretation, including the verified (2\pi) geometry factor and CGS conversion.
- Electrons use the existing drift-kinetic Heyn/Markl coefficients. Ions use the new KIM integral tensor.
- (B^r) and (B_\parallel) are prescribed periodic drives in the first implementation and share one complex shielding amplitude.
- The production local-to-global transition is smooth and compact. Do not differentiate a tapered potential.
- No production code is written before its focused test has failed for the expected missing behavior.
- Every accepted implementation slice gets a focused commit and an issue reference in the commit message where practical.

## Baseline and build setup

### Task 0.1: Record the clean implementation baseline

**Files:**

- Create: `docs/plans/2026-08-11-periodic-kim-ql-balance-integral.md`
- Read: `AGENTS.md`, `CONTEXT-MAP.md`, `KIM/CONTEXT.md`, `QL-Balance/CONTEXT.md`, `KiLCA/CONTEXT.md`

**Steps:**

1. Confirm the implementation worktree is based on `roadmap/periodic-kim-ql-balance-gyrokinetic`.
2. Record the dirty main checkout as out of scope.
3. Confirm no build directory exists in the implementation worktree.
4. Locate the developer's dependency/build instructions and the existing KIM/QL-Balance CTest registration.

**Verification:**

```bash
git status --short --branch
git worktree list
rg -n "KIM_lib|ql-balance_lib|add_fortran_test|test_kim_adapter" KIM QL-Balance
```

Expected: implementation branch is clean before code changes; source/test registration points are identified.

### Task 0.2: Configure a reproducible debug build

**Files:**

- No source changes.
- Build directory: `build-periodic-kim-ql-balance-integral/`

**Steps:**

1. Configure with the repository's documented CMake command and Debug checks enabled.
2. Build the smallest existing KIM and QL-Balance targets.
3. Run the existing KIM periodic tests and `test_kim_adapter` if dependencies are available.
4. Save the exact configure/build/test commands and any pre-existing failures in the implementation notes.

**Verification:**

```bash
cmake -S . -B build-periodic-kim-ql-balance-integral -DCMAKE_BUILD_TYPE=Debug
cmake --build build-periodic-kim-ql-balance-integral --target KIM_lib ql-balance_lib
ctest --test-dir build-periodic-kim-ql-balance-integral --output-on-failure
```

If dependency acquisition is unavailable, stop at the first reproducible external failure and retain the command/output; do not silently claim a baseline pass.

## #285: state-isolated periodic local-response contract

### Task 1.1: Specify the failing result-contract test

**Files:**

- Create: `KIM/tests/test_periodic_result_contract.f90`
- Modify: `KIM/tests/CMakeLists.txt`
- Read/modify target: `KIM/src/general/kim_solver_m.f90`

**Steps:**

1. Add a test that initializes `kim_solver_t` with `run_type='electrostatic_periodic'` and a valid in-memory profile fixture.
2. Request one signed resonant mode.
3. Assert that the result exposes the local radial grid, potential, physical electric fields, radial magnetic field, misalignment component, and species currents with consistent bounds.
4. Run the new test before extending the result type.

**Verification:**

```bash
cmake --build build-periodic-kim-ql-balance-integral --target test_periodic_result_contract
ctest --test-dir build-periodic-kim-ql-balance-integral -R test_periodic_result_contract --output-on-failure
```

Expected RED: the periodic result lacks allocated fields and/or the requested result members.

### Task 1.2: Add typed periodic metadata and (B_\parallel) result storage

**Files:**

- Modify: `KIM/src/general/kim_solver_m.f90`
- Modify: `KIM/src/electrostatic_poisson/fields_mod.f90`

**Steps:**

1. Extend the plain-data result with periodic-window metadata: resonance radius, as-is width, transition width, period, Fourier bounds, and trusted-core bounds.
2. Add `Bparallel` and the misalignment-field component to the result contract.
3. Preserve deep-copy semantics in `solver_results`.
4. Define allocation/empty-result behavior for global and periodic run types explicitly.
5. Re-run the test and confirm the new failure is now in periodic packing rather than the public type.

**Verification:** focused result-contract test passes its shape/allocation assertions while field-value assertions remain RED until packing is implemented.

### Task 1.3: Pack periodic (Phi), physical fields, (B^r), misalignment, and current

**Files:**

- Modify: `KIM/src/electrostatic_poisson/poisson_periodic.f90`
- Modify: `KIM/src/electrostatic_poisson/fields_mod.f90`
- Modify: `KIM/src/electrostatic_poisson/periodic_solve.f90`
- Modify: `KIM/src/general/kim_solver_m.f90`

**Steps:**

1. Preserve the endpoint-exclusive periodic grid semantics already used by `compute_periodic_delta_phi`.
2. Add the spectral radial derivative of the periodic potential and construct physical cylindrical and flux-surface components.
3. Pack prescribed (B^r), zero/explicitly configured (B_\parallel), the misalignment field, total current, and species current into `EBdat` and `kim_results_t`.
4. Ensure the periodic result never reads global-only arrays that are not allocated.
5. Add a single-Fourier-mode analytic derivative fixture and run it RED before implementation, then GREEN after.

**Verification:** `test_periodic_result_contract` passes; the single-mode derivative agrees within the declared Fourier tolerance.

### Task 1.4: Make repeated periodic solves restore global state

**Files:**

- Modify: `KIM/src/general/kim_solver_m.f90`
- Modify: `KIM/src/background_equilibrium/periodic_background.f90`
- Create: `KIM/tests/test_periodic_state_restore.f90`
- Modify: `KIM/tests/CMakeLists.txt`

**Steps:**

1. Write a failing A-B-A test that snapshots global grid/plasma dimensions and representative profile values.
2. Add an injected-failure path if the periodic core exposes one; otherwise use a controlled invalid request.
3. Snapshot and restore all state mutated by periodic background construction around each solve.
4. Confirm mode changes and profile updates trigger the intended recomputation without leaking the local window.
5. Verify A-B-A equality and failure cleanup.

**Verification:**

```bash
ctest --test-dir build-periodic-kim-ql-balance-integral -R 'test_periodic_(result_contract|state_restore)' --output-on-failure
```

### Task 1.5: Verify linear magnetic-drive scaling

**Files:**

- Modify: `KIM/tests/test_periodic_result_contract.f90`

**Steps:**

1. Add a failing test using a nontrivial complex scale (s).
2. Solve with unit (B^r), then with (sB^r).
3. Assert linear scaling of (Phi), (E), misalignment field, and current.
4. Keep tensor scaling for later transport slices.

**Verification:** focused scaling test passes after the periodic response is complete.

## #286: full moments and Phi-only integral-ion D11

### Task 2.1: Add a failing full-moment storage test

**Files:**

- Modify: `KIM/src/background_equilibrium/species_mod.f90`
- Create: `KIM/tests/test_integral_moments.f90`
- Modify: `KIM/tests/CMakeLists.txt`

**Steps:**

1. Test access to every (mathcal I_{\ell\sigma}^{p,q}) for (p,q=0,\ldots,3), each species, radial point, and retained harmonic.
2. Verify the test fails because only selected named moments are retained.
3. Replace the selected-only storage with a bounded dense tensor/accessor while retaining compatibility aliases until callers migrate.
4. Verify moment indexing and conjugation/symmetry identities.

**Verification:** moment test passes for both cell-centred and boundary grids and for at least two harmonic indices.

### Task 2.2: Import and pin the Mathematica transport oracle

**Files:**

- Create: `KIM/mathematica/verify_quasilinear_integral_transport.wl`
- Create: `KIM/mathematica/README.md`
- Create: `KIM/src/transport/generated/README.md`
- Create: `KIM/tests/test_quasilinear_integral_algebra.f90`

**Steps:**

1. Copy the verified symbolic script from the KIM documentation repository.
2. Define a deterministic export format for coefficient-polynomial tables and a generator/content hash.
3. Add a Fortran test with exact/rational and high-precision sample vectors.
4. Run the test before generated tables exist and confirm RED.
5. Generate checked data without requiring Mathematica at runtime.

**Verification:** the test reports a deliberate missing-table failure first, then passes against the generated oracle export.

### Task 2.3: Implement the Phi-only ion D11 contraction

**Files:**

- Create: `KIM/src/transport/quasilinear_field_operators_m.f90`
- Create: `KIM/src/transport/quasilinear_integral_m.f90`
- Modify: `KIM/src/general/kim_solver_m.f90`
- Modify: `KIM/src/electrostatic_poisson/poisson_periodic.f90`
- Modify: `KIM/CMakeSources.in`
- Modify: `KIM/tests/CMakeLists.txt`

**Steps:**

1. Add a failing Phi-only D11 test using synthetic local fields and moments.
2. Implement the Gaussian/Bessel field operator and the (Phi)-(Phi) contraction over retained harmonics and full moments.
3. Return real ion D11 on the local KIM grid.
4. Check finite values and preserve the complex field phase convention.
5. Run the exact oracle and synthetic contraction tests.

**Verification:** Phi-only D11 tests pass and no Mathematica runtime is required.

### Task 2.4: Add the D11 drift-limit test

**Files:**

- Create: `KIM/tests/test_integral_drift_limit.f90`
- Modify: `KIM/tests/CMakeLists.txt`
- Modify: `KIM/src/transport/quasilinear_integral_m.f90`

**Steps:**

1. Write the failing zero-FLR, (ell=0), (B_\parallel=0), local-diagonal test against the old D11 expression.
2. Match collision frequency, thermal speed, charge, phase, and CGS units explicitly.
3. Implement any minimal normalization correction needed by the derived formalism.
4. Record separate algebraic and asymptotic tolerances.

**Verification:** D11 agrees with the old ion drift-kinetic result at the predicted limit.

## #287: complete the Phi/Br ion tensor

### Task 3.1: Add failing ordered-field-pair tests

**Files:**

- Modify: `KIM/tests/test_quasilinear_integral_algebra.f90`
- Create: `KIM/tests/test_quasilinear_integral_properties.f90`

**Steps:**

1. Add independent Phi-Phi, Phi-Br, Br-Phi, and Br-Br fixtures for all four tensor entries.
2. Assert D12/D21 symmetry, phase invariance, and positive semidefiniteness.
3. Run and confirm failure for missing cross channels.

### Task 3.2: Implement and expose the Phi/Br tensor

**Files:**

- Modify: `KIM/src/transport/quasilinear_field_operators_m.f90`
- Modify: `KIM/src/transport/quasilinear_integral_m.f90`
- Modify: `KIM/src/general/kim_solver_m.f90`
- Modify: `KIM/src/electrostatic_poisson/poisson_periodic.f90`

**Steps:**

1. Add all four tensor entries and the four ordered Phi/Br field pairs.
2. Hermitianize the two-wave-number quadratic form before taking its real part.
3. Return `D_ion(1:2,1:2,:)` in the periodic result.
4. Add finite, Onsager, and positive-semidefinite diagnostics with test-only strictness.
5. Run algebra, property, and drift-limit tests.

**Verification:** all Phi/Br tests pass; no B-parallel terms are accidentally included.

## #288: carry B-parallel from KiLCA through KIM

### Task 4.1: Specify the KiLCA-to-KIM B-parallel mapping test

**Files:**

- Create: `QL-Balance/tests/test_kim_bparallel_mapping.f90`
- Modify: `QL-Balance/src/base/kim_wave_code_adapter.f90`
- Modify: `QL-Balance/src/base/wave_code_data_mod.f90`

**Steps:**

1. Add a failing mapping fixture with a known KiLCA `Bp` component, signed mode, phase, units, and resonance.
2. Assert the mapped B-parallel shape is independent of Br and preserves the complex scale.
3. Run RED before adding storage and mapping.

### Task 4.2: Add B-parallel to KIM fields and the periodic request

**Files:**

- Modify: `KIM/src/electrostatic_poisson/fields_mod.f90`
- Modify: `KIM/src/general/kim_solver_m.f90`
- Modify: `KIM/src/electrostatic_poisson/poisson_periodic.f90`
- Modify: `QL-Balance/src/base/kim_wave_code_adapter.f90`
- Modify: `QL-Balance/src/base/wave_code_data_mod.f90`
- Modify: `QL-Balance/src/base/wave_code_data_64bit.f90`

**Steps:**

1. Add independent B-parallel request/result storage and provenance.
2. Extract and normalize the KiLCA vacuum Bp shape at the same resonance as Br.
3. Preserve signed mode and phase conventions at the adapter boundary.
4. Implement charge/current B-parallel response for the supported periodic kernel.
5. Add an explicit benchmark-only B-parallel-zero switch; production cannot silently zero it.

**Verification:** mapping, B-parallel-only, zero-B-parallel, and complex scaling tests pass.

## #289: complete B-parallel transport channels

### Task 5.1: Add failing B-parallel tensor-channel tests

**Files:**

- Modify: `KIM/tests/test_quasilinear_integral_algebra.f90`
- Modify: `KIM/tests/test_quasilinear_integral_properties.f90`

**Steps:**

1. Add Bparallel-Bparallel and both ordered cross pairs with Phi and Br.
2. Assert that zero B-parallel removes exactly the five associated contributions.
3. Add the analytic long-wavelength magnetic-moment fixture.
4. Run RED before adding the channels.

### Task 5.2: Implement all nine field-pair contractions

**Files:**

- Modify: `KIM/src/transport/quasilinear_field_operators_m.f90`
- Modify: `KIM/src/transport/quasilinear_integral_m.f90`
- Modify: `KIM/src/general/kim_solver_m.f90`

**Steps:**

1. Add the five B-parallel-containing ordered pairs to every D entry.
2. Match the long-wavelength (v_m^r) phase and species-sign convention.
3. Generate and validate every coefficient-polynomial family against Mathematica.
4. Re-run Hermiticity, Onsager, positivity, phase, harmonic, and finite-FLR tests.

**Verification:** all nine ordered field pairs pass independent and combined tests.

## #290: couple one periodic mode into QL-Balance

### Task 6.1: Add failing local-to-global embedding tests

**Files:**

- Create: `QL-Balance/tests/test_local_response_embedding.f90`
- Modify: `QL-Balance/src/test/CMakeLists.txt`
- Create: `QL-Balance/src/base/local_response_embedding_mod.f90`

**Steps:**

1. Test a compact smooth weight that is one in the as-is core and zero outside support.
2. Assert endpoint smoothness and (D_\mathrm{global}=w^2D_\mathrm{local}).
3. Add a fixture proving no derivative-of-window electric field is created.
4. Add a fixture proving misalignment components retain a common transition and phase.
5. Run RED before implementing the embedding helper.

### Task 6.2: Replace the global-KIM single-mode path

**Files:**

- Modify: `QL-Balance/src/base/kim_wave_code_adapter.f90`
- Modify: `QL-Balance/src/base/get_dql.f90`
- Modify: `QL-Balance/src/base/wave_code_data_mod.f90`
- Modify: `QL-Balance/src/base/wave_code_data_64bit.f90`
- Modify: `QL-Balance/src/base/diff_coeffs.f90`
- Modify: `QL-Balance/src/base/local_response_embedding_mod.f90`

**Steps:**

1. Select `electrostatic_periodic` for the coupled KIM mode.
2. Return and store local window fields/tensors instead of interpolating global-only arrays.
3. Apply the compact physical-field transition.
4. Call the legacy drift-kinetic routine for electrons and KIM's integral tensor for ions.
5. Add explicit electron/ion model validation; remove hard-coded common-model branches.
6. Keep hard-zero/extrapolate modes diagnostic-only.

**Verification:** one stationary periodic mode produces finite localized electron/ion tensors and no transport outside support.

## #291: drift-kinetic limiting benchmark

### Task 7.1: Add QL-Balance dual-evaluation output

**Files:**

- Modify: `QL-Balance/src/base/get_dql.f90`
- Modify: `QL-Balance/src/base/writeData.f90`
- Modify: `QL-Balance/src/base/wave_code_data_mod.f90`
- Create: `QL-Balance/tests/test_integral_drift_limit.f90`

**Steps:**

1. Add a benchmark switch that computes old and new ion tensors but evolves with the selected model.
2. Record all four absolute and relative residuals.
3. Keep electron output and evolution on the old drift-kinetic path.
4. Add the zero-FLR/harmonic/B-parallel/local-diagonal fixture.
5. Add a constant-profile periodic-window fixture.

**Verification:** benchmark residuals converge at separate algebraic, DFT, and asymptotic tolerances; electron tensors are unchanged.

### Task 7.2: Add the end-to-end limiting golden case

**Files:**

- Create: `test/golden/cases/kim_periodic_integral_transport_limit/`
- Modify: `test/golden/README.md` as required

**Steps:**

1. Add a minimal periodic KIM/QL-Balance input with B-parallel disabled explicitly.
2. Run old/new ion dual evaluation and freeze the residual output.
3. Ensure the golden case is independent of Mathematica at runtime.

**Verification:** golden harness passes and reports all four tensor residuals.

## #292: target-current normalization

### Task 8.1: Add failing amplitude-scaling unit tests

**Files:**

- Create: `QL-Balance/src/test/test_shielding_amplitude.f90`
- Modify: `QL-Balance/src/test/CMakeLists.txt`
- Modify: `QL-Balance/src/base/transp_coeffs_mod.f90`

**Steps:**

1. Add unit-response, target-current, double-current, current-floor, and double-scaling fixtures.
2. Verify RED because only the squared `antenna_factor` exists.
3. Add an explicit complex per-mode amplitude state and linear/quadratic scaling helpers.
4. Preserve `antenna_factor=abs(s)^2` only as compatibility output and prevent duplicate scaling.

**Verification:** amplitude tests pass without changing the legacy path when periodic integral coupling is disabled.

### Task 8.2: Integrate the first-step constant-psi target-current update

**Files:**

- Modify: `QL-Balance/src/base/time_evolution.f90`
- Modify: `QL-Balance/src/base/SingleStep.f90`
- Modify: `QL-Balance/src/base/paramscan.f90`
- Modify: `QL-Balance/src/base/transp_coeffs_mod.f90`
- Modify: `QL-Balance/src/base/wave_code_data_mod.f90`

**Steps:**

1. Initialize a unit constant-psi Br drive and the common normalized B-parallel shape.
2. Integrate shielding current over the trusted local layer using verified (2\pi)/CGS conventions.
3. Update the amplitude with target-current normalization, relaxation, current floor, finite checks, and scale bounds.
4. Scale fields/current by s and tensors by abs(s)^2 exactly once.
5. Emit target, calculated current, residual, amplitude, and guard diagnostics.

**Verification:** one-step target-current and scaling tests pass through SingleStep, TimeEvolution, and ParameterScan entry points.

## #293: multi-mode profile feedback

### Task 9.1: Add failing mode-order and profile-refresh tests

**Files:**

- Modify: `QL-Balance/src/test/test_kim_adapter.f90`
- Create: `QL-Balance/src/test/test_periodic_multimode_feedback.f90`
- Modify: `QL-Balance/src/test/CMakeLists.txt`

**Steps:**

1. Add two signed modes with separated and overlapping supports.
2. Assert order-independent resonance locations and accumulated tensors.
3. Assert updated QL profiles are passed before each periodic response refresh.
4. Assert KIM state is not reused across modes incorrectly.
5. Run RED before adapter multi-mode changes.

### Task 9.2: Implement mode-resolved periodic feedback

**Files:**

- Modify: `QL-Balance/src/base/kim_wave_code_adapter.f90`
- Modify: `QL-Balance/src/base/get_dql.f90`
- Modify: `QL-Balance/src/base/gengrid.f90` only where needed to preserve signed mode tuples
- Modify: `QL-Balance/src/base/wave_code_data_mod.f90`

**Steps:**

1. Refresh KIM from accepted global profiles for every mode and time state.
2. Keep per-mode local fields, windows, currents, tensors, and amplitudes separate.
3. Accumulate global diffusion incoherently after local embedding.
4. Preserve nonresonant-mode diagnostics and signed mode provenance.

**Verification:** separated/overlapping two-mode tests pass and close the periodic portion of #18.

## #294: time evolution and restart

### Task 10.1: Add failing accepted/trial/restart tests

**Files:**

- Create: `QL-Balance/src/test/test_shielding_restart.f90`
- Modify: `QL-Balance/src/test/CMakeLists.txt`
- Modify: restart/checkpoint implementation identified by `rg -n "restart|checkpoint|h5.*restart" QL-Balance/src`

**Steps:**

1. Test accepted and rejected steps with distinct profile and amplitude states.
2. Test uninterrupted versus restarted continuation.
3. Test nonfinite/current-floor/scale-bound rollback.
4. Run RED before persistence and rollback changes.

### Task 10.2: Integrate amplitude into accepted time evolution

**Files:**

- Modify: `QL-Balance/src/base/time_evolution.f90`
- Modify: `QL-Balance/src/base/SingleStep.f90`
- Modify: `QL-Balance/src/base/paramscan.f90`
- Modify: restart/checkpoint implementation
- Modify: `QL-Balance/src/base/writeData.f90`

**Steps:**

1. Use accepted profiles to calculate unit response/current, update amplitudes, embed fields/tensors, and advance profiles.
2. Commit amplitudes only with accepted time steps.
3. Restore amplitudes and profiles on rejection.
4. Persist phase policy, target/current residuals, relaxation, and normalization version.
5. Ensure SingleStep, ParameterScan, and TimeEvolution share one normalization path.

**Verification:** multi-step, rejection, and restart tests pass deterministically.

## #295: production workflow and validation

### Task 11.1: Expose validated configuration

**Files:**

- Modify: `KIM/src/setup/config_mod.f90`
- Modify: KIM configuration reader and `KIM/nmls/KIM_config.nml`
- Modify: QL-Balance control/read configuration modules
- Modify: `QL-Balance/namelists/balance_conf.nml`
- Modify: `python/balance_interface/balance_conf.py`

**Steps:**

1. Add explicit periodic run type, transport models, harmonic bounds, B-parallel source/zero switch, transition controls, benchmark switch, target-current guards, and restart settings.
2. Reject unsupported combinations during configuration validation.
3. Add minimal configuration fixtures for periodic production and drift-limit benchmark.

**Verification:** Python and Fortran configuration tests reject invalid combinations and accept the two supported golden configurations.

### Task 11.2: Add output provenance and diagnostics

**Files:**

- Modify: `QL-Balance/src/base/writeData.f90`
- Modify: KIM periodic HDF5 output routines
- Modify: relevant `Data.md`/README documentation

**Steps:**

1. Write local-window/transition metadata, Br, B-parallel, physical E, species currents, unscaled/scaled tensors, mode amplitudes, and target-current history.
2. Record model names, Fourier convention, units, harmonic cutoffs, transition policy, and Mathematica-table hash.
3. Add output readback assertions to the golden harness.

**Verification:** output is sufficient to reproduce normalization, local embedding, and benchmark residuals.

### Task 11.3: Add production golden, convergence, and performance validation

**Files:**

- Create: `test/golden/cases/kim_periodic_integral_transport/`
- Modify: `test/golden/README.md`
- Create/modify validation scripts under `test/ql-balance/`

**Steps:**

1. Add the coupled periodic production golden case with B-parallel, target-current evolution, profile feedback, and restart.
2. Add Fourier/radial/harmonic/FLR/time-step convergence scans.
3. Add transition-width sensitivity and central-layer periodic/global comparison.
4. Measure wall time and memory against global KIM for representative workloads.
5. Run existing KiLCA, global-KIM, and default QL-Balance regressions.
6. Produce a human-review report with tolerances, residuals, sensitivity, convergence, and performance.

**Verification:** all automated tests and golden cases pass; #295 receives the documented scientific sign-off separately from implementation completion.

## Commit and issue sequence

Use focused commits in this order, preserving green tests after each:

1. `feat(KIM): expose periodic local response` (#285)
2. `feat(KIM): add integral ion Phi D11` (#286)
3. `feat(KIM): complete Phi Br ion tensor` (#287)
4. `feat(KIM): carry Bparallel through periodic response` (#288)
5. `feat(KIM): complete Bparallel transport channels` (#289)
6. `feat(QL): couple periodic KIM local response` (#290)
7. `test: benchmark integral ions in drift limit` (#291)
8. `feat(QL): normalize periodic response to target current` (#292)
9. `feat(QL): add periodic multimode feedback` (#293)
10. `feat(QL): persist shielding amplitudes through time evolution` (#294)
11. `test: validate production periodic integral workflow` (#295)

Before claiming completion, run the verification-before-completion checklist: inspect the final diff, run focused and full relevant CTests, run golden comparisons, verify clean status, and report any unavailable external dependency instead of masking it.
