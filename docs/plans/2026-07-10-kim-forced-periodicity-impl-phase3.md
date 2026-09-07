# KIM Forced-Periodicity — Phase 3 (Strategy-A Convergence + Golden) Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (or
> subagent-driven-development) to implement this plan task-by-task. Each code task is TDD.

**Goal:** Validate the `electrostatic_periodic` solver on its OWN terms — three-knob
convergence (`N_rg`, `M`, `Δr`) reporting the self-consistent periodization-deformation error
bar — and freeze a golden `(6,2)` case so the numbers are pinned in CI.

**Architecture:** Refactor the run-type's core (build → assemble → solve → reconstruct) into one
reusable routine `compute_periodic_delta_phi(rm, dx_asis, dx_tr, M, n_rg, r_out) -> δΦ` that both
the run-type and the convergence harness call. Convergence tests evaluate `δΦ` on a FIXED
diagnostic resonant-layer grid across parameter refinements and report the Cauchy relative
residual `‖δΦ_i − δΦ_{i−1}‖ / ‖δΦ_i‖` (L2 + max). `N_rg` (quadrature) and `M`/`k_max` (spectral)
must converge toward zero; the `Δr` scan's plateau residual IS the reported periodization error.

**Tech Stack:** Fortran, the Phase-2 modules (`periodic_background_m`, `periodic_assembly_m`,
`periodic_solve_m`, `rt_electrostatic_periodic_m`), KIM `error stop` CTest harness, the golden
harness under `test/golden/kim/`.

**Baseline:** all KIM ctest pass except the pre-existing `test_rhs_balance`.

---

## Open item carried into Phase 3 (do NOT try to fix here)

The strategy-B cross-check (periodic vs global) shows **~86% relative difference** on (6,2)
(`max|Φ_global|=0.121` vs `max|Φ_periodic|=0.027`, ~4.5× magnitude gap; DC-removed larger).
Leading suspects: **assembly/solve normalization** (the dominant ~4.5× factor), the **`k_s²`
convention** (`kern_include_ks2` default `.false.`), and the **`D_m=−k_m²`** operator. Phase 3
deliberately validates the periodic solver's *internal numerical convergence* independent of the
global reference; the cross-check magnitude is a separate physics investigation (recorded in
memory `kim-global-electrostatic-umfpack-bug` follow-on note). Strategy-A convergence is still
meaningful even if an overall normalization factor is later corrected (a constant scale cancels in
the *relative* Cauchy residual).

---

## Ground truth (Phase-2 pieces to reuse)

- `periodic_background_m::build_periodic_plasma(rm, dx_asis, dx_tr, n_rg)` — redirects `rg_grid`
  to the window, populates the periodic derived plasma (with the `u0_seed` B0 fix).
- `periodic_assembly_m::assemble_periodic_matrices(plasma, L, M, Kphi, KB)`.
- `periodic_solve_m::solve_periodic(Kphi, KB, L, M, Br_const, Phi_m, info)` and
  `reconstruct_delta_phi(Phi_m, L, M, r_out)`.
- `rt_electrostatic_periodic_m::run_electrostatic_periodic` (`KIM/src/electrostatic_poisson/poisson_periodic.f90`)
  currently inlines: `prepare_resonances` → rm → interp `rho_L(rm)` → (`dx_asis, dx_tr, L,
  k_max, M = ceil(k_max·L/2π), n_rg`) → build → assemble → solve → reconstruct on `rg_grid%xb` →
  pack `EBdat%Phi`.
- `Br_const = cmplx(Br_boundary_re, Br_boundary_im, dp)` (default (1,0)).
- Golden harness: `test/golden/kim/cases/<case>/` with `KIM_config.nml` + profile `.dat`; runner
  auto-discovers via `GR_INPUT=KIM_config.nml`; compared numerically (`gr_numcompare.py`, rtol
  1e-7/atol 1e-12); runs in the GitHub `golden` job, not `make test`. Existing case:
  `test/golden/kim/cases/electromagnetic/`.
- In-memory (6,2) test setup pattern: `KIM/tests/test_kim_solver_periodic.f90` /
  `test_periodic_assembly.f90` (profiles → namelist `m_mode=-6,n_mode=2,type_br_field=12,
  collision_model='FokkerPlanck'` → `kim_init` → `set_profiles_from_arrays` → factory → `%init()`).

---

## Task 1: Extract a reusable `compute_periodic_delta_phi` core (refactor, no behavior change)

Pull the run-type's build→assemble→solve→reconstruct into one routine taking EXPLICIT parameters,
so both the run-type and the convergence harness share it.

**Files:**
- Modify: `KIM/src/electrostatic_poisson/poisson_periodic.f90` — add public
  `compute_periodic_delta_phi(rm, dx_asis, dx_tr, M, n_rg, Br_const, r_out, dPhi, info)` and make
  `run_electrostatic_periodic` call it.
- Test: `KIM/tests/test_kim_solver_periodic.f90` (existing) must still pass unchanged.

**Step 1 — Write the failing test** (in a NEW `KIM/tests/test_periodic_convergence.f90`, which
Task 2 extends): after the (6,2) setup + `%init()`, compute `rm` (`prepare_resonances`),
`rhoL_rm` (interp), pick `dx_asis=5*rhoL_rm, dx_tr=10*rhoL_rm, M=24, n_rg=96`, a diagnostic grid
`r_out` of 41 points over `[rm-3*rhoL_rm, rm+3*rhoL_rm]`, then
`call compute_periodic_delta_phi(rm, dx_asis, dx_tr, 24, 96, (1.0_dp,0.0_dp), r_out, dPhi, info)`
and assert `info==0`, `dPhi` finite, `maxval(abs(dPhi)) > 0`. `error stop` on failure.

**Step 2 — Run red** (routine absent): `cmake -S . -B build && cmake --build build --target
test_periodic_convergence.x && ctest --test-dir build -R test_periodic_convergence --output-on-failure`.

**Step 3 — Implement** `compute_periodic_delta_phi`: body = the current `run_electrostatic_periodic`
core (`build_periodic_plasma(rm,dx_asis,dx_tr,n_rg)` → `L=2*(dx_asis+dx_tr)` →
`assemble_periodic_matrices(plasma,L,M,Kphi,KB)` → `solve_periodic(Kphi,KB,L,M,Br_const,Phi_m,info)`
→ `dPhi = reconstruct_delta_phi(Phi_m,L,M,r_out)`), returning `dPhi(size(r_out))` and `info`.
Refactor `run_electrostatic_periodic` to compute its params then call this with
`r_out = rg_grid%xb` and pack the result into `EBdat%Phi` (behavior preserved — the existing
end-to-end test must still pass).

**Step 4 — Run green** + full suite (only `test_rhs_balance` fails).

**Step 5 — Commit** (explicit `git add`; NEVER `git add -A`):
`git commit -m "refactor(kim): extract compute_periodic_delta_phi core for convergence harness"`.

---

## Task 2: `N_rg` (quadrature) convergence test

**Files:** Modify `KIM/tests/test_periodic_convergence.f90` (add `test_nrg_convergence`).

**Step 1 — Write the test.** Fix `M=32` and `dx_asis=5*rhoL_rm, dx_tr=10*rhoL_rm`; fixed
diagnostic grid `r_out` (41 pts over the resonant layer). For `n_rg ∈ {48, 96, 192, 384}` call
`compute_periodic_delta_phi(...)` and collect `dPhi_k`. Compute the Cauchy relative residual
`res_k = ‖dPhi_k − dPhi_{k-1}‖_2 / ‖dPhi_k‖_2` for k=2..4 (and the max-norm version). PRINT all
residuals. ASSERT: the residual sequence is DECREASING and the finest residual `res_4 < 1e-2`
(quadrature converges as `N_rg` grows). `error stop` if not monotonically converging or if
`res_4` exceeds the threshold. (State the threshold; loosen only with justification if the
periodization non-periodicity floor — the ~2.7e-5 `rho_L` deformation — sets a higher plateau.)

**Step 2 — red / 3 — (already implemented core; run) / 4 — green** (this is a validation test; it
may pass on first run). If it does NOT converge, that is a real finding — investigate the
quadrature (seam handling, equidistant weight), don't loosen blindly.

**Step 5 — Commit** `git commit -m "test(kim): N_rg quadrature convergence of periodic solver"`.

---

## Task 3: `M`/`k_max` (spectral) convergence test

**Files:** Modify `KIM/tests/test_periodic_convergence.f90` (add `test_M_convergence`).

**Step 1 — Write the test.** Fix `n_rg=192` (converged from Task 2) and the window
`dx_asis/dx_tr`; fixed `r_out`. For `M ∈ {8, 16, 24, 32}` call `compute_periodic_delta_phi(...)`;
compute the Cauchy residuals `res_k` (L2 + max). PRINT them. ASSERT the residual decreases and
`res(last) < 1e-2` (the spectral basis converges as `M` grows). `error stop` on non-convergence.

**Step 2-4** as Task 2 (validation; investigate real non-convergence). **Step 5 — Commit**
`git commit -m "test(kim): M spectral convergence of periodic solver"`.

---

## Task 4: `Δr` deformation scan (the reported error bar)

**Files:** Modify `KIM/tests/test_periodic_convergence.f90` (add `test_dr_deformation_scan`).

**Step 1 — Write the test.** Fix `M=32`, `n_rg=192` (both converged). Keep the transition ratio
`dx_tr = 2*dx_asis` (memo Fig.1). Scan `dx_asis ∈ {3, 5, 8, 12} * rhoL_rm` (so the window grows).
Evaluate `δΦ` on a FIXED diagnostic resonant-layer grid `r_out` = 41 pts over
`[rm - 3*rhoL_rm, rm + 3*rhoL_rm]` (the SAME physical grid for every `dx_asis`, i.e. the innermost
layer common to all windows). Compute the Cauchy residual `res_k` between successive `dx_asis`.
PRINT them prominently as the **periodization-deformation error bar** (this is the reported
number, design §5.1). ASSERT only that all `δΦ` are finite and `res_k` is finite and does NOT
GROW without bound (a converging/plateauing sequence is the physical expectation; a growing
sequence is a finding to investigate). This task's output is a NUMBER (the plateau residual at the
largest affordable `Δr`), logged — soft gate on finiteness, not a tight threshold.

**Step 2-4** (validation). **Step 5 — Commit**
`git commit -m "test(kim): Dr periodization-deformation error scan"`.

---

## Task 5: Freeze a golden `(6,2)` periodic case

Pin the periodic solver's numbers against regression via the golden harness (CI-only).

**Files:**
- Create: `test/golden/kim/cases/electrostatic_periodic/KIM_config.nml` (mirror the
  `electromagnetic` golden case's config, but `type_of_run='electrostatic_periodic'`,
  `collision_model='FokkerPlanck'`, `m_mode=-6, n_mode=2, type_br_field=12`, and the
  `&KIM_PERIODIC` scales `periodic_dr_asis_scale=5, periodic_dr_tr_scale=10, periodic_kmax_scale=5,
  periodic_n_rg=96`).
- Create: the profile `.dat` files the case needs (`n.dat, q.dat, Te.dat, Ti.dat` — copy/adapt
  from `test/golden/kim/cases/electromagnetic/`, ensuring the profile makes q cross 3 for the
  (6,2) resonance).

**Step 1 — Author the case.** Copy `test/golden/kim/cases/electromagnetic/` to
`.../electrostatic_periodic/`; edit `KIM_config.nml` per above; keep/adjust the profiles so the
run resolves the (6,2) resonance and produces the localized `δΦ`.

**Step 2 — Generate the baseline.** Follow `test/golden/kim/config.sh` / the golden README to
build and run the case with the runner, producing the reference output (the `δΦ`/fields the
comparator will pin). Confirm the run exits 0 and writes a nonzero `Phi`. (The golden runner uses
`GR_EXE=KIM.x`, `GR_INPUT=KIM_config.nml`.)

**Step 3 — Verify it re-compares clean.** Re-run the golden comparison (`gr_numcompare.py`, rtol
1e-7/atol 1e-12) against the just-frozen baseline → identical. Confirm the case is discovered by
`gr_dispatch.sh`.

**Step 4 — Commit** the case + baseline (explicit `git add` of the case dir):
`git commit -m "test(golden): freeze electrostatic_periodic (6,2) case"`. NOTE: golden `.dat`
baselines are large/numeric — stage ONLY the new case directory; do not sweep other untracked
golden trees (e.g. `test/ql-balance/`).

---

## Phase 3 exit criteria

- `test_periodic_convergence` green: `N_rg` and `M` residuals converge (< 1e-2); the `Δr`
  deformation scan reports a finite, plateauing error bar (the periodization error number).
- A frozen `electrostatic_periodic` golden case in `test/golden/kim/cases/`, re-comparing clean.
- Full suite unchanged (only `test_rhs_balance` fails).
- **The periodic solver is self-validated** (internal numerical convergence + regression-pinned),
  with the periodization-deformation error quantified.

## Deferred / open

- **Strategy-B cross-check discrepancy (~86%)** — the periodic-vs-global magnitude/shape mismatch;
  investigate normalization (~4.5× factor), `kern_include_ks2`, and `D_m` separately.
- R1 deviation formulation; current-diagnostic kernels; electromagnetic + toroidal extensions;
  wiring `kern_include_ks2` to config.
