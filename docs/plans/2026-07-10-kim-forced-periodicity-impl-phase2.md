# KIM Forced-Periodicity — Phase 2 (Periodization + Assembly + Solve) Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (or
> subagent-driven-development) to implement this plan task-by-task. Each code task is TDD.

**Goal:** Turn the Phase-1 fused FP kernel into a working forced-periodicity electrostatic
solver: periodize the background over one period `L = 2(Δr + Δr_tr)`, assemble the dense
Fourier matrix `K_{m,m'} = (2π/L) ∫ G dr_g`, solve `(D + 4π K^{ρΦ})Φ = −4π K^{ρB}Br` with
LAPACK `zgesv`, reconstruct `δΦ(r)`, and expose it through a new `electrostatic_periodic`
run-type.

**Architecture:** Reuse KIM's equilibrium/background machinery on a **redirected grid**: build
an equidistant window grid `[rm−L/2, rm+L/2]`, populate the plasma primitives (n, Te, Ti, Eᵣ, q)
with their *periodized* values (via the vendored `make_periodic`), then re-run
`calculate_equil`/`set_plasma_quantities`/`calculate_thermodynamic_forces_and_susc` so every
derived quantity (`B0, kp, om_E, ρ_L, λ_D, I00…`) is periodic because it is *derived* from the
periodic primitives. The Phase-1 kernel `hatG_rho_phi`/`hatG_rho_B` is reused unchanged.

**Tech Stack:** Fortran (gfortran), LAPACK `zgesv`, the vendored `make_periodic`/`localizer`
(from `~/1_projects/KIM_forced_periodicity/FORCED_PERIODICITY/SRC/`), KIM `plasma_t`/`grid_type`,
KIM custom `error stop` CTest harness. Build/test as in Phase 1.

---

## DESIGN DECISIONS (baked in; flag at review if any should change)

1. **Redirect `rg_grid`, don't refactor.** `calculate_thermodynamic_forces_and_susc` is
   hardcoded to the global `rg_grid` (`species_mod.f90:326-445`). The periodic run-type is a
   standalone run, so it **owns** `rg_grid`: it rebuilds `rg_grid` as the equidistant window and
   runs the recompute chain on it. (Alternative — parameterize that routine by grid — is more
   invasive and deferred.)
2. **Periodize primitives, derive the rest.** Only n, Te, Ti, Eᵣ, q are periodized (via
   `make_periodic`); `B0, kp, om_E, ρ_L, λ_D`, and the susceptibility functions follow from
   KIM's existing recompute chain, so they are automatically periodic and mutually consistent.
   Gradients `dndr, dTdr, dqdr` are recomputed on the window from the periodized profiles.
3. **Direct `Φ`, no R1.** Periodization makes the solution genuinely `L`-periodic, so the
   periodic-basis edge mismatch that motivated the R1 deviation formulation (design §3.5) does
   not arise. Phase 2 solves for the full `Φ` directly. (R1 stays deferred unless a measured
   Gibbs problem reappears.)
4. **Expose `δΦ` via `kim_results_t`.** Add a `complex(dp), allocatable :: Phi(:)` field to
   `kim_results_t` (`kim_solver_m.f90:39-53`) and copy it in `copy_results_from_globals`
   (`kim_solver_m.f90:256-300`) — needed for cross-check B and golden records.
5. **ρ_L-anchored defaults, memo-ratio** (config knobs, overridable): keep the memo Fig.1 ratio
   `Δr : Δr_tr = 1 : 2` at ρ_L-anchored absolute size — `Δr = 5·ρ_L(rm)`, `Δr_tr = 10·ρ_L(rm)`
   ⇒ `L = 2(Δr+Δr_tr) = 30·ρ_L(rm)`; `k_max = 5/ρ_L(rm)` ⇒
   `M = ceil(k_max·L/(2π)) = ceil(150/2π) ≈ 24` (~49 modes, `2M+1`); `N_rg ≳ 4M ≈ 100` equidistant
   `r_g` points. A trivial ~49×49 dense complex solve — comfortably below the memo's "~1000
   modes, minutes" ceiling. (These are the Phase-3 strategy-A convergence knobs.)
6. **Benchmark case:** `(m_mode,n_mode) = (−6, 2)`, resonant at `q = 3`, `type_br_field = 12`
   (constant `Br`), matching the design's strategy-B case.

---

## Ground truth (from Phase-2 investigation; so the executor has zero re-discovery)

**Periodization source** (`~/1_projects/KIM_forced_periodicity/FORCED_PERIODICITY/SRC/`):
- `make_periodic(aperfuns, nfuns, x_mid, dx_asis, dx_tr, x, funs_per)` — `x_mid=rm`,
  `dx_asis=Δr`, `dx_tr=Δr_tr`, period `L=2(Δr+Δr_tr)`; keeps `|x−rm|≤Δr` true, blends the
  transition via `localizer`. `aperfuns(nfuns,x,funs)` is the true-background callback.
- `localizer(sig, x1, x2, x, weight)` — C∞ transition `exp(-2π/(1-t)·exp(-√2/t))`.
- Reference output `fort.3`: with `aperfuns=(x, x²)`, `rm=0.5, Δr=0.25, Δr_tr=0.5`, over
  `x=-1+0.01i` (i=1..300). E.g. `x=-0.99 → (0.51, 0.2601)`, `x=2.0 → (0.5, 0.25)`. Use as the
  Task-1 regression anchor.

**KIM background machinery** (`KIM/src/background_equilibrium/species_mod.f90`):
- Primitives: `plasma%spec(sp)%n/T/dndr/dTdr`, `plasma%q/Er/r_grid`, `plasma%grid_size`.
- Evaluate a profile at arbitrary radius: `binsrc(r_grid,1,n,x,ir)` + `plag_coeff(4,nder,x,
  r_grid(ibeg:iend),coef)` then `sum(coef(0,:)*f(ibeg:iend))` (`util/plag_coeff.f90`;
  example `calculate_equil.f90:169-180`).
- Recompute chain `set_plasma_quantities(plasma)` (`:250-276`): `calculate_plasma_backs(plasma)`
  (grid-size-aware; `:447-561`) → `interpolate_plasma_backs(plasma, grid(:))` (grid-parameterized;
  `:734-799`) → `calculate_thermodynamic_forces_and_susc(plasma)` (**hardcoded `rg_grid`**;
  `:326-445`, `getIfunc` for `I00…`).
- Grid `grid_type` (`grid/grid_mod.f90:44-63`): `%xb(:), %xc(:), %npts_b, %npts_c, %min_val,
  %max_val`; build equidistant via `rg_grid%grid_init_equidistant(n, rmin, rmax, 'rg')` +
  `grid_generate_equidistant()` (`grid/gengrid.f90:14-15`).
- Resonant surface: `prepare_resonances` sets `kim_resonances_m::r_res` where `q(r_res)=|m/n|`
  (`grid/prepare_resonances.f90`).

**Run-type + solve + results:**
- `kim_t` base (`general/KIM_base.f90`): deferred `init(this)`, `run(this)`.
- Mirror `electrostatic_t` (`electrostatic_poisson/poisson.f90:9-48`): `init` =
  create_output_directories → generate_grids → calculate_equil(.true.) → set_plasma_quantities →
  interpolate_equil(rg_grid%xb). Simpler example: `WKB_dispersion_t` (`dispersion/wkb_dispersion.f90:36-109`).
- Factory `KIM_mod.f90:20-33`: add `case("electrostatic_periodic")`.
- `zgesv`: `call zgesv(n, 1, A, n, ipiv, b, n, info)`, solution overwrites `b`
  (`electromagnetic/electromagnetic_solver.f90:241-250`). LAPACK already linked to KIM_lib.
- `kim_results_t` (`general/kim_solver_m.f90:39-53`) — no `Phi`; copy logic
  `copy_results_from_globals` (`:256-300`). `EBdat%Phi` exists (`electrostatic_poisson/fields_mod.f90:7-29`).
- `Br` drive: `type_br_field=12` → constant `Br=(1,0)`, RHS `=−4π·K^{ρB}·Br`
  (`solve_poisson.f90:98-131`, `fields_mod.f90:35-70`). Mode `(m,n,ω)` in `setup_m`.
- Build registration: `KIM/src/CMakeSources.in`; tests in `KIM/tests/` + `KIM/tests/CMakeLists.txt`.

---

## Task 1: Vendor `make_periodic` + `localizer` into KIM (regression vs reference)

**Files:**
- Create: `KIM/src/background_equilibrium/periodization.f90` (module `periodization_m`).
- Modify: `KIM/src/CMakeSources.in` (add to `KIM_background_equilibrium`).
- Create: `KIM/tests/test_periodization.f90`; register in `KIM/tests/CMakeLists.txt`.

**Step 1 — Port the code.** Translate `localizer.f90` and `make_periodic.f90` into
`periodization_m` with `use KIM_kinds_m, only: dp`, `real(dp)`, and an ABSTRACT INTERFACE for
the callback so `make_periodic` takes a procedure argument:
```fortran
abstract interface
    subroutine aperfuns_i(nfuns, x, funs)
        import :: dp
        integer, intent(in) :: nfuns
        real(dp), intent(in) :: x
        real(dp), intent(out) :: funs(nfuns)
    end subroutine
end interface
! public :: make_periodic, localizer
subroutine make_periodic(aperfuns, nfuns, x_mid, dx_asis, dx_tr, x, funs_per)
    procedure(aperfuns_i) :: aperfuns
    ...  ! body transcribed verbatim from make_periodic.f90 (double precision -> real(dp))
```
Keep `localizer` as a module procedure. Transcribe the math EXACTLY (the `modulo` wrap, the two
`localizer`-weighted blends).

**Step 2 — Write the failing regression test** `test_periodization.f90`: define a module
`test_aperfuns` (or a contained subroutine matching `aperfuns_i`) returning `funs=(x, x**2)`.
For `rm=0.5, Δr=0.25, Δr_tr=0.5`, `x = -1 + 0.01*i` (i=1..300), call `make_periodic` and assert
the results match a hardcoded subset of `fort.3` (embed ~6 reference rows: e.g.
`x=-0.99 → (0.51, 0.2601)`, `x=0.30 → (0.30,0.09)` [as-is region], `x=2.0 → (0.5,0.25)`) to
`< 1e-12`. Also assert periodicity: `make_periodic(...,x) == make_periodic(...,x+L)` for a few x
(`L=1.5`), and that in the as-is region `|x-rm|≤Δr` the result equals `(x, x²)` exactly.

**Step 3 — Run red / Step 4 — port until green / Step 5 — full suite / Step 6 — commit**
`git commit` (explicit `git add` of the 4 files — NOT `git add -A`):
`git commit -m "feat(kim): vendor make_periodic/localizer periodization utility"`.

---

## Task 2: KIM background `aperfuns` + periodic primitive sampling

Provide the callback that evaluates KIM's TRUE primitives at arbitrary radius, and a routine
that samples the PERIODIZED primitives on a grid.

**Files:** Modify `periodization.f90` (or a new `periodic_background.f90` module); test in
`test_periodization.f90` (or a new `test_periodic_background.f90`).

**Step 1 — Failing test.** Using the Phase-1 in-memory plasma setup (profiles → namelist →
`kim_init` → `set_profiles_from_arrays`; see `KIM/tests/test_flr2_fourier_kernel.f90`), obtain a
populated `plasma`. Compute `rm` (call `prepare_resonances`; read `kim_resonances_m::r_res`).
Choose `Δr, Δr_tr` (a few `ρ_L`), build an equidistant window grid over `[rm−L/2, rm+L/2]`
(`N_rg` points). Call the new `sample_periodic_primitives(rm, Δr, Δr_tr, rg, prof_out)` for each
grid point and assert: (a) in `|rg−rm|≤Δr` the sampled n/T/Er/q equal the true-profile
interpolation to `<1e-10`; (b) periodicity `prof(rg) == prof(rg+L)`; (c) all finite.

**Step 2 — Run red. Step 3 — Implement.** `sample_periodic_primitives` builds an `aperfuns`
callback (a module procedure with access to the global `plasma`) that returns
`funs = [n_e, Te, Ti, Er, q]` evaluated at `x` via `binsrc`+`plag_coeff` on `plasma%r_grid`
(handle out-of-range x by clamping to the profile domain — `make_periodic` queries
`x_inperiod ± L`, which may exceed the raw domain, so `aperfuns` must extrapolate/clamp
sensibly; document the choice). Then call `make_periodic(aperfuns, 5, rm, Δr, Δr_tr, rg, funs)`.
**Step 4 — green. Step 5 — commit** (explicit add).

> Subtlety to resolve here: `make_periodic` evaluates `aperfuns` at `x_inperiod−L`, `x_inperiod`,
> `x_inperiod+L`, which can fall outside the raw profile domain `[r_min, r_plas]`. Decide and
> test the extrapolation policy (clamp to endpoint value is the safe default, since the ideal
> region is flat-ish). This is why Task 2 exists separately from Task 3.

---

## Task 3 (CRUX): Periodic-background builder — recompute derived plasma on the window

Redirect `rg_grid` to the window and produce a fully-populated PERIODIC `plasma` whose derived
quantities (`ρ_L, λ_D, B0, kp, om_E, I00…`) the Phase-1 kernel can read.

**Files:** New `periodic_background.f90` (`build_periodic_plasma(rm, Δr, Δr_tr, N_rg)`); test
in `test_periodic_background.f90`.

**Step 1 — Failing test.** After setup + `build_periodic_plasma(...)`, assert on the window
grid: `plasma%spec(0)%rho_L`, `lambda_D`, `I00(:,0)` are finite and > 0 (rho_L, lambda_D);
periodicity of a derived quantity (`rho_L(j)` at `rg` vs `rg+L` — pick indices one period apart);
`kp ≈ 0` at the grid point nearest `rm` (resonance condition); and that in the as-is region the
derived `rho_L` matches the non-periodic `rho_L` (recompute the global solver's value at the same
radius, tol ~1e-6).

**Step 2 — Run red. Step 3 — Implement `build_periodic_plasma`:**
1. `prepare_resonances` → `rm`.
2. `L = 2(Δr+Δr_tr)`; rebuild `rg_grid` as equidistant `[rm−L/2, rm+L/2]`, `N_rg` boundary pts
   (`rg_grid%grid_init_equidistant` + `grid_generate_equidistant`).
3. Overwrite `plasma%r_grid` = window `xb`, `plasma%grid_size = N_rg`; fill
   `plasma%spec(sp)%n/T`, `plasma%Er`, `plasma%q` from `sample_periodic_primitives` (Task 2);
   recompute gradients `dndr/dTdr/dqdr` on the window (finite-difference or the routine KIM
   already uses — TRACE how the raw profiles' derivatives are computed in profile setup and reuse
   it).
4. Run `calculate_equil` (recompute `B0, kp, om_E` from periodized `q`, `btor`, `R0`) →
   `calculate_plasma_backs(plasma)` → `calculate_thermodynamic_forces_and_susc(plasma)` on the
   redirected `rg_grid`.
**Step 4 — green. Step 5 — full suite. Step 6 — commit** (explicit add).

> This is the highest-risk task. If `calculate_equil`/the derivative recompute cannot be cleanly
> re-driven on the redirected grid, STOP and report what blocks it (the fallback is to
> parameterize `calculate_thermodynamic_forces_and_susc` by grid — a bigger refactor). Do a small
> spike first if needed. Do not fake derived values.

---

## Task 4: Discrete k-grid + dense matrix assembly `K_{m,m'}`

**Files:** New `periodic_assembly.f90` (`assemble_periodic_matrices(...) -> Kphi(:,:), KB(:,:)`);
test `test_periodic_assembly.f90`.

**Step 1 — Failing test.** With a periodic plasma built (Task 3), for small `M` assemble
`Kphi`, `KB` and assert: shape `(2M+1)×(2M+1)`; the `(m,m')` element equals a brute-force
reference `(2π/L)·Σ_g hatG_rho_phi(plasma, k_m, k_{m'}, g)·Δrg` computed inline in the test
(non-tautological only if the production routine is structured differently — otherwise this is a
characterization test pinning the assembly); all finite. Check the `m=m'` diagonal is nonzero.

**Step 2 — red. Step 3 — Implement:** `k_m = 2π m / L`, `m=−M..M`; equidistant `r_g` = window
`xb`; `Kphi(m,m') = (2π/L)·Σ_g hatG_rho_phi(plasma, k_m, k_{m'}, g)·Δrg` (trapezoid/midpoint —
equidistant over one period is spectrally exact), likewise `KB` with `hatG_rho_B`. Loop cost
`~(2M+1)²·N_rg`. **Step 4 — green. Step 5 — commit** (explicit add).

---

## Task 5: Dense `zgesv` solve + inverse DFT → `δΦ(r)`

**Files:** New `periodic_solve.f90` (`solve_periodic(Kphi, KB, ...) -> Phi_m(:)`,
`reconstruct_delta_phi(Phi_m, r_out) -> dPhi(:)`); test `test_periodic_solve.f90`.

**Step 1 — Failing test (manufactured):** build a small system with a KNOWN solution — e.g. set
`Kphi = 0` so the operator is just the Fourier-Laplacian `D_m = −k_m²` and a chosen RHS, and
assert `Phi_m` matches `−4π KB Br_m / D_m` (skip `m=0` or regularize). Assert `zgesv` `info==0`
and the reconstructed `δΦ(r) = Σ_m Phi_m e^{i k_m r}` matches `Σ (analytic Phi_m) e^{i k_m r}`
to `1e-10`.

**Step 2 — red. Step 3 — Implement:** form `A = D_m δ_{m,m'} + 4π Kphi` (`D_m = −k_m²` up to the
radial-operator convention — verify sign against `solve_poisson`), `b = −4π KB · Br_m` (for
`type_br_field=12`, `Br_m` is the Fourier coefficient of a constant over the window ⇒ `m=0` only);
`call zgesv(2M+1, 1, A, 2M+1, ipiv, b, 2M+1, info)`; check `info`; `Phi_m = b`; reconstruct on an
output grid. Report the condition number (optional: `zgecon`). **Step 4 — green. Step 5 — commit**
(explicit add).

---

## Task 6: New run-type `electrostatic_periodic_t` + config + results `Phi`

**Files:**
- New `KIM/src/electrostatic_poisson/poisson_periodic.f90` (module `rt_electrostatic_periodic_m`).
- Modify `KIM/src/general/KIM_mod.f90` (factory `case`).
- Modify `KIM/src/setup/config_mod.f90` + `general/read_config.f90` (knobs `dr_asis_scale`,
  `dr_tr_scale`, `kmax_scale`, `n_rg_periodic` — or absolute `dr_asis`, `dr_tr`).
- Modify `KIM/src/general/kim_solver_m.f90` (add `Phi(:)` to `kim_results_t` + copy).
- Modify `KIM/src/CMakeSources.in`.
- Test `KIM/tests/test_kim_solver_periodic.f90`.

**Step 1 — Failing test (end-to-end):** in-memory (6,2) plasma; `type_of_run =
'electrostatic_periodic'`; drive the `kim_solver_t` lifecycle (init/set_profiles/solve/results)
OR `from_kim_factory_get_kim('electrostatic_periodic', kim); kim%init(); kim%run()`. Assert
`results%Phi` (or `EBdat%Phi`) is allocated, finite, non-zero, and localized (|Phi| decays from
the resonant layer toward the window edges).

**Step 2 — red. Step 3 — Implement** `init_electrostatic_periodic` (mirror `electrostatic_t`
init but call `build_periodic_plasma` with the config knobs instead of the global grid setup) and
`run_electrostatic_periodic` (assemble → solve → reconstruct → pack `EBdat%Phi` + copy to
`results%Phi`). Register in factory + CMake; add the `Phi` results field + copy in
`copy_results_from_globals`. **Step 4 — green. Step 5 — full suite. Step 6 — commit** (explicit add).

---

## Task 7: Cross-check B — periodic vs global `electrostatic` on (6,2)

**Files:** Test `KIM/tests/test_periodic_vs_global.f90` (or a golden case under
`test/golden/kim/cases/electrostatic_periodic/`).

**Step 1 — Test:** run the existing `electrostatic` (global) and the new
`electrostatic_periodic` on the SAME (6,2) in-memory case; interpolate both `δΦ(r)` onto a common
resonant-layer grid `|r−rm| ≤ Δr`; report `‖δΦ_periodic − δΦ_global‖ / ‖δΦ_global‖` in L2 and
max-norm. Assert it is below a **loose** informational threshold (e.g. 20%) — this is
approx-vs-approx (design §5.2), so treat a large deviation as a diagnostic to investigate, not a
hard pass/fail initially. Log the number. **Commit** (explicit add).

---

## Phase 2 exit criteria

- `electrostatic_periodic` run-type produces a finite, localized `δΦ(r)` on the (6,2) case.
- `test_periodization`, `test_periodic_background`, `test_periodic_assembly`,
  `test_periodic_solve`, `test_kim_solver_periodic` all green; full suite has no new failures
  (only `test_rhs_balance` pre-existing).
- Cross-check B number logged for (6,2).
- **Checkpoint before Phase 3** (strategy-A `Δr`/`N_rg`/`M` convergence harness + golden record).

## Deferred to later phases

- **Strategy A** (Δr / Δr_tr / N_rg / M convergence → periodization-deformation error bar) +
  frozen golden case — Phase 3.
- **R1 deviation** — only if a Gibbs problem is measured (periodization is expected to remove it).
- **Current-diagnostic kernels** `K^{jΦ}, K^{jB}`; electromagnetic extension; toroidal geometry.
- Wiring `kern_include_ks2` to config; A/B of the two k_s² paths.
