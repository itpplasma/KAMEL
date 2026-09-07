# KIM Forced-Periodicity — Phase 1 (Fused FP Fourier Kernel) Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan
> task-by-task. Each code task is TDD (red → green → commit).

**Goal:** Deliver a tested Fortran function giving the fused Fokker–Planck off-diagonal
Fourier-space kernel integrand `Ĝ^{ρΦ}(k_r, k'_r, r_g)` (and `Ĝ^{ρB}`), which collapses
**exactly** to the existing diagonal source-of-truth `calc_hatK_Phi_in_Fourier` when
`k_r = k'_r`. Pure physics; no solver, no run-type yet (those are Phase 2+).

**Architecture:** Extract the per-species diagonal integrand from
`calc_hatK_Phi_in_Fourier` into a shared, reusable function, then generalize it to two
wavenumbers by (a) replacing the FLR argument `b = k_r² ρ_L²` with the symmetric pair
`b₊ = ρ_L²(k_r²+k'_r²)/2`, `b× = ρ_L²|k_r k'_r|`, and (b) multiplying by the Fourier phase
`exp(i(k_r−k'_r) r_g)`. Because the generalized form reduces to the diagonal expression
term-by-term at `k_r=k'_r` (phase → 1, `b₊=b×=b`), collapse is exact **by construction**;
the collapse test is a regression guard. Everything reuses the KIM `plasma_t`/`species_t`
backgrounds sampled on the guiding-centre grid `rg_grid%xb`.

**Tech Stack:** Fortran (gfortran 14, `-O3 -fcheck=all`), `fortnum_special::bessel_in`
(unscaled modified Bessel Iₙ), KIM custom `error stop`/`stop 1` unit-test harness under
`KIM/tests/`, CTest. Build: `cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release &&
cmake --build build`; test: `ctest --test-dir build --output-on-failure`.

---

## DECISION (ks² convention — RESOLVED 2026-07-10: keep k_s² behind a toggle)

The diagonal source-of-truth `calc_hatK_Phi_in_Fourier`
(`KIM/src/asymptotics/FLR2_asymptotics.f90:195-199`) deliberately uses the FLR argument

```
b = k_r² · ρ_L²          ! k_s² DROPPED — comment: "k_s diverges at r=0"
```

whereas the design doc (§3.1) and the deprecated Krook off-diagonal integrand
(`bfa9d70c:src/deprecated/integrands.f90`) use the **full perpendicular** argument
(`k_s²` KEPT). These disagree off-diagonal, and at the diagonal (`k_r=k'_r`) they agree
**only if `k_s=0`**.

**RESOLUTION (author decision): implement BOTH behind a toggle** `kern_include_ks2`
(module variable in `flr2_fourier_kernel_m`, `logical, save, public`, **default `.false.`**):

```
kern_include_ks2 = .false.  (DEFAULT — drop k_s²; exact collapse to calc_hatK_Phi_in_Fourier)
    b₊ = ρ_L² (k_r² + k'_r²)/2 ,          b× = ρ_L² |k_r k'_r|

kern_include_ks2 = .true.   (design §3.1 full form)
    b₊ = ρ_L² (2k_s² + k_r² + k'_r²)/2 ,  b× = ρ_L² √(k_s²+k_r²)·√(k_s²+k'_r²)
```

Both forms are expressed in terms of the diagonal's `ρ_L²` (which carries
`ion_flr_scale_factor`, default `1.0`) — so the `k_s²` term is also scaled by that factor;
acceptable and default-neutral. `k_s = plasma%ks(j)`.

**Testing consequence:**
- Default path (`.false.`): **exact** diagonal collapse to `calc_hatK_Phi_in_Fourier`
  (Task 3). This is the primary R5 guard.
- Toggle path (`.true.`): collapse cannot match the ks-dropped diagonal; instead assert (new
  Task 3b) that at `k_r=k'_r` it equals the diagonal FP expression fed `b_full=(k_s²+k_r²)ρ_L²`
  (a structural exact check via the shared `core_*_sp(bplus,bcross,j)` helper — no need to
  manufacture `k_s=0`).

---

## Reconciliation with GitHub issue #175 (RESOLVED 2026-07-10)

Issue #175 ("revive Fourier-space solver via enforced radial periodicity") covers the same
goal from the Kasilov memo's angle. Two forks were decided by the author:

- **Kernel flavor → stay Fokker–Planck** (design decision #4), NOT restore the issue's Krook
  kernel. The recoverable `1ae0aeeb~1:KIM/src/kernels/kernel_mod.f90` (module `kernels_m`,
  `kernel_rho_phi_of_kr_krp_rg`/`kernel_rho_B_of_kr_krp_rg`) is **Krook/plasma-Z**, normalized
  `/(2³π²)`. We use it as a **reference only**, not restored:
  - it **validates our `b₊/b×`**: `eval_bp = ρ_L²/2·(kperp²+kperpp²)`,
    `eval_bt = ρ_L²·kperp·kperpp` — identical shape to `hatG_rho_phi`;
  - it carries the **large-argument Bessel asymptotics** (`bessel_large_arg_limit` branch,
    `exp(bt−bp)/√(2π bt)` etc.) — lift this for Task 6 (adapt `gsl_sf_bessel_In` → `bessel_in`);
  - its `fill_rho_kernels` is the ready-made **Phase-2 assembly inner loop** (3-D
    `(kr,krp,rg)` grid).
- **Background → periodize** (REVISED 2026-07-10, superseding the earlier truncate decision).
  Adopt the author's existing periodization code
  (`~/1_projects/KIM_forced_periodicity/FORCED_PERIODICITY/SRC/`: `make_periodic.f90`,
  `localizer.f90`) — sample the background *periodically* over one period
  `L = 2(Δr + Δr_tr)` (`Δr` = as-is half-width `dx_asis`, `Δr_tr` = transition half-width
  `dx_tr`), with `localizer` (a C∞ transition, `exp(-2π/(1-t)·exp(-√2/t))`) blending the
  as-is layer into the periodic images so all derivatives stay continuous.

  **Key consequence (memo eq. `matelrhoPhi_discr`):** the matrix element uses the SAME `G` we
  built — `K^{ρΦ}_{m,m'} = (2π/L) ∫_{rm−L/2}^{rm+L/2} dr_g G(k_m, k_{m'}, r_g)` — and because
  the periodized background makes `G(k_m,k_{m'},r_g+L)=G(k_m,k_{m'},r_g)` periodic, the `r_g`
  integral over one period is spectrally exact with equidistant (lowest-order) quadrature.
  **Phase-1 kernel work is therefore fully reused, unchanged.** What changes is Phase 2
  (background sampling → periodic via `make_periodic`; equidistant one-period `r_g` grid) and
  Phase 5 (error to quantify = periodization *deformation* controlled by `Δr`, not truncation;
  strategy A scans `Δr`, and `Δr_tr` controls smoothness).

**`k_s²` note:** the canonical `kernel_mod.f90` sets `kperp = √(ks²+kr²)` — i.e. it **includes
`k_s²`** (= our `kern_include_ks2=.true.` path). Our FP diagonal source-of-truth drops `k_s²`
(an r=0 diagnostic hack). We keep the toggle **default `.false.`** (exact diagonal collapse) as
decided, but record that the physically-canonical form is `.true.`.

**Adopt from the issue regardless:** the orphaned `test_kernel_rho_phi.f90` analytic Debye
value `−1/(2³π²λ_D²)` is a concrete anchor — but that is the **Krook `/(2³π²)` normalization**;
our FP kernel is `/(4π)` with a different adiabatic term, so the FP Debye limit differs. Our
Task 4 (`b→0` closed form) is the FP analogue; a full FP Debye cross-check is a Phase-5 item.
The Vaclavik cross-benchmark (issue) → Phase 5.

## Ground truth captured during planning (so the executor has zero guesswork)

- **Diagonal ρΦ integrand** — `FLR2_asymptotics.f90:189-249`, module `flr2_asymptotics_m`.
  Per species `sp`, per rg-grid index `j`, with `b = kr²·plasma%spec(sp)%rho_L(j)²`:
  - Debye term (if `artificial_debye_case <= 1`): `-1/λ_D(j)²`.
  - FP term (if `artificial_debye_case == 0 .or. == 2`):
    `(1/λ_D²)·i·vT²·ks/(ω_c·ν)·exp(-b)·[ I00(j,0)·(I₀(b)·(A1+A2·(1-b)) + A2·b·I₋₁(b)) + ½·I20(j,0)·A2·I₀(b) ]`
    where `I₀=bessel_in(0,b)`, `I₋₁=bessel_in(-1,b)`, all `plasma%spec(sp)%…(j)`.
  - Species-summed, then `kernel_phi = (Σ_sp …)/(4π)`.
- **Diagonal ρB integrand** — same block, `FLR2_asymptotics.f90:216-224`: identical shape but
  prefactor `-1/λ_D²·vT³/(ω_c·ν·sol)`, and susceptibility functions `I01(j,0)`, `I21(j,0)`
  instead of `I00`, `I20`. Also `/(4π)`. (No `com_unit`, no `ks` factor, and negative sign.)
- **`bessel_in`** is the *unscaled* modified Bessel Iₙ (diagonal multiplies by `exp(-b)`
  explicitly). `bessel_in(-1,b) = bessel_in(1,b)`. Source: linked `fortnum_special`.
- **`ρ_L = |vT/ω_c|·ion_flr_scale_factor`** (`species_mod.f90:482-484`);
  `ion_flr_scale_factor` default `1.0` (`KIM/nmls/KIM_config.nml:14`).
- **Susceptibility funcs** `I00,I20,I01,I21,I11,I13` are `complex(dp)`, live in `species_t`
  (`species_mod.f90:49-57`), indexed `(rg-grid j, mphi)`; `mphi=0` is the default harmonic.
  Populated by `calculate_thermodynamic_forces_and_susc` → `getIfunc(x1,x2,symbI)`, then
  interpolated to `rg_grid%xb` via `interpolate_plasma_backs`.
- **Callers of the diagnostic** — `poisson.f90:220`, `flr2_benchmark.f90:161` (both write
  `hatK_*` files). Only the electromagnetic golden case runs in CI today, so a
  behavior-preserving refactor of `calc_hatK_Phi_in_Fourier` disturbs no golden record.
- **Test harness** — CTest programs in `KIM/tests/`, each a `program` that ends `error stop`
  / `stop 1` on failure, registered in `KIM/tests/CMakeLists.txt` with
  `target_link_libraries(<t> KIM_lib kilca_lib lapack cerf ddeabm slatec)`. In-memory plasma
  setup recipe: `test_kim_diagnostics.f90:146-167`
  (`make_test_profiles` → `write_test_namelist` → `nml_config_path=…` →
  `profiles_in_memory=.true.` → `kim_init()` → `set_profiles_from_arrays(...)` →
  `from_kim_factory_get_kim('electrostatic', kim)` → `kim%init()`), which leaves the global
  `species_m::plasma` populated on the rg grid. The config must set
  `collision_model = 'FokkerPlanck'` so the susceptibility functions are computed.
- **New source registration** — add the new file to `set(KIM_asymptotics …)` in
  `KIM/src/CMakeSources.in` (currently one line: `asymptotics/FLR2_asymptotics.f90`).

---

## Task 1: De-risk the test harness — new module + a test that gets a populated plasma

Prove we can (a) create a new asymptotics module and link a new CTest program, and (b) obtain
a fully-populated `plasma` (with rg-grid susceptibility functions) inside a unit test —
BEFORE writing any kernel math.

**Files:**
- Create: `KIM/src/asymptotics/flr2_fourier_kernel.f90`
- Modify: `KIM/src/CMakeSources.in` (append to `KIM_asymptotics`)
- Create: `KIM/tests/test_flr2_fourier_kernel.f90`
- Modify: `KIM/tests/CMakeLists.txt` (register the new test)

**Step 1 — Create the module stub.**

```fortran
module flr2_fourier_kernel_m
    ! Fused Fokker-Planck off-diagonal Fourier-space kernel integrand
    ! G(k_r, k'_r, r_g) for the forced-periodicity solver. Collapses to the
    ! diagonal source-of-truth calc_hatK_Phi_in_Fourier at k_r = k'_r.
    use KIM_kinds_m, only: dp
    use species_m, only: plasma_t
    implicit none
    private
    public :: hatG_rho_phi, hatG_rho_B
contains
    complex(dp) function hatG_rho_phi(plasma_in, kr, krp, j) result(G)
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr, krp
        integer, intent(in) :: j
        G = (0.0_dp, 0.0_dp)   ! stub — replaced in Task 3
    end function hatG_rho_phi

    complex(dp) function hatG_rho_B(plasma_in, kr, krp, j) result(G)
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr, krp
        integer, intent(in) :: j
        G = (0.0_dp, 0.0_dp)   ! stub — replaced in Task 5
    end function hatG_rho_B
end module flr2_fourier_kernel_m
```

**Step 2 — Register the source.** In `KIM/src/CMakeSources.in`, change
`set(KIM_asymptotics asymptotics/FLR2_asymptotics.f90)` to also list
`asymptotics/flr2_fourier_kernel.f90`.

**Step 3 — Write the harness test** `KIM/tests/test_flr2_fourier_kernel.f90`. Copy the
in-memory setup helpers `make_test_profiles` and `write_test_namelist` from
`test_kim_diagnostics.f90` (config must have `collision_model='FokkerPlanck'`,
`type_of_run='electrostatic'`). The `main` program:
1. runs the setup sequence (profiles → namelist → `kim_init` → `set_profiles_from_arrays` →
   factory `'electrostatic'` → `kim%init()`);
2. `use species_m, only: plasma`; picks an interior rg-grid index `j` (e.g. `rg_grid%npts_b/2`);
3. asserts `plasma%spec(0)%I00(j,0)` is finite and `plasma%spec(0)%rho_L(j) > 0` and
   `plasma%spec(0)%lambda_D(j) > 0` (i.e. the background really is populated);
4. calls `hatG_rho_phi(plasma, 1.0_dp, 1.0_dp, j)` and prints it;
5. prints `All tests PASSED`; `stop 0`.

> If `kim%init()` alone does not populate `plasma%spec%I00` on the rg grid, trace the
> electrostatic setup (`poisson.f90` init/run before line 220) to find the routine that runs
> `calculate_thermodynamic_forces_and_susc` + `interpolate_plasma_backs`, and call up to that
> point. Resolving this is the whole point of Task 1 — do not proceed to Task 2 until a unit
> test can read a populated `plasma%spec(0)%I00(j,0)`.

**Step 4 — Register the test** in `KIM/tests/CMakeLists.txt` (mirror an existing block):

```cmake
add_executable(test_flr2_fourier_kernel ${CMAKE_SOURCE_DIR}/KIM/tests/test_flr2_fourier_kernel.f90)
set_target_properties(test_flr2_fourier_kernel PROPERTIES
                        OUTPUT_NAME test_flr2_fourier_kernel.x
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")
target_link_libraries(test_flr2_fourier_kernel KIM_lib kilca_lib lapack cerf ddeabm slatec)
add_test(NAME test_flr2_fourier_kernel
         COMMAND ${CMAKE_BINARY_DIR}/tests/test_flr2_fourier_kernel.x)
```

**Step 5 — Build & run.**
Run: `cmake --build build --target test_flr2_fourier_kernel.x && ctest --test-dir build -R test_flr2_fourier_kernel --output-on-failure`
Expected: builds, PASS, prints a finite `I00` and a `(0,0)` kernel stub value.

**Step 6 — Commit.**
`git add KIM/src/asymptotics/flr2_fourier_kernel.f90 KIM/src/CMakeSources.in KIM/tests/test_flr2_fourier_kernel.f90 KIM/tests/CMakeLists.txt`
`git commit -m "test(kim): scaffold fused FLR2 Fourier kernel + harness"`

---

## Task 2: Extract the diagonal integrand into a shared function (behavior-preserving)

Make `calc_hatK_Phi_in_Fourier` delegate its per-`(sp,j,kr)` ρΦ and ρB integrand to new
functions in `flr2_fourier_kernel_m`, so the diagonal math has exactly one home. This is a
pure refactor — outputs must be byte-identical.

**Files:**
- Modify: `KIM/src/asymptotics/flr2_fourier_kernel.f90` (add `hatG_rho_phi_diag_sp`,
  `hatG_rho_B_diag_sp` — per-species, at a single wavenumber `kr`, no phase, no `/(4π)`)
- Modify: `KIM/src/asymptotics/FLR2_asymptotics.f90:189-249` (call the new functions)
- Modify: `KIM/tests/test_flr2_fourier_kernel.f90` (add the characterization assertion)

**Step 1 — Write the failing characterization test.** In the test program, after setup,
add a subroutine `test_diagonal_matches_inline()` that, for `kr ∈ {0.1,1.0,5.0}` and a few
`j`, computes the diagonal value **two ways**: (i) an inline copy of the exact expression
from `FLR2_asymptotics.f90:200-248` (Debye + FP term, summed over species), and (ii)
`sum over sp of hatG_rho_phi_diag_sp(plasma, sp, kr, j)`. Assert
`abs(inline - fromfunc) < 1e-12 * (1 + abs(inline))`. (The inline copy is the anchor; it is
deleted in Step 5 once the function is the single source.)

**Step 2 — Run red.**
Run: `ctest --test-dir build -R test_flr2_fourier_kernel --output-on-failure`
Expected: FAIL (function returns 0 stub ≠ inline value).

**Step 3 — Implement `hatG_rho_phi_diag_sp` / `hatG_rho_B_diag_sp`** in the module, copying
the exact per-species Debye + FP expression (ρΦ: `I00`,`I20`; ρB: `I01`,`I21`, its own
prefactor/sign), honoring `artificial_debye_case`, `turn_off_ions`, `turn_off_electrons`
from `config_m`. `b = kr²·plasma%spec(sp)%rho_L(j)²`; `I₀=bessel_in(0,b)`,
`I₋₁=bessel_in(-1,b)`. Do NOT apply `/(4π)` or phase here (caller does).

**Step 4 — Run green.** Expected: PASS.

**Step 5 — Refactor the caller.** Replace the inline bodies at
`FLR2_asymptotics.f90:200-248` with
`kernel_phi_temp = kernel_phi_temp + hatG_rho_phi_diag_sp(plasma_in, sp, kr, j)` and the ρB
analogue (keep the `/(4π)` scaling and file writes untouched). Verify no golden/diagnostic
change: build `KIM.x`, or if a quick check is unavailable, rely on the byte-identical math
(the expression moved verbatim).

**Step 6 — Run full KIM test suite.**
Run: `ctest --test-dir build --output-on-failure`
Expected: all previously-passing tests still PASS.

**Step 7 — Commit.**
`git commit -am "refactor(kim): extract diagonal FLR2 Fourier integrand to shared fn"`

---

## Task 3: Generalize to two wavenumbers (b₊, b×, phase) — default (drop k_s²) path

**Prerequisite:** ks² convention resolved (keep-k_s²-behind-toggle; see top).

**Files:**
- Modify: `KIM/src/asymptotics/flr2_fourier_kernel.f90` (add `kern_include_ks2` module var,
  private `core_rho_phi_sp(plasma,sp,bplus,bcross,j)`, implement `hatG_rho_phi`; refactor
  `hatG_rho_phi_diag_sp` to call `core_rho_phi_sp` with `bplus=bcross=b`)
- Modify: `KIM/tests/test_flr2_fourier_kernel.f90` (collapse + phase tests)

This task implements ONLY the default path (`kern_include_ks2 = .false.`). Leave the variable
declared and defaulted `.false.`; the `.true.` branch lands in Task 3b.

**Step 1 — Write the failing collapse test.** `test_collapse_rho_phi()`
(with `kern_include_ks2 = .false.`): for `kr ∈ {0.1,1.0,5.0}` and several `j`, assert
`abs(hatG_rho_phi(plasma,kr,kr,j) - dref) < 1e-12*(1+abs(dref))`, where
`dref = (Σ_sp hatG_rho_phi_diag_sp(plasma,sp,kr,j))/(4π)` (the diagonal value **with** the
`/(4π)`, phase = 1). Also `test_phase_rho_phi()`: for `kr≠krp`, assert the two-wavenumber
core is phase-only-different, i.e.
`abs( hatG_rho_phi(plasma,kr,krp,j)*exp(-ci*(kr-krp)*rg) - hatG_rho_phi(plasma,krp,kr,j)*exp(-ci*(krp-kr)*rg) ) < 1e-12`
(symmetric `b₊,b×` ⇒ equal cores), with `rg = rg_grid%xb(j)`, `ci=(0,1)`.

**Step 2 — Run red.** Expected: FAIL (stub returns 0).

**Step 3 — Implement the default path.** Introduce the private helper
`core_rho_phi_sp(plasma,sp,bplus,bcross,j)` holding the per-species FP+Debye expression with
`exp(-bplus)`, `bessel_in(0,bcross)`, `bessel_in(-1,bcross)`, and `(1-bplus)`. Refactor
`hatG_rho_phi_diag_sp` (Task 2) to call it with `bplus=bcross = kr²·rho_L²`. Implement
`hatG_rho_phi`: for each `sp`, compute (default branch)
`b₊ = rho_L²·(kr²+krp²)/2`, `b× = rho_L²·abs(kr*krp)`, sum `core_rho_phi_sp(...,b₊,b×,j)`,
then multiply the species sum by `exp(ci*(kr-krp)*rg_grid%xb(j))/(4π)`. Collapse is exact by
construction (at `kr=krp`, `b₊=b×=kr²rho_L²`, phase=1).

**Step 4 — Run green.** Expected: collapse + phase tests PASS.

**Step 5 — Commit.**
`git commit -am "feat(kim): fused two-wavenumber FP kernel G^rhoPhi (drop-k_s2 path)"`

---

## Task 3b: k_s²-inclusive path behind `kern_include_ks2` + structural collapse test

**Files:**
- Modify: `KIM/src/asymptotics/flr2_fourier_kernel.f90` (`hatG_rho_phi` `.true.` branch)
- Modify: `KIM/tests/test_flr2_fourier_kernel.f90` (`test_collapse_rho_phi_ks2()`)

**Step 1 — Write the failing test** `test_collapse_rho_phi_ks2()`. Set
`kern_include_ks2 = .true.`. For `kr ∈ {0.1,1.0,5.0}` and several `j`, build the reference by
calling the shared core at the **full** diagonal argument:
`ref = (Σ_sp core_rho_phi_sp(plasma,sp,bf,bf,j))/(4π)` with `bf = (plasma%ks(j)²+kr²)*rho_L(j)²`
per species. Assert `abs(hatG_rho_phi(plasma,kr,kr,j) - ref) < 1e-12*(1+abs(ref))`.
(`core_rho_phi_sp` must be accessible to the test — either make it `public`, or expose a thin
public test shim `core_rho_phi_sp_pub`.) Reset `kern_include_ks2 = .false.` at test end so
later tests keep the default. **Guard:** the default-path tests from Task 3 must still pass
with the toggle back off.

**Step 2 — Run red** (the `.true.` branch not yet implemented → returns default-path value ≠ ref).

**Step 3 — Implement the `.true.` branch** in `hatG_rho_phi`: with `ks = plasma%ks(j)`,
`b₊ = rho_L²·(2*ks²+kr²+krp²)/2`, `b× = rho_L²·sqrt(ks²+kr²)*sqrt(ks²+krp²)`; same
`core`/phase/`(4π)` assembly.

**Step 4 — Run green.** Expected: `test_collapse_rho_phi_ks2` PASS and Task 3 tests still PASS.

**Step 5 — Commit.**
`git commit -am "feat(kim): k_s2-inclusive FP kernel path behind kern_include_ks2 toggle"`

---

## Task 4: Analytic b→0 / homogeneous limit test

Validate the Bessel/exp assembly against a closed form, independent of the diagonal code.

**Files:** Modify `KIM/tests/test_flr2_fourier_kernel.f90`.

**Step 1 — Write the failing test** `test_zero_wavenumber_limit()`. At `kr=krp=0`:
`b₊=b×=0`, `exp(0)=1`, `I₀(0)=1`, `I₋₁(0)=0`, phase `=1`. So the closed form per species is
`debye + (1/λ_D²)·ci·vT²·ks/(ω_c·ν)·[ I00(j,0)·(A1+A2) + ½·I20(j,0)·A2 ]`, summed, `/(4π)`.
Compute this inline from `plasma%spec` and assert
`abs(hatG_rho_phi(plasma,0,0,j) - closed) < 1e-12*(1+abs(closed))` at several `j`.

**Step 2 — Run red / Step 3 — (already implemented; should pass) / Step 4 — green.**
If it fails, the bug is in the Bessel-argument generalization; fix in Task 3's code.

**Step 5 — Commit.** `git commit -am "test(kim): b->0 analytic limit for fused FP kernel"`

---

## Task 5: ρB kernel `hatG_rho_B` — same generalization + tests

**Files:** Modify `flr2_fourier_kernel.f90` (`hatG_rho_B` via shared `core_rho_B_sp`),
`test_flr2_fourier_kernel.f90` (collapse + b→0 tests for ρB, mirroring Tasks 3–4).

Steps mirror Task 3 **and Task 3b** (default drop-k_s² collapse test → implement default
branch; then `kern_include_ks2=.true.` structural collapse test → implement `.true.` branch;
reusing the shared `core_rho_B_sp(plasma,sp,bplus,bcross,j)` helper) and Task 4 (b→0:
`I₋₁(0)=0`, `I₀(0)=1`). `hatG_rho_B` honors the same `kern_include_ks2` toggle. Commit:
`git commit -am "feat(kim): fused two-wavenumber FP kernel G^rhoB + tests"`

---

## Task 6: Large-`b×` overflow guard (numerical robustness)

Unscaled `bessel_in(0,b×)` overflows for large `b×` (large `k_r`), even though `exp(-b₊)`
tames the product (`b₊ ≥ b×` by AM–GM). Lift the scaled/asymptotic branch from the canonical
`1ae0aeeb~1:KIM/src/kernels/kernel_mod.f90` (`bessel_large_arg_limit` branch:
`eval_besselI0 = exp(bt-bp)/√(2π bt)`, `eval_besselIm1 = exp(-bp + asinh(-1/bt) + bt√(1+1/bt²))/√(2π bt√(1+1/bt²))`)
so `exp(-b₊)·I₀(b×)` and `exp(-b₊)·I₋₁(b×)` are evaluated without intermediate overflow. (Same
form appears in `bfa9d70c:src/deprecated/integrands.f90:157-165`.)

**Files:** Modify `flr2_fourier_kernel.f90` (`core_*_sp` helpers), `test_flr2_fourier_kernel.f90`.

**Step 1 — Write the failing test** `test_large_kr_finite()`: assert
`hatG_rho_phi(plasma, 50.0_dp, 50.0_dp, j)` is finite (`x == x` and `abs < huge`). With naive
unscaled Bessel this is NaN/Inf → red.

**Step 2 — Run red.** **Step 3 — Implement** the `if (b× > ~10) then <asymptotic> else
<bessel_in>` branch computing the combined `exp(-b₊)·Iₙ(b×)` directly. **Step 4 — green.**
Re-run the collapse test (Task 3) to confirm the branch still matches the diagonal in its
overlap. **Step 5 — Commit** `git commit -am "fix(kim): overflow-safe Bessel in fused FP kernel"`

---

## Phase 1 exit criteria

- `ctest --test-dir build -R test_flr2_fourier_kernel --output-on-failure` — all green.
- `ctest --test-dir build --output-on-failure` — no regressions in the existing suite.
- `hatG_rho_phi`, `hatG_rho_B` public in `flr2_fourier_kernel_m`; diagonal math has one home
  (shared with `calc_hatK_Phi_in_Fourier`); default-path collapse exact to 1e-12;
  `kern_include_ks2=.true.` path structurally verified; b→0 limit verified; overflow-safe.
- `kern_include_ks2` module toggle present, default `.false.`, both paths tested.
- **Checkpoint — report and wait before Phase 2** (assembly + dense `zgesv` solve + new
  `electrostatic_periodic_t` run-type).

## Deferred to later phases (recorded so Phase 1 stays scoped)

- A/B convergence comparison of the two `kern_include_ks2` paths — Phase 5 (strategy A).
- Wiring `kern_include_ks2` to the config namelist (module-var only for now) — Phase 2.
- Interpolating variant of the kernel onto a window grid `[rm−W, rm+W]` — Phase 2.
- `mφ ≠ 0` cyclotron harmonics — off by default (design §3.2); switch-on for verification only.
- Current-diagnostic kernels `K^{jΦ}, K^{jB}` — out of scope (design §9).
