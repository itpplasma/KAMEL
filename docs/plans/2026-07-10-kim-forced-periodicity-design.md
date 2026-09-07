# KIM Forced-Periodicity Solver — Design

- **Date:** 2026-07-10
- **Branch:** `feature/kim-forced-periodicity`
- **Status:** Design (brainstormed + grilled); ready to turn into an implementation plan
- **Primary spec:** `~/1_projects/KIM_forced_periodicity/forced_periodicity.pdf`
- **Supporting theory:** `~/1_projects/PhD_Thesis.pdf` (M. Markl, 2024), Part III — esp. Ch. 13 (Krook kernels), Ch. 14 (Fokker–Planck kernels), Ch. 15 (variational implementation)

---

## 1. Motivation

KIM currently solves the plasma-response problem *globally*: hat-function basis over the
whole radius, the integral kernel built in Fourier space then transformed to the hat basis
(thesis §15.2), giving a sparse-but-approximated matrix. In practice the observed solutions
are **highly localized** at the resonant surface (localization scale ~ ion Larmor radius
`ρ_L`), so the sparsity advantage is lost and the basis transform carries avoidable
approximations.

The **forced-periodicity** approach exploits the *locality* of the response: a few
centimetres from the resonant surface, in the ideal region, the potential is the known
aligned-MHD solution plus a small FLR2 correction, and it is decoupled from the resonant
layer. This lets us:

- represent the perturbation in a **periodic Fourier-harmonic basis** `e^{i k_m r}`,
  `k_m = 2πm/L`, over a **finite radial window** centred on the resonant surface `rm`;
- build the matrix elements directly from the analytic Fourier-space kernel
  `G(k_r, k'_r, r_g)` (spec Eq. 15–16), **avoiding the hat-basis transform** entirely;
- solve a small **dense** system instead of a large sparse one.

### 1.1 Periodize the background (per the spec; REVISED 2026-07-10)

**Earlier this design deviated from the spec by keeping the true, non-periodic background and
truncating the integration (single knob `W = L/2`). That decision is REVERSED: we follow the
spec and periodize the background**, reusing the author's existing periodization code
(`~/1_projects/KIM_forced_periodicity/FORCED_PERIODICITY/SRC/make_periodic.f90`,
`localizer.f90`).

The background is made periodic over period `L = 2(Δr + Δr_tr)`: the **as-is** region
`[rm − Δr, rm + Δr]` (`Δr = dx_asis`) keeps the true background, and a **transition** region
of half-width `Δr_tr = dx_tr` on each side uses `localizer` — a C∞ transition
`exp(−2π/(1−t)·exp(−√2/t))` — to blend smoothly into the periodic images so all derivatives
stay continuous. Consequences:

- `G(k_r, k'_r, r_g)` **is** periodic in `r_g` (spec §; memo eq. `matelrhoPhi_discr`), so the
  matrix element `K_{m,m'} = (2π/L) ∫_{rm−L/2}^{rm+L/2} dr_g G(k_m,k_{m'},r_g)` is evaluated
  with **cheap lowest-order equidistant quadrature** that is spectrally exact over one period.
- Two geometric knobs: `Δr` (as-is half-width — controls how much of the resonant layer is
  kept true) and `Δr_tr` (transition half-width — controls smoothness of the periodic
  extension); together `L = 2(Δr + Δr_tr)`.
- The error to quantify is now the **periodization deformation** — how much periodizing the
  background outside the as-is layer perturbs the resonant-layer solution — controlled by
  `Δr` (larger `Δr` ⇒ less deformation). This replaces the truncation-error interpretation.

The FP kernel `G` is unchanged by this choice (the memo assembles the matrix element from the
same `G`); periodization affects only how the background is **sampled** and how the `r_g`
integral is discretized.

---

## 2. Resolved design decisions

| # | Topic | Decision |
|---|-------|----------|
| 1 | Deliverable | Working finite-periodic Fourier solver **and** a quantified periodization error. |
| 2 | Integration | New `electrostatic_periodic_t` run-type behind the `kim_solver_t` seam; reuses config/grid/species/equilibrium and the FP susceptibility functions. |
| 3 | System | **Electrostatic**, mirroring the existing `electrostatic` run-type: `(Δ + 4π K^{ρΦ}) Φ = −4π K^{ρB} Br`. |
| 4 | Kernel | **Full FP** kernel `G(k_r,k'_r,r_g)` (thesis Eq. 14.3–14.4), evaluated at **`m_φ = 0` by default**. |
| 5 | Kernel source | Fuse the FP susceptibility assembly (`calc_hatK_Phi_in_Fourier`) with the off-diagonal two-wavenumber structure (deprecated k-space integrand). |
| 6 | Basis / window | Periodic Fourier harmonics on the finite window; `ρ_L`-anchored resolution; three-knob convergence. |
| 7 | R1 | Solve for the deviation `ψ = δΦ − Φ_MA`; first-cut naive `Φ_MA` projection, then measure. |
| 8 | Error strategy | **A** window-convergence + **B** cross-check vs the `electrostatic`-global solver; metric = resonant-layer `δΦ(r)`. |
| 9 | Validation floor | Diagonal-collapse + homogeneous kernel unit tests; a frozen golden case. Approx-vs-approx for cross-check (no analytic anchor initially). |

---

## 3. Physics core

### 3.1 The Fourier-space kernel `G(k_r, k'_r, r_g)`

The charge-density response is `δρ(r) = ∫ dr' K^{ρΦ}(r,r') δΦ(r')` (spec Eq. 10). The kernel
has the continuous Fourier-space representation (spec Eq. 14–15)

```
K^{ρΦ}_{k_r,k'_r} = ∫ dr_g G(k_r, k'_r, r_g),
```

so the finite-window matrix element in the periodic basis is (spec Eq. 16, but over the
**true** non-periodic background):

```
K^{ρΦ}_{m,m'} = (2π/L) ∫_{rm−W}^{rm+W} dr_g  G(k_m, k_{m'}, r_g).
```

**Where `G` lives in the code (verified):**

- **Diagonal `G(k_r, k_r, r_g)`** is already implemented in
  `KIM/src/asymptotics/FLR2_asymptotics.f90 :: calc_hatK_Phi_in_Fourier`. It assembles all
  four kernels (ρΦ, ρB, jΦ, jB) from the FP susceptibility functions
  (`A1, A2, I00, I20, I01, I21, I11, I13, ν, vₜ, ω_c, λ_D`) at a *single* wavenumber
  (`b = k_r² ρ_L²`). It is currently a **diagnostic** (loops over hardcoded
  `k_r ∈ {0.1, 1, 5}`, writes to files), not a matrix builder.
- **Off-diagonal two-wavenumber structure** (`b_+`, `b_×`, and the phase
  `exp(i(k_r − k'_r) r_g)`) exists in git history in `src/deprecated/integrands.f90`
  (`integrand_K_rho_phi_krp`, developed across commits
  `9d78ee34 → 292eeb68 → 4f6b1adf → bfa9d70c`), but in the **plasma-Z (Krook)** flavour.

**Neither source alone is the FP off-diagonal `G`.** The revival is a *fusion*: lift the FP
susceptibility assembly (from `calc_hatK_Phi_in_Fourier`) to two wavenumbers using the
off-diagonal FLR structure

```
b_+ = vₜ²/(2ω_c²) (2k_s² + k_r² + k'_r²)
b_× = vₜ²/ω_c²  √(k_s² + k_r²) · √(k_s² + k'_r²)
```

with `exp(−b_+)`, `I₀(b_×)`, `I₋₁(b_×)` (full Bessel functions, plus the large-argument
asymptotics already in the deprecated code), and the Fourier phase `exp(i(k_r − k'_r) r_g)`.
In the deprecated code the two-wavenumber dependence enters *only* through these FLR factors
and the phase, while the susceptibility coefficients stay local at `r_g` — so the lift is
structurally well-defined.

### 3.2 Cyclotron-harmonic truncation (`m_φ = 0`)

Thesis Eq. 14.3–14.6 give the FP kernels with an **explicit, unresolved sum over cyclotron
harmonics** `Σ_{m_φ}`, each carrying `I_{m_φ}(b_×)`, `I_{m_φ−1}(b_×)` and a harmonic-shifted
collisionality `x₂ = −(ω_E + m_φ ω_c − ω)/ν_σ` (Eq. 14.2). The thesis leaves resolving this
sum as future work.

**Decision:** target the full FP kernel structure but evaluate **`m_φ = 0` only by default**.
For `m_φ ≠ 0`, `x₂ ≈ m_φ ω_c / ν` is very large, and the FP susceptibility functions
`I^{kl}(x₁, x₂)` are strongly suppressed (the cyclotron resonances are far off-resonance and
collisionally damped). Keep `m_φ` as a structurally-present loop (default range `{0}`) so a
harmonic can be switched on to *verify* its smallness, but it never runs by default.

### 3.3 Gap-2 sign (resolved)

At the diagonal (`m_φ = 0`, `b_+ = b_× = b`), `calc_hatK_Phi_in_Fourier` uses the coefficient
`A₂·(1 − b)`, whereas thesis Eq. 14.3 reads `A₂·(1 + b_+)`. **The code sign is authoritative**
(user-confirmed); `calc_hatK_Phi_in_Fourier` is the kernel source of truth, and the thesis
`(1 + b_+)` is a documentation discrepancy — do **not** propagate it.

### 3.4 The linear system (electrostatic)

Mirroring the `electrostatic` run-type (`electrostatic_poisson/poisson.f90` +
`solve_poisson.f90`), in the Fourier basis:

```
[ D_m δ_{m,m'} + 4π K^{ρΦ}_{m,m'} ] Φ_{m'} = −4π K^{ρB}_{m,m'} Br_{m'},
```

with `D_m` the Fourier-space Poisson operator (`−k_m²` up to the radial-operator convention)
and `Br_{m'}` the Fourier coefficients of the external radial-field drive over the window
(trivial for the default `type_br_field = 12`, a constant `Br` → `m = 0` only). Dense complex
solve (LAPACK `zgesv`) for `Φ_m`, then reconstruct `δΦ(r) = Σ_m Φ_m e^{i k_m r}`.

Only `K^{ρΦ}` and `K^{ρB}` are needed for the electrostatic solve. The current-diagnostic
kernels `K^{jΦ}, K^{jB}` (Eq. 14.5–14.6) are deferred.

### 3.5 R1 — solve for the deviation from the aligned potential

A periodic Fourier basis implicitly forces `δΦ(rm−W) = δΦ(rm+W)`, but the true `δΦ` tends to
the (non-periodic) aligned potential at the edges → Gibbs oscillations. Mitigation: solve for
the **deviation** `ψ = δΦ − Φ_MA`, which is localized to the resonant layer and decays to ~0
at the edges (periodic-basis-friendly). `Φ_MA` is already computed by
`calc_flr2_asymptotic_Phi_MA` (output `EBdat%Phi_MA_asymptotic`).

Substituting `Φ = Φ_MA + ψ` into `L Φ = g` (with `L = Δ + 4π K^{ρΦ}`,
`g = −4π K^{ρB} Br`):

```
L ψ = g − L Φ_MA  ≡  s.
```

The source `s` is physically **localized** (Φ_MA solves the equation in the ideal region, so
`s ≈ 0` there), hence its Fourier projection is clean. The only difficulty is computing
`L Φ_MA = ΔΦ_MA + 4π K^{ρΦ} Φ_MA`: the `K^{ρΦ} Φ_MA` term needs `Φ_MA`'s Fourier
coefficients, and `Φ_MA` is non-periodic at the edges.

**Approach (tracer-bullet):** first-cut **naive projection** of `Φ_MA` with ample modes;
get the pipeline running; use strategy A to measure whether the residual Gibbs actually
degrades `ψ`-convergence for the (6,2) case. Escalate to (a) an analytic FLR2 localized-source
derivation, or (b) polynomial/ramp subtraction, **only if measured necessary**.

**Note:** R1 is new machinery — the existing `electrostatic` solver solves for the full `Φ`
directly and uses `Φ_MA_asymptotic` only as a diagnostic.

---

## 4. Numerics

### 4.1 Window and spectral resolution — `ρ_L`-anchored, three-knob

Anchor length scale: the ion Larmor radius `ρ_L` at the resonant surface `rm`.

- `k_max ≈ c_k / ρ_L` (default `c_k ≈ 5`) — resolve `ρ_L`-scale structure.
- `W = L/2 ≈ c_W · ρ_L` (default `c_W ≈ 15`) — contain the decayed deviation `ψ`.
- `M = k_max · W / π ≈ 24` → ~50 harmonics → a trivial dense `50×50` complex solve
  (validates the write-up's "minutes" claim; the hypothetical "1000" basis functions is a
  comfortable ceiling).
- `N_rg` (r_g quadrature points): the window is exactly one period `L`, so the phase
  `exp(i 2π(m−m') r_g / L)` completes whole cycles and equidistant sampling is *spectrally*
  accurate for the smooth part → `N_rg ≳ 4M`.

**Three convergence knobs** (used by strategy A):

1. `M` / `k_max` — spectral/basis error → converge to negligible.
2. `N_rg` — quadrature error → converge to negligible.
3. `W` — physical domain-truncation error → **the periodization error knob**; its plateau
   residual *is* the reported error bar.

**Subtlety:** because the true susceptibility is not periodic over the window, the `r_g`
quadrature residual is entangled with the physical truncation error — both stem from the same
non-periodicity. Report `M` and `N_rg` convergence separately from the `W` scan so the
physical error is not conflated with numerical error.

### 4.2 Assembly

For each `(m, m')`, evaluate `G(k_m, k_{m'}, r_g)` on the `N_rg` grid over the window (FP
susceptibility coefficients local at `r_g` × two-wavenumber FLR factors × phase) and quadrature
to `K^{ρΦ}_{m,m'}` / `K^{ρB}_{m,m'}`. Cost `~ M² · N_rg · n_species` susceptibility
evaluations — negligible at `M ~ 50`. The matrix is dense and **not** Toeplitz (the FLR
factors depend on both `k_m` and `k_{m'}`, not just their difference).

### 4.3 Boundary conditions

None imposed explicitly: the periodic Fourier basis is inherently periodic, and the deviation
formulation (`ψ → 0` at the edges) provides well-posedness. No Dirichlet plumbing.

---

## 5. Error quantification

### 5.1 Strategy A — window convergence (self-consistent error bar)

Scan the three knobs. With `M` and `N_rg` converged, grow `W`, extract `δΦ(r)` on the
resonant layer `|r − rm| ≤ Δr` (`Δr` = a few `ρ_L`, diagnostic-only), and report the Cauchy
residual `‖δΦ_i − δΦ_{i−1}‖ / ‖δΦ_i‖` in both L2 and max-norm. The plateau value at the
largest affordable `W` is the self-consistent periodization error bar.

### 5.2 Strategy B — cross-check vs the `electrostatic`-global solver

Run the existing `electrostatic` (global, hat-basis) run-type and the periodic run-type on the
**same** (6,2) case (`m_mode = −6`, `n_mode = 2`, resonant at `q = 3`, constant `Br`);
interpolate both `δΦ(r)` to a common resonant-layer grid; report
`‖δΦ_periodic − δΦ_global‖ / ‖δΦ_global‖`. Same physics on both sides → isolates the
finite-window + Fourier approximation.

**Caveat (accepted):** the global solver carries its own approximations (the basis transform
we are avoiding), so B is approx-vs-approx. An analytic anchor (homogeneous manufactured case)
is deferred; a small `δΦ_periodic − δΦ_global` on the localized (6,2) case is nonetheless
strong evidence.

---

## 6. Architecture

New run-type `electrostatic_periodic_t` (module `rt_electrostatic_periodic_m`) implementing
`class(kim_t)`, registered in the factory (`KIM/src/general/KIM_mod.f90 ::
from_kim_factory_get_kim`) under `type_of_run = "electrostatic_periodic"`. It reuses,
unchanged: config parsing, `grid_m`, `species_m`, `equilibrium_m`, and the FP susceptibility
functions on the profile grid. Output `δΦ(r)` is packed into the same `kim_results_t`
(`Es, Er, …, Phi`) so cross-check B and golden-record CI come almost for free.

**Solve pipeline** for a given `rm` and `(m, n)`:

1. Sample the **true** background on an equidistant `r_g` grid over `[rm−W, rm+W]`.
2. Build the Fourier basis `k_m = 2πm/L`, `m = −M…M`.
3. Assemble the dense `K^{ρΦ}_{m,m'}`, `K^{ρB}_{m,m'}` (§4.2).
4. Form the R1 deviation source `s` (§3.5, naive-projection first-cut).
5. Dense `zgesv` solve for `ψ_m`; reconstruct `ψ(r)`; add back `Φ_MA` → `δΦ(r)`.
6. Pack into `kim_results_t`.

---

## 7. Testing & CI

1. **Kernel unit tests.** The fused off-diagonal `G(k_r,k'_r,r_g)` must collapse to
   `calc_hatK_Phi_in_Fourier` on the diagonal `k_r = k'_r` (exact regression anchor) and to
   analytic limits (homogeneous background, `b → 0`).
2. **Assembly test.** Matrix elements vs a brute-force quadrature reference.
3. **Solver integration test.** Homogeneous case with a known `δΦ`.
4. **Physics validation.** Strategies A and B as CI-tracked cases.
5. **Golden record.** Freeze one periodic (6,2) case in `test/golden/kim` so the numbers are
   pinned against regressions.

Robustness: reuse the seam's status codes (`KIM_SOLVE_FAILED`, `KIM_GRID_COVERAGE`); report the
dense-solve condition number; fail soft, never halt the host; handle resonance-near-edge and
singular-matrix degeneracies explicitly.

---

## 8. Risks

- **R1 — periodic-basis edge mismatch.** Mitigated by the deviation formulation; residual
  handling escalates only if the first-cut naive projection is measured insufficient (§3.5).
- **R2 — oscillatory `r_g` quadrature.** The phase oscillates fast for large `|k_r − k'_r|`;
  managed by `N_rg ≳ 4M` (spectral accuracy over one period) plus the entanglement caveat
  (§4.1).
- **R3 — the magnetic drive.** `Br_{m'}` specification/normalization; trivial for the default
  constant `Br`, needs care for richer field models.
- **R4 — cost / conditioning** of the dense `(2M+1)²` solve as `M ∝ W`; expected trivial at
  `M ~ 50`, monitored via condition number.
- **R5 — kernel fusion correctness.** The FP off-diagonal `G` is derived by fusion, not copied;
  the diagonal-collapse unit test is the primary guard.

---

## 9. Baked-in assumptions

- **Kernel scope:** only `K^{ρΦ}` + `K^{ρB}`; current-diagnostic kernels `K^{jΦ}, K^{jB}`
  deferred.
- **Boundary conditions:** none explicit (periodic basis + `ψ → 0`).
- **Drive:** `Br` projected onto the Fourier basis over the window; trivial for
  `type_br_field = 12`.

---

## 10. Implementation phasing (tracer-bullet)

1. **Kernel `G`** — fuse FP susceptibility assembly + off-diagonal structure; diagonal-collapse
   unit test. Pure physics, no solver.
2. **Assembly + dense solve** — `electrostatic_periodic_t`, direct-`Φ` first (no R1);
   homogeneous sanity test.
3. **Cross-check B** — wire the `δΦ(r)` comparison vs `electrostatic`-global on (6,2).
4. **R1 deviation** — add the `ψ` formulation; measure Gibbs impact.
5. **Strategy A + golden** — three-knob convergence harness; freeze the golden case; CI.

---

## 11. Out of scope / future work

- **Toroidal geometry** (spec §C): the sub-integrand `G` is no longer a 1D integral; requires
  numerical guiding-centre orbit integrals (POTATO infrastructure).
- **Electromagnetic extension:** the same finite-window Fourier assembly applied to the coupled
  `(Φ, A_par)` block (mirroring `electromagnetic_t`).
- **Higher cyclotron harmonics** (`m_φ ≠ 0`): switch-on path exists for verification; not run by
  default.
- **Analytic anchor:** a homogeneous manufactured true-reference for strategy B.
- **Full-FP harmonic sum:** resolving `Σ_{m_φ}` analytically (the thesis' open problem).

---

## 12. Key references (code + theory)

- **Fourier kernel (diagonal, FP):** `KIM/src/asymptotics/FLR2_asymptotics.f90 ::
  calc_hatK_Phi_in_Fourier`.
- **Aligned potential `Φ_MA`:** same file, `calc_flr2_asymptotic_Phi_MA`.
- **Off-diagonal k-space integrand (git history, Krook):** `src/deprecated/integrands.f90 ::
  integrand_K_rho_phi_krp` @ `bfa9d70c` (and `9d78ee34 → 292eeb68 → 4f6b1adf`).
- **Electrostatic solve to mirror:** `KIM/src/electrostatic_poisson/poisson.f90`,
  `solve_poisson.f90`.
- **Solver seam:** `KIM/src/general/kim_solver_m.f90`; factory `KIM/src/general/KIM_mod.f90`.
- **Theory:** spec Eq. 9–16; thesis Ch. 13 (Krook kernels, §13.7 dispersion), Ch. 14 (FP
  kernels Eq. 14.1–14.6), Ch. 15 (variational implementation).
