# FLRE golden-record divergence — investigation log

Reference case: `test/golden/kilca/cases/flre_m6n2` (m=6, n=2).
Oracle: C++ KiLCA built at the port branch's merge-base with `main`
(`aa92bfb5`, equivalently the pre-port tree `3c364ab8`), which is what the
golden record was generated with. Comparison tool: `test/golden/bin/gr_numcompare.py`
(rtol=1e-7, atol=1e-12).

## Status — RESOLVED (2026-07-02)

The EB.dat/poy divergence was not irreducible 1F1 cancellation (the earlier
conclusion below is superseded). Four discrete bugs caused it; all are fixed
on this branch:

1. **Electron mass, 2 ULP** (`fix(kilca): mirror the C++ oracle's electron-mass
   literal in ported settings`). `constants.h` has `m_e = 9.10938185917485e-28`;
   the Fortran constants module computes `me = mp/1836.1526675`, 2 ULP away.
   Every ported electron-parameter path (`omc_`, `vT_`) needs the C++ literal,
   like the speed-of-light split above.
2. **libmvec vector cos** (`fix(kilca): keep scalar libm cos in the ported
   adaptive-grid unit`). gfortran has no math-errno semantics and vectorizes
   the Chebyshev-node loop to `_ZGVbN2v_cos` at `-O3`; g++ keeps scalar `cos`.
   Vector and scalar cos differ in last bits, so the adaptive sampling grid
   nodes sat 2 ULP off. Diagnostic: `nm -D <exe> | grep ZGV` must show the
   same set for port and oracle.
3. **1F1m via fortnum on both sides** (`port(kilca): route the 1F1m hot path
   through fortnum hyperg_1f1m_a1`). The modified form's low bits are
   codegen-dependent (analysis below still holds); the fix is one shared
   compiled kernel, not reproduction of cancellation garbage. Pairs with the
   C++-side swap in #146.
4. **DOMINANT: basis-reconstruction off-by-one** (`fix(kilca): start basis
   reconstruction at the last renorm's stored entries`). `solver_m.f90` passed
   post-advance indices to `renorm_basis_vecs`; the oracle passes
   `rdata-1, ydata-Neq, taudata-2*Nfs`. With at least one QR renormalization
   active, the first segment crossing inverted the orthonormal Q slot instead
   of the R factor and shifted every later segment lookup. The result was a
   local field defect at renorm radii: Poynting self-consistency 1.01 at
   r=5.03 against the healthy 2.9e-2, and 34 extra rows kept by the output
   thinning (2682 vs 2648). Harmless when no renorm fires, which is why unit
   tests and 56/58 files passed.

Verified against an oracle built as current `main` + #145 (fortnum Bessel
shim) + #146 (fortnum 1F1m): **58/58 kilca files pass** (EB.dat max_rel 0.0,
poy residual equal at 2.867e-2), and the QL-Balance `f_6_2` golden passes
**33/33 quantities at rtol 1e-8** including the `LinearProfiles` Br and D^ql
series. Against the pre-swap baseline the only residual is the deliberate
#145/#146 kernel-value change (EB.dat max_rel 9.7e-6), which lands with those
PRs and one golden-baseline re-bless.

The bisection method that found these: identical bit-pattern dump probes
(hex via `transfer`) inserted into legacy files present verbatim in both
trees, walked stage by stage (cti/cte at the RHS, K-matrix samples, W2 inputs,
adaptive-grid nodes, thinning input rows).

## Fixed and committed

### Speed-of-light precision split (`fix(kilca): mirror oracle's split speed-of-light precision`)

The C++ oracle carried two speed-of-light constants:
`constants.h` `const double c = 29979245800.0` (full double, used by C++ code)
and the Fortran `constants` module `c = 29979245800.0` (a bare literal gfortran
parses at single precision → stored as 29979246592.0, used by legacy Fortran
physics). The golden record encodes both: C++-path quantities (dPhi0, j0th) use
the exact value; Fortran-path quantities (`conduct_parameters` → Vth_m, j0th_m)
use the truncated one.

Fix: keep the Fortran `constants` module `c` at the legacy single-precision
value (`real(29979245800.0_sp, dp)`, bit-identical to the bare literal) so every
unported Fortran consumer matches with zero change; give the three
ported-from-C++ modules (`background_data_m`, `flre_zone_m`, `transforms_m`) a
local double-precision `c` mirroring `constants.h`. This moved dPhi0, j0th,
Vth_m, j0th_m from FAIL to PASS (54/58 → 56/58). ctest 36/36.

## Ruled out (with trace evidence, at r = 3.0 unless noted)

All checks instrumented identically in the oracle worktree and the port, run on
the same fixture, compared bit-for-bit.

- **Adaptive-grid parameters identical.** `DBGSGP`: pre-thinning grid dimension
  `dim_old = 92754` and `r_sum = 3.7931986556236963e+06` are **bit-identical**;
  `deg=5`, `step=0.1` identical. Only the post-thinning `dimnew` (648 vs 682)
  and the field-magnitude-scaled `sold_sum` differ. → the grid inputs are
  correct; the fields feeding the thinning differ.

- **Background pitch/field state bit-identical.** `DBGCOEF`/`back_data`: `ht_`,
  `hz_`, `dht_`, `dhz_`, `ddht_`, `ddhz_`, `dddh`, `r_`, `omega_`, and all 12
  spec-independent `back_data` slots (B, dB, ddB, dPhi0 + derivs) match to the
  bit.

- **Spec-dependent background state bit-identical.** `f0_x`: `n_`, `dn_`,
  `ddn_`, `Vp_`, `dVp_`, `vT_`, `dvT_`, `ddvT_`, `nu_`, `dnu_`, `ddnu_`, `omc_`,
  `domc_` all match to the bit. → the ported background splines feed the
  conductivity identical inputs.

- **Conductivity adaptive sample grid identical.** `DBGCOND`: the polynomial
  adaptive grid places `dimx = 911` sample points with `xt_sum` **bit-identical**
  (5.1867800045448821e+04). The sampled conductivity `yt` differs (below).

- **The confluent-hypergeometric 1F1 kernel is not the dominant driver.**
  `calc_Imn_array` → `hypergeometric1f1_cont_fract_1_modified_0_ada` computes the
  small modified `1F1m(1,b,z)` with `b ≈ z`, so `1F1 - 1 - z/b` cancels to
  ~1e-11 of its operands and the low bits depend on the exact rounding of each
  complex division. The port's `f_re` differs from the oracle by ~2e-12 for the
  first call. Two rounds of experiments:
    - Rewriting the division naively (explicit libstdc++-style `cdiv`/`cmul`, or
      `-fcx-limited-range`) made the conductivity `yt` mismatch **worse**
      (33231/262368 vs 25433 with the plain intrinsic). This is because
      gfortran's default intrinsic complex division and libstdc++'s
      `std::complex operator/` **both** lower to the same libgcc `__divdc3`
      range-reduced routine, so the plain intrinsic is already the closest match;
      naive division is a regression. Reverted.
    - Even with matching division, the ill-conditioned result carries
      codegen-dependent low bits (a small standalone unit inlines the division
      and matches the oracle for a given input; the full library object emits the
      `__divdc3` call and shifts the cancelled bits). Isolating the kernel made
      `f_re` bit-exact for the probed input but left `f_im` 1 ULP off and did
      **not** change EB.dat/poy — so 1F1 is at most a minor contributor.

## SUPERSEDED: ill-conditioned 1F1 as sole root cause

The 1F1 low-bit analysis below is correct and motivated routing the modified
form through fortnum on both sides. Its conclusion (that the 1F1 bits alone
set the row count and that exclusion was the only option) was wrong: with the
1F1 unified, the row-count gap persisted until the electron-mass, vector-cos,
and renorm-reconstruction fixes above landed.

Bisecting the conductivity pipeline (multi-agent workflow + manual, dumping each
stage in the port and a fresh oracle worktree at 3c364ab8, hex, bit-for-bit):

- Conductivity **samples** feeding the spline (`cp%K`) already differ. Walking
  up: `W2` (calc_W2_array_gen) differs at the first sample; its inputs
  `collmod/x1/x2/nu` are bit-identical but `Imn` (calc_Imn_array) differs; and
  inside `calc_Imn_array` the 1F1 **inputs are bit-identical** while its **output
  differs**:
  - `b = 40D7BCF4D4A45489 (−26889.34i)`, `z = 40D7BCB4D4A45489` — identical.
  - oracle `F_re = 3F1A70E4EA8CC000`, port `F_re = 3F1A70E4EA8D0000` (~2e-12).

- The call is `hypergeometric1f1_cont_fract_1_modified_0_ada` with `b ≈ z`
  (|z/b| ≥ 0.1). The raw continued fraction `S2` is bit-identical port/oracle;
  the divergence is entirely in the modified form
  `1F1m = (S2 − 1 − z/b)(b/z)((b+1)/z) − 1`, whose subtraction cancels to ~1e-11
  of its operands. High-precision (mpmath) shows the true value is
  `≈ 3F1A70E4EA8C8840`; both the oracle (`8CC000`) and the port (`8D0000`) are
  ~10^4 ULP of cancellation garbage away from it, in different directions.

- The exact garbage bits are **codegen-dependent** and not reproducible across
  compilers. Verified: gfortran intrinsic in a small unit → `8CC000`, in the
  library object → `8D0000`; explicit libstdc++-formula `cdiv`/`cmul` →
  `8CC000` (f_re) but `f_im` 1 ULP off; even the **C++ standalone** (`...B58F`)
  differs from the **C++ oracle full build** (`...B591`). No arithmetic form,
  `-fcx-limited-range`, `-ffp-contract`, `-O` level, module split, or forced
  inlining reproduces the oracle full-build's exact result.

- This ~2e-12 amplifies through the `Imn` formula's own denominator
  cancellation to ~1e-9 in the conductivity tensor `cti/cte`, then through the
  CVODE integration to ~5e-6 near the resonant layer, changing the adaptive
  output-grid row count (2682 vs 2648). **The physical values agree to ~1e-9
  (well inside rtol=1e-7); the golden fails only on shape (row count).**

Conclusion: reproducing the pre-port C++ output for `EB.dat`/`poy_test_err`
byte-for-byte requires reproducing a specific compiler's rounding of a
catastrophically ill-conditioned hypergeometric evaluation — a fundamental
cross-compiler floating-point limit, not a porting bug.

### Tested: keeping the 1F1 kernel as C++ does NOT close it

Compiling just the three routines (`modified_0_ada`/`kummer_modified_0_ada`/
`cont_fract_1_inv_ada`) as a small C++ file, verbatim from the pre-port source,
made the FIRST 1F1 call byte-exact (`8CC000 B591`) but left `EB.dat` at 2682
rows. Reason: the ill-conditioned result depends on the whole compilation unit's
codegen, and a 3-function C++ file differs from the oracle's full
`hyper1F1.cpp`. Dumping all calls: of 3411 identical-input 1F1 calls, **79 still
diverge** (the most ill-conditioned), which propagate. Reproducing the oracle's
exact codegen needs the *entire* `hyper1F1.cpp` compilation unit — which depends
on **GSL** (`gsl_sum_levin_u` in the dead `..._accel` variant, `gsl_integration`
in `..._quad`), a dependency the port removed. So the C++-kernel route cannot
close it without restoring GSL and ~400 lines of C++.

### RETRACTED: exclusion of the two outputs

An earlier revision excluded `EB.dat` and `zone_0_poy_test_err.dat` from the
comparison via a `GR_EXCLUDE` mechanism, arguing the divergence was a
pre-existing chaotic property of the case (control experiment: a 1e-12
perturbation of the oracle's own 1F1m moved EB.dat by 1.6e-6). The control
experiment was sound but measured the wrong thing: it proved the case
amplifies kernel-value changes, not that the port's row-count difference was
unfixable. The actual defect was the renorm-reconstruction off-by-one (see
Status). The exclusion and its mechanism are removed; the comparison now runs
main's harness unmodified (`zone_*_poy_test_err.dat` is checked against an
absolute self-consistency bar there, per #172, because it amplifies last-bit
noise ~1e8x between any two runs).

## Earlier framing — the ~1e-9 divergence (superseded by the root cause above)

The system matrix `Dmat` from `calc_diff_sys_matrix` (legacy) already differs at
the first RHS evaluation (r=3.0): oracle `Dmat_sum = 6.8086399267114894e+03` vs
port `6.8086399267037505e+03` (~1.1e-9 rel). The background inputs are
bit-identical (above); the difference is in the conductivity tensor `cti`/`cte`
(`DBGCOEF3`): oracle `cti_sum = 8.3454239222024362e+14` vs port
`8.3454239221962688e+14` (~7e-11 rel). This ~1e-9 RHS difference integrates
through CVODE to a ~5e-6 field difference near the resonant layer, setting the
different adaptive-grid row count and the ~7.77e-9 Poynting-theorem residual
(oracle 1.9e-12).

The raw conductivity K-matrix samples `yt` differ (~25k of 262k components
> 1e-12, first at the boundary grid point) even with the plain intrinsic build,
i.e. the residual is **not** the complex-division algorithm. The sampling
routines (`eval_a_matrix`, `calc_dem_djmi_arrays`, `eval_fgi_arrays`,
`calc_w2_array`, `calc_d_array`, `calc_k_matrices`, `calc_Imn_array`) are
unchanged legacy Fortran compiled by the same gfortran in both builds and are
fed bit-identical background inputs, so — apart from the 1F1 kernel (ruled
minor) — they should be bit-identical. The next lead is therefore any remaining
**ported** helper reachable from the FLRE conductivity sampling (or a subtle
compiled-once-vs-recompiled ordering difference between the oracle's legacy
objects and the port's), to be isolated by dumping `calc_Imn_array` and
`calc_d_array` outputs at the first sample against the oracle.
