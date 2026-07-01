# FLRE golden-record divergence — investigation log

Reference case: `test/golden/kilca/cases/flre_m6n2` (m=6, n=2).
Oracle: C++ KiLCA built at the port branch's merge-base with `main`
(`aa92bfb5`, equivalently the pre-port tree `3c364ab8`), which is what the
golden record was generated with. Comparison tool: `test/golden/bin/gr_numcompare.py`
(rtol=1e-7, atol=1e-12).

## Status

58 files compared. **56 pass byte/float-exact.** 2 remain:

- `linear-data/m_6_n_2_flab_[1,0]/EB.dat` — DIFFER (shape: 2682 rows vs oracle 2648)
- `linear-data/m_6_n_2_flab_[1,0]/zone_0_poy_test_err.dat` — DIFFER (shape: 681 vs 647)

Both are FLRE zone-0 outputs. The row-count difference is downstream of the
adaptive output grid (`sparse_grid_polynom`), which keeps points based on the
integrated basis fields; the fields differ at the ~1e-9 level near the resonant
layer, so the thinning keeps a slightly different number of points.

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

## Open — the remaining ~1e-9 divergence

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
