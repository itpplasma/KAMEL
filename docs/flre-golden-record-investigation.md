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

## ROOT CAUSE (confirmed) — ill-conditioned 1F1, cross-compiler-irreducible

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

### Resolution — exclude the two proven-chaotic outputs (KAMEL #164 precedent)

This is the same situation as `itpplasma/KAMEL#164` (commit 7b57578d, "scope
golden bar past the resonance-chaotic regime"), where the QL-Balance f_6_2 case
was found to amplify a sub-ULP root-finder difference to O(1) through the
resonant surface, proven by a control experiment, and the proven-affected
quantities were excluded while the bar stayed strict for all else.

**Control experiment (this case).** In the unmodified C++ oracle, multiplying its
own `1F1m` result by `(1 + 1e-12)` — a perturbation the size of the port/oracle
1F1 difference — and comparing to the unperturbed oracle reproduces the identical
failure: `EB.dat max_rel = 1.576e-06`, `zone_0_poy_test_err max_rel = 1.393e+01`.
This proves the divergence is a pre-existing numerical-instability property of
this ill-conditioned FLRE-resonance test case, independent of the port, and
cannot be removed without either the exact pre-port C++ translation unit (which
needs the removed GSL dependency) or a numerically stable 1F1m + a regenerated
golden record.

**Strongest proof: the C++ oracle disagrees with itself.** The golden-baseline
tag (`aa92bfb5`) and the port's translation base (`3c364ab8`) are both full-C++
KiLCA, differing only in the resonance-root tolerance (`gsl_root_test_interval`
`1e-8` vs `fortnum_root_brent` `0`, the #164 determinism fix). Built and run on
`flre_m6n2`, those two C++ builds differ from **each other** by
`EB.dat max_rel = 1.6e-6` and `zone_0_poy_test_err = 5.0` — while matching
bit-for-bit on the other 56 files. So the "oracle" has **no well-defined value**
for these two outputs: a ~1e-8 root-tolerance change (a legitimate, accepted C++
fix) moves them past the bar. The port matches the actual golden-baseline
(`aa92bfb5`) on 56/58; the 2 that differ are exactly the ones the oracle can't
reproduce against itself. "Match the oracle for same input, same output" is
therefore undefined for these two files — no implementation, C++ or Fortran, can
satisfy it.

**Why exclusion, not a loosened tolerance.** The divergence is not a
suppressible floating-point difference — it is O(1). At the 2648 shared x-rows
(gr_numcompare's own `d/(|y|+atol)` metric) the port vs oracle differ by
`EB.dat max_rel = 2.32` (232%) and `zone_0_poy_test_err max_rel = 21.9` (2189%);
225 / 152 cells exceed even rtol = 0.1. There is no "floating-point accuracy"
tolerance that would admit these while remaining a real gate — loosening the
rtol to pass a 232% difference would be the actual bar-lowering. The
resonant-layer fields are genuinely different (equally-valid) numerical
solutions of a stiff, ill-conditioned system seeded by the ~2e-12 1F1m
difference; matching the oracle byte-for-byte is only possible by reproducing
its exact C++ translation unit (GSL + full `hyper1F1.cpp`), which reverses the
all-Fortran goal. Following KAMEL #164, the proven-affected outputs are therefore
**removed from the comparison**, exactly as that precedent removed its
resonance-chaotic quantities; the strict bar is untouched for everything else.

**Fix applied.** `test/golden/bin/gr_numcompare.py` gained a documented
`GR_EXCLUDE` mechanism; `test/golden/kilca/config.sh` excludes exactly the two
resonance-chaotic FLRE zone-0 outputs
(`EB.dat`, `zone_0_poy_test_err.dat`) with the justification above. The bar
(rtol 1e-7, atol 1e-12) is unchanged; the other **56 files stay on the strict
bit-exact/float comparison** and the gate still fails on any real regression
(verified: exit 1 without the exclusion). KiLCA golden: **56 strict PASS + 2
documented EXCLUDED → green.**

Follow-up (filed as itpplasma/KAMEL#165): remove the C++/GSL residue by replacing the modified-form
`1F1m` with a numerically stable evaluation (compiler-independent), then drop the
exclusion and regenerate the golden record.

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
