# KIM #291: drift-kinetic limiting benchmark

The benchmark compares the integral ion tensor with the established
Heyn/Markl drift-kinetic tensor. It is diagnostic-only: QL-Balance continues
to use the explicitly selected production model. Electrons remain on the
drift-kinetic path in all cases.

## Limit

The oracle takes one ordered wave pair with

\[
\ell=0,\qquad k_r=k_r'=0,\qquad k_s\rho_i,k_s'\rho_i\to0,
\qquad B_\parallel=0.
\]

For the Gaussian--Bessel moments,

\[
W_0\to1,\qquad W_1\to1,\qquad W_2\to2.
\]

The potential insertion is converted to the physical perpendicular field by
\(\Phi=iE_s/(c k_s)\) before taking the limit. The resulting four entries
are exactly the existing drift-kinetic ledger, including the antisymmetric
phase term in \(D_{12}-D_{21}\). A nonzero \(B_\parallel\) regression remains
in `test_periodic_ion_tensor` and is intentionally outside this limiting
case.

## Automated coverage

- `KIM/tests/test_dqli_limit_benchmark.f90` checks the exact zero-FLR branch
  and a finite small-FLR point.
- `KIM/mathematica/verify_quasilinear_integral_transport.wl` proves the
  Gaussian--Bessel limiting moments and records the reduced four-entry ledger.
- `kim_qldiff_m::calc_dqli_limit_benchmark` exposes a side-by-side diagnostic
  result with absolute and normalized residuals; it does not alter evolution.

The tolerances are separated by source: the exact algebraic branch is tested
at `2e-12`, the finite-FLR point at `2e-7` (quadrature/discretization), and
the existing Mathematica fixture comparisons retain their `1e-13`--`1e-35`
symbolic/high-precision gates. This keeps finite-FLR asymptotic error distinct
from algebraic and numerical integration error.
