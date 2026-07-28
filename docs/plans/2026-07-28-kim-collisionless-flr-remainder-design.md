# Collisionless Finite-FLR Remainder Design

## Problem

The periodic collisionless ion `rho-Phi` core evaluates only the
`m_phi = 0` contribution from Markl, Eq. (13.58). At finite Larmor radius,
that contribution contains `-exp(-b_plus) I_0(b_cross) / lambda_D^2` in the
homogeneous static limit. The nonzero cyclotron harmonics supply the analytic
remainder derived in Eqs. (13.125)--(13.127):

```text
[exp(-b_plus) I_0(b_cross)
 - exp(-rho_L^2 (kr-krp)^2 / 2)] / lambda_D^2.
```

Without it, the collisionless kernel disagrees with both the global kernel and
the periodic Fokker--Planck convention even before gradients or resonances are
introduced. Existing tests set `rho_L = 0`, where the missing term vanishes.

## Design

Add the closed-form nonzero-harmonic remainder directly to `rho_phi` after the
zeroth-harmonic moment is evaluated. Reuse the already computed scaled Bessel
factor and compute the Gaussian from `rho_L`, `kr`, and `krp`. Keep the common
`1/(8*pi^2)` normalization explicit. Do not alter `rho_B`, `j_phi`, or `j_B`:
the same derivation shows the higher-harmonic charge response driven by `B_r`
reduces to the zeroth harmonic, and this issue concerns only `rho-Phi`.

## Regression Strategy

Replace the degenerate zero-FLR-only static check with independent finite-FLR
tests covering:

1. diagonal `kr = krp`, where the complete static response must be exactly
   `-1/lambda_D^2` for arbitrary nonzero `ks*rho_L`;
2. off-diagonal `kr != krp`, where the response must be the canonical Gaussian;
3. zero FLR, retained as a limiting-case guard;
4. non-static manufactured-oracle values, updated from an independent formula
   that includes the analytic remainder;
5. magnetic and current outputs unchanged by the correction.

The finite-FLR diagonal case is the essential mutation test: removing the
remainder changes the result by `1-exp(-b)I_0(b)` and must fail decisively.
