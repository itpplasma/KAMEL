# KIM Radial-Current Fourier Kernels Design

## Goal

Implement the three documented one-dimensional radial-current response kernels

\[
(B^r,B_\parallel,\Phi)\longmapsto j^r
\]

and use them to measure convergence with cyclotron-harmonic cutoff.

The source of truth is
`~/2_areas/KIM_documentation/documentation/chapters/radial_current_1d.tex`,
including its audited Fourier convention, output-first susceptibility indices,
signed cyclotron frequency, and zero-perpendicular-wavenumber limits.

## Scope

Add a dedicated `radial_current_fourier_kernel_m` module beside
`flr2_fourier_kernel_m`. It will expose:

- per-species, per-harmonic contributions for
  `jrad-Phi`, `jrad-Br`, and `jrad-Bparallel`;
- species-summed kernels with an optional harmonic cutoff;
- phase-stripped radial FLR coefficients for direct unit testing;
- symmetric harmonic accumulation from `-cutoff` through `+cutoff`.

The periodic solver will assemble and reconstruct

\[
j^r=K^{j^rB^r}B^r+K^{j^r\Phi}\Phi,
\]

with the presently unavailable `B_parallel` field set to zero. The
`jrad-Bparallel` kernel remains a tested public API for future magnetic-compression
support. The global real-space FEM kernel is not changed in this feature.

## Data model

The documented `jrad-Br` kernel needs susceptibility moments `I01` and `I03`.
KIM already computes a complete `symbI(0:3,0:3)` table but does not retain
`I03`. Add boundary and cell-centred `I03` arrays to `species_t`, allocate,
populate, write, periodize, and deallocate them wherever the existing moments
are managed.

The harmonic index uses KIM's existing signed `omega_c` convention and

\[
x_{2\ell}=-(\omega_E+\ell\Omega-\omega)/\nu.
\]

No absolute value may be applied to `omega_c` in the resonance or in the
remaining `jrad-Bparallel` prefactor. `rho_L` continues to use the positive
magnitude already stored by KIM.

## Radial FLR coefficients

For each species and harmonic, form the documented coefficients

\[
A_{0\ell}=\Gamma_\ell[A_1+(1-b_+-\ell)A_2]+\Lambda_\ell A_2,
\qquad
A_{2\ell}=\tfrac12 A_2\Gamma_\ell,
\]

then evaluate

\[
O_{j\ell}=D_\ell^o A_{j\ell},\qquad
W_{j\ell}=D_\ell^oD_\ell^s A_{j\ell}.
\]

The implementation uses analytic first and mixed derivatives of the scaled
modified-Bessel products. Neighboring orders provide derivatives of `I_l` and
avoid finite differencing. A small-`k_perp` branch uses the adjacent-order
limit required by the documentation so coordinate factors proportional to
`1/k_perp**2` never produce `NaN` or `Inf`.

## Kernel API and phase

The public summed functions accept source `kr`, response `krp`, grid index,
and optional cutoff, matching the argument order already used by
`flr2_fourier_kernel_m`. The radial Fourier phase therefore follows the
repository's existing transform orientation. The harmonic angular phase is
added consistently with the documented source/response gyroaverages after
mapping them to that argument order.

Each function validates that `0 <= cutoff <= mphi_max`; an invalid cutoff is a
programming error and terminates with a descriptive `error stop`. Species
turn-off flags and `artificial_debye_case` follow the existing Fourier kernels.

## Verification

Tests are independent of the implementation formulas wherever practical:

1. retained `I03` equals `symbI(0,3)` for every allocated harmonic;
2. analytic radial coefficients agree with centred numerical derivatives away
   from coordinate singularities;
3. at `k_perp=k_perp'=0`, one radial insertion has no ordinary zero-harmonic
   contribution and the double insertion is carried by `ell=+/-1`;
4. signed-frequency reversal has the documented species parity;
5. per-harmonic contributions sum exactly to each cutoff kernel;
6. all three kernels remain finite at zero and large FLR arguments;
7. existing Fourier and periodic tests remain green;
8. a representative KIM case reports shell-by-shell Frobenius norms and the
   change in reconstructed `jrad`.
