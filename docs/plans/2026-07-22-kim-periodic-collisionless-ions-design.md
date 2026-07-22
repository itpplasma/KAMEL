# Collisionless Ions in the Forced-Periodicity Solver

## Objective

Extend `electrostatic_periodic` so `ion_collision_model='collisionless'`
selects the same hybrid kinetic model as the global electrostatic solver:
Fokker--Planck electrons and collisionless Krook/Hamiltonian ions. Unlike the
current global implementation, the periodic implementation will include the
collisionless-ion contribution to both charge and parallel current. The
cyclotron-harmonic truncation remains `m_phi=0`, matching the global charge
model.

## Theory contract

The mathematical source is `~/1_projects/PhD_thesis.pdf`, Chapter 13,
especially equations (13.53)--(13.59), (13.63), (13.67), (13.151), and
(13.153). The implementation uses the established KIM radial Fourier pair and
the periodic phase convention already fixed for the FP kernels. The four
continuous Fourier integrands are evaluated locally at each periodized
gyrocentre-grid point and assembled with the existing endpoint-exclusive
periodic trapezoidal rule.

For each ion species, define

- `b_plus` and `b_cross` with the existing two-wavenumber FLR convention;
- `sI0 = exp(-b_plus) I_0(b_cross)` and
  `sIm1 = exp(-b_plus) I_{-1}(b_cross)` using the existing overflow-safe
  implementation;
- `k_abs = sqrt(k_parallel**2 + epsilon**2)`;
- `k_pole = k_parallel + i epsilon`;
- `z0 = -(omega_E-omega)/(sqrt(2) v_T k_abs)`.

The charge-potential response uses the even denominator `k_abs`, consistent
with the existing global collisionless-ion prefactors. Signed inverse-
`k_parallel` factors in the magnetic drive and both current moments use the
causal pole `k_pole`. This keeps `z0` on the real axis, avoiding the unphysical
lower-half-plane exponential growth of the plasma dispersion function, while
retaining causal signed response factors.

## Architecture

Add `collisionless_fourier_kernel_m` beside the existing FP Fourier module.
It owns pure/local per-species `m_phi=0` collisionless cores for
`rho-Phi`, `rho-B`, `j-Phi`, and `j-B`, plus hybrid wrapper functions that:

1. retain the current FP electron core;
2. select FP or collisionless ion cores from `ion_collision_model`;
3. respect `turn_off_electrons` and `turn_off_ions`;
4. apply the existing radial Fourier phase and continuous-kernel
   normalization exactly once.

The new module reuses the public FLR argument and scaled-Bessel helpers from
`flr2_fourier_kernel_m`. This is a one-way dependency; the existing FP module
is not modified to depend on the collisionless implementation. The periodic
assembler calls the configured hybrid wrappers. The default
`ion_collision_model='FokkerPlanck'` path must reproduce the current matrices
without numerical change.

No collisionless susceptibility arrays are stored in `plasma_t`. The
collisionless kernels consume local primitive/derived quantities (`vT`,
`omega_c`, `lambda_D`, `A1`, `A2`, `ks`, `kp`, and `om_E`) and evaluate `z0`
and `Z(z0)` directly. Consequently, `ion_fp_collision_scale` and the computed
ion collision frequency cannot leak into collisionless ion kernels.

## Current-density behavior

For collisionless runs, reconstructed parallel current becomes

`j_parallel = j_parallel,e(FP) + sum_i j_parallel,i(collisionless)`.

Both ion-current kernels come from the first parallel-velocity moment of the
same `m_phi=0` distribution used for the ion charge kernels. In particular,
the Mathematica reduction must verify

`G_jB = sqrt(2) vT z0 G_rhoB`

for every ion species and Fourier pair. The implementation must not infer the
potential-driven current from a charge-kernel ratio; it evaluates the thesis
first-moment expression directly.

## Verification

The derivation will be retained as both a Wolfram Language source file and a
standalone LaTeX document. Tests will cover:

- Mathematica-generated complex oracle values for all four kernels;
- homogeneous, static, zero-FLR Debye and zero-current anchors;
- the exact `jB/rhoB` recurrence;
- finite values at `k_parallel=0` for positive epsilon;
- epsilon sensitivity;
- FP-path preservation and explicit hybrid species decomposition;
- independence from `ion_fp_collision_scale` in collisionless mode;
- finite, nonzero ion current in a low-cost periodic assembly and solve;
- an end-to-end `electrostatic_periodic` collisionless run.

The clean baseline has one unrelated, documented failure in
`test_rhs_balance` (`gamma_a_e` stale-`forces_nl` mismatch). That failure is
not part of this feature and will be reported separately from new regressions.
