# KIM Periodic Poisson Deviation Design

## Decision

Use the one-dimensional radial Poisson operator already defined by the KIM
integral-model hat formulation, not an empirical `-k_s^2` correction. The full
cylindrical harmonic Laplacian contains radial metric and helical terms, but
chapter 15 reduces the numerical integral model to the local radial coordinate
`x` and the hat solver implements that same reduction. The periodic symbol is
therefore `D_m=-k_m^2`. The Mathematica gate records both formulas and the
declared reduction so a future real-geometry extension cannot silently inherit
the slab approximation.

## Deviation Formulation

Let `L = d^2/dr^2 + 4 pi K^{rho Phi}` and
`g = -4 pi K^{rho B} B_r`. Choose a finite, smooth reference field
`Phi_ref(r)` that agrees with global `bc_type=3` at the two physical window
edges. A quintic C2 lift connects those values without evaluating singular
`1/k_parallel` expressions at resonance. Solve

`L psi = g - L Phi_ref`

in the periodic Fourier space and reconstruct the physical field as

`delta Phi(r) = Phi_ref(r) + psi(r)`.

`Phi_ref` must remain a physical-space lifting function during reconstruction.
If it is projected into the same finite Fourier space and added there, exact
linear algebra cancels the shift and reproduces the original direct-periodic
solution. Tests must include this no-op identity as a negative design check.

For the tracer bullet, the reference and its analytic radial second derivative
are sampled on the endpoint-exclusive periodic grid. `K Phi_ref` is evaluated
by the same Fourier kernel after projection. Unequal physical edge values
necessarily create a Fourier jump, so the declared residual measures Gibbs
contamination in the untouched central plasma region, not at the artificial
period boundary; residuals above 5% are rejected. The default benchmark gives
3.16%. Excessive residual is a configuration error, not a curve-fit tolerance.

## Operator and Gauge Policy

Screening makes the zero Fourier mode nonsingular. In a true vacuum/null-space
case, the mode-zero source must satisfy the compatibility condition and an
explicit mean-zero gauge is required; otherwise the solver returns a nonzero
status. No unconditional deletion of mode zero is allowed.

## Verification

The independent oracle covers Fourier measure and phase, hat overlap and weak
derivative matrices, radial versus full cylindrical symbols, `4 pi` source
sign, constant/nonconstant magnetic-drive coefficients, screened mode zero,
manufactured complex `Phi_ref`/`psi`, electric-field reconstruction, and
rho-B/j-Phi/j-B response signs. Fortran tests consume all 14 committed rows
without requiring Mathematica, recover a manufactured complex deviation,
exercise screened and vacuum zero modes, and reject an under-resolved lift.
The independent mutation gate checks source sign, transform factor, transpose,
boundary term, gauge removal, and `k_s^2`.

The first benchmark after implementation showed that this boundary correction
is not the dominant curve discrepancy: DC-removed periodic/global relative L2
remained 0.72443. That result is evidence for the next diagnostic step rather
than a reason to tune the lift to the target curve.
