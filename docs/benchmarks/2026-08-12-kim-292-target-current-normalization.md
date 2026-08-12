# KIM #292: target-current normalization

For periodic KIM, each mode is first solved with the unit constant-\(\psi\)
drive. On the trusted current layer, the adapter evaluates

\[
I_{\rm unit}=\int J_\parallel r\,dr,
\qquad
s=\frac{c I_{\parallel,\rm tor}}{2\pi I_{\rm unit}}.
\]

The linear response fields \((B_r,B_\parallel,E,J_\parallel)\) are multiplied
by the complex amplitude \(s\); the ion tensor is multiplied once by
\(|s|^2\). The legacy antenna-factor rescaling is disabled for this path, so
it cannot apply a second \(|s|^2\) factor. A zero target keeps the explicit
unit-response benchmark behavior.

`periodic_current_normalization_m` owns the trusted-layer quadrature, the
2-\(\pi\) convention, relaxation, current-floor, non-finite, and maximum
scale-ratio guards. The adapter records per-mode unit current, complex scale,
and guard status. `test_periodic_current_normalization` covers geometry,
target doubling, and the current-floor guard.
