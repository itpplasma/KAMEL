# Collisionless Causal Pole Design

The ion-only collisionless charge kernel must retain the rational-surface
singularity without letting the background-grid placement select an accidental
regularization.  The current implementation mixes an endpoint average of
`z0` with a denominator formed from the signed cell average of `k_parallel`.
When a background cell straddles the rational surface, the latter nearly
cancels while the former remains on the endpoint scale.

The collisionless branch defines two cell-centred forms with the same positive
width:

\[
\widetilde{k}_{\parallel,j} = k^{cc}_{\parallel,j} + i\epsilon_\parallel,
\qquad
|k|_{\epsilon,j}=\sqrt{(k^{cc}_{\parallel,j})^2+\epsilon_\parallel^2}.
\]

The causal complex pole is used only where the analytical expression contains
a signed `1/k_parallel`, currently the magnetic `kappa_rho_B` prefactor.  Terms
whose collisionless expression contains `1/abs(k_parallel)`, including `z0`
and `G1/G2/G3 rho_phi`, use the broadened magnitude.  Putting the complex pole
inside `z0` would send `z0` deep into the lower half-plane on one side of the
resonance; the analytic continuation of the plasma-dispersion function then
contains an exponentially growing contribution and overflows.  The split
prescription preserves the sign information exactly where it belongs without
changing the original magnitude terms into signed poles.

`epsilon_parallel` is a required positive namelist value for collisionless-ion
runs, expressed in inverse centimetres.  The FP and collisional Krook paths
retain their current real/absolute-value conventions.

The regulator does not remove the collisionless singularity.  At fixed
positive epsilon it defines the causal pole and makes the problem numerically
resolvable.  Decreasing epsilon recovers the singular response.  A meaningful
epsilon scan must increase radial resolution so the scale over which
`k_parallel` changes by epsilon is resolved.

Verification has three levels.  Unit tests check the pole sign, broadened
magnitude, finite plasma-dispersion evaluation at the AUG resonance,
cell-centred `z0`, and use by the appropriate collisionless prefactors.  Configuration tests
check that collisionless runs reject a non-positive epsilon.  The AUG benchmark
checks that, at fixed resolved epsilon, a small grid translation no longer
changes the kernel by orders of magnitude; an epsilon scan checks that the
singular response grows as epsilon decreases.  Electron kernels and FP outputs
remain unchanged by construction.
