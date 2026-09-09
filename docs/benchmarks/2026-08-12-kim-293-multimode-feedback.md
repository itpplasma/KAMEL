# KIM #293: multi-mode periodic feedback

The periodic adapter retains signed mode identity (`m`, `n`), solve status,
resonant radius, current-normalization metadata, and compact-transition state
per configured mode. `kim_run_for_all_modes` resets per-mode storage before
each batch and stores each response independently. Call `kim_update_profiles`
first to supply updated QL density, temperatures, safety factor, and radial
electric field; production `get_dql` performs those calls in that order when
refreshing the wave response. Rotation continues to enter the established QL
force balance and radial electric field; this change adds no kinetic rotation model.

The cylindrical Fourier convention has a rational surface at `q = -m/n`.
The locator supports increasing and decreasing safety factor, exact endpoints,
and selects the innermost crossing when several exist. It returns no resonance
for wrong-sign or out-of-range modes and for `n = 0`; periodic KIM rejects those
modes explicitly. Invalid or non-finite profile geometry is rejected. The legacy
nonperiodic point-charge placement remains separate from periodic localization.
A failed solve terminates the batch; nonperiodic resonance metadata stays zero.

`get_dql` adds the independently computed mode diffusion tensors incoherently
on the global radial grid. Local fields use one compact transition; quadratic
transport uses its square. The resolved transition stays within sampled support.
No kinetic kernels or compression response have been added. The existing
restriction to `B_parallel = 0` remains.

Regression coverage includes:

- Production-locator tests for signed modes, decreasing q, endpoints, multiple
  crossings, absent roots, and invalid geometry, plus adapter rejection tests.
- The original electromagnetic density-feedback regression, alongside periodic
  tests using solved electric fields and currents rather than prescribed Br.
- Two distinct periodic modes with separated and overlapping supports, reordered
  and resized mode lists, and separate updates of density, Te, Ti, Er, and q.
  Refreshed responses are compared with clean initialization from those profiles.
- Actual `get_dql` assembly of two overlapping modes compared with the sum of
  independent single-mode runs, including order and overall phase invariance.

These are stationary refresh and transport-assembly regressions, not a validation
of long time evolution or of every rational surface in a reversed-shear plasma.
