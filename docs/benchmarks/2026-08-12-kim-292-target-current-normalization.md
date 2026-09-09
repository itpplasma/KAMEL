# KIM #292: target-current normalization

For a positive `I_par_toroidal` target, the QL-Balance adapter solves every
periodic KIM mode with `Br=1+0i` across the local window. It restores the
configured `Br_boundary_re/im` immediately after the solve. Zero or complex
configured drive amplitudes therefore cannot contaminate this unit response.

On the native trusted layer `[r_resonance-dx_asis, r_resonance+dx_asis]`,

\[
I_{\rm unit}=\int J_\parallel r\,dr,\qquad
s_{\rm target}=\frac{c I_{\parallel,\rm tor}}{2\pi I_{\rm unit}}.
\]

The integration interpolates the cylindrical integrand `r*J_parallel` at
clipped cell boundaries. It includes every overlapping cell and requires
finite, strictly increasing coordinates and an interval contained in the
sampled grid. This is the native KIM core integral, before interpolation or
tapering; it differs from the legacy Gaussian-fit current diagnostic on the
global QL grid. The `c/(2*pi)` convention matches the established target-current
conversion. Targets use c=1 CGS; the raw integral uses full CGS without `2*pi`.

All linear response quantities, including potential, electric and magnetic
fields, species parallel currents and radial current, scale by the same
complex `s`. The ion tensor scales once by `abs(s)**2`; the electron and
selectable drift-kinetic ion coefficients are calculated from the scaled
fields. For this positive-target path, the later antenna factor is one.

When `I_par_toroidal <= 0`, the configured drive and manual antenna factor
retain their previous behavior. A unit-response benchmark in that mode
requires explicitly setting the drive to one and the antenna factor to one.

## Relaxation and guards

The following optional `BALANCENML` settings have defaults:

```fortran
kim_current_floor = 1.0d-30
kim_current_max_scale = 1.0d12
kim_current_relaxation = 1.0d0
```

Relaxation describes one stationary update from the unit amplitude:
`s = (1-alpha) + alpha*s_target`. It is not an update from a previous time
step. With `alpha < 1`, the achieved current generally differs from the
target; the actual residual is saved. `kim_current_max_scale` limits `abs(s)`
relative to the unit drive.

The numerical helper reports status 0 for success, 1 for invalid finite
configuration, 2 for the current floor, and 3 for non-finite or excessive
response. Non-finite configuration is rejected on input. Guarded response
failure directly assigns zero to every linear response array and the ion
tensor, rather than multiplying invalid values by zero. The background and
radial geometry remain unchanged. A warning and the saved guard status
identify suppression; it is not reported as successful normalization.

## Saved diagnostics

With HDF5 output enabled, the normal saved-profile cadence writes
`/<mode_group>/CurrentNormalization/<time_index>/mode_<index>/`.
Each record contains:

- `target_current`, complex `unit_current`, `achieved_current`, `residual`,
  `scale`, and scalar `relative_residual`;
- `relaxation`, `current_floor`, `max_scale_ratio`, `status`, `m`, `n`, and
  `core_bounds`;
- native `r`, `unit_jpar`, and `normalized_jpar` profiles for independent
  integration of the reported current.

Complex scalars use `[real,imag]`; current profiles have Fortran shape
`(nrad,2)`. The achieved current includes `2*pi/c`, so its units match the
target. Relative residual is `abs(achieved-target)/max(abs(achieved),abs(target))`.
Guarded failures record zero achieved current and a nonzero residual, along
with their status; a failed raw unit profile may contain NaNs. Manual-drive
runs create no new target-normalization records. The time group has an
`active` flag: one for a positive-target snapshot, zero if a manual solve
rewrites that index. Repeated writes clear all old per-mode datasets first;
empty groups may remain because the HDF5 wrapper cannot delete groups.
Fresh solves also discard stale in-memory records.

## Regression coverage and limits

`test_periodic_current_normalization` checks clipped and within-cell bounds,
complex current, second-order quadrature convergence, geometry rejection,
and finite-value/scale guards. `test_periodic_response_normalization` checks
all response quantities, phase, doubling, relaxation, and clean suppression
of NaN/Infinity and multiplication overflow.

`test_periodic_target_current` runs the actual solver for two modes. It
checks target doubling (linear quantities double, both species' tensors
quadruple), zero/complex configured-drive invariance, restored configuration,
manual antenna behavior, independent integration of native HDF5 currents,
and persisted current-floor failure. The native integration tolerance is
`1e-10` relative; a separate `1e-2` tolerance covers interpolation onto the
401-point global QL grid. This is a stationary solver/transport integration
test, not a long time-evolution validation.

The supported periodic solver currently has `B_parallel=0`. It rejects a
nonzero compression drive because its linear charge and parallel-current
response kernels are not implemented. This PR does not enable a
transport-only magnetic insertion or claim to complete issue #292's vacuum
compression-response requirement. That physics requirement remains open.
