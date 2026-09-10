# KiLCA FLR mathematics verification

This document tracks executable checks of the differential KiLCA plasma-response model. A
passing golden record is useful for detecting change, but is not evidence that an equation is
correct.

Primary reference: Markus Markl, *Kinetic investigation of the bifurcation of resonant magnetic
perturbations in magnetically confined plasmas and development of an integral response model*
(2024), [doi:10.3217/efp2p-0x485](https://doi.org/10.3217/efp2p-0x485).

## Compiled model

The default CMake configuration selects:

- `MAGNETIC_DRIFT=MDNO`, so the compiled conductivity sources have no `_drift` suffix;
- `COLLISION_MODEL=FPGEN`, so the compiled collision source is
  `calc_W2_array_gen.f90`;
- `flre_order` values 1, 3, and 5, dispatched explicitly by
  `KiLCA/flre/conductivity/calc_D_array.f90`.

The corresponding generated coefficient sources are `calc_D_array_1_fac.f90`,
`calc_D_array_3_fac.f90`, and `calc_D_array_5_fac.f90`. Other similarly named source files are
not part of the default target and are not evidence for the compiled model.

## Verification status

### Verified

- Time/phase convention: thesis (4.1), with real fields represented by
  `Re(F exp(-i omega t))`. The complex amplitudes in `KiLCA/flre/quants/` are covered by the
  Lorentz-force CTest.
- Lorentz-force density: the full phasor average
  `1/2 Re(q n E* + j x B*/c)`. Thesis (6.17) gives the magnetic-force torque reduction. The
  production symbols are `calc_time_averaged_lorentz_force` and
  `calc_lorentz_torque_density`.
- Cylindrical torque arms: `T_theta = r F_theta` and `T_z = R0 F_z`, implemented by
  `calc_cylindrical_torque_density`.
- Radial integration: the signed numerical weighting for the straight-cylinder volume element
  `4 pi^2 R0 r dr`, implemented by `flre_quants::vol_fac`,
  `integrate_over_cylinder`, and `calc_lorentz_torque_on_cylinder`.

### Pending

- FLR expansion: thesis (5.5)-(5.7), a Taylor expansion in radial Larmor displacement through
  order `N`. The small parameter is locally `rho/L`, or `k_perp rho` for a Fourier scale.
  Coefficient identities and quantitative remainder/domain bounds remain to be checked in
  `calc_D_array_{1,3,5}_fac.f90`, `kmatrices.f90`, and `conductivity.f90`.
- Constant-psi reduction: thesis (6.18) adds assumptions beyond numerical radial integration.
- Resonant-layer width: `flre_zone::D` is read from zone settings and controls conductivity-grid
  refinement. No derivation is attached to this input, so it is not a computed physical layer
  width.
- Island width: no direct consuming KiLCA symbol was identified in this audit. Response-field
  output does not by itself verify an island-width formula.

Here `rho = v_T/omega_c`. The statement `rho/L << 1` is the asymptotic condition for the
Taylor expansion, not a repository-enforced validity threshold. No numerical bound for the
neglected `O((rho/L)^(N+1))` contribution has yet been established for the generated
coefficient files.

## Lorentz torque convention

For complex amplitudes with the thesis convention `exp(-i omega t)`, the period average of two
real harmonics is

`<Re(a exp(-i omega t)) Re(b exp(-i omega t))> = 1/2 Re(a b*)`.

KiLCA therefore computes, in CGS units,

`f = 1/2 Re(q n E* + j x B*/c)`.

The CTest uses deliberately complex amplitudes so that removing either conjugation changes the
answer. It also chooses basis-aligned `j` and `B` vectors so that a cross-product sign mutation
fails independently of the electric term.

After transforming `(r,s,p)` components to cylinder components `(r,theta,z)`, the local output
slots contain:

- radial force density `f_r` (retained for historical output compatibility);
- poloidal torque density `r f_theta`;
- toroidal torque density `R0 f_z`.

The radial slot is not integrated: `calc_lorentz_torque_on_cylinder` passes a zero volume factor
for that component. Tangential components use

`Integral T dV = 4 pi^2 R0 Integral r T(r) dr`

with a signed trapezoidal rule. Ion and electron torque densities are calculated separately and
then summed.

The focused test is registered as `kilca_lorentz_torque` in CTest.

## Constant-psi boundary

Thesis equation (6.18) is not merely the numerical radial integration above. It assumes one
perturbation mode and replaces all radially varying quantities except the parallel current by
their values at the resonant surface. Equation (6.19) additionally assumes strong,
electron-only shielding and neglects the perpendicular electric field; the thesis explicitly
states that this limit is invalid close to the gyrocenter resonance.

No production helper currently exposes that analytic reduction as an independently testable
function. Until such a consumer is identified or introduced, the constant-psi part of issue
#195 remains unverified.

## Reproduction

```bash
cmake -S . -B build -G Ninja
cmake --build build --target test_kilca_lorentz_torque
ctest --test-dir build -R kilca_lorentz_torque --output-on-failure
```
