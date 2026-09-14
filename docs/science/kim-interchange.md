# KIM interchange conventions

**Status: evidence draft for A06 maintainer review.** This document records what
the implementation currently does. It is not approval for a collaborator
importer, a unit conversion, or a comparison claim. Source references below
refer to the checkout at commit `131b635b` unless a later source revision is
named.

## Scope and run types

The Fortran factory accepts `electrostatic`, `electrostatic_periodic`,
`electromagnetic`, `flr2`, `flr2_benchmark`, and `WKB_dispersion`
(`KIM/src/general/KIM_mod.f90:7-39`). The stable Python configuration exposes
only `electrostatic_periodic`, `electrostatic`, and `flr2`
(`KIM/python/src/kim/config.py:33-39`). The convention status below therefore
distinguishes the three Python-supported run types from Fortran run types that
are outside the current Python contract.

## Fourier and phase conventions

The forced-periodicity implementation defines the radial Fourier grid as

\[
 k_\ell = {2\pi\ell\over L},\qquad \ell=-M,\ldots,M,
\]

and reconstructs a field with

\[
 \delta\Phi(r)=\sum_{\ell=-M}^{M}\Phi_\ell\exp(+i k_\ell r).
\]

Both statements are executable conventions: `periodic_solve_m::solve_periodic`
sets `k_m = 2*pi*m/L` and the diagonal operator to `-k_m**2`
(`KIM/src/electrostatic_poisson/periodic_solve.f90:63-73`), while
`reconstruct_delta_phi` uses `exp(com_unit*k_m*r_out)`
(`KIM/src/electrostatic_poisson/periodic_solve.f90:200-222`). The derivative is consequently

\[
 {d\delta\Phi\over dr}=\sum_\ell (+i k_\ell)\Phi_\ell e^{+i k_\ell r}
\]

(`KIM/src/electrostatic_poisson/periodic_solve.f90:156-173`).

The two-wavenumber kinetic kernel uses the corresponding source/observation
phase

\[
 \exp[-i(k_r-k'_r)r_g]/(8\pi^2),
\]

in `hatG_rho_phi`, `hatG_rho_B`, `hatG_j_phi`, and `hatG_j_B`
(`KIM/src/asymptotics/flr2_fourier_kernel.f90:96-123`, `KIM/src/asymptotics/flr2_fourier_kernel.f90:125-151`,
`KIM/src/asymptotics/flr2_fourier_kernel.f90:153-185`, `KIM/src/asymptotics/flr2_fourier_kernel.f90:187-217`). The periodic assembly integrates this kernel over one
period with `K=(2*pi/L) int G dr_g`, implemented as an endpoint-exclusive,
equal-weight sum `2*pi/N` (`KIM/src/electrostatic_poisson/periodic_assembly.f90:8-29`,
`125-131`). A comparison implementation must preserve the endpoint convention
and the `1/(8*pi^2)` kernel normalization; neither may be inferred from an
FFT library default.

For the cylindrical field representation, the periodic solver derives

\[
 E_r=-\partial_r\Phi,\qquad E_\theta=-i\,m\Phi/r,\qquad
 E_z=-i\,n\Phi/R_0
\]

from `EBdat%Phi` (`KIM/src/electrostatic_poisson/poisson_periodic.f90:580-590`).
These signs are evidence for the positive spatial phase convention. The source
does not, in the files inspected here, state a complete time-domain sentence
such as `exp(-i omega t)`; frequency sign must therefore be carried as the
signed input `omega` until a maintainer approves a time convention.

## Mode geometry, frequency, and resonance

`parallel_wavenumber` and `perpendicular_wavenumber` define

\[
 k_\parallel=h_\theta {m\over r}+h_z {n\over R_0},\qquad
 k_s=h_z {m\over r}-h_\theta {n\over R_0},
\]

(`KIM/src/background_equilibrium/wavenumber_geometry.f90:15-29`). Here
`h_z=B_{0z}/B_0` and `h_theta=B_{0theta}/B_0`; those ratios are computed in
`calculate_equil` (`KIM/src/background_equilibrium/calculate_equil.f90:146-164`).
The same routine computes

\[
 \omega_E=-c E_r k_s/B_0
\]

(`KIM/src/background_equilibrium/wavenumber_geometry.f90:31-35`, with `c=sol` from
`KIM/src/util/constants_mod.f90:7-15`). Thus a resonance of the cylindrical
Fourier mode is \(k_\parallel=0\), which reduces to

\[
 q(r_\mathrm{res})=-m/n
\]

when the equilibrium relation used by `calculate_equil` is applied. The
Fortran resolver computes exactly `qres = -real(m_mode)/real(n_mode)` and
linearly interpolates the first radial crossing in radial order
(`KIM/src/grid/prepare_resonances.f90:33-56`). For `type_br_field == 2`, the
legacy nonperiodic solvers then replace the root by half the profile radius;
this override is explicitly skipped for `electrostatic_periodic`
(`KIM/src/grid/prepare_resonances.f90:58-61`). A comparison case must record whether it
is comparing the physical root or this legacy point-charge placement.

Frequency enters the collisionless response through

\[
 z_0=-{\omega_E-\omega\over
 \sqrt{2}\,v_T\sqrt{k_\parallel^2+\epsilon^2}},
\]

as implemented by `Krook_collisionless_z0`
(`KIM/src/kernels/Krook_kernel_plasma_prefacs.f90:34-51`). The radial-current
harmonic detuning is `om_E + ell*omega_c - omega`
(`KIM/src/asymptotics/radial_current_fourier_kernel.f90:267-305`). The signed
pole is `k_parallel+i epsilon`, while the magnitude regularizes even factors
(`KIM/src/kernels/Krook_kernel_plasma_prefacs.f90:5-32`, `KIM/src/asymptotics/collisionless_fourier_kernel.f90:292-349`).
Replacing signed `omega`, `k_parallel`, or `omega_c` by an absolute value is
therefore not an approved interchange transformation.

## Units and coordinate handling

The numerical constants are CGS: speed of light `sol=29979245800`, electron
mass and proton mass in grams, charge in statcoulomb, and `ev` in erg/eV
(`KIM/src/util/constants_mod.f90:7-13`). The equilibrium writes `B0`, `B0z`,
and `B0theta` in G (`KIM/src/background_equilibrium/calculate_equil.f90:214-223`)
and writes density in `1/cm^3` (`KIM/src/background_equilibrium/species_mod.f90:789-795`).
The force-balance implementation documents and computes

\[
 E_r={T_i eV\over e n}{dn\over dr}+{eV\over e}{dT_i\over dr}
       +{r B_0 V_z\over c q R_0},
\]

with \(T_i\) in eV and output in statV/cm
(`KIM/src/background_equilibrium/profile_input_m.f90:380-414`). Missing
`Er.dat` is filled by this `k=0`, no-poloidal-rotation path
(`KIM/src/background_equilibrium/profile_input_m.f90:252-271`, `KIM/src/background_equilibrium/profile_input_m.f90:327-378`).

The stable Python model declares `btor` in G, `major_radius` in cm, `frequency`
in `1/s`, boundary `Br` in G, and the signed mode numbers as dimensionless
(`KIM/python/src/kim/config.py:259-289`). Its profile validator labels radius
as cm, density as `1/cm^3`, temperatures as eV, `q` as dimensionless, and
`Er` as statV/cm (`KIM/python/src/kim/profiles.py:166-184`). The Fortran reader
consumes two columns and associates rows by position; it imposes quasineutral
ion densities from electron density and configured ion charge
(`KIM/src/background_equilibrium/species_mod.f90:1358-1444`). These are
implementation facts, not permission to convert an external source.

`coord_type='r_eff'` uses effective-radius profiles. `coord_type='sqrt_psiN'`
requires an equilibrium file or GEQDSK and runs the profile preprocessor
(`KIM/src/background_equilibrium/profile_input_m.f90:26-57`, `119-183`).
`auto` classifies a maximum first-column value above 2 as `r_eff`, otherwise
`sqrt_psiN` (`KIM/src/background_equilibrium/profile_input_m.f90:59-117`). This heuristic is unsuitable as
an external-format contract until a maintainer approves it.

## Python resonance check: evidence and status

The current stable Python validator does **not** use an absolute value. It sets
`target = -m_mode / n_mode` and finds a sign-changing or exact crossing in
`q-target` (`KIM/python/src/kim/profiles.py:88-123`, `224-251`). This agrees
with `kim_prepare_resonances` for the three Python-supported run types when
the run uses the physical-root path and has a scalar signed `q` profile.

An absolute-value check does exist in the legacy profile generator:
`python/utility/create_parabolic_profiles.py:16-17` interpolates `m_mode/n_mode`
against `abs(q)`. That helper is not the stable KIM Python validator and does
not establish a valid interchange convention. The Fortran source also uses
`abs` for the collisionless magnitude and for some numerical guards, but keeps
the signed `q=-m/n` root and signed pole. Consequently:

| Run type | Resonance evidence | Status of `abs(q)` / `abs(m/n)` check |
| --- | --- | --- |
| `electrostatic` | Cylindrical `k_parallel` and `kim_prepare_resonances`; legacy `type_br_field=2` override applies | **Incomplete/unsafe**: can locate a root with the wrong sign and cannot model the point-charge override. Use signed `q=-m/n`; maintainer must approve any broader rule. |
| `electrostatic_periodic` | Same signed root; periodic window is centered on `r_res`, then scaled by active-species `rho_L(r_res)` (`KIM/src/electrostatic_poisson/poisson_periodic.f90:500-549`) | **Wrong as a general rule**: the periodic window requires the physical signed root; absolute-value matching can center the window on a different surface. |
| `flr2` | Run calls the same initialization and uses signed `m`, `n`, `B0z`, and `omega` in `solve_flr2_response` (`KIM/src/flr2/flr2_run_type.f90:106-116`) | **Insufficient evidence for an absolute-value rule**: the inspected FLR2 response does not define an alternate absolute-value resonance criterion. Retain signed inputs pending maintainer review. |

The stable API currently requires nonzero `m_mode` and `n_mode`
(`KIM/python/src/kim/config.py:329-335`) and rejects a missing signed crossing
(`KIM/python/src/kim/profiles.py:241-251`). Whether nonresonant runs should be
supported, and how reversed shear or multiple crossings should be selected,
remain scientific decisions.

## Interchange decisions still requiring approval

The following are intentionally unresolved and must be approved before A07–A11
add conversion or comparison behavior:

- the physical time-phase convention and whether external frequency values may
  be sign-transformed;
- whether an external collaborator and dataset may be named, redistributed,
  or only accessed through restricted, maintainer-supplied fixtures;
- one permitted source format and its exact field names, coordinate declaration,
  unit declaration, and equilibrium provenance;
- whether only `r_eff` is accepted initially, or whether the
  `sqrt_psiN`/GEQDSK preprocessing path is part of the contract;
- treatment of multiple signed (q=-m/n) crossings, reversed shear, absent
  crossings, and `type_br_field=2` legacy placement;
- approved signed phase, field, and unit transformations, with independent
  reference values;
- whether the comparison uses complex `Phi`, complex `jpar`, a derived
  observable, or an integrated quantity, and the matching radial domain.

No collaborator, dataset, redistribution permission, conversion, tolerance, or
scientific equivalence is approved by this evidence draft.
