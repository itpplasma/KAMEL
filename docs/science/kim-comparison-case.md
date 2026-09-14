# KIM comparison case

**Status: proposed comparison contract for A06 review.** No collaborator,
external dataset, redistribution permission, reference result, tolerance, or
scientific equivalence has been approved. This file separates implementation
evidence from decisions that require a scientific maintainer and, where input
rights are involved, the data owner.

## What KIM can currently expose

The stable Python API has a typed view for forced-periodicity output. It reads
`backs/e/r`, `fields/Phi`, `fields/jpar`, and
`setup/periodic_scale/dx_asis` (`KIM/python/src/kim/results.py:92-105`). The
Fortran periodic run writes complex `Phi`, total `jpar`, electron and ion
`jpar`, and radial `jrad` under `fields/`, with units statV, statA/cm², and
related CGS units (`KIM/src/electrostatic_poisson/poisson_periodic.f90:626-657`).
The Python reader decodes the HDF5 real/imag compound representation into
complex arrays (`KIM/python/src/kim/results.py:219-230`).

For a periodic result, the currently implemented derived observable is

\[
 I_\parallel(\mathcal R)=2\pi\int_\mathcal R r\,j_\parallel(r)\,dr,
\]

evaluated by a trapezoidal rule (`KIM/python/src/kim/results.py:149-168`). The
Fortran diagnostics routine uses the same expression over the grid supplied
by its caller (`KIM/src/diagnostics/kim_diagnostics_mod.f90:14-33`). The typed
Python view offers the regions `as_is` (the points within `dx_asis` of the
reconstructed resonance radius) and `full_window`; these labels describe code
behavior, not a scientific comparison choice.

The other stable run types write fields as well. `electrostatic` writes
potential and charge/current fields through
`KIM/src/electrostatic_poisson/poisson.f90:114-128` and
`KIM/src/electrostatic_poisson/poisson.f90:179-216`; `flr2` writes `Br`, `Phi`, and `jpar` in
`KIM/src/flr2/flr2_run_type.f90:118-123`. Their output grids, field equations, and
normalizations are not interchangeable with the periodic Fourier spectrum.
No comparison observable is selected for them here.

## Candidate case (not approved)

The most constrained candidate for maintainer review is one
`electrostatic_periodic` run with `type_br_field=12`, a single signed mode
pair `(m,n)`, a supplied `r_eff` profile set, and `Bparallel=0`. The reason for
this proposal is traceability: the periodic solver documents the constant
radial drive and its Fourier column (`KIM/src/electrostatic_poisson/periodic_solve.f90:43-89`),
and stores the reconstructed field on the endpoint-exclusive window
(`KIM/src/electrostatic_poisson/poisson_periodic.f90:542-590`). This is a
candidate boundary for discussion, not permission to use a dataset.

The case record would need these values before implementation is authorized:

| Item | Current evidence | Approval needed |
| --- | --- | --- |
| Mode/resonance | Signed `m_mode`, `n_mode`; `q=-m/n` | Sign and crossing policy |
| Radius | `r_eff` cm or `sqrt_psiN` preprocessing | Coordinate and equilibrium source |
| Profiles | Two-column KIM profile roles | Names, units, grid and missing fields |
| Equilibrium | `B0`, `B0z`, `B0theta`, `h_z`, `h_theta` | External-to-reduced mapping |
| Excitation | Complex constant `Br` in G; optional `Bparallel` | Complex phase/amplitude rule |
| Frequency | Signed `omega` in 1/s; FLR2 requires zero | Time phase and sign |
| Observable | Complex field or integrated current | One observable and domain |
| Reference | None selected | Collaborator, dataset and data rights |

The requested source formats are therefore open. A maintainer must select and
approve either (a) already permitted two-column profile text plus an approved
equilibrium artifact, (b) one named experimental format with a read-only
fixture, or (c) native KIM HDF5 as a solver-to-analysis reference only. No
external format, collaborator, or redistribution license is inferred from the
repository. If inputs cannot be redistributed, the approved process must say
how a maintainer-supplied private fixture is identified, staged, hashed, and
run without copying it into the repository.

## Evidence-backed reduction boundary

The only mappings that can be stated from current code are:

1. `r_eff` profile rows are read as radius in cm and values in KIM's CGS
   conventions. The Python validator requires matching grids for `n`, `Te`,
   `Ti`, `q`, and present `Vz` (`KIM/python/src/kim/profiles.py:88-123`).
2. The Fortran reader implements `n_i=n_e Z_i/sum(Z)`, despite the source
   comment calling this quasineutrality. For mixed charges, the charge sum
   generally differs from `n_e` (`KIM/src/background_equilibrium/species_mod.f90:1434-1444`).
   External multi-species mapping requires scientific review.
3. The equilibrium ODE integrates pressure work, computes signed `B0z`,
   `B0theta`, `B0`, `h_z`, `h_theta`, `k_s`, `k_parallel`, and `omega_E`
   (`KIM/src/background_equilibrium/calculate_equil.f90:96-168`).
4. For a periodic run, the resonant radius is found from signed `q=-m/n`, the
   active-species Larmor radius sets the window scales, and the solution is
   reconstructed with the positive Fourier phase (`KIM/src/grid/prepare_resonances.f90:33-61`,
   `KIM/src/electrostatic_poisson/poisson_periodic.f90:500-559`,
   `KIM/src/electrostatic_poisson/periodic_solve.f90:200-222`).

No source-to-equilibrium reduction for a collaborator's experimental
coordinate, magnetic geometry, perturbation phase, or reference observable is
established by those routines. A07–A11 must refuse undeclared conventions and
unsupported mappings rather than silently applying a presumed conversion.

## Missing physics and unsupported comparisons

The following limitations are part of the case boundary until a maintainer
changes it deliberately:

- The stable API does not expose Fortran-only `electromagnetic`,
  `flr2_benchmark`, or `WKB_dispersion` run types
  (`KIM/src/general/KIM_mod.f90:22-38`, `KIM/python/src/kim/config.py:33-39`).
- `flr2` currently requires `omega=0`, one ion species, and Fokker–Planck
  collisions (`KIM/src/flr2/flr2_run_type.f90:57-70`); it is not a frequency-
  sweep comparison case.
- The periodic solver's full kernel includes both radial and perpendicular
  wavenumber contributions; the optional global-matching approximation drops
  only the `k_s^2` contribution from the Bessel arguments. The kernel still
  depends on `k_r` and `k'_r`; it also sets electron FLR to zero
  (`KIM/src/asymptotics/flr2_fourier_kernel.f90:19-40`,
  `KIM/src/setup/config_mod.f90:62-65`). Which model is desired is open.
- Collisionless ions use a causal pole and a magnitude in different factors;
  an external model that supplies only `|k_parallel|` cannot be declared
  equivalent (`KIM/src/kernels/Krook_kernel_plasma_prefacs.f90:5-51`).
- KIM's profiles and equilibrium are reduced cylindrical quantities. The
  repository does not establish a general flux-surface, 3-D equilibrium, or
  experimental magnetic-perturbation reduction.
- The periodic current integral depends on its selected radial region, grid,
  complex phase, and (2\pi r) weighting. Comparing magnitudes only would
  discard information and requires explicit approval.

## Proposed analytic vectors (open test material)

These vectors test representation and comparison plumbing only. They are not
approved physics reference values and must not be used to bless a scientific
tolerance.

For the documented periodic basis, choose (L=2\pi), (M=1), ordered modes
((-1,0,+1)), and

\[
 \Phi_{-1}=1+2i,\quad \Phi_0=3-i,\quad \Phi_{+1}=-2+0.5i.
\]

The documented reconstruction predicts the exact algebraic checks

\[
 \Phi(0)=2+1.5i,\qquad \Phi(\pi)=4-3.5i,
\]

and

\[
 \partial_r\Phi(0)=1.5-3i.
\]

These values exercise both Fourier signs and the mode ordering. A second
plumbing vector may use (M=0), any positive (L), and
\(\Phi_0=2+3i\); reconstruction must return the constant (2+3i) at every
radius. An HDF5 fixture, its compound-complex encoding, and expected detached
array dtypes still require approval.

For a proposed observable test, use a synthetic complex current
\(j_\parallel(r)=J_0\) on a declared grid and compare the implementation's
trapezoidal (2\pi\int rj_\parallel dr) with the analytically integrated value
\(\pi J_0(r_b^2-r_a^2)\). The grid, region, and (J_0) are test choices, not
the external reference case.

## Tolerance questions

Before A07 or A11 records a pass/fail result, the scientific reviewer must
answer:

- Is agreement required for complex real and imaginary parts, magnitude and
  phase, or a derived real observable?
- Is error absolute, relative, or mixed near a zero crossing? What denominator
  and floor are approved?
- Are tolerances applied pointwise, in an L2 norm, or to the selected integral?
- Must the radial domains and grids match exactly, or may one approved
  interpolation be used? If interpolation is allowed, which direction and
  order?
- Are Fourier truncation, radial window, collision model, FLR switches,
  `omega`, and profile interpolation held fixed between reference and KIM?
- What evidence distinguishes a parsing/conversion error from a model
  disagreement?

## Proposed A07–A11 expansion after approval

These are concrete follow-on questions for the local implementation plan; they
are not authorization to implement them in A06.

**A07 — conventions and bounded conversions.** Define a versioned metadata
record containing source coordinate, radius unit, each profile unit, magnetic
field sign, mode signs, frequency sign/time phase, perturbation phase, and
equilibrium provenance. Approve only listed conversions. Required tests should
include the Fourier vectors above, one independently calculated CGS conversion,
missing metadata errors, and refusal of dimensionally incompatible fields.

**A08 — one input reader.** Select exactly one approved source format and one
fixture delivery method. Validate structure, required quantities, finite values,
monotone grids, and source metadata. Test malformed/truncated input and exact
known fixture values without downloading at test time.

**A09 — preparation.** Implement only the approved equilibrium/profile mapping,
excitation mapping, and domain restrictions. Stage source files, converted
profiles, request, hashes, and a conversion report in a fresh directory; refuse
overwrite and refuse extrapolation unless explicitly approved. Do not launch
KIM from preparation.

**A10 — public CLI and walkthrough.** Expose preparation and validation through
the approved API/CLI, preserve relative paths and report fields, and exercise
one approved case with exact commands. The walkthrough must show how each
conversion-report field is interpreted and how restricted inputs are supplied.

**A11 — comparison result.** Read only the approved reference representation,
select the approved observable and domain, report complex/real normalization,
and preserve source identity and model settings. Test incompatible metadata,
detached arrays, analytic values, and a known approved comparison result. Do
not add a generic dataset mapper or claim magnetic/torque agreement where KIM
does not expose that observable.
