# Roadmap: periodic KIM, integral ion transport, and QL-Balance coupling

Date: 2026-08-11

Branch: `roadmap/periodic-kim-ql-balance-gyrokinetic`

Status: design and implementation plan; no production implementation yet

## 1. Objective

Replace the global-KIM path used by QL-Balance with the faster local, forced-periodicity KIM calculation and add the quasilinear transport model derived in KIM's integral formalism.

The intended species split is:

- electrons: retain the existing Heyn/Markl drift-kinetic coefficients;
- ions: use the new gyrokinetic integral coefficients;
- magnetic perturbations: include both radial \(B^r\) and parallel \(B_\parallel\) channels;
- coupling: KIM computes the local resonant response and ion coefficients, while QL-Balance embeds those local results into its global radial profiles and advances them in time.

The first implementation must also provide an exact limiting case in which the new ion coefficients reduce to, and are numerically benchmarked against, the old drift-kinetic coefficients.

## 2. Executive design decision

Treat periodic KIM as a **local response service**, not as a global wave solver.

For each mode and resonant surface, periodic KIM receives equilibrium profiles plus prescribed complex magnetic drives \(B^r\) and \(B_\parallel\). It returns the local window, electrostatic response, physical electric-field components, species-resolved parallel current, and the ion \(2\times2\) quasilinear transport tensor. QL-Balance then applies a smooth compact transition to place the result on its global radial grid, computes the electron coefficients with the existing drift-kinetic model, sums independent mode contributions, and advances the global profiles.

This separates three responsibilities that are currently mixed together:

```text
KAMEL configuration and orchestration
                  |
                  v
QL-Balance global profiles and time integration
       |          ^
       | local    | local fields, current, ion D_ij
       v          |
periodic KIM response at each (m,n) resonance
       |
       +-- Phi, E, B^r, B_parallel, j_parallel,sigma
       +-- gyrokinetic ion D11, D12, D21, D22
```

The production local-to-global policy will be a smooth compact taper to zero. Raw extrapolation remains an optional diagnostic only. A hard cutoff is not acceptable because it creates discontinuities in the transport profiles.

## 3. Current state of the code

### 3.1 Repository and active work

The repository is a monorepo containing KIM, KiLCA, QL-Balance, common Fortran code, Python interfaces, and golden tests. This roadmap branch was created from `origin/main` at `b60f1127` in an isolated worktree, leaving the user's dirty main checkout untouched.

Relevant pre-existing development must be reconciled rather than independently reimplemented:

- `feat/kim-periodic-cyclotron-harmonics` contains work toward a finite ion cyclotron-harmonic sum in the periodic Fourier kernels;
- `feature/kim-periodic-dqle22` contains a partial periodic diagnostic for the old electron \(D_{22}\), not the requested ion tensor;
- the user's main checkout contains uncommitted periodic/global test changes and must not be overwritten or used as an integration base without review.

### 3.2 KIM periodic mode

The periodic calculation in `KIM/src/electrostatic_poisson/poisson_periodic.f90` currently:

- constructs a local Fourier window around one resonant surface;
- accepts a single prescribed complex constant \(B^r\) drive;
- solves the forced-periodic Poisson equation for \(\Phi\);
- computes total and species parallel-current diagnostics;
- does not return physical electric-field components through the public solver result;
- does not return \(B^r\) through the public solver result;
- does not represent or calculate \(B_\parallel\);
- does not calculate the requested four ion diffusion coefficients.

`KIM/src/general/kim_solver_m.f90` exposes a result type designed around the global solver. It has fields for \(E\), \(B^r\), \(j_\parallel\), and \(\Phi\), but lacks \(B_\parallel\), the transport tensor, and periodic-window metadata. The QL adapter would currently dereference periodic result arrays that are not allocated.

`KIM/src/background_equilibrium/periodic_background.f90` caches the global profiles and mutates process-global plasma/grid state to the periodic window. This is unsafe for repeated in-process solves over several modes and time steps unless every solve restores the prior state. State isolation is therefore a prerequisite for QL-Balance coupling.

The susceptibility calculation temporarily has moments \(\mathcal I^{p,q}\) for \(0\le p,q\le3\), but only a selected subset is retained. The integral transport tensor requires the full moment set for every retained cyclotron harmonic.

### 3.3 Existing KIM--QL-Balance adapter

`QL-Balance/src/base/kim_wave_code_adapter.f90` currently hard-codes electromagnetic global KIM, interpolates global KIM fields to the full QL grid, and applies a vacuum continuation outside the KIM plasma boundary. That continuation solves a different problem from embedding a local periodic window and cannot be reused as-is.

The adapter already supports profile feedback from QL-Balance to KIM. This is a useful seam and should be retained. However:

- it assumes global `r_field`/`r_plasma` semantics;
- it sets `Bp`, the equilibrium-parallel perturbation component, to zero;
- it loads only vacuum \(B^r\) from the KiLCA vacuum calculation;
- it has no local-window validity or transition metadata;
- it has no path for KIM-computed ion diffusion coefficients.

`QL-Balance/src/base/wave_code_data_64bit.f90` already runs a KiLCA vacuum calculation when KIM is selected. This can initially provide normalized magnetic drive shapes for both \(B^r\) and \(B_\parallel\), but the latter must be captured rather than discarded.

### 3.4 Existing quasilinear coefficients

`QL-Balance/src/base/diff_coeffs.f90::calc_transport_coeffs_ornuhl` implements the existing local drift-kinetic Heyn/Markl coefficients. `QL-Balance/src/base/get_dql.f90` currently calls this same routine for electrons and ions behind a hard-coded branch. The more elaborate routine in `QL-Balance/src/base/Wmn.f90` is dormant.

The current model uses \(E_s\) and \(B^r\), but not \(B_\parallel\). A misalignment-field term is constructed in `get_dql.f90`, although its addition to the coefficient is presently disabled. The new coupling must make the treatment of the misalignment contribution explicit and testable.

`KIM/src/diagnostics/kim_qldiff_mod.f90` is only a scalar port of the old electron \(D_{22}\) expression. It is not a foundation for the full ion integral tensor except as a regression reference.

### 3.5 Existing time evolution and shielding normalization

`QL-Balance/src/base/time_evolution.f90` already follows the broad sequence required by the new model:

1. compute wave response and diffusion coefficients;
2. derive an antenna/current normalization;
3. rescale diffusion coefficients;
4. update profiles;
5. repeat with the new profiles.

`QL-Balance/src/base/transp_coeffs_mod.f90::compute_antenna_factor_from_Ipar` reduces the current normalization immediately to a non-negative squared factor and applies it only to \(D_{ij}\). It therefore loses the signed or complex field amplitude and never updates \(B^r\), \(B_\parallel\), \(E\), or \(j_\parallel\) consistently. The new implementation needs an explicit per-mode magnetic-amplitude state.

## 4. Formal model to implement

The derivation is documented in the companion mono-manuscript chapter `documentation/chapters/quasilinear_integral_transport.tex`, with symbolic checks in `mathematica/verify_quasilinear_integral_transport.wl` in the KIM documentation repository. The implementation must preserve the master quadratic form

\[
D_{\sigma,ij}(r_g)=\operatorname{Re}\sum_{a,b}
\int\frac{dk'_r}{2\pi}\frac{dk_r}{2\pi}
\left(F^a_{k'_r}\right)^*
\mathscr D^{ab}_{\sigma,ij}(r_g;k'_r,k_r)F^b_{k_r},
\]

with

\[
F=(\Phi,B^r,B_\parallel)
\]

and

\[
\mathscr D^{ab}_{\sigma,ij}=
\frac{e^{i(k_r-k'_r)r_g}}{2\nu_\sigma}
\sum_\ell\sum_{p,q=0}^{3}
P^{ab}_{\ell,ij,pq}\,\mathcal I^{p,q}_{\ell\sigma}.
\]

The production implementation therefore contains four transport coefficients, nine ordered field-pair channels, the configured cyclotron-harmonic sum, and the complete required moment set. It must not collapse \(B_\parallel\) into an effective electrostatic term prematurely.

In the long-wavelength representation, the magnetic-moment contribution appears through

\[
v_m^r=\frac{cE_{m\perp}+uB_m^r
-i k_s v_\perp^2 B_{m\parallel}/(2\Omega_\sigma)}{B_0}.
\]

This expression fixes the relative phase and species-sign conventions for the \(B_\parallel\) channel. The Mathematica notebook/script is the algebraic oracle; generated Fortran coefficient tables must be checked against it and reviewed rather than transcribed manually.

## 5. Alternatives considered

### A. Make periodic KIM impersonate the global solver

KIM would extend every local field to the complete QL grid and return the existing global result type.

This minimizes adapter changes but hides the validity domain, encourages accidental use of periodic edge artifacts, and gives KIM responsibility for QL's global transport policy. Do not use this design.

### B. Return local fields and let QL-Balance calculate both species' coefficients

This keeps all transport in one executable, but it moves the full two-wave-number integral kernel and KIM susceptibility moments across the library boundary. It duplicates KIM conventions and produces a large, fragile interface.

This is acceptable only as a temporary independent reference implementation, not as production architecture.

### C. Return a typed local response and ion tensor from KIM

KIM owns the periodic Fourier representation, susceptibility moments, and gyrokinetic ion contraction. QL-Balance owns the global embedding, electron drift-kinetic coefficients, mode accumulation, and time integration.

This is the recommended design. It keeps each calculation next to its native data, makes the species policy explicit, and permits a small stable coupling contract.

## 6. Target interfaces

Introduce a periodic request/result contract at the KIM library boundary. Names below are illustrative; implementation may follow the repository's naming conventions.

```fortran
type :: kim_periodic_request_t
    integer :: m, n
    integer :: ell_min, ell_max
    complex(dp) :: Br_amplitude
    complex(dp) :: Bparallel_amplitude
    character(len=...) :: Bparallel_shape
    ! Window, transition, and benchmark controls.
end type

type :: kim_periodic_result_t
    real(dp), allocatable :: r(:)
    real(dp), allocatable :: taper_coordinate(:)
    logical, allocatable :: trusted_core(:)
    complex(dp), allocatable :: Phi(:), Phi_m(:)
    complex(dp), allocatable :: Er(:), Etheta(:), Ez(:), Es(:), Ep(:)
    complex(dp), allocatable :: Br(:), Bparallel(:)
    complex(dp), allocatable :: jpar_species(:,:)
    real(dp), allocatable :: D_ion(:,:,:) ! (2,2,n_r)
    real(dp) :: r_resonance, dx_asis, dx_transition, period
    integer :: status
end type
```

Required contract properties:

- all complex quantities use one documented \(e^{i(m\theta+n\varphi-\omega t)}\) convention;
- grid endpoint semantics and Fourier normalization are explicit;
- the trusted as-is core is distinct from the periodic transition zone;
- no module-global KIM state remains changed when a library call returns;
- an error or rejected solve cannot leave a partially updated QL mode;
- \(B^r\) and \(B_\parallel\) amplitudes and normalization are reported in output metadata;
- the tensor indices are unambiguous: \(D_{11},D_{12},D_{21},D_{22}\).

## 7. Local fields, global profiles, and misalignment

### 7.1 Construct physical local fields first

Derive \(E_r\) spectrally from the periodic \(\Phi_m\) coefficients. Derive angular/toroidal components algebraically with the same mode convention, then transform to \(E_s\) and \(E_p\). Do not reuse the global-grid derivative routine without first proving that its basis and endpoint convention match the periodic DFT.

Apply gauge checks to \(\Phi\), but couple QL-Balance through physical fields and transport coefficients rather than through an extrapolated potential.

### 7.2 Production local-to-global transition

For every resonance define a common real weight \(w(r)\):

- \(w=1\) in the periodic as-is core;
- \(w\) is a compact \(C^\infty\), or at minimum \(C^2\), transition in the trusted buffer;
- \(w=0\) outside the local support.

The initial implementation should use a standard bump/smoothstep with an analytically fixed endpoint behavior. The transition width is a configuration value bounded by the KIM window metadata.

The taper is applied to already-derived physical perturbation components. Do **not** form \(w\Phi\) and then differentiate, because the artificial \(w'\Phi\) term would create a spurious radial electric field. Use the same weight for the cancelling components of the misalignment velocity so that their phase relation is not altered.

Since the transport tensor is quadratic in perturbation amplitude, either recompute it from the tapered fields or, for a local KIM tensor, use

\[
D^{\mathrm{global}}_{ij}(r)=w(r)^2D^{\mathrm{local}}_{ij}(r).
\]

This preserves positive-semidefinite structure under the embedding. Independent mode contributions are accumulated incoherently as in the current QL-Balance model.

### 7.3 Diagnostic alternatives

Provide two non-default diagnostic modes:

- `hard_zero`, solely to expose cutoff sensitivity;
- `extrapolate`, with a bounded extension distance and explicit warnings.

Neither is an acceptance path for production. A transition-width/extension sensitivity scan is required before the periodic coupling is considered validated.

## 8. Magnetic drive and shielding-current evolution

### 8.1 First-step constant-psi drive

At the first accepted time step, initialize each mode with a normalized complex amplitude \(s_{mn}=1\) and prescribe \(B^r\) as constant across the local KIM window. This is the constant-\(\psi\) approximation requested for the local calculation.

The initial \(B_\parallel\) shape should come from the existing KiLCA vacuum calculation at the same resonance and be normalized consistently with \(B^r\). Both components are multiplied by the same shielding amplitude unless a later physical closure establishes an independent response. A configuration option may set \(B_\parallel=0\) for the required drift-limit benchmark, but production must not silently discard it.

The periodic electrostatic model does not yet solve Ampere's equation, so neither magnetic component is claimed to be self-consistent inside KIM. They are prescribed drives whose amplitude is closed by the shielding-current update.

### 8.2 Explicit amplitude state

Replace the present diffusion-only `antenna_factor` concept by a persistent per-mode amplitude \(s_{mn}\). Given a unit-amplitude local solve,

\[
B^r=s_{mn}\widehat B^r,\qquad
B_\parallel=s_{mn}\widehat B_\parallel,\qquad
E=s_{mn}\widehat E,\qquad
j_\parallel=s_{mn}\widehat j_\parallel,
\]

and

\[
D_{ij}=|s_{mn}|^2\widehat D_{ij}.
\]

Linear scaling permits one KIM solve per mode and profile state. A second KIM call after normalization is unnecessary unless a future nonlinear model is introduced.

### 8.3 Current closure per accepted time step

For mode \((m,n)\):

1. restore the latest accepted QL profiles;
2. build the KIM local equilibrium and unit magnetic drives;
3. solve periodic Poisson and calculate local \(E\), \(j_\parallel\), and ion \(D_{ij}\);
4. integrate the shielding current over the trusted current layer with documented geometry and CGS/SI conversions;
5. update \(s_{mn}\) from the new shielding current;
6. scale all fields and currents by \(s_{mn}\), and both species' coefficients by \(|s_{mn}|^2\);
7. taper local quantities onto the global QL grid;
8. advance QL profiles and accept or reject the time step;
9. commit the new amplitude only when the time step is accepted.

Under the target-current interpretation already present in QL-Balance, the unit-response update is

\[
|s_{mn}^{\star}|=
\frac{|I_{\parallel,mn}^{\mathrm{target}}|}
{2\pi\,|\widehat I_{\parallel,mn}^{\mathrm{shield}}|},
\]

with the exact \(2\pi\) and unit conversion verified against the current implementation. For an iterative/non-unit evaluation, use

\[
s_{mn}^{k+1}=s_{mn}^{k}
\frac{I_{\parallel,mn}^{\mathrm{target}}}
{I_{\parallel,mn}^{\mathrm{shield},k}}
\]

with a documented phase policy. Apply under-relaxation, a current floor, finite-value checks, and a configurable maximum scale ratio. On failure, retain the last accepted amplitude and reject or reduce the QL time step.

The old scalar `antenna_factor` may remain as a backward-compatible output equal to \(|s_{mn}|^2\), but it must not cause a second coefficient scaling.

### 8.4 Deferred electromagnetic-response alternative

If “rescale \(B^r\) according to the shielding current” instead means an electromagnetic response relation such as

\[
B^r=B^r_{\mathrm{vac}}+\mathcal G I_\parallel^{\mathrm{shield}},
\]

then the response operator \(\mathcal G\), its phase, and its radial/mode normalization must be supplied or derived. This is different from the current `I_par_toroidal` target-current normalization and is explicitly deferred from the first implementation.

## 9. Species policy

Make the model selection explicit in configuration and output:

```text
electron_transport_model = drift_kinetic
ion_transport_model      = kim_integral
```

The production path must reject unsupported combinations rather than fall through hard-coded `.true.`/`.false.` branches.

Electrons use the existing Heyn/Markl drift-kinetic routine evaluated from the same tapered physical local fields. Ions use the KIM integral tensor. In benchmark mode, QL-Balance also evaluates the old ion drift-kinetic tensor but does not use it for evolution unless explicitly requested.

## 10. Implementation milestones

Each milestone ends with a small reviewable commit and must pass its listed gate before the next dependent milestone starts.

### M0. Freeze conventions and reconcile prerequisites

Files:

- `roadmap.md`
- `KIM/docs/` or the repository's chosen design-document location
- `KIM/src/asymptotics/flr2_fourier_kernel.f90`
- relevant changes from `feat/kim-periodic-cyclotron-harmonics`

Tasks:

- document phase, Fourier, radial-wavenumber, charge, \(\Omega_\sigma\), and field-component conventions;
- compare the harmonic branch against current main and extract a clean tested change;
- decide the supported \(\ell\) truncation and convergence control;
- record the exact mapping from KiLCA `Bp` to KIM/QL \(B_\parallel\);
- decide whether target-current normalization is the intended shielding closure.

Gate:

- convention tests on one analytic Fourier mode;
- harmonic \(\ell=0\) result unchanged from main;
- no dependence on a dirty working tree or uncommitted branch state.

Suggested commit: `docs: freeze periodic transport conventions`

### M1. Isolate periodic KIM state and expose the local response

Files:

- `KIM/src/general/kim_solver_m.f90`
- `KIM/src/background_equilibrium/periodic_background.f90`
- `KIM/src/electrostatic_poisson/poisson_periodic.f90`
- `KIM/src/electrostatic_poisson/periodic_solve.f90`
- `KIM/src/electrostatic_poisson/fields_mod.f90`
- `KIM/tests/CMakeLists.txt`
- new `KIM/tests/test_periodic_result_contract.f90`
- new `KIM/tests/test_periodic_state_restore.f90`

Tasks:

- add typed request/result data or extend `kim_results_t` without ambiguous allocation states;
- return local window metadata, \(\Phi\), Fourier coefficients, and physical \(E\) components;
- return the prescribed \(B^r\) and \(B_\parallel\) profiles;
- snapshot/restore the global plasma and grid state around every periodic solve;
- support sequential solves for different modes and updated profiles in one process;
- ensure failure paths restore state.

Tests first:

- a periodic result contains all required arrays with consistent bounds;
- spectral differentiation of a single Fourier mode is exact to tolerance;
- solve A, solve B, then solve A gives the original A result;
- the global background is byte/numerically unchanged after success and injected failure;
- scaling both magnetic drives by a complex scalar scales \(E\) and \(j_\parallel\) linearly.

Gate: the public KIM library can run a multi-mode, multi-profile periodic sequence without leaked state or unallocated result access.

Suggested commits:

1. `test: specify periodic KIM result contract`
2. `feat: expose isolated periodic KIM responses`

### M2. Retain the full susceptibility moment set

Files:

- `KIM/src/background_equilibrium/species_mod.f90`
- `KIM/src/asymptotics/flr2_fourier_kernel.f90`
- periodic susceptibility assembly modules
- new `KIM/tests/test_integral_moments.f90`

Tasks:

- store or expose \(\mathcal I_{\ell\sigma}^{p,q}\) for \(p,q=0,\ldots,3\);
- index the values by species, radial point, and retained harmonic;
- avoid adding dozens of individually named scalar arrays; use a dense bounded tensor or accessor;
- retain symmetry/conjugation identities as runtime-debug assertions and unit tests;
- measure memory use for production window and harmonic counts.

Gate: every moment used by the Mathematica-derived polynomial tables is available with tested index and normalization conventions.

Suggested commit: `feat: expose periodic susceptibility moments`

### M3. Implement the gyrokinetic integral ion tensor, including \(B_\parallel\)

Files:

- new `KIM/src/transport/quasilinear_integral_m.f90`
- new `KIM/src/transport/quasilinear_field_operators_m.f90`
- generated/read-only coefficient data under `KIM/src/transport/generated/`
- copy of `mathematica/verify_quasilinear_integral_transport.wl` under `KIM/mathematica/`
- a generation/export script and checked output manifest
- `KIM/CMakeLists.txt` and `KIM/tests/CMakeLists.txt`
- new `KIM/tests/test_quasilinear_integral_algebra.f90`
- new `KIM/tests/test_quasilinear_integral_properties.f90`

Tasks:

- construct the Gaussian--Bessel weights/operators for the retained harmonics;
- implement all four \(D_{ij}\) entries and all nine ordered pairs of \(\Phi,B^r,B_\parallel\);
- contract the field vectors with the complete moment set;
- explicitly Hermitianize the quadratic form before taking the real part;
- expose per-channel diagnostics in debug/test builds;
- keep Mathematica as a build-time/development oracle, not a runtime dependency;
- attach a generator version/hash to the produced coefficient data.

Tests first:

- generated tables exactly match the Mathematica export;
- each of the 36 coefficient-polynomial families is checked at rational/random high-precision points;
- every one of the nine field-pair insertions is activated alone and compared with a direct reference contraction;
- \(D_{12}=D_{21}\) within the derived Onsager tolerance;
- the symmetric tensor is positive semidefinite within roundoff for admissible inputs;
- phase rotations of the complete field vector leave \(D\) invariant;
- setting \(B_\parallel=0\) removes only its five associated ordered/cross contributions;
- the analytic long-wavelength \(B_\parallel\)-only result is recovered.

Gate: the ion tensor passes algebraic, Hermiticity, Onsager, positivity, and harmonic-convergence checks independently of QL-Balance.

Suggested commits:

1. `test: add Mathematica transport oracle`
2. `feat: implement integral ion diffusion tensor`
3. `feat: include parallel magnetic transport channel`

### M4. Establish the old-model limiting benchmark

Files:

- new `KIM/tests/test_integral_drift_limit.f90`
- new cross-code driver under `test/ql-balance/`
- `QL-Balance/src/base/diff_coeffs.f90` only if a pure/reference interface is needed
- benchmark fixtures under `test/golden/cases/`

Limit to impose:

- \(k_\perp\rho_i\rightarrow0\);
- retain \(\ell=0\);
- set \(B_\parallel=0\);
- use the diagonal/local radial-wave-number limit \(k'_r=k_r\);
- match collision frequency, thermal-speed, field, phase, and unit conventions exactly.

Benchmark ladder:

1. symbolic Mathematica reduction of all four new coefficients to the old formulas;
2. one-radius Fortran comparison using synthetic moments and fields;
3. constant-profile periodic-window comparison point by point;
4. a QL-Balance dual-evaluation run that writes old/new ion tensors and relative residuals;
5. an end-to-end limiting golden case in which evolution uses the new ion tensor while electrons remain on the old model.

Acceptance tolerances must be separated into algebraic roundoff, quadrature/DFT error, and finite-FLR truncation error. Do not hide discrepancies in one broad tolerance.

Gate: every \(D_{ij}\) converges to the old ion result at the predicted order, with a checked report stored as test output/artifact.

Suggested commit: `test: benchmark integral tensor in drift limit`

### M5. Couple periodic KIM to QL-Balance and globalize the local response

Files:

- `QL-Balance/src/base/kim_wave_code_adapter.f90`
- `QL-Balance/src/base/wave_code_data_mod.f90`
- `QL-Balance/src/base/wave_code_data_64bit.f90`
- `QL-Balance/src/base/get_dql.f90`
- new `QL-Balance/src/base/local_response_embedding_mod.f90`
- `QL-Balance/src/base/diff_coeffs.f90`
- `QL-Balance/tests/test_kim_wave_code_adapter.f90`
- new `QL-Balance/tests/test_local_response_embedding.f90`

Tasks:

- select `electrostatic_periodic` rather than electromagnetic/global KIM;
- send each mode's local request and current QL profiles to KIM;
- capture vacuum `Bp` as the initial \(B_\parallel\) shape alongside \(B^r\);
- interpolate only within the local support, preserving complex phase;
- implement the shared smooth transition weight and \(w^2D\) embedding;
- calculate old drift-kinetic electron coefficients from the same tapered local fields;
- receive and use the integral KIM ion tensor;
- replace hard-coded transport branches by validated species configuration;
- expose optional old/new ion dual evaluation without changing the selected evolution model;
- keep per-mode quantities separate until the established incoherent accumulation point.

Tests first:

- partition checks: \(w=1\) in the core, all endpoint derivatives required by the chosen smoothness vanish, and \(w=0\) outside support;
- a constant local tensor maps to \(w^2D\);
- no derivative-of-window electric field appears;
- overlapping resonances add diffusion coefficients but do not coherently add unrelated field phases;
- electrons call only the drift-kinetic model and ions call only the integral model in production configuration;
- the adapter handles several modes and profile refreshes without stale KIM state;
- hard-zero and extrapolation modes require explicit diagnostic configuration.

Gate: a stationary QL-Balance solve uses periodic KIM for all configured modes and produces finite, localized electron and ion coefficient profiles including \(B_\parallel\).

Suggested commits:

1. `test: specify local-to-global response embedding`
2. `feat: couple periodic KIM responses to QL-Balance`
3. `feat: split electron and ion transport models`

### M6. Add constant-psi shielding-amplitude evolution

Files:

- `QL-Balance/src/base/transp_coeffs_mod.f90`
- `QL-Balance/src/base/time_evolution.f90`
- `QL-Balance/src/base/SingleStep.f90`
- `QL-Balance/src/base/paramscan.f90`
- restart/checkpoint modules identified during implementation
- `QL-Balance/src/base/writeData.f90`
- new `QL-Balance/tests/test_shielding_amplitude.f90`
- new `QL-Balance/tests/test_shielding_restart.f90`

Tasks:

- add accepted/trial per-mode complex amplitude state;
- initialize \(B^r\) with constant-\(\psi\) unit amplitude at the first step;
- integrate the updated shielding current only over the trusted local layer;
- update and under-relax the amplitude with current-floor and ratio guards;
- scale \(B^r,B_\parallel,E,j_\parallel\) by \(s\) and \(D\) by \(|s|^2\);
- remove/detect double application of the old antenna factor;
- roll amplitude back when a time step is rejected;
- serialize amplitude, phase policy, and normalization metadata in restart files;
- report target current, calculated shielding current, residual, amplitude, relaxation, and guard activation per mode/time step.

Tests first:

- known linear unit response gives the analytic amplitude in one update;
- doubling target current doubles fields/current and quadruples \(D\);
- a rejected step restores both profiles and amplitudes;
- restart continuation is identical to uninterrupted evolution;
- current below the floor fails safely without `NaN` or an unbounded field;
- relaxation converges monotonically for a synthetic linear response;
- `SingleStep`, parameter scans, and full time evolution use the same normalization path.

Gate: a multi-step run updates profiles, recomputes shielding current, rescales both magnetic components consistently, and restarts reproducibly.

Suggested commits:

1. `test: specify shielding amplitude evolution`
2. `feat: evolve constant-psi magnetic amplitude`
3. `feat: persist shielding state in restarts`

### M7. Complete KAMEL configuration, output, and golden coverage

Files:

- `KIM/src/setup/config_mod.f90`
- KIM configuration reader and `KIM/nmls/KIM_config.nml`
- `KIM/nmls/README.md`
- QL-Balance control/configuration readers
- `QL-Balance/namelists/balance_conf.nml`
- `python/balance_interface/balance_conf.py`
- `QL-Balance/src/base/writeData.f90`
- top-level `Data.md` and relevant README files
- new `test/golden/cases/kim_periodic_integral_transport/`

Configuration to expose:

- periodic run type and window controls;
- harmonic bounds/convergence policy;
- electron and ion transport models;
- \(B_\parallel\) source/shape and benchmark-zero switch;
- local-to-global transition kind and width;
- shielding-current closure, relaxation, current floor, and scale bounds;
- benchmark/dual-evaluation switches.

Output/provenance to add:

- local-window and transition metadata per mode;
- \(B^r\), \(B_\parallel\), and relevant \(E\) components;
- species-resolved \(j_\parallel\);
- unscaled and scaled electron/ion \(D_{ij}\);
- per-channel integral diagnostics when requested;
- old/new benchmark residuals;
- shielding target/current/amplitude history;
- model names, convention version, Mathematica-table hash, and KIM revision.

Gate: a Python-configured end-to-end golden run is reproducible and its output is sufficient to diagnose every normalization step.

Suggested commits:

1. `feat: configure periodic integral transport`
2. `feat: record periodic transport provenance`
3. `test: add periodic KIM QL-Balance golden case`

### M8. Scientific validation and performance acceptance

Validation matrix:

- Fourier grid convergence \(M\) and radial grid convergence;
- periodic as-is width and transition-width sensitivity;
- cyclotron-harmonic convergence;
- finite-Larmor-radius to drift-limit convergence;
- \(B_\parallel=0\), \(B^r=0\), \(\Phi=0\), and every two-field cross channel;
- positive-semidefinite and Onsager residuals across a representative profile scan;
- comparison of periodic versus global KIM in the central trusted layer for cases where both are valid;
- branch-cut/causal-prescription checks near resonance;
- shielding-amplitude convergence and time-step refinement;
- multi-resonance overlap behavior;
- electron results unchanged relative to the established drift-kinetic path;
- wall-clock and memory comparison against global KIM.

Acceptance criteria:

- all formal limiting identities pass at their specified tolerances;
- transition choices do not materially affect core observables within an agreed uncertainty band;
- no transport is produced outside compact support except by an explicitly selected diagnostic mode;
- periodic coupling is demonstrably faster than global KIM for the target workloads;
- output contains enough provenance to reproduce the result;
- no regressions occur in existing KiLCA and global-KIM configurations.

Suggested commit: `test: validate periodic integral transport workflow`

## 11. Test strategy and continuous integration

Use a red--green--refactor sequence for each milestone. The test pyramid should contain:

- **symbolic oracle tests:** Mathematica export and exact identities;
- **pure Fortran unit tests:** weights, operators, moments, contractions, amplitude updates;
- **KIM integration tests:** local field/current/tensor response and state restoration;
- **QL adapter tests:** interpolation, tapering, species dispatch, mode accumulation;
- **restart/time tests:** accepted/rejected state and reproducibility;
- **golden tests:** one small end-to-end periodic case and one drift-limit case;
- **scientific convergence jobs:** larger non-blocking or scheduled tests with stored reports.

Never use the global-periodic comparison alone as proof of correctness: the formulations and valid domains differ. Use it as one diagnostic after algebraic and limiting tests have passed.

## 12. Risks and mitigations

| Risk | Consequence | Mitigation |
|---|---|---|
| Periodic KIM mutates global plasma state | mode/time-order-dependent results | snapshot/restore immediately; later refactor toward immutable local context |
| Inconsistent Fourier/phase convention | wrong cross terms, especially \(B_\parallel\) | freeze conventions, analytic one-mode tests, generated algebra oracle |
| Taper creates a false electric field | artificial transport at window edge | derive physical fields before tapering; never differentiate tapered \(\Phi\) |
| Current normalization is applied twice | diffusion scales as \(|s|^4\) | single explicit amplitude owner and scaling-invariance tests |
| Near-zero shielding current | unbounded \(B^r\) update | floor, bounded update, relaxation, rollback, diagnostic failure |
| Lost complex phase in current closure | incorrect shielding response | preserve complex amplitudes until closure; document magnitude-only option |
| \(B_\parallel\) convention differs between KiLCA and KIM | incorrect magnetic-moment transport | explicit component mapping and vacuum-field unit test |
| Full moment tensor increases memory | slower local solver | bounded dense layout, window-local allocation, memory benchmark |
| Feature branches contain uncommitted or overlapping work | silently lost fixes | extract reviewed commits, do not merge dirty worktrees wholesale |
| Benchmark “agreement” obtained with loose tolerance | hidden algebra error | separate roundoff, discretization, and asymptotic tolerances |

## 13. Definition of done

The project is complete when:

- QL-Balance uses periodic KIM rather than global KIM for the configured coupled path;
- periodic KIM safely supports repeated modes and profile updates in one process;
- KIM returns local \(E\), \(B^r\), \(B_\parallel\), species current, window metadata, and the complete ion \(D_{ij}\);
- all nine \((\Phi,B^r,B_\parallel)\) field-pair channels are included and symbolically verified;
- electrons use the established drift-kinetic coefficients and ions use the new integral coefficients;
- local results are embedded into global QL profiles with a validated smooth compact transition;
- the first-step constant-\(\psi\) drive and subsequent shielding-current amplitude update work across accepted steps and restarts;
- the new ion tensor reproduces the old ion tensor in the documented drift limit;
- golden, restart, convergence, and regression tests pass;
- configuration and output fully identify the model, normalization, transition, and algebra version;
- the periodic coupled workflow meets an agreed speedup target relative to global KIM.

## 14. Resolved shielding-closure decision

Target-current normalization was approved on 2026-08-11. M6 will preserve the existing `I_par_toroidal` interpretation and choose the common \(B^r\)/\(B_\parallel\) amplitude so that the integrated shielding current matches that target.

An electromagnetic response closure of the form \(B^r=B^r_{\mathrm{vac}}+\mathcal G I_\parallel^{\mathrm{shield}}\) is a possible later extension. It is not part of this implementation chain.

## 15. Published implementation chain

The roadmap is split into AFK-ready GitHub issues in dependency order:

1. Existing prerequisites: #231 (signed resonances) and #258 (cyclotron harmonics).
2. #285 returns a state-isolated periodic KIM local response.
3. #286 delivers a \(\Phi\)-only integral-ion \(D_{11}\) tracer bullet.
4. #287 completes the \((\Phi,B^r)\) ion tensor.
5. #288 carries the KiLCA \(B_\parallel\) vacuum drive through periodic KIM.
6. #289 completes all \(B_\parallel\) transport channels.
7. #290 couples one periodic mode to QL-Balance with compact global embedding.
8. #291 benchmarks the integral ion tensor against the drift-kinetic limit.
9. #292 normalizes one mode to the target shielding current.
10. #293 adds multi-mode periodic response with profile feedback.
11. #294 evolves and restarts the shielding amplitudes.
12. #295 ships configuration, provenance, golden coverage, scientific validation, and performance reporting.

Issues #285--#295 carry the `ready-for-agent` label. The final production acceptance in #295 remains a human scientific sign-off after the automated implementation and validation report are complete.
