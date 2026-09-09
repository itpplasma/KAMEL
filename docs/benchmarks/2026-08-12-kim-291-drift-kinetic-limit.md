# KIM #291: drift-kinetic limiting benchmark

This diagnostic compares the finite-Larmor-radius ion tensor with the existing
Heyn/Markl drift-kinetic tensor using the same physical electric and magnetic
fields. It does not replace the configured ion evolution model or the existing
drift-kinetic electron coefficients.

## Physical limit and comparison contract

The local limiting fixture takes one diagonal radial wave pair with

\[
\ell=0,\qquad k_r=k_r'=0,\qquad k_s\rho_i,k_s'\rho_i\to0,
\qquad B_\parallel=0.
\]

The electric field is in statV/cm and the magnetic perturbation in Gauss.
Since \(E_s=-i k_s\Phi\), the benchmark supplies
\(\Phi=iE_s/k_s\). The factor \(c\) belongs to the radial electric drift:
the legacy coefficients contain \(c^2|E_s|^2\) and electric/magnetic cross
terms proportional to \(c v_{Ti}\). The local API requires nonnegative
\(k_s\); at exactly zero, a nonzero supplied \(E_s\) is rejected because its
potential representation is singular. The exact zero-wave-number fixture is
therefore magnetic-only; the electric limit is approached at nonzero \(k_s\).

The Gaussian–Bessel moments reduce to

\[
W_0\to1,\qquad W_1\to1,\qquad W_2\to2.
\]

The current integral kernel symmetrizes its real diagonal-wave tensor. It
recovers the **symmetric part** of the legacy tensor in this limit. It does
not in general recover the legacy antisymmetric phase term. Define

\[
a={1\over 2\nu_i B_0^2}
  {2c v_{Ti}\,\operatorname{Im}(E_s^*B_r)\over4}
  \operatorname{Im}(I_{21}-I_{30}).
\]

The signed full-tensor difference at zero FLR is

\[
D^{\mathrm{integral}}-D^{\mathrm{legacy}}
=\begin{pmatrix}0&-a\\a&0\end{pmatrix}.
\]

The diagnostic retains this discrepancy. It does not project the legacy
reference onto its symmetric part before reporting the four coefficient
residuals. Full equality requires the phase term to vanish. Balanced complex
drives with number-only and momentum-conserving susceptibility models expose
a finite discrepancy even when the FLR parameter tends to zero; the
energy-conserving and in-phase fixtures test the symmetric case. This is a
limitation of the existing integral model, not a change to its evolution.

## Residuals

`kim_qldiff_m::calc_dqli_limit_benchmark` returns both complete 2-by-2 tensors
and componentwise residuals. For each coefficient,

\[
\Delta_{ij}=|D^{\mathrm{new}}_{ij}-D^{\mathrm{old}}_{ij}|,\qquad
R_{ij}={\Delta_{ij}\over
\max(|D^{\mathrm{old}}_{ij}|,|D^{\mathrm{new}}_{ij}|)}.
\]

When both coefficients are zero, the relative residual is zero. When only
one is zero it is one; opposite signs can give a value up to two. There is
no dimensionful floor of one. Scaling both field amplitudes by a common
factor scales the tensors and absolute errors by its square, while leaving
the relative errors unchanged. The returned tensors preserve the signs
needed to reconstruct the signed difference and its antisymmetric part.

## QL-Balance runtime diagnostic

Enable the diagnostic in `BALANCENML` with `kim_transport_benchmark = .true.`.
It requires `wave_code = 'KIM'`, `kim_run_type = 'electrostatic_periodic'`,
and `ihdf5IO = 1`. Keep `kim_ion_transport_model` set to the model intended
for evolution. Both choices retain their usual selected coefficients; the
extra evaluation is diagnostic. Electrons continue to use the existing
drift-kinetic path.

At the saved-profile cadence the output group is
`/<h5_mode_groupname>/TransportBenchmark/<time_index>/mode_<1-based-index>/`.
It contains the full `drift_kinetic`, `finite_larmor_radius`,
`absolute_residual`, `relative_residual`, and `drift_antisymmetric` arrays
with Fortran shape `(2,2,nrad)`, plus `r`, `m`, `n`, `transition_weight`, and
`embedding_bounds_and_widths`. The last array records core bounds and the
requested/effective transition widths. `drift_antisymmetric` is half the
difference between the drift tensor and its transpose.

The time group records `selected_ion_model`, `kim_ion_conservation_model`,
`drift_energy_conservation`, and `drift_cutoff_halfwidth` so the comparison
can be interpreted. The runtime reference is the actual QL-Balance
`calc_transport_coeffs_ornuhl` result, including its cutoff at
`r_resonant +/- 2*gg_width` and its global energy-conservation boolean. The
integral result retains KIM's per-species conservation model, finite radial
spectrum, and compact embedding. These policies can differ. Runtime
residuals therefore combine model, background, cutoff, embedding, and
discretization differences; they are not automatically an algebraic
zero-FLR acceptance test. Inspect the trusted core and recorded settings
before attributing an observed difference to a particular source.

## Symbolic and numerical coverage

`KIM/mathematica/verify_quasilinear_flr_transport.wl` limits the actual
field-contracted transport polynomials with arbitrary complex fields and
symmetric complex susceptibility moments. It checks all four symmetric
coefficients against separately written legacy formulas and proves the
explicit full-tensor residual above. It also checks that the phase term can
be nonzero. This is stronger than recording the three limiting perpendicular
moments alone.

`KIM/tests/test_dqli_limit_benchmark.f90` compares the physical electric
channel directly with the established drift implementation, exercises
balanced mixed fields and several phases for all supported conservation
models, and checks decreasing FLR parameters. It tests quadratic amplitude
scaling, scale-independent residuals, small coefficients, sign reversal,
and zero-reference cases. Tests of the symmetric limiting identity coexist
with tests that require the full antisymmetric mismatch to remain visible.

`KIM/tests/test_dqli_periodic_limit.f90` uses constant plasma parameters and
a nontrivial periodic field spectrum. It compares the reconstructed
pointwise tensor with the drift reference, independently varying Fourier
cutoff, sampling count, and the physical FLR parameter.

`QL-Balance/src/test/test_periodic_kim_coupling.f90` exercises a real KIM
solve through QL transport selection and HDF5 readback. With the diagnostic
enabled and disabled, both selectable ion models return identical evolution
coefficients and leave physical fields unchanged. The saved residuals are
checked against the recorded tensors, and disabling the diagnostic prevents
stale benchmark output.

Nonzero \(B_\parallel\) field-channel coverage remains in the low-level FLR
algebra/oracle tests. The production electrostatic periodic solver rejects a
nonzero parallel magnetic drive; `test_periodic_ion_tensor` checks that
rejection. These low-level extension tests are not evidence of a working
production periodic compression response.

## Error budget

The error sources have different meanings and must be evaluated separately:

- Algebraic identities use exact symbolic simplification; the numerical
  physical-field reference comparison uses a `2e-13` relative tolerance for
  floating-point evaluation. Existing oracle fixtures retain their separate
  high-precision checks.
- Finite-FLR model error is measured relative to the symmetric legacy
  tensor. The local test decreases `ks` from `1e-2` through `1e-3` to `1e-4`,
  requires at least a factor of 50 reduction per decade, and bounds the last
  error by `2e-7`. This is an asymptotic error check, not a radial quadrature
  or DFT accuracy estimate.
- Periodic reconstruction is tested separately. At fixed physical
  parameters, increasing the Fourier cutoff through 1, 2, 4, and 6 bounds
  truncation error by `1e-5`; sampling counts 10, 20, and 40 at fixed cutoff
  isolate aliasing, with the resolved error below `1e-11`. A separate FLR
  scan halves tangential and radial wave numbers together, requires the
  error ratio to lie between 0.15 and 0.35, and bounds the final error by
  `1e-3`. These gates apply to this smooth synthetic periodic fixture,
  not to arbitrary plasma profiles or to the full QL evolution solver.
- A nonzero antisymmetric legacy term is a model discrepancy. It must not
  be described as discretization error or expected to vanish with grid
  refinement.
