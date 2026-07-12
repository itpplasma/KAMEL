# Electrostatic no-op benchmark for the periodic-Fourier solver

Acceptance slice for issue #175: with the periodization layer wide enough to
cover the localized response region, the enforced-periodicity Fourier solve
(`type_of_run = 'fourier_periodic'`) must reproduce the hat-basis
electrostatic solve (`type_of_run = 'electrostatic'`) inside the flat part of
the layer, where the periodized background is identical to the original one.
Both solves use the same constant unit Br drive (`type_br_field = 12`), the
same Krook background, and the same profiles.

Status: the benchmark currently FAILS, and neither solve is converged on
grids that fit the intended test runtime. The measured numbers are recorded
below; this case is committed as the reproducible starting point, not as a
registered ctest. Do not register a pass/fail test from these settings.

## Case definition

- Profiles `n.dat`, `Te.dat`, `Ti.dat`, `q.dat`: copied unchanged from
  `test/golden/kim/cases/electromagnetic/` (r_eff range [10, 40] cm).
  `Er.dat` is regenerated deterministically by KIM from force balance on the
  first run and is therefore not committed.
- `KIM_config.nml`: hat-basis reference, derived from `KIM/nmls/KIM_config.nml`
  with `type_of_run = 'electrostatic'`, `bc_type = 3`, `type_br_field = 12`,
  `m_mode = -3`, `n_mode = 2` (resonance q = -m/n = 1.5 at r = 20 cm, kp zero
  crossing at r = 19.96 cm on the run grid), `omega = 0`, domain
  r in [10, 35.5] cm, `l_space_dim = rg_space_dim = 128`.
- `KIM_config_periodic.nml`: identical except `type_of_run =
  'fourier_periodic'` plus `fp_r_res = 19.995`, `fp_dr_layer = 2.3`,
  `fp_dr_transition = 0.9`, `fp_n_modes = 128`, `fp_n_quad = 512`,
  `fp_grid_points = 2048`.
- Collision model: `Krook` is forced, not a speed choice. The periodic solver
  evaluates the restored continuum kernels of `KIM/src/kernels/kernel_mod.f90`,
  which are built on `plasma_Z(z0)` with the Krook z0; there is no
  periodic-Fourier counterpart of the FokkerPlanck hat kernels. Note that
  `KIM/nmls/KIM_config.nml` labels Krook "(unstable)", consistent with the
  non-convergence observed below.
- Layer sizing: the periodization needs profile coverage of
  fp_r_res +- 3 (fp_dr_layer + fp_dr_transition). With profiles on [10, 40]
  and r_res near 20, the half-width is capped at about 3.33 cm; 3.2 cm
  (2.3 + 0.9) is close to that maximum, so the flat layer [17.7, 22.3] cm is
  the largest no-op window this profile set allows. The hat solution is
  peaked inside that window but retains about 5 percent of its peak
  amplitude at r = 12 cm, i.e. it is only approximately localized.

## Procedure (manual)

Run each solve in its own scratch directory containing the four profile
`.dat` files and the corresponding namelist as `KIM_config.nml` (or pass the
namelist path as the single command-line argument of `KIM.x`):

    KIM.x                 # hat:      out/m-3_n2/fields/Phi.dat
    KIM.x                 # periodic: out/m-3_n2/fields/Phi_periodic.dat

then compare on the flat layer:

    python3 compare_phi.py <hat>/out/m-3_n2/fields/Phi.dat \
        <periodic>/out/m-3_n2/fields/Phi_periodic.dat

`compare_phi.py` interpolates the hat solution onto the periodic output grid
restricted to |r - fp_r_res| <= fp_dr_layer, removes the k = 0 gauge
constant (the hat solve pins boundary values via bc_type = 3, the periodic
solve fixes the mean over one period), and reports relative L2 and Linf of
the difference normalized to the hat solution. Acceptance threshold:
offset-removed relative L2 <= 0.05.

## Measured results (2026-07-12, Release build)

Runtimes on a 32-core desktop: hat 128 grid 39 s, hat 256 grid 241 s,
periodic 128 modes 90 s, periodic 256 modes 829 s.

Comparison of the committed settings (hat 128 grid vs periodic 128 modes) on
the flat layer, recorded in `hat_window.txt` / `periodic_window.txt`:

- max |Phi| hat: 1.741e10 statV; max |Phi| periodic: 2.132e-1 statV
  (amplitude ratio about 8e10).
- raw relative L2 = 1.000; offset-removed relative L2 = 0.975,
  relative Linf = 0.880. FAIL against the 0.05 threshold.
- Allowing an optimal complex rescaling of the periodic solution still
  leaves a residual of 0.976, so the two solutions are uncorrelated in
  shape; the disagreement is not a normalization constant.

Convergence checks (required < 1e-6 relative L2 change under one doubling;
both fail):

- periodic, fp_n_modes 128 -> 256 with fp_n_quad 512 -> 1024:
  relative L2 change 1.072 (max |Phi| 0.242 -> 0.346 statV).
- hat, l_space_dim = rg_space_dim 128 -> 256: relative L2 change 1.000 in
  the window, peak amplitude grows by a factor 3.6e3
  (1.74e10 -> 6.24e13 statV). The Krook hat reference diverges under grid
  refinement at omega = 0, so it is not a usable reference at these
  settings.

## Consequences

- No ctest is registered: with both solves unconverged and disagreeing by
  eleven orders of magnitude, any threshold either always fails or encodes
  the bug.
- The disagreement is reported as a finding, not fixed here (no solver
  edits in this change). Follow-up needs (a) a hat-basis reference that is
  stable under grid refinement for a Krook background (or a FokkerPlanck
  variant of the continuum kernels), and (b) a mode/quadrature convergence
  study of the periodic solve once such a reference exists.
