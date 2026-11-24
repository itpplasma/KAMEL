# Configuration of KIM
KIM is configured via the namelist file KIM_config.nml containing multiple namelists. The following gives an overview of the purpose of each variable.

`KIM.x` takes the namelist path as its only argument (default `./KIM_config.nml`).
The groups `KIM_CONFIG`, `WKB_DISPERSION`, `KIM_IO`, `KIM_SETUP`, `KIM_GRID` and
`KIM_PROFILES` are required and are read sequentially in this order, so they
must appear in the file in this order (other groups may be interleaved).
`WKB_DISPERSION` must be present even when `type_of_run` is not
`WKB_dispersion`; an empty `&WKB_DISPERSION /` is sufficient. `KIM_PERIODIC`
and `KIM_FLR2` are optional. `KIM_SPECIES` is read only when
`read_species_from_namelist = .true.`. An option that the code does not read
aborts the run, so remove options that are no longer listed here.

Options listed without a default have no default assigned in the source and
must be set explicitly.

Each run writes into the subdirectory `m<m_mode>_n<n_mode>/` of `output_path`.

The complete option set is checked against the source by
`test/test_kim_namelist_docs.py`.

## KIM_CONFIG
- number_of_ion_species ... integer, number of ion species
- read_species_from_namelist ... boolean, if true, read the ion species from `&KIM_SPECIES`; if false, initialize a single-ion plasma selected by `plasma_type`
- plasma_type ... string, ion species used when `read_species_from_namelist = .false.`: 'D' (deuterium) or 'H' (hydrogen). Any other value aborts the run
- type_of_run ... string, specifies the type of run, options include: 'electrostatic', 'electrostatic_periodic', 'electromagnetic', 'WKB_dispersion', 'flr2', 'flr2_benchmark'. `flr2` runs the imported KiLCA-FLR2 local second-order response using KIM's prepared background; `flr2_benchmark` runs KIM's non-local FLR2 kernel benchmark.
- collision_model ... string, collision model to use, options are: 'Krook', 'FokkerPlanck'
- ion_collision_model ... string, ion model used with `collision_model='FokkerPlanck'`: `FokkerPlanck` (default) or `collisionless`. In `electrostatic_periodic`, collisionless mode keeps FP electrons and uses the $m_\phi=0$ collisionless ion kernels for both charge and parallel current.
- collisionless_kpar_epsilon ... real, positive causal-pole width in 1/cm required by `ion_collision_model='collisionless'`. KIM uses $k_{\parallel}+i\epsilon$ for signed inverse factors and $\sqrt{k_{\parallel}^2+\epsilon^2}$ for even charge-response factors.
- ion_fp_collision_scale ... real, positive multiplier applied only to ion FP collision frequencies. It has no effect on collisionless ion kernels.
- collision_frequency_scale ... real, positive multiplier applied to the calculated electron and ion collision frequencies before the Krook argument z0 and the Fokker-Planck susceptibility inputs are assembled; default 1.0 leaves the collision formulas unchanged
- electron_ifunc_conservation_model ... integer, electron FPGEN conservation model: `0` particle number, `1` number and energy, `2` number and parallel momentum, `3` number, energy, and parallel momentum. Default `-1` inherits `boole_energy_conservation`. Use `1` for resistive calculations; electron momentum conservation removes the resistive response.
- ion_ifunc_conservation_model ... integer, ion FPGEN conservation model with the same numbering. Use `1` for the tokamak-like flow-braking model and `3` for the momentum-conserving cylindrical model.
- boole_energy_conservation ... deprecated boolean fallback used only for a species whose integer model is `-1`; `.true.` resolves to model `1`, `.false.` to model `0`
- artificial_debye_case ... integer, if 0: full calculation of kernel, if 1: Debye case, if 2: exclude Debye case
- turn_off_electrons ... boolean, if true, the electron (species 0) contributions are skipped when the response kernels are assembled. The `flr2` run type also drops the electron potential and current terms and rejects disabling both electrons and ions
- turn_off_ions ... boolean, if true, the ion contributions are skipped when the response kernels are assembled; in the global electrostatic run the per-ion potentials `Phi_m_<species>` are written as zero
- rescale_density ... boolean, if true, multiply the density and its radial derivative of every species (electrons and ions) by `number_density_rescale` and recompute the Debye lengths. The rescaling is applied after the collision frequencies are calculated, so these use the unscaled density
- number_density_rescale ... real, density multiplier used when `rescale_density = .true.`
- ion_flr_scale_factor ... real, multiplier applied to the ion thermal Larmor radius and to the ion FLR terms of the Fokker-Planck kernels; electron terms are never scaled. Use `1.0` for the physical Larmor radius
- ion_temperature_gradient_model ... string, diagnostic control of the ion thermodynamic-force terms (default `'full'`): `'full'` keeps both temperature-gradient terms; `'zero_A2'` sets the ion `A2` to zero; `'zero_Tprime'` sets both the temperature-gradient part of the ion `A1` and the ion `A2` to zero. Only the ion forces change; the Ti profile and every other quantity derived from it are unchanged. Electron forces are never modified

Collisionless forced periodicity currently requires `artificial_debye_case=0`,
`theta_integration='GaussLegendre'`, and a positive
`collisionless_kpar_epsilon`. Higher cyclotron harmonics are not included in
this ion model; the implemented derivation uses $m_\phi=0$.

For the cylindrical comparison requested in the Er scan, use:

```fortran
electron_ifunc_conservation_model = 1  ! N+E: retain electron resistivity
ion_ifunc_conservation_model = 3       ! N+E+P: conserve ion parallel momentum
```

## KIM_IO
- profile_location ... string, directory of the prepared `r_eff` profiles (`n.dat`, `Te.dat`, `Ti.dat`, `q.dat`, `Er.dat`). The profile reader appends the file names without a separator, so the path must end with `/`
- output_path ... string, path to the output directory; results are written to `<output_path>/m<m_mode>_n<n_mode>/`
- hdf5_input ... boolean, must be `.false.`. Reading profiles from HDF5 is not implemented: with `.true.` KIM only logs a message and reads no profiles
- hdf5_output ... boolean, if true writes to a hdf5 output file
- h5_out_file ... string, name of the HDF5 output file inside the run output directory; required when `hdf5_output = .true.`. An existing file of that name is deleted at startup
- log_level ... integer, console log verbosity (default 3): -1 silent, 0 results only, 1 errors, 2 warnings, 3 info, 4 debug, 5 trace; values outside this range are clamped
- data_verbosity ... integer, amount of auxiliary data written (default 1). Values of 3 or more additionally write the kernel prefactors and a sparse-matrix check file of the Poisson system
- write_diagnostics_dat ... boolean, write the scalar resonance diagnostics (`dqle22`, `Ipar`, `Ipar_e`) to `kim_diagnostics.dat` in the run output directory, independent of `hdf5_output` (default false). The file is not written when the resonant surface lies outside the radial grid
- calculate_asymptotics ... boolean, if true calculates asymptotic quantities (FLR2 asymptotic potential and Fourier space kernel)

## KIM_SETUP
- btor ... double, toroidal magnetic field at the magnetic axis in Gauss
- r0 ... double, major radius of the magnetic axis in cm
- m_mode ... integer, poloidal mode number of the resonant surface
- n_mode ... integer, toroidal mode number of the resonant surface
- omega ... double, perturbation frequency (rad/s)
- spline_base ... integer, specifies which spline base to use in the FEM procedure, options: 1 - hat functions
- type_br_field ... integer, specifies which type of Br field (array) is used. The `flr2` run type supports 11 (read `./inp/Br_in.dat`) and 12 (constant complex field from `Br_boundary_re/im`).
- collisions_off ... boolean, if true, sets collision frequencies to zero
- set_profiles_constant ... integer, if 1, sets profiles constant to the core value
- bc_type ... integer, boundary condition type. 0: None; 1: zero Neuman left, zero Dirichlet right; 2: Dirichlet left and right; 3: zero-misalignment
- mphi_max ... integer, maximum number of cyclotron harmonics to include in kernel calculations
- Br_boundary_re ... double, real part of Br at right boundary (default 1.0), used by the `electromagnetic` run type and as the constant Br amplitude for `flr2` and `electrostatic_periodic`
- Br_boundary_im ... double, imaginary part of Br at right boundary (default 0.0), used by the `electromagnetic` run type and as the constant Br amplitude for `flr2` and `electrostatic_periodic`

## WKB_DISPERSION

Used by `type_of_run = 'WKB_dispersion'`, which solves the local dispersion
relation D = 0 for the complex radial wavenumber at each point of the r_g grid.
The group must be present for every run type. Results are written to
`<run output directory>/dispersion/`.

- WKB_dispersion_mode ... string, dispersion function (default `'KIM'`): `'KIM'` uses the full Bessel-function form, `'FLRE'` the finite-Larmor-radius expansion
- WKB_dispersion_solver ... string, root finder (default `'Muller'`): `'Muller'` tracks one root from grid point to grid point, using the previous root as initial guess (recommended); `'region_roots'` uses contour integration with Newton refinement (fortnum `complex_region_roots`) and tracks several branches
- WKB_solve_for_kr_squared ... boolean, if true the unknown is kr^2 and D(kr^2) = 0 is solved; if false (default) the unknown is kr
- WKB_use_initial_root ... boolean, `Muller` only: seed the first grid point with `WKB_initial_root_real + i WKB_initial_root_imag` (default false)
- WKB_initial_root_real ... real, real part of the initial root, in the same variable as the unknown (kr or kr^2; default 0)
- WKB_initial_root_imag ... real, imaginary part of the initial root (default 0)
- WKB_initial_root_profile ... string, `Muller` only: optional file with one line `r, Re(root), Im(root)` per r_g grid point. Its radii must match the r_g grid to 1e-8; when set, it seeds every grid point and takes precedence over the previous-root guess
- WKB_max_tracked_branches ... integer, `region_roots` only: maximum number of simultaneously tracked branches (default 4)
- WKB_branch_search_halfwidth ... real, `region_roots` only: half-width of the square search region in the complex plane of the unknown around each tracked branch (default 1.5)
- WKB_broad_search_halfwidth ... real, `region_roots` only: half-width of the square region used for branch discovery (default 5.0)
- WKB_broad_search_interval ... integer, `region_roots` only: repeat the broad search every N grid points; 0 (default) searches broadly only at the first point and whenever all branches are lost
- WKB_root_tolerance ... real, `region_roots` only: roots with \|D\| at or above this value are discarded (default 1e-6)
- WKB_verbose ... boolean, `region_roots` only: print per-grid-point search details (default false)

## KIM_GRID
- r_plas ... double, minor radius of the plasma in cm
- r_min ... double, minimum radius for the computational domain in cm
- width_res ... double, width of the resonant layer
- ampl_res ... double, amplitude factor for the resonant layer
- hrmax_scaling ... double, scaling factor for maximum grid spacing
- l_space_dim ... integer, dimension of the spline space
- rg_space_dim ... integer, number of r_g grid points (cell boundaries)
   With "equidistant" spacing the grid has exactly this many points, from r_min to r_plas
   inclusive, so the spacing is (r_plas - r_min)/(rg_space_dim - 1); the same holds for
   l_space_dim.
 - grid_spacing_rg ... string, r_g grid spacing: "equidistant", "non-equidistant", or "adaptive"
 - grid_spacing_xl ... string, x_l grid spacing: "equidistant", "non-equidistant", or "adaptive"
- theta_integration ... string, angular integration selector:
  - "GaussLegendre": fixed-order Gauss–Legendre quadrature
  - "RKF45" or "QUADPACK": adaptive path (the actual adaptive integrator is selected by `theta_integration_method`)
- theta_integration_method ... string, adaptive integrator; normally omitted, because KIM sets it from `theta_integration` ("RKF45" or "QUADPACK") and overwrites any value given here in that case:
  - "RKF45": embedded Runge–Kutta–Fehlberg (adaptive)
  - "QUADPACK": QUADPACK-style adaptive Gauss–Kronrod integration (QAG/QAGS) provided by fortnum
- quadpack_algorithm ... string, QUADPACK algorithm for finite-interval θ integrals: "QAG" or "QAGS"
  - 'QAG': General-purpose adaptive integration (recommended for most cases)
  - 'QAGS': Adaptive integration with singularities at endpoints
- quadpack_key ... integer, Gauss–Kronrod rule key for QAG (1–6 for 15–61 points; 6 = 61‑point); QAGS uses its own fixed rule
    - 1: QK15 (7-15 points) - Fast but less accurate
    - 2: QK21 (10-21 points)
    - 3: QK31 (15-31 points)
    - 4: QK41 (20-41 points)
    - 5: QK51 (25-51 points)
    - 6: QK61 (30-61 points) - Slow but most accurate (recommended)
- quadpack_limit ... integer, maximum number of subintervals for QAG/QAGS
    - Increase if integration fails to converge
    - Higher values allow better handling of difficult integrands
- quadpack_epsabs ... double, absolute tolerance for QAG/QAGS
- quadpack_epsrel ... double, relative tolerance for QAG/QAGS
- quadpack_use_u_substitution ... boolean, use u = sin(θ/2) mapping to improve stability near endpoints
- rkf45_atol ... double, absolute tolerance for RKF45 adaptive θ-integration
- rkf45_rtol ... double, relative tolerance for RKF45 adaptive θ-integration
- kernel_taper_skip_threshold ... double, threshold for skipping a kernel contribution when the distance-based taper weight falls below this value; must lie in (0, 1)
- Larmor_skip_factor ... double, scaling of the distance-based taper exp(-(d/(alpha*rhoT))^2); larger values reduce skipping by broadening support
- gauss_int_nodes_Nx ... integer, number of Gauss integration nodes in x direction (should differ from Nxp)
- gauss_int_nodes_Nxp ... integer, number of Gauss integration nodes in x' direction
- gauss_int_nodes_Ntheta ... integer, number of Gauss integration nodes in theta direction (used only when theta_integration = "GaussLegendre")

Notes:
- The QAG/QAGS integrators come from fortnum, which CMake fetches with the other pinned dependencies; no separate QUADPACK or SLATEC sources are needed.
- QUADPACK supports only finite-interval θ integrals here; "QAGI" (infinite intervals) is not used by the kernels.

## KIM_SPECIES

Read only when `read_species_from_namelist = .true.`; the first
`number_of_ion_species` entries are used (at most 16).

- zi ... integer array, ion charge numbers (e.g., zi = 1, 2 for H+ and He++)
- ai ... integer array, ion mass numbers (e.g., ai = 2, 4 for D and He)

## KIM_PROFILES

Controls how the profiles in `profile_location` are prepared before the solver
reads them. With `coord_type = 'sqrt_psiN'`, KIM transforms raw profiles from
`input_profile_dir` to `r_eff` and writes them to `profile_location`. If
`Er_file` is missing from `profile_location`, Er is computed from radial force
balance without poloidal rotation. KIM also checks `r_min` and `r_plas` against
the radial range of the density profile, and compares `btor` and `R0` with
`btor_rbig.dat` in `profile_location` when that file exists.

The solver itself reads `n.dat`, `Te.dat`, `Ti.dat`, `q.dat` and `Er.dat` from
`profile_location`; the `*_file` options below rename only the files used by
the preparation step.

- coord_type ... string, radial coordinate of the input profiles (default `'auto'`): `'r_eff'` effective radius in cm, used as is; `'sqrt_psiN'` normalized poloidal-flux radius, transformed to `r_eff` with the equilibrium; `'auto'` selects `'r_eff'` if the largest radius in the density input exceeds 2, otherwise `'sqrt_psiN'`
- input_profile_dir ... string, directory of the raw input profiles (default `'./'`). For `'r_eff'` input in a different directory than `profile_location`, only `Er_file` is copied from here
- equil_file ... string, equilibrium table `equil_r_q_psi.dat` used for the `sqrt_psiN` transformation (default empty: `<input_profile_dir>/equil_r_q_psi.dat`)
- geqdsk_file ... string, GEQDSK file used to compute the equilibrium when the equilibrium table does not exist (default empty). This path requires the environment variable `$CODE` to locate `$CODE/KAMEL/common/equil/convexwall/convexwall.asdex`, an ASDEX Upgrade wall contour, and writes `field_divB0.inp` into the working directory
- n_input_file ... string, density input file in `sqrt_psiN` coordinates (default `'n_of_psiN.dat'`); also the first file tried by `coord_type = 'auto'`
- Te_input_file ... string, electron-temperature input file in `sqrt_psiN` coordinates (default `'Te_of_psiN.dat'`)
- Ti_input_file ... string, ion-temperature input file in `sqrt_psiN` coordinates (default `'Ti_of_psiN.dat'`)
- Vz_input_file ... string, toroidal-rotation input file in `sqrt_psiN` coordinates (default `'Vz_of_psiN.dat'`)
- n_file ... string, prepared density file name (default `'n.dat'`); used by the `auto` detection fallback (in `input_profile_dir`), the radial-range check and the Er force-balance calculation
- Te_file ... string, prepared electron-temperature file name (default `'Te.dat'`); accepted but currently not used
- Ti_file ... string, prepared ion-temperature file name (default `'Ti.dat'`); used by the Er force-balance calculation
- Vz_file ... string, prepared toroidal-rotation file name (default `'Vz.dat'`); used by the Er force-balance calculation, which assumes Vz = 0 if the file is missing
- Er_file ... string, prepared radial-electric-field file name (default `'Er.dat'`); if it is missing from `profile_location`, Er is computed from force balance
- q_file ... string, prepared safety-factor file name (default `'q.dat'`); used by the Er force-balance calculation
- geometry_file ... string, optional flux-surface geometry from a host code (default `''`, KIM's screw pinch); a relative path is taken relative to `profile_location`. Columns `r B0 kp_m kp_n ks_m ks_n radial_metric` give the field magnitude, $k_\parallel = m\,k_{p,m} + n\,k_{p,n}$, $k_s = m\,k_{s,m} + n\,k_{s,n}$ and the radial metric $M$ (physical radial wavenumber $\sqrt{M}\,k_r$). The table replaces B0, kp and ks (and hence om_E) of every run type; a metric other than 1 is supported only by `electrostatic_periodic`, where it enters the Larmor-radius arguments and the radial Poisson operator $\partial_r (M \partial_r)$

## KIM_PERIODIC

- periodic_dr_asis_scale ... real, half-width of the unchanged profile region in units of the selected reference Larmor radius
- periodic_dr_tr_scale ... real, width of the transition buffer on each side in units of the selected reference Larmor radius
- periodic_kmax_scale ... real, radial Fourier cutoff multiplied by the reference Larmor radius; the solver derives the integer harmonic count $M$ from this cutoff and the period
- periodic_n_rg ... integer, number of endpoint-exclusive guiding-centre samples over one period
- periodic_calculate_radial_current ... boolean, assemble and reconstruct the periodic radial-current response; default false. Set true when `/fields/jrad` is needed
- periodic_calculate_ion_tensor ... boolean, calculate and return the periodic ion diffusion tensor; default false. Set true for QL-Balance transport or other runs that need `D_ion`; field/current/diagnostic-only runs leave it unavailable
- periodic_Bparallel_ratio ... complex, ratio of the constant parallel magnetic-field drive to the constant radial drive: the periodic solve uses `Bparallel = periodic_Bparallel_ratio * (Br_boundary_re + i Br_boundary_im)`. Default `(0.0, 0.0)` disables the parallel drive (purely radial benchmark); set it explicitly for coupled runs
- periodic_match_global_kernel_approximations ... boolean, if true drops $k_s^2$ from the FLR arguments and takes electrons in the zero-FLR limit for comparison with the simplified global model; default false

For `ion_collision_model='collisionless'`, these settings control the same
periodic window and Fourier basis as in FP mode. The assembled response is FP
electrons plus collisionless ion charge and current. Both `Phi` and `jpar` are
written through the standard forced-periodicity output path.

## KIM_FLR2

This optional group controls only terms in the KiLCA-FLR2 response. KIM's
existing profile, equilibrium, collision-frequency, density-rescaling, and
energy-conservation settings remain authoritative.

- flr2_electron_flr ... boolean, include the electron second-order FLR terms
- flr2_ion_flr ... boolean, include the ion second-order FLR terms
- flr2_electron_potential ... boolean, include electrons in the potential equation
- flr2_ion_potential ... boolean, include ions in the potential equation
- flr2_electron_current ... boolean, include electrons in the parallel current
- flr2_ion_current ... boolean, include ions in the parallel current
- flr2_include_potential_in_current ... boolean, include the solved potential contribution in the current; if false, `Phi` is set to zero and only the magnetic-drive current is returned

All seven switches default to `.true.` when the group is omitted. The `flr2`
run type requires `collision_model='FokkerPlanck'`, exactly one ion species,
`omega=0`, nonzero radial electric field over the solve domain, and
`type_br_field` 11 or 12. Existing `turn_off_electrons` and `turn_off_ions`
settings also gate the corresponding FLR2 potential and current contributions.
It writes `Br`, `Phi`, and `jpar` on KIM's guiding-centre boundary grid.
