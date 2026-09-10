# Configuration of KIM
KIM is configured via the namelist file KIM_config.nml containing multiple namelists. The following gives an overview of the purpose of each variable.

## KIM_CONFIG
- number_of_ion_species ... integer, number of ion species
- read_species_from_namelist ... boolean, if false, initialize deuterium plasma
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
- profile_location ... string, path to the profiles
- output_path ... string, path to the output directory
- hdf5_input ... boolean, if true reads from a hdf5 file
- hdf5_output ... boolean, if true writes to a hdf5 output file
- fdebug ... integer, flag for debugging (0=off, 1=basic, 2=detailed)
- fstatus ... integer, flag for giving status updates during the execution of the code
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
- Br_boundary_re ... double, real part of Br at right boundary (default 1.0), used by the `electromagnetic` run type and as the constant Br amplitude for `flr2`
- Br_boundary_im ... double, imaginary part of Br at right boundary (default 0.0), used by the `electromagnetic` run type and as the constant Br amplitude for `flr2`

## KIM_GRID
- r_plas ... double, minor radius of the plasma in cm
- r_min ... double, minimum radius for the computational domain in cm
- width_res ... double, width of the resonant layer
- ampl_res ... double, amplitude factor for the resonant layer
- hrmax_scaling ... double, scaling factor for maximum grid spacing
- k_space_dim ... integer, dimension of the k space
- l_space_dim ... integer, dimension of the spline space
- rg_space_dim ... integer, number of r_g grid points (cell boundaries)
 - grid_spacing_rg ... string, r_g grid spacing: "equidistant", "non-equidistant", or "adaptive"
 - grid_spacing_xl ... string, x_l grid spacing: "equidistant", "non-equidistant", or "adaptive"
- num_gengrid_points ... integer, minimal number of grid points in the l grid
- kr_grid_width_res ... double, width parameter for k-space grid near resonance
- kr_grid_ampl_res ... double, amplitude parameter for k-space grid near resonance
- theta_integration ... string, angular integration selector:
  - "GaussLegendre": fixed-order Gauss–Legendre quadrature
  - "RKF45" or "QUADPACK": adaptive path (the actual adaptive integrator is selected by `theta_integration_method`)
- theta_integration_method ... string, adaptive integrator selection when `theta_integration` is "RKF45" or "QUADPACK"; options:
  - "RKF45": embedded Runge–Kutta–Fehlberg (adaptive)
  - "QUADPACK": Netlib QUADPACK (adaptive, global)
- quadpack_algorithm ... string, QUADPACK algorithm for finite-interval θ integrals: "QAG" or "QAGS"
  - 'QAG': General-purpose adaptive integration (recommended for most cases)
  - 'QAGS': Adaptive integration with singularities at endpoints
- quadpack_key ... integer, Gauss–Kronrod rule key for QAG/QAGS (1–6 for 15–61 points; 6 = 61‑point)
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
- kernel_taper_skip_threshold ... double, threshold for skipping a kernel contribution when the distance-based taper weight falls below this value
- Larmor_skip_factor ... double, scaling of the distance-based taper exp(-(d/(alpha*rhoT))^2); larger values reduce skipping by broadening support
- gauss_int_nodes_Nx ... integer, number of Gauss integration nodes in x direction (should differ from Nxp)
- gauss_int_nodes_Nxp ... integer, number of Gauss integration nodes in x' direction
- gauss_int_nodes_Ntheta ... integer, number of Gauss integration nodes in theta direction (used only when theta_integration = "GaussLegendre")

Notes:
- QUADPACK integration fetches Netlib sources at CMake configure time; network access is required unless cached.
- QUADPACK supports only finite-interval θ integrals here; "QAGI" (infinite intervals) is not used by the kernels.
- A fallback `xerror` (SLATEC error handler) can be enabled via CMake option `KIM_XERROR_FALLBACK` (ON by default) when SLATEC does not provide it.

## KIM_SPECIES
- zi ... integer array, ion charge numbers (e.g., zi = 1, 2 for H+ and He++)
- ai ... integer array, ion mass numbers (e.g., ai = 2, 4 for D and He)

## KIM_PERIODIC

- periodic_dr_asis_scale ... real, half-width of the unchanged profile region in units of the selected reference Larmor radius
- periodic_dr_tr_scale ... real, width of the transition buffer on each side in units of the selected reference Larmor radius
- periodic_kmax_scale ... real, radial Fourier cutoff multiplied by the reference Larmor radius; the solver derives the integer harmonic count $M$ from this cutoff and the period
- periodic_n_rg ... integer, number of endpoint-exclusive guiding-centre samples over one period
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
