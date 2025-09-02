# Configuration of KIM
KIM is configured via the namelist file KIM_config.nml containing multiple namelists. The following gives an overview of the purpose of each variable. 

## KIM_CONFIG
- number_of_ion_species ... integer, number of ion species
- read_species_from_namelist ... boolean, if false, initialize deuterium plasma
- type_of_run ... string, specifies the type of run, options are: 'electrostatic'
- collision_model ... string, collision model to use, options are: 'Krook', 'FokkerPlanck'
- artificial_debye_case ... boolean, if true, only considers the term in the rho phi kernel that describes Debye shielding

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
- type_br_field ... integer, specifies which type of Br field (array) is used
- collisions_off ... boolean, if true, sets collision frequencies to zero
- eps_reg ... double, regularization variable for numerical stability
- set_profiles_constant ... integer, if 1, sets profiles constant to the core value

## KIM_GRID
- r_plas ... double, minor radius of the plasma in cm
- r_min ... double, minimum radius for the computational domain in cm
- width_res ... double, width of the resonant layer
- ampl_res ... double, amplitude factor for the resonant layer
- hrmax_scaling ... double, scaling factor for maximum grid spacing
- k_space_dim ... integer, dimension of the k space
- l_space_dim ... integer, dimension of the spline space
- reduce_r ... boolean, if true reduces r dim of the input profiles
- reduced_rg_dim ... integer, dimension of reduced r dim
- grid_spacing ... integer, choose type of grid spacing (1=equidistant, 2=non-equidistant, 3=adaptive)
- num_gengrid_points ... integer, minimal number of grid points in the l grid
- kr_grid_width_res ... double, width parameter for k-space grid near resonance
- kr_grid_ampl_res ... double, amplitude parameter for k-space grid near resonance
- theta_integration ... string, angular integration method: "RKF45" (adaptive) or "GaussLegendre" (fixed)
- rkf45_atol ... double, absolute tolerance for RKF45 adaptive θ-integration
- rkf45_rtol ... double, relative tolerance for RKF45 adaptive θ-integration
- kernel_taper_skip_threshold ... double, threshold for skipping a kernel contribution when the distance-based taper weight falls below this value
- Larmor_skip_factor ... double, scaling of the distance-based taper exp(-(d/(alpha*rhoT))^2); larger values reduce skipping by broadening support
- gauss_int_nodes_Nx ... integer, number of Gauss integration nodes in x direction (should differ from Nxp)
- gauss_int_nodes_Nxp ... integer, number of Gauss integration nodes in x' direction
- gauss_int_nodes_Ntheta ... integer, number of Gauss integration nodes in theta direction (used only when theta_integration = "GaussLegendre")

## KIM_SPECIES
- zi ... integer array, ion charge numbers (e.g., zi = 1, 2 for H+ and He++)
- ai ... integer array, ion mass numbers (e.g., ai = 2, 4 for D and He)
