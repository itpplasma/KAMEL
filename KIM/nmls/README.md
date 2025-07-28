# Configuration of KIM
KIM is configured via the namelist file KIM_config.nml containing multiple namelists. The following gives an overview of the purpose of each variable. 

## KIM_CONFIG
- profile_location ... string, path to the profiles
- output_path ... string, path to the output directory
- hdf5_input ... boolean, if true reads from a hdf5 file
- hdf5_output ... boolean, if true writes to a hdf5 output file
- fdebug ... integer, flag for debugging
- fstatus ... integer, flag for giving status updates during the execution of the code
- number_of_ion_species ... integer, number of ion species
- type_of_run ... string, specifies the type of run, options are: electrostatic
- artificial_debye_case ... boolean, if true, only considers the term in the rho phi kernel that describes Debye shielding

## KIM_SETUP
- btor ... double, toroidal magnetic field at the magnetic axis in Gauss
- R0 ... double, major radius of the magnetic axis in cm
- m_mode ... integer, poloidal mode number of the resonant surface
- n_mode ... integer, toroidal mode number of the resonant surface
- Zi ... integer, ion charge number; if number_of_ion_species > 1 write multiple separated by spaces, i.e. Zi = 1 2 4
- Ai ... integer, ion mass number; if ispeceis > 1 write multiple separated by spaces, i.e. Ai = 2 4 8
- r_plas ... double, minor radius of the plasma
- omega ... double, perturbation frequency
- spline_base ... integer, specifies which spline base to use in the FEM procedure, options: 1 - hat functions
- cut_off_fac ... double, factor to configure cut-off of off-diagonal elements of kernels in spline space, i.e. roughly
- type_br_field ... integer, specifies which type of Br field (array) is used
- collisions_off ... boolean, if true, sets collision frequencies to zero
- eps_reg ... double, regularization variable
- set_profiles_constant ... integer, if 1, sets profiles constant to the core value

## KIM_GRID
- k_space_dim ... integer, dimension of the k space
- l_space_dim ... integer, dimension of the spline space
- reduce_r ... boolean, if true reduces r dim of the input profiles
- reduced_rg_dim ... integer, dimension of reduced r dim
- grid_spacing ... integer, choose type of grid spacing (1=equidistant, 2=non-equidistant)
- num_gengrid_points ... integer, minimal number of grid points in the l grid
- kr_grid_width_res ... double
- kr_grid_ampl_res ... double
- delta_l_max ... integer, kernel matrix elements further apart than delta_l_max are not calculated
