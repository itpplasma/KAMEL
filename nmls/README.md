# Configuration of KIM
KIM is configured via namelist the file KIM_config.nml containing multiple namelists. The following gives an overview of the purpose of each variable. Note that the variables in the namelist file migth be in a different lower/upper case.

## KIM_CONFIG
- profile_location ... string, path to the profiles
- output_path ... string, path to the output directory
- hdf5_input ... boolean, if true reads from a hdf5 file
- hdf5_output ... boolean, if true writes to a hdf5 output file
- fdebug ... integer, flag for debugging
- fstatus ... integer, flag for giving status updates during the execution of the code
- ispecies ... integer, number of ion species

## KIM_SETUP
- btor ... double precision, toroidal magnetic field at the magnetic axis in Gauss
- R0 ... double precision, major radius of the magnetic axis in cm
- m_mode ... integer, poloidal mode number of the resonant surface
- n_mode ... integer, toroidal mode number of the resonant surface
- Zi ... integer, ion charge number; if ispecies > 1 write multiple separated by spaces, i.e. Zi = 1 2 4
- Ai ... integer, ion mass number; if ispeceis > 1 write multiple separated by spaces, i.e. Ai = 2 4 8
- r_plas ... double precision, minor radius of the plasma
- k_space_dim ... integer, dimension of the k space
- l_space_dim ... integer, dimension of the spline space
- reduce_r ... boolean, if true reduces r dim of the input profiles
- reduced_r_dim ... integer, dimension of reduced r dim
- omega ... double precision, perturbation frequency
- spline_base ... integer, choose the spline basis functions (1=hat functions)
- grid_spacing ... integer, choose type of grid spacing (1=equidistant, 2=non-equidistant)
- cut_off_fac ... double precision
- num_gengrid_points ... integer, minimal number of grid points in the l grid