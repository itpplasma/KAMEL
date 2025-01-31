subroutine read_config

    use plasma_parameter
    use config
    use constants
    use setup
    use grid
    use cut_off_integration
    use equilibrium, only: calculate_equil
    use unit_tests, only: test_all, test_sparse_solver
    use poisson_solver, only: solve_poisson

    implicit none

    namelist /KIM_CONFIG/ profile_location, hdf5_input, hdf5_output, &
                        fdebug, fstatus, number_of_ion_species, output_path, artificial_debye_case, &
                        kernel_debye_case, type_of_run

    namelist /KIM_SETUP/ btor, R0, m_mode, n_mode, Zi, Ai, omega, spline_base, &
                        cut_off_fac, kr_cut_off_fac, r_plas, type_br_field, collisions_off, eps_reg

    namelist /KIM_GRID/ reduce_r, grid_spacing, l_space_dim, num_gengrid_points, &
                        reduced_rg_dim, kr_grid_width_res, kr_grid_ampl_res, k_space_dim

    open(unit = 77, file = './KIM_config.nml')
    read(unit = 77, nml = KIM_CONFIG)
    allocate(Zi(number_of_ion_species), Ai(number_of_ion_species))
    read(unit = 77, nml = KIM_SETUP)
    read(unit = 77, nml = KIM_GRID)
    close(unit = 77)

    write(*,*) '+ + + + + + + + KIM + + + + + + + + + + + + + + + +'
    write(*,*) ' type of run = ', type_of_run
    write(*,*) ' - - - - - - - - - - - - - - - - - - - - - - - -'
    write(*,*) 'Configuration namelist'
    write(*,*) '  profile_location = ', trim(profile_location)
    write(*,*) '  output_path      = ', trim(output_path)
    write(*,*) '  hdf5_input       = ', hdf5_input
    write(*,*) '  hdf5_output      = ', hdf5_output
    write(*,*) '  fdebug           = ', fdebug
    write(*,*) '  fstatus          = ', fstatus
    write(*,*) '  number_of_ion_species         = ', number_of_ion_species
    write(*,*) '  artificial_debye_case= ', artificial_debye_case
    write(*,*) '  kernel_debye_case= ', kernel_debye_case
    write(*,*) ' - - - - - - - - - - - - - - - - - - - - - - - -'
    write(*,*) 'Setup namelist'
    write(*,*) '  btor             = ', btor
    write(*,*) '  R0               = ', R0
    write(*,*) '  r_plas           = ', r_plas
    write(*,*) '  m_mode           = ', m_mode
    write(*,*) '  n_mode           = ', n_mode
    write(*,*) '  Zi               = ', Zi
    write(*,*) '  Ai               = ', Ai
    write(*,*) '  omega            = ', omega
    write(*,*) '  spline_base      = ', spline_base
    write(*,*) '  cut_off_fac      = ', cut_off_fac
    write(*,*) '  kr_cut_off_fac   = ', kr_cut_off_fac 
    write(*,*) '  type_br_field    = ', type_br_field
    write(*,*) '  collisions_off   = ', collisions_off
    write(*,*) '  eps_reg          = ', eps_reg
    write(*,*) ' - - - - - - - - - - - - - - - - - - - - - - - -'
    write(*,*) 'Grid namelist'
    write(*,*) '  k_space_dim      = ', k_space_dim
    write(*,*) '  l_space_dim      = ', l_space_dim
    write(*,*) '  reduce_r         = ', reduce_r
    write(*,*) '  reduced_rg_dim    = ', reduced_rg_dim
    write(*,*) '  num_gengrid_points = ', num_gengrid_points
    write(*,*) '  grid_spacing     = ', grid_spacing
    write(*,*) '  kr_grid_width_res = ', kr_grid_width_res
    write(*,*) '  kr_grid_ampl_res = ', kr_grid_ampl_res
    write(*,*) ' - - - - - - - - - - - - - - - - - - - - - - - -'
    write(*,*) ' - - - - - - - - - - - - - - - - - - - - - - - -'

end subroutine