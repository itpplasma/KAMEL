program kim_main

    use plas_parameter
    use config
    use constants
    use setup
    use use_libcerf
    use grid

    implicit none

    double complex :: z= (1d0, 1d0)
    double complex :: plasma_Z

    namelist /KIM_CONFIG/ profile_location, hdf5_input, hdf5_output, &
                          fdebug, fstatus, ispecies, output_path

    namelist /KIM_SETUP/ btor, R0, m_mode, n_mode, Zi, Ai, k_space_dim, &
                        reduce_r, reduced_r_dim, omega, spline_base, &
                        grid_spacing, l_space_dim

    open(unit = 77, file = './nmls/KIM_config.nml')
    read(unit = 77, nml = KIM_CONFIG)
    allocate(Zi(ispecies), Ai(ispecies))
    read(unit = 77, nml = KIM_SETUP)
    close(unit = 77)

    write(*,*) ''
    write(*,*) ' - - - - - - - - - - - - - - - - - -'
    write(*,*) 'Configuration namelist'
    write(*,*) '  profile_location = ', trim(profile_location)
    write(*,*) '  output_path      = ', trim(output_path)
    write(*,*) '  hdf5_input       = ', hdf5_input
    write(*,*) '  hdf5_output      = ', hdf5_output
    write(*,*) '  fdebug           = ', fdebug
    write(*,*) '  fstatus          = ', fstatus
    write(*,*) '  ispecies         = ', ispecies
    write(*,*) ' - - - - - - - - - - - - - - - - - -'

    write(*,*) ''
    write(*,*) ' - - - - - - - - - - - - - - - - - -'
    write(*,*) 'Setup namelist'
    write(*,*) '  btor             = ', btor
    write(*,*) '  R0               = ', R0
    write(*,*) '  m_mode           = ', m_mode
    write(*,*) '  n_mode           = ', n_mode
    write(*,*) '  Zi               = ', Zi
    write(*,*) '  Ai               = ', Ai
    write(*,*) '  k_space_dim      = ', k_space_dim
    write(*,*) '  l_space_dim      = ', l_space_dim
    write(*,*) '  reduce_r         = ', reduce_r
    write(*,*) '  reduced_r_dim    = ', reduced_r_dim
    write(*,*) '  omega            = ', omega
    write(*,*) '  spline_base      = ', spline_base
    write(*,*) '  grid_spacing     = ', grid_spacing
    write(*,*) ' - - - - - - - - - - - - - - - - - -'

    ! for the moment:
    l_space_dim = reduced_r_dim

    call generate_k_space_grid(.false.)
    call read_profiles(reduce_r)

    ! calculate equilibrium B field and J
    call calculate_equil(.true.)

    ! calculate quantities used for the kernels, e.g. A1, A2, dndr, omega_c,...
    call calculate_backs(.true.)
    
    !call generate_l_space_grid
    ! calculate kernels
    call kernel_rho_phi(.true.)
    call kernel_rho_B(.true.)

    !call calculate_fourier_trans_spline_funcs(.true.)

    !call basis_transform_kernel(.true.)

end program