program kim_main

    use plas_parameter
    use config
    use constants
    use setup
    use use_libcerf
    use grid

    implicit none

    real :: t_start, t_finish

    !double complex :: z= (1d0, 1d0)
    !double complex :: plasma_Z

    namelist /KIM_CONFIG/ profile_location, hdf5_input, hdf5_output, &
                          fdebug, fstatus, ispecies, output_path

    namelist /KIM_SETUP/ btor, R0, m_mode, n_mode, Zi, Ai, k_space_dim, &
                        reduce_r, reduced_r_dim, omega, spline_base, &
                        grid_spacing, l_space_dim, cut_off_fac, num_gengrid_points

    open(unit = 77, file = './nmls/KIM_config.nml')
    read(unit = 77, nml = KIM_CONFIG)
    allocate(Zi(ispecies), Ai(ispecies))
    read(unit = 77, nml = KIM_SETUP)
    close(unit = 77)

    write(*,*) '+ + + + + + + + KIM + + + + + + + + + +'
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
    write(*,*) '  cut_off_fac      = ', cut_off_fac
    write(*,*) '  num_gengrid_points = ', num_gengrid_points
    write(*,*) ' - - - - - - - - - - - - - - - - - -'

    ! for the moment:
    !l_space_dim = reduced_r_dim

    call cpu_time(t_start)

    call generate_k_space_grid(.true.)
    call read_profiles(reduce_r)

    ! calculate equilibrium B field and J
    call calculate_equil(.true.)

    ! calculate quantities used for the kernels, e.g. A1, A2, dndr, omega_c,...
    call calculate_backs(.true.)
    
    ! generate non-equidistant grid for spline functions (i.e. real space)
    call gengrid(num_gengrid_points, .true.)

    ! calculate kernels
    call kernel_phi(.true.)
    call kernel_B(.true.)

    ! transform the kernels from k space to spline space
    call basis_transform_kernel(.true.)

    ! solve poisson's equation with spline solver
    call solve_poisson

    call cpu_time(t_finish)

    write(*,*) ' Time: ', (t_finish - t_start), ' s'

end program