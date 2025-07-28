subroutine read_config

    use config
    use constants
    use setup
    use grid
    use poisson_solver, only: solve_poisson

    implicit none

    character(len=256), dimension(:), allocatable :: args
    integer :: ix, num_args
    logical :: ex

    namelist /KIM_CONFIG/ profile_location, hdf5_input, hdf5_output, &
                        fdebug, fstatus, number_of_ion_species, output_path, artificial_debye_case, &
                        kernel_debye_case, type_of_run, collision_model, read_species_from_namelist 

    namelist /KIM_SETUP/ btor, R0, m_mode, n_mode, omega, spline_base, &
                        type_br_field, collisions_off, eps_reg, &
                        set_profiles_constant

    namelist /KIM_GRID/ reduce_r, grid_spacing, l_space_dim, num_gengrid_points, &
                        reduced_rg_dim, kr_grid_width_res, kr_grid_ampl_res, k_space_dim, &
                        delta_l_max, gauss_int_nodes_Ntheta, gauss_int_nodes_Nx, gauss_int_nodes_Nxp, &
                        r_plas, r_min, width_res, ampl_res, hrmax_scaling

    num_args = command_argument_count()
    if (num_args > 1) then
        write(*,*) 'Too many arguments'
        stop
    else if (num_args == 1) then
        allocate(args(num_args))  ! I've omitted checking the return status of the allocation 
        do ix = 1, num_args
            call get_command_argument(ix,args(ix))
            print *, 'Argument ', ix, ': ', args(ix)
        end do
        nml_config_path = trim(args(1))
        write(*,*) 'Namelist path provided: ', nml_config_path
    end if

    open(unit = 77, file = trim(nml_config_path))
    read(unit = 77, nml = KIM_CONFIG)
    read(unit = 77, nml = KIM_SETUP)
    read(unit = 77, nml = KIM_GRID)
    close(unit = 77)

    write(output_path, '(A,A,I0,A,I0,A)') trim(output_path), '/m', m_mode, '_n', n_mode, '/'
    inquire(file=trim(output_path), exist=ex)
    if (.not. ex) then
        call system('mkdir -p '//trim(output_path))
    end if

    write(*,*) '+ + + + + + + + KIM + + + + + + + + + + + + + + + +'
    write(*,*) ' type of run = ', type_of_run
    write(*,*) ' - - - - - - - - - - - - - - - - - - - - - - - -'
    write(*,*) 'Configuration namelist'
    write(*,*) '  profile_location           = ', trim(profile_location)
    write(*,*) '  output_path                = ', trim(output_path)
    write(*,*) '  hdf5_input                 = ', hdf5_input
    write(*,*) '  hdf5_output                = ', hdf5_output
    write(*,*) '  fdebug                     = ', fdebug
    write(*,*) '  fstatus                    = ', fstatus
    write(*,*) '  number_of_ion_species      = ', number_of_ion_species
    write(*,*) '  read_species_from_namelist = ', read_species_from_namelist
    write(*,*) '  artificial_debye_case      = ', artificial_debye_case
    write(*,*) '  kernel_debye_case          = ', kernel_debye_case
    write(*,*) '  collision_model            = ', collision_model
    write(*,*) ' - - - - - - - - - - - - - - - - - - - - - - - -'
    write(*,*) 'Setup namelist'
    write(*,*) '  btor                  = ', btor
    write(*,*) '  R0                    = ', R0
    write(*,*) '  m_mode                = ', m_mode
    write(*,*) '  n_mode                = ', n_mode
    write(*,*) '  omega                 = ', omega
    write(*,*) '  spline_base           = ', spline_base
    write(*,*) '  type_br_field         = ', type_br_field
    write(*,*) '  collisions_off        = ', collisions_off
    write(*,*) '  eps_reg               = ', eps_reg
    write(*,*) '  set_profiles_constant = ', set_profiles_constant
    write(*,*) ' - - - - - - - - - - - - - - - - - - - - - - - -'
    write(*,*) 'Grid namelist'
    write(*,*) '  r_min                 = ', r_min
    write(*,*) '  r_plas                = ', r_plas
    write(*,*) '  width_res             = ', width_res
    write(*,*) '  ampl_res              = ', ampl_res
    write(*,*) '  hrmax_scaling         = ', hrmax_scaling
    write(*,*) '  k_space_dim           = ', k_space_dim
    write(*,*) '  l_space_dim           = ', l_space_dim
    write(*,*) '  reduce_r              = ', reduce_r
    write(*,*) '  reduced_rg_dim        = ', reduced_rg_dim
    write(*,*) '  num_gengrid_points    = ', num_gengrid_points
    write(*,*) '  grid_spacing          = ', grid_spacing
    write(*,*) '  kr_grid_width_res     = ', kr_grid_width_res
    write(*,*) '  kr_grid_ampl_res      = ', kr_grid_ampl_res
    write(*,*) '  delta_l_max           = ', delta_l_max
    write(*,*) '  gauss_int_nodes_Ntheta= ', gauss_int_nodes_Ntheta
    write(*,*) '  gauss_int_nodes_Nx    = ', gauss_int_nodes_Nx
    write(*,*) '  gauss_int_nodes_Nxp   = ', gauss_int_nodes_Nxp
    write(*,*) ' - - - - - - - - - - - - - - - - - - - - - - - -'
    write(*,*) ' - - - - - - - - - - - - - - - - - - - - - - - -'

end subroutine