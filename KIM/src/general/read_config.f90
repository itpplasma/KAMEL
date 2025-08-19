subroutine read_config

    use config
    use constants
    use setup
    use grid
    use poisson_solver, only: solve_poisson
    use config_display, only: display_kim_configuration

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

    namelist /KIM_GRID/ reduce_r, grid_spacing, l_space_dim, theta_integration, &
                        reduced_rg_dim, kr_grid_width_res, kr_grid_ampl_res, k_space_dim, &
                        Larmor_skip_factor, gauss_int_nodes_Ntheta, gauss_int_nodes_Nx, gauss_int_nodes_Nxp, &
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

    if (collisions_off .and. collision_model == "FokkerPlanck") then
        write(*,*) 'Error: collision_model is set to "FokkerPlanck" but collisions_off is true.'
        write(*,*) 'Please set collisions_off to false or change collision_model.'
        stop
    end if

    write(output_path, '(A,A,I0,A,I0,A)') trim(output_path), '/m', m_mode, '_n', n_mode, '/'
    inquire(file=trim(output_path), exist=ex)
    if (.not. ex) then
        call system('mkdir -p '//trim(output_path))
    end if

    ! Display formatted configuration
    call display_kim_configuration()

end subroutine
