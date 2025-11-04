module IO_collection_m

    use KIM_kinds_m, only: dp
    use KAMEL_hdf5_tools, only: HID_T

    implicit none

    integer(HID_T) :: h5id

    contains

    subroutine initialize_hdf5_output()

        use KIM_kinds_m, only: dp
        use KAMEL_hdf5_tools, only: h5_init, h5_open_rw, h5overwrite, h5_create
        use config_m, only: output_path, h5_out_file

        implicit none

        logical :: ex

        inquire(file=trim(output_path), exist=ex)
        if (.not. ex) then
            call system('mkdir -p '//trim(output_path))
        end if

        ! remove existing HDF5 output file
        inquire(file=trim(output_path)//trim(h5_out_file), exist=ex)
        if (ex) then
            call system('rm '//trim(output_path)//trim(h5_out_file))
            call h5_create(trim(output_path)//trim(h5_out_file), h5id)
        end if
        
        h5overwrite = .true.
        CALL h5_init()
        CALL h5_open_rw(trim(output_path)//trim(h5_out_file), h5id)

    end subroutine initialize_hdf5_output

    subroutine deinitialize_hdf5_output()

        use KAMEL_hdf5_tools, only: h5_close

        implicit none

        call h5_close(h5id)

    end subroutine deinitialize_hdf5_output

    subroutine write_KIM_namelist_to_hdf5()

        implicit none

        call write_config_namelist_to_hdf5()
        call write_io_namelist_to_hdf5()
        call write_setup_namelist_to_hdf5()
        call write_grid_namelist_to_hdf5()

    end subroutine write_KIM_namelist_to_hdf5

    subroutine write_config_namelist_to_hdf5()

        use KAMEL_hdf5_tools, only: HID_T, h5_define_group, h5_obj_exists, h5_add, h5_close
        use config_m

        implicit none

        logical :: ex
        integer(HID_T) :: h5grpid

        call h5_obj_exists(h5id, 'config/', ex)
        if (.not. ex) then
            call h5_define_group(h5id, 'config/', h5grpid)
        end if

        call h5_add(h5grpid, 'config/number_of_ion_species', number_of_ion_species, &
            'Number of ion species in the simulation', '1')
        call h5_add(h5grpid, 'config/artificial_debye_case', artificial_debye_case, &
            'Switch for Debye case only. Deactivates everything else. Was used for benchmarking.', '1')
        call h5_add(h5grpid, 'config/type_of_run', trim(type_of_run), &
            'Type of run: electrostatic, FLR2_benchmark, etc.', 'str')
        call h5_add(h5grpid, 'config/collision_model', trim(collision_model), &
            'Type of collision model used in the simulation.', 'str')
        call h5_add(h5grpid, 'config/read_species_from_namelist', read_species_from_namelist, &
            'Logical switch to read species from namelist or use default deuterium plasma.', 'true/false')
        call h5_add(h5grpid, 'config/turn_off_ions', turn_off_ions, &
            'If true, only the first species (electrons) is considered in calculations.', 'true/false')
        call h5_add(h5grpid, 'config/turn_off_electrons', turn_off_electrons, &
            'If true, ions only simulation.', 'true/false')
        call h5_add(h5grpid, 'config/plasma_type', trim(plasma_type), &
            'Type of plasma: H for hydrogen, D for deuterium.', 'str')
        call h5_add(h5grpid, 'config/rescale_density', rescale_density, &
            'Logical switch to rescale density.', 'true/false')
        call h5_add(h5grpid, 'config/number_density_rescale', number_density_rescale, &
            'Factor by which to rescale number density.', 'float')

        call h5_close(h5grpid)

    end subroutine write_config_namelist_to_hdf5


    subroutine write_io_namelist_to_hdf5()

        use KAMEL_hdf5_tools, only: HID_T, h5_define_group, h5_obj_exists, h5_add, h5_close
        use config_m

        implicit none

        logical :: ex
        integer(HID_T) :: h5grpid

        call h5_obj_exists(h5id, 'io/', ex)
        if (.not. ex) then
            call h5_define_group(h5id, 'io/', h5grpid)
        end if

        call h5_add(h5grpid, 'io/profile_location', trim(profile_location), &
            'Location of profile data.', 'str')
        call h5_add(h5grpid, 'io/hdf5_input', hdf5_input, &
            'Logical switch for HDF5 input.', 'true/false')
        call h5_add(h5grpid, 'io/hdf5_output', hdf5_output, &
            'Logical switch for HDF5 output.', 'true/false')
        call h5_add(h5grpid, 'io/fdebug', fdebug, &
            'Logical switch for debug output.', 'true/false')
        call h5_add(h5grpid, 'io/fstatus', fstatus, &
            'Logical switch for status output.', 'true/false')
        call h5_add(h5grpid, 'io/output_path', trim(output_path), &
            'Path for output files.', 'str')
        call h5_add(h5grpid, 'io/calculate_asymptotics', calculate_asymptotics, &
            'Logical switch to calculate asymptotics.', 'true/false')
        call h5_add(h5grpid, 'io/fdiagnostics', fdiagnostics, &
            'Diagnostics output level.', 'integer')
        call h5_add(h5grpid, 'io/h5_out_file', trim(h5_out_file), &
            'Name of the HDF5 output file.', 'str')

        call h5_close(h5grpid)

    end subroutine write_io_namelist_to_hdf5 

    subroutine write_setup_namelist_to_hdf5()

        use KAMEL_hdf5_tools, only: HID_T, h5_define_group, h5_obj_exists, h5_add, h5_close
        use setup_m

        implicit none

        logical :: ex
        integer(HID_T) :: h5grpid

        call h5_obj_exists(h5id, 'setup/', ex)
        if (.not. ex) then
            call h5_define_group(h5id, 'setup/', h5grpid)
        end if

        call h5_add(h5grpid, 'setup/btor', btor, &
            'Toroidal magnetic field at major radius R0.', 'float')
        call h5_add(h5grpid, 'setup/R0', R0, &
            'Major radius of the magnetic axis.', 'float')
        call h5_add(h5grpid, 'setup/m_mode', m_mode, &
            'Poloidal mode number.', 'integer')
        call h5_add(h5grpid, 'setup/n_mode', n_mode, &
            'Toroidal mode number.', 'integer')
        call h5_add(h5grpid, 'setup/omega', omega, &
            'Angular frequency of the perturbation mode.', 'float')
        call h5_add(h5grpid, 'setup/spline_base', spline_base, &
            'Base for FEM basis functions.', 'integer')
        call h5_add(h5grpid, 'setup/type_br_field', type_br_field, &
            'Integer type of delta Br.', '1')
        call h5_add(h5grpid, 'setup/collisions_off', collisions_off, &
            'Logical switch to turn off collisions.', 'true/false')
        call h5_add(h5grpid, 'setup/set_profiles_constant', set_profiles_constant, &
            'Integer switch for setting (some) profiles constant.', '1')
        call h5_add(h5grpid, 'setup/bc_type', bc_type, &
            'Integer type of boundary condition.', '1')

        call h5_close(h5grpid)

    end subroutine write_setup_namelist_to_hdf5

    subroutine write_grid_namelist_to_hdf5()

        use KAMEL_hdf5_tools, only: HID_T, h5_define_group, h5_obj_exists, h5_add, h5_close
        use grid_m

        implicit none

        logical :: ex
        integer(HID_T) :: h5grpid

        call h5_obj_exists(h5id, 'grid/', ex)
        if (.not. ex) then
            call h5_define_group(h5id, 'grid/', h5grpid)
        end if

        call h5_add(h5grpid, 'grid/grid_spacing_rg', trim(grid_spacing_rg), &
            'Grid spacing mode for rg grid.', 'str')
        call h5_add(h5grpid, 'grid/grid_spacing_xl', trim(grid_spacing_xl), &
            'Grid spacing mode for xl grid.', 'str')
        call h5_add(h5grpid, 'grid/l_space_dim', l_space_dim, &
            'Dimension of spline grid.', 'integer')
        call h5_add(h5grpid, 'grid/rg_space_dim', rg_space_dim, &
            'Dimension of rg grid.', 'integer')
        call h5_add(h5grpid, 'grid/theta_integration', trim(theta_integration), &
            'Theta integration method.', 'str')
        call h5_add(h5grpid, 'grid/theta_integration_method', trim(theta_integration_method), &
            'Theta integration method details.', 'str')
        call h5_add(h5grpid, 'grid/Larmor_skip_factor', Larmor_skip_factor, &
            'Larmor skip factor.', 'float')
        call h5_add(h5grpid, 'grid/gauss_int_nodes_Ntheta', gauss_int_nodes_Ntheta, &
            'Number of Gauss integration nodes in theta.', 'integer')
        call h5_add(h5grpid, 'grid/gauss_int_nodes_Nx', gauss_int_nodes_Nx, &
            'Number of Gauss integration nodes in x.', 'integer')
        call h5_add(h5grpid, 'grid/gauss_int_nodes_Nxp', gauss_int_nodes_Nxp, &
            'Number of Gauss integration nodes in x prime.', 'integer')
        call h5_add(h5grpid, 'grid/r_plas', r_plas, &
            'Plasma radius.', 'float')
        call h5_add(h5grpid, 'grid/r_min', r_min, &
            'Minimum radius.', 'float')
        call h5_add(h5grpid, 'grid/width_res', width_res, &
            'Width of resolution region.', 'float')
        call h5_add(h5grpid, 'grid/ampl_res', ampl_res, &
            'Amplitude of resolution region.', 'float')
        call h5_add(h5grpid, 'grid/hrmax_scaling', hrmax_scaling, &
            'Scaling factor for maximum grid spacing.', 'float')
        call h5_add(h5grpid, 'grid/rkf45_atol', rkf45_atol, &
            'Absolute tolerance for RKF45 integration.', 'float')
        call h5_add(h5grpid, 'grid/rkf45_rtol', rkf45_rtol, &
            'Relative tolerance for RKF45 integration.', 'float')
        call h5_add(h5grpid, 'grid/kernel_taper_skip_threshold', kernel_taper_skip_threshold, &
            'Threshold for kernel tapering skip.', 'float')
        call h5_add(h5grpid, 'grid/quadpack_algorithm', quadpack_algorithm, &
            'Algorithm used in Quadpack integration.', 'integer')
        call h5_add(h5grpid, 'grid/quadpack_key', quadpack_key, &
            'Key for Quadpack integration.', 'integer')
        call h5_add(h5grpid, 'grid/quadpack_limit', quadpack_limit, &
            'Limit for Quadpack integration.', 'integer')
        call h5_add(h5grpid, 'grid/quadpack_epsabs', quadpack_epsabs, &
            'Absolute error tolerance for Quadpack integration.', 'float')
        call h5_add(h5grpid, 'grid/quadpack_epsrel', quadpack_epsrel, &
            'Relative error tolerance for Quadpack integration.', 'float')
        call h5_add(h5grpid, 'grid/quadpack_use_u_substitution', quadpack_use_u_substitution, &
            'Logical switch for using u substitution in Quadpack integration.', 'true/false')

        call h5_close(h5grpid)

    end subroutine write_grid_namelist_to_hdf5


    subroutine write_matrix(filename, A, nx, ny, comment, unit)
        
        use KIM_kinds_m, only: dp
        use config_m, only: output_path, hdf5_output
        use KAMEL_hdf5_tools

        implicit none

        character(*), intent(in) :: filename
        character(*), optional :: comment, unit
        real(dp), intent(in) :: A(nx, ny)
        integer, intent(in) :: nx, ny

        integer :: i, j

        if (hdf5_output) then

            if (.not. present(comment)) then
                comment = ''
            end if

            if (.not. present(unit)) then
                comment = ''
            end if

            call h5_add(h5id, filename, A, [1, 1], [nx, ny], comment, unit)
            return

        end if

        open(unit=10, file=trim(output_path)//filename//'.dat', status='replace')

        do i = 1, nx
            write(10, *) (A(i,j), j = 1, ny)
        end do

        close(10)

    end subroutine write_matrix


    subroutine write_profile(x, y, n, filename, comment, unit)

        use KIM_kinds_m, only: dp
        use config_m, only: output_path, hdf5_output
        use KAMEL_hdf5_tools, only: h5_add

        implicit none

        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), y(n)
        character(len=*), intent(in) :: filename
        character(*), optional :: comment, unit

        integer :: i

        if (hdf5_output) then

            if (.not. present(comment)) then
                comment = ''
            end if

            if (.not. present(unit)) then
                comment = ''
            end if

            call h5_add(h5id, filename, y, [1], [n], comment, unit)
            return
        end if

        open(unit=10, file=trim(output_path)//filename//'.dat', status='replace', action='write')

        do i = 1, n
            write(10, *) x(i), y(i)
        end do

        close(10)

    end subroutine write_profile


    subroutine write_complex_profile(x, y, n, filename, comment, unit)

        use KIM_kinds_m, only: dp
        use config_m, only: output_path, hdf5_output
        use KAMEL_hdf5_tools, only: h5_add

        implicit none

        integer, intent(in) :: n
        real(dp), intent(in) :: x(n)
        complex(dp), intent(in) :: y(n)
        character(len=*), intent(in) :: filename
        character(*), optional :: comment, unit

        integer :: i

        if (hdf5_output) then

            if (.not. present(comment)) then
                comment = ''
            end if

            if (.not. present(unit)) then
                comment = ''
            end if

            call h5_add(h5id, filename, y, [1], [n], comment, unit)
            return
        end if


        open(unit=10, file=trim(output_path)//filename//'.dat', status='replace', action='write')

        ! Write data as two columns
        do i = 1, n
            write(10, *) x(i), real(y(i)), dimag(y(i))
        end do

        close(10)

    end subroutine write_complex_profile

    subroutine write_complex_profile_abs(x, y, n, filename, comment, unit)

        use KIM_kinds_m, only: dp
        use config_m, only: output_path, hdf5_output

        implicit none

        integer, intent(in) :: n
        real(dp), intent(in) :: x(n)
        complex(dp), intent(in) :: y(n)
        character(len=*), intent(in) :: filename
        character(*), optional :: comment, unit

        integer :: i

        if (hdf5_output) then
            call write_complex_profile(x, y, n, filename, comment, unit)
        end if

        open(unit=10, file=trim(output_path)//filename//'.dat', status='replace', action='write')

        ! Write data as two columns
        do i = 1, n
            write(10, *) x(i), real(y(i)), dimag(y(i)), abs(y(i))
        end do

        close(10)

    end subroutine write_complex_profile_abs


    subroutine plot_2D(datafile)

        implicit none

        character(*), intent(in) :: datafile
        character(len=200) :: cmd

        ! Create command string
        write(cmd, '(A)') 'gnuplot -persist -e "splot '''//trim(datafile)//''' with lines"'

        call execute_command_line(trim(cmd))

    end subroutine plot_2D

    subroutine plot_1D(datafile)

        implicit none

        character(*), intent(in) :: datafile
        character(len=256) :: cmd

        write(cmd, '(A)') 'gnuplot -persist -e "plot '''//trim(datafile)//''' using 1:2 with lines"'

        call execute_command_line(trim(cmd))

    end subroutine

    subroutine plot_profile(x,y)

        use KIM_kinds_m, only: dp
        use config_m, only: hdf5_output

        implicit none

        real(dp), intent(in) :: x(:), y(:)

        if (.not. hdf5_output) then
            call write_profile(x, y, size(x), 'profile.dat')
            call plot_1D('profile.dat')
            call remove_file('profile.dat')
        else 
            print *, 'Plotting not supported for HDF5 output mode.'
        end if

    end subroutine

    subroutine plot_complex_1D(datafile)

        implicit none

        character(*), intent(in) :: datafile
        character(len=256) :: cmd

        write(cmd, '(A)') 'gnuplot -persist -e "plot '''//trim(datafile)//''' using 1:2 with lines; plot '''//trim(datafile)//''' using 1:3 with lines"'

        call execute_command_line(trim(cmd))

    end subroutine

    subroutine plot_matrix(datafile)
        implicit none
        character(*), intent(in) :: datafile
        character(len=300) :: cmd

        write(cmd, '(A)') 'gnuplot -persist -e "plot '''//trim(datafile)//''' matrix with image"'
        call execute_command_line(trim(cmd))
    end subroutine plot_matrix

    subroutine plot_1D_labeled(datafile, xlabel, ylabel, title)

        implicit none

        character(*), intent(in) :: datafile
        character(*), intent(in) :: xlabel, ylabel, title
        character(len=1024) :: cmd
        character(len=1024) :: plot_cmd

        plot_cmd = 'set xlabel ''' // trim(xlabel) // '''; ' // &
                'set ylabel ''' // trim(ylabel) // '''; ' // &
                'set title '''  // trim(title)  // '''; ' // &
                'plot ''' // trim(datafile) // ''' using 1:2 with lines;'

        write(cmd, '(A)') 'gnuplot -persist -e "' // trim(plot_cmd) // '"'

        call execute_command_line(trim(cmd))

    end subroutine

    subroutine remove_file(filename)
    
        implicit none

        character(*), intent(in) :: filename
        character(len=300) :: cmd

        ! Build the command to remove the file
        write(cmd, '(A)') 'rm -f "' // trim(filename) // '"'

        ! Execute the command
        call execute_command_line(trim(cmd))

    end subroutine

    subroutine create_output_directories

        use config_m, only: output_path, hdf5_output
        use KAMEL_hdf5_tools, only: h5_obj_exists, h5_define_group, HID_T

        implicit none

        logical :: ex
        integer(HID_T) :: h5grpid

        if (hdf5_output) then

            call h5_obj_exists(h5id, 'fields/', ex)
            if (.not. ex) then
                call h5_define_group(h5id, 'fields/', h5grpid)
            end if

            call h5_obj_exists(h5id, 'kernel/', ex)
            if (.not. ex) then
                call h5_define_group(h5id, 'kernel/', h5grpid)
            end if

            call h5_obj_exists(h5id, 'backs/', ex)
            if (.not. ex) then
                call h5_define_group(h5id, 'backs/', h5grpid)
            end if

            call h5_obj_exists(h5id, 'grid/', ex)
            if (.not. ex) then
                call h5_define_group(h5id, 'grid/', h5grpid)
            end if

            return
        end if

        inquire(file=trim(output_path)//'fields', exist=ex)
        if (.not. ex) then
            call system('mkdir -p '//trim(output_path)//'fields')
        end if
        inquire(file=trim(output_path)//'kernel', exist=ex)
        if (.not. ex) then
            call system('mkdir -p '//trim(output_path)//'kernel')
        end if
        inquire(file=trim(output_path)//'backs', exist=ex)
        if (.not. ex) then
            call system('mkdir -p '//trim(output_path)//'backs')
        end if
        inquire(file = trim(output_path)//'grid', exist = ex)
        if (.not. ex) then
            call system('mkdir -p '//trim(output_path)//'grid')
        end if

    end subroutine


end module
