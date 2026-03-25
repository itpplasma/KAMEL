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

        use KAMEL_hdf5_tools, only: HID_T, h5_define_group, h5_obj_exists, h5_add, h5_close_group
        use config_m

        implicit none

        logical :: ex
        integer(HID_T) :: h5grpid

        call h5_obj_exists(h5id, 'config/', ex)
        if (.not. ex) then
            call h5_define_group(h5id, 'config/', h5grpid)
        end if

        call h5_add(h5grpid, 'number_of_ion_species', number_of_ion_species, &
            'Number of ion species in the simulation', '1')
        call h5_add(h5grpid, 'artificial_debye_case', artificial_debye_case, &
            'Switch for Debye case only. Deactivates everything else. Was used for benchmarking.', '1')
        call h5_add(h5grpid, 'type_of_run', trim(type_of_run), &
            'Type of run: electrostatic, FLR2_benchmark, etc.', 'str')
        call h5_add(h5grpid, 'collision_model', trim(collision_model), &
            'Type of collision model used in the simulation.', 'str')
        call h5_add(h5grpid, 'read_species_from_namelist', read_species_from_namelist, &
            'Logical switch to read species from namelist or use default deuterium plasma.', 'true/false')
        call h5_add(h5grpid, 'turn_off_ions', turn_off_ions, &
            'If true, only the first species (electrons) is considered in calculations.', 'true/false')
        call h5_add(h5grpid, 'turn_off_electrons', turn_off_electrons, &
            'If true, ions only simulation.', 'true/false')
        call h5_add(h5grpid, 'plasma_type', trim(plasma_type), &
            'Type of plasma: H for hydrogen, D for deuterium.', 'str')
        call h5_add(h5grpid, 'rescale_density', rescale_density, &
            'Logical switch to rescale density.', 'true/false')
        call h5_add(h5grpid, 'number_density_rescale', number_density_rescale, &
            'Factor by which to rescale number density.', 'float')

        call h5_close_group(h5grpid)

    end subroutine write_config_namelist_to_hdf5


    subroutine write_io_namelist_to_hdf5()

        use KAMEL_hdf5_tools, only: HID_T, h5_define_group, h5_obj_exists, h5_add, h5_close_group
        use config_m

        implicit none

        logical :: ex
        integer(HID_T) :: h5grpid

        call h5_obj_exists(h5id, 'io/', ex)
        if (.not. ex) then
            call h5_define_group(h5id, 'io/', h5grpid)
        end if

        call h5_add(h5grpid, 'profile_location', trim(profile_location), &
            'Location of profile data.', 'str')
        call h5_add(h5grpid, 'hdf5_input', hdf5_input, &
            'Logical switch for HDF5 input.', 'true/false')
        call h5_add(h5grpid, 'hdf5_output', hdf5_output, &
            'Logical switch for HDF5 output.', 'true/false')
        call h5_add(h5grpid, 'log_level', log_level, &
            'Log verbosity level', 'i')
        call h5_add(h5grpid, 'data_verbosity', data_verbosity, &
            'Data output verbosity', 'i')
        call h5_add(h5grpid, 'output_path', trim(output_path), &
            'Path for output files.', 'str')
        call h5_add(h5grpid, 'calculate_asymptotics', calculate_asymptotics, &
            'Logical switch to calculate asymptotics.', 'true/false')
        call h5_add(h5grpid, 'h5_out_file', trim(h5_out_file), &
            'Name of the HDF5 output file.', 'str')

        call h5_close_group(h5grpid)

    end subroutine write_io_namelist_to_hdf5

    subroutine write_setup_namelist_to_hdf5()

        use KAMEL_hdf5_tools, only: HID_T, h5_define_group, h5_obj_exists, h5_add, h5_close_group
        use setup_m

        implicit none

        logical :: ex
        integer(HID_T) :: h5grpid

        call h5_obj_exists(h5id, 'setup/', ex)
        if (.not. ex) then
            call h5_define_group(h5id, 'setup/', h5grpid)
        end if

        call h5_add(h5grpid, 'btor', btor, &
            'Toroidal magnetic field at major radius R0.', 'float')
        call h5_add(h5grpid, 'B_ref', btor, &
            'Reference magnetic field at magnetic axis, used for effective radius r = sqrt(2 psi_tor / B_ref).', 'G')
        call h5_add(h5grpid, 'R0', R0, &
            'Major radius of the magnetic axis.', 'float')
        call h5_add(h5grpid, 'm_mode', m_mode, &
            'Poloidal mode number.', 'integer')
        call h5_add(h5grpid, 'n_mode', n_mode, &
            'Toroidal mode number.', 'integer')
        call h5_add(h5grpid, 'omega', omega, &
            'Angular frequency of the perturbation mode.', 'float')
        call h5_add(h5grpid, 'spline_base', spline_base, &
            'Base for FEM basis functions.', 'integer')
        call h5_add(h5grpid, 'type_br_field', type_br_field, &
            'Integer type of delta Br.', '1')
        call h5_add(h5grpid, 'collisions_off', collisions_off, &
            'Logical switch to turn off collisions.', 'true/false')
        call h5_add(h5grpid, 'set_profiles_constant', set_profiles_constant, &
            'Integer switch for setting (some) profiles constant.', '1')
        call h5_add(h5grpid, 'bc_type', bc_type, &
            'Integer type of boundary condition.', '1')
        call h5_add(h5grpid, 'mphi_max', mphi_max, &
            'Max. number of cyclotron harmonics used', '1')

        call h5_close_group(h5grpid)

    end subroutine write_setup_namelist_to_hdf5

    subroutine write_grid_namelist_to_hdf5()

        use KAMEL_hdf5_tools, only: HID_T, h5_define_group, h5_obj_exists, h5_add, h5_close_group
        use grid_m

        implicit none

        logical :: ex
        integer(HID_T) :: h5grpid

        call h5_obj_exists(h5id, 'grid/', ex)
        if (.not. ex) then
            call h5_define_group(h5id, 'grid/', h5grpid)
        end if

        call h5_add(h5grpid, 'grid_spacing_rg', trim(grid_spacing_rg), &
            'Grid spacing mode for rg grid.', 'str')
        call h5_add(h5grpid, 'grid_spacing_xl', trim(grid_spacing_xl), &
            'Grid spacing mode for xl grid.', 'str')
        call h5_add(h5grpid, 'l_space_dim', l_space_dim, &
            'Dimension of spline grid.', 'integer')
        call h5_add(h5grpid, 'rg_space_dim', rg_space_dim, &
            'Dimension of rg grid.', 'integer')
        call h5_add(h5grpid, 'theta_integration', trim(theta_integration), &
            'Theta integration method.', 'str')
        call h5_add(h5grpid, 'theta_integration_method', trim(theta_integration_method), &
            'Theta integration method details.', 'str')
        call h5_add(h5grpid, 'Larmor_skip_factor', Larmor_skip_factor, &
            'Larmor skip factor.', 'float')
        call h5_add(h5grpid, 'gauss_int_nodes_Ntheta', gauss_int_nodes_Ntheta, &
            'Number of Gauss integration nodes in theta.', 'integer')
        call h5_add(h5grpid, 'gauss_int_nodes_Nx', gauss_int_nodes_Nx, &
            'Number of Gauss integration nodes in x.', 'integer')
        call h5_add(h5grpid, 'gauss_int_nodes_Nxp', gauss_int_nodes_Nxp, &
            'Number of Gauss integration nodes in x prime.', 'integer')
        call h5_add(h5grpid, 'r_plas', r_plas, &
            'Plasma radius.', 'float')
        call h5_add(h5grpid, 'r_min', r_min, &
            'Minimum radius.', 'float')
        call h5_add(h5grpid, 'width_res', width_res, &
            'Width of resolution region.', 'float')
        call h5_add(h5grpid, 'ampl_res', ampl_res, &
            'Amplitude of resolution region.', 'float')
        call h5_add(h5grpid, 'hrmax_scaling', hrmax_scaling, &
            'Scaling factor for maximum grid spacing.', 'float')
        call h5_add(h5grpid, 'rkf45_atol', rkf45_atol, &
            'Absolute tolerance for RKF45 integration.', 'float')
        call h5_add(h5grpid, 'rkf45_rtol', rkf45_rtol, &
            'Relative tolerance for RKF45 integration.', 'float')
        call h5_add(h5grpid, 'kernel_taper_skip_threshold', kernel_taper_skip_threshold, &
            'Threshold for kernel tapering skip.', 'float')
        call h5_add(h5grpid, 'quadpack_algorithm', quadpack_algorithm, &
            'Algorithm used in Quadpack integration.', 'integer')
        call h5_add(h5grpid, 'quadpack_key', quadpack_key, &
            'Key for Quadpack integration.', 'integer')
        call h5_add(h5grpid, 'quadpack_limit', quadpack_limit, &
            'Limit for Quadpack integration.', 'integer')
        call h5_add(h5grpid, 'quadpack_epsabs', quadpack_epsabs, &
            'Absolute error tolerance for Quadpack integration.', 'float')
        call h5_add(h5grpid, 'quadpack_epsrel', quadpack_epsrel, &
            'Relative error tolerance for Quadpack integration.', 'float')
        call h5_add(h5grpid, 'quadpack_use_u_substitution', quadpack_use_u_substitution, &
            'Logical switch for using u substitution in Quadpack integration.', 'true/false')

        call h5_close_group(h5grpid)

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

    subroutine write_complex_profile_abs(x, y, n, filename, comment, unit, output_dir)

        use KIM_kinds_m, only: dp
        use config_m, only: output_path, hdf5_output

        implicit none

        integer, intent(in) :: n
        real(dp), intent(in) :: x(n)
        complex(dp), intent(in) :: y(n)
        character(len=*), intent(in) :: filename
        character(*), optional :: comment, unit
        character(len=*), intent(in), optional :: output_dir

        integer :: i
        character(len=512) :: out_dir

        if (hdf5_output) then
            call write_complex_profile(x, y, n, filename, comment, unit)
            return
        end if

        ! Use custom output directory if provided, otherwise use default
        if (present(output_dir)) then
            out_dir = output_dir
        else
            out_dir = output_path
        end if

        open(unit=10, file=trim(out_dir)//filename//'.dat', status='replace', action='write')

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

    subroutine ensure_dispersion_dir_exists()
        !> Ensure the dispersion output directory exists, creating it if necessary

        use config_m, only: dispersion_output_path

        implicit none

        logical :: ex

        inquire(file=trim(dispersion_output_path), exist=ex)
        if (.not. ex) then
            call system('mkdir -p '//trim(dispersion_output_path))
        end if

    end subroutine ensure_dispersion_dir_exists

    function itoa(i) result(res)

        implicit none

        character(:),allocatable :: res
        integer,intent(in) :: i
        character(range(i)+2) :: tmp

        write(tmp,'(i0)') i

        res = trim(tmp)

    end function


    subroutine write_roots_with_multiplicities(x_grid, n_grid, zeros, fzeros, &
            multiplicities, n_roots_per_point, filename, comment)
        !> Write complex roots with their multiplicities at each grid point.
        !> Supports both text file and HDF5 output.
        !>
        !> For text output: writes a file with columns:
        !>   grid_index, x_position, Re(zero), Im(zero), Re(f(zero)), Im(f(zero)), multiplicity
        !>
        !> For HDF5 output: writes three datasets:
        !>   - filename/zeros: complex array (max_roots, n_grid)
        !>   - filename/fzeros: complex array (max_roots, n_grid)
        !>   - filename/multiplicities: integer array (max_roots, n_grid)
        !>   - filename/n_roots: integer array (n_grid) with actual count per point

        use KIM_kinds_m, only: dp
        use config_m, only: output_path, hdf5_output
        use KAMEL_hdf5_tools, only: h5_add, h5_define_group, h5_obj_exists, h5_close_group, HID_T

        implicit none

        integer, intent(in) :: n_grid
        real(dp), intent(in) :: x_grid(n_grid)
        complex(dp), intent(in) :: zeros(:,:)        ! (max_roots, n_grid)
        complex(dp), intent(in) :: fzeros(:,:)       ! (max_roots, n_grid)
        integer, intent(in) :: multiplicities(:,:)   ! (max_roots, n_grid)
        integer, intent(in) :: n_roots_per_point(n_grid)
        character(len=*), intent(in) :: filename
        character(len=*), intent(in), optional :: comment

        integer :: i, j, max_roots
        integer(HID_T) :: h5grpid
        logical :: ex
        character(len=256) :: comment_str

        max_roots = size(zeros, 1)

        if (present(comment)) then
            comment_str = comment
        else
            comment_str = 'ZEAL root finding results'
        end if

        if (hdf5_output) then
            ! Create group for this output
            call h5_obj_exists(h5id, trim(filename)//'/', ex)
            if (.not. ex) then
                call h5_define_group(h5id, trim(filename)//'/', h5grpid)
            end if

            call h5_add(h5grpid, 'zeros', zeros, [1, 1], [max_roots, n_grid], &
                trim(comment_str)//' - complex zeros', 'cm^{-1}')
            call h5_add(h5grpid, 'fzeros', fzeros, [1, 1], [max_roots, n_grid], &
                trim(comment_str)//' - function values at zeros', '1')
            call h5_add(h5grpid, 'multiplicities', multiplicities, [1, 1], [max_roots, n_grid], &
                trim(comment_str)//' - multiplicities', '1')
            call h5_add(h5grpid, 'n_roots', n_roots_per_point, [1], [n_grid], &
                'Number of distinct roots at each grid point', '1')
            call h5_add(h5grpid, 'x_grid', x_grid, [1], [n_grid], &
                'Radial grid positions', 'cm')

            call h5_close_group(h5grpid)
            return
        end if

        ! Text file output
        open(unit=10, file=trim(output_path)//trim(filename)//'.dat', status='replace', action='write')

        ! Write header
        write(10, '(A)') '# ' // trim(comment_str)
        write(10, '(A)') '# grid_idx    x_position    Re(zero)    Im(zero)    ' // &
            'Re(f(zero))    Im(f(zero))    multiplicity'

        do j = 1, n_grid
            do i = 1, n_roots_per_point(j)
                write(10, '(I8, 5(ES22.14), I8)') j, x_grid(j), &
                    real(zeros(i,j), dp), aimag(zeros(i,j)), &
                    real(fzeros(i,j), dp), aimag(fzeros(i,j)), &
                    multiplicities(i,j)
            end do
        end do

        close(10)

    end subroutine write_roots_with_multiplicities


    subroutine track_root_branches(zeros, n_roots_per_point, n_grid, max_roots, &
            branch_id, n_branches)
        !> Track roots across grid points by matching closest roots between adjacent points.
        !> Uses greedy nearest-neighbor algorithm to follow root branches.
        !>
        !> Output:
        !>   branch_id(i, j) = branch index for root i at grid point j
        !>                     0 means no root at this position
        !>   n_branches = total number of distinct branches found

        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: n_grid, max_roots
        complex(dp), intent(in) :: zeros(max_roots, n_grid)
        integer, intent(in) :: n_roots_per_point(n_grid)
        integer, intent(out) :: branch_id(max_roots, n_grid)
        integer, intent(out) :: n_branches

        integer :: j, i, k, best_match, current_branch
        real(dp) :: dist, best_dist
        logical :: matched(max_roots)
        integer :: prev_branch(max_roots)  ! branch IDs at previous grid point

        ! Initialize
        branch_id = 0
        n_branches = 0

        ! Handle empty case
        if (n_grid == 0) return
        if (n_roots_per_point(1) == 0) return

        ! Initialize first grid point: each root gets a unique branch ID
        do i = 1, n_roots_per_point(1)
            n_branches = n_branches + 1
            branch_id(i, 1) = n_branches
        end do

        ! Track through remaining grid points
        do j = 2, n_grid
            if (n_roots_per_point(j) == 0) cycle

            ! Store branch IDs from previous point for matching
            prev_branch = 0
            do i = 1, n_roots_per_point(j-1)
                prev_branch(i) = branch_id(i, j-1)
            end do

            matched = .false.

            ! For each root at current grid point, find closest match at previous point
            do i = 1, n_roots_per_point(j)
                best_match = 0
                best_dist = huge(1.0_dp)

                ! Search for closest unmatched root at previous grid point
                do k = 1, n_roots_per_point(j-1)
                    if (matched(k)) cycle

                    dist = abs(zeros(i, j) - zeros(k, j-1))
                    if (dist < best_dist) then
                        best_dist = dist
                        best_match = k
                    end if
                end do

                if (best_match > 0) then
                    ! Continue existing branch
                    branch_id(i, j) = prev_branch(best_match)
                    matched(best_match) = .true.
                else
                    ! Start new branch (more roots than previous point)
                    n_branches = n_branches + 1
                    branch_id(i, j) = n_branches
                end if
            end do
        end do

    end subroutine track_root_branches


    subroutine write_tracked_roots(x_grid, n_grid, zeros, fzeros, multiplicities, &
            n_roots_per_point, branch_id, n_branches, filename, comment, output_dir)
        !> Write roots organized by branch (tracked across grid points).
        !> Each branch represents a continuous root that can be followed radially.
        !>
        !> For text output: writes one file per branch with columns:
        !>   x_position, Re(zero), Im(zero), Re(f(zero)), Im(f(zero)), multiplicity
        !>
        !> For HDF5 output: writes datasets organized by branch index

        use KIM_kinds_m, only: dp
        use config_m, only: output_path, hdf5_output
        use KAMEL_hdf5_tools, only: h5_add, h5_define_group, h5_obj_exists, h5_close_group, HID_T

        implicit none

        integer, intent(in) :: n_grid, n_branches
        real(dp), intent(in) :: x_grid(n_grid)
        complex(dp), intent(in) :: zeros(:,:)
        complex(dp), intent(in) :: fzeros(:,:)
        integer, intent(in) :: multiplicities(:,:)
        integer, intent(in) :: n_roots_per_point(n_grid)
        integer, intent(in) :: branch_id(:,:)
        character(len=*), intent(in) :: filename
        character(len=*), intent(in), optional :: comment
        character(len=*), intent(in), optional :: output_dir

        integer :: i, j, b, max_roots, n_points_in_branch
        integer(HID_T) :: h5grpid, h5branchid
        logical :: ex
        character(len=256) :: comment_str, branch_filename
        character(len=16) :: branch_name
        character(len=512) :: out_dir

        ! Arrays for single branch data
        real(dp), allocatable :: branch_x(:)
        complex(dp), allocatable :: branch_zeros(:), branch_fzeros(:)
        integer, allocatable :: branch_mult(:)

        max_roots = size(zeros, 1)

        ! Use custom output directory if provided, otherwise use default
        if (present(output_dir)) then
            out_dir = output_dir
        else
            out_dir = output_path
        end if

        if (present(comment)) then
            comment_str = comment
        else
            comment_str = 'ZEAL tracked root branches'
        end if

        ! Allocate temporary arrays for branch data
        allocate(branch_x(n_grid))
        allocate(branch_zeros(n_grid))
        allocate(branch_fzeros(n_grid))
        allocate(branch_mult(n_grid))

        if (hdf5_output) then
            ! Create main group
            call h5_obj_exists(h5id, trim(filename)//'/', ex)
            if (.not. ex) then
                call h5_define_group(h5id, trim(filename)//'/', h5grpid)
            end if

            call h5_add(h5grpid, 'n_branches', n_branches, &
                'Total number of tracked branches', '1')
            call h5_add(h5grpid, 'x_grid', x_grid, [1], [n_grid], &
                'Radial grid positions', 'cm')

            ! Write each branch
            do b = 1, n_branches
                write(branch_name, '(A,I0)') 'branch_', b

                ! Collect data for this branch
                n_points_in_branch = 0
                do j = 1, n_grid
                    do i = 1, n_roots_per_point(j)
                        if (branch_id(i, j) == b) then
                            n_points_in_branch = n_points_in_branch + 1
                            branch_x(n_points_in_branch) = x_grid(j)
                            branch_zeros(n_points_in_branch) = zeros(i, j)
                            branch_fzeros(n_points_in_branch) = fzeros(i, j)
                            branch_mult(n_points_in_branch) = multiplicities(i, j)
                            exit  ! Only one root per branch per grid point
                        end if
                    end do
                end do

                if (n_points_in_branch > 0) then
                    call h5_define_group(h5grpid, trim(branch_name)//'/', h5branchid)
                    call h5_add(h5branchid, 'x', branch_x(1:n_points_in_branch), &
                        [1], [n_points_in_branch], 'Radial positions for this branch', 'cm')
                    call h5_add(h5branchid, 'zeros', branch_zeros(1:n_points_in_branch), &
                        [1], [n_points_in_branch], 'Complex zeros along branch', 'cm^{-1}')
                    call h5_add(h5branchid, 'fzeros', branch_fzeros(1:n_points_in_branch), &
                        [1], [n_points_in_branch], 'Function values at zeros', '1')
                    call h5_add(h5branchid, 'multiplicities', branch_mult(1:n_points_in_branch), &
                        [1], [n_points_in_branch], 'Multiplicities along branch', '1')
                    call h5_add(h5branchid, 'n_points', n_points_in_branch, &
                        'Number of grid points in this branch', '1')
                    call h5_close_group(h5branchid)
                end if
            end do

            call h5_close_group(h5grpid)

        else
            ! Text file output: one file per branch
            do b = 1, n_branches
                write(branch_filename, '(A,A,I0,A)') trim(filename), '_branch_', b, '.dat'

                ! Collect data for this branch
                n_points_in_branch = 0
                do j = 1, n_grid
                    do i = 1, n_roots_per_point(j)
                        if (branch_id(i, j) == b) then
                            n_points_in_branch = n_points_in_branch + 1
                            branch_x(n_points_in_branch) = x_grid(j)
                            branch_zeros(n_points_in_branch) = zeros(i, j)
                            branch_fzeros(n_points_in_branch) = fzeros(i, j)
                            branch_mult(n_points_in_branch) = multiplicities(i, j)
                            exit
                        end if
                    end do
                end do

                if (n_points_in_branch > 0) then
                    open(unit=10, file=trim(out_dir)//trim(branch_filename), &
                        status='replace', action='write')

                    write(10, '(A,I0)') '# ' // trim(comment_str) // ' - Branch ', b
                    write(10, '(A)') '# x_position    Re(zero)    Im(zero)    ' // &
                        'Re(f(zero))    Im(f(zero))    multiplicity'

                    do i = 1, n_points_in_branch
                        write(10, '(5(ES22.14), I8)') branch_x(i), &
                            real(branch_zeros(i), dp), aimag(branch_zeros(i)), &
                            real(branch_fzeros(i), dp), aimag(branch_fzeros(i)), &
                            branch_mult(i)
                    end do

                    close(10)
                end if
            end do

            ! Also write a summary file
            open(unit=10, file=trim(out_dir)//trim(filename)//'_summary.dat', &
                status='replace', action='write')
            write(10, '(A)') '# ' // trim(comment_str)
            write(10, '(A,I0)') '# Total branches found: ', n_branches
            write(10, '(A)') '# branch_id    n_points    start_x    end_x'

            do b = 1, n_branches
                n_points_in_branch = 0
                do j = 1, n_grid
                    do i = 1, n_roots_per_point(j)
                        if (branch_id(i, j) == b) then
                            n_points_in_branch = n_points_in_branch + 1
                            if (n_points_in_branch == 1) then
                                branch_x(1) = x_grid(j)  ! start position
                            end if
                            branch_x(2) = x_grid(j)  ! end position (updated each time)
                            exit
                        end if
                    end do
                end do
                if (n_points_in_branch > 0) then
                    write(10, '(I8, I12, 2(ES16.8))') b, n_points_in_branch, &
                        branch_x(1), branch_x(2)
                end if
            end do
            close(10)
        end if

        deallocate(branch_x, branch_zeros, branch_fzeros, branch_mult)

    end subroutine write_tracked_roots


end module
