program test_radial_current_fourier_kernel
    use KIM_kinds_m, only: dp

    implicit none

    call test_i03_is_retained()

    print *, 'All radial-current Fourier-kernel tests PASSED'
    stop 0

contains

    subroutine test_i03_is_retained()
        use config_m, only: nml_config_path, profiles_in_memory
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim
        use species_m, only: evaluate_susceptibility, plasma, set_profiles_from_arrays

        integer, parameter :: npts = 101
        real(dp), parameter :: tolerance = 1.0e-12_dp
        real(dp) :: r(npts), density(npts), te(npts), ti(npts), q(npts), er(npts)
        complex(dp) :: moments(0:3, 0:3)
        class(kim_t), allocatable :: kim
        integer :: ell, j, sp

        call make_test_profiles(r, density, te, ti, q, er)
        call write_test_namelist('KIM_config_radial_current_test.nml')
        nml_config_path = 'KIM_config_radial_current_test.nml'
        profiles_in_memory = .true.

        call kim_init()
        call set_profiles_from_arrays(r, density, te, ti, q, er, npts)
        call from_kim_factory_get_kim('electrostatic', kim)
        call kim%init()

        j = rg_grid%npts_b / 2
        do sp = 0, plasma%n_species - 1
            if (.not. allocated(plasma%spec(sp)%I03)) then
                error stop 'I03 boundary profile was not allocated'
            end if
            if (.not. allocated(plasma%spec(sp)%I03_cc)) then
                error stop 'I03 cell-centred profile was not allocated'
            end if

            do ell = -2, 2
                call evaluate_susceptibility(plasma%spec(sp)%x1(j), &
                    plasma%spec(sp)%x2(j, ell), moments)
                if (abs(plasma%spec(sp)%I03(j, ell) - moments(0, 3)) > tolerance) then
                    error stop 'retained boundary I03 does not match symbI(0,3)'
                end if

                call evaluate_susceptibility(plasma%spec(sp)%x1_cc(j), &
                    plasma%spec(sp)%x2_cc(j, ell), moments)
                if (abs(plasma%spec(sp)%I03_cc(j, ell) - moments(0, 3)) > tolerance) then
                    error stop 'retained cell-centred I03 does not match symbI(0,3)'
                end if
            end do
        end do

        print *, 'PASS: I03 retained for every configured harmonic'
    end subroutine test_i03_is_retained

    subroutine make_test_profiles(r, density, te, ti, q, er)
        real(dp), intent(out) :: r(:), density(:), te(:), ti(:), q(:), er(:)
        integer :: i
        real(dp) :: fraction

        do i = 1, size(r)
            fraction = real(i - 1, dp) / real(size(r) - 1, dp)
            r(i) = 2.0_dp + 58.0_dp * fraction
            density(i) = 2.0e13_dp * (1.1_dp - r(i) / 100.0_dp)
            te(i) = 1.0e3_dp * (1.2_dp - r(i) / 100.0_dp)
            ti(i) = te(i)
            q(i) = 1.5_dp + fraction
            er(i) = -0.5_dp
        end do
    end subroutine make_test_profiles

    subroutine write_test_namelist(path)
        character(len=*), intent(in) :: path
        integer :: unit

        open(newunit=unit, file=path, status='replace', action='write')
        write(unit, '(A)') '&KIM_CONFIG'
        write(unit, '(A)') ' number_of_ion_species = 1'
        write(unit, '(A)') ' artificial_debye_case = 0'
        write(unit, '(A)') " type_of_run = 'electrostatic'"
        write(unit, '(A)') " collision_model = 'FokkerPlanck'"
        write(unit, '(A)') ' read_species_from_namelist = .false.'
        write(unit, '(A)') " plasma_type = 'D'"
        write(unit, '(A)') '/'
        write(unit, '(A)') '&WKB_DISPERSION'
        write(unit, '(A)') '/'
        write(unit, '(A)') '&KIM_IO'
        write(unit, '(A)') " profile_location = './'"
        write(unit, '(A)') " output_path = './out_radial_current_test/'"
        write(unit, '(A)') ' hdf5_input = .false.'
        write(unit, '(A)') ' hdf5_output = .false.'
        write(unit, '(A)') ' data_verbosity = 0'
        write(unit, '(A)') ' calculate_asymptotics = .false.'
        write(unit, '(A)') '/'
        write(unit, '(A)') '&KIM_SETUP'
        write(unit, '(A)') ' btor = -18000.0'
        write(unit, '(A)') ' R0 = 165.0'
        write(unit, '(A)') ' m_mode = -2'
        write(unit, '(A)') ' n_mode = 1'
        write(unit, '(A)') ' omega = 0.0'
        write(unit, '(A)') ' spline_base = 1'
        write(unit, '(A)') ' type_br_field = 12'
        write(unit, '(A)') ' collisions_off = .false.'
        write(unit, '(A)') ' set_profiles_constant = 0'
        write(unit, '(A)') ' bc_type = 3'
        write(unit, '(A)') ' mphi_max = 2'
        write(unit, '(A)') '/'
        write(unit, '(A)') '&KIM_GRID'
        write(unit, '(A)') " grid_spacing_rg = 'equidistant'"
        write(unit, '(A)') " grid_spacing_xl = 'equidistant'"
        write(unit, '(A)') ' l_space_dim = 32'
        write(unit, '(A)') ' rg_space_dim = 32'
        write(unit, '(A)') " theta_integration = 'GaussLegendre'"
        write(unit, '(A)') ' Larmor_skip_factor = 5'
        write(unit, '(A)') ' gauss_int_nodes_Ntheta = 17'
        write(unit, '(A)') ' gauss_int_nodes_Nx = 31'
        write(unit, '(A)') ' gauss_int_nodes_Nxp = 30'
        write(unit, '(A)') ' r_plas = 50.0'
        write(unit, '(A)') ' r_min = 10.0'
        write(unit, '(A)') ' width_res = 0.5'
        write(unit, '(A)') ' ampl_res = 15.0'
        write(unit, '(A)') ' hrmax_scaling = 1.0'
        write(unit, '(A)') ' kernel_taper_skip_threshold = 1.0e-6'
        write(unit, '(A)') '/'
        write(unit, '(A)') '&KIM_PROFILES'
        write(unit, '(A)') '/'
        close(unit)
    end subroutine write_test_namelist

end program test_radial_current_fourier_kernel
