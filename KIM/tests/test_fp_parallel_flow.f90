program test_fp_parallel_flow
    use KIM_kinds_m, only: dp
    use periodic_background_m, only: build_periodic_plasma, &
                                     sample_periodic_species_profiles
    use profile_input_m, only: read_parallel_flow_profile, FLOW_PROFILE_OK, &
                               FLOW_PROFILE_AMBIGUOUS_CONVENTION, &
                               FLOW_PROFILE_INCOMPATIBLE_COVERAGE

    implicit none

    call test_global_and_periodic_flow()
    call test_file_input_contract()

    print *, 'All tests PASSED'
    stop 0

contains

    subroutine test_global_and_periodic_flow()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use config_m, only: profiles_in_memory, nml_config_path
        use species_m, only: plasma, set_profiles_from_arrays
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim
        use kim_resonances_m, only: r_res
        use grid_m, only: rg_grid
        use setup_m, only: omega

        integer, parameter :: npts = 201
        real(dp) :: r(npts), ne(npts), Te(npts), Ti(npts), q(npts), Er(npts)
        real(dp) :: Vpar(npts), rm, rho_ref, dx_asis, dx_tr, period_length
        real(dp) :: r_asis, Vpar_true, Vpar_built, expected, delta_expected
        real(dp), allocatable :: x2_zero_e(:, :), x2_zero_i(:, :)
        complex(dp), allocatable :: I00_zero_i(:, :)
        real(dp) :: n_per(0:1), T_per(0:1), Vpar_per(0:1)
        real(dp) :: n_shift(0:1), T_shift(0:1), Vpar_shift(0:1)
        real(dp) :: q_per, Er_per, q_shift, Er_shift
        class(kim_t), allocatable :: kim_instance
        integer :: j, mphi, n_rg

        call make_profiles(r, ne, Te, Ti, q, Er)
        call write_test_namelist('./KIM_config_fp_parallel_flow.nml')
        nml_config_path = './KIM_config_fp_parallel_flow.nml'
        profiles_in_memory = .true.
        call kim_init()

        ! Baseline: an explicitly zero field-parallel ion flow must retain the
        ! exact historical arithmetic for x2 and its susceptibility tables.
        Vpar = 0.0_dp
        call set_profiles_from_arrays(r, ne, Te, Ti, q, Er, npts, Vpar_in=Vpar)
        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        allocate(x2_zero_e(size(plasma%spec(0)%x2, 1), -1:1))
        allocate(x2_zero_i(size(plasma%spec(1)%x2, 1), -1:1))
        allocate(I00_zero_i(size(plasma%spec(1)%I00, 1), -1:1))
        x2_zero_e = plasma%spec(0)%x2
        x2_zero_i = plasma%spec(1)%x2
        I00_zero_i = plasma%spec(1)%I00

        do j = 1, rg_grid%npts_b
            do mphi = -1, 1
                expected = -(plasma%om_E(j) + &
                    real(mphi, dp)*plasma%spec(1)%omega_c(j) - omega) / &
                    plasma%spec(1)%nu(j)
                call assert_bitwise_equal('zero-flow x2', &
                                          plasma%spec(1)%x2(j, mphi), expected)
            end do
        end do
        print *, 'PASS: zero flow is bitwise equivalent to the historical x2 formula'

        ! Finite constant field-parallel flow: all ion harmonics receive the
        ! same -kp*Vpar/nu shift, while the non-drifting electrons are unchanged.
        Vpar = 2.5e6_dp
        call set_profiles_from_arrays(r, ne, Te, Ti, q, Er, npts, Vpar_in=Vpar)
        call kim_instance%init()

        do j = 1, rg_grid%npts_b
            if (abs(plasma%spec(0)%Vpar(j)) > tiny(1.0_dp)) then
                print *, 'FAIL: electron Maxwellian was enabled for ion flow'
                error stop
            end if
            call assert_relative_close('constant ion Vpar', plasma%spec(1)%Vpar(j), &
                                       Vpar(1), 10.0_dp*epsilon(1.0_dp))
            do mphi = -1, 1
                expected = -(plasma%om_E(j) + plasma%kp(j)*plasma%spec(1)%Vpar(j) + &
                    real(mphi, dp)*plasma%spec(1)%omega_c(j) - omega) / &
                    plasma%spec(1)%nu(j)
                call assert_relative_close('finite-flow x2', &
                                           plasma%spec(1)%x2(j, mphi), expected, 5.0e-13_dp)
                delta_expected = -plasma%kp(j)*plasma%spec(1)%Vpar(j) / &
                                 plasma%spec(1)%nu(j)
                call assert_relative_close('finite-flow x2 shift', &
                    plasma%spec(1)%x2(j, mphi) - x2_zero_i(j, mphi), &
                    delta_expected, 1.0e-10_dp)
                call assert_bitwise_equal('stationary-electron x2', &
                                          plasma%spec(0)%x2(j, mphi), &
                                          x2_zero_e(j, mphi))
            end do
        end do
        if (maxval(abs(plasma%spec(1)%I00 - I00_zero_i)) <= 100.0_dp*epsilon(1.0_dp)) then
            print *, 'FAIL: finite flow did not change ion susceptibility'
            error stop
        end if
        print *, 'PASS: finite flow shifts every ion harmonic with the correct sign'

        ! Preserve the same field-parallel profile through forced periodization.
        call prepare_resonances
        rm = r_res
        if (.not. (rm > 0.0_dp)) error stop 'parallel-flow test resonance not found'
        j = minloc(abs(plasma%r_grid - rm), dim=1)
        rho_ref = plasma%spec(1)%rho_L(j)
        dx_asis = 5.0_dp*rho_ref
        dx_tr = 10.0_dp*rho_ref
        n_rg = 96
        r_asis = rm + 0.5_dp*dx_asis
        Vpar_true = interp_lagrange4(plasma%r_grid, plasma%spec(1)%Vpar, &
                                     plasma%grid_size, r_asis)

        call build_periodic_plasma(rm, dx_asis, dx_tr, n_rg)
        Vpar_built = interp_lagrange4(rg_grid%xb, plasma%spec(1)%Vpar, &
                                      plasma%grid_size, r_asis)
        call assert_relative_close('periodic as-is Vpar', Vpar_built, Vpar_true, 1.0e-10_dp)

        period_length = 2.0_dp*(dx_asis + dx_tr)
        call sample_periodic_species_profiles(rm, dx_asis, dx_tr, &
            rm - 0.5_dp*period_length, n_per, T_per, q_per, Er_per, Vpar_per)
        call sample_periodic_species_profiles(rm, dx_asis, dx_tr, &
            rm + 0.5_dp*period_length, n_shift, T_shift, q_shift, Er_shift, Vpar_shift)
        call assert_relative_close('periodic seam Vpar', Vpar_shift(1), &
                                   Vpar_per(1), 1.0e-12_dp)

        do j = 1, rg_grid%npts_b
            do mphi = -1, 1
                expected = -(plasma%om_E(j) + plasma%kp(j)*plasma%spec(1)%Vpar(j) + &
                    real(mphi, dp)*plasma%spec(1)%omega_c(j) - omega) / &
                    plasma%spec(1)%nu(j)
                if (.not. ieee_is_finite(plasma%spec(1)%x2(j, mphi))) then
                    print *, 'FAIL: non-finite periodic x2 at ', j, mphi
                    error stop
                end if
                call assert_relative_close('periodic finite-flow x2', &
                                           plasma%spec(1)%x2(j, mphi), expected, 5.0e-13_dp)
            end do
        end do
        print *, 'PASS: finite flow survives periodization and susceptibility rebuild'
    end subroutine test_global_and_periodic_flow

    subroutine test_file_input_contract()
        use config_m, only: profile_location, Vz_file, parallel_flow_convention
        use setup_m, only: btor, R0

        real(dp) :: r(5), q(5), velocity(5), Vpar(5), expected_hz
        integer :: info, i, iunit

        call execute_command_line('mkdir -p ./flow_profile_contract')
        profile_location = './flow_profile_contract'
        Vz_file = 'Vz.dat'
        r = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp]
        q = 3.0_dp
        velocity = 1.0e6_dp
        call write_flow_file(trim(profile_location)//'/'//trim(Vz_file), r, velocity)

        parallel_flow_convention = 'none'
        call read_parallel_flow_profile(r, q, Vpar, info)
        if (info /= FLOW_PROFILE_AMBIGUOUS_CONVENTION) then
            print *, 'FAIL: undeclared finite velocity convention status = ', info
            error stop
        end if

        parallel_flow_convention = 'toroidal'
        call read_parallel_flow_profile(r, q, Vpar, info)
        if (info /= FLOW_PROFILE_OK) then
            print *, 'FAIL: valid toroidal velocity profile status = ', info
            error stop
        end if
        do i = 1, size(r)
            expected_hz = sign(1.0_dp, btor) / &
                          sqrt(1.0_dp + (r(i)/(q(i)*R0))**2)
            call assert_relative_close('toroidal-to-parallel projection', Vpar(i), &
                                       velocity(i)/expected_hz, 5.0e-15_dp)
        end do

        call write_flow_file(trim(profile_location)//'/'//trim(Vz_file), &
                             r(2:4), velocity(2:4))
        call read_parallel_flow_profile(r, q, Vpar, info)
        if (info /= FLOW_PROFILE_INCOMPATIBLE_COVERAGE) then
            print *, 'FAIL: short velocity profile coverage status = ', info
            error stop
        end if

        open(newunit=iunit, file=trim(profile_location)//'/'//trim(Vz_file), &
             status='old')
        close(iunit, status='delete')
        print *, 'PASS: velocity convention and radial coverage are validated'
    end subroutine test_file_input_contract

    subroutine make_profiles(r, ne, Te, Ti, q, Er)
        real(dp), intent(out) :: r(:), ne(:), Te(:), Ti(:), q(:), Er(:)
        integer :: i

        do i = 1, size(r)
            r(i) = 2.0_dp + 58.0_dp*real(i - 1, dp)/real(size(r) - 1, dp)
            ne(i) = 2.0e13_dp*(1.1_dp - r(i)/100.0_dp)
            Te(i) = 1.0e3_dp*(1.2_dp - r(i)/100.0_dp)
            Ti(i) = 0.8e3_dp*(1.15_dp - r(i)/120.0_dp)
            q(i) = 1.5_dp + 2.5_dp*(r(i) - 2.0_dp)/58.0_dp
            Er(i) = -0.5_dp
        end do
    end subroutine make_profiles

    subroutine write_flow_file(path, r, velocity)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: r(:), velocity(:)
        integer :: i, iunit

        open(newunit=iunit, file=path, status='replace', action='write')
        do i = 1, size(r)
            write(iunit, '(2(ES24.16,1X))') r(i), velocity(i)
        end do
        close(iunit)
    end subroutine write_flow_file

    subroutine assert_bitwise_equal(label, actual, expected)
        use, intrinsic :: iso_fortran_env, only: int64

        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual, expected

        if (transfer(actual, 0_int64) /= transfer(expected, 0_int64)) then
            print *, 'FAIL: ', trim(label), ' is not bitwise equal'
            print *, '  actual = ', actual, ' expected = ', expected
            error stop
        end if
    end subroutine assert_bitwise_equal

    subroutine assert_relative_close(label, actual, expected, tolerance)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual, expected, tolerance
        real(dp) :: relative_error

        relative_error = abs(actual - expected) / max(abs(expected), tiny(1.0_dp))
        if (relative_error > tolerance) then
            print *, 'FAIL: ', trim(label)
            print *, '  actual = ', actual, ' expected = ', expected
            print *, '  relative error = ', relative_error, ' tolerance = ', tolerance
            error stop
        end if
    end subroutine assert_relative_close

    real(dp) function interp_lagrange4(grid, vals, ngrid, x) result(fx)
        integer, intent(in) :: ngrid
        real(dp), intent(in) :: grid(ngrid), vals(ngrid), x
        integer, parameter :: nlagr = 4, nder = 0
        real(dp) :: coef(0:nder, nlagr)
        integer :: ir, ibeg, iend

        call binsrc(grid, 1, ngrid, x, ir)
        ibeg = max(1, ir - nlagr/2)
        iend = ibeg + nlagr - 1
        if (iend > ngrid) then
            iend = ngrid
            ibeg = iend - nlagr + 1
        end if
        call plag_coeff(nlagr, nder, x, grid(ibeg:iend), coef)
        fx = sum(coef(0, :)*vals(ibeg:iend))
    end function interp_lagrange4

    subroutine write_test_namelist(path)
        character(len=*), intent(in) :: path
        integer :: iunit

        open(newunit=iunit, file=path, status='replace', action='write')
        write(iunit, '(A)') '&KIM_CONFIG'
        write(iunit, '(A)') ' number_of_ion_species = 1'
        write(iunit, '(A)') ' artificial_debye_case = 0'
        write(iunit, '(A)') " type_of_run = 'electrostatic'"
        write(iunit, '(A)') " collision_model = 'FokkerPlanck'"
        write(iunit, '(A)') ' read_species_from_namelist = .false.'
        write(iunit, '(A)') ' turn_off_ions = .false.'
        write(iunit, '(A)') ' turn_off_electrons = .false.'
        write(iunit, '(A)') " plasma_type = 'D'"
        write(iunit, '(A)') ' rescale_density = .false.'
        write(iunit, '(A)') ' number_density_rescale = 1.0'
        write(iunit, '(A)') ' ion_flr_scale_factor = 1.0'
        write(iunit, '(A)') '/'
        write(iunit, '(A)') '&WKB_DISPERSION'
        write(iunit, '(A)') '/'
        write(iunit, '(A)') '&KIM_IO'
        write(iunit, '(A)') " profile_location = './'"
        write(iunit, '(A)') " output_path = './out_fp_parallel_flow/'"
        write(iunit, '(A)') ' hdf5_input = .false.'
        write(iunit, '(A)') ' hdf5_output = .false.'
        write(iunit, '(A)') ' log_level = 3'
        write(iunit, '(A)') ' data_verbosity = 0'
        write(iunit, '(A)') ' calculate_asymptotics = .false.'
        write(iunit, '(A)') " h5_out_file = ''"
        write(iunit, '(A)') ' write_diagnostics_dat = .false.'
        write(iunit, '(A)') '/'
        write(iunit, '(A)') '&KIM_SETUP'
        write(iunit, '(A)') ' btor = -18000.0'
        write(iunit, '(A)') ' R0 = 165.0'
        write(iunit, '(A)') ' m_mode = -6'
        write(iunit, '(A)') ' n_mode = 2'
        write(iunit, '(A)') ' omega = 0.0'
        write(iunit, '(A)') ' spline_base = 1'
        write(iunit, '(A)') ' type_br_field = 12'
        write(iunit, '(A)') ' collisions_off = .false.'
        write(iunit, '(A)') ' set_profiles_constant = 0'
        write(iunit, '(A)') ' bc_type = 3'
        write(iunit, '(A)') ' mphi_max = 1'
        write(iunit, '(A)') ' Br_boundary_re = 1.0'
        write(iunit, '(A)') ' Br_boundary_im = 0.0'
        write(iunit, '(A)') '/'
        write(iunit, '(A)') '&KIM_GRID'
        write(iunit, '(A)') " grid_spacing_rg = 'equidistant'"
        write(iunit, '(A)') " grid_spacing_xl = 'equidistant'"
        write(iunit, '(A)') ' l_space_dim = 32'
        write(iunit, '(A)') ' rg_space_dim = 32'
        write(iunit, '(A)') " theta_integration = 'GaussLegendre'"
        write(iunit, '(A)') " theta_integration_method = ''"
        write(iunit, '(A)') ' Larmor_skip_factor = 5'
        write(iunit, '(A)') ' gauss_int_nodes_Ntheta = 17'
        write(iunit, '(A)') ' gauss_int_nodes_Nx = 31'
        write(iunit, '(A)') ' gauss_int_nodes_Nxp = 30'
        write(iunit, '(A)') ' r_plas = 50.0'
        write(iunit, '(A)') ' r_min = 10.0'
        write(iunit, '(A)') ' width_res = 0.5'
        write(iunit, '(A)') ' ampl_res = 15.0'
        write(iunit, '(A)') ' hrmax_scaling = 1.0'
        write(iunit, '(A)') ' rkf45_atol = 1.0e-9'
        write(iunit, '(A)') ' rkf45_rtol = 1.0e-6'
        write(iunit, '(A)') ' kernel_taper_skip_threshold = 1.0e-6'
        write(iunit, '(A)') " quadpack_algorithm = 'QAG'"
        write(iunit, '(A)') ' quadpack_key = 6'
        write(iunit, '(A)') ' quadpack_limit = 500'
        write(iunit, '(A)') ' quadpack_epsabs = 1.0e-10'
        write(iunit, '(A)') ' quadpack_epsrel = 1.0e-10'
        write(iunit, '(A)') ' quadpack_use_u_substitution = .true.'
        write(iunit, '(A)') '/'
        write(iunit, '(A)') '&KIM_PROFILES'
        write(iunit, '(A)') " parallel_flow_convention = 'none'"
        write(iunit, '(A)') '/'
        close(iunit)
    end subroutine write_test_namelist

end program test_fp_parallel_flow
