program test_kim_diagnostics

    use KIM_kinds_m, only: dp
    use kim_diagnostics_m, only: integrate_Ipar
    use kim_qldiff_m, only: calc_dqle22

    implicit none

    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
    real(dp), parameter :: r1 = 1.5_dp, r2 = 4.0_dp

    call test_constant_jpar()
    call test_linear_jpar()
    call test_dqle22_resonant('Dqle22 at resonance, omE/nue = 0.5', 0.5_dp)
    call test_dqle22_resonant('Dqle22 at resonance, omE/nue = 2.0', 2.0_dp)
    call test_em_solve_writes_diagnostics_file()
    call test_em_solve_no_resonance_skips_file()

    print *, 'All kim_diagnostics tests passed.'

contains

    subroutine test_constant_jpar()
        ! jpar = c on [r1, r2]: I_par = 2*pi*int(r*c dr) = pi*c*(r2^2 - r1^2).
        ! Integrand is linear in r, so trapezoid is exact even on a
        ! non-equidistant grid (catches wrong trapezoid weights).
        integer, parameter :: npts = 37
        complex(dp), parameter :: c = (2.0_dp, -3.0_dp)
        real(dp) :: r(npts)
        complex(dp) :: jpar(npts), Ipar, Ipar_exact

        call make_nonequidistant_grid(npts, r)
        jpar = c

        Ipar = integrate_Ipar(npts, r, jpar)
        Ipar_exact = pi * c * (r2**2 - r1**2)

        call assert_close('constant jpar', Ipar, Ipar_exact, 1.0e-10_dp)
    end subroutine

    subroutine test_linear_jpar()
        ! jpar = c*r on [r1, r2]: I_par = 2*pi*c*(r2^3 - r1^3)/3.
        ! Trapezoid is 2nd order; 200 points suffice for 1e-3 relative error.
        integer, parameter :: npts = 200
        complex(dp), parameter :: c = (-1.0_dp, 4.0_dp)
        real(dp) :: r(npts)
        complex(dp) :: jpar(npts), Ipar, Ipar_exact
        integer :: i

        call make_nonequidistant_grid(npts, r)
        do i = 1, npts
            jpar(i) = c * r(i)
        end do

        Ipar = integrate_Ipar(npts, r, jpar)
        Ipar_exact = 2.0_dp * pi * c * (r2**3 - r1**3) / 3.0_dp

        call assert_close('linear jpar', Ipar, Ipar_exact, 1.0e-3_dp)
    end subroutine

    subroutine test_dqle22_resonant(label, omE_over_nue)
        ! At the resonant surface (k_par -> 0, i.e. x1 -> 0) with a pure
        ! radial perturbation field |Br| (Es = 0), the susceptibility-based
        ! D_ql,e22 reduces to the arbitrary-collisionality closed form of
        ! Markl et al 2023 NF 63 126007, Eq. 31:
        !   D_qle22 = (1/8) vTe^2 nue (47 omE^2 + 279 nue^2)
        !             / (omE^4 + 10 omE^2 nue^2 + 9 nue^4) |Br|^2 / B0^2
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: omE_over_nue

        ! Physically plausible CGS magnitudes at a tokamak resonant surface
        real(dp), parameter :: vTe = 1.0e9_dp     ! cm/s
        real(dp), parameter :: nue = 1.0e5_dp     ! 1/s
        real(dp), parameter :: B0 = 1.8e4_dp      ! G
        complex(dp), parameter :: Br = (1.0_dp, 0.0_dp)  ! G
        complex(dp), parameter :: Es = (0.0_dp, 0.0_dp)
        real(dp), parameter :: kpar = 1.0e-10_dp  ! 1/cm -> x1 = 1e-6

        real(dp) :: omE, dqle22, dqle22_exact

        omE = omE_over_nue * nue

        dqle22 = calc_dqle22(vTe, nue, omE, B0, kpar, Es, Br)

        dqle22_exact = 0.125_dp * vTe**2 * nue &
                       * (47.0_dp * omE**2 + 279.0_dp * nue**2) &
                       / (omE**4 + 10.0_dp * omE**2 * nue**2 + 9.0_dp * nue**4) &
                       * abs(Br)**2 / B0**2

        call assert_close(label, cmplx(dqle22, 0.0_dp, dp), &
                          cmplx(dqle22_exact, 0.0_dp, dp), 1.0e-2_dp)
    end subroutine

    subroutine test_em_solve_writes_diagnostics_file()
        ! Integration test: run the full electromagnetic solve on a small
        ! in-memory plasma (same injection path as the QL-Balance adapter)
        ! with write_diagnostics_dat enabled and assert that
        ! kim_diagnostics.dat is written with five finite values and
        ! dqle22 > 0. HDF5 output is off, mirroring the scan driver.
        use config_m, only: profiles_in_memory, nml_config_path
        use species_m, only: set_profiles_from_arrays
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101
        character(len=*), parameter :: out_file = &
            './out_diag_test/m-2_n1/kim_diagnostics.dat'

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        real(dp) :: vals(5)
        class(kim_t), allocatable :: kim_instance
        character(len=256) :: header
        integer :: i, iunit, ios
        logical :: ex

        ! Analytic profiles on r in [2, 60] cm; q crosses |m/n| = 2 at
        ! r ~ 31 cm, inside the solver domain [r_min, r_plas] = [10, 50].
        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_diag_test.nml', -2)
        nml_config_path = './KIM_config_diag_test.nml'

        ! Remove stale output so an old file cannot mask a failure
        inquire(file=out_file, exist=ex)
        if (ex) then
            open(newunit=iunit, file=out_file, status='old')
            close(iunit, status='delete')
        end if

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electromagnetic', kim_instance)
        call kim_instance%init()
        call kim_instance%run()

        inquire(file=out_file, exist=ex)
        if (.not. ex) then
            print *, 'FAIL: EM solve did not write ', out_file
            error stop
        end if

        open(newunit=iunit, file=out_file, status='old', action='read')
        read(iunit, '(A)') header
        read(iunit, *, iostat=ios) vals
        close(iunit)

        if (header(1:1) /= '#') then
            print *, 'FAIL: kim_diagnostics.dat missing # header line'
            error stop
        end if
        if (ios /= 0) then
            print *, 'FAIL: could not parse 5 values from kim_diagnostics.dat'
            error stop
        end if
        do i = 1, 5
            if (vals(i) /= vals(i) .or. abs(vals(i)) > huge(1.0_dp)) then
                print *, 'FAIL: non-finite diagnostics value at column ', i
                print *, '  values = ', vals
                error stop
            end if
        end do
        if (vals(1) <= 0.0_dp) then
            print *, 'FAIL: dqle22 must be positive, got ', vals(1)
            error stop
        end if

        print *, 'PASS: EM solve wrote kim_diagnostics.dat'
        print *, '  dqle22 = ', vals(1)
        print *, '  Ipar   = ', vals(2), vals(3)
        print *, '  Ipar_e = ', vals(4), vals(5)
    end subroutine

    subroutine test_em_solve_no_resonance_skips_file()
        ! Guard test: with m chosen so that |m/n| = 4 lies outside the
        ! q range [1.5, 2.5] of the analytic test profiles, no resonant
        ! surface exists. The diagnostics must then NOT be written --
        ! a missing kim_diagnostics.dat is the scan driver's failure
        ! signal, while plausible-looking values evaluated at a bogus
        ! radius would be silently machine-read.
        use config_m, only: profiles_in_memory, nml_config_path
        use species_m, only: set_profiles_from_arrays
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim
        use fields_m, only: EBdat, EBdat_t

        integer, parameter :: npts = 101
        character(len=*), parameter :: out_file = &
            './out_diag_test/m-4_n1/kim_diagnostics.dat'

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        integer :: iunit
        logical :: ex

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_diag_test_nores.nml', -4)
        nml_config_path = './KIM_config_diag_test_nores.nml'

        ! Remove stale output so an old file cannot mask a failure
        inquire(file=out_file, exist=ex)
        if (ex) then
            open(newunit=iunit, file=out_file, status='old')
            close(iunit, status='delete')
        end if

        profiles_in_memory = .true.
        call kim_init()
        ! Resets the lazy resonance-detection flag, so the stale r_res
        ! from the previous test cannot leak into this run.
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)
        ! The EM solver allocates EBdat components without re-entry
        ! handling; reset them for the second solve in this process.
        EBdat = EBdat_t()

        call from_kim_factory_get_kim('electromagnetic', kim_instance)
        call kim_instance%init()
        call kim_instance%run()

        inquire(file=out_file, exist=ex)
        if (ex) then
            print *, 'FAIL: EM solve without resonance wrote ', out_file
            print *, '  (values would be evaluated at a bogus radius)'
            error stop
        end if

        print *, 'PASS: no resonance -> kim_diagnostics.dat not written'
    end subroutine

    subroutine make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                  q_prof, Er_prof)
        ! Analytic profiles on r in [2, 60] cm with q in [1.5, 2.5];
        ! q crosses |m/n| = 2 at r ~ 31 cm, inside the solver domain
        ! [r_min, r_plas] = [10, 50].
        integer, intent(in) :: npts
        real(dp), intent(out) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp), intent(out) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        integer :: i

        do i = 1, npts
            r_prof(i) = 2.0_dp + 58.0_dp * real(i - 1, dp) / real(npts - 1, dp)
            n_prof(i) = 2.0e13_dp * (1.1_dp - r_prof(i) / 100.0_dp)
            Te_prof(i) = 1.0e3_dp * (1.2_dp - r_prof(i) / 100.0_dp)
            Ti_prof(i) = Te_prof(i)
            q_prof(i) = 1.5_dp + (r_prof(i) - 2.0_dp) / 58.0_dp
            Er_prof(i) = -0.5_dp
        end do
    end subroutine

    subroutine write_test_namelist(path, m_mode_val)
        ! Minimal electromagnetic configuration on a coarse grid, written
        ! in the exact group order kim_read_config reads them.
        character(len=*), intent(in) :: path
        integer, intent(in) :: m_mode_val
        integer :: iunit

        open(newunit=iunit, file=path, status='replace', action='write')
        write(iunit, '(A)') '&KIM_CONFIG'
        write(iunit, '(A)') ' number_of_ion_species = 1'
        write(iunit, '(A)') ' artificial_debye_case = 0'
        write(iunit, '(A)') " type_of_run = 'electromagnetic'"
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
        write(iunit, '(A)') " output_path = './out_diag_test/'"
        write(iunit, '(A)') ' hdf5_input = .false.'
        write(iunit, '(A)') ' hdf5_output = .false.'
        write(iunit, '(A)') ' log_level = 3'
        write(iunit, '(A)') ' data_verbosity = 0'
        write(iunit, '(A)') ' calculate_asymptotics = .false.'
        write(iunit, '(A)') " h5_out_file = ''"
        write(iunit, '(A)') ' write_diagnostics_dat = .true.'
        write(iunit, '(A)') '/'
        write(iunit, '(A)') '&KIM_SETUP'
        write(iunit, '(A)') ' btor = -18000.0'
        write(iunit, '(A)') ' R0 = 165.0'
        write(iunit, '(A,I0)') ' m_mode = ', m_mode_val
        write(iunit, '(A)') ' n_mode = 1'
        write(iunit, '(A)') ' omega = 0.0'
        write(iunit, '(A)') ' spline_base = 1'
        write(iunit, '(A)') ' type_br_field = 12'
        write(iunit, '(A)') ' collisions_off = .false.'
        write(iunit, '(A)') ' set_profiles_constant = 0'
        write(iunit, '(A)') ' bc_type = 3'
        write(iunit, '(A)') ' mphi_max = 0'
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
        write(iunit, '(A)') '/'
        close(iunit)
    end subroutine

    subroutine make_nonequidistant_grid(npts, r)
        integer, intent(in) :: npts
        real(dp), intent(out) :: r(npts)
        real(dp) :: s
        integer :: i

        do i = 1, npts
            s = real(i - 1, dp) / real(npts - 1, dp)
            r(i) = r1 + (r2 - r1) * s**2
        end do
    end subroutine

    subroutine assert_close(label, actual, expected, rel_tol)
        character(len=*), intent(in) :: label
        complex(dp), intent(in) :: actual, expected
        real(dp), intent(in) :: rel_tol
        real(dp) :: err_re, err_im, denom

        ! Normalize by the full complex magnitude so that zero real or
        ! imaginary parts of the expected value cannot cause division
        ! by zero; tiny() guards the all-zero case.
        denom = max(abs(expected), tiny(1.0_dp))
        err_re = abs(real(actual, dp) - real(expected, dp)) / denom
        err_im = abs(aimag(actual) - aimag(expected)) / denom

        if (err_re > rel_tol .or. err_im > rel_tol) then
            print *, 'FAIL: ', label
            print *, '  actual   = ', actual
            print *, '  expected = ', expected
            print *, '  rel err (re, im) = ', err_re, err_im
            error stop
        end if

        print *, 'PASS: ', label, ' rel err (re, im) = ', err_re, err_im
    end subroutine

end program
