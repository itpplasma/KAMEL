program test_radial_current_fourier_kernel
    use KIM_kinds_m, only: dp
    use radial_current_fourier_kernel_m, only: radial_flr_coefficients, &
        scaled_bessel_harmonic

    implicit none

    call test_i03_is_retained()
    call test_scaled_bessel_harmonics()
    call test_radial_flr_coefficients()
    call test_zero_wavenumber_radial_flr_coefficients()

    print *, 'All radial-current Fourier-kernel tests PASSED'
    stop 0

contains

    subroutine test_radial_flr_coefficients()
        real(dp), parameter :: derivative_step = 2.0e-5_dp
        real(dp), parameter :: tolerance = 2.0e-7_dp
        real(dp), parameter :: ks = 0.71_dp, kr = -0.37_dp
        real(dp), parameter :: ksp = -0.46_dp, krp = 0.29_dp
        real(dp), parameter :: rho_l = 0.83_dp
        real(dp), parameter :: a1 = -0.24_dp, a2 = 0.17_dp
        complex(dp) :: o0, o2, w0, w2
        complex(dp) :: expected_o0, expected_o2, expected_w0, expected_w2
        complex(dp) :: phase
        integer :: ell

        do ell = -3, 3
            call radial_flr_coefficients(ell, ks, kr, ksp, krp, rho_l, a1, a2, &
                o0, o2, w0, w2)

            phase = harmonic_phase(ell, ks, kr, ksp, krp)
            expected_o0 = centred_response_derivative(0, ell, ks, kr, ksp, krp, &
                rho_l, a1, a2, derivative_step) / phase
            expected_o2 = centred_response_derivative(2, ell, ks, kr, ksp, krp, &
                rho_l, a1, a2, derivative_step) / phase
            expected_w0 = centred_mixed_derivative(0, ell, ks, kr, ksp, krp, &
                rho_l, a1, a2, derivative_step) / phase
            expected_w2 = centred_mixed_derivative(2, ell, ks, kr, ksp, krp, &
                rho_l, a1, a2, derivative_step) / phase

            call assert_complex_close(o0, expected_o0, tolerance, &
                'O0 disagrees with independent response derivative')
            call assert_complex_close(o2, expected_o2, tolerance, &
                'O2 disagrees with independent response derivative')
            call assert_complex_close(w0, expected_w0, tolerance, &
                'W0 disagrees with independent mixed derivative')
            call assert_complex_close(w2, expected_w2, tolerance, &
                'W2 disagrees with independent mixed derivative')
        end do

        print *, 'PASS: radial FLR coefficients match centred derivatives'

    end subroutine test_radial_flr_coefficients

    complex(dp) function centred_response_derivative(moment, ell, ks, kr, ksp, krp, &
        rho_l, a1, a2, step) result(value)
        integer, intent(in) :: moment, ell
        real(dp), intent(in) :: ks, kr, ksp, krp, rho_l, a1, a2, step

        value = (ordinary_flr_coefficient(moment, ell, ks, kr, ksp + step, krp, &
            rho_l, a1, a2) - ordinary_flr_coefficient(moment, ell, ks, kr, &
            ksp - step, krp, rho_l, a1, a2)) / (2.0_dp * step)
    end function centred_response_derivative

    complex(dp) function centred_mixed_derivative(moment, ell, ks, kr, ksp, krp, &
        rho_l, a1, a2, step) result(value)
        integer, intent(in) :: moment, ell
        real(dp), intent(in) :: ks, kr, ksp, krp, rho_l, a1, a2, step

        value = (ordinary_flr_coefficient(moment, ell, ks + step, kr, ksp + step, &
            krp, rho_l, a1, a2) - ordinary_flr_coefficient(moment, ell, ks + step, &
            kr, ksp - step, krp, rho_l, a1, a2) &
            - ordinary_flr_coefficient(moment, ell, ks - step, kr, ksp + step, krp, &
            rho_l, a1, a2) + ordinary_flr_coefficient(moment, ell, ks - step, kr, &
            ksp - step, krp, rho_l, a1, a2)) / (4.0_dp * step**2)
    end function centred_mixed_derivative

    subroutine test_zero_wavenumber_radial_flr_coefficients()
        real(dp), parameter :: tolerance = 1.0e-14_dp
        real(dp), parameter :: derivative_tolerance = 2.0e-6_dp
        real(dp), parameter :: derivative_step = 2.0e-4_dp
        real(dp), parameter :: rho_l = 0.64_dp
        real(dp), parameter :: a1 = -0.31_dp, a2 = 0.23_dp
        complex(dp) :: o0, o2, w0, w2
        complex(dp) :: phase
        real(dp) :: expected_w0, expected_w2
        integer :: ell

        expected_w0 = rho_l**2 * (0.5_dp * a1 + a2)
        expected_w2 = 0.25_dp * rho_l**2 * a2

        do ell = -3, 3
            call radial_flr_coefficients(ell, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                rho_l, a1, a2, o0, o2, w0, w2)

            call assert_complex_close(o0, (0.0_dp, 0.0_dp), tolerance, &
                'single radial insertion must vanish at zero wavenumber')
            call assert_complex_close(o2, (0.0_dp, 0.0_dp), tolerance, &
                'single radial insertion must vanish at zero wavenumber')
            if (abs(ell) == 1) then
                call assert_complex_close(w0, cmplx(expected_w0, 0.0_dp, dp), tolerance, &
                    'W0 zero-wavenumber anchor must be carried by ell=+/-1')
                call assert_complex_close(w2, cmplx(expected_w2, 0.0_dp, dp), tolerance, &
                    'W2 zero-wavenumber anchor must be carried by ell=+/-1')
            else
                call assert_complex_close(w0, (0.0_dp, 0.0_dp), tolerance, &
                    'W0 zero-wavenumber anchor leaked to the wrong harmonic')
                call assert_complex_close(w2, (0.0_dp, 0.0_dp), tolerance, &
                    'W2 zero-wavenumber anchor leaked to the wrong harmonic')
            end if
        end do

        do ell = -2, 2
            call radial_flr_coefficients(ell, 0.73_dp, -0.28_dp, 0.0_dp, 0.0_dp, &
                rho_l, a1, a2, o0, o2, w0, w2)
            phase = harmonic_phase(ell, 0.73_dp, -0.28_dp, 0.0_dp, 0.0_dp)
            call assert_complex_close(o0, centred_response_derivative(0, ell, &
                0.73_dp, -0.28_dp, 0.0_dp, 0.0_dp, rho_l, a1, a2, &
                derivative_step) / phase, derivative_tolerance, &
                'O0 response-zero adjacent-order limit failed')
            call assert_complex_close(w0, centred_mixed_derivative(0, ell, &
                0.73_dp, -0.28_dp, 0.0_dp, 0.0_dp, rho_l, a1, a2, &
                derivative_step) / phase, derivative_tolerance, &
                'W0 response-zero adjacent-order limit failed')

            call radial_flr_coefficients(ell, 0.0_dp, 0.0_dp, -0.61_dp, 0.34_dp, &
                rho_l, a1, a2, o0, o2, w0, w2)
            phase = harmonic_phase(ell, 0.0_dp, 0.0_dp, -0.61_dp, 0.34_dp)
            call assert_complex_close(o0, centred_response_derivative(0, ell, &
                0.0_dp, 0.0_dp, -0.61_dp, 0.34_dp, rho_l, a1, a2, &
                derivative_step) / phase, derivative_tolerance, &
                'O0 source-zero adjacent-order limit failed')
            call assert_complex_close(w0, centred_mixed_derivative(0, ell, &
                0.0_dp, 0.0_dp, -0.61_dp, 0.34_dp, rho_l, a1, a2, &
                derivative_step) / phase, derivative_tolerance, &
                'W0 source-zero adjacent-order limit failed')
        end do

        print *, 'PASS: zero-wavenumber anchors require ell=+/-1'
    end subroutine test_zero_wavenumber_radial_flr_coefficients

    complex(dp) function ordinary_flr_coefficient(moment, ell, ks, kr, ksp, krp, &
        rho_l, a1, a2) result(value)
        use fortnum_special, only: bessel_in

        integer, intent(in) :: moment, ell
        real(dp), intent(in) :: ks, kr, ksp, krp, rho_l, a1, a2
        real(dp) :: bplus, bcross, gamma, lambda, kperp, kperp_p

        kperp = hypot(ks, kr)
        kperp_p = hypot(ksp, krp)
        bplus = 0.5_dp * rho_l**2 * (kperp**2 + kperp_p**2)
        bcross = rho_l**2 * kperp * kperp_p
        gamma = exp(-bplus) * bessel_in(ell, bcross)
        lambda = bcross * exp(-bplus) * bessel_in(ell - 1, bcross)

        select case (moment)
        case (0)
            value = gamma * (a1 + (1.0_dp - bplus - real(ell, dp)) * a2) &
                + lambda * a2
        case (2)
            value = 0.5_dp * a2 * gamma
        case default
            error stop 'unsupported ordinary FLR moment in test'
        end select
        value = value * harmonic_phase(ell, ks, kr, ksp, krp)
    end function ordinary_flr_coefficient

    complex(dp) function harmonic_phase(ell, ks, kr, ksp, krp) result(value)
        integer, intent(in) :: ell
        real(dp), intent(in) :: ks, kr, ksp, krp
        real(dp) :: alpha, alpha_p

        alpha = atan2(kr, ks)
        alpha_p = atan2(krp, ksp)
        value = exp(cmplx(0.0_dp, real(ell, dp) * (alpha_p - alpha), dp))
    end function harmonic_phase

    subroutine assert_complex_close(actual, expected, tolerance, message)
        complex(dp), intent(in) :: actual, expected
        real(dp), intent(in) :: tolerance
        character(len=*), intent(in) :: message

        if (abs(actual - expected) > tolerance * (1.0_dp + abs(expected))) then
            print *, 'actual:  ', actual
            print *, 'expected:', expected
            error stop message
        end if
    end subroutine assert_complex_close

    subroutine test_scaled_bessel_harmonics()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use fortnum_special, only: bessel_in

        real(dp), parameter :: tolerance = 1.0e-13_dp
        real(dp) :: actual, expected
        integer :: ell

        do ell = -4, 4
            actual = scaled_bessel_harmonic(ell, 1.25_dp, 0.75_dp)
            expected = exp(-1.25_dp) * bessel_in(ell, 0.75_dp)
            if (abs(actual - expected) > tolerance * (1.0_dp + abs(expected))) then
                error stop 'scaled harmonic Bessel mismatch at moderate argument'
            end if

            if (scaled_bessel_harmonic(ell, 1.25_dp, 0.75_dp) /= &
                scaled_bessel_harmonic(-ell, 1.25_dp, 0.75_dp)) then
                error stop 'integer-order scaled Bessel symmetry failed'
            end if
        end do

        if (scaled_bessel_harmonic(0, 0.0_dp, 0.0_dp) /= 1.0_dp) then
            error stop 'scaled I0 has wrong zero-argument limit'
        end if
        do ell = 1, 4
            if (scaled_bessel_harmonic(ell, 0.0_dp, 0.0_dp) /= 0.0_dp) then
                error stop 'nonzero-order scaled Bessel has wrong zero-argument limit'
            end if
        end do

        do ell = -4, 4
            actual = scaled_bessel_harmonic(ell, 1000.5_dp, 1000.0_dp)
            if (.not. ieee_is_finite(actual)) then
                error stop 'scaled harmonic Bessel is not finite at large argument'
            end if
            if (actual < 0.0_dp) then
                error stop 'scaled harmonic Bessel must be nonnegative for positive argument'
            end if
        end do

        print *, 'PASS: arbitrary-order scaled Bessel values are finite and symmetric'
    end subroutine test_scaled_bessel_harmonics

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
