program test_periodic_ion_tensor
    use KIM_kinds_m, only: dp
    use constants_m, only: pi, sol
    use config_m, only: resolved_ion_ifunc_conservation_model, periodic_Bparallel_ratio
    use setup_m, only: mphi_max
    use species_m, only: evaluate_susceptibility
    use quasilinear_flr_m, only: calc_ion_flr_harmonic
    use rt_electrostatic_periodic_m, only: compute_periodic_ion_tensor, &
        compute_periodic_ion_tensor_spectrum, electrostatic_periodic_t
    use fortnum_special, only: bessel_in
    implicit none

    real(dp), parameter :: ks = 0.7_dp, kpar = 0.3_dp, vti = 1.0_dp
    real(dp), parameter :: nui = 0.4_dp, omega_ci = 1.0_dp, om_e = 0.7_dp
    real(dp), parameter :: b0 = 1.0_dp, period = pi, radius = 0.37_dp
    complex(dp) :: moments(0:3,0:3)
    integer :: failures
    character(len=32) :: argument
    type(electrostatic_periodic_t) :: solver

    call get_command_argument(1, argument)
    if (trim(argument) == 'reject_bparallel') then
        periodic_Bparallel_ratio = (1.0_dp, 0.0_dp)
        call solver%run()
        stop 0
    end if

    failures = 0
    mphi_max = 0
    call evaluate_susceptibility(kpar*vti/nui, -om_e/nui, &
        resolved_ion_ifunc_conservation_model, moments)
    call test_single_radial_wave()
    call test_radial_interference()
    call test_constant_drive()
    call test_negative_wave_harmonics()
    call test_zero_padded_spectrum()
    if (failures /= 0) error stop 'periodic ion spectral tensor regression failed'
    print *, 'periodic ion spectral tensor tests passed'

contains

    subroutine evaluate_spectrum(phi_m, br, x, tensor)
        complex(dp), intent(in) :: phi_m(:), br
        real(dp), intent(in) :: x
        real(dp), intent(out) :: tensor(2,2)
        call compute_periodic_ion_tensor_spectrum(phi_m, br, period, x, ks, kpar, &
            vti, nui, omega_ci, 0.0_dp, om_e, b0, tensor)
    end subroutine evaluate_spectrum

    subroutine test_single_radial_wave()
        complex(dp) :: phi_m(3), fields(3)
        real(dp) :: got(2,2), expected(2,2)

        phi_m = (0.0_dp, 0.0_dp)
        phi_m(3) = cmplx(0.8_dp, -0.3_dp, dp)/sol
        fields = [phi_m(3), (0.0_dp,0.0_dp), (0.0_dp,0.0_dp)]
        call calc_ion_flr_harmonic(0, ks, 2.0_dp*pi/period, ks, 2.0_dp*pi/period, &
            vti, omega_ci, omega_ci, sol, b0, nui, fields, fields, moments, expected)
        call evaluate_spectrum(phi_m, (0.0_dp,0.0_dp), radius, got)
        call check_tensor('single radial wave retains its FLR argument', got, expected)
    end subroutine test_single_radial_wave

    subroutine test_radial_interference()
        complex(dp) :: phi_m(3), rotated(3), phase, product
        real(dp) :: got(2,2), other(2,2), expected11, prefactor, bfirst, bsecond, bcross
        real(dp) :: wfirst, wsecond, wcross
        integer :: radial_sign, wave_index

        do radial_sign = -1, 1, 2
            wave_index = 2 + radial_sign
            phi_m = (0.0_dp, 0.0_dp)
            phi_m(2) = cmplx(0.7_dp, 0.2_dp, dp)/sol
            phi_m(wave_index) = cmplx(-0.4_dp, 0.6_dp, dp)/sol
            bfirst = (vti/omega_ci)**2*ks**2
            bsecond = (vti/omega_ci)**2*(ks**2+(2.0_dp*pi/period)**2)
            bcross = sqrt(bfirst*bsecond)
            wfirst = exp(-bfirst)*bessel_in(0,bfirst)
            wsecond = exp(-bsecond)*bessel_in(0,bsecond)
            wcross = exp(-0.5_dp*(bfirst+bsecond))*bessel_in(0,bcross)
            prefactor = sol**2*ks**2*real(moments(0,0),dp)/(2.0_dp*nui*b0**2)
            phase = exp(cmplx(0.0_dp, real(radial_sign,dp)*2.0_dp*pi*radius/period, dp))
            product = conjg(phi_m(2))*phi_m(wave_index)*phase
            expected11 = prefactor*(abs(phi_m(2))**2*wfirst + abs(phi_m(wave_index))**2*wsecond &
                + 2.0_dp*real(product,dp)*wcross)
            call evaluate_spectrum(phi_m, (0.0_dp,0.0_dp), radius, got)
            call check_scalar('two-wave D11 includes complex spatial interference', &
                got(1,1), expected11)

            rotated = phi_m*cmplx(0.6_dp,0.8_dp,dp)
            call evaluate_spectrum(rotated, (0.0_dp,0.0_dp), radius, other)
            call check_tensor('common phase invariance', other, got)
            call evaluate_spectrum(2.0_dp*phi_m, (0.0_dp,0.0_dp), radius, other)
            call check_tensor('quadratic amplitude scaling', other, 4.0_dp*got)

            phi_m(wave_index) = -phi_m(wave_index)
            expected11 = prefactor*(abs(phi_m(2))**2*wfirst + abs(phi_m(wave_index))**2*wsecond &
                - 2.0_dp*real(product,dp)*wcross)
            call evaluate_spectrum(phi_m, (0.0_dp,0.0_dp), radius, other)
            call check_scalar('relative phase changes interference sign', other(1,1), expected11)
            if (maxval(abs(got-transpose(got))) > 1.0e-12_dp) then
                failures = failures + 1
                print *, 'FAIL: spectral tensor is not symmetric'
            end if
            if (min(got(1,1),got(2,2)) < -1.0e-12_dp .or. &
                    got(1,1)*got(2,2)-got(1,2)*got(2,1) < -1.0e-12_dp) then
                failures = failures + 1
                print *, 'FAIL: admissible spectral fixture lost positive semidefiniteness'
            end if
        end do
    end subroutine test_radial_interference

    subroutine test_constant_drive()
        complex(dp) :: phi_m(3), fields(3), br
        real(dp) :: got(2,2), expected(2,2)

        phi_m = (0.0_dp, 0.0_dp)
        phi_m(2) = cmplx(0.8_dp,0.3_dp,dp)/sol
        br = cmplx(0.2_dp,-0.1_dp,dp)
        fields = [phi_m(2), br, (0.0_dp,0.0_dp)]
        call compute_periodic_ion_tensor(fields, ks, 0.0_dp, kpar, vti, nui, omega_ci, &
            0.0_dp, om_e, b0, expected)
        call evaluate_spectrum(phi_m, br, radius, got)
        call check_tensor('constant drive retains zero-mode normalization', got, expected)
    end subroutine test_constant_drive

    subroutine test_negative_wave_harmonics()
        complex(dp) :: phi_m(3), fields(3), harmonic_moments(0:3,0:3)
        real(dp) :: got(2,2), expected(2,2), term(2,2), zero_harmonic(2,2)
        real(dp), parameter :: omega_mode = 0.23_dp, kr_negative = -2.0_dp*pi/period
        integer :: ell

        phi_m = (0.0_dp, 0.0_dp)
        phi_m(1) = cmplx(-0.5_dp,0.9_dp,dp)/sol
        fields = [phi_m(1), (0.0_dp,0.0_dp), (0.0_dp,0.0_dp)]
        expected = 0.0_dp
        mphi_max = 2
        do ell = -2, 2
            call evaluate_susceptibility(kpar*vti/nui, &
                -(om_e+real(ell,dp)*omega_ci-omega_mode)/nui, &
                resolved_ion_ifunc_conservation_model, harmonic_moments)
            call calc_ion_flr_harmonic(ell, ks, kr_negative, ks, kr_negative, &
                vti, omega_ci, omega_ci, sol, b0, nui, fields, fields, harmonic_moments, term)
            expected = expected + term
            if (ell == 0) zero_harmonic = term
        end do
        call compute_periodic_ion_tensor_spectrum(phi_m, (0.0_dp,0.0_dp), period, radius, &
            ks, kpar, vti, nui, omega_ci, omega_mode, om_e, b0, got)
        call check_tensor('negative radial wave sums every configured cyclotron harmonic', &
            got, expected)
        if (maxval(abs(expected-zero_harmonic)) <= 1.0e-8_dp*maxval(abs(expected))) then
            failures = failures + 1
            print *, 'FAIL: harmonic fixture does not distinguish ell=0 truncation'
        end if
        mphi_max = 0
    end subroutine test_negative_wave_harmonics

    subroutine test_zero_padded_spectrum()
        complex(dp) :: phi_m(3), padded(7), br
        real(dp) :: original(2,2), enlarged(2,2)
        real(dp), parameter :: sample_radii(3) = [0.0_dp, radius, -0.23_dp]
        integer :: i

        phi_m = [cmplx(0.2_dp,-0.5_dp,dp), cmplx(0.7_dp,0.1_dp,dp), &
            cmplx(-0.3_dp,0.4_dp,dp)]/sol
        padded = (0.0_dp, 0.0_dp)
        ! At fixed period, indices 3:5 of -3:3 retain the same -1:1 wave numbers.
        padded(3:5) = phi_m
        br = cmplx(0.15_dp,-0.07_dp,dp)
        mphi_max = 2
        do i = 1, size(sample_radii)
            call evaluate_spectrum(phi_m, br, sample_radii(i), original)
            call evaluate_spectrum(padded, br, sample_radii(i), enlarged)
            call check_tensor('zero padding preserves radial waves and constant Br normalization', &
                enlarged, original)
        end do
        mphi_max = 0
    end subroutine test_zero_padded_spectrum

    subroutine check_scalar(label, got, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, expected
        if (abs(got-expected) > 2.0e-11_dp*max(abs(expected),1.0e-14_dp)) then
            failures = failures + 1
            print *, 'FAIL: ', label, ' got=',got,' expected=',expected
        end if
    end subroutine check_scalar

    subroutine check_tensor(label, got, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got(2,2), expected(2,2)
        if (maxval(abs(got-expected)) > 2.0e-11_dp*max(maxval(abs(expected)),1.0e-14_dp)) then
            failures = failures + 1
            print *, 'FAIL: ', label
            print *, 'got=',got,' expected=',expected
        end if
    end subroutine check_tensor
end program test_periodic_ion_tensor
