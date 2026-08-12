program test_quasilinear_integral_algebra
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use KIM_kinds_m, only: dp
    use quasilinear_integral_m, only: calc_ion_integral_harmonic, &
        calc_ion_integral_harmonic_debug, build_transverse_moments, &
        QL_PHI, QL_BR, QL_BPAR, &
        QL_INTEGRAL_ALGEBRA_VERSION
    implicit none

    call test_mathematica_oracles()
    call test_zero_bparallel_removes_five_pairs()
    call test_bparallel_long_wavelength()
    call test_exact_zero_kperp()
    call test_large_flr_is_finite()
    call test_species_sign_convention()
    call test_explicit_hermitianization()
    call test_reciprocity_and_tensor_properties()
    call test_phase_invariance()
    call test_harmonic_convergence()
    call test_flr_convergence()

    print *, 'quasilinear integral algebra tests passed'

contains

    subroutine set_oracle_inputs(fields_s, fields_o, ifunc)
        complex(dp), intent(out) :: fields_s(3), fields_o(3)
        complex(dp), intent(out) :: ifunc(0:3,0:3)
        integer :: p, q

        fields_s(QL_PHI) = cmplx(3.0_dp/5.0_dp, -2.0_dp/7.0_dp, dp)
        fields_s(QL_BR) = cmplx(-4.0_dp/9.0_dp, 1.0_dp/6.0_dp, dp)
        fields_s(QL_BPAR) = cmplx(5.0_dp/8.0_dp, 3.0_dp/11.0_dp, dp)
        fields_o(QL_PHI) = cmplx(-2.0_dp/3.0_dp, 1.0_dp/5.0_dp, dp)
        fields_o(QL_BR) = cmplx(7.0_dp/10.0_dp, 2.0_dp/9.0_dp, dp)
        fields_o(QL_BPAR) = cmplx(-3.0_dp/7.0_dp, -1.0_dp/4.0_dp, dp)
        do p = 0, 3
            do q = 0, 3
                ifunc(p,q) = cmplx(real(p+1,dp)/real(q+2,dp), &
                    real(p-q,dp)/13.0_dp, dp)
            end do
        end do
    end subroutine set_oracle_inputs

    subroutine test_mathematica_oracles()
        complex(dp) :: fields_s(3), fields_o(3), ifunc(0:3,0:3)
        complex(dp) :: transverse(0:2,3,3), channels(2,2,3,3), reference
        real(dp) :: tensor(2,2), production_tensor(2,2)
        real(dp) :: re_reference, im_reference
        character(len=160) :: header
        integer :: unit, ios, row, q, observation, source, i, j

        call build_transverse_moments(2, 0.7_dp, -0.4_dp, 1.1_dp, 0.3_dp, &
            1.2_dp, 1.8_dp, -1.8_dp, 1.3_dp, 0.9_dp, transverse)
        open(newunit=unit, file='quasilinear_transverse_oracle.dat', &
            status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'cannot open Mathematica transverse oracle'
        read(unit,'(A)') header
        call assert_oracle_version(header)
        read(unit,'(A)') header
        do row = 1, 27
            read(unit,*,iostat=ios) q, observation, source, re_reference, im_reference
            if (ios /= 0) error stop 'invalid Mathematica transverse oracle'
            reference = cmplx(re_reference, im_reference, dp)
            call assert_complex('Mathematica transverse coefficient', &
                transverse(q,observation,source), reference, 3.0e-13_dp)
        end do
        close(unit)

        call set_oracle_inputs(fields_s, fields_o, ifunc)
        call calc_ion_integral_harmonic_debug(2, 0.7_dp, -0.4_dp, 1.1_dp, 0.3_dp, &
            1.2_dp, 1.8_dp, -1.8_dp, 1.3_dp, 0.9_dp, 1.4_dp, fields_s, &
            fields_o, ifunc, tensor, channels)
        open(newunit=unit, file='quasilinear_channel_oracle.dat', &
            status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'cannot open Mathematica channel oracle'
        read(unit,'(A)') header
        call assert_oracle_version(header)
        read(unit,'(A)') header
        do row = 1, 36
            read(unit,*,iostat=ios) i, j, observation, source, &
                re_reference, im_reference
            if (ios /= 0) error stop 'invalid Mathematica channel oracle'
            reference = cmplx(re_reference, im_reference, dp)
            call assert_complex('Mathematica contracted channel', &
                channels(i,j,observation,source), reference, 5.0e-13_dp)
        end do
        close(unit)

        call calc_ion_integral_harmonic(2, 0.7_dp, -0.4_dp, 1.1_dp, 0.3_dp, &
            1.2_dp, 1.8_dp, -1.8_dp, 1.3_dp, 0.9_dp, 1.4_dp, fields_s, &
            fields_o, ifunc, production_tensor)
        call assert_array_real('optional debug channels leave production tensor', &
            production_tensor, tensor, 3.0e-15_dp)
        do i = 1, 2
            do j = 1, 2
                call assert_relative('channel sum is production tensor', tensor(i,j), &
                    real(sum(channels(i,j,:,:)),dp), 3.0e-15_dp)
            end do
        end do
    end subroutine test_mathematica_oracles

    subroutine test_zero_bparallel_removes_five_pairs()
        complex(dp) :: fields(3), fields_o(3), ifunc(0:3,0:3)
        complex(dp) :: channels_all(2,2,3,3), channels_no(2,2,3,3), removed
        real(dp) :: tensor_all(2,2), tensor_no(2,2)
        integer :: i, j

        call set_oracle_inputs(fields, fields_o, ifunc)
        fields_o = fields
        call calc_ion_integral_harmonic_debug(1, 0.7_dp, -0.4_dp, 1.1_dp, 0.3_dp, &
            1.2_dp, 1.8_dp, -1.8_dp, 1.3_dp, 0.9_dp, 1.4_dp, fields, &
            fields_o, ifunc, tensor_all, channels_all)
        fields(QL_BPAR) = (0.0_dp, 0.0_dp)
        fields_o(QL_BPAR) = (0.0_dp, 0.0_dp)
        call calc_ion_integral_harmonic_debug(1, 0.7_dp, -0.4_dp, 1.1_dp, 0.3_dp, &
            1.2_dp, 1.8_dp, -1.8_dp, 1.3_dp, 0.9_dp, 1.4_dp, fields, &
            fields_o, ifunc, tensor_no, channels_no)

        if (maxval(abs(channels_no(:,:,QL_BPAR,:))) > 0.0_dp) &
            error stop 'zero Bparallel left observation channels'
        if (maxval(abs(channels_no(:,:,:,QL_BPAR))) > 0.0_dp) &
            error stop 'zero Bparallel left source channels'
        do i = 1, 2
            do j = 1, 2
                removed = sum(channels_all(i,j,QL_BPAR,:)) &
                    + sum(channels_all(i,j,:,QL_BPAR)) &
                    - channels_all(i,j,QL_BPAR,QL_BPAR)
                call assert_relative('exactly five Bparallel pairs removed', &
                    tensor_all(i,j)-tensor_no(i,j), real(removed,dp), 3.0e-14_dp)
            end do
        end do
    end subroutine test_zero_bparallel_removes_five_pairs

    subroutine test_bparallel_long_wavelength()
        complex(dp) :: fields(3), ifunc(0:3,0:3), channels(2,2,3,3)
        real(dp) :: tensor(2,2), expected, expected_phi, expected_cross
        real(dp), parameter :: vT = 3.0e7_dp, omega_c = 4.0e7_dp
        real(dp), parameter :: c_light = 2.99792458e10_dp
        real(dp), parameter :: B0 = 1.8e4_dp, nu = 2.0e5_dp
        real(dp), parameter :: ks = 1.0e-10_dp

        fields = (0.0_dp, 0.0_dp)
        ifunc = (0.0_dp, 0.0_dp)
        ifunc(0,0) = (1.0_dp, 0.0_dp)
        fields(QL_BPAR) = cmplx(0.7_dp, -0.2_dp, dp)
        call calc_ion_integral_harmonic_debug(0, ks, 0.0_dp, ks, 0.0_dp, vT, &
            omega_c, -omega_c, c_light, B0, nu, fields, fields, ifunc, &
            tensor, channels)
        expected = vT**4*ks**2*abs(fields(QL_BPAR))**2 &
            /(nu*B0**2*omega_c**2)
        call assert_relative('Bparallel-only long-wave D11', tensor(1,1), &
            expected, 2.0e-8_dp)
        call assert_relative('isolated Bparallel long-wave channel', &
            real(channels(1,1,QL_BPAR,QL_BPAR),dp), expected, 2.0e-8_dp)

        ! Add Phi with a non-trivial relative phase.  The two ordered
        ! Phi-Bparallel channels combine to the signed real interference term.
        fields(QL_PHI) = cmplx(5.0e-4_dp,-3.0e-4_dp,dp)
        call calc_ion_integral_harmonic_debug(0, ks, 0.0_dp, ks, 0.0_dp, vT, &
            omega_c, -omega_c, c_light, B0, nu, fields, fields, ifunc, &
            tensor, channels)
        expected_phi = c_light**2*ks**2*abs(fields(QL_PHI))**2 &
            /(2.0_dp*nu*B0**2)
        expected_cross = -c_light*vT**2*ks**2 &
            *real(conjg(fields(QL_PHI))*fields(QL_BPAR),dp) &
            /(nu*B0**2*omega_c)
        call assert_relative('signed relative-phase long-wave D11', tensor(1,1), &
            expected+expected_phi+expected_cross, 3.0e-8_dp)
    end subroutine test_bparallel_long_wavelength

    subroutine test_exact_zero_kperp()
        complex(dp) :: transverse(0:2,3,3)
        real(dp), parameter :: vT = 1.2_dp, B0 = 0.9_dp

        call build_transverse_moments(0, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            vT, 1.8_dp, -1.8_dp, 1.3_dp, B0, transverse)
        call assert_complex('zero-kperp Br-Br W0', transverse(0,QL_BR,QL_BR), &
            cmplx(vT**2/B0**2,0.0_dp,dp), 3.0e-14_dp)
        call assert_complex('zero-kperp Br-Br W1', transverse(1,QL_BR,QL_BR), &
            cmplx(vT**2/B0**2,0.0_dp,dp), 3.0e-14_dp)
        call assert_complex('zero-kperp Br-Br W2', transverse(2,QL_BR,QL_BR), &
            cmplx(2.0_dp*vT**2/B0**2,0.0_dp,dp), 3.0e-14_dp)
        if (maxval(abs(transverse(:,QL_BPAR,:))) > 0.0_dp .or. &
            maxval(abs(transverse(:,: ,QL_BPAR))) > 0.0_dp) &
            error stop 'zero-kperp ell=0 left Bparallel derivative channels'

        ! For |ell|=1 the base gyro-factor vanishes at the origin while its
        ! source/observation mixed derivative is finite.  This is the branch
        ! most easily lost by evaluating polar angles or 0/0 ratios.
        call build_transverse_moments(1, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            vT, 1.8_dp, -1.8_dp, 1.3_dp, B0, transverse)
        call assert_complex('zero-kperp ell=1 Bparallel W0', &
            transverse(0,QL_BPAR,QL_BPAR), &
            cmplx(0.5_dp*vT**2/B0**2,0.0_dp,dp), 3.0e-14_dp)
        call assert_complex('zero-kperp ell=1 Bparallel W1', &
            transverse(1,QL_BPAR,QL_BPAR), &
            cmplx(vT**2/B0**2,0.0_dp,dp), 3.0e-14_dp)
        call assert_complex('zero-kperp ell=1 Bparallel W2', &
            transverse(2,QL_BPAR,QL_BPAR), &
            cmplx(3.0_dp*vT**2/B0**2,0.0_dp,dp), 3.0e-14_dp)
    end subroutine test_exact_zero_kperp

    subroutine test_large_flr_is_finite()
        use fortnum_special, only: bessel_in
        complex(dp) :: transverse(0:2,3,3)
        real(dp) :: leading_asymptotic, expected, wave_number

        call build_transverse_moments(0, 30.0_dp, 0.0_dp, 30.0_dp, 0.0_dp, &
            1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, transverse)
        if (.not. all(ieee_is_finite(real(transverse,dp))) .or. &
            .not. all(ieee_is_finite(aimag(transverse)))) &
            error stop 'large finite-Larmor-radius block is non-finite'
        leading_asymptotic = 1.0_dp/sqrt(2.0_dp*acos(-1.0_dp)*900.0_dp)
        call assert_relative('large-FLR scaled Bessel product', &
            real(transverse(0,QL_BR,QL_BR),dp), leading_asymptotic, 2.0e-4_dp)

        wave_number = sqrt(501.0_dp)
        call build_transverse_moments(50, wave_number, 0.0_dp, wave_number, &
            0.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, transverse)
        expected = exp(-501.0_dp)*bessel_in(50,501.0_dp)
        call assert_relative('large-FLR high-harmonic scaled Bessel product', &
            real(transverse(0,QL_BR,QL_BR),dp), expected, 2.0e-11_dp)
    end subroutine test_large_flr_is_finite

    subroutine test_species_sign_convention()
        complex(dp) :: fields_s(3), fields_o(3), ifunc(0:3,0:3)
        complex(dp) :: minus_channels(2,2,3,3), plus_channels(2,2,3,3)
        real(dp) :: minus_tensor(2,2), plus_tensor(2,2)
        integer :: i, j, observation, source, bpar_count

        call set_oracle_inputs(fields_s, fields_o, ifunc)
        call calc_ion_integral_harmonic_debug(2, 0.7_dp, -0.4_dp, 1.1_dp, 0.3_dp, &
            1.2_dp, 1.8_dp, -1.8_dp, 1.3_dp, 0.9_dp, 1.4_dp, fields_s, &
            fields_o, ifunc, minus_tensor, minus_channels)
        call calc_ion_integral_harmonic_debug(2, 0.7_dp, -0.4_dp, 1.1_dp, 0.3_dp, &
            1.2_dp, 1.8_dp, 1.8_dp, 1.3_dp, 0.9_dp, 1.4_dp, fields_s, &
            fields_o, ifunc, plus_tensor, plus_channels)
        do i = 1, 2
            do j = 1, 2
                do observation = 1, 3
                    do source = 1, 3
                        bpar_count = merge(1,0,observation == QL_BPAR) &
                            + merge(1,0,source == QL_BPAR)
                        if (mod(bpar_count,2) == 0) then
                            call assert_complex('even Bparallel species sign', &
                                plus_channels(i,j,observation,source), &
                                minus_channels(i,j,observation,source), 2.0e-14_dp)
                        else
                            call assert_complex('odd Bparallel species sign', &
                                plus_channels(i,j,observation,source), &
                                -minus_channels(i,j,observation,source), 2.0e-14_dp)
                        end if
                    end do
                end do
            end do
        end do
    end subroutine test_species_sign_convention

    subroutine test_explicit_hermitianization()
        complex(dp) :: fields(3), unused_fields(3), ifunc(0:3,0:3)
        complex(dp) :: channels(2,2,3,3)
        real(dp) :: tensor(2,2)
        integer :: i, j, observation, source

        ! set_oracle_inputs deliberately returns a non-Hermitian I(p,q).
        ! Equal source/observation waves and fields isolate the explicit
        ! (K + K^dagger)/2 construction inside the production routine.
        call set_oracle_inputs(fields, unused_fields, ifunc)
        call calc_ion_integral_harmonic_debug(2, 0.7_dp, -0.4_dp, 0.7_dp, -0.4_dp, &
            1.2_dp, 1.8_dp, -1.8_dp, 1.3_dp, 0.9_dp, 1.4_dp, fields, &
            fields, ifunc, tensor, channels)
        do i = 1, 2
            do j = 1, 2
                do observation = 1, 3
                    do source = 1, 3
                        call assert_complex('explicit quadratic Hermitianization', &
                            channels(i,j,observation,source), &
                            conjg(channels(j,i,source,observation)), 3.0e-13_dp)
                    end do
                end do
            end do
        end do
    end subroutine test_explicit_hermitianization

    subroutine test_reciprocity_and_tensor_properties()
        complex(dp) :: first(0:2,3,3), swapped(0:2,3,3)
        complex(dp) :: fields(3), ifunc(0:3,0:3), channels(2,2,3,3)
        real(dp) :: tensor(2,2), scale, xi
        integer :: q, observation, source, p, r, i, j

        call build_transverse_moments(2, 0.7_dp, -0.4_dp, 1.1_dp, 0.3_dp, &
            1.2_dp, 1.8_dp, -1.8_dp, 1.3_dp, 0.9_dp, first)
        call build_transverse_moments(2, 1.1_dp, 0.3_dp, 0.7_dp, -0.4_dp, &
            1.2_dp, 1.8_dp, -1.8_dp, 1.3_dp, 0.9_dp, swapped)
        do q = 0, 2
            do observation = 1, 3
                do source = 1, 3
                    call assert_complex('source-observation reciprocity', &
                        first(q,observation,source), &
                        conjg(swapped(q,source,observation)), 2.0e-13_dp)
                end do
            end do
        end do

        fields(QL_PHI) = cmplx(0.4_dp,-0.2_dp,dp)
        fields(QL_BR) = cmplx(-0.3_dp,0.1_dp,dp)
        fields(QL_BPAR) = cmplx(0.2_dp,0.35_dp,dp)
        xi = 0.6_dp
        do p = 0, 3
            do r = 0, 3
                ifunc(p,r) = cmplx(xi**(p+r),0.0_dp,dp)
            end do
        end do
        call calc_ion_integral_harmonic_debug(1, 0.7_dp, -0.4_dp, 0.7_dp, &
            -0.4_dp, 1.2_dp, 1.8_dp, -1.8_dp, 1.3_dp, 0.9_dp, 1.4_dp, &
            fields, fields, ifunc, tensor, channels)
        do i = 1, 2
            do j = 1, 2
                do observation = 1, 3
                    do source = 1, 3
                        call assert_complex('channel Hermiticity', &
                            channels(i,j,observation,source), &
                            conjg(channels(j,i,source,observation)), 3.0e-13_dp)
                    end do
                end do
            end do
        end do
        call assert_relative('Onsager D12-D21', tensor(1,2), tensor(2,1), &
            3.0e-13_dp)
        scale = max(1.0_dp,maxval(abs(tensor)))
        if (tensor(1,1) < -3.0e-13_dp*scale .or. &
            tensor(2,2) < -3.0e-13_dp*scale) error stop 'negative tensor diagonal'
        if (tensor(1,1)*tensor(2,2)-tensor(1,2)*tensor(2,1) &
            < -3.0e-13_dp*scale**2) error stop 'ion tensor is not positive semidefinite'
    end subroutine test_reciprocity_and_tensor_properties

    subroutine test_phase_invariance()
        complex(dp) :: fields_s(3), fields_o(3), phased_s(3), phased_o(3)
        complex(dp) :: ifunc(0:3,0:3), channels(2,2,3,3), phased_channels(2,2,3,3)
        complex(dp) :: phase
        real(dp) :: tensor(2,2), phased_tensor(2,2)

        call set_oracle_inputs(fields_s, fields_o, ifunc)
        phase = exp(cmplx(0.0_dp,0.731_dp,dp))
        phased_s = phase*fields_s
        phased_o = phase*fields_o
        call calc_ion_integral_harmonic_debug(2, 0.7_dp, -0.4_dp, 1.1_dp, 0.3_dp, &
            1.2_dp, 1.8_dp, -1.8_dp, 1.3_dp, 0.9_dp, 1.4_dp, fields_s, &
            fields_o, ifunc, tensor, channels)
        call calc_ion_integral_harmonic_debug(2, 0.7_dp, -0.4_dp, 1.1_dp, 0.3_dp, &
            1.2_dp, 1.8_dp, -1.8_dp, 1.3_dp, 0.9_dp, 1.4_dp, phased_s, &
            phased_o, ifunc, phased_tensor, phased_channels)
        call assert_array_complex('common field phase', phased_channels, channels, &
            3.0e-14_dp)
        call assert_array_real('common phase tensor', phased_tensor, tensor, 3.0e-14_dp)
    end subroutine test_phase_invariance

    subroutine test_harmonic_convergence()
        real(dp) :: sum8(2,2), sum10(2,2)
        call harmonic_sum(8, sum8)
        call harmonic_sum(10, sum10)
        call assert_array_real('nonzero-Bparallel harmonic convergence', &
            sum10, sum8, 2.0e-9_dp)
    end subroutine test_harmonic_convergence

    subroutine harmonic_sum(limit, tensor_sum)
        integer, intent(in) :: limit
        real(dp), intent(out) :: tensor_sum(2,2)
        complex(dp) :: fields(3), ifunc(0:3,0:3), channels(2,2,3,3)
        real(dp) :: tensor(2,2)
        integer :: harmonic, p, q

        fields(QL_PHI) = cmplx(0.4_dp,-0.2_dp,dp)
        fields(QL_BR) = cmplx(-0.3_dp,0.1_dp,dp)
        fields(QL_BPAR) = cmplx(0.2_dp,0.35_dp,dp)
        do p = 0, 3
            do q = 0, 3
                ifunc(p,q) = cmplx(0.55_dp**(p+q),0.0_dp,dp)
            end do
        end do
        tensor_sum = 0.0_dp
        do harmonic = -limit, limit
            call calc_ion_integral_harmonic_debug(harmonic, 0.7_dp, -0.4_dp, &
                0.7_dp, -0.4_dp, 1.2_dp, 1.8_dp, -1.8_dp, 1.3_dp, 0.9_dp, &
                1.4_dp, fields, fields, ifunc, tensor, channels)
            tensor_sum = tensor_sum + tensor
        end do
    end subroutine harmonic_sum

    subroutine test_flr_convergence()
        real(dp) :: error_large, error_small
        error_large = bparallel_longwave_error(0.2_dp)
        error_small = bparallel_longwave_error(0.1_dp)
        if (.not. (error_small < 0.35_dp*error_large)) &
            error stop 'Bparallel finite-Larmor-radius limit did not converge'
    end subroutine test_flr_convergence

    function bparallel_longwave_error(ks) result(relative_error)
        real(dp), intent(in) :: ks
        real(dp) :: relative_error, tensor(2,2), expected
        complex(dp) :: fields(3), ifunc(0:3,0:3), channels(2,2,3,3)

        fields = (0.0_dp,0.0_dp)
        fields(QL_BPAR) = cmplx(0.7_dp,-0.2_dp,dp)
        ifunc = (0.0_dp,0.0_dp)
        ifunc(0,0) = (1.0_dp,0.0_dp)
        call calc_ion_integral_harmonic_debug(0, ks, 0.0_dp, ks, 0.0_dp, &
            1.2_dp, 1.8_dp, -1.8_dp, 1.3_dp, 0.9_dp, 1.4_dp, fields, &
            fields, ifunc, tensor, channels)
        expected = 1.2_dp**4*ks**2*abs(fields(QL_BPAR))**2 &
            /(1.4_dp*0.9_dp**2*1.8_dp**2)
        relative_error = abs(tensor(1,1)/expected-1.0_dp)
    end function bparallel_longwave_error

    subroutine assert_oracle_version(header)
        character(len=*), intent(in) :: header
        character(len=64) :: expected
        write(expected,'(a,i0)') '# algebra_version ', QL_INTEGRAL_ALGEBRA_VERSION
        if (trim(header) /= trim(expected)) then
            print *, 'FAIL: incompatible Mathematica oracle: ', trim(header)
            error stop
        end if
    end subroutine assert_oracle_version

    subroutine assert_complex(label, actual, reference, tolerance)
        character(len=*), intent(in) :: label
        complex(dp), intent(in) :: actual, reference
        real(dp), intent(in) :: tolerance
        if (abs(actual-reference) > &
            tolerance*max(abs(reference),sqrt(tiny(1.0_dp)))) then
            print *, 'FAIL: ', trim(label), actual, reference
            error stop
        end if
    end subroutine assert_complex

    subroutine assert_relative(label, actual, reference, tolerance)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual, reference, tolerance
        if (abs(actual-reference) > tolerance*max(abs(reference),tiny(1.0_dp))) then
            print *, 'FAIL: ', trim(label), actual, reference
            error stop
        end if
    end subroutine assert_relative

    subroutine assert_array_complex(label, actual, reference, tolerance)
        character(len=*), intent(in) :: label
        complex(dp), intent(in) :: actual(:,:,:,:), reference(:,:,:,:)
        real(dp), intent(in) :: tolerance
        if (maxval(abs(actual-reference)) > &
            tolerance*max(maxval(abs(reference)),sqrt(tiny(1.0_dp)))) then
            print *, 'FAIL: ', trim(label), maxval(abs(actual-reference))
            error stop
        end if
    end subroutine assert_array_complex

    subroutine assert_array_real(label, actual, reference, tolerance)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual(:,:), reference(:,:), tolerance
        if (maxval(abs(actual-reference)) > &
            tolerance*max(maxval(abs(reference)),sqrt(tiny(1.0_dp)))) then
            print *, 'FAIL: ', trim(label), maxval(abs(actual-reference))
            error stop
        end if
    end subroutine assert_array_real

end program test_quasilinear_integral_algebra
