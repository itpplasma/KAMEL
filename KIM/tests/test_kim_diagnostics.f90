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
