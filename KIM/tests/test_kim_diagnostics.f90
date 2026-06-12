program test_kim_diagnostics

    use KIM_kinds_m, only: dp
    use kim_diagnostics_m, only: integrate_Ipar

    implicit none

    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
    real(dp), parameter :: r1 = 1.5_dp, r2 = 4.0_dp

    call test_constant_jpar()
    call test_linear_jpar()

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
        real(dp) :: err_re, err_im

        err_re = abs(real(actual, dp) - real(expected, dp)) / abs(real(expected, dp))
        err_im = abs(aimag(actual) - aimag(expected)) / abs(aimag(expected))

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
