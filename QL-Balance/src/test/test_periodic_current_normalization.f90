program test_periodic_current_normalization
    use QLBalance_kinds, only: dp
    use periodic_current_normalization_m, only: integrate_trusted_current, periodic_drive_scale
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, &
                                                ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    implicit none
    integer :: failures
    real(dp) :: bad(3)
    character(len=32) :: argument

    failures = 0
    bad = [ieee_value(0.0_dp, ieee_quiet_nan), ieee_value(0.0_dp, ieee_positive_inf), &
           ieee_value(0.0_dp, ieee_negative_inf)]
    call get_command_argument(1, argument)
    if (len_trim(argument) > 0) then
        call invalid_geometry(trim(argument))
        stop 0
    end if
    call test_clipped_integrals()
    call test_grid_refinement()
    call test_complex_scaling()
    call test_guards()
    if (failures /= 0) then
        print *, 'Failed current normalization checks:', failures
        error stop 'periodic current normalization regression'
    end if
    print *, 'periodic current normalization tests passed'

contains

    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, 'FAIL: ', label
            failures = failures + 1
        end if
    end subroutine check

    subroutine check_complex(got, expected, label)
        complex(dp), intent(in) :: got, expected
        character(*), intent(in) :: label
        call check(ieee_is_finite(real(got)) .and. ieee_is_finite(aimag(got)), label//' finite')
        call check(abs(got - expected) <= 1.0e-12_dp * max(abs(expected), 1.0e-30_dp), label)
    end subroutine check_complex

    subroutine test_clipped_integrals()
        real(dp), parameter :: r(5) = [0.0_dp, 0.8_dp, 1.4_dp, 2.2_dp, 4.0_dp]
        complex(dp), parameter :: phase = cmplx(0.6_dp, -0.8_dp, dp)
        complex(dp) :: j(5), current
        integer :: i
        j = phase
        current = integrate_trusted_current(r, j, 1.1_dp, 2.1_dp)
        call check_complex(current, 1.6_dp * phase, 'off-grid complex cylindrical integral')
        current = integrate_trusted_current(r, j, 1.6_dp, 1.9_dp)
        call check_complex(current, 0.525_dp * phase, 'core strictly inside one cell')
        current = integrate_trusted_current(r, j, r(1), r(5))
        call check_complex(current, 8.0_dp * phase, 'full sampled interval')
        do i = 1, size(bad)
            j = phase
            j(3) = cmplx(bad(i), 0.0_dp, dp)
            current = integrate_trusted_current(r, j, 1.1_dp, 2.1_dp)
            call check(ieee_is_nan(real(current)), 'nonfinite real current returns NaN')
            j(3) = cmplx(0.0_dp, bad(i), dp)
            current = integrate_trusted_current(r, j, 1.1_dp, 2.1_dp)
            call check(ieee_is_nan(real(current)), 'nonfinite imaginary current returns NaN')
        end do
    end subroutine test_clipped_integrals

    subroutine test_grid_refinement()
        real(dp), allocatable :: r(:)
        complex(dp), allocatable :: j(:)
        real(dp) :: errors(3)
        complex(dp), parameter :: phase = cmplx(1.0_dp, 0.4_dp, dp)
        complex(dp) :: current, exact
        integer :: level, n, i
        exact = phase * (2.31_dp**3 - 0.73_dp**3) / 3.0_dp
        do level = 1, 3
            n = 16 * 2**(level - 1)
            allocate (r(n + 1), j(n + 1))
            r = [(3.0_dp * real(i, dp) / real(n, dp), i=0, n)]
            j = phase * r
            current = integrate_trusted_current(r, j, 0.73_dp, 2.31_dp)
            errors(level) = abs(current - exact)
            deallocate (r, j)
        end do
        print *, 'Complex r*j quadrature refinement errors:', errors
        call check(all(errors(1:2) > 3.0_dp * errors(2:3)), 'quadrature converges at second order')
        call check(all(errors(1:2) < 5.0_dp * errors(2:3)), 'quadrature resolves nonzero error')
    end subroutine test_grid_refinement

    subroutine test_complex_scaling()
        real(dp), parameter :: c_light = 2.99792458e10_dp, pi = acos(-1.0_dp)
        complex(dp), parameter :: current = cmplx(3.0_dp, -4.0_dp, dp)
        complex(dp) :: scale, doubled, relaxed
        integer :: status
        call periodic_drive_scale(8.0_dp, current, c_light, 1.0e-20_dp, &
                                  1.0e12_dp, 1.0_dp, scale, status)
        call check(status == 0, 'complex physical target accepted')
        call check_complex(2.0_dp * pi * scale * current / c_light, cmplx(8.0_dp, 0.0_dp, dp), &
                           'normalized complex current reaches real target in c=1 CGS')
        call periodic_drive_scale(16.0_dp, current, c_light, 1.0e-20_dp, &
                                  1.0e12_dp, 1.0_dp, doubled, status)
        call check(status == 0, 'doubled target accepted')
        call check_complex(doubled, 2.0_dp * scale, 'target doubling doubles complex amplitude')
        call check(abs(abs(doubled)**2 / abs(scale)**2 - 4.0_dp) < 1.0e-12_dp, &
                   'target doubling quadruples quadratic response')
        call periodic_drive_scale(8.0_dp, current, c_light, 1.0e-20_dp, &
                                  1.0e12_dp, 0.5_dp, relaxed, status)
       call check_complex(relaxed, 0.5_dp * (scale + 1.0_dp), 'one-shot relaxation from unit scale')
    end subroutine test_complex_scaling

    subroutine test_guards()
        real(dp), parameter :: valid(5) = [8.0_dp, 1.0_dp, 1.0e-20_dp, 1.0e6_dp, 1.0_dp]
        real(dp) :: config(5)
        complex(dp) :: scale, current
        integer :: i, input, status
        do input = 1, size(valid)
            config = valid
            config(input) = 0.0_dp
            call evaluate(config, (4.0_dp, 0.0_dp), scale, status)
            call check(status == 1 .and. scale == (0.0_dp, 0.0_dp), &
                       'finite invalid configuration suppresses response')
            do i = 1, size(bad)
                config = valid
                config(input) = bad(i)
                call evaluate(config, (4.0_dp, 0.0_dp), scale, status)
                call check(status == 3 .and. scale == (0.0_dp, 0.0_dp), &
                           'every nonfinite configuration input suppresses response')
            end do
        end do
        config = valid
        config(5) = 1.1_dp
        call evaluate(config, (4.0_dp, 0.0_dp), scale, status)
        call check(status == 1 .and. scale == (0.0_dp, 0.0_dp), 'relaxation above one rejected')
        do i = 1, size(bad)
            current = cmplx(bad(i), 0.0_dp, dp)
            call evaluate(valid, current, scale, status)
            call check(status == 3 .and. scale == (0.0_dp, 0.0_dp), 'nonfinite real unit current')
            current = cmplx(0.0_dp, bad(i), dp)
            call evaluate(valid, current, scale, status)
         call check(status == 3 .and. scale == (0.0_dp, 0.0_dp), 'nonfinite imaginary unit current')
        end do
        call evaluate(valid, (0.0_dp, 0.0_dp), scale, status)
        call check(status == 2 .and. scale == (0.0_dp, 0.0_dp), 'zero current floor')
        call evaluate(valid, cmplx(valid(3), 0.0_dp, dp), scale, status)
        call check(status == 2 .and. scale == (0.0_dp, 0.0_dp), 'current at floor')
        config = valid
        config(4) = 1.0e-2_dp
        call evaluate(config, (4.0_dp, 0.0_dp), scale, status)
        call check(status == 3 .and. scale == (0.0_dp, 0.0_dp), 'excessive scale suppressed')
    end subroutine test_guards

    subroutine evaluate(config, current, scale, status)
        real(dp), intent(in) :: config(5)
        complex(dp), intent(in) :: current
        complex(dp), intent(out) :: scale
        integer, intent(out) :: status
        call periodic_drive_scale(config(1), current, config(2), config(3), &
                                  config(4), config(5), scale, status)
    end subroutine evaluate

    subroutine invalid_geometry(which)
        character(*), intent(in) :: which
        real(dp) :: r(3), lo, hi
        complex(dp) :: current, j(3)
        r = [0.0_dp, 1.0_dp, 2.0_dp]
        lo = 0.2_dp
        hi = 1.8_dp
        j = (1.0_dp, 0.0_dp)
        select case (which)
        case ('geometry_duplicate')
            r(2) = r(1)
        case ('geometry_descending')
            r = [2.0_dp, 1.0_dp, 0.0_dp]
        case ('geometry_nan')
            r(2) = bad(1)
        case ('geometry_inf')
            r(3) = bad(2)
        case ('bounds_nan')
            lo = bad(1)
        case ('bounds_inf')
            hi = bad(2)
        case ('bounds_reversed')
            lo = 1.5_dp
            hi = 0.5_dp
        case ('bounds_empty')
            lo = hi
        case ('bounds_outside')
            lo = -0.1_dp
        case ('grid_shape')
            current = integrate_trusted_current(r, j(:2), lo, hi)
            return
        case default
            error stop 'unknown geometry test'
        end select
        current = integrate_trusted_current(r, j, lo, hi)
        print *, 'Invalid geometry was accepted:', which, current
    end subroutine invalid_geometry
end program test_periodic_current_normalization
