program test_periodize
    ! Behavioral checks of the radial periodization: untouched resonant
    ! layer, exact L-periodicity, pass-through of an already-periodic
    ! profile, and seam smoothness.

    use KIM_kinds_m, only: dp
    use constants_m, only: pi
    use periodize_m, only: periodized_profile_value

    implicit none

    integer, parameter :: npts = 381
    real(dp), parameter :: r_lo = -10.0d0
    real(dp), parameter :: r_hi = 85.0d0
    real(dp), parameter :: r_mid = 36.5d0
    real(dp), parameter :: dr_layer = 5.0d0
    real(dp), parameter :: dr_transition = 10.0d0
    real(dp), parameter :: period = 2.0d0*(dr_layer + dr_transition)

    real(dp) :: r_grid(npts), quadratic(npts), harmonic(npts)
    real(dp) :: x, got, want, slope_left, slope_right, h
    integer :: i

    do i = 1, npts
        r_grid(i) = r_lo + (r_hi - r_lo)*(i - 1)/(npts - 1)
    end do
    quadratic = r_grid**2
    harmonic = sin(2.0d0*pi*(r_grid - r_mid)/period)

    ! Inside |x - r_mid| <= dr_layer the profile is unmodified;
    ! four-point Lagrange interpolation is exact for a quadratic.
    do i = 0, 10
        x = r_mid - dr_layer + 2.0d0*dr_layer*i/10.0d0
        got = periodized_profile_value(r_grid, quadratic, r_mid, dr_layer, &
                                       dr_transition, x)
        call require_close("resonant layer is unmodified", got, x**2, 1.0d-12)
    end do

    ! Exact L-periodicity across several periods.
    do i = 0, 12
        x = r_mid - 2.5d0*period + 5.0d0*period*i/12.0d0
        got = periodized_profile_value(r_grid, quadratic, r_mid, dr_layer, &
                                       dr_transition, x + period)
        want = periodized_profile_value(r_grid, quadratic, r_mid, dr_layer, &
                                        dr_transition, x)
        call require_close("construction is L-periodic", got, want, 1.0d-12)
    end do

    ! A profile that is already L-periodic passes through unchanged.
    do i = 0, 12
        x = r_mid - 1.4d0*period + 2.8d0*period*i/12.0d0
        got = periodized_profile_value(r_grid, harmonic, r_mid, dr_layer, &
                                       dr_transition, x)
        want = sin(2.0d0*pi*(x - r_mid)/period)
        call require_close("periodic profile passes through", got, want, &
                           1.0d-8)
    end do

    ! The seam is smooth: one-sided slopes agree to the finite-difference
    ! error, far below the O(1) jump an unblended wrap would give.
    h = 1.0d-3
    x = r_mid - 0.5d0*period
    slope_left = (periodized_profile_value(r_grid, quadratic, r_mid, &
                                           dr_layer, dr_transition, x - h) &
                  - periodized_profile_value(r_grid, quadratic, r_mid, &
                                             dr_layer, dr_transition, x - 3.0d0*h)) &
                 /(2.0d0*h)
    slope_right = (periodized_profile_value(r_grid, quadratic, r_mid, &
                                            dr_layer, dr_transition, x + 3.0d0*h) &
                   - periodized_profile_value(r_grid, quadratic, r_mid, &
                                              dr_layer, dr_transition, x + h)) &
                  /(2.0d0*h)
    call require_close("seam slopes agree", slope_left, slope_right, 1.0d-3)

    ! The layer edge starts the transition zone: a localizer without
    ! vanishing edge derivatives would kink the profile here.
    x = r_mid + dr_layer
    slope_left = (periodized_profile_value(r_grid, quadratic, r_mid, &
                                           dr_layer, dr_transition, x - h) &
                  - periodized_profile_value(r_grid, quadratic, r_mid, &
                                             dr_layer, dr_transition, x - 3.0d0*h)) &
                 /(2.0d0*h)
    slope_right = (periodized_profile_value(r_grid, quadratic, r_mid, &
                                            dr_layer, dr_transition, x + 3.0d0*h) &
                   - periodized_profile_value(r_grid, quadratic, r_mid, &
                                              dr_layer, dr_transition, x + h)) &
                  /(2.0d0*h)
    call require_close("transition edge slopes agree", slope_left, &
                       slope_right, 1.0d-3)

    print *, "PASS all periodize checks"

contains

    subroutine require_close(name, got, want, rel_tol)

        implicit none

        character(*), intent(in) :: name
        real(dp), intent(in) :: got, want, rel_tol

        real(dp) :: scale

        scale = max(abs(want), 1.0d0)
        if (abs(got - want) > rel_tol*scale) then
            print *, "FAIL ", name
            print *, "  got  = ", got
            print *, "  want = ", want
            error stop
        end if

    end subroutine require_close

end program test_periodize
