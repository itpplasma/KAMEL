module periodize_m
    ! Radial periodization for the enforced-periodicity solver (issue #175,
    ! Kasilov note 2026-07-04): background profiles stay untouched inside
    ! the resonant layer |r - r_m| <= dr_layer and are blended with their
    ! period-shifted copies inside the transition zones so that the result
    ! is periodic with L = 2 (dr_layer + dr_transition) and keeps all
    ! derivatives continuous. Gradients must be computed from the
    ! periodized profile afterwards; periodizing a gradient is wrong in
    ! the transition zone (mathematica/27 in the programme repository).

    use KIM_kinds_m, only: dp
    use constants_m, only: pi

    implicit none

contains

    ! C-infinity transition weight of the prototype localizer: 1 below x1,
    ! 0 above x2, essential-singularity decay at both ends.
    real(dp) function localizer_weight(x1, x2, x)

        implicit none

        real(dp), intent(in) :: x1, x2, x

        real(dp) :: t

        t = (x - x1)/(x2 - x1)
        if (t <= 0.0d0) then
            localizer_weight = 1.0d0
        else if (t >= 1.0d0) then
            localizer_weight = 0.0d0
        else
            localizer_weight = exp(-2.0d0*pi/(1.0d0 - t)*exp(-sqrt(2.0d0)/t))
        end if

    end function localizer_weight

    ! Periodized profile value at x from a tabulation (r_grid, f) of the
    ! original profile. The tabulation must cover three half-periods on
    ! both sides of r_mid because the blend samples the shifted copies.
    real(dp) function periodized_profile_value(r_grid, f, r_mid, dr_layer, &
                                               dr_transition, x)

        implicit none

        real(dp), intent(in) :: r_grid(:), f(:)
        real(dp), intent(in) :: r_mid, dr_layer, dr_transition, x

        real(dp) :: half_period, period, x_leftboundary, x_inperiod
        real(dp) :: value_in, value_left, value_right, weight

        if (dr_layer <= 0.0d0 .or. dr_transition <= 0.0d0) then
            error stop "periodize_m: layer and transition widths must be positive"
        end if
        half_period = dr_layer + dr_transition
        period = 2.0d0*half_period
        if (r_grid(1) > r_mid - 3.0d0*half_period) then
            error stop "periodize_m: tabulation does not cover the left copies"
        end if
        if (r_grid(size(r_grid)) < r_mid + 3.0d0*half_period) then
            error stop "periodize_m: tabulation does not cover the right copies"
        end if

        x_leftboundary = r_mid - half_period
        x_inperiod = x_leftboundary + modulo(x - x_leftboundary, period)

        value_in = tabulated_value(r_grid, f, x_inperiod)
        value_left = tabulated_value(r_grid, f, x_inperiod - period)
        value_right = tabulated_value(r_grid, f, x_inperiod + period)

        weight = localizer_weight(r_mid + dr_layer, &
                                  r_mid + dr_layer + 2.0d0*dr_transition, &
                                  x_inperiod)
        periodized_profile_value = value_in*weight + value_left*(1.0d0 - weight)

        weight = localizer_weight(r_mid - dr_layer - 2.0d0*dr_transition, &
                                  r_mid - dr_layer, x_inperiod)
        periodized_profile_value = periodized_profile_value*(1.0d0 - weight) &
                                   + value_right*weight

    end function periodized_profile_value

    real(dp) function tabulated_value(r_grid, f, x)

        implicit none

        real(dp), intent(in) :: r_grid(:), f(:), x

        integer, parameter :: nlagr = 4
        integer, parameter :: nder = 0
        real(dp) :: coef(0:nder, nlagr)
        integer :: ir, ibeg, iend

        call binsrc(r_grid, 1, size(r_grid), x, ir)
        ibeg = max(1, ir - nlagr/2)
        iend = ibeg + nlagr - 1
        if (iend > size(r_grid)) then
            iend = size(r_grid)
            ibeg = iend - nlagr + 1
        end if
        call plag_coeff(nlagr, nder, x, r_grid(ibeg:iend), coef)
        tabulated_value = sum(coef(0, :)*f(ibeg:iend))

    end function tabulated_value

end module periodize_m
