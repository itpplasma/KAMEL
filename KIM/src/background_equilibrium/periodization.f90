module periodization_m
    !> Forced-periodicity utility vendored from the FORCED_PERIODICITY reference
    !> (make_periodic.f90, localizer.f90). Wraps an aperiodic profile into an
    !> L-periodic one by mirroring neighbouring periods and blending them with a
    !> C-infinity localizer weight across a transition band.
    !>
    !> The layout, for x_mid the period centre, is
    !>   half_period = dx_asis + dx_tr,  period L = 2*half_period,
    !> with an "as-is" core of half-width dx_asis (where the periodized functions
    !> equal the aperiodic ones exactly) flanked by transition bands of width
    !> dx_tr that smoothly hand over to the mirrored neighbouring periods.

    use KIM_kinds_m, only: dp
    implicit none
    private

    public :: make_periodic, localizer

    abstract interface
        !> Aperiodic sample functions evaluated at a single point x.
        subroutine aperfuns_i(nfuns, x, funs)
            import :: dp
            integer, intent(in) :: nfuns
            real(dp), intent(in) :: x
            real(dp), intent(out) :: funs(nfuns)
        end subroutine aperfuns_i
    end interface

contains

    !> Wrap aperiodic sample functions into an L-periodic set at point x.
    subroutine make_periodic(aperfuns, nfuns, x_mid, dx_asis, dx_tr, x, funs_per)
        procedure(aperfuns_i) :: aperfuns
        integer, intent(in) :: nfuns
        real(dp), intent(in) :: x_mid, dx_asis, dx_tr, x
        real(dp), intent(out) :: funs_per(nfuns)

        real(dp) :: funs_aper(nfuns)
        real(dp) :: period, half_period, x_inperiod, x_leftboundary, weight

        half_period = dx_asis + dx_tr
        period = 2.0_dp*half_period
        x_leftboundary = x_mid - half_period

        x_inperiod = x_leftboundary + modulo(x - x_leftboundary, period)

        call aperfuns(nfuns, x_inperiod, funs_per)

        call aperfuns(nfuns, x_inperiod - period, funs_aper)
        call localizer(1.0_dp, x_mid + dx_asis, x_mid + dx_asis + 2.0_dp*dx_tr, x_inperiod, weight)
        funs_per = funs_per*weight + funs_aper*(1.0_dp - weight)

        call aperfuns(nfuns, x_inperiod + period, funs_aper)
        call localizer(1.0_dp, x_mid - dx_asis - 2.0_dp*dx_tr, x_mid - dx_asis, x_inperiod, weight)
        funs_per = funs_per*(1.0_dp - weight) + funs_aper*weight
    end subroutine make_periodic

    !> C-infinity transition weight over [x1, x2]. sig>0 ramps from 1 to 0,
    !> otherwise from 0 to 1. The literals are 2*pi and sqrt(2).
    subroutine localizer(sig, x1, x2, x, weight)
        real(dp), intent(in) :: sig, x1, x2, x
        real(dp), intent(out) :: weight

        real(dp) :: t

        if (sig .gt. 0.0_dp) then
            t = (x - x1)/(x2 - x1)      ! from 1 to 0
        else
            t = (x2 - x)/(x2 - x1)      ! from 0 to 1
        end if

        if (t .le. 0.0_dp) then
            weight = 1.0_dp
        else if (t .ge. 1.0_dp) then
            weight = 0.0_dp
        else
            weight = exp(-6.283185307179586d0/(1.0_dp - t)*exp(-1.414213562373095d0/t))
        end if
    end subroutine localizer

end module periodization_m
