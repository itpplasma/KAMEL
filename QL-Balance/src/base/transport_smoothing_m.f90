module transport_smoothing_m
    use QLBalance_kinds, only: dp
    implicit none
    private
    public :: smooth_transport_profile

contains

    subroutine smooth_transport_profile(values, periodic)
        !! Periodic responses already carry the physical compact transition.
        !! A second spatial filter would alter the core and spread transport
        !! outside that support. Keep the established filter for other solvers.
        real(dp), intent(inout) :: values(:)
        logical, intent(in) :: periodic
        real(dp) :: smoothed(size(values))

        if (periodic) return
        call smooth_array_gauss(size(values), 30, values, smoothed)
        values = smoothed
    end subroutine smooth_transport_profile

end module transport_smoothing_m
