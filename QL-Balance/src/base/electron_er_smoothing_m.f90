module electron_er_smoothing_m
    !! Optional radial smoothing of the E x B frequency seen by the electron
    !! drift-kinetic coefficients (calc_transport_coeffs_ornuhl). Near the
    !! resonant surface k_par -> 0 and the electron D has a sharp peak at
    !! om_E = 0; with om_E taken pointwise from the local gradients this lets
    !! grid-scale kinks grow. Smoothing only this input is a diagnostic switch:
    !! the Er passed to the wave code is unchanged.
    use QLBalance_kinds, only: dp
    implicit none
    private
    public :: gauss_smooth_cells

contains

    subroutine gauss_smooth_cells(values, sigma_cells)
        !! Gaussian smoothing in grid index: weights exp(-i^2 / (2 sigma^2)),
        !! |i| <= ceiling(4 sigma), renormalised where the window is cut by the
        !! array ends. sigma_cells <= 0 leaves the values unchanged.
        real(dp), intent(inout) :: values(:)
        real(dp), intent(in) :: sigma_cells
        real(dp) :: smoothed(size(values))
        real(dp), allocatable :: w(:)
        integer :: n, half, i, lo, hi

        if (sigma_cells <= 0.0_dp) return
        n = size(values)
        half = ceiling(4.0_dp*sigma_cells)
        allocate (w(-half:half))
        do i = -half, half
            w(i) = exp(-0.5_dp*(real(i, dp)/sigma_cells)**2)
        end do
        do i = 1, n
            lo = max(1, i - half)
            hi = min(n, i + half)
            smoothed(i) = sum(w(lo - i:hi - i)*values(lo:hi))/sum(w(lo - i:hi - i))
        end do
        values = smoothed
    end subroutine gauss_smooth_cells

end module electron_er_smoothing_m
