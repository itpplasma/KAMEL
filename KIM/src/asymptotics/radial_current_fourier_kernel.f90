module radial_current_fourier_kernel_m
    use KIM_kinds_m, only: dp

    implicit none
    private

    public :: scaled_bessel_harmonic

contains

    real(dp) function scaled_bessel_harmonic(ell, bplus, bcross) result(value)
        !> Return exp(-bplus) I_ell(bcross) without overflowing the unscaled
        !> modified Bessel function. Integer-order symmetry is applied before
        !> the large-argument expansion.
        use constants_m, only: pi
        use fortnum_special, only: bessel_in

        integer, intent(in) :: ell
        real(dp), intent(in) :: bplus, bcross
        real(dp), parameter :: asymptotic_threshold = 500.0_dp
        real(dp) :: mu, inverse_eight_x, series
        integer :: order

        order = abs(ell)
        if (bcross == 0.0_dp) then
            if (order == 0) then
                value = exp(-bplus)
            else
                value = 0.0_dp
            end if
            return
        end if

        if (bcross <= asymptotic_threshold) then
            value = exp(-bplus) * bessel_in(order, bcross)
            return
        end if

        mu = 4.0_dp * real(order * order, dp)
        inverse_eight_x = 1.0_dp / (8.0_dp * bcross)
        series = 1.0_dp &
            - (mu - 1.0_dp) * inverse_eight_x &
            + 0.5_dp * (mu - 1.0_dp) * (mu - 9.0_dp) * inverse_eight_x**2 &
            - (mu - 1.0_dp) * (mu - 9.0_dp) * (mu - 25.0_dp) * &
            inverse_eight_x**3 / 6.0_dp
        value = exp(bcross - bplus) * series / sqrt(2.0_dp * pi * bcross)
    end function scaled_bessel_harmonic

end module radial_current_fourier_kernel_m
