module numerics_utils_m

    use KIM_kinds_m, only: dp
    use constants_m, only: pi
    use gsl_mod, only: erf => gsl_sf_erf

    implicit none

    contains

    function erf_diff(z1, z2) result(d)
        ! Numerically stable difference of error functions.
        ! Uses linearized form for small deltas to avoid cancellation.
        real(dp), intent(in) :: z1, z2
        real(dp) :: d
        real(dp) :: dz, zm

        dz = z1 - z2

        if (abs(dz) < 1.0e-6_dp) then
            zm = 0.5d0 * (z1 + z2)
            d = 2.0d0 / sqrt(pi) * exp(-zm*zm) * dz
        else
            d = erf(z1) - erf(z2)
        end if
    end function erf_diff

end module numerics_utils_m
