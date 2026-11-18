module periodic_coil_ramp_m
    !! Coil-current ramp for the periodic KIM TimeEvolution.
    !!
    !! In periodic mode every response is normalised to the target current, so
    !! the legacy antenna factor (a field-squared multiplier on D) cannot ramp
    !! the drive. With kim_periodic_coil_ramp the ramp routines (ramp_up_mode)
    !! still set antenna_factor exactly as for KiLCA, and the target current
    !! follows the coil current:
    !!     I_target = I_full * sqrt(antenna_factor / antenna_factor_max),
    !! where I_full is the namelist I_par_toroidal (GPEC current at the full
    !! coil current) and antenna_factor_max the namelist antenna_factor. The
    !! transport coefficients are not scaled again (they carry |s|^2).
    use QLBalance_kinds, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: coil_ramp_target

contains

    function coil_ramp_target(I_full, antenna_factor, antenna_factor_max) result(I_target)
        real(dp), intent(in) :: I_full, antenna_factor, antenna_factor_max
        real(dp) :: I_target

        if (.not. (ieee_is_finite(I_full) .and. ieee_is_finite(antenna_factor) .and. &
                   ieee_is_finite(antenna_factor_max))) &
            error stop 'coil ramp: non-finite input'
        if (I_full <= 0.0_dp) error stop 'coil ramp: I_par_toroidal must be > 0'
        if (antenna_factor_max <= 0.0_dp) error stop 'coil ramp: antenna_factor (max) must be > 0'
        if (antenna_factor <= 0.0_dp) error stop 'coil ramp: ramp antenna factor must be > 0'
        I_target = I_full*sqrt(antenna_factor/antenna_factor_max)
    end function coil_ramp_target

end module periodic_coil_ramp_m
