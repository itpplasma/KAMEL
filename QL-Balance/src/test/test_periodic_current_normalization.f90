program test_periodic_current_normalization
    use QLBalance_kinds, only: dp
    use periodic_current_normalization_m, only: integrate_trusted_current, periodic_drive_scale
    implicit none
    real(dp), parameter :: r(5) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    complex(dp), parameter :: j(5) = [(1.0_dp,0.0_dp), (1.0_dp,0.0_dp), &
                                      (1.0_dp,0.0_dp), (1.0_dp,0.0_dp), (1.0_dp,0.0_dp)]
    complex(dp) :: unit_current, scale
    integer :: status
    unit_current = integrate_trusted_current(r, j, 1.0_dp, 3.0_dp)
    if (abs(unit_current - 4.0_dp) > 1.0e-12_dp) &
        error stop 'trusted current integral geometry'
    call periodic_drive_scale(8.0_dp, unit_current, 1.0_dp, 1.0e-20_dp, 1.0e6_dp, 1.0_dp, scale, status)
    if (status /= 0 .or. abs(scale - cmplx(1.0_dp/acos(-1.0_dp),0.0_dp,dp)) > 1.0e-12_dp) &
        error stop 'target current scale'
    call periodic_drive_scale(16.0_dp, unit_current, 1.0_dp, 1.0e-20_dp, 1.0e6_dp, 1.0_dp, scale, status)
    if (status /= 0 .or. abs(scale - 2.0_dp*cmplx(1.0_dp/acos(-1.0_dp),0.0_dp,dp)) > 1.0e-12_dp) &
        error stop 'target current doubling'
    call periodic_drive_scale(8.0_dp, (0.0_dp,0.0_dp), 1.0_dp, 1.0e-20_dp, 1.0e6_dp, 1.0_dp, scale, status)
    if (status /= 2 .or. abs(scale) > 0.0_dp) error stop 'current floor guard'
    call periodic_drive_scale(8.0_dp, unit_current, 1.0_dp, 1.0e-20_dp, 1.0e-2_dp, 1.0_dp, scale, status)
    if (status /= 3 .or. abs(scale) > 0.0_dp) error stop 'maximum scale guard'
    call periodic_drive_scale(-1.0_dp, unit_current, 1.0_dp, 1.0e-20_dp, 1.0e6_dp, 1.0_dp, scale, status)
    if (status /= 1 .or. abs(scale-cmplx(1.0_dp,0.0_dp,dp)) > 0.0_dp) &
        error stop 'invalid-input guard'
    call periodic_drive_scale(8.0_dp, unit_current, 1.0_dp, 1.0e-20_dp, 1.0e6_dp, 0.5_dp, scale, status)
    if (status /= 0 .or. abs(scale-cmplx(0.5_dp*(1.0_dp+1.0_dp/acos(-1.0_dp)),0.0_dp,dp)) &
            > 1.0e-12_dp) error stop 'normalization relaxation'
    print *, 'periodic current normalization tests passed'
end program test_periodic_current_normalization
