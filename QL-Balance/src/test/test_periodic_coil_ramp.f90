program test_periodic_coil_ramp
    use QLBalance_kinds, only: dp
    use periodic_coil_ramp_m, only: coil_ramp_target
    implicit none

    real(dp), parameter :: I_full = 160.8456083_dp, af_max = 7.910800579_dp

    ! Full coil current at af = af_max.
    if (abs(coil_ramp_target(I_full, af_max, af_max) - I_full) > 1.0e-12_dp*I_full) &
        error stop 'target must equal the full current at af = af_max'
    ! The target current scales with the coil current, i.e. sqrt(af).
    if (abs(coil_ramp_target(I_full, 0.25_dp*af_max, af_max) - 0.5_dp*I_full) &
            > 1.0e-12_dp*I_full) error stop 'target must scale with sqrt(af)'
    ! Ramp start of ramp_up_mode 0/2 (af = 1e-4).
    if (abs(coil_ramp_target(I_full, 1.0e-4_dp, af_max) &
            - I_full*sqrt(1.0e-4_dp/af_max)) > 1.0e-12_dp*I_full) &
        error stop 'wrong target at the ramp start'
    ! The paper stops at af = antenna_max_stopping * af_max = 2 af_max.
    if (abs(coil_ramp_target(I_full, 2.0_dp*af_max, af_max) - sqrt(2.0_dp)*I_full) &
            > 1.0e-12_dp*I_full) error stop 'wrong target beyond af_max'

    write (*, *) 'test_periodic_coil_ramp passed'
end program test_periodic_coil_ramp
