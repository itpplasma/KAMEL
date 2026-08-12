program test_periodic_amplitude_state
    use QLBalance_kinds, only: dp
    use periodic_amplitude_state_m, only: periodic_amplitude_state_t
    implicit none
    type(periodic_amplitude_state_t) :: state
    complex(dp) :: initial(2), trial(2), unit_initial(2), unit_trial(2)
    complex(dp) :: residual_initial(2), residual_trial(2)
    integer :: status_initial(2), status_trial(2)
    initial = [cmplx(1.0_dp,0.0_dp,dp), cmplx(2.0_dp,0.0_dp,dp)]
    trial = [cmplx(3.0_dp,1.0_dp,dp), cmplx(4.0_dp,-1.0_dp,dp)]
    unit_initial = [cmplx(5.0_dp,0.0_dp,dp), cmplx(6.0_dp,0.0_dp,dp)]
    unit_trial = [cmplx(7.0_dp,1.0_dp,dp), cmplx(8.0_dp,-1.0_dp,dp)]
    residual_initial = [cmplx(0.1_dp,0.2_dp,dp), cmplx(0.3_dp,0.4_dp,dp)]
    residual_trial = [cmplx(0.5_dp,0.6_dp,dp), cmplx(0.7_dp,0.8_dp,dp)]
    status_initial = [0, 2]
    status_trial = [1, 3]
    call state%initialize(initial, unit_initial, residual_initial, status_initial, 9.0_dp, 0.75_dp)
    call state%begin_trial(trial, unit_trial, residual_trial, status_trial, 10.0_dp, 0.5_dp)
    call state%reject()
    if (maxval(abs(state%trial-initial)) > 1.0e-14_dp) error stop 'amplitude rollback'
    if (maxval(abs(state%trial_current_unit-unit_initial)) > 1.0e-14_dp) error stop 'current rollback'
    if (any(state%trial_status /= status_initial)) error stop 'status rollback'
    if (abs(state%target_current-9.0_dp) > 1.0e-14_dp) error stop 'target metadata'
    call state%begin_trial(trial)
    call state%accept()
    if (maxval(abs(state%accepted-trial)) > 1.0e-14_dp) error stop 'amplitude commit'
    if (maxval(abs(state%accepted_current_unit-unit_initial)) > 1.0e-14_dp) error stop 'optional metadata'
    print *, 'periodic amplitude state tests passed'
end program test_periodic_amplitude_state
