program test_br_ramp_residual
    !! Penetration stop of the TimeEvolution ramp (ramp_up_mode = 2: amplitude linear in time,
    !! antenna factor quadratic). br_abs is stored per unit amplitude, so a linear (shielded)
    !! response keeps it constant during the ramp; the residual against the fitted line through
    !! the origin must then stay near zero, and a jump of the response must exceed the 0.1 limit.
    use QLBalance_kinds, only: dp
    use time_evolution, only: br_ramp_residual
    implicit none

    integer, parameter :: n = 60
    real(dp), parameter :: af_exp = 65.64_dp, t_ramp = 1.0_dp, br_unit_lin = 0.1468_dp
    real(dp) :: t(n), af(n), br_unit(n), beta, delta
    integer :: i

    do i = 1, n
        t(i) = 1.3_dp*real(i - 1, dp)/real(n - 1, dp)
    end do
    af = af_exp*(t/t_ramp)**2
    af(1) = 1.0e-4_dp
    br_unit = br_unit_lin

    beta = 0.0_dp
    call br_ramp_residual(br_unit, af, t, beta, delta)
    print *, 'shielded ramp: beta =', beta, ' residual =', delta
    if (delta > 1.0e-3_dp) error stop 'shielded linear ramp reported as departure from linear'
    if (abs(beta - br_unit_lin*sqrt(af_exp)/t_ramp) > 1.0e-2_dp*beta) &
        error stop 'slope differs from the physical ramp rate'

    br_unit(n) = 2.0_dp*br_unit_lin
    call br_ramp_residual(br_unit, af, t, beta, delta)
    print *, 'response doubled at the last step: residual =', delta
    if (delta <= 0.1_dp) error stop 'jump of the response not detected'
end program test_br_ramp_residual
