program test_equilibrium_sinks
    !! Equilibrium source rates with fixed (explicit) and state-dependent
    !! (implicit) sinks.
    !!
    !! Asserted:
    !!   - at the reference state both treatments give the same rates, so the
    !!     equilibrium stays stationary;
    !!   - away from it, implicit sinks scale negative source densities with
    !!     n / n_0 (particles) and n T / (n_0 T_0) (pressure), positive ones stay;
    !!   - the analytic linearisation used in the Jacobian agrees with a central
    !!     finite difference of the rates, for sources of either sign.
    use QLBalance_kinds, only: dp
    use rhs_balance_m, only: compute_equilibrium_source_rates, &
                             linearize_equilibrium_source_rates

    implicit none

    real(dp), parameter :: reference(4) = [2.0e13_dp, 1.0e4_dp, 1.6e-9_dp, 1.2e-9_dp]
    real(dp), parameter :: rel_tol = 1.0e-6_dp
    integer :: num_failed

    num_failed = 0
    call check_reference_state([3.0e13_dp, 5.0_dp, 4.0e4_dp, 2.0e4_dp])
    call check_reference_state([-3.0e13_dp, 5.0_dp, -4.0e4_dp, -2.0e4_dp])
    call check_scaling()
    call check_linearization([3.0e13_dp, 5.0_dp, 4.0e4_dp, 2.0e4_dp])
    call check_linearization([-3.0e13_dp, 5.0_dp, -4.0e4_dp, 2.0e4_dp])
    call check_linearization([3.0e13_dp, 5.0_dp, 4.0e4_dp, -2.0e4_dp])

    if (num_failed > 0) then
        print *, 'TESTS FAILED: ', num_failed
        error stop 1
    end if
    print *, 'PASS: equilibrium sinks'

contains

    subroutine check_reference_state(sources)
        real(dp), intent(in) :: sources(4)
        real(dp) :: explicit_rates(4), implicit_rates(4)

        call compute_equilibrium_source_rates(sources, reference, reference, .false., &
                                              explicit_rates)
        call compute_equilibrium_source_rates(sources, reference, reference, .true., &
                                              implicit_rates)
        call require(all(abs(implicit_rates - explicit_rates) <= &
                         rel_tol * abs(explicit_rates)), &
                     'implicit sinks change the rates at the reference state')
    end subroutine check_reference_state

    subroutine check_scaling()
        real(dp), parameter :: sources(4) = [-3.0e13_dp, 0.0_dp, -4.0e4_dp, 2.0e4_dp]
        real(dp) :: state(4), rates(4), expected(4), effective(4)

        state = reference * [1.3_dp, 1.0_dp, 0.8_dp, 1.1_dp]
        effective(1) = sources(1) * state(1) / reference(1)
        effective(2) = sources(2)
        effective(3) = sources(3) * state(1) * state(3) / (reference(1) * reference(3))
        effective(4) = sources(4)
        expected(1:2) = effective(1:2)
        expected(3:4) = (effective(3:4) - state(3:4) * effective(1)) / state(1)

        call compute_equilibrium_source_rates(sources, state, reference, .true., rates)
        call require(all(abs(rates - expected) <= rel_tol * abs(expected) + tiny(1.0_dp)), &
                     'implicit sinks are not scaled with the local state')
    end subroutine check_scaling

    subroutine check_linearization(sources)
        real(dp), intent(in) :: sources(4)
        real(dp) :: state(4), perturbation(4), rates_plus(4), rates_minus(4)
        real(dp) :: finite_difference(4), analytic(4)
        logical :: implicit_sinks
        integer :: imode, ieq
        real(dp), parameter :: step = 1.0e-4_dp

        state = reference * [1.05_dp, 1.0_dp, 0.97_dp, 1.02_dp]
        do imode = 0, 1
            implicit_sinks = imode == 1
            do ieq = 1, 4
                if (ieq == 2) cycle
                perturbation = 0.0_dp
                perturbation(ieq) = step * state(ieq)
                call compute_equilibrium_source_rates(sources, state + perturbation, &
                                                      reference, implicit_sinks, rates_plus)
                call compute_equilibrium_source_rates(sources, state - perturbation, &
                                                      reference, implicit_sinks, rates_minus)
                finite_difference = 0.5_dp * (rates_plus - rates_minus)
                call linearize_equilibrium_source_rates(sources, state, reference, &
                                                        implicit_sinks, perturbation, analytic)
                call require(all(abs(analytic - finite_difference) <= &
                                 1.0e-5_dp * maxval(abs(finite_difference)) + tiny(1.0_dp)), &
                             'linearised source rates disagree with the finite difference')
            end do
        end do
    end subroutine check_linearization

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, 'FAIL: ', message
            num_failed = num_failed + 1
        end if
    end subroutine require

end program test_equilibrium_sinks
