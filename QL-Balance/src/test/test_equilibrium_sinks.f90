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
    !!     finite difference of the rates, for sources of either sign;
    !!   - for Z_i = 2 the ion pressure source is (2/3) S_w,i = d(n_i T_i)/dt with
    !!     n_i = n_e / Z_i, and converting it back gives the original rates.
    use QLBalance_kinds, only: dp
    use rhs_balance_m, only: compute_equilibrium_source_rates, &
                             linearize_equilibrium_source_rates, &
                             convert_rates_to_equilibrium_sources

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
    call check_linearization([-3.0e13_dp, 5.0_dp, -4.0e4_dp, -2.0e4_dp], 2.0_dp)
    call check_ion_charge()

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
                                              1.0_dp, explicit_rates)
        call compute_equilibrium_source_rates(sources, reference, reference, .true., &
                                              1.0_dp, implicit_rates)
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

        call compute_equilibrium_source_rates(sources, state, reference, .true., 1.0_dp, &
                                              rates)
        call require(all(abs(rates - expected) <= rel_tol * abs(expected) + tiny(1.0_dp)), &
                     'implicit sinks are not scaled with the local state')
    end subroutine check_scaling

    subroutine check_linearization(sources, Z_i)
        real(dp), intent(in) :: sources(4)
        real(dp), intent(in), optional :: Z_i
        real(dp) :: Z
        real(dp) :: state(4), perturbation(4), rates_plus(4), rates_minus(4)
        real(dp) :: finite_difference(4), analytic(4)
        logical :: implicit_sinks
        integer :: imode, ieq
        real(dp), parameter :: step = 1.0e-4_dp

        Z = 1.0_dp
        if (present(Z_i)) Z = Z_i
        state = reference * [1.05_dp, 1.0_dp, 0.97_dp, 1.02_dp]
        do imode = 0, 1
            implicit_sinks = imode == 1
            do ieq = 1, 4
                if (ieq == 2) cycle
                perturbation = 0.0_dp
                perturbation(ieq) = step * state(ieq)
                call compute_equilibrium_source_rates(sources, state + perturbation, &
                                                      reference, implicit_sinks, Z, rates_plus)
                call compute_equilibrium_source_rates(sources, state - perturbation, &
                                                      reference, implicit_sinks, Z, rates_minus)
                finite_difference = 0.5_dp * (rates_plus - rates_minus)
                call linearize_equilibrium_source_rates(sources, state, reference, &
                                                        implicit_sinks, Z, perturbation, analytic)
                call require(all(abs(analytic - finite_difference) <= &
                                 1.0e-5_dp * maxval(abs(finite_difference)) + tiny(1.0_dp)), &
                             'linearised source rates disagree with the finite difference')
            end do
        end do
    end subroutine check_linearization

    subroutine check_ion_charge()
        ! Rates that keep reference stationary, with n_i = n_e / 2.
        real(dp), parameter :: Z = 2.0_dp
        real(dp), parameter :: rates(4) = [1.0e12_dp, 3.0_dp, -2.0e-11_dp, 4.0e-11_dp]
        real(dp) :: sources(4), back(4), expected_ion

        call convert_rates_to_equilibrium_sources(rates, reference, Z, sources)
        expected_ion = (reference(1) * rates(4) + reference(4) * rates(1)) / Z
        call require(abs(sources(4) - expected_ion) <= rel_tol * abs(expected_ion), &
                     'ion pressure source is not d(n_i T_i)/dt for Z_i = 2')
        call require(abs(sources(3) - (reference(1) * rates(3) + reference(3) * rates(1))) <= &
                     rel_tol * abs(reference(1) * rates(3)), &
                     'electron pressure source depends on Z_i')

        call compute_equilibrium_source_rates(sources, reference, reference, .false., Z, back)
        call require(all(abs(back - rates) <= rel_tol * abs(rates)), &
                     'equilibrium sources do not reproduce their rates for Z_i = 2')
    end subroutine check_ion_charge

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, 'FAIL: ', message
            num_failed = num_failed + 1
        end if
    end subroutine require

end program test_equilibrium_sinks
