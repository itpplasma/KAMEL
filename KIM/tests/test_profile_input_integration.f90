program test_profile_input_integration
    !> Integration tests for profile_input_m: exercise the production
    !> calculate_derivative and compute_er_force_balance routines directly
    !> (no inline reimplementation of the force-balance formula).
    use KIM_kinds_m, only: dp
    use profile_input_m, only: calculate_derivative, compute_er_force_balance
    implicit none

    logical :: all_passed
    all_passed = .true.

    call test_derivative_linear(all_passed)
    call test_er_sign_decreasing(all_passed)
    call test_er_zero_for_flat(all_passed)

    if (all_passed) then
        print *, 'All integration tests PASSED'
        stop 0
    else
        print *, 'Some integration tests FAILED'
        stop 1
    end if

contains

    subroutine test_derivative_linear(passed)
        !> Forward/central/backward differences are exact for a linear profile:
        !> y = 3 r + 1  =>  dy/dr = 3 at every node.
        logical, intent(inout) :: passed
        real(dp) :: r(4), y(4), dydx(4)
        real(dp), parameter :: tol = 1.0e-10_dp

        r = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp]
        y = 3.0_dp * r + 1.0_dp
        call calculate_derivative(r, y, dydx, 4)

        call report('derivative exact for linear profile', &
                    maxval(abs(dydx - 3.0_dp)) < tol, passed)
    end subroutine test_derivative_linear

    subroutine test_er_sign_decreasing(passed)
        !> Decreasing n and Ti with Vz = 0  =>  Er < 0 at interior nodes
        !> (both pressure-gradient terms are negative).
        logical, intent(inout) :: passed
        real(dp) :: r(4), n(4), Ti(4), Vz(4), q(4), Er(4)

        r  = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp]
        n  = [1.0e19_dp, 8.0e18_dp, 5.0e18_dp, 2.0e18_dp]
        Ti = [1800.0_dp, 1600.0_dp, 1000.0_dp, 500.0_dp]
        Vz = 0.0_dp
        q  = 1.0_dp
        call compute_er_force_balance(r, n, Ti, Vz, q, 2.5_dp, 165.0_dp, Er)

        call report('Er < 0 for decreasing profiles (interior node)', &
                    Er(2) < 0.0_dp, passed)
    end subroutine test_er_sign_decreasing

    subroutine test_er_zero_for_flat(passed)
        !> Constant n and Ti with Vz = 0  =>  Er = 0 everywhere.
        logical, intent(inout) :: passed
        real(dp) :: r(4), n(4), Ti(4), Vz(4), q(4), Er(4)
        real(dp), parameter :: tol = 1.0e-12_dp

        r  = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp]
        n  = 5.0e18_dp
        Ti = 1000.0_dp
        Vz = 0.0_dp
        q  = 1.0_dp
        call compute_er_force_balance(r, n, Ti, Vz, q, 2.5_dp, 165.0_dp, Er)

        call report('Er = 0 for flat profiles', maxval(abs(Er)) < tol, passed)
    end subroutine test_er_zero_for_flat

    subroutine report(name, ok, passed)
        character(*), intent(in) :: name
        logical, intent(in) :: ok
        logical, intent(inout) :: passed
        if (ok) then
            print *, 'PASS: ', name
        else
            print *, 'FAIL: ', name
            passed = .false.
        end if
    end subroutine report

end program test_profile_input_integration
