!> @file test_rk4_integrator.f90
!> @brief Unit tests for the RK4 integrator module
!>
!> Tests the RK4 integrator using a simple harmonic oscillator ODE
!> with known analytical solution: y'' + y = 0
!> Solution: y(t) = cos(t), y'(t) = -sin(t) for y(0) = 1, y'(0) = 0

program test_rk4_integrator
    use rk4_integrator_m, only: rk4_step

    implicit none

    double precision, parameter :: pi = 3.14159265358979d0
    double precision, parameter :: tolerance = 1.0d-6

    logical :: all_passed

    all_passed = .true.

    print *, "========================================"
    print *, "Testing RK4 integrator"
    print *, "========================================"

    call test_harmonic_oscillator(all_passed)
    call test_exponential_decay(all_passed)

    print *, ""
    if (all_passed) then
        print *, "All RK4 integrator tests PASSED"
        stop 0
    else
        print *, "Some RK4 integrator tests FAILED"
        stop 1
    end if

contains

    !> Test RK4 with harmonic oscillator: y'' + y = 0
    !> State: y(1) = position, y(2) = velocity
    !> Analytical solution: y(t) = cos(t), v(t) = -sin(t)
    subroutine test_harmonic_oscillator(passed)
        logical, intent(inout) :: passed

        integer, parameter :: ndim = 2
        integer, parameter :: nsteps = 1000
        double precision, parameter :: t_final = 2.d0 * pi  ! One full period

        double precision :: y(ndim), t, h
        double precision :: y_exact, v_exact, error_y, error_v
        integer :: i

        external :: harmonic_rhs

        print *, ""
        print *, "Test: Harmonic oscillator"

        ! Initial conditions: y(0) = 1, v(0) = 0
        y(1) = 1.d0
        y(2) = 0.d0
        t = 0.d0
        h = t_final / nsteps

        ! Integrate for one full period
        do i = 1, nsteps
            call rk4_step(y, ndim, t, h, harmonic_rhs)
        end do

        ! Check against analytical solution at t = 2*pi
        y_exact = cos(t)
        v_exact = -sin(t)

        error_y = abs(y(1) - y_exact)
        error_v = abs(y(2) - v_exact)

        print *, "  Final t = ", t
        print *, "  y(t):    computed = ", y(1), "  exact = ", y_exact, "  error = ", error_y
        print *, "  v(t):    computed = ", y(2), "  exact = ", v_exact, "  error = ", error_v

        if (error_y < tolerance .and. error_v < tolerance) then
            print *, "  PASSED"
        else
            print *, "  FAILED: errors exceed tolerance ", tolerance
            passed = .false.
        end if

    end subroutine test_harmonic_oscillator

    !> Test RK4 with exponential decay: y' = -y
    !> Analytical solution: y(t) = exp(-t)
    subroutine test_exponential_decay(passed)
        logical, intent(inout) :: passed

        integer, parameter :: ndim = 1
        integer, parameter :: nsteps = 100
        double precision, parameter :: t_final = 5.d0

        double precision :: y(ndim), t, h
        double precision :: y_exact, error
        integer :: i

        external :: decay_rhs

        print *, ""
        print *, "Test: Exponential decay"

        ! Initial condition: y(0) = 1
        y(1) = 1.d0
        t = 0.d0
        h = t_final / nsteps

        ! Integrate
        do i = 1, nsteps
            call rk4_step(y, ndim, t, h, decay_rhs)
        end do

        ! Check against analytical solution
        y_exact = exp(-t)
        error = abs(y(1) - y_exact)

        print *, "  Final t = ", t
        print *, "  y(t):    computed = ", y(1), "  exact = ", y_exact, "  error = ", error

        if (error < tolerance) then
            print *, "  PASSED"
        else
            print *, "  FAILED: error exceeds tolerance ", tolerance
            passed = .false.
        end if

    end subroutine test_exponential_decay

end program test_rk4_integrator


!> RHS for harmonic oscillator: dy/dt = v, dv/dt = -y
subroutine harmonic_rhs(ndim, t, y, dydt)
    implicit none
    integer, intent(in) :: ndim
    double precision, intent(in) :: t, y(ndim)
    double precision, intent(out) :: dydt(ndim)

    dydt(1) = y(2)       ! dy/dt = v
    dydt(2) = -y(1)      ! dv/dt = -y
end subroutine harmonic_rhs


!> RHS for exponential decay: dy/dt = -y
subroutine decay_rhs(ndim, t, y, dydt)
    implicit none
    integer, intent(in) :: ndim
    double precision, intent(in) :: t, y(ndim)
    double precision, intent(out) :: dydt(ndim)

    dydt(1) = -y(1)      ! dy/dt = -y
end subroutine decay_rhs
