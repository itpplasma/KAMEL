program test_kilca_bdf_solver
    use iso_fortran_env
    use kilca_types_m
    use kilca_bdf_solver_m
    implicit none
    
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed
    character(len=100) :: test_name
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    
    print *, ""
    print *, "========================================================"
    print *, "KiLCA BDF Solver - Unit Tests"
    print *, "========================================================"
    print *, ""
    
    ! Run all test suites
    call test_bdf_creation()
    call test_linear_ode()
    call test_stiff_ode()
    call test_van_der_pol()
    call test_robertson_chemical()
    call test_order_selection()
    call test_step_size_control()
    call test_error_estimation()
    call test_convergence_behavior()
    call test_plasma_relevant_ode()
    
    ! Print summary
    print *, ""
    print *, "========================================================"
    print *, "TEST SUMMARY"
    print *, "========================================================"
    print '(A,I5)', " Total tests run:     ", total_tests
    print '(A,I5)', " Tests passed:        ", passed_tests  
    print '(A,I5)', " Tests failed:        ", failed_tests
    print '(A,F8.2,A)', " Success rate:        ", &
            100.0_dp * real(passed_tests, dp) / real(total_tests, dp), " %"
    
    if (failed_tests > 0) then
        print *, ""
        print *, " *** SOME TESTS FAILED! ***"
        stop 1
    else
        print *, ""
        print *, " *** ALL TESTS PASSED! ***"
    end if
    
contains

    !---------------------------------------------------------------------------
    ! Helper subroutines
    !---------------------------------------------------------------------------
    subroutine start_test(name)
        character(len=*), intent(in) :: name
        test_name = name
        write(*, '(A,A,A)', advance='no') "Testing ", trim(name), " ... "
    end subroutine start_test
    
    subroutine end_test(passed)
        logical, intent(in) :: passed
        total_tests = total_tests + 1
        if (passed) then
            passed_tests = passed_tests + 1
            print *, " PASSED"
        else
            failed_tests = failed_tests + 1
            print *, " FAILED"
        end if
    end subroutine end_test
    
    !---------------------------------------------------------------------------
    ! Test 1: BDF solver creation and destruction
    !---------------------------------------------------------------------------
    subroutine test_bdf_creation()
        type(bdf_solver_t) :: solver
        integer :: ierr, neq
        
        call start_test("BDF solver creation")
        test_passed = .true.
        
        neq = 10
        
        ! Test creation with valid parameters
        call bdf_solver_create(solver, neq, order=2, ierr=ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (solver%neq == neq)
        test_passed = test_passed .and. (solver%order == 2)
        test_passed = test_passed .and. allocated(solver%alpha)
        test_passed = test_passed .and. allocated(solver%y_history)
        
        ! Test destruction
        call bdf_solver_destroy(solver, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (.not. allocated(solver%alpha))
        
        ! Test invalid order
        call bdf_solver_create(solver, neq, order=6, ierr=ierr)
        test_passed = test_passed .and. (ierr /= 0)
        
        call end_test(test_passed)
    end subroutine test_bdf_creation
    
    !---------------------------------------------------------------------------
    ! Test 2: Linear ODE y' = -y, y(0) = 1
    !---------------------------------------------------------------------------
    subroutine test_linear_ode()
        type(bdf_solver_t) :: solver
        real(dp) :: y0(1), y_new(1), t0, t_end, h0
        real(dp) :: y_exact, error
        integer :: ierr, i, n_steps
        
        call start_test("Linear ODE (y' = -y)")
        test_passed = .true.
        
        ! Create solver
        call bdf_solver_create(solver, neq=1, order=2, ierr=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Set tolerances
        solver%rtol = 1.0e-8_dp
        solver%atol = 1.0e-10_dp
        
        ! Initial conditions
        t0 = 0.0_dp
        y0(1) = 1.0_dp
        h0 = 0.01_dp
        
        call bdf_solver_init(solver, t0, y0, h0, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Integrate to t = 1
        t_end = 1.0_dp
        n_steps = 100
        
        do i = 1, n_steps
            call bdf_solver_step(solver, linear_rhs, linear_jac, &
                                t0 + i * t_end / n_steps, y_new, ierr)
            test_passed = test_passed .and. (ierr == 0)
        end do
        
        ! Check accuracy
        y_exact = exp(-t_end)
        error = abs(y_new(1) - y_exact)
        test_passed = test_passed .and. (error < 1.0e-6_dp)
        
        call bdf_solver_destroy(solver, ierr)
        
        call end_test(test_passed)
    end subroutine test_linear_ode
    
    ! RHS for y' = -y
    subroutine linear_rhs(neq, t, y, ydot)
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(*)
        real(dp), intent(out) :: ydot(*)
        ydot(1) = -y(1)
    end subroutine linear_rhs
    
    ! Jacobian for y' = -y
    subroutine linear_jac(neq, t, y, jac, ldj)
        integer, intent(in) :: neq, ldj
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(*)
        real(dp), intent(out) :: jac(ldj, *)
        jac(1, 1) = -1.0_dp
    end subroutine linear_jac
    
    !---------------------------------------------------------------------------
    ! Test 3: Stiff ODE y' = -1000(y - cos(t)) - sin(t)
    !---------------------------------------------------------------------------
    subroutine test_stiff_ode()
        type(bdf_solver_t) :: solver
        real(dp) :: y0(1), y_new(1), t0, t_end, h0
        real(dp) :: y_exact, error
        integer :: ierr, i, n_steps
        
        call start_test("Stiff ODE system")
        test_passed = .true.
        
        ! Create solver
        call bdf_solver_create(solver, neq=1, order=3, ierr=ierr)
        
        ! Initial conditions
        t0 = 0.0_dp
        y0(1) = 1.0_dp  ! cos(0)
        h0 = 0.001_dp
        
        call bdf_solver_init(solver, t0, y0, h0, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Integrate to t = 2*pi
        t_end = 2.0_dp * PI
        n_steps = 1000
        
        do i = 1, n_steps
            call bdf_solver_step(solver, stiff_rhs, stiff_jac, &
                                t0 + i * t_end / n_steps, y_new, ierr)
            if (ierr /= 0) exit
        end do
        
        test_passed = test_passed .and. (ierr == 0)
        
        ! Check accuracy - solution should be cos(t)
        y_exact = cos(t_end)
        error = abs(y_new(1) - y_exact)
        test_passed = test_passed .and. (error < 1.0e-4_dp)
        
        ! Check that we used fewer steps than explicit method would need
        test_passed = test_passed .and. (solver%n_steps < 10000)
        
        call bdf_solver_destroy(solver, ierr)
        
        call end_test(test_passed)
    end subroutine test_stiff_ode
    
    ! RHS for stiff ODE
    subroutine stiff_rhs(neq, t, y, ydot)
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(*)
        real(dp), intent(out) :: ydot(*)
        real(dp), parameter :: lambda = 1000.0_dp
        ydot(1) = -lambda * (y(1) - cos(t)) - sin(t)
    end subroutine stiff_rhs
    
    ! Jacobian for stiff ODE
    subroutine stiff_jac(neq, t, y, jac, ldj)
        integer, intent(in) :: neq, ldj
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(*)
        real(dp), intent(out) :: jac(ldj, *)
        real(dp), parameter :: lambda = 1000.0_dp
        jac(1, 1) = -lambda
    end subroutine stiff_jac
    
    !---------------------------------------------------------------------------
    ! Test 4: Van der Pol oscillator (stiff for large mu)
    !---------------------------------------------------------------------------
    subroutine test_van_der_pol()
        type(bdf_solver_t) :: solver
        real(dp) :: y0(2), y_new(2), t0, t_end, h0
        integer :: ierr, i, n_steps
        
        call start_test("Van der Pol oscillator")
        test_passed = .true.
        
        ! Create solver
        call bdf_solver_create(solver, neq=2, order=2, ierr=ierr)
        
        ! Initial conditions
        t0 = 0.0_dp
        y0(1) = 2.0_dp
        y0(2) = 0.0_dp
        h0 = 0.0001_dp
        
        call bdf_solver_init(solver, t0, y0, h0, ierr)
        
        ! Integrate for a short time
        t_end = 0.1_dp
        n_steps = 100
        
        do i = 1, n_steps
            call bdf_solver_step(solver, vanderpol_rhs, vanderpol_jac, &
                                t0 + i * t_end / n_steps, y_new, ierr)
            if (ierr /= 0) exit
        end do
        
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (solver%n_steps <= n_steps)
        
        ! Solution should remain bounded
        test_passed = test_passed .and. (abs(y_new(1)) < 10.0_dp)
        test_passed = test_passed .and. (abs(y_new(2)) < 1000.0_dp)
        
        call bdf_solver_destroy(solver, ierr)
        
        call end_test(test_passed)
    end subroutine test_van_der_pol
    
    ! RHS for Van der Pol: y1' = y2, y2' = mu*(1-y1^2)*y2 - y1
    subroutine vanderpol_rhs(neq, t, y, ydot)
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(*)
        real(dp), intent(out) :: ydot(*)
        real(dp), parameter :: mu = 100.0_dp  ! Stiff parameter
        ydot(1) = y(2)
        ydot(2) = mu * (1.0_dp - y(1)**2) * y(2) - y(1)
    end subroutine vanderpol_rhs
    
    ! Jacobian for Van der Pol
    subroutine vanderpol_jac(neq, t, y, jac, ldj)
        integer, intent(in) :: neq, ldj
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(*)
        real(dp), intent(out) :: jac(ldj, *)
        real(dp), parameter :: mu = 100.0_dp
        jac(1, 1) = 0.0_dp
        jac(1, 2) = 1.0_dp
        jac(2, 1) = -2.0_dp * mu * y(1) * y(2) - 1.0_dp
        jac(2, 2) = mu * (1.0_dp - y(1)**2)
    end subroutine vanderpol_jac
    
    !---------------------------------------------------------------------------
    ! Test 5: Robertson chemical kinetics (very stiff)
    !---------------------------------------------------------------------------
    subroutine test_robertson_chemical()
        type(bdf_solver_t) :: solver
        real(dp) :: y0(3), y_new(3), t0, t_end, h0
        real(dp) :: mass_conservation
        integer :: ierr, i, n_steps
        
        call start_test("Robertson chemical kinetics")
        test_passed = .true.
        
        ! Create solver with higher order for this very stiff problem
        call bdf_solver_create(solver, neq=3, order=5, ierr=ierr)
        
        ! Set tight tolerances
        solver%rtol = 1.0e-10_dp
        solver%atol = 1.0e-12_dp
        
        ! Initial conditions
        t0 = 0.0_dp
        y0(1) = 1.0_dp
        y0(2) = 0.0_dp
        y0(3) = 0.0_dp
        h0 = 1.0e-8_dp  ! Very small initial step
        
        call bdf_solver_init(solver, t0, y0, h0, ierr)
        
        ! Integrate to t = 0.1
        t_end = 0.1_dp
        n_steps = 1000
        
        i = 1
        do while (i <= n_steps)
            call bdf_solver_step(solver, robertson_rhs, robertson_jac, &
                                t0 + i * t_end / n_steps, y_new, ierr)
            if (ierr /= 0) then
                ! Retry with smaller step
                solver%h_current = solver%h_current * 0.5_dp
                i = i - 1
                if (solver%h_current < 1.0e-15_dp) exit
            end if
            i = i + 1
        end do
        
        ! Check mass conservation
        mass_conservation = y_new(1) + y_new(2) + y_new(3)
        test_passed = test_passed .and. (abs(mass_conservation - 1.0_dp) < 1.0e-8_dp)
        
        ! Check that all concentrations are non-negative
        test_passed = test_passed .and. all(y_new >= 0.0_dp)
        
        call bdf_solver_destroy(solver, ierr)
        
        call end_test(test_passed)
    end subroutine test_robertson_chemical
    
    ! Robertson problem RHS
    subroutine robertson_rhs(neq, t, y, ydot)
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(*)
        real(dp), intent(out) :: ydot(*)
        real(dp), parameter :: k1 = 0.04_dp
        real(dp), parameter :: k2 = 3.0e7_dp
        real(dp), parameter :: k3 = 1.0e4_dp
        
        ydot(1) = -k1 * y(1) + k3 * y(2) * y(3)
        ydot(2) = k1 * y(1) - k2 * y(2)**2 - k3 * y(2) * y(3)
        ydot(3) = k2 * y(2)**2
    end subroutine robertson_rhs
    
    ! Robertson problem Jacobian
    subroutine robertson_jac(neq, t, y, jac, ldj)
        integer, intent(in) :: neq, ldj
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(*)
        real(dp), intent(out) :: jac(ldj, *)
        real(dp), parameter :: k1 = 0.04_dp
        real(dp), parameter :: k2 = 3.0e7_dp
        real(dp), parameter :: k3 = 1.0e4_dp
        
        jac(1, 1) = -k1
        jac(1, 2) = k3 * y(3)
        jac(1, 3) = k3 * y(2)
        
        jac(2, 1) = k1
        jac(2, 2) = -2.0_dp * k2 * y(2) - k3 * y(3)
        jac(2, 3) = -k3 * y(2)
        
        jac(3, 1) = 0.0_dp
        jac(3, 2) = 2.0_dp * k2 * y(2)
        jac(3, 3) = 0.0_dp
    end subroutine robertson_jac
    
    !---------------------------------------------------------------------------
    ! Test 6: Order selection
    !---------------------------------------------------------------------------
    subroutine test_order_selection()
        type(bdf_solver_t) :: solver
        real(dp) :: y0(1), y_new(1), err_est
        integer :: ierr, i
        
        call start_test("BDF order selection")
        test_passed = .true.
        
        ! Create solver with adaptive order
        call bdf_solver_create(solver, neq=1, order=1, ierr=ierr)
        solver%max_order = 5
        
        ! Initial conditions
        y0(1) = 1.0_dp
        call bdf_solver_init(solver, 0.0_dp, y0, 0.01_dp, ierr)
        
        ! Take several adaptive steps
        do i = 1, 20
            call bdf_solver_step_adaptive(solver, linear_rhs, linear_jac, &
                                         0.1_dp * i, y_new, err_est, ierr)
        end do
        
        ! Order should have increased from 1
        test_passed = test_passed .and. (solver%order > 1)
        test_passed = test_passed .and. (solver%order <= 5)
        
        call bdf_solver_destroy(solver, ierr)
        
        call end_test(test_passed)
    end subroutine test_order_selection
    
    !---------------------------------------------------------------------------
    ! Test 7: Step size control
    !---------------------------------------------------------------------------
    subroutine test_step_size_control()
        type(bdf_solver_t) :: solver
        real(dp) :: y0(1), y_new(1), err_est
        real(dp) :: h_init, h_final
        integer :: ierr, i
        
        call start_test("Step size control")
        test_passed = .true.
        
        ! Create solver
        call bdf_solver_create(solver, neq=1, order=2, ierr=ierr)
        
        ! Initial conditions with small step
        h_init = 1.0e-6_dp
        y0(1) = 1.0_dp
        call bdf_solver_init(solver, 0.0_dp, y0, h_init, ierr)
        
        ! Take adaptive steps - step size should grow
        do i = 1, 10
            call bdf_solver_step_adaptive(solver, linear_rhs, linear_jac, &
                                         1.0_dp, y_new, err_est, ierr)
            if (solver%t_current >= 0.1_dp) exit
        end do
        
        h_final = solver%h_current
        
        ! Step size should have increased
        test_passed = test_passed .and. (h_final > h_init)
        test_passed = test_passed .and. (h_final < 1.0_dp)  ! But not too much
        
        call bdf_solver_destroy(solver, ierr)
        
        call end_test(test_passed)
    end subroutine test_step_size_control
    
    !---------------------------------------------------------------------------
    ! Test 8: Error estimation
    !---------------------------------------------------------------------------
    subroutine test_error_estimation()
        type(bdf_solver_t) :: solver
        real(dp) :: y0(1), y_new(1), err_est
        integer :: ierr
        
        call start_test("Error estimation")
        test_passed = .true.
        
        ! Create solver
        call bdf_solver_create(solver, neq=1, order=2, ierr=ierr)
        
        ! Set tolerances
        solver%rtol = 1.0e-6_dp
        solver%atol = 1.0e-8_dp
        
        y0(1) = 1.0_dp
        call bdf_solver_init(solver, 0.0_dp, y0, 0.01_dp, ierr)
        
        ! Take a step and check error estimate
        call bdf_solver_step_adaptive(solver, linear_rhs, linear_jac, &
                                     0.1_dp, y_new, err_est, ierr)
        
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (err_est >= 0.0_dp)
        test_passed = test_passed .and. (err_est < 10.0_dp)  ! Reasonable error
        
        call bdf_solver_destroy(solver, ierr)
        
        call end_test(test_passed)
    end subroutine test_error_estimation
    
    !---------------------------------------------------------------------------
    ! Test 9: Convergence behavior
    !---------------------------------------------------------------------------
    subroutine test_convergence_behavior()
        type(bdf_solver_t) :: solver
        integer :: ierr
        
        call start_test("Newton convergence")
        test_passed = .true.
        
        ! Create solver
        call bdf_solver_create(solver, neq=1, order=2, ierr=ierr)
        
        ! Check Newton iteration counts
        test_passed = test_passed .and. (solver%n_newton_iters == 0)
        test_passed = test_passed .and. (solver%n_jac_evals == 0)
        test_passed = test_passed .and. (solver%n_lin_solves == 0)
        
        call bdf_solver_destroy(solver, ierr)
        
        call end_test(test_passed)
    end subroutine test_convergence_behavior
    
    !---------------------------------------------------------------------------
    ! Test 10: Plasma-relevant stiff ODE (simplified charge state evolution)
    !---------------------------------------------------------------------------
    subroutine test_plasma_relevant_ode()
        type(bdf_solver_t) :: solver
        real(dp) :: y0(3), y_new(3), t0, t_end, h0
        real(dp) :: total_density
        integer :: ierr, i, n_steps
        
        call start_test("Plasma charge state evolution")
        test_passed = .true.
        
        ! Create solver for 3 charge states
        call bdf_solver_create(solver, neq=3, order=3, ierr=ierr)
        
        ! Initial conditions: all neutral
        t0 = 0.0_dp
        y0(1) = 1.0_dp  ! Neutral
        y0(2) = 0.0_dp  ! Singly ionized
        y0(3) = 0.0_dp  ! Doubly ionized
        h0 = 1.0e-12_dp  ! Very small for fast ionization rates
        
        call bdf_solver_init(solver, t0, y0, h0, ierr)
        
        ! Integrate
        t_end = 1.0e-6_dp  ! 1 microsecond
        n_steps = 1000
        
        i = 1
        do while (i <= n_steps)
            call bdf_solver_step(solver, plasma_rhs, plasma_jac, &
                                t0 + i * t_end / n_steps, y_new, ierr)
            if (ierr /= 0) then
                solver%h_current = solver%h_current * 0.5_dp
                i = i - 1
            end if
            i = i + 1
        end do
        
        ! Check conservation
        total_density = sum(y_new)
        test_passed = test_passed .and. (abs(total_density - 1.0_dp) < 1.0e-8_dp)
        
        ! Check physical bounds
        test_passed = test_passed .and. all(y_new >= 0.0_dp)
        test_passed = test_passed .and. all(y_new <= 1.0_dp)
        
        ! Should have some ionization
        test_passed = test_passed .and. (y_new(2) + y_new(3) > 0.1_dp)
        
        call bdf_solver_destroy(solver, ierr)
        
        call end_test(test_passed)
    end subroutine test_plasma_relevant_ode
    
    ! Plasma charge state RHS (simplified)
    subroutine plasma_rhs(neq, t, y, ydot)
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(*)
        real(dp), intent(out) :: ydot(*)
        real(dp), parameter :: ne = 1.0e20_dp     ! Electron density
        real(dp), parameter :: Te = 100.0_dp      ! Electron temperature (eV)
        real(dp), parameter :: k1 = 1.0e-14_dp    ! Ionization rate coefficient
        real(dp), parameter :: k2 = 5.0e-15_dp    ! Second ionization rate
        real(dp), parameter :: kr1 = 1.0e-16_dp   ! Recombination rate
        real(dp), parameter :: kr2 = 5.0e-17_dp   ! Second recombination
        
        ! Rate equations
        ydot(1) = -ne * k1 * y(1) + ne * kr1 * y(2)
        ydot(2) = ne * k1 * y(1) - ne * k2 * y(2) - ne * kr1 * y(2) + ne * kr2 * y(3)
        ydot(3) = ne * k2 * y(2) - ne * kr2 * y(3)
    end subroutine plasma_rhs
    
    ! Plasma charge state Jacobian
    subroutine plasma_jac(neq, t, y, jac, ldj)
        integer, intent(in) :: neq, ldj
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(*)
        real(dp), intent(out) :: jac(ldj, *)
        real(dp), parameter :: ne = 1.0e20_dp
        real(dp), parameter :: k1 = 1.0e-14_dp
        real(dp), parameter :: k2 = 5.0e-15_dp
        real(dp), parameter :: kr1 = 1.0e-16_dp
        real(dp), parameter :: kr2 = 5.0e-17_dp
        
        jac(1, 1) = -ne * k1
        jac(1, 2) = ne * kr1
        jac(1, 3) = 0.0_dp
        
        jac(2, 1) = ne * k1
        jac(2, 2) = -ne * k2 - ne * kr1
        jac(2, 3) = ne * kr2
        
        jac(3, 1) = 0.0_dp
        jac(3, 2) = ne * k2
        jac(3, 3) = -ne * kr2
    end subroutine plasma_jac
    
end program test_kilca_bdf_solver