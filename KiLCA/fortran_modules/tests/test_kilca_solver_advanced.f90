program test_kilca_solver_advanced
    !---------------------------------------------------------------------------
    ! KiLCA Advanced Solver System - Unit Tests
    !
    ! Tests advanced features of the KiLCA solver system including method 
    ! selection, performance optimization, and integration with RHS functions.
    !
    ! Author: Claude (Anthropic)
    ! Date: 2024
    !
    ! Testing: kilca_solver_m.f90 (advanced features)
    !---------------------------------------------------------------------------
    
    use iso_fortran_env, only: real64, int32, error_unit
    use kilca_solver_m
    use kilca_rhs_func_m
    use kilca_types_m
    implicit none
    
    ! Test counters
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed
    
    ! External procedure declarations
    external :: test_rhs_wrapper, simple_decay_rhs, dummy_rhs, exponential_decay_rhs
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    
    write(*,'(A)') " ========================================================"
    write(*,'(A)') " KiLCA Advanced Solver System - Unit Tests"
    write(*,'(A)') " ========================================================"
    write(*,'(A)') " "
    
    ! Run all tests
    call test_solver_method_selection()
    call test_solver_rhs_integration()
    call test_solver_performance_settings()
    call test_solver_memory_management()
    call test_solver_error_recovery()
    call test_solver_convergence_criteria()
    
    ! Print summary
    write(*,'(A)') " "
    write(*,'(A)') " ========================================================"
    write(*,'(A)') " TEST SUMMARY"
    write(*,'(A)') " ========================================================"
    write(*,'(A,I15)') " Total tests run:              ", total_tests
    write(*,'(A,I15)') " Tests passed:                 ", passed_tests
    write(*,'(A,I15)') " Tests failed:                 ", failed_tests
    write(*,'(A,F15.7,A)') " Success rate:         ", &
        100.0_real64 * real(passed_tests, real64) / real(max(total_tests,1), real64), " %"
    write(*,'(A)') " "
    
    if (failed_tests > 0) then
        write(*,'(A)') " *** SOME TESTS FAILED! ***"
        stop 1
    else
        write(*,'(A)') " *** ALL TESTS PASSED! ***"
    end if

contains

    !---------------------------------------------------------------------------
    ! Test 1: Solver method selection and configuration
    !---------------------------------------------------------------------------
    subroutine test_solver_method_selection()
        type(solver_settings_t) :: settings
        integer :: ierr
        
        call start_test("Solver method selection and configuration")
        test_passed = .true.
        
        ! Test creating settings with different method configurations
        call solver_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            ! Test NORMAL_EXACT method configuration
            call solver_settings_set(settings, 5, 1.0e-8_real64, 1.0e-10_real64, 1000.0_real64, 0, ierr)
            test_passed = test_passed .and. (ierr == 0)
            test_passed = test_passed .and. (settings%Nort == 5)
            test_passed = test_passed .and. (abs(settings%eps_rel - 1.0e-8_real64) < 1.0e-12_real64)
            
            ! Test different accuracy settings
            call solver_settings_set(settings, 10, 1.0e-6_real64, 1.0e-8_real64, 500.0_real64, 1, ierr)
            test_passed = test_passed .and. (ierr == 0)
            test_passed = test_passed .and. (settings%debug == 1)
            
            ! Test boundary values
            call solver_settings_set(settings, 1, 1.0e-12_real64, 1.0e-15_real64, 10.0_real64, 0, ierr)
            test_passed = test_passed .and. (ierr == 0)
        end if
        
        call end_test(test_passed)
    end subroutine test_solver_method_selection
    
    !---------------------------------------------------------------------------
    ! Test 2: Integration with RHS function system
    !---------------------------------------------------------------------------
    subroutine test_solver_rhs_integration()
        type(solver_settings_t) :: solver_settings
        type(rhs_func_params_t) :: rhs_params
        real(real64), allocatable :: rvec(:), y(:)
        integer :: ierr, i
        integer, parameter :: Nfs = 1, Nw = 2, dim = 5
        
        call start_test("Solver integration with RHS function system")
        test_passed = .true.
        
        ! Create solver settings
        call solver_settings_create(solver_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            call solver_settings_set(solver_settings, 3, 1.0e-6_real64, 1.0e-8_real64, 100.0_real64, 0, ierr)
            test_passed = test_passed .and. (ierr == 0)
        end if
        
        ! Create RHS parameters
        call rhs_func_params_create(Nw, 2, Nfs, rhs_params, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            ! Set up integration grid and initial conditions
            allocate(rvec(dim), y(2*Nfs*Nw))
            do i = 1, dim
                rvec(i) = real(i-1, real64) * 0.2_real64  ! [0, 0.2, 0.4, 0.6, 0.8]
            end do
            
            ! Initial conditions: [1.0, 0.0, 0.0, 1.0] for 2 waves, 1 fundamental solution
            y = [1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64]
            
            ! Test ODE integration with RHS function
            call solver_integrate_ode(test_rhs_wrapper, Nfs, Nw, dim, rvec, y, solver_settings, ierr)
            test_passed = test_passed .and. (ierr == 0)
            
            ! Verify solution magnitude is reasonable (not infinity/NaN)
            test_passed = test_passed .and. all(abs(y) < 1000.0_real64)
            test_passed = test_passed .and. all(y == y)  ! NaN check
            
            deallocate(rvec, y)
            call rhs_func_params_destroy(rhs_params, ierr)
        end if
        
        call end_test(test_passed)
    end subroutine test_solver_rhs_integration
    
    !---------------------------------------------------------------------------
    ! Test 3: Performance settings and optimization
    !---------------------------------------------------------------------------
    subroutine test_solver_performance_settings()
        type(solver_settings_t) :: settings_fast, settings_accurate
        real(real64), allocatable :: rvec(:), y1(:), y2(:)
        integer :: ierr, i
        integer, parameter :: dim = 10
        
        call start_test("Solver performance settings and optimization")
        test_passed = .true.
        
        ! Create fast settings (lower accuracy, fewer orthonormalizations)
        call solver_settings_create(settings_fast, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            call solver_settings_set(settings_fast, 2, 1.0e-4_real64, 1.0e-6_real64, 50.0_real64, 0, ierr)
            test_passed = test_passed .and. (ierr == 0)
        end if
        
        ! Create accurate settings (higher accuracy, more orthonormalizations)
        call solver_settings_create(settings_accurate, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            call solver_settings_set(settings_accurate, 10, 1.0e-10_real64, 1.0e-12_real64, 1000.0_real64, 0, ierr)
            test_passed = test_passed .and. (ierr == 0)
        end if
        
        if (ierr == 0) then
            ! Set up test problem
            allocate(rvec(dim), y1(4), y2(4))
            do i = 1, dim
                rvec(i) = real(i-1, real64) * 0.1_real64
            end do
            y1 = [1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64]
            y2 = y1
            
            ! Solve with both settings
            call solver_integrate_ode(simple_decay_rhs, 1, 2, dim, rvec, y1, settings_fast, ierr)
            test_passed = test_passed .and. (ierr == 0)
            
            call solver_integrate_ode(simple_decay_rhs, 1, 2, dim, rvec, y2, settings_accurate, ierr)
            test_passed = test_passed .and. (ierr == 0)
            
            ! Both should give reasonable results
            test_passed = test_passed .and. all(abs(y1) < 10.0_real64)
            test_passed = test_passed .and. all(abs(y2) < 10.0_real64)
            
            deallocate(rvec, y1, y2)
        end if
        
        call end_test(test_passed)
    end subroutine test_solver_performance_settings
    
    !---------------------------------------------------------------------------
    ! Test 4: Memory management and vector operations
    !---------------------------------------------------------------------------
    subroutine test_solver_memory_management()
        real(real64), allocatable :: vectors_in(:,:), vectors_out(:,:)
        real(real64), allocatable :: basis(:,:), coeffs(:), result(:)
        integer :: ierr
        integer, parameter :: n = 6, nvec = 3
        
        call start_test("Solver memory management and vector operations")
        test_passed = .true.
        
        ! Test orthonormalization memory management
        allocate(vectors_in(n, nvec), vectors_out(n, nvec))
        
        ! Create test vectors (linearly independent)
        vectors_in(:,1) = [1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, 0.0_real64, 0.0_real64]
        vectors_in(:,2) = [0.0_real64, 1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, 0.0_real64]  
        vectors_in(:,3) = [0.0_real64, 0.0_real64, 1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64]
        
        call solver_orthonormalize_vectors(n, nvec, vectors_in, vectors_out, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test superposition memory management
        allocate(basis(n, nvec), coeffs(nvec), result(n))
        basis = vectors_out  ! Use orthonormalized vectors as basis
        coeffs = [0.5_real64, 0.3_real64, 0.2_real64]  ! Linear combination coefficients
        
        call solver_superpose_vectors(n, nvec, basis, coeffs, result, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Verify result magnitude is reasonable
        test_passed = test_passed .and. (sqrt(sum(result**2)) > 0.1_real64)
        test_passed = test_passed .and. (sqrt(sum(result**2)) < 10.0_real64)
        
        deallocate(vectors_in, vectors_out, basis, coeffs, result)
        
        call end_test(test_passed)
    end subroutine test_solver_memory_management
    
    !---------------------------------------------------------------------------
    ! Test 5: Error recovery and robustness
    !---------------------------------------------------------------------------
    subroutine test_solver_error_recovery()
        type(solver_settings_t) :: settings
        real(real64), allocatable :: rvec(:), y(:)
        integer :: ierr
        
        call start_test("Solver error recovery and robustness")
        test_passed = .true.
        
        ! Test invalid settings detection
        call solver_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            ! Test invalid Nort
            call solver_settings_set(settings, -1, 1.0e-6_real64, 1.0e-8_real64, 100.0_real64, 0, ierr)
            test_passed = test_passed .and. (ierr /= 0)  ! Should fail
            
            ! Test invalid tolerances
            call solver_settings_set(settings, 5, -1.0e-6_real64, 1.0e-8_real64, 100.0_real64, 0, ierr)
            test_passed = test_passed .and. (ierr /= 0)  ! Should fail
            
            call solver_settings_set(settings, 5, 1.0e-6_real64, -1.0e-8_real64, 100.0_real64, 0, ierr)
            test_passed = test_passed .and. (ierr /= 0)  ! Should fail
        
            ! Test invalid grid
            allocate(rvec(0), y(4))  ! Empty grid
            call solver_integrate_ode(dummy_rhs, 1, 2, 0, rvec, y, settings, ierr)
            test_passed = test_passed .and. (ierr /= 0)  ! Should fail
            
            deallocate(rvec, y)
        end if
        
        call end_test(test_passed)
    end subroutine test_solver_error_recovery
    
    !---------------------------------------------------------------------------
    ! Test 6: Convergence criteria and accuracy
    !---------------------------------------------------------------------------
    subroutine test_solver_convergence_criteria()
        type(solver_settings_t) :: settings
        real(real64), allocatable :: rvec(:), y(:)
        real(real64) :: analytical_solution
        integer :: ierr, i
        integer, parameter :: dim = 11
        
        call start_test("Solver convergence criteria and accuracy")
        test_passed = .true.
        
        call solver_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            ! Use tight tolerances for accuracy test
            call solver_settings_set(settings, 5, 1.0e-10_real64, 1.0e-12_real64, 1000.0_real64, 0, ierr)
            test_passed = test_passed .and. (ierr == 0)
        end if
        
        if (ierr == 0) then
            ! Test exponential decay: dy/dr = -y, y(0) = 1, analytical: y(r) = exp(-r)
            allocate(rvec(dim), y(2))
            do i = 1, dim
                rvec(i) = real(i-1, real64) * 0.1_real64  ! [0, 0.1, ..., 1.0]
            end do
            y = [1.0_real64, 0.0_real64]  ! Initial condition: y = 1 + 0i
            
            call solver_integrate_ode(exponential_decay_rhs, 1, 1, dim, rvec, y, settings, ierr)
            test_passed = test_passed .and. (ierr == 0)
            
            ! Compare with analytical solution at final point
            analytical_solution = exp(-rvec(dim))
            test_passed = test_passed .and. (abs(y(1) - analytical_solution) < 1.0e-6_real64)
            test_passed = test_passed .and. (abs(y(2)) < 1.0e-6_real64)  ! Imaginary part should remain ~0
            
            deallocate(rvec, y)
        end if
        
        call end_test(test_passed)
    end subroutine test_solver_convergence_criteria
    
    !---------------------------------------------------------------------------
    ! Test utilities  
    !---------------------------------------------------------------------------
    subroutine start_test(name)
        character(len=*), intent(in) :: name
        total_tests = total_tests + 1
        write(*,'(A,A)', advance='no') "Testing ", name
        write(*,'(A)', advance='no') " ... "
    end subroutine start_test
    
    subroutine end_test(passed)
        logical, intent(in) :: passed
        if (passed) then
            print *, "PASSED"
            passed_tests = passed_tests + 1
        else
            print *, "FAILED"
            failed_tests = failed_tests + 1
        end if
    end subroutine end_test

end program test_kilca_solver_advanced

!---------------------------------------------------------------------------
! Helper RHS functions for testing
!---------------------------------------------------------------------------
subroutine test_rhs_wrapper(r, y_in, ydot)
    use iso_fortran_env, only: real64
    implicit none
    real(real64), intent(in) :: r
    real(real64), intent(in) :: y_in(*)
    real(real64), intent(out) :: ydot(*)
    
    ! Simple test RHS: dy/dr = -0.1 * y (exponential decay)
    ydot(1) = -0.1_real64 * y_in(1)  ! Real part wave 1
    ydot(2) = -0.1_real64 * y_in(2)  ! Imag part wave 1
    ydot(3) = -0.1_real64 * y_in(3)  ! Real part wave 2
    ydot(4) = -0.1_real64 * y_in(4)  ! Imag part wave 2
end subroutine test_rhs_wrapper

subroutine simple_decay_rhs(r, y, ydot)
    use iso_fortran_env, only: real64
    implicit none
    real(real64), intent(in) :: r, y(*)
    real(real64), intent(out) :: ydot(*)
    ydot(1:4) = -0.2_real64 * y(1:4)
end subroutine simple_decay_rhs

subroutine dummy_rhs(r, y, ydot)
    use iso_fortran_env, only: real64
    implicit none
    real(real64), intent(in) :: r, y(*)
    real(real64), intent(out) :: ydot(*)
    ydot(1:4) = 0.0_real64
end subroutine dummy_rhs

subroutine exponential_decay_rhs(r, y, ydot)
    use iso_fortran_env, only: real64
    implicit none
    real(real64), intent(in) :: r, y(*)
    real(real64), intent(out) :: ydot(*)
    ydot(1) = -y(1)  ! Real part
    ydot(2) = -y(2)  ! Imaginary part
end subroutine exponential_decay_rhs