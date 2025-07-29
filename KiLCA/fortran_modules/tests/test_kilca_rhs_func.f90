program test_kilca_rhs_func
    !---------------------------------------------------------------------------
    ! KiLCA RHS Function System - Unit Tests
    !
    ! Tests the right-hand side function system used by the ODE solver,
    ! including system matrix evaluation and Jacobian computation.
    !
    ! Author: Claude (Anthropic)
    ! Date: 2024
    !
    ! Testing: kilca_rhs_func_m.f90
    !---------------------------------------------------------------------------
    
    use iso_fortran_env, only: real64, int32, int64, error_unit
    use kilca_rhs_func_m
    use kilca_types_m
    implicit none
    
    ! Test counters
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    
    print *, "========================================================"
    print *, "KiLCA RHS Function System - Unit Tests"
    print *, "========================================================"
    print *, ""
    
    ! Test 1: RHS function parameters structure
    call test_rhs_func_params()
    
    ! Test 2: Basic RHS function evaluation
    call test_rhs_func_basic()
    
    ! Test 3: RHS function with different system sizes
    call test_rhs_func_scaling()
    
    ! Test 4: System matrix evaluation (placeholder functionality)
    call test_system_matrix_evaluation()
    
    ! Test 5: Jacobian computation
    call test_jacobian_computation()
    
    ! Test 6: Error handling and validation
    call test_rhs_func_error_handling()
    
    ! Final summary
    print *, ""
    print *, "========================================================"
    print *, "TEST SUMMARY"
    print *, "========================================================"
    print *, "Total tests run:    ", total_tests
    print *, "Tests passed:       ", passed_tests
    print *, "Tests failed:       ", failed_tests
    print *, "Success rate:       ", real(passed_tests)/real(total_tests)*100.0, "%"
    print *, ""
    
    if (failed_tests == 0) then
        print *, "*** ALL TESTS PASSED! ***"
    else
        print *, "*** SOME TESTS FAILED! ***"
        stop 1
    end if
    
contains

    !---------------------------------------------------------------------------
    ! Test 1: RHS function parameters structure
    !---------------------------------------------------------------------------
    subroutine test_rhs_func_params()
        type(rhs_func_params_t) :: params
        integer :: ierr
        
        call start_test("RHS function parameters structure")
        
        ! Test parameter initialization
        call rhs_func_params_create(2, 1, 2, params, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            test_passed = (params%Nwaves == 2)
            test_passed = test_passed .and. (params%Nphys == 1)
            test_passed = test_passed .and. (params%Nfs == 2)
            test_passed = test_passed .and. allocated(params%Dmat)
        end if
        
        if (test_passed) then
            call rhs_func_params_destroy(params, ierr)
            test_passed = (ierr == 0)
        end if
        
        call end_test(test_passed)
    end subroutine test_rhs_func_params
    
    !---------------------------------------------------------------------------
    ! Test 2: Basic RHS function evaluation
    !---------------------------------------------------------------------------
    subroutine test_rhs_func_basic()
        type(rhs_func_params_t) :: params
        real(real64), allocatable :: y(:), ydot(:)
        real(real64) :: r = 0.1_real64
        integer :: Nw, Nfs, neq, ierr
        
        call start_test("Basic RHS function evaluation")
        
        ! Setup test problem
        Nw = 2
        Nfs = 2
        neq = 2 * Nw * Nfs  ! 8 equations
        
        call rhs_func_params_create(Nw, 1, Nfs, params, ierr)
        test_passed = (ierr == 0)
        
        if (test_passed) then
            allocate(y(neq), ydot(neq))
            
            ! Initialize test solution vector
            y = [1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, &
                 1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64]
            ydot = 0.0_real64
            
            ! Evaluate RHS function
            call rhs_func(r, y, ydot, params, ierr)
            
            test_passed = (ierr == 0)
            if (test_passed) then
                ! Check that ydot was modified (not all zeros)
                test_passed = (maxval(abs(ydot)) > 1.0e-15_real64)
            end if
            
            deallocate(y, ydot)
        end if
        
        if (allocated(params%Dmat)) then
            call rhs_func_params_destroy(params, ierr)
        end if
        
        call end_test(test_passed)
    end subroutine test_rhs_func_basic
    
    !---------------------------------------------------------------------------
    ! Test 3: RHS function with different system sizes
    !---------------------------------------------------------------------------
    subroutine test_rhs_func_scaling()
        type(rhs_func_params_t) :: params
        real(real64), allocatable :: y(:), ydot(:)
        real(real64) :: r = 0.2_real64
        integer :: Nw, Nfs, neq, ierr, i
        
        call start_test("RHS function scaling with system size")
        
        ! Test with larger system
        Nw = 3
        Nfs = 3
        neq = 2 * Nw * Nfs  ! 18 equations
        
        call rhs_func_params_create(Nw, 1, Nfs, params, ierr)
        test_passed = (ierr == 0)
        
        if (test_passed) then
            allocate(y(neq), ydot(neq))
            
            ! Initialize with more complex pattern
            do i = 1, neq
                y(i) = sin(real(i, real64) * 0.1_real64)
            end do
            ydot = 0.0_real64
            
            call rhs_func(r, y, ydot, params, ierr)
            
            test_passed = (ierr == 0)
            if (test_passed) then
                ! Check that all components were affected
                test_passed = (maxval(abs(ydot)) > 1.0e-15_real64)
                test_passed = test_passed .and. (count(abs(ydot) > 1.0e-15_real64) > neq/2)
            end if
            
            deallocate(y, ydot)
        end if
        
        if (allocated(params%Dmat)) then
            call rhs_func_params_destroy(params, ierr)
        end if
        
        call end_test(test_passed)
    end subroutine test_rhs_func_scaling
    
    !---------------------------------------------------------------------------
    ! Test 4: System matrix evaluation
    !---------------------------------------------------------------------------
    subroutine test_system_matrix_evaluation()
        real(real64), allocatable :: Dmat(:)
        real(real64) :: r = 0.3_real64
        integer :: Nw = 2, ierr
        
        call start_test("System matrix evaluation")
        
        allocate(Dmat(2*Nw*Nw))  ! Complex matrix storage
        
        ! Test system matrix evaluation
        call eval_diff_sys_matrix(r, 1, Dmat, Nw, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            ! Check that matrix contains reasonable values
            test_passed = (maxval(abs(Dmat)) > 1.0e-15_real64)
            test_passed = test_passed .and. (maxval(abs(Dmat)) < 1.0e10_real64)
        end if
        
        deallocate(Dmat)
        call end_test(test_passed)
    end subroutine test_system_matrix_evaluation
    
    !---------------------------------------------------------------------------
    ! Test 5: Jacobian computation
    !---------------------------------------------------------------------------
    subroutine test_jacobian_computation()
        type(rhs_func_params_t) :: params
        real(real64), allocatable :: y(:), fy(:), jac(:,:)
        real(real64) :: r = 0.15_real64
        integer :: Nw, Nfs, neq, ierr
        
        call start_test("Jacobian computation")
        
        ! Setup test problem
        Nw = 2
        Nfs = 2
        neq = 2 * Nw * Nfs
        
        call rhs_func_params_create(Nw, 1, Nfs, params, ierr)
        test_passed = (ierr == 0)
        
        if (test_passed) then
            allocate(y(neq), fy(neq), jac(neq, neq))
            
            ! Initialize vectors
            y = [1.0_real64, 0.5_real64, -0.5_real64, 1.0_real64, &
                 0.0_real64, 1.0_real64, 1.0_real64, 0.0_real64]
            fy = 0.0_real64
            jac = 0.0_real64
            
            ! Compute Jacobian
            call jacobian_func(r, y, fy, jac, params, ierr)
            
            test_passed = (ierr == 0)
            if (test_passed) then
                ! Check that Jacobian contains non-zero values
                test_passed = (maxval(abs(jac)) > 1.0e-15_real64)
                ! Check matrix structure (should be sparse but not empty)
                test_passed = test_passed .and. (count(abs(jac) > 1.0e-15_real64) > 0)
            end if
            
            deallocate(y, fy, jac)
        end if
        
        if (allocated(params%Dmat)) then
            call rhs_func_params_destroy(params, ierr)
        end if
        
        call end_test(test_passed)
    end subroutine test_jacobian_computation
    
    !---------------------------------------------------------------------------
    ! Test 6: Error handling and validation
    !---------------------------------------------------------------------------
    subroutine test_rhs_func_error_handling()
        type(rhs_func_params_t) :: params
        real(real64), allocatable :: y(:), ydot(:)
        integer :: ierr
        
        call start_test("RHS function error handling")
        
        ! Test with invalid dimensions
        call rhs_func_params_create(0, 1, 1, params, ierr)
        test_passed = (ierr /= 0)  ! Should fail
        
        ! Test with valid parameters but invalid inputs
        call rhs_func_params_create(1, 1, 1, params, ierr)
        if (ierr == 0) then
            allocate(y(2), ydot(2))
            y = [1.0_real64, 0.0_real64]
            ydot = 0.0_real64
            
            ! Test with NaN input
            y(1) = transfer(int(z'7FF8000000000000', int64), 1.0_real64)  ! NaN
            call rhs_func(0.1_real64, y, ydot, params, ierr)
            test_passed = test_passed .and. (ierr /= 0 .or. .not. any(ydot /= ydot))  ! Should handle NaN
            
            deallocate(y, ydot)
            call rhs_func_params_destroy(params, ierr)
        end if
        
        call end_test(test_passed)
    end subroutine test_rhs_func_error_handling
    
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

end program test_kilca_rhs_func