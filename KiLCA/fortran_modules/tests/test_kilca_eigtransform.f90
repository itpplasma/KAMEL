program test_kilca_eigtransform
    use iso_fortran_env, only: real64, int32, error_unit
    use kilca_eigtransform_m
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
    print *, "KiLCA Eigenvalue Transformation - Unit Tests"
    print *, "========================================================"
    print *, ""
    
    ! Test 1: Coefficient starting values computation
    call test_coeff_start_vals()
    
    ! Test 2: Coefficients to solution transformation
    call test_coeffs_to_solution()
    
    ! Test 3: Eigenvalue matrix evaluation (placeholder functionality)
    call test_eig_mat_eval()
    
    ! Test 4: Error handling and validation
    call test_eigtransform_error_handling()
    
    ! Test 5: Complex number operations
    call test_complex_matrix_operations()
    
    ! Test 6: Multiple wave/fundamental solution support
    call test_multiple_modes()
    
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
    ! Test 1: Coefficient starting values computation
    !---------------------------------------------------------------------------
    subroutine test_coeff_start_vals()
        integer, parameter :: Nw = 2, Nfs = 2
        real(real64), parameter :: r = 0.1_real64
        real(real64), allocatable :: zstart(:)
        integer :: ierr
        
        call start_test("Coefficient starting values computation")
        
        ! Allocate complex array (2*Nw*Nfs real values)
        allocate(zstart(2*Nw*Nfs))
        
        ! Initialize with test values
        zstart = [1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, &
                  1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64]
        
        ! Test coeff_start_vals function
        call coeff_start_vals(Nw, Nfs, r, zstart, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            ! Check that the solution was modified (not exactly the initial values)
            test_passed = (abs(zstart(1) - 1.0_real64) > 1.0e-10_real64 .or. &
                          abs(zstart(2) - 0.0_real64) > 1.0e-10_real64)
        end if
        
        if (.not. test_passed) then
            print *, "  Coefficient starting values failed, ierr =", ierr
        end if
        
        deallocate(zstart)
        call end_test(test_passed)
    end subroutine test_coeff_start_vals
    
    !---------------------------------------------------------------------------
    ! Test 2: Coefficients to solution transformation  
    !---------------------------------------------------------------------------
    subroutine test_coeffs_to_solution()
        integer, parameter :: Nw = 2, Nfs = 2, dim = 3
        real(real64), allocatable :: rgrid(:), sol(:)
        real(real64) :: initial_sol(2*Nw*Nfs*dim)
        integer :: ierr, i
        
        call start_test("Coefficients to solution transformation")
        
        ! Setup test problem
        allocate(rgrid(dim), sol(2*Nw*Nfs*dim))
        
        ! Create radial grid
        do i = 1, dim
            rgrid(i) = real(i-1, real64) * 0.1_real64
        end do
        
        ! Initialize coefficient matrix (identity-like pattern)
        sol = 0.0_real64
        do i = 1, dim
            ! Set diagonal elements to 1 for each grid point
            sol((i-1)*2*Nw*Nfs + 1) = 1.0_real64  ! (1,1) real part
            sol((i-1)*2*Nw*Nfs + 3) = 1.0_real64  ! (2,1) real part  
            sol((i-1)*2*Nw*Nfs + 5) = 1.0_real64  ! (1,2) real part
            sol((i-1)*2*Nw*Nfs + 7) = 1.0_real64  ! (2,2) real part
        end do
        
        ! Store initial values for comparison
        initial_sol = sol
        
        ! Transform coefficients to solution
        call coeffs_to_solution(Nw, Nfs, dim, rgrid, sol, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            ! Check that solution was modified
            test_passed = (maxval(abs(sol - initial_sol)) > 1.0e-10_real64)
        end if
        
        if (.not. test_passed) then
            print *, "  Coefficients to solution failed, ierr =", ierr
        end if
        
        deallocate(rgrid, sol)
        call end_test(test_passed)
    end subroutine test_coeffs_to_solution
    
    !---------------------------------------------------------------------------
    ! Test 3: Eigenvalue matrix evaluation
    !---------------------------------------------------------------------------
    subroutine test_eig_mat_eval()
        integer, parameter :: Nw = 2
        real(real64), parameter :: r = 0.5_real64
        real(real64), allocatable :: eigmat(:)
        integer :: ierr
        
        call start_test("Eigenvalue matrix evaluation")
        
        ! Allocate complex matrix (2*Nw*Nw real values)
        allocate(eigmat(2*Nw*Nw))
        
        ! Test eigenvalue matrix evaluation
        call eig_mat_eval(r, Nw, eigmat, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            ! Check that matrix contains reasonable values
            test_passed = (maxval(abs(eigmat)) > 1.0e-15_real64)
        end if
        
        if (.not. test_passed) then
            print *, "  Eigenvalue matrix evaluation failed, ierr =", ierr
        end if
        
        deallocate(eigmat)
        call end_test(test_passed)
    end subroutine test_eig_mat_eval
    
    !---------------------------------------------------------------------------
    ! Test 4: Error handling and validation
    !---------------------------------------------------------------------------
    subroutine test_eigtransform_error_handling()
        real(real64), allocatable :: zstart(:), eigmat(:)
        integer :: ierr
        
        call start_test("Eigenvalue transformation error handling")
        
        ! Test with invalid dimensions
        allocate(zstart(4), eigmat(4))
        
        ! Test coeff_start_vals with invalid Nw
        call coeff_start_vals(0, 2, 0.1_real64, zstart, ierr)
        test_passed = (ierr /= 0)  ! Should fail
        
        ! Test eig_mat_eval with invalid Nw
        call eig_mat_eval(0.1_real64, 0, eigmat, ierr)
        test_passed = test_passed .and. (ierr /= 0)  ! Should fail
        
        deallocate(zstart, eigmat)
        call end_test(test_passed)
    end subroutine test_eigtransform_error_handling
    
    !---------------------------------------------------------------------------
    ! Test 5: Complex number operations
    !---------------------------------------------------------------------------
    subroutine test_complex_matrix_operations()
        integer, parameter :: Nw = 2, Nfs = 1
        real(real64), allocatable :: zstart(:)
        real(real64) :: r = 0.2_real64
        integer :: ierr
        
        call start_test("Complex matrix operations")
        
        allocate(zstart(2*Nw*Nfs))
        
        ! Initialize with complex values: z1 = 1+2i, z2 = 3-i
        zstart(1) = 1.0_real64   ! Re(z1)
        zstart(2) = 2.0_real64   ! Im(z1)
        zstart(3) = 3.0_real64   ! Re(z2)
        zstart(4) = -1.0_real64  ! Im(z2)
        
        call coeff_start_vals(Nw, Nfs, r, zstart, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            ! Check that complex arithmetic was performed correctly
            ! (result should be different from input)
            test_passed = (abs(zstart(1) - 1.0_real64) > 1.0e-10_real64 .or. &
                          abs(zstart(2) - 2.0_real64) > 1.0e-10_real64)
        end if
        
        deallocate(zstart)
        call end_test(test_passed)
    end subroutine test_complex_matrix_operations
    
    !---------------------------------------------------------------------------
    ! Test 6: Multiple wave/fundamental solution support
    !---------------------------------------------------------------------------
    subroutine test_multiple_modes()
        integer, parameter :: Nw = 3, Nfs = 3, dim = 2
        real(real64), allocatable :: rgrid(:), sol(:)
        integer :: ierr, i
        
        call start_test("Multiple wave/fundamental solution support")
        
        allocate(rgrid(dim), sol(2*Nw*Nfs*dim))
        
        ! Setup grid
        rgrid = [0.0_real64, 0.5_real64]
        
        ! Initialize with identity-like pattern for multiple modes
        sol = 0.0_real64
        do i = 1, min(Nw, Nfs)
            sol(2*(i-1)*(Nw+1) + 1) = 1.0_real64  ! Diagonal elements
        end do
        
        call coeffs_to_solution(Nw, Nfs, dim, rgrid, sol, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            ! Basic sanity check - solution should contain finite values
            test_passed = all(abs(sol) < 1.0e10_real64) .and. any(abs(sol) > 1.0e-15_real64)
        end if
        
        deallocate(rgrid, sol)
        call end_test(test_passed)
    end subroutine test_multiple_modes
    
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

end program test_kilca_eigtransform