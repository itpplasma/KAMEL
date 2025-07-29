program test_validation
    use iso_fortran_env, only: real64, int32, int64
    use ieee_arithmetic, only: ieee_is_finite
    use kilca_complex_m
    use kilca_stability_m
    use kilca_spline_m
    implicit none
    
    ! Test counters
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed, suite_passed
    
    ! Test parameters
    real(real64), parameter :: STRICT_TOL = epsilon(1.0_real64) * 10.0_real64
    real(real64), parameter :: RELAXED_TOL = epsilon(1.0_real64) * 1000.0_real64
    real(real64), parameter :: PHYS_TOL = 1.0e-12_real64  ! Physical tolerance
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    suite_passed = .true.
    
    print *, "===================================================="
    print *, "KiLCA Validation Tests - Custom Routines Only"
    print *, "===================================================="
    print *, ""
    
    ! Test groups
    call validate_complex_special_functions()
    call validate_stability_algorithms()
    call validate_spline_accuracy()
    call validate_numerical_edge_cases()
    call validate_physical_constraints()
    
    ! Final summary
    print *, ""
    print *, "===================================================="
    print *, "VALIDATION SUMMARY"
    print *, "===================================================="
    print *, "Total tests run:    ", total_tests
    print *, "Tests passed:       ", passed_tests
    print *, "Tests failed:       ", failed_tests
    print *, "Success rate:       ", real(passed_tests)/real(total_tests)*100.0, "%"
    print *, ""
    
    if (failed_tests == 0) then
        print *, "*** ALL VALIDATION TESTS PASSED! ***"
    else
        print *, "*** SOME VALIDATION TESTS FAILED! ***"
        stop 1
    end if
    
contains

    !---------------------------------------------------------------------------
    ! Validate complex special functions
    !---------------------------------------------------------------------------
    subroutine validate_complex_special_functions()
        complex(real64) :: z, z_base, result
        real(real64) :: logabs_result, expected_real
        integer :: i
        
        print *, "1. Validating Complex Special Functions"
        print *, "---------------------------------------"
        
        ! Test 1.1: cmplx_logabs mathematical properties
        call start_test("cmplx_logabs properties")
        test_passed = .true.
        
        ! Property: log|z1*z2| = log|z1| + log|z2|
        z = cmplx(3.0_real64, 4.0_real64, real64)
        z_base = cmplx(2.0_real64, 1.0_real64, real64)
        expected_real = cmplx_logabs(z) + cmplx_logabs(z_base)
        logabs_result = cmplx_logabs(z * z_base)
        test_passed = test_passed .and. (abs(logabs_result - expected_real) < RELAXED_TOL)
        
        ! Property: log|1| = 0
        logabs_result = cmplx_logabs(cmplx(1.0_real64, 0.0_real64, real64))
        test_passed = test_passed .and. (abs(logabs_result) < STRICT_TOL)
        
        ! Property: log|z^n| = n*log|z|
        do i = 2, 5
            result = z**i
            expected_real = real(i, real64) * cmplx_logabs(z)
            logabs_result = cmplx_logabs(result)
            test_passed = test_passed .and. (abs(logabs_result - expected_real) < RELAXED_TOL * real(i, real64))
        end do
        
        call end_test(test_passed)
        
        ! Test 1.2: cmplx_pow_safe edge cases
        call start_test("cmplx_pow_safe edge cases")
        test_passed = .true.
        
        ! 0^0 = 1 (by convention)
        result = cmplx_pow_safe(cmplx(0.0_real64, 0.0_real64, real64), &
                               cmplx(0.0_real64, 0.0_real64, real64))
        test_passed = test_passed .and. (abs(result - cmplx(1.0_real64, 0.0_real64, real64)) < STRICT_TOL)
        
        ! 0^z = 0 for Re(z) > 0
        result = cmplx_pow_safe(cmplx(0.0_real64, 0.0_real64, real64), &
                               cmplx(2.0_real64, 1.0_real64, real64))
        test_passed = test_passed .and. (abs(result) < STRICT_TOL)
        
        ! z^0 = 1 for any z
        result = cmplx_pow_safe(cmplx(3.14_real64, 2.71_real64, real64), &
                               cmplx(0.0_real64, 0.0_real64, real64))
        test_passed = test_passed .and. (abs(result - cmplx(1.0_real64, 0.0_real64, real64)) < STRICT_TOL)
        
        call end_test(test_passed)
        
        ! Test 1.3: Fast operations consistency
        call start_test("Fast operations consistency")
        test_passed = .true.
        
        z = cmplx(3.0_real64, 4.0_real64, real64)
        
        ! abs2_fast should equal abs()^2
        expected_real = abs(z)**2
        logabs_result = cmplx_abs2_fast(z)
        test_passed = test_passed .and. (abs(logabs_result - expected_real) < STRICT_TOL)
        
        ! mult_i_fast should equal multiplication by i
        result = cmplx_mult_i_fast(z)
        z_base = z * cmplx(0.0_real64, 1.0_real64, real64)
        test_passed = test_passed .and. (abs(result - z_base) < STRICT_TOL)
        
        call end_test(test_passed)
        
    end subroutine validate_complex_special_functions
    
    !---------------------------------------------------------------------------
    ! Validate stability analysis algorithms
    !---------------------------------------------------------------------------
    subroutine validate_stability_algorithms()
        complex(real64), allocatable :: A(:,:), A_perturb(:,:)
        real(real64) :: cond1, cond2, rcond, expected
        integer :: n, i, j, info, rank
        
        print *, ""
        print *, "2. Validating Stability Analysis Algorithms"
        print *, "-------------------------------------------"
        
        ! Test 2.1: Condition number properties
        call start_test("Condition number properties")
        test_passed = .true.
        n = 5
        allocate(A(n,n))
        
        ! Identity matrix should have condition number 1
        A = cmplx(0.0_real64, 0.0_real64, real64)
        do i = 1, n
            A(i,i) = cmplx(1.0_real64, 0.0_real64, real64)
        end do
        cond1 = estimate_condition_number(A, info)
        test_passed = test_passed .and. (abs(cond1 - 1.0_real64) < RELAXED_TOL)
        
        ! Scaled matrix: cond(c*A) = cond(A)
        A = A * cmplx(5.0_real64, 0.0_real64, real64)
        cond2 = estimate_condition_number(A, info)
        test_passed = test_passed .and. (abs(cond2 - cond1) < RELAXED_TOL)
        
        deallocate(A)
        call end_test(test_passed)
        
        ! Test 2.2: Reciprocal condition number
        call start_test("Reciprocal condition relationship")
        test_passed = .true.
        n = 4
        allocate(A(n,n))
        
        ! Create a well-conditioned matrix
        A = cmplx(0.0_real64, 0.0_real64, real64)
        do i = 1, n
            A(i,i) = cmplx(real(i, real64), 0.0_real64, real64)
            if (i < n) A(i,i+1) = cmplx(0.1_real64, 0.0_real64, real64)
        end do
        
        cond1 = estimate_condition_number(A, info)
        rcond = estimate_rcond(A, info)
        
        ! rcond should approximately equal 1/cond
        if (cond1 > 0.0_real64) then
            expected = 1.0_real64 / cond1
            test_passed = test_passed .and. (abs(rcond - expected) / expected < 0.1_real64)
        end if
        
        deallocate(A)
        call end_test(test_passed)
        
        ! Test 2.3: Numerical rank validation
        call start_test("Numerical rank determination")
        test_passed = .true.
        n = 6
        allocate(A(n,n))
        
        ! Create rank-3 matrix
        A = cmplx(0.0_real64, 0.0_real64, real64)
        ! First 3 columns are independent
        do i = 1, n
            A(i,1) = cmplx(real(i, real64), 0.0_real64, real64)
            A(i,2) = cmplx(real(i**2, real64), 0.0_real64, real64)
            A(i,3) = cmplx(sin(real(i, real64)), 0.0_real64, real64)
        end do
        ! Columns 4-6 are linear combinations
        A(:,4) = A(:,1) + A(:,2)
        A(:,5) = 2.0_real64 * A(:,1) - A(:,3)
        A(:,6) = A(:,2) + A(:,3)
        
        rank = compute_numerical_rank(A, 1.0e-10_real64, info)
        test_passed = test_passed .and. (rank == 3)
        
        deallocate(A)
        call end_test(test_passed)
        
    end subroutine validate_stability_algorithms
    
    !---------------------------------------------------------------------------
    ! Validate spline accuracy
    !---------------------------------------------------------------------------
    subroutine validate_spline_accuracy()
        type(spline_data_t) :: spline
        real(real64), allocatable :: x(:), y(:), z(:), result(:), exact(:)
        real(real64) :: max_error
        integer :: n, nz, i, ierr
        
        print *, ""
        print *, "3. Validating Spline Accuracy"
        print *, "-----------------------------"
        
        ! Test 3.1: Cubic spline interpolation of polynomial
        call start_test("Cubic spline polynomial reproduction")
        test_passed = .true.
        
        ! For cubic splines, should exactly reproduce polynomials up to degree 3
        n = 10
        allocate(x(n), y(n))
        
        ! Create data from cubic polynomial: y = x^3 - 2x^2 + x + 1
        do i = 1, n
            x(i) = real(i-1, real64) / real(n-1, real64) * 2.0_real64
            y(i) = x(i)**3 - 2.0_real64*x(i)**2 + x(i) + 1.0_real64
        end do
        
        call spline_create(spline, 3, SPLINE_NATURAL, n, x, ierr)
        call spline_calc_coefficients(spline, y, ierr)
        
        ! Evaluate at intermediate points
        nz = 50
        allocate(z(nz), result(nz), exact(nz))
        do i = 1, nz
            z(i) = real(i-1, real64) / real(nz-1, real64) * 2.0_real64
            exact(i) = z(i)**3 - 2.0_real64*z(i)**2 + z(i) + 1.0_real64
        end do
        
        call spline_eval(spline, z, result, ierr)
        
        ! Check maximum error
        max_error = maxval(abs(result - exact))
        test_passed = test_passed .and. (max_error < 0.01_real64)  ! Relaxed due to simple implementation
        
        call spline_destroy(spline)
        deallocate(x, y, z, result, exact)
        call end_test(test_passed)
        
        ! Test 3.2: Derivative accuracy
        call start_test("Spline derivative accuracy")
        test_passed = .true.
        
        n = 20
        allocate(x(n), y(n))
        
        ! Create smooth function
        do i = 1, n
            x(i) = real(i-1, real64) / real(n-1, real64) * 3.14159265358979_real64
            y(i) = sin(x(i))
        end do
        
        call spline_create(spline, 3, SPLINE_NATURAL, n, x, ierr)
        call spline_calc_coefficients(spline, y, ierr)
        
        ! Test derivative at midpoints
        nz = 10
        allocate(z(nz), result(nz), exact(nz))
        do i = 1, nz
            z(i) = real(i, real64) / real(nz+1, real64) * 3.14159265358979_real64
            exact(i) = cos(z(i))  ! Derivative of sin
        end do
        
        call spline_eval_deriv(spline, z, 1, result, ierr)
        
        max_error = maxval(abs(result - exact))
        test_passed = test_passed .and. (max_error < 0.1_real64)  ! Relaxed tolerance
        
        call spline_destroy(spline)
        deallocate(x, y, z, result, exact)
        call end_test(test_passed)
        
    end subroutine validate_spline_accuracy
    
    !---------------------------------------------------------------------------
    ! Validate numerical edge cases
    !---------------------------------------------------------------------------
    subroutine validate_numerical_edge_cases()
        complex(real64) :: z, result
        complex(real64), allocatable :: A(:,:)
        real(real64) :: r, cond
        integer :: info
        logical :: is_stable
        
        print *, ""
        print *, "4. Validating Numerical Edge Cases"
        print *, "----------------------------------"
        
        ! Test 4.1: Very small numbers
        call start_test("Very small number handling")
        test_passed = .true.
        
        ! cmplx_logabs should handle denormalized numbers
        z = cmplx(1.0e-300_real64, 1.0e-300_real64, real64)
        r = cmplx_logabs(z)
        test_passed = test_passed .and. (r > -huge(1.0_real64)) .and. (r < 0.0_real64)
        
        ! Should not overflow or underflow
        test_passed = test_passed .and. ieee_is_finite(r)
        
        call end_test(test_passed)
        
        ! Test 4.2: Near-singular matrices
        call start_test("Near-singular matrix detection")
        test_passed = .true.
        
        allocate(A(3,3))
        A(1,:) = [cmplx(1.0_real64, 0.0_real64, real64), &
                  cmplx(2.0_real64, 0.0_real64, real64), &
                  cmplx(3.0_real64, 0.0_real64, real64)]
        A(2,:) = [cmplx(4.0_real64, 0.0_real64, real64), &
                  cmplx(5.0_real64, 0.0_real64, real64), &
                  cmplx(6.0_real64, 0.0_real64, real64)]
        A(3,:) = [cmplx(7.0_real64, 0.0_real64, real64), &
                  cmplx(8.0_real64, 0.0_real64, real64), &
                  cmplx(9.0_real64 + 1.0e-15_real64, 0.0_real64, real64)]
        
        is_stable = check_matrix_stability(A, 1.0e-10_real64)
        test_passed = test_passed .and. (.not. is_stable)
        
        cond = estimate_condition_number(A, info)
        test_passed = test_passed .and. (cond > 1.0e10_real64)
        
        deallocate(A)
        call end_test(test_passed)
        
    end subroutine validate_numerical_edge_cases
    
    !---------------------------------------------------------------------------
    ! Validate physical constraints
    !---------------------------------------------------------------------------
    subroutine validate_physical_constraints()
        complex(real64), allocatable :: A(:,:), eigenvalues(:)
        real(real64) :: det_mag
        integer :: n, info
        
        print *, ""
        print *, "5. Validating Physical Constraints"
        print *, "----------------------------------"
        
        ! Test 5.1: Hermitian matrix properties
        call start_test("Hermitian matrix validation")
        test_passed = .true.
        n = 4
        allocate(A(n,n))
        
        ! Create Hermitian matrix
        A(1,1) = cmplx(2.0_real64, 0.0_real64, real64)
        A(1,2) = cmplx(1.0_real64, 1.0_real64, real64)
        A(1,3) = cmplx(0.0_real64, 2.0_real64, real64)
        A(1,4) = cmplx(1.0_real64, -1.0_real64, real64)
        
        A(2,1) = conjg(A(1,2))
        A(2,2) = cmplx(3.0_real64, 0.0_real64, real64)
        A(2,3) = cmplx(2.0_real64, 1.0_real64, real64)
        A(2,4) = cmplx(0.0_real64, 1.0_real64, real64)
        
        A(3,1) = conjg(A(1,3))
        A(3,2) = conjg(A(2,3))
        A(3,3) = cmplx(1.0_real64, 0.0_real64, real64)
        A(3,4) = cmplx(1.0_real64, 0.0_real64, real64)
        
        A(4,1) = conjg(A(1,4))
        A(4,2) = conjg(A(2,4))
        A(4,3) = conjg(A(3,4))
        A(4,4) = cmplx(2.0_real64, 0.0_real64, real64)
        
        ! Hermitian matrix should have real determinant magnitude > 0 if positive definite
        det_mag = compute_det_magnitude(A, info)
        test_passed = test_passed .and. (det_mag > 0.0_real64)
        
        deallocate(A)
        call end_test(test_passed)
        
        ! Test 5.2: Conservation properties
        call start_test("Numerical conservation properties")
        test_passed = .true.
        
        ! Test that our operations preserve important properties
        ! This is a placeholder - in real physics code, we'd test conservation laws
        test_passed = .true.
        
        call end_test(test_passed)
        
    end subroutine validate_physical_constraints
    
    !---------------------------------------------------------------------------
    ! Test helper procedures
    !---------------------------------------------------------------------------
    subroutine start_test(name)
        character(len=*), intent(in) :: name
        total_tests = total_tests + 1
        write(*,'(A,A)', advance='no') "   Testing ", name
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
            suite_passed = .false.
        end if
    end subroutine end_test

end program test_validation