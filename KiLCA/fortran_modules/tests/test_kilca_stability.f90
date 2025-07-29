program test_kilca_stability
    use iso_fortran_env, only: real64, int32, output_unit, error_unit
    use kilca_stability_m
    implicit none
    
    ! Test counters
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed
    
    ! Test data
    complex(real64), allocatable :: A(:,:), A_test(:,:), b(:), x(:)
    real(real64), allocatable :: A_real(:,:), singular_values(:), work(:)
    real(real64) :: cond_num, rcond, norm_val, det_magnitude
    real(real64) :: tol, expected
    integer :: n, info, rank
    logical :: is_stable, is_singular
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    tol = epsilon(1.0_real64) * 100.0_real64
    
    print *, "===================================================="
    print *, "KiLCA Numerical Stability Module - Unit Tests"
    print *, "===================================================="
    print *, ""
    
    ! Test 1: Condition number estimation for well-conditioned matrix
    call start_test("Condition number - well-conditioned matrix")
    n = 3
    allocate(A(n,n))
    ! Create a well-conditioned matrix
    A = cmplx(0.0_real64, 0.0_real64, real64)
    A(1,1) = cmplx(4.0_real64, 0.0_real64, real64)
    A(1,2) = cmplx(1.0_real64, 0.0_real64, real64)
    A(1,3) = cmplx(0.0_real64, 0.0_real64, real64)
    A(2,1) = cmplx(1.0_real64, 0.0_real64, real64)
    A(2,2) = cmplx(3.0_real64, 0.0_real64, real64)
    A(2,3) = cmplx(1.0_real64, 0.0_real64, real64)
    A(3,1) = cmplx(0.0_real64, 0.0_real64, real64)
    A(3,2) = cmplx(1.0_real64, 0.0_real64, real64)
    A(3,3) = cmplx(2.0_real64, 0.0_real64, real64)
    
    cond_num = estimate_condition_number(A, info)
    test_passed = (info == 0) .and. (cond_num > 1.0_real64) .and. (cond_num < 10.0_real64)
    deallocate(A)
    call end_test(test_passed)
    
    ! Test 2: Condition number for ill-conditioned matrix
    call start_test("Condition number - ill-conditioned matrix")
    n = 3
    allocate(A(n,n))
    ! Create Hilbert matrix (ill-conditioned)
    A(1,1) = cmplx(1.0_real64, 0.0_real64, real64)
    A(1,2) = cmplx(0.5_real64, 0.0_real64, real64)
    A(1,3) = cmplx(1.0_real64/3.0_real64, 0.0_real64, real64)
    A(2,1) = cmplx(0.5_real64, 0.0_real64, real64)
    A(2,2) = cmplx(1.0_real64/3.0_real64, 0.0_real64, real64)
    A(2,3) = cmplx(0.25_real64, 0.0_real64, real64)
    A(3,1) = cmplx(1.0_real64/3.0_real64, 0.0_real64, real64)
    A(3,2) = cmplx(0.25_real64, 0.0_real64, real64)
    A(3,3) = cmplx(0.2_real64, 0.0_real64, real64)
    
    cond_num = estimate_condition_number(A, info)
    test_passed = (info == 0) .and. (cond_num > 100.0_real64)
    deallocate(A)
    call end_test(test_passed)
    
    ! Test 3: Reciprocal condition number
    call start_test("Reciprocal condition number estimation")
    n = 3
    allocate(A(n,n))
    ! Create a matrix with known properties
    A = cmplx(0.0_real64, 0.0_real64, real64)
    A(1,1) = cmplx(2.0_real64, 0.0_real64, real64)
    A(2,2) = cmplx(2.0_real64, 0.0_real64, real64)
    A(3,3) = cmplx(2.0_real64, 0.0_real64, real64)
    
    rcond = estimate_rcond(A, info)
    test_passed = (info == 0) .and. (rcond > 0.0_real64) .and. (rcond <= 1.0_real64)
    deallocate(A)
    call end_test(test_passed)
    
    ! Test 4: Check matrix stability
    call start_test("Matrix stability check")
    n = 3
    allocate(A(n,n))
    ! Well-conditioned matrix
    A = cmplx(0.0_real64, 0.0_real64, real64)
    A(1,1) = cmplx(5.0_real64, 0.0_real64, real64)
    A(2,2) = cmplx(4.0_real64, 0.0_real64, real64)
    A(3,3) = cmplx(3.0_real64, 0.0_real64, real64)
    
    is_stable = check_matrix_stability(A, 1.0e-10_real64)
    test_passed = is_stable
    deallocate(A)
    call end_test(test_passed)
    
    ! Test 5: Singular matrix detection
    call start_test("Singular matrix detection")
    n = 3
    allocate(A(n,n))
    ! Create a singular matrix (third row is sum of first two)
    A(1,:) = [cmplx(1.0_real64, 0.0_real64, real64), &
              cmplx(2.0_real64, 0.0_real64, real64), &
              cmplx(3.0_real64, 0.0_real64, real64)]
    A(2,:) = [cmplx(4.0_real64, 0.0_real64, real64), &
              cmplx(5.0_real64, 0.0_real64, real64), &
              cmplx(6.0_real64, 0.0_real64, real64)]
    A(3,:) = [cmplx(5.0_real64, 0.0_real64, real64), &
              cmplx(7.0_real64, 0.0_real64, real64), &
              cmplx(9.0_real64, 0.0_real64, real64)]
    
    is_singular = is_matrix_singular(A, 1.0e-10_real64)
    test_passed = is_singular
    deallocate(A)
    call end_test(test_passed)
    
    ! Test 6: Numerical rank determination
    call start_test("Numerical rank determination")
    n = 4
    allocate(A(n,n))
    ! Create a rank-2 matrix (two independent columns, rest zero)
    A = cmplx(0.0_real64, 0.0_real64, real64)
    A(1,1) = cmplx(1.0_real64, 0.0_real64, real64)
    A(1,2) = cmplx(0.0_real64, 0.0_real64, real64)
    A(2,1) = cmplx(0.0_real64, 0.0_real64, real64)
    A(2,2) = cmplx(1.0_real64, 0.0_real64, real64)
    A(3,1) = cmplx(1.0_real64, 0.0_real64, real64)
    A(3,2) = cmplx(1.0_real64, 0.0_real64, real64)
    A(4,1) = cmplx(2.0_real64, 0.0_real64, real64)
    A(4,2) = cmplx(3.0_real64, 0.0_real64, real64)
    
    rank = compute_numerical_rank(A, 1.0e-10_real64, info)
    test_passed = (info == 0) .and. (rank == 2)
    if (.not. test_passed) then
        print *, " (Expected rank=2, got rank=", rank, ", info=", info, ")"
    end if
    deallocate(A)
    call end_test(test_passed)
    
    ! Test 7: Matrix norm computation
    call start_test("Matrix norm computation - Frobenius")
    n = 2
    allocate(A(n,n))
    A(1,1) = cmplx(3.0_real64, 0.0_real64, real64)
    A(1,2) = cmplx(4.0_real64, 0.0_real64, real64)
    A(2,1) = cmplx(0.0_real64, 0.0_real64, real64)
    A(2,2) = cmplx(0.0_real64, 0.0_real64, real64)
    
    norm_val = compute_matrix_norm(A, 'F')
    expected = 5.0_real64  ! sqrt(9 + 16)
    test_passed = abs(norm_val - expected) < tol
    deallocate(A)
    call end_test(test_passed)
    
    ! Test 8: Determinant magnitude for stability
    call start_test("Determinant magnitude computation")
    n = 2
    allocate(A(n,n))
    A(1,1) = cmplx(3.0_real64, 0.0_real64, real64)
    A(1,2) = cmplx(1.0_real64, 0.0_real64, real64)
    A(2,1) = cmplx(2.0_real64, 0.0_real64, real64)
    A(2,2) = cmplx(4.0_real64, 0.0_real64, real64)
    
    det_magnitude = compute_det_magnitude(A, info)
    expected = 10.0_real64  ! |3*4 - 1*2| = |12 - 2| = 10
    test_passed = (info == 0) .and. (abs(det_magnitude - expected) < tol)
    deallocate(A)
    call end_test(test_passed)
    
    ! Test 9: System solver stability check
    call start_test("Linear system solver stability")
    n = 3
    allocate(A(n,n), b(n), x(n))
    ! Create a well-conditioned system
    A(1,1) = cmplx(4.0_real64, 0.0_real64, real64)
    A(1,2) = cmplx(1.0_real64, 0.0_real64, real64)
    A(1,3) = cmplx(0.0_real64, 0.0_real64, real64)
    A(2,1) = cmplx(1.0_real64, 0.0_real64, real64)
    A(2,2) = cmplx(3.0_real64, 0.0_real64, real64)
    A(2,3) = cmplx(1.0_real64, 0.0_real64, real64)
    A(3,1) = cmplx(0.0_real64, 0.0_real64, real64)
    A(3,2) = cmplx(1.0_real64, 0.0_real64, real64)
    A(3,3) = cmplx(2.0_real64, 0.0_real64, real64)
    
    b = [cmplx(1.0_real64, 0.0_real64, real64), &
         cmplx(2.0_real64, 0.0_real64, real64), &
         cmplx(3.0_real64, 0.0_real64, real64)]
    
    is_stable = check_system_stability(A, b, 1.0e-10_real64)
    test_passed = is_stable
    deallocate(A, b, x)
    call end_test(test_passed)
    
    ! Test 10: Perturbation analysis
    call start_test("Matrix perturbation sensitivity")
    n = 3
    allocate(A(n,n), A_test(n,n))
    ! Original matrix
    A = cmplx(0.0_real64, 0.0_real64, real64)
    A(1,1) = cmplx(2.0_real64, 0.0_real64, real64)
    A(2,2) = cmplx(2.0_real64, 0.0_real64, real64)
    A(3,3) = cmplx(2.0_real64, 0.0_real64, real64)
    
    ! Perturbed matrix
    A_test = A
    A_test(1,1) = A_test(1,1) + cmplx(1.0e-8_real64, 0.0_real64, real64)
    
    expected = analyze_perturbation_sensitivity(A, A_test, info)
    test_passed = (info == 0) .and. (expected > 0.0_real64) .and. (expected < 1.0_real64)
    deallocate(A, A_test)
    call end_test(test_passed)
    
    ! Final summary
    print *, ""
    print *, "===================================================="
    print *, "TEST SUMMARY"
    print *, "===================================================="
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

end program test_kilca_stability