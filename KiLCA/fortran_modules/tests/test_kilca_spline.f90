program test_kilca_spline
    use iso_fortran_env, only: real64, int32, int64
    use kilca_spline_m
    implicit none
    
    ! Test counters
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed
    
    ! Test data
    type(spline_data_t) :: spline
    real(real64), allocatable :: x(:), y(:), z(:), result(:)
    real(real64), allocatable :: coeff(:,:)
    real(real64) :: expected, tol
    integer :: n, i, ierr, degree
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    tol = epsilon(1.0_real64) * 1000.0_real64
    
    print *, "===================================================="
    print *, "KiLCA Spline Module - Unit Tests"
    print *, "===================================================="
    print *, ""
    
    ! Test 1: Create and destroy spline
    call start_test("Spline creation and destruction")
    n = 10
    degree = 3  ! Cubic spline
    allocate(x(n), y(n))
    do i = 1, n
        x(i) = real(i-1, real64)
        y(i) = sin(x(i))
    end do
    
    call spline_create(spline, degree, SPLINE_NATURAL, n, x, ierr)
    test_passed = (ierr == 0) .and. associated(spline%x) .and. (spline%N == degree)
    
    if (test_passed) then
        call spline_destroy(spline)
        test_passed = .not. associated(spline%x)
    end if
    
    deallocate(x, y)
    call end_test(test_passed)
    
    ! Test 2: Calculate spline coefficients
    call start_test("Spline coefficient calculation")
    n = 5
    degree = 3
    allocate(x(n), y(n))
    x = [0.0_real64, 1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    y = [0.0_real64, 1.0_real64, 4.0_real64, 9.0_real64, 16.0_real64]  ! y = x^2
    
    call spline_create(spline, degree, SPLINE_NATURAL, n, x, ierr)
    call spline_calc_coefficients(spline, y, ierr)
    test_passed = (ierr == 0)
    
    deallocate(x, y)
    call spline_destroy(spline)
    call end_test(test_passed)
    
    ! Test 3: Evaluate spline at points
    call start_test("Spline evaluation")
    n = 5
    degree = 3
    allocate(x(n), y(n), z(3), result(3))
    x = [0.0_real64, 1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    y = [0.0_real64, 1.0_real64, 4.0_real64, 9.0_real64, 16.0_real64]
    z = [0.5_real64, 1.5_real64, 2.5_real64]
    
    call spline_create(spline, degree, SPLINE_NATURAL, n, x, ierr)
    call spline_calc_coefficients(spline, y, ierr)
    call spline_eval(spline, z, result, ierr)
    
    ! Check interpolation at known points
    test_passed = (ierr == 0)
    if (test_passed) then
        ! For quadratic function, cubic spline should be quite accurate
        test_passed = test_passed .and. (abs(result(1) - 0.25_real64) < 0.1_real64)
        test_passed = test_passed .and. (abs(result(2) - 2.25_real64) < 0.1_real64)
        test_passed = test_passed .and. (abs(result(3) - 6.25_real64) < 0.1_real64)
    end if
    
    deallocate(x, y, z, result)
    call spline_destroy(spline)
    call end_test(test_passed)
    
    ! Test 4: Evaluate spline derivatives
    call start_test("Spline derivative evaluation")
    n = 5
    degree = 3
    allocate(x(n), y(n), z(3), result(3))
    x = [0.0_real64, 1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    y = [0.0_real64, 1.0_real64, 4.0_real64, 9.0_real64, 16.0_real64]  ! y = x^2, dy/dx = 2x
    z = [1.0_real64, 2.0_real64, 3.0_real64]
    
    call spline_create(spline, degree, SPLINE_NATURAL, n, x, ierr)
    call spline_calc_coefficients(spline, y, ierr)
    call spline_eval_deriv(spline, z, 1, result, ierr)
    
    test_passed = (ierr == 0)
    if (test_passed) then
        ! Derivative of x^2 is 2x
        test_passed = test_passed .and. (abs(result(1) - 2.0_real64) < 0.1_real64)
        test_passed = test_passed .and. (abs(result(2) - 4.0_real64) < 0.1_real64)
        test_passed = test_passed .and. (abs(result(3) - 6.0_real64) < 0.1_real64)
    end if
    
    deallocate(x, y, z, result)
    call spline_destroy(spline)
    call end_test(test_passed)
    
    ! Test 5: Test periodic boundary conditions
    call start_test("Periodic boundary conditions")
    n = 9
    degree = 3
    allocate(x(n), y(n), z(2), result(2))
    do i = 1, n
        x(i) = 2.0_real64 * 3.14159265358979_real64 * real(i-1, real64) / real(n-1, real64)
        y(i) = sin(x(i))
    end do
    z = [0.1_real64, 6.2_real64]  ! Near boundaries
    
    call spline_create(spline, degree, SPLINE_PERIODIC, n, x, ierr)
    call spline_calc_coefficients(spline, y, ierr)
    call spline_eval(spline, z, result, ierr)
    
    test_passed = (ierr == 0)
    if (test_passed) then
        test_passed = test_passed .and. (abs(result(1) - sin(z(1))) < 0.01_real64)
        test_passed = test_passed .and. (abs(result(2) - sin(z(2))) < 0.01_real64)
    end if
    
    deallocate(x, y, z, result)
    call spline_destroy(spline)
    call end_test(test_passed)
    
    ! Test 6: Higher degree splines
    call start_test("Higher degree splines (degree 5)")
    n = 10
    degree = 5
    allocate(x(n), y(n), z(1), result(1))
    do i = 1, n
        x(i) = real(i-1, real64)
        y(i) = exp(-x(i))
    end do
    z = [4.5_real64]
    
    call spline_create(spline, degree, SPLINE_NATURAL, n, x, ierr)
    test_passed = (ierr == 0) .and. (spline%N == degree)
    
    if (test_passed) then
        call spline_calc_coefficients(spline, y, ierr)
        call spline_eval(spline, z, result, ierr)
        test_passed = (ierr == 0) .and. (abs(result(1) - exp(-z(1))) < 0.001_real64)
    end if
    
    deallocate(x, y, z, result)
    call spline_destroy(spline)
    call end_test(test_passed)
    
    ! Test 7: Error handling - invalid degree
    call start_test("Error handling - invalid degree")
    n = 5
    degree = 4  ! Even degree not allowed
    allocate(x(n), y(n))
    x = [0.0_real64, 1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    y = [0.0_real64, 1.0_real64, 4.0_real64, 9.0_real64, 16.0_real64]
    
    call spline_create(spline, degree, SPLINE_NATURAL, n, x, ierr)
    test_passed = (ierr /= 0)  ! Should fail
    
    deallocate(x, y)
    call end_test(test_passed)
    
    ! Test 8: Extrapolation beyond boundaries
    call start_test("Extrapolation handling")
    n = 4
    degree = 3
    allocate(x(n), y(n), z(2), result(2))
    x = [0.0_real64, 1.0_real64, 2.0_real64, 3.0_real64]
    y = [0.0_real64, 1.0_real64, 0.0_real64, -1.0_real64]
    z = [-0.5_real64, 3.5_real64]  ! Outside boundaries
    
    call spline_create(spline, degree, SPLINE_NATURAL, n, x, ierr)
    call spline_calc_coefficients(spline, y, ierr)
    call spline_eval(spline, z, result, ierr)
    
    test_passed = (ierr == 0)  ! Should handle extrapolation gracefully
    
    deallocate(x, y, z, result)
    call spline_destroy(spline)
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

end program test_kilca_spline