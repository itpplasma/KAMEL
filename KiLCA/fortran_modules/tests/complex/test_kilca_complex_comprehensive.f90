program test_kilca_complex_comprehensive
    use iso_fortran_env, only: real64, int32, output_unit, error_unit
    use kilca_complex_m
    use kilca_constants_m, only: pi
    implicit none
    
    ! Test counters
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed, suite_passed
    
    ! Test data
    complex(real64) :: z1, z2, z3, result, expected
    complex(real64), allocatable :: array1(:), array2(:), result_array(:)
    complex(real64), allocatable :: matrix(:,:)
    real(real64) :: tol, r, theta
    integer :: i, n
    character(len=256) :: test_name
    character(len=100) :: str
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    suite_passed = .true.
    tol = epsilon(1.0_real64) * 100.0_real64
    
    print *, "===================================================="
    print *, "KiLCA Complex Number Module - Comprehensive Unit Tests"
    print *, "===================================================="
    print *, ""
    
    ! ============= Test Group 1: Basic Creation and Access =============
    call start_test_group("Basic Creation and Access Functions")
    
    ! Test 1.1: Complex number creation
    call start_test("cmplx_make")
    z1 = cmplx_make(3.0_real64, 4.0_real64)
    test_passed = (real(z1) == 3.0_real64) .and. (aimag(z1) == 4.0_real64)
    call end_test(test_passed)
    
    ! Test 1.2: Polar creation
    call start_test("cmplx_polar")
    z1 = cmplx_polar(5.0_real64, atan2(4.0_real64, 3.0_real64))
    test_passed = abs(real(z1) - 3.0_real64) < tol .and. abs(aimag(z1) - 4.0_real64) < tol
    call end_test(test_passed)
    
    ! Test 1.3: Real and imaginary part extraction
    call start_test("cmplx_real and cmplx_imag")
    z1 = cmplx(3.0_real64, 4.0_real64, real64)
    test_passed = cmplx_real(z1) == 3.0_real64 .and. cmplx_imag(z1) == 4.0_real64
    call end_test(test_passed)
    
    ! Test 1.4: Magnitude and argument
    call start_test("cmplx_abs and cmplx_arg")
    z1 = cmplx(3.0_real64, 4.0_real64, real64)
    test_passed = abs(cmplx_abs(z1) - 5.0_real64) < tol .and. &
                  abs(cmplx_arg(z1) - atan2(4.0_real64, 3.0_real64)) < tol
    call end_test(test_passed)
    
    ! ============= Test Group 2: Basic Arithmetic =============
    call start_test_group("Basic Arithmetic Operations")
    
    ! Test 2.1: Addition
    call start_test("cmplx_add")
    z1 = cmplx(1.0_real64, 2.0_real64, real64)
    z2 = cmplx(3.0_real64, 4.0_real64, real64)
    result = cmplx_add(z1, z2)
    expected = cmplx(4.0_real64, 6.0_real64, real64)
    test_passed = abs(result - expected) < tol
    call end_test(test_passed)
    
    ! Test 2.2: Subtraction
    call start_test("cmplx_sub")
    result = cmplx_sub(z1, z2)
    expected = cmplx(-2.0_real64, -2.0_real64, real64)
    test_passed = abs(result - expected) < tol
    call end_test(test_passed)
    
    ! Test 2.3: Multiplication
    call start_test("cmplx_mult")
    result = cmplx_mult(z1, z2)
    expected = cmplx(-5.0_real64, 10.0_real64, real64)
    test_passed = abs(result - expected) < tol
    call end_test(test_passed)
    
    ! Test 2.4: Division
    call start_test("cmplx_div")
    result = cmplx_div(z2, z1)
    expected = cmplx(2.2_real64, -0.4_real64, real64)
    test_passed = abs(result - expected) < tol
    call end_test(test_passed)
    
    ! ============= Test Group 3: Transcendental Functions =============
    call start_test_group("Transcendental Functions")
    
    ! Test 3.1: Exponential
    call start_test("cmplx_exp")
    z1 = cmplx(1.0_real64, pi/2.0_real64, real64)
    result = cmplx_exp(z1)
    expected = cmplx(0.0_real64, exp(1.0_real64), real64)
    test_passed = abs(result - expected) < tol * 10.0_real64
    call end_test(test_passed)
    
    ! Test 3.2: Logarithm
    call start_test("cmplx_log")
    z1 = cmplx(1.0_real64, 1.0_real64, real64)
    result = cmplx_log(z1)
    expected = cmplx(log(sqrt(2.0_real64)), pi/4.0_real64, real64)
    test_passed = abs(result - expected) < tol
    call end_test(test_passed)
    
    ! Test 3.3: Square root
    call start_test("cmplx_sqrt")
    z1 = cmplx(-1.0_real64, 0.0_real64, real64)
    result = cmplx_sqrt(z1)
    expected = cmplx(0.0_real64, 1.0_real64, real64)
    test_passed = abs(result - expected) < tol
    call end_test(test_passed)
    
    ! Test 3.4: Trigonometric functions
    call start_test("cmplx_sin, cmplx_cos, cmplx_tan")
    z1 = cmplx(0.0_real64, 0.0_real64, real64)
    test_passed = abs(cmplx_sin(z1)) < tol .and. &
                  abs(cmplx_cos(z1) - cmplx_E) < tol .and. &
                  abs(cmplx_tan(z1)) < tol
    call end_test(test_passed)
    
    ! ============= Test Group 4: Power Functions =============
    call start_test_group("Power Functions")
    
    ! Test 4.1: Integer power
    call start_test("cmplx_pow_int")
    z1 = cmplx(2.0_real64, 1.0_real64, real64)
    result = cmplx_pow_int(z1, 3)
    expected = cmplx(2.0_real64, 11.0_real64, real64)
    test_passed = abs(result - expected) < tol
    call end_test(test_passed)
    
    ! Test 4.2: Real power
    call start_test("cmplx_pow_real")
    z1 = cmplx(4.0_real64, 0.0_real64, real64)
    result = cmplx_pow_real(z1, 0.5_real64)
    expected = cmplx(2.0_real64, 0.0_real64, real64)
    test_passed = abs(result - expected) < tol
    call end_test(test_passed)
    
    ! Test 4.3: Complex power
    call start_test("cmplx_pow_complex")
    z1 = cmplx(exp(1.0_real64), 0.0_real64, real64)
    z2 = cmplx(0.0_real64, pi, real64)
    result = cmplx_pow_complex(z1, z2)
    expected = cmplx(-1.0_real64, 0.0_real64, real64)
    test_passed = abs(result - expected) < tol * 10.0_real64
    call end_test(test_passed)
    
    ! ============= Test Group 5: Special Functions =============
    call start_test_group("Special Functions")
    
    ! Test 5.1: Safe logarithm of absolute value
    call start_test("cmplx_logabs")
    z1 = cmplx(1.0e-200_real64, 1.0e-200_real64, real64)
    r = cmplx_logabs(z1)
    expected = log(sqrt(2.0_real64)) - 200.0_real64 * log(10.0_real64)
    test_passed = abs(r - expected) < tol * abs(expected)
    call end_test(test_passed)
    
    ! Test 5.2: Safe power function
    call start_test("cmplx_pow_safe")
    z1 = cmplx(0.0_real64, 0.0_real64, real64)
    z2 = cmplx(0.0_real64, 0.0_real64, real64)
    result = cmplx_pow_safe(z1, z2)
    expected = cmplx(1.0_real64, 0.0_real64, real64)
    test_passed = abs(result - expected) < tol
    call end_test(test_passed)
    
    ! ============= Test Group 6: Validation Functions =============
    call start_test_group("Validation and Comparison Functions")
    
    ! Test 6.1: Finite check
    call start_test("cmplx_is_finite")
    z1 = cmplx(1.0_real64, 2.0_real64, real64)
    test_passed = cmplx_is_finite(z1)
    call end_test(test_passed)
    
    ! Test 6.2: Zero check
    call start_test("cmplx_is_zero")
    z1 = cmplx(0.0_real64, 0.0_real64, real64)
    z2 = cmplx(1.0e-10_real64, 0.0_real64, real64)
    ! Test that exact zero is detected and small value is not (with default tolerance)
    test_passed = cmplx_is_zero(z1) .and. .not. cmplx_is_zero(z2)
    call end_test(test_passed)
    
    ! Test 6.3: Real/Imaginary checks
    call start_test("cmplx_is_real and cmplx_is_imaginary")
    z1 = cmplx(3.0_real64, 0.0_real64, real64)
    z2 = cmplx(0.0_real64, 4.0_real64, real64)
    test_passed = cmplx_is_real(z1) .and. cmplx_is_imaginary(z2)
    call end_test(test_passed)
    
    ! Test 6.4: Approximate equality
    call start_test("cmplx_equals")
    z1 = cmplx(1.0_real64, 2.0_real64, real64)
    z2 = cmplx(1.0_real64 + 1.0e-15_real64, 2.0_real64, real64)
    test_passed = cmplx_equals(z1, z2, 1.0e-14_real64)
    call end_test(test_passed)
    
    ! ============= Test Group 7: Performance Optimizations =============
    call start_test_group("Performance Optimization Functions")
    
    ! Test 7.1: Fast operations
    call start_test("Fast arithmetic operations")
    z1 = cmplx(3.0_real64, 4.0_real64, real64)
    r = cmplx_abs2_fast(z1)
    test_passed = abs(r - 25.0_real64) < tol
    result = cmplx_div_real_fast(z1, 2.0_real64)
    test_passed = test_passed .and. abs(result - cmplx(1.5_real64, 2.0_real64, real64)) < tol
    result = cmplx_mult_i_fast(z1)
    test_passed = test_passed .and. abs(result - cmplx(-4.0_real64, 3.0_real64, real64)) < tol
    call end_test(test_passed)
    
    ! Test 7.2: In-place operations
    call start_test("In-place operations")
    z1 = cmplx(1.0_real64, 2.0_real64, real64)
    z2 = cmplx(3.0_real64, 4.0_real64, real64)
    call cmplx_add_inplace(z1, z2)
    test_passed = abs(z1 - cmplx(4.0_real64, 6.0_real64, real64)) < tol
    call end_test(test_passed)
    
    ! ============= Test Group 8: Coordinate Transformations =============
    call start_test_group("Coordinate Transformation Functions")
    
    ! Test 8.1: Polar conversions
    call start_test("Polar coordinate conversions")
    z1 = cmplx(3.0_real64, 4.0_real64, real64)
    call cmplx_to_polar(z1, r, theta)
    test_passed = abs(r - 5.0_real64) < tol
    z2 = cmplx_from_polar(r, theta)
    test_passed = test_passed .and. abs(z2 - z1) < tol
    call end_test(test_passed)
    
    ! Test 8.2: Complex rotation
    call start_test("Complex number rotation")
    z1 = cmplx(1.0_real64, 0.0_real64, real64)
    z2 = cmplx_rotate(z1, pi/2.0_real64)
    expected = cmplx(0.0_real64, 1.0_real64, real64)
    test_passed = abs(z2 - expected) < tol
    call end_test(test_passed)
    
    ! ============= Test Group 9: Array Operations =============
    call start_test_group("Array Operations")
    
    ! Test 9.1: Array allocation
    call start_test("Array allocation and deallocation")
    n = 10
    call cmplx_allocate_1d(array1, n)
    test_passed = allocated(array1) .and. size(array1) == n
    call cmplx_deallocate_1d(array1)
    test_passed = test_passed .and. .not. allocated(array1)
    call end_test(test_passed)
    
    ! Test 9.2: Array operations
    call start_test("Array arithmetic operations")
    n = 5
    allocate(array1(n), array2(n), result_array(n))
    array1 = [(cmplx(real(i, real64), 0.0_real64, real64), i=1,n)]
    array2 = [(cmplx(0.0_real64, real(i, real64), real64), i=1,n)]
    call cmplx_array_add_fast(array1, array2, result_array)
    test_passed = .true.
    do i = 1, n
        expected = cmplx(real(i, real64), real(i, real64), real64)
        if (abs(result_array(i) - expected) > tol) test_passed = .false.
    end do
    deallocate(array1, array2, result_array)
    call end_test(test_passed)
    
    ! ============= Test Group 10: Matrix Operations =============
    call start_test_group("Matrix Operations")
    
    ! Test 10.1: Matrix trace
    call start_test("Matrix trace calculation")
    allocate(matrix(3,3))
    matrix = cmplx(0.0_real64, 0.0_real64, real64)
    matrix(1,1) = cmplx(1.0_real64, 1.0_real64, real64)
    matrix(2,2) = cmplx(2.0_real64, 2.0_real64, real64)
    matrix(3,3) = cmplx(3.0_real64, 3.0_real64, real64)
    result = cmplx_matrix_trace(matrix)
    expected = cmplx(6.0_real64, 6.0_real64, real64)
    test_passed = abs(result - expected) < tol
    deallocate(matrix)
    call end_test(test_passed)
    
    ! ============= Test Group 11: I/O Functions =============
    call start_test_group("Input/Output Functions")
    
    ! Test 11.1: String conversion
    call start_test("String conversion")
    z1 = cmplx(3.14_real64, -2.71_real64, real64)
    str = cmplx_to_string(z1)
    test_passed = len_trim(str) > 0  ! Basic check that string is not empty
    call end_test(test_passed)
    
    ! Test 11.2: String parsing
    call start_test("String parsing")
    str = "(1.5,2.5)"
    z1 = cmplx_parse_string(str)
    expected = cmplx(1.5_real64, 2.5_real64, real64)
    test_passed = abs(z1 - expected) < tol
    call end_test(test_passed)
    
    ! ============= Final Summary =============
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

    subroutine start_test_group(group_name)
        character(len=*), intent(in) :: group_name
        print *, ""
        print *, "---------- ", trim(group_name), " ----------"
    end subroutine start_test_group
    
    subroutine start_test(name)
        character(len=*), intent(in) :: name
        test_name = name
        total_tests = total_tests + 1
    end subroutine start_test
    
    subroutine end_test(passed)
        logical, intent(in) :: passed
        if (passed) then
            print *, "[PASS] ", trim(test_name)
            passed_tests = passed_tests + 1
        else
            print *, "[FAIL] ", trim(test_name)
            failed_tests = failed_tests + 1
            suite_passed = .false.
        end if
    end subroutine end_test
    
end program test_kilca_complex_comprehensive