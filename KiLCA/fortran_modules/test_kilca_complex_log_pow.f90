program test_kilca_complex_log_pow
    use iso_fortran_env, only: real64, int32
    use kilca_complex_m
    use kilca_constants_m, only: pi
    implicit none
    
    logical :: test_passed
    complex(real64) :: z1, z2, z3, expected
    real(real64) :: tol
    integer :: n
    
    test_passed = .true.
    tol = epsilon(1.0_real64) * 1000.0_real64
    
    print *, "Testing complex logarithmic and power functions..."
    print *, ""
    
    ! Test 1: Power functions with integer exponents
    print *, "Testing power functions with integer exponents..."
    
    ! Test z^0 = 1
    z1 = cmplx_make(2.0_real64, 3.0_real64)
    z2 = cmplx_pow_int(z1, 0)
    if (abs(z2 - cmplx_E) > tol) then
        print *, "FAIL: z^0 != 1"
        test_passed = .false.
    else
        print *, "PASS: z^0 = 1"
    end if
    
    ! Test z^1 = z
    z2 = cmplx_pow_int(z1, 1)
    if (abs(z2 - z1) > tol) then
        print *, "FAIL: z^1 != z"
        test_passed = .false.
    else
        print *, "PASS: z^1 = z"
    end if
    
    ! Test z^(-1) = 1/z
    z2 = cmplx_pow_int(z1, -1)
    z3 = cmplx_E / z1
    if (abs(z2 - z3) > tol) then
        print *, "FAIL: z^(-1) != 1/z"
        test_passed = .false.
    else
        print *, "PASS: z^(-1) = 1/z"
    end if
    
    ! Test (z^2)^3 = z^6
    z2 = cmplx_pow_int(cmplx_pow_int(z1, 2), 3)
    z3 = cmplx_pow_int(z1, 6)
    if (abs(z2 - z3) > tol * 10.0_real64) then
        print *, "FAIL: (z^2)^3 != z^6"
        test_passed = .false.
    else
        print *, "PASS: (z^2)^3 = z^6"
    end if
    
    ! Test 2: Power functions with real exponents
    print *, ""
    print *, "Testing power functions with real exponents..."
    
    ! Test z^0.5 * z^0.5 = z
    z1 = cmplx_make(3.0_real64, 4.0_real64)
    z2 = cmplx_pow_real(z1, 0.5_real64)
    z3 = z2 * z2
    if (abs(z3 - z1) > tol * 10.0_real64) then
        print *, "FAIL: z^0.5 * z^0.5 != z"
        test_passed = .false.
    else
        print *, "PASS: z^0.5 * z^0.5 = z"
    end if
    
    ! Test i^pi
    z3 = cmplx_pow_real(cmplx_I, pi)
    ! i^pi = exp(pi * log(i)) = exp(pi * i*pi/2) = exp(i*pi^2/2)
    ! |exp(i*theta)| = 1 for any real theta
    if (abs(abs(z3) - 1.0_real64) < tol) then
        print *, "PASS: |i^pi| = 1"
    else
        print *, "FAIL: |i^pi| != 1"
        test_passed = .false.
    end if
    
    ! Test (-1)^0.5 = i or -i
    z1 = cmplx_pow_real(cmplx(-1.0_real64, 0.0_real64, real64), 0.5_real64)
    if (abs(abs(z1) - 1.0_real64) < tol .and. abs(cmplx_real(z1)) < tol) then
        print *, "PASS: (-1)^0.5 has correct form"
    else
        print *, "FAIL: (-1)^0.5 incorrect"
        test_passed = .false.
    end if
    
    ! Test 3: Power functions with complex exponents
    print *, ""
    print *, "Testing power functions with complex exponents..."
    
    ! Test z^(1+0i) = z
    z1 = cmplx_make(2.0_real64, 1.0_real64)
    z2 = cmplx_pow_complex(z1, cmplx_E)
    if (abs(z2 - z1) > tol) then
        print *, "FAIL: z^(1+0i) != z"
        test_passed = .false.
    else
        print *, "PASS: z^(1+0i) = z"
    end if
    
    ! Test z^(a+bi) * z^(c+di) = z^((a+c)+(b+d)i)
    z2 = cmplx_make(0.5_real64, 0.3_real64)
    z3 = cmplx_make(0.7_real64, -0.2_real64)
    expected = cmplx_pow_complex(z1, z2) * cmplx_pow_complex(z1, z3)
    z3 = cmplx_pow_complex(z1, z2 + z3)
    if (abs(expected - z3) / abs(z3) > tol * 100.0_real64) then
        print *, "FAIL: z^w1 * z^w2 != z^(w1+w2)"
        test_passed = .false.
    else
        print *, "PASS: z^w1 * z^w2 = z^(w1+w2)"
    end if
    
    ! Test 4: Logarithm functions
    print *, ""
    print *, "Testing logarithm functions..."
    
    ! Test log(e^z) = z (principal branch)
    z1 = cmplx_make(1.0_real64, 0.5_real64)
    z2 = cmplx_log(cmplx_exp(z1))
    if (abs(z2 - z1) > tol) then
        print *, "FAIL: log(exp(z)) != z"
        test_passed = .false.
    else
        print *, "PASS: log(exp(z)) = z"
    end if
    
    ! Test log(z1*z2) = log(z1) + log(z2) (mod 2πi)
    z1 = cmplx_make(2.0_real64, 1.0_real64)
    z2 = cmplx_make(3.0_real64, -2.0_real64)
    z3 = cmplx_log(z1 * z2)
    expected = cmplx_log(z1) + cmplx_log(z2)
    ! Check if they differ by a multiple of 2πi
    if (abs(real(z3 - expected, real64)) > tol .and. &
        abs(aimag(z3 - expected) - 2.0_real64*pi*nint(aimag(z3 - expected)/(2.0_real64*pi))) > tol) then
        print *, "FAIL: log(z1*z2) != log(z1) + log(z2)"
        test_passed = .false.
    else
        print *, "PASS: log(z1*z2) = log(z1) + log(z2)"
    end if
    
    ! Test log10
    z1 = cmplx(1000.0_real64, 0.0_real64, real64)
    z2 = cmplx_log10(z1)
    expected = cmplx(3.0_real64, 0.0_real64, real64)
    if (abs(z2 - expected) > tol) then
        print *, "FAIL: log10(1000) != 3"
        test_passed = .false.
    else
        print *, "PASS: log10(1000) = 3"
    end if
    
    ! Test log2
    z1 = cmplx(16.0_real64, 0.0_real64, real64)
    z2 = cmplx_log2(z1)
    expected = cmplx(4.0_real64, 0.0_real64, real64)
    if (abs(z2 - expected) > tol) then
        print *, "FAIL: log2(16) != 4"
        test_passed = .false.
    else
        print *, "PASS: log2(16) = 4"
    end if
    
    ! Test 5: Special cases and edge conditions
    print *, ""
    print *, "Testing special cases and edge conditions..."
    
    ! Test 0^0 (should give 1 by convention in many systems)
    z1 = cmplx_pow_int(cmplx_O, 0)
    if (abs(z1 - cmplx_E) > tol) then
        print *, "FAIL: 0^0 != 1"
        test_passed = .false.
    else
        print *, "PASS: 0^0 = 1"
    end if
    
    ! Test log(1) = 0
    z1 = cmplx_log(cmplx_E)
    if (abs(z1) > tol) then
        print *, "FAIL: log(1) != 0"
        test_passed = .false.
    else
        print *, "PASS: log(1) = 0"
    end if
    
    ! Test log(-1) = i*pi
    z1 = cmplx_log(cmplx(-1.0_real64, 0.0_real64, real64))
    expected = cmplx(0.0_real64, pi, real64)
    if (abs(z1 - expected) > tol) then
        print *, "FAIL: log(-1) != i*pi"
        test_passed = .false.
    else
        print *, "PASS: log(-1) = i*pi"
    end if
    
    ! Test log(i) = i*pi/2
    z1 = cmplx_log(cmplx_I)
    expected = cmplx(0.0_real64, pi/2.0_real64, real64)
    if (abs(z1 - expected) > tol) then
        print *, "FAIL: log(i) != i*pi/2"
        test_passed = .false.
    else
        print *, "PASS: log(i) = i*pi/2"
    end if
    
    ! Test 6: Branch cut behavior
    print *, ""
    print *, "Testing branch cut behavior..."
    
    ! Test that log is continuous away from negative real axis
    z1 = cmplx(-1.0_real64, 1.0e-10_real64, real64)
    z2 = cmplx(-1.0_real64, -1.0e-10_real64, real64)
    z3 = cmplx_log(z1) - cmplx_log(z2)
    ! The imaginary parts should differ by approximately 2*pi
    if (abs(aimag(z3) - 2.0_real64*pi) > tol * 1000.0_real64) then
        print *, "WARN: Branch cut discontinuity larger than expected"
        print *, "  Difference in imag part: ", aimag(z3)
    else
        print *, "PASS: Branch cut behavior as expected"
    end if
    
    print *, ""
    if (test_passed) then
        print *, "All logarithmic and power function tests PASSED!"
    else
        print *, "Some logarithmic and power function tests FAILED!"
        stop 1
    end if
    
end program test_kilca_complex_log_pow