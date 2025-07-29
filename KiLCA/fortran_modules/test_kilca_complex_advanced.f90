program test_kilca_complex_advanced
    use iso_fortran_env, only: real64, int32
    use kilca_complex_m
    use kilca_constants_m, only: pi
    implicit none
    
    logical :: test_passed
    complex(real64) :: z1, z2, z3, expected
    real(real64) :: tol, r1, r2
    
    test_passed = .true.
    tol = epsilon(1.0_real64) * 1000.0_real64
    
    print *, "Testing advanced complex transcendental functions..."
    print *, ""
    
    ! Test 1: Inverse hyperbolic functions
    print *, "Testing inverse hyperbolic functions..."
    
    ! Test sinh(asinh(z)) = z
    z1 = cmplx_make(1.5_real64, 2.0_real64)
    z2 = cmplx_sinh(cmplx_asinh(z1))
    if (abs(z2 - z1) > tol) then
        print *, "FAIL: sinh(asinh(z)) != z"
        print *, "  z = ", z1
        print *, "  Result = ", z2
        test_passed = .false.
    else
        print *, "PASS: sinh(asinh(z)) = z"
    end if
    
    ! Test cosh(acosh(z)) = z
    z1 = cmplx_make(2.0_real64, 1.0_real64)
    z2 = cmplx_cosh(cmplx_acosh(z1))
    if (abs(z2 - z1) > tol) then
        print *, "FAIL: cosh(acosh(z)) != z"
        test_passed = .false.
    else
        print *, "PASS: cosh(acosh(z)) = z"
    end if
    
    ! Test tanh(atanh(z)) = z for |z| < 1
    z1 = cmplx_make(0.5_real64, 0.3_real64)
    z2 = cmplx_tanh(cmplx_atanh(z1))
    if (abs(z2 - z1) > tol) then
        print *, "FAIL: tanh(atanh(z)) != z"
        test_passed = .false.
    else
        print *, "PASS: tanh(atanh(z)) = z"
    end if
    
    ! Test 2: Base-2 functions
    print *, ""
    print *, "Testing base-2 functions..."
    
    ! Test exp2(log2(z)) = z
    z1 = cmplx_make(3.0_real64, 4.0_real64)
    z2 = cmplx_exp2(cmplx_log2(z1))
    if (abs(z2 - z1) > tol) then
        print *, "FAIL: exp2(log2(z)) != z"
        test_passed = .false.
    else
        print *, "PASS: exp2(log2(z)) = z"
    end if
    
    ! Test exp2(3) = 8
    z1 = cmplx(3.0_real64, 0.0_real64, real64)
    z2 = cmplx_exp2(z1)
    expected = cmplx(8.0_real64, 0.0_real64, real64)
    if (abs(z2 - expected) > tol) then
        print *, "FAIL: exp2(3) != 8"
        test_passed = .false.
    else
        print *, "PASS: exp2(3) = 8"
    end if
    
    ! Test log2(8) = 3
    z1 = cmplx(8.0_real64, 0.0_real64, real64)
    z2 = cmplx_log2(z1)
    expected = cmplx(3.0_real64, 0.0_real64, real64)
    if (abs(z2 - expected) > tol) then
        print *, "FAIL: log2(8) != 3"
        test_passed = .false.
    else
        print *, "PASS: log2(8) = 3"
    end if
    
    ! Test 3: Functions for small arguments
    print *, ""
    print *, "Testing functions for small arguments..."
    
    ! Test expm1 for small z
    z1 = cmplx(1.0e-8_real64, 1.0e-8_real64, real64)
    z2 = cmplx_expm1(z1)
    z3 = cmplx_exp(z1) - cmplx_E
    if (abs(z2 - z1) > tol * 10.0_real64) then  ! For small z, expm1(z) ≈ z
        print *, "FAIL: expm1(small z) != z"
        print *, "  z = ", z1
        print *, "  expm1(z) = ", z2
        test_passed = .false.
    else
        print *, "PASS: expm1(small z) ≈ z"
    end if
    
    ! Test log1p for small z
    z1 = cmplx(1.0e-8_real64, 1.0e-8_real64, real64)
    z2 = cmplx_log1p(z1)
    if (abs(z2 - z1) > tol * 10.0_real64) then  ! For small z, log1p(z) ≈ z
        print *, "FAIL: log1p(small z) != z"
        test_passed = .false.
    else
        print *, "PASS: log1p(small z) ≈ z"
    end if
    
    ! Test expm1 and log1p are inverses
    z1 = cmplx_make(0.1_real64, 0.1_real64)
    z2 = cmplx_expm1(cmplx_log1p(z1))
    if (abs(z2 - z1) > tol) then
        print *, "FAIL: expm1(log1p(z)) != z"
        test_passed = .false.
    else
        print *, "PASS: expm1(log1p(z)) = z"
    end if
    
    ! Test 4: Cube root
    print *, ""
    print *, "Testing cube root function..."
    
    ! Test cbrt(z)^3 = z
    z1 = cmplx_make(2.0_real64, 3.0_real64)
    z2 = cmplx_cbrt(z1)**3
    if (abs(z2 - z1) > tol * 10.0_real64) then
        print *, "FAIL: cbrt(z)^3 != z"
        test_passed = .false.
    else
        print *, "PASS: cbrt(z)^3 = z"
    end if
    
    ! Test cbrt(8) = 2
    z1 = cmplx(8.0_real64, 0.0_real64, real64)
    z2 = cmplx_cbrt(z1)
    ! Note: There are three cube roots, we check the principal value
    r1 = cmplx_real(z2)
    if (abs(r1 - 2.0_real64) > tol .or. abs(cmplx_imag(z2)) > tol) then
        print *, "FAIL: cbrt(8) != 2"
        print *, "  Result = ", z2
        test_passed = .false.
    else
        print *, "PASS: cbrt(8) = 2"
    end if
    
    ! Test cbrt(-8)
    z1 = cmplx(-8.0_real64, 0.0_real64, real64)
    z2 = cmplx_cbrt(z1)
    z3 = z2**3
    if (abs(z3 - z1) > tol * 10.0_real64) then
        print *, "FAIL: cbrt(-8)^3 != -8"
        test_passed = .false.
    else
        print *, "PASS: cbrt(-8)^3 = -8"
    end if
    
    ! Test 5: Hypot function
    print *, ""
    print *, "Testing hypot function..."
    
    ! Test hypot(3, 4) = 5
    r1 = cmplx_hypot(3.0_real64, 4.0_real64)
    if (abs(r1 - 5.0_real64) > tol) then
        print *, "FAIL: hypot(3, 4) != 5"
        test_passed = .false.
    else
        print *, "PASS: hypot(3, 4) = 5"
    end if
    
    ! Test hypot with large values (overflow protection)
    r1 = cmplx_hypot(1.0e200_real64, 1.0e200_real64)
    r2 = 1.0e200_real64 * sqrt(2.0_real64)
    if (abs(r1 - r2) / r2 > tol) then
        print *, "FAIL: hypot with large values"
        test_passed = .false.
    else
        print *, "PASS: hypot with large values"
    end if
    
    ! Test hypot(0, x) = |x|
    r1 = cmplx_hypot(0.0_real64, -7.5_real64)
    if (abs(r1 - 7.5_real64) > tol) then
        print *, "FAIL: hypot(0, x) != |x|"
        test_passed = .false.
    else
        print *, "PASS: hypot(0, x) = |x|"
    end if
    
    ! Test 6: Complex atan2
    print *, ""
    print *, "Testing complex atan2 function..."
    
    ! Test basic atan2
    z1 = cmplx_E
    z2 = cmplx_I
    z3 = cmplx_atan2(z2, z1)
    ! For real arguments, this should give pi/4
    if (abs(cmplx_real(z3) - pi/4.0_real64) > tol .or. &
        abs(cmplx_imag(z3)) > tol) then
        print *, "FAIL: atan2(i, 1) != pi/4"
        print *, "  Result = ", z3
        test_passed = .false.
    else
        print *, "PASS: atan2(i, 1) ≈ pi/4"
    end if
    
    print *, ""
    if (test_passed) then
        print *, "All advanced transcendental function tests PASSED!"
    else
        print *, "Some advanced transcendental function tests FAILED!"
        stop 1
    end if
    
end program test_kilca_complex_advanced