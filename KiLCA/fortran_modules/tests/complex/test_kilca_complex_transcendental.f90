program test_kilca_complex_transcendental
    use iso_fortran_env, only: real64, int32
    use kilca_complex_m
    use kilca_constants_m, only: pi
    implicit none
    
    logical :: test_passed
    complex(real64) :: z1, z2, z3, expected
    real(real64) :: tol, re, im
    integer :: i
    
    test_passed = .true.
    tol = epsilon(1.0_real64) * 100.0_real64
    
    print *, "Testing kilca_complex_m transcendental functions..."
    print *, ""
    
    ! Test 1: Exponential and logarithm
    print *, "Testing exponential and logarithm functions..."
    
    ! Test exp(ln(z)) = z
    z1 = cmplx_make(2.0_real64, 3.0_real64)
    z2 = cmplx_exp(cmplx_log(z1))
    if (abs(z2 - z1) > tol) then
        print *, "FAIL: exp(log(z)) != z"
        test_passed = .false.
    else
        print *, "PASS: exp(log(z)) = z"
    end if
    
    ! Test Euler's identity: exp(i*pi) = -1
    z1 = cmplx_I * pi
    z2 = cmplx_exp(z1)
    expected = cmplx(-1.0_real64, 0.0_real64, real64)
    if (abs(z2 - expected) > tol) then
        print *, "FAIL: exp(i*pi) != -1"
        test_passed = .false.
    else
        print *, "PASS: exp(i*pi) = -1"
    end if
    
    ! Test log10
    z1 = cmplx(100.0_real64, 0.0_real64, real64)
    z2 = cmplx_log10(z1)
    expected = cmplx(2.0_real64, 0.0_real64, real64)
    if (abs(z2 - expected) > tol) then
        print *, "FAIL: log10(100) != 2"
        test_passed = .false.
    else
        print *, "PASS: log10(100) = 2"
    end if
    
    ! Test 2: Trigonometric functions
    print *, ""
    print *, "Testing trigonometric functions..."
    
    ! Test sin^2 + cos^2 = 1
    z1 = cmplx_make(1.5_real64, 2.0_real64)
    z2 = cmplx_sin(z1)**2 + cmplx_cos(z1)**2
    expected = cmplx_E
    if (abs(z2 - expected) > tol * 10.0_real64) then
        print *, "FAIL: sin^2 + cos^2 != 1"
        print *, "  Result: ", z2
        test_passed = .false.
    else
        print *, "PASS: sin^2 + cos^2 = 1"
    end if
    
    ! Test tan = sin/cos
    z2 = cmplx_tan(z1)
    z3 = cmplx_sin(z1) / cmplx_cos(z1)
    if (abs(z2 - z3) > tol * 10.0_real64) then
        print *, "FAIL: tan != sin/cos"
        test_passed = .false.
    else
        print *, "PASS: tan = sin/cos"
    end if
    
    ! Test periodicity: sin(z + 2*pi) = sin(z)
    z2 = cmplx_sin(z1)
    z3 = cmplx_sin(z1 + 2.0_real64 * pi)
    if (abs(z2 - z3) > tol * 100.0_real64) then
        print *, "FAIL: sin(z + 2*pi) != sin(z)"
        test_passed = .false.
    else
        print *, "PASS: sin(z + 2*pi) = sin(z)"
    end if
    
    ! Test 3: Hyperbolic functions
    print *, ""
    print *, "Testing hyperbolic functions..."
    
    ! Test cosh^2 - sinh^2 = 1
    z1 = cmplx_make(0.5_real64, 1.0_real64)
    z2 = cmplx_cosh(z1)**2 - cmplx_sinh(z1)**2
    expected = cmplx_E
    if (abs(z2 - expected) > tol * 10.0_real64) then
        print *, "FAIL: cosh^2 - sinh^2 != 1"
        print *, "  Result: ", z2
        test_passed = .false.
    else
        print *, "PASS: cosh^2 - sinh^2 = 1"
    end if
    
    ! Test tanh = sinh/cosh
    z2 = cmplx_tanh(z1)
    z3 = cmplx_sinh(z1) / cmplx_cosh(z1)
    if (abs(z2 - z3) > tol * 10.0_real64) then
        print *, "FAIL: tanh != sinh/cosh"
        test_passed = .false.
    else
        print *, "PASS: tanh = sinh/cosh"
    end if
    
    ! Test relation: sinh(z) = -i*sin(i*z)
    z2 = cmplx_sinh(z1)
    z3 = -cmplx_I * cmplx_sin(cmplx_I * z1)
    if (abs(z2 - z3) > tol * 10.0_real64) then
        print *, "FAIL: sinh(z) != -i*sin(i*z)"
        test_passed = .false.
    else
        print *, "PASS: sinh(z) = -i*sin(i*z)"
    end if
    
    ! Test 4: Inverse trigonometric functions
    print *, ""
    print *, "Testing inverse trigonometric functions..."
    
    ! Test sin(asin(z)) = z for |z| < 1
    z1 = cmplx_make(0.5_real64, 0.3_real64)
    z2 = cmplx_sin(cmplx_asin(z1))
    if (abs(z2 - z1) > tol * 100.0_real64) then
        print *, "FAIL: sin(asin(z)) != z"
        test_passed = .false.
    else
        print *, "PASS: sin(asin(z)) = z"
    end if
    
    ! Test cos(acos(z)) = z
    z2 = cmplx_cos(cmplx_acos(z1))
    if (abs(z2 - z1) > tol * 100.0_real64) then
        print *, "FAIL: cos(acos(z)) != z"
        test_passed = .false.
    else
        print *, "PASS: cos(acos(z)) = z"
    end if
    
    ! Test tan(atan(z)) = z
    z2 = cmplx_tan(cmplx_atan(z1))
    if (abs(z2 - z1) > tol * 100.0_real64) then
        print *, "FAIL: tan(atan(z)) != z"
        test_passed = .false.
    else
        print *, "PASS: tan(atan(z)) = z"
    end if
    
    ! Test 5: Square root
    print *, ""
    print *, "Testing square root function..."
    
    ! Test sqrt(z)^2 = z
    z1 = cmplx_make(3.0_real64, 4.0_real64)
    z2 = cmplx_sqrt(z1)**2
    if (abs(z2 - z1) > tol * 10.0_real64) then
        print *, "FAIL: sqrt(z)^2 != z"
        test_passed = .false.
    else
        print *, "PASS: sqrt(z)^2 = z"
    end if
    
    ! Test sqrt(-1) = i
    z1 = cmplx(-1.0_real64, 0.0_real64, real64)
    z2 = cmplx_sqrt(z1)
    ! Note: There are two square roots, we check if it's ±i
    if (abs(abs(z2) - 1.0_real64) > tol .or. &
        abs(abs(cmplx_imag(z2)) - 1.0_real64) > tol) then
        print *, "FAIL: sqrt(-1) != ±i"
        test_passed = .false.
    else
        print *, "PASS: sqrt(-1) = ±i"
    end if
    
    ! Test 6: Edge cases
    print *, ""
    print *, "Testing edge cases..."
    
    ! Test functions at zero
    z1 = cmplx_O
    
    ! sin(0) = 0
    z2 = cmplx_sin(z1)
    if (abs(z2) > tol) then
        print *, "FAIL: sin(0) != 0"
        test_passed = .false.
    else
        print *, "PASS: sin(0) = 0"
    end if
    
    ! cos(0) = 1
    z2 = cmplx_cos(z1)
    if (abs(z2 - cmplx_E) > tol) then
        print *, "FAIL: cos(0) != 1"
        test_passed = .false.
    else
        print *, "PASS: cos(0) = 1"
    end if
    
    ! exp(0) = 1
    z2 = cmplx_exp(z1)
    if (abs(z2 - cmplx_E) > tol) then
        print *, "FAIL: exp(0) != 1"
        test_passed = .false.
    else
        print *, "PASS: exp(0) = 1"
    end if
    
    ! Test 7: Large argument behavior
    print *, ""
    print *, "Testing large argument behavior..."
    
    ! Test exp with large imaginary part (should oscillate)
    z1 = cmplx(0.0_real64, 100.0_real64, real64)
    z2 = cmplx_exp(z1)
    re = cmplx_real(z2)
    im = cmplx_imag(z2)
    if (abs(abs(z2) - 1.0_real64) > tol * 10.0_real64) then
        print *, "FAIL: |exp(100i)| != 1"
        test_passed = .false.
    else
        print *, "PASS: |exp(100i)| = 1"
    end if
    
    print *, ""
    if (test_passed) then
        print *, "All transcendental function tests PASSED!"
    else
        print *, "Some transcendental function tests FAILED!"
        stop 1
    end if
    
end program test_kilca_complex_transcendental