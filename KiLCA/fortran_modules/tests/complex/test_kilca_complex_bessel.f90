program test_kilca_complex_bessel
    use iso_fortran_env, only: real64, int32
    use kilca_complex_m
    use kilca_constants_m, only: pi
    implicit none
    
    logical :: test_passed
    complex(real64) :: z1, z2, z3, expected
    real(real64) :: tol, x
    integer :: nu, n
    
    test_passed = .true.
    tol = epsilon(1.0_real64) * 10000.0_real64  ! More relaxed tolerance for special functions
    
    print *, "Testing complex Bessel function interfaces..."
    print *, ""
    
    ! Test 1: Bessel J functions
    print *, "Testing Bessel J functions..."
    
    ! Test J_0(0) = 1
    z1 = cmplx_O
    z2 = cmplx_besselj(0, z1)
    expected = cmplx_E
    if (abs(z2 - expected) > tol) then
        print *, "FAIL: J_0(0) != 1"
        print *, "  Got: ", z2
        test_passed = .false.
    else
        print *, "PASS: J_0(0) = 1"
    end if
    
    ! Test J_1(0) = 0
    z2 = cmplx_besselj(1, z1)
    if (abs(z2) > tol) then
        print *, "FAIL: J_1(0) != 0"
        test_passed = .false.
    else
        print *, "PASS: J_1(0) = 0"
    end if
    
    ! Test J_0(x) for real x (should match known values)
    x = 1.0_real64
    z1 = cmplx(x, 0.0_real64, real64)
    z2 = cmplx_besselj(0, z1)
    ! J_0(1) ≈ 0.765197686557966551
    expected = cmplx(0.765197686557966551_real64, 0.0_real64, real64)
    if (abs(z2 - expected) > tol) then
        print *, "FAIL: J_0(1) doesn't match expected value"
        print *, "  Expected: ", expected
        print *, "  Got: ", z2
        test_passed = .false.
    else
        print *, "PASS: J_0(1) matches expected value"
    end if
    
    ! Test recurrence relation: J_{n-1}(z) + J_{n+1}(z) = (2n/z) * J_n(z)
    z1 = cmplx(2.0_real64, 1.0_real64, real64)
    nu = 1
    z2 = cmplx_besselj(nu-1, z1) + cmplx_besselj(nu+1, z1)
    z3 = (2.0_real64 * nu / z1) * cmplx_besselj(nu, z1)
    if (abs(z2 - z3) / abs(z3) > tol * 10.0_real64) then
        print *, "FAIL: Bessel J recurrence relation"
        print *, "  LHS: ", z2
        print *, "  RHS: ", z3
        test_passed = .false.
    else
        print *, "PASS: Bessel J recurrence relation"
    end if
    
    ! Test 2: Modified Bessel I functions
    print *, ""
    print *, "Testing modified Bessel I functions..."
    
    ! Test I_0(0) = 1
    z1 = cmplx_O
    z2 = cmplx_besseli(0, z1)
    expected = cmplx_E
    if (abs(z2 - expected) > tol) then
        print *, "FAIL: I_0(0) != 1"
        test_passed = .false.
    else
        print *, "PASS: I_0(0) = 1"
    end if
    
    ! Test I_1(0) = 0
    z2 = cmplx_besseli(1, z1)
    if (abs(z2) > tol) then
        print *, "FAIL: I_1(0) != 0"
        test_passed = .false.
    else
        print *, "PASS: I_1(0) = 0"
    end if
    
    ! Test relation: I_n(z) = (-i)^n J_n(iz)
    z1 = cmplx(1.5_real64, 0.5_real64, real64)
    nu = 1
    z2 = cmplx_besseli(nu, z1)
    z3 = (-cmplx_I)**nu * cmplx_besselj(nu, cmplx_I * z1)
    if (abs(z2 - z3) / abs(z2) > tol * 100.0_real64) then
        print *, "FAIL: I_n(z) = (-i)^n J_n(iz) relation"
        print *, "  I_n(z): ", z2
        print *, "  (-i)^n J_n(iz): ", z3
        test_passed = .false.
    else
        print *, "PASS: I_n(z) = (-i)^n J_n(iz) relation"
    end if
    
    ! Test 3: Bessel Y functions
    print *, ""
    print *, "Testing Bessel Y functions..."
    
    ! Test Y_0(1) (known value)
    x = 1.0_real64
    z1 = cmplx(x, 0.0_real64, real64)
    z2 = cmplx_bessely(0, z1)
    ! Y_0(1) ≈ 0.08825696421567696
    expected = cmplx(0.08825696421567696_real64, 0.0_real64, real64)
    if (abs(z2 - expected) > tol * 10.0_real64) then
        print *, "FAIL: Y_0(1) doesn't match expected value"
        print *, "  Expected: ", expected
        print *, "  Got: ", z2
        test_passed = .false.
    else
        print *, "PASS: Y_0(1) matches expected value"
    end if
    
    ! Test 4: Modified Bessel K functions
    print *, ""
    print *, "Testing modified Bessel K functions..."
    
    ! Test K_0(1) (known value)
    x = 1.0_real64
    z1 = cmplx(x, 0.0_real64, real64)
    z2 = cmplx_besselk(0, z1)
    ! K_0(1) ≈ 0.421024438240708333
    expected = cmplx(0.421024438240708333_real64, 0.0_real64, real64)
    if (abs(z2 - expected) > tol * 10.0_real64) then
        print *, "FAIL: K_0(1) doesn't match expected value"
        print *, "  Expected: ", expected
        print *, "  Got: ", z2
        test_passed = .false.
    else
        print *, "PASS: K_0(1) matches expected value"
    end if
    
    ! Test 5: Derivatives
    print *, ""
    print *, "Testing Bessel function derivatives..."
    
    ! Test J_0'(z) = -J_1(z)
    z1 = cmplx(1.5_real64, 0.8_real64, real64)
    z2 = cmplx_besselj_derivative(0, z1, 1)
    z3 = -cmplx_besselj(1, z1)
    if (abs(z2 - z3) / abs(z3) > tol * 10.0_real64) then
        print *, "FAIL: J_0'(z) != -J_1(z)"
        test_passed = .false.
    else
        print *, "PASS: J_0'(z) = -J_1(z)"
    end if
    
    print *, ""
    if (test_passed) then
        print *, "All Bessel function interface tests PASSED!"
    else
        print *, "Some Bessel function interface tests FAILED!"
        stop 1
    end if
    
end program test_kilca_complex_bessel