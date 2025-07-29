program test_kilca_complex_safe
    use iso_fortran_env, only: real64, int32
    use kilca_complex_m
    use kilca_constants_m, only: pi
    implicit none
    
    logical :: test_passed
    complex(real64) :: z1, z2, z3, z4
    real(real64) :: r1, r2, tol
    
    test_passed = .true.
    tol = epsilon(1.0_real64) * 1000.0_real64
    
    print *, "Testing safe complex logarithmic and power functions..."
    print *, ""
    
    ! Test 1: cmplx_logabs function
    print *, "Testing cmplx_logabs function..."
    
    ! Test normal case
    z1 = cmplx_make(3.0_real64, 4.0_real64)
    r1 = cmplx_logabs(z1)
    r2 = log(abs(z1))
    if (abs(r1 - r2) > tol) then
        print *, "FAIL: cmplx_logabs(z) != log(abs(z))"
        test_passed = .false.
    else
        print *, "PASS: cmplx_logabs gives correct result"
    end if
    
    ! Test with very large numbers
    z1 = cmplx(1.0e200_real64, 1.0e200_real64, real64)
    r1 = cmplx_logabs(z1)
    r2 = log(sqrt(2.0_real64)) + log(1.0e200_real64)
    if (abs(r1 - r2) / abs(r2) > tol) then
        print *, "FAIL: cmplx_logabs fails for large numbers"
        test_passed = .false.
    else
        print *, "PASS: cmplx_logabs handles large numbers"
    end if
    
    ! Test with very small numbers
    z1 = cmplx(1.0e-200_real64, 1.0e-200_real64, real64)
    r1 = cmplx_logabs(z1)
    r2 = log(sqrt(2.0_real64)) + log(1.0e-200_real64)
    if (abs(r1 - r2) / abs(r2) > tol) then
        print *, "FAIL: cmplx_logabs fails for small numbers"
        test_passed = .false.
    else
        print *, "PASS: cmplx_logabs handles small numbers"
    end if
    
    ! Test with zero
    z1 = cmplx_O
    r1 = cmplx_logabs(z1)
    if (r1 /= -huge(1.0_real64)) then
        print *, "FAIL: cmplx_logabs(0) != -infinity"
        test_passed = .false.
    else
        print *, "PASS: cmplx_logabs(0) = -infinity"
    end if
    
    ! Test 2: cmplx_pow_safe function
    print *, ""
    print *, "Testing cmplx_pow_safe function..."
    
    ! Test normal case
    z1 = cmplx_make(2.0_real64, 1.0_real64)
    z2 = cmplx_make(3.0_real64, 0.5_real64)
    z3 = cmplx_pow_safe(z1, z2)
    z4 = cmplx_pow_complex(z1, z2)
    if (abs(z3 - z4) / abs(z4) > tol) then
        print *, "FAIL: cmplx_pow_safe != cmplx_pow_complex for normal values"
        test_passed = .false.
    else
        print *, "PASS: cmplx_pow_safe agrees with standard power"
    end if
    
    ! Test large base
    z1 = cmplx(1.0e100_real64, 0.0_real64, real64)
    z2 = cmplx(2.0_real64, 0.0_real64, real64)
    z3 = cmplx_pow_safe(z1, z2)
    ! Should not overflow since 1e100^2 = 1e200 < huge
    if (abs(z3) < 1.0e199_real64 .or. abs(z3) > 1.0e201_real64) then
        print *, "FAIL: cmplx_pow_safe fails for large base"
        test_passed = .false.
    else
        print *, "PASS: cmplx_pow_safe handles large base"
    end if
    
    ! Test potential overflow
    z1 = cmplx(1.0e200_real64, 0.0_real64, real64)
    z2 = cmplx(10.0_real64, 0.0_real64, real64)
    z3 = cmplx_pow_safe(z1, z2)
    ! Should saturate to huge
    if (real(z3, real64) /= huge(1.0_real64)) then
        print *, "FAIL: cmplx_pow_safe doesn't handle overflow correctly"
        test_passed = .false.
    else
        print *, "PASS: cmplx_pow_safe handles overflow"
    end if
    
    ! Test potential underflow
    z1 = cmplx(1.0e-200_real64, 0.0_real64, real64)
    z2 = cmplx(10.0_real64, 0.0_real64, real64)
    z3 = cmplx_pow_safe(z1, z2)
    ! Should underflow to zero
    if (abs(z3) /= 0.0_real64) then
        print *, "FAIL: cmplx_pow_safe doesn't handle underflow correctly"
        test_passed = .false.
    else
        print *, "PASS: cmplx_pow_safe handles underflow"
    end if
    
    ! Test 0^0
    z1 = cmplx_O
    z2 = cmplx_O
    z3 = cmplx_pow_safe(z1, z2)
    if (z3 /= cmplx_E) then
        print *, "FAIL: cmplx_pow_safe: 0^0 != 1"
        test_passed = .false.
    else
        print *, "PASS: cmplx_pow_safe: 0^0 = 1"
    end if
    
    ! Test 0^positive
    z1 = cmplx_O
    z2 = cmplx(2.5_real64, 1.0_real64, real64)
    z3 = cmplx_pow_safe(z1, z2)
    if (z3 /= cmplx_O) then
        print *, "FAIL: cmplx_pow_safe: 0^positive != 0"
        test_passed = .false.
    else
        print *, "PASS: cmplx_pow_safe: 0^positive = 0"
    end if
    
    ! Test 0^negative
    z1 = cmplx_O
    z2 = cmplx(-1.0_real64, 0.0_real64, real64)
    z3 = cmplx_pow_safe(z1, z2)
    if (real(z3, real64) /= huge(1.0_real64)) then
        print *, "FAIL: cmplx_pow_safe: 0^negative != infinity"
        test_passed = .false.
    else
        print *, "PASS: cmplx_pow_safe: 0^negative = infinity"
    end if
    
    print *, ""
    if (test_passed) then
        print *, "All safe function tests PASSED!"
    else
        print *, "Some safe function tests FAILED!"
        stop 1
    end if
    
end program test_kilca_complex_safe