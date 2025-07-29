program test_kilca_complex_hypergeometric
    use iso_fortran_env, only: real64, int32
    use kilca_complex_m
    use kilca_constants_m, only: pi
    implicit none
    
    logical :: test_passed
    complex(real64) :: z1, z2, z3, expected
    real(real64) :: tol, a_val, b_val
    
    test_passed = .true.
    tol = epsilon(1.0_real64) * 1000.0_real64  ! Relaxed tolerance for special functions
    
    print *, "Testing complex hypergeometric function interfaces..."
    print *, ""
    
    ! Test 1: Confluent hypergeometric function 1F1(a,b,z)
    print *, "Testing confluent hypergeometric function 1F1(a,b,z)..."
    
    ! Test 1F1(1,1,z) = exp(z)
    z1 = cmplx(0.5_real64, 0.0_real64, real64)
    z2 = cmplx_hyp1f1(1.0_real64, 1.0_real64, z1)
    expected = cmplx_exp(z1)
    if (abs(z2 - expected) > tol) then
        print *, "FAIL: 1F1(1,1,z) != exp(z)"
        print *, "  Expected: ", expected
        print *, "  Got: ", z2
        test_passed = .false.
    else
        print *, "PASS: 1F1(1,1,z) = exp(z)"
    end if
    
    ! Test 1F1(1,2,z) = (exp(z) - 1) / z
    z1 = cmplx(0.5_real64, 0.0_real64, real64)
    z2 = cmplx_hyp1f1(1.0_real64, 2.0_real64, z1)
    expected = (cmplx_exp(z1) - cmplx_E) / z1
    if (abs(z2 - expected) > tol) then
        print *, "FAIL: 1F1(1,2,z) != (exp(z)-1)/z"
        print *, "  Expected: ", expected
        print *, "  Got: ", z2
        test_passed = .false.
    else
        print *, "PASS: 1F1(1,2,z) = (exp(z)-1)/z"
    end if
    
    ! Test 1F1(0,b,z) = 1 for any b != 0
    z1 = cmplx(2.0_real64, 1.0_real64, real64)
    z2 = cmplx_hyp1f1(0.0_real64, 3.0_real64, z1)
    expected = cmplx_E
    if (abs(z2 - expected) > tol) then
        print *, "FAIL: 1F1(0,b,z) != 1"
        print *, "  Expected: ", expected
        print *, "  Got: ", z2
        test_passed = .false.
    else
        print *, "PASS: 1F1(0,b,z) = 1"
    end if
    
    ! Test 1F1(a,a,z) = exp(z) for any a != 0
    a_val = 2.5_real64
    z1 = cmplx(0.3_real64, 0.2_real64, real64)
    z2 = cmplx_hyp1f1(a_val, a_val, z1)
    expected = cmplx_exp(z1)
    if (abs(z2 - expected) > tol * 10.0_real64) then
        print *, "FAIL: 1F1(a,a,z) != exp(z)"
        print *, "  Expected: ", expected
        print *, "  Got: ", z2
        test_passed = .false.
    else
        print *, "PASS: 1F1(a,a,z) = exp(z)"
    end if
    
    ! Test relation: d/dz 1F1(a,b,z) = (a/b) * 1F1(a+1,b+1,z)
    a_val = 1.5_real64
    b_val = 2.5_real64
    z1 = cmplx(0.5_real64, 0.3_real64, real64)
    z2 = cmplx_hyp1f1_derivative(a_val, b_val, z1, 1)
    z3 = (a_val / b_val) * cmplx_hyp1f1(a_val + 1.0_real64, b_val + 1.0_real64, z1)
    if (abs(z2 - z3) / abs(z3) > tol * 100.0_real64) then
        print *, "FAIL: d/dz 1F1(a,b,z) != (a/b)*1F1(a+1,b+1,z)"
        print *, "  Derivative: ", z2
        print *, "  Relation: ", z3
        test_passed = .false.
    else
        print *, "PASS: d/dz 1F1(a,b,z) = (a/b)*1F1(a+1,b+1,z)"
    end if
    
    ! Test 2: Gauss hypergeometric function 2F1(a,b,c,z)
    print *, ""
    print *, "Testing Gauss hypergeometric function 2F1(a,b,c,z)..."
    
    ! Test 2F1(a,b,b,z) = (1-z)^(-a)
    a_val = 1.5_real64
    b_val = 2.0_real64
    z1 = cmplx(0.3_real64, 0.0_real64, real64)
    z2 = cmplx_hyp2f1(a_val, b_val, b_val, z1)
    expected = cmplx_pow_real(cmplx_E - z1, -a_val)
    if (abs(z2 - expected) > tol * 10.0_real64) then
        print *, "FAIL: 2F1(a,b,b,z) != (1-z)^(-a)"
        print *, "  Expected: ", expected
        print *, "  Got: ", z2
        test_passed = .false.
    else
        print *, "PASS: 2F1(a,b,b,z) = (1-z)^(-a)"
    end if
    
    ! Test 2F1(1,1,2,z) = -log(1-z)/z
    z1 = cmplx(0.5_real64, 0.0_real64, real64)
    z2 = cmplx_hyp2f1(1.0_real64, 1.0_real64, 2.0_real64, z1)
    expected = -cmplx_log(cmplx_E - z1) / z1
    if (abs(z2 - expected) > tol * 10.0_real64) then
        print *, "FAIL: 2F1(1,1,2,z) != -log(1-z)/z"
        print *, "  Expected: ", expected
        print *, "  Got: ", z2
        test_passed = .false.
    else
        print *, "PASS: 2F1(1,1,2,z) = -log(1-z)/z"
    end if
    
    ! Test 2F1(0,b,c,z) = 1 for any b,c
    z1 = cmplx(0.7_real64, 0.2_real64, real64)
    z2 = cmplx_hyp2f1(0.0_real64, 1.5_real64, 2.5_real64, z1)
    expected = cmplx_E
    if (abs(z2 - expected) > tol) then
        print *, "FAIL: 2F1(0,b,c,z) != 1"
        print *, "  Expected: ", expected
        print *, "  Got: ", z2
        test_passed = .false.
    else
        print *, "PASS: 2F1(0,b,c,z) = 1"
    end if
    
    print *, ""
    if (test_passed) then
        print *, "All hypergeometric function interface tests PASSED!"
    else
        print *, "Some hypergeometric function interface tests FAILED!"
        stop 1
    end if
    
end program test_kilca_complex_hypergeometric