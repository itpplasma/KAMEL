program test_kilca_complex_validation
    use iso_fortran_env, only: real64, int32
    use kilca_complex_m
    use kilca_constants_m, only: pi
    implicit none
    
    logical :: test_passed
    complex(real64) :: z1, z2, z3
    real(real64) :: tol, rel_tol
    logical :: is_valid
    
    test_passed = .true.
    tol = epsilon(1.0_real64) * 10.0_real64
    rel_tol = 1.0e-12_real64
    
    print *, "Testing complex number validation and testing utilities..."
    print *, ""
    
    ! Test 1: Complex number validation functions
    print *, "Testing complex number validation functions..."
    
    ! Test cmplx_is_finite
    z1 = cmplx(1.0_real64, 2.0_real64, real64)
    is_valid = cmplx_is_finite(z1)
    if (.not. is_valid) then
        print *, "FAIL: cmplx_is_finite failed for finite number"
        test_passed = .false.
    else
        print *, "PASS: cmplx_is_finite correctly identifies finite number"
    end if
    
    ! Test cmplx_is_finite with infinity
    z1 = cmplx(huge(1.0_real64), 1.0_real64, real64)
    is_valid = cmplx_is_finite(z1)
    if (is_valid) then
        print *, "FAIL: cmplx_is_finite should fail for infinite number"
        test_passed = .false.
    else
        print *, "PASS: cmplx_is_finite correctly identifies infinite number"
    end if
    
    ! Test cmplx_is_nan
    z1 = cmplx(1.0_real64, 2.0_real64, real64)
    is_valid = cmplx_is_nan(z1)
    if (is_valid) then
        print *, "FAIL: cmplx_is_nan should be false for normal number"
        test_passed = .false.
    else
        print *, "PASS: cmplx_is_nan correctly identifies normal number"
    end if
    
    ! Test cmplx_is_normal
    z1 = cmplx(1.0_real64, 2.0_real64, real64)
    is_valid = cmplx_is_normal(z1)
    if (.not. is_valid) then
        print *, "FAIL: cmplx_is_normal failed for normal number"
        test_passed = .false.
    else
        print *, "PASS: cmplx_is_normal correctly identifies normal number"
    end if
    
    ! Test cmplx_is_zero
    z1 = cmplx_O
    is_valid = cmplx_is_zero(z1)
    if (.not. is_valid) then
        print *, "FAIL: cmplx_is_zero failed for zero"
        test_passed = .false.
    else
        print *, "PASS: cmplx_is_zero correctly identifies zero"
    end if
    
    z1 = cmplx(1.0_real64, 0.0_real64, real64)
    is_valid = cmplx_is_zero(z1)
    if (is_valid) then
        print *, "FAIL: cmplx_is_zero should be false for non-zero"
        test_passed = .false.
    else
        print *, "PASS: cmplx_is_zero correctly identifies non-zero"
    end if
    
    ! Test 2: Complex number comparison utilities
    print *, ""
    print *, "Testing complex number comparison utilities..."
    
    ! Test cmplx_equals with exact match
    z1 = cmplx(1.0_real64, 2.0_real64, real64)
    z2 = cmplx(1.0_real64, 2.0_real64, real64)
    is_valid = cmplx_equals(z1, z2, tol)
    if (.not. is_valid) then
        print *, "FAIL: cmplx_equals failed for exact match"
        test_passed = .false.
    else
        print *, "PASS: cmplx_equals correctly identifies exact match"
    end if
    
    ! Test cmplx_equals with tolerance
    z1 = cmplx(1.0_real64, 2.0_real64, real64)
    z2 = cmplx(1.0_real64 + tol/2.0_real64, 2.0_real64, real64)
    is_valid = cmplx_equals(z1, z2, tol)
    if (.not. is_valid) then
        print *, "FAIL: cmplx_equals failed within tolerance"
        test_passed = .false.
    else
        print *, "PASS: cmplx_equals correctly works within tolerance"
    end if
    
    ! Test cmplx_equals outside tolerance
    z1 = cmplx(1.0_real64, 2.0_real64, real64)
    z2 = cmplx(1.0_real64 + tol*10.0_real64, 2.0_real64, real64)
    is_valid = cmplx_equals(z1, z2, tol)
    if (is_valid) then
        print *, "FAIL: cmplx_equals should fail outside tolerance"
        test_passed = .false.
    else
        print *, "PASS: cmplx_equals correctly identifies difference outside tolerance"
    end if
    
    ! Test cmplx_approx_zero
    z1 = cmplx(tol/2.0_real64, tol/2.0_real64, real64)
    is_valid = cmplx_approx_zero(z1, tol)
    if (.not. is_valid) then
        print *, "FAIL: cmplx_approx_zero failed for small number"
        test_passed = .false.
    else
        print *, "PASS: cmplx_approx_zero correctly identifies small number"
    end if
    
    ! Test 3: Complex number property testing
    print *, ""
    print *, "Testing complex number property testing..."
    
    ! Test cmplx_is_real
    z1 = cmplx(5.0_real64, 0.0_real64, real64)
    is_valid = cmplx_is_real(z1, tol)
    if (.not. is_valid) then
        print *, "FAIL: cmplx_is_real failed for real number"
        test_passed = .false.
    else
        print *, "PASS: cmplx_is_real correctly identifies real number"
    end if
    
    z1 = cmplx(5.0_real64, 1.0_real64, real64)
    is_valid = cmplx_is_real(z1, tol)
    if (is_valid) then
        print *, "FAIL: cmplx_is_real should fail for complex number"
        test_passed = .false.
    else
        print *, "PASS: cmplx_is_real correctly identifies complex number"
    end if
    
    ! Test cmplx_is_imaginary
    z1 = cmplx(0.0_real64, 5.0_real64, real64)
    is_valid = cmplx_is_imaginary(z1, tol)
    if (.not. is_valid) then
        print *, "FAIL: cmplx_is_imaginary failed for pure imaginary"
        test_passed = .false.
    else
        print *, "PASS: cmplx_is_imaginary correctly identifies pure imaginary"
    end if
    
    ! Test cmplx_is_unit
    z1 = cmplx_E
    is_valid = cmplx_is_unit(z1, tol)
    if (.not. is_valid) then
        print *, "FAIL: cmplx_is_unit failed for unit magnitude"
        test_passed = .false.
    else
        print *, "PASS: cmplx_is_unit correctly identifies unit magnitude"
    end if
    
    z1 = cmplx(1.0_real64/sqrt(2.0_real64), 1.0_real64/sqrt(2.0_real64), real64)
    is_valid = cmplx_is_unit(z1, tol)
    if (.not. is_valid) then
        print *, "FAIL: cmplx_is_unit failed for unit magnitude complex"
        test_passed = .false.
    else
        print *, "PASS: cmplx_is_unit correctly identifies unit magnitude complex"
    end if
    
    ! Test 4: Mathematical relationship validation
    print *, ""
    print *, "Testing mathematical relationship validation..."
    
    ! Test cmplx_verify_conjugate
    z1 = cmplx(3.0_real64, 4.0_real64, real64)
    z2 = cmplx(3.0_real64, -4.0_real64, real64)
    is_valid = cmplx_verify_conjugate(z1, z2, tol)
    if (.not. is_valid) then
        print *, "FAIL: cmplx_verify_conjugate failed for conjugate pair"
        test_passed = .false.
    else
        print *, "PASS: cmplx_verify_conjugate correctly identifies conjugate pair"
    end if
    
    ! Test cmplx_verify_identity
    z1 = cmplx(2.0_real64, 3.0_real64, real64)
    z2 = cmplx_exp(cmplx_log(z1))
    is_valid = cmplx_verify_identity(z1, z2, rel_tol)
    if (.not. is_valid) then
        print *, "FAIL: cmplx_verify_identity failed for exp(log(z)) = z"
        test_passed = .false.
    else
        print *, "PASS: cmplx_verify_identity correctly verifies exp(log(z)) = z"
    end if
    
    ! Test cmplx_verify_euler
    z1 = cmplx(0.0_real64, pi, real64)
    z2 = cmplx_exp(z1) + cmplx_E
    is_valid = cmplx_verify_euler(z2, tol)
    if (.not. is_valid) then
        print *, "FAIL: cmplx_verify_euler failed for Euler's identity"
        test_passed = .false.
    else
        print *, "PASS: cmplx_verify_euler correctly verifies Euler's identity"
    end if
    
    ! Test 5: Error analysis utilities
    print *, ""
    print *, "Testing error analysis utilities..."
    
    ! Test cmplx_relative_error
    z1 = cmplx(100.0_real64, 0.0_real64, real64)
    z2 = cmplx(101.0_real64, 0.0_real64, real64)
    tol = cmplx_relative_error(z1, z2)
    if (abs(tol - 0.01_real64) > 1.0e-10_real64) then
        print *, "FAIL: cmplx_relative_error calculation incorrect"
        test_passed = .false.
    else
        print *, "PASS: cmplx_relative_error correctly calculated"
    end if
    
    ! Test cmplx_absolute_error
    z1 = cmplx(100.0_real64, 0.0_real64, real64)
    z2 = cmplx(101.0_real64, 0.0_real64, real64)
    tol = cmplx_absolute_error(z1, z2)
    if (abs(tol - 1.0_real64) > 1.0e-10_real64) then
        print *, "FAIL: cmplx_absolute_error calculation incorrect"
        test_passed = .false.
    else
        print *, "PASS: cmplx_absolute_error correctly calculated"
    end if
    
    print *, ""
    if (test_passed) then
        print *, "All complex number validation and testing utility tests PASSED!"
    else
        print *, "Some complex number validation and testing utility tests FAILED!"
        stop 1
    end if
    
end program test_kilca_complex_validation