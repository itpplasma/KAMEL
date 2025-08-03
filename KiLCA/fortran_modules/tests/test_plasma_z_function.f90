program test_plasma_z_function
    use iso_fortran_env, only: real64
    use kilca_plasma_physics_m
    implicit none
    
    integer :: test_status = 0
    real(real64), parameter :: tolerance = 1.0e-6_real64
    
    print *, "==========================================="
    print *, "Testing Plasma Dispersion Z-Function [RED PHASE - SHOULD FAIL]"
    print *, "==========================================="
    
    call test_z_function_known_values()
    call test_z_function_properties()
    call test_z_function_asymptotic_limits()
    call test_z_function_derivatives()
    call test_z_function_series_expansion()
    
    if (test_status == 0) then
        print *, ""
        print *, "Plasma Z-function tests PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "Plasma Z-function tests FAILED:", test_status, "test(s) failed (expected in RED phase)"
        stop 1
    end if

contains

    !> Test Z-function against known tabulated values
    subroutine test_z_function_known_values()
        complex(real64) :: zeta, z_val, z_expected
        real(real64) :: error
        
        print *, ""
        print *, "Testing Z-function known values..."
        print *, ""
        
        ! Test 1: Z(0) = i*sqrt(pi)
        zeta = cmplx(0.0_real64, 0.0_real64, real64)
        z_val = plasma_z_function(zeta)
        z_expected = cmplx(0.0_real64, sqrt(3.14159265358979323846_real64), real64)
        error = abs(z_val - z_expected)
        
        if (error > tolerance) then
            print *, "FAIL: Z(0) incorrect. Got:", z_val, "Expected:", z_expected
            test_status = test_status + 1
        else
            print *, "PASS: Z(0) = i*sqrt(pi)"
        end if
        
        ! Test 2: Z(i) ≈ 0.757872 + 0.394936i (tabulated value)
        zeta = cmplx(0.0_real64, 1.0_real64, real64)
        z_val = plasma_z_function(zeta)
        z_expected = cmplx(0.757872_real64, 0.394936_real64, real64)
        error = abs(z_val - z_expected)
        
        if (error > 0.001_real64) then  ! Relaxed tolerance for complex value
            print *, "FAIL: Z(i) incorrect. Got:", z_val, "Expected:", z_expected
            test_status = test_status + 1
        else
            print *, "PASS: Z(i) matches tabulated value"
        end if
        
        ! Test 3: Z(1) ≈ -1.07689 + 0.628573i (tabulated value)
        zeta = cmplx(1.0_real64, 0.0_real64, real64)
        z_val = plasma_z_function(zeta)
        z_expected = cmplx(-1.07689_real64, 0.628573_real64, real64)
        error = abs(z_val - z_expected)
        
        if (error > 0.001_real64) then
            print *, "FAIL: Z(1) incorrect. Got:", z_val, "Expected:", z_expected
            test_status = test_status + 1
        else
            print *, "PASS: Z(1) matches tabulated value"
        end if
        
        ! Test 4: Z(2+2i) ≈ -0.0161043 - 0.116506i (tabulated value)
        zeta = cmplx(2.0_real64, 2.0_real64, real64)
        z_val = plasma_z_function(zeta)
        z_expected = cmplx(-0.0161043_real64, -0.116506_real64, real64)
        error = abs(z_val - z_expected)
        
        if (error > 0.001_real64) then
            print *, "FAIL: Z(2+2i) incorrect. Got:", z_val, "Expected:", z_expected
            test_status = test_status + 1
        else
            print *, "PASS: Z(2+2i) matches tabulated value"
        end if
        
    end subroutine test_z_function_known_values
    
    !> Test Z-function mathematical properties
    subroutine test_z_function_properties()
        complex(real64) :: zeta, z_val, z_conj, z_val_conj
        real(real64) :: error
        
        print *, ""
        print *, "Testing Z-function properties..."
        print *, ""
        
        ! Test 1: Symmetry property: Z(-ζ*) = -Z(ζ)*
        zeta = cmplx(1.5_real64, 0.5_real64, real64)
        z_val = plasma_z_function(zeta)
        z_conj = plasma_z_function(-conjg(zeta))
        z_val_conj = -conjg(z_val)
        error = abs(z_conj - z_val_conj)
        
        if (error > tolerance) then
            print *, "FAIL: Symmetry property Z(-ζ*) = -Z(ζ)* not satisfied"
            test_status = test_status + 1
        else
            print *, "PASS: Symmetry property satisfied"
        end if
        
        ! Test 2: Real axis property: Im[Z(x)] > 0 for real x
        zeta = cmplx(1.0_real64, 0.0_real64, real64)
        z_val = plasma_z_function(zeta)
        
        if (aimag(z_val) <= 0.0_real64) then
            print *, "FAIL: Im[Z(x)] should be positive for real x. Got:", aimag(z_val)
            test_status = test_status + 1
        else
            print *, "PASS: Im[Z(x)] > 0 for real x"
        end if
        
        ! Test 3: Recurrence relation: Z'(ζ) = -2(1 + ζZ(ζ))
        zeta = cmplx(0.5_real64, 0.5_real64, real64)
        z_val = plasma_z_function(zeta)
        z_val_conj = plasma_z_function_derivative(zeta)
        z_val_conj = -2.0_real64 * (1.0_real64 + zeta * z_val)  ! Expected value
        error = abs(plasma_z_function_derivative(zeta) - z_val_conj)
        
        if (error > tolerance) then
            print *, "FAIL: Recurrence relation Z'(ζ) = -2(1 + ζZ(ζ)) not satisfied"
            test_status = test_status + 1
        else
            print *, "PASS: Recurrence relation satisfied"
        end if
        
    end subroutine test_z_function_properties
    
    !> Test Z-function asymptotic limits
    subroutine test_z_function_asymptotic_limits()
        complex(real64) :: zeta, z_val, z_asymptotic
        real(real64) :: error
        
        print *, ""
        print *, "Testing Z-function asymptotic limits..."
        print *, ""
        
        ! Test 1: Large argument limit |ζ| >> 1: Z(ζ) ≈ -1/ζ - 1/(2ζ³)
        zeta = cmplx(10.0_real64, 0.0_real64, real64)
        z_val = plasma_z_function(zeta)
        z_asymptotic = -1.0_real64/zeta - 1.0_real64/(2.0_real64*zeta**3)
        error = abs(z_val - z_asymptotic) / abs(z_asymptotic)
        
        if (error > 0.01_real64) then  ! 1% relative error for asymptotic
            print *, "FAIL: Large argument asymptotic incorrect. Error:", error
            test_status = test_status + 1
        else
            print *, "PASS: Large argument asymptotic correct"
        end if
        
        ! Test 2: Small imaginary part limit
        zeta = cmplx(0.0_real64, 0.001_real64, real64)
        z_val = plasma_z_function(zeta)
        
        ! For small Im(ζ), should be close to i*sqrt(pi)
        z_asymptotic = cmplx(0.0_real64, sqrt(3.14159265358979323846_real64), real64)
        error = abs(aimag(z_val) - aimag(z_asymptotic)) / abs(aimag(z_asymptotic))
        
        if (error > 0.01_real64) then
            print *, "FAIL: Small imaginary part limit incorrect"
            test_status = test_status + 1
        else
            print *, "PASS: Small imaginary part limit correct"
        end if
        
        ! Test 3: Large negative real part
        zeta = cmplx(-10.0_real64, 0.01_real64, real64)
        z_val = plasma_z_function(zeta)
        
        ! Should approach 2*exp(-ζ²) for large negative Re(ζ)
        z_asymptotic = 2.0_real64 * exp(-zeta**2)
        error = abs(z_val - z_asymptotic)
        
        ! This is a weak test due to exponential decay
        if (abs(real(z_val)) > 1.0_real64) then
            print *, "FAIL: Large negative real part limit incorrect"
            test_status = test_status + 1
        else
            print *, "PASS: Large negative real part limit reasonable"
        end if
        
    end subroutine test_z_function_asymptotic_limits
    
    !> Test Z-function derivatives
    subroutine test_z_function_derivatives()
        complex(real64) :: zeta, z_val, z_prime, z_prime_expected
        real(real64) :: error
        
        print *, ""
        print *, "Testing Z-function derivatives..."
        print *, ""
        
        ! Test 1: First derivative at ζ = 0.5
        zeta = cmplx(0.5_real64, 0.0_real64, real64)
        z_val = plasma_z_function(zeta)
        z_prime = plasma_z_function_derivative(zeta)
        z_prime_expected = -2.0_real64 * (1.0_real64 + zeta * z_val)
        error = abs(z_prime - z_prime_expected)
        
        if (error > tolerance) then
            print *, "FAIL: First derivative incorrect. Got:", z_prime, "Expected:", z_prime_expected
            test_status = test_status + 1
        else
            print *, "PASS: First derivative Z'(ζ) correct"
        end if
        
        ! Test 2: Second derivative relation: Z''(ζ) = -2(Z(ζ) + ζZ'(ζ))
        z_val = plasma_z_function(zeta)
        z_prime = plasma_z_function_derivative(zeta)
        z_prime_expected = plasma_z_function_second_derivative(zeta)
        z_prime_expected = -2.0_real64 * (z_val + zeta * z_prime)  ! Expected value
        error = abs(plasma_z_function_second_derivative(zeta) - z_prime_expected)
        
        if (error > tolerance) then
            print *, "FAIL: Second derivative relation incorrect"
            test_status = test_status + 1
        else
            print *, "PASS: Second derivative Z''(ζ) correct"
        end if
        
    end subroutine test_z_function_derivatives
    
    !> Test Z-function series expansion
    subroutine test_z_function_series_expansion()
        complex(real64) :: zeta, z_val, z_series
        real(real64) :: error
        integer :: n
        
        print *, ""
        print *, "Testing Z-function series expansion..."
        print *, ""
        
        ! Test 1: Taylor series for small |ζ|
        zeta = cmplx(0.1_real64, 0.1_real64, real64)
        z_val = plasma_z_function(zeta)
        
        ! Taylor series: Z(ζ) ≈ i√π exp(-ζ²) - 2ζ(1 - 2ζ²/3 + 4ζ⁴/15 - ...)
        z_series = cmplx(0.0_real64, sqrt(3.14159265358979323846_real64), real64) * exp(-zeta**2) &
                 - 2.0_real64 * zeta * (1.0_real64 - 2.0_real64*zeta**2/3.0_real64 &
                 + 4.0_real64*zeta**4/15.0_real64)
        
        error = abs(z_val - z_series)
        
        if (error > 0.01_real64) then
            print *, "FAIL: Taylor series expansion incorrect for small argument"
            test_status = test_status + 1
        else
            print *, "PASS: Taylor series expansion correct for small |ζ|"
        end if
        
        ! Test 2: Asymptotic series for large |ζ|
        zeta = cmplx(5.0_real64, 0.0_real64, real64)
        z_val = plasma_z_function(zeta)
        
        ! Asymptotic series: Z(ζ) ≈ -1/ζ(1 + 1/(2ζ²) + 3/(4ζ⁴) + ...)
        z_series = -1.0_real64/zeta * (1.0_real64 + 1.0_real64/(2.0_real64*zeta**2) &
                 + 3.0_real64/(4.0_real64*zeta**4))
        
        error = abs(z_val - z_series) / abs(z_series)
        
        if (error > 0.01_real64) then
            print *, "FAIL: Asymptotic series expansion incorrect for large argument"
            test_status = test_status + 1
        else
            print *, "PASS: Asymptotic series expansion correct for large |ζ|"
        end if
        
    end subroutine test_z_function_series_expansion

end program test_plasma_z_function