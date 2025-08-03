program test_bessel_functions
    use iso_fortran_env, only: real64
    use bessel_gsl_m, only: besselj, besseli
    implicit none
    
    integer :: test_status = 0
    real(real64), parameter :: tolerance = 1.0e-8_real64
    
    print *, "===========================================" 
    print *, "Testing Bessel Function Implementations [GREEN PHASE]"
    print *, "==========================================="
    
    call test_bessel_function_accuracy()
    
    if (test_status == 0) then
        print *, ""
        print *, "Bessel function tests PASSED - GREEN phase successful"
        stop 0
    else
        print *, ""
        print *, "Bessel function tests FAILED:", test_status, "test(s) failed"
        stop 1
    end if

contains

    !> Test Bessel function accuracy against known values
    subroutine test_bessel_function_accuracy()
        complex(real64) :: result, expected
        real(real64) :: error
        
        print *, "Testing Bessel function accuracy against known mathematical values..."
        print *, ""
        
        ! Test J_0 at specific points with known exact values
        print *, "Testing J_0 (Bessel function of first kind, order 0):"
        
        ! J_0(0) = 1 (exact)
        result = besselj(0, cmplx(0.0_real64, 0.0_real64, real64), 0)
        expected = cmplx(1.0_real64, 0.0_real64, real64)
        call check_result("J_0(0)", result, expected)
        
        ! J_0(1) ≈ 0.7651976865579666 (high precision)
        result = besselj(0, cmplx(1.0_real64, 0.0_real64, real64), 0)
        expected = cmplx(0.7651976865579666_real64, 0.0_real64, real64)
        call check_result("J_0(1)", result, expected)
        
        ! J_0(5) ≈ -0.17759677131433830 (high precision)
        result = besselj(0, cmplx(5.0_real64, 0.0_real64, real64), 0)
        expected = cmplx(-0.17759677131433830_real64, 0.0_real64, real64)
        call check_result("J_0(5)", result, expected)
        
        print *, ""
        print *, "Testing J_1 (Bessel function of first kind, order 1):"
        
        ! J_1(0) = 0 (exact)
        result = besselj(1, cmplx(0.0_real64, 0.0_real64, real64), 0)
        expected = cmplx(0.0_real64, 0.0_real64, real64)
        call check_result("J_1(0)", result, expected)
        
        ! J_1(1) ≈ 0.4400505857449335 (high precision)
        result = besselj(1, cmplx(1.0_real64, 0.0_real64, real64), 0)
        expected = cmplx(0.4400505857449335_real64, 0.0_real64, real64)
        call check_result("J_1(1)", result, expected)
        
        ! J_1(2) ≈ 0.5767248077568734 (high precision)
        result = besselj(1, cmplx(2.0_real64, 0.0_real64, real64), 0)
        expected = cmplx(0.5767248077568734_real64, 0.0_real64, real64)
        call check_result("J_1(2)", result, expected)
        
        print *, ""
        print *, "Testing I_0 (Modified Bessel function of first kind, order 0):"
        
        ! I_0(0) = 1 (exact)
        result = besseli(0, cmplx(0.0_real64, 0.0_real64, real64), 0)
        expected = cmplx(1.0_real64, 0.0_real64, real64)
        call check_result("I_0(0)", result, expected)
        
        ! I_0(1) ≈ 1.2660658777520083 (high precision)
        result = besseli(0, cmplx(1.0_real64, 0.0_real64, real64), 0)
        expected = cmplx(1.2660658777520083_real64, 0.0_real64, real64)
        call check_result("I_0(1)", result, expected)
        
        ! I_0(2) ≈ 2.2795853023360673 (high precision)
        result = besseli(0, cmplx(2.0_real64, 0.0_real64, real64), 0)
        expected = cmplx(2.2795853023360673_real64, 0.0_real64, real64)
        call check_result("I_0(2)", result, expected)
        
        print *, ""
        print *, "Testing I_1 (Modified Bessel function of first kind, order 1):"
        
        ! I_1(0) = 0 (exact)
        result = besseli(1, cmplx(0.0_real64, 0.0_real64, real64), 0)
        expected = cmplx(0.0_real64, 0.0_real64, real64)
        call check_result("I_1(0)", result, expected)
        
        ! I_1(1) ≈ 0.5651591039924851 (high precision)
        result = besseli(1, cmplx(1.0_real64, 0.0_real64, real64), 0)
        expected = cmplx(0.5651591039924851_real64, 0.0_real64, real64)
        call check_result("I_1(1)", result, expected)
        
        ! I_1(2) ≈ 1.5906368546373291 (high precision)
        result = besseli(1, cmplx(2.0_real64, 0.0_real64, real64), 0)
        expected = cmplx(1.5906368546373291_real64, 0.0_real64, real64)
        call check_result("I_1(2)", result, expected)
        
        print *, ""
        print *, "Testing complex arguments:"
        
        ! J_0(i) = I_0(1) ≈ 1.2660658777520083 (mathematical identity)
        result = besselj(0, cmplx(0.0_real64, 1.0_real64, real64), 0)
        expected = cmplx(1.2660658777520083_real64, 0.0_real64, real64)
        call check_result("J_0(i)", result, expected)
        
        ! Test higher order functions
        print *, ""
        print *, "Testing higher order Bessel functions:"
        
        ! J_2(1) ≈ 0.11490348493190047 (high precision)
        result = besselj(2, cmplx(1.0_real64, 0.0_real64, real64), 0)
        expected = cmplx(0.11490348493190047_real64, 0.0_real64, real64)
        call check_result("J_2(1)", result, expected)
        
        ! Test derivatives (n = 1 means first derivative)
        print *, ""
        print *, "Testing Bessel function derivatives:"
        
        ! d/dz J_0(z) = -J_1(z), so J_0'(1) = -J_1(1) ≈ -0.4400505857449335
        result = besselj(0, cmplx(1.0_real64, 0.0_real64, real64), 1)
        expected = cmplx(-0.4400505857449335_real64, 0.0_real64, real64)
        call check_result("J_0'(1)", result, expected)
        
        print *, ""
        print *, "Summary: All tests should PASS in GREEN phase with proper GSL implementation"
        
    end subroutine test_bessel_function_accuracy
    
    !> Check if computed result matches expected value within tolerance
    subroutine check_result(test_name, result, expected)
        character(len=*), intent(in) :: test_name
        complex(real64), intent(in) :: result, expected
        real(real64) :: error
        
        error = abs(result - expected)
        
        if (error < tolerance) then
            print *, "PASS: ", test_name, " = ", result, " (expected ", expected, ")"
        else
            print *, "FAIL: ", test_name, " = ", result, " (expected ", expected, ", error = ", error, ")"
            test_status = test_status + 1
        end if
        
    end subroutine check_result

end program test_bessel_functions