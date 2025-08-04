program test_bessel_derivatives_correct
    ! Test program for CORRECT Bessel function higher derivatives
    ! This will initially FAIL with current broken implementation
    
    use iso_fortran_env, only: real64, error_unit
    use bessel_gsl_m
    implicit none
    
    integer, parameter :: N_TESTS = 10
    integer :: test_count = 0
    integer :: pass_count = 0
    
    write(*, '(A)') "Testing CORRECT Bessel Function Higher Derivatives"
    write(*, '(A)') "================================================="
    
    ! Test 1: Second derivatives are different from first derivatives
    call test_second_derivative_different()
    
    ! Test 2: Third derivatives are different from first derivatives  
    call test_third_derivative_different()
    
    ! Test 3: Known second derivative values for J_0 and J_1
    call test_known_second_derivatives_j()
    
    ! Test 4: Known second derivative values for I_0 and I_1
    call test_known_second_derivatives_i()
    
    ! Test 5: Derivative recurrence relations for J_n
    call test_j_derivative_recurrence()
    
    ! Test 6: Derivative recurrence relations for I_n
    call test_i_derivative_recurrence()
    
    ! Test 7: Complex argument higher derivatives
    call test_complex_higher_derivatives()
    
    ! Test 8: Derivative chain rule consistency
    call test_derivative_chain_rule()
    
    ! Test 9: Mathematical identities for higher derivatives
    call test_higher_derivative_identities()
    
    ! Test 10: Numerical vs analytical derivative comparison
    call test_numerical_vs_analytical()
    
    ! Summary
    write(*, '(A)') ""
    write(*, '(A, I0, A, I0, A)') "Results: ", pass_count, "/", test_count, " tests passed"
    
    if (pass_count == test_count) then
        write(*, '(A)') "SUCCESS: All Bessel derivative tests passed!"
        stop 0
    else
        write(*, '(A)') "FAILURE: Bessel derivative implementation has shortcuts!"
        stop 1
    end if
    
contains

    subroutine test_second_derivative_different()
        ! Test that second derivatives are NOT equal to first derivatives
        complex(real64) :: z_test = (2.0_real64, 0.5_real64)
        complex(real64) :: j_prime, j_double_prime, i_prime, i_double_prime
        
        call increment_test("Second derivatives different from first")
        
        ! Get first and second derivatives
        j_prime = besselj(1, z_test, 1)          ! First derivative of J_1
        j_double_prime = besselj(1, z_test, 2)   ! Second derivative of J_1
        
        i_prime = besseli(1, z_test, 1)          ! First derivative of I_1  
        i_double_prime = besseli(1, z_test, 2)   ! Second derivative of I_1
        
        ! Current implementation returns same value - this will FAIL
        if (abs(j_prime - j_double_prime) < 1.0e-14_real64) then
            write(*, '(A)') "FAIL: J_1 second derivative returns first derivative (shortcut!)"
            return
        end if
        
        if (abs(i_prime - i_double_prime) < 1.0e-14_real64) then
            write(*, '(A)') "FAIL: I_1 second derivative returns first derivative (shortcut!)"
            return
        end if
        
        call increment_pass("Second derivatives different")
    end subroutine

    subroutine test_third_derivative_different()
        ! Test that third derivatives are NOT equal to first derivatives
        complex(real64) :: z_test = (1.5_real64, 0.0_real64)
        complex(real64) :: j_prime, j_triple_prime
        
        call increment_test("Third derivatives different from first")
        
        j_prime = besselj(0, z_test, 1)          ! First derivative of J_0
        j_triple_prime = besselj(0, z_test, 3)   ! Third derivative of J_0
        
        ! Current implementation returns same value - this will FAIL
        if (abs(j_prime - j_triple_prime) < 1.0e-14_real64) then
            write(*, '(A)') "FAIL: J_0 third derivative returns first derivative (major shortcut!)"
            return
        end if
        
        call increment_pass("Third derivatives different")
    end subroutine

    subroutine test_known_second_derivatives_j()
        ! Test known second derivative values for J_n
        real(real64) :: x
        complex(real64) :: z
        complex(real64) :: j0_val, j_minus1, j_plus1, j0_second_expected, j0_second_actual, j2_val
        
        x = 2.0_real64
        z = cmplx(x, 0.0_real64, real64)
        
        call increment_test("Known J_n second derivative values")
        
        ! J_0''(x) = [J_{-1}(x) - 2*J_0(x) + J_1(x)]/1 = [-J_1(x) - 2*J_0(x) + J_1(x)]/1 = -2*J_0(x)
        ! Actually: J_0''(x) = [J_{-2}(x) - 2*J_0(x) + J_2(x)]/4
        j0_val = besselj(0, z, 0)                  ! J_0(x)
        j_minus1 = -besselj(1, z, 0)              ! J_{-1}(x) = -J_1(x)
        j_plus1 = besselj(1, z, 0)                ! J_1(x)
        
        ! For J_0: J_0''(x) = [J_{-2}(x) - 2*J_0(x) + J_2(x)]/4
        ! But J_{-2}(x) = J_2(x), so: J_0''(x) = [2*J_2(x) - 2*J_0(x)]/4 = [J_2(x) - J_0(x)]/2
        j2_val = besselj(2, z, 0)  ! J_2(x)
        j0_second_expected = (j2_val - j0_val) / 2.0_real64
        
        j0_second_actual = besselj(0, z, 2)        ! Current implementation (wrong)
        
        ! This should be different but current implementation makes them equal
        if (abs(j0_second_expected - j0_second_actual) < 1.0e-12_real64) then
            call increment_pass("J_n second derivatives")
        else
            write(*, '(A)') "FAIL: J_0 second derivative formula incorrect"
            return
        end if
    end subroutine

    subroutine test_known_second_derivatives_i()
        ! Test known second derivative values for I_n  
        real(real64) :: x
        complex(real64) :: z
        complex(real64) :: i0_val, i2_val, i0_second_expected, i0_second_actual
        
        x = 1.5_real64
        z = cmplx(x, 0.0_real64, real64)
        
        call increment_test("Known I_n second derivative values")
        
        ! For I_0: I_0''(x) = [I_{-2}(x) + 2*I_0(x) + I_2(x)]/4 = [I_2(x) + 2*I_0(x) + I_2(x)]/4 = [I_0(x) + I_2(x)]/2
        i0_val = besseli(0, z, 0)                  ! I_0(x)
        i2_val = besseli(2, z, 0)                  ! I_2(x)
        i0_second_expected = (i0_val + i2_val) / 2.0_real64
        
        i0_second_actual = besseli(0, z, 2)        ! Current implementation (wrong)
        
        ! Current implementation will return first derivative, not second
        if (abs(i0_second_expected - i0_second_actual) > 1.0e-12_real64) then
            write(*, '(A)') "FAIL: I_0 second derivative uses wrong formula"
            return
        end if
        
        call increment_pass("I_n second derivatives")
    end subroutine

    subroutine test_j_derivative_recurrence()
        ! Test J_n derivative recurrence relations
        real(real64) :: x
        complex(real64) :: z
        complex(real64) :: j1_prime_formula, j1_prime_direct
        
        x = 2.5_real64
        z = cmplx(x, 0.0_real64, real64)
        
        call increment_test("J_n derivative recurrence relations")
        
        ! J_1'(x) = [J_0(x) - J_2(x)]/2
        j1_prime_formula = (besselj(0, z, 0) - besselj(2, z, 0)) / 2.0_real64
        j1_prime_direct = besselj(1, z, 1)
        
        if (abs(j1_prime_formula - j1_prime_direct) > 1.0e-12_real64) then
            write(*, '(A)') "FAIL: J_1 first derivative recurrence relation incorrect"
            return
        end if
        
        call increment_pass("J_n derivative recurrence")
    end subroutine

    subroutine test_i_derivative_recurrence()
        ! Test I_n derivative recurrence relations
        real(real64) :: x
        complex(real64) :: z
        complex(real64) :: i1_prime_formula, i1_prime_direct
        
        x = 1.8_real64
        z = cmplx(x, 0.0_real64, real64)
        
        call increment_test("I_n derivative recurrence relations")
        
        ! I_1'(x) = [I_0(x) + I_2(x)]/2
        i1_prime_formula = (besseli(0, z, 0) + besseli(2, z, 0)) / 2.0_real64
        i1_prime_direct = besseli(1, z, 1)
        
        if (abs(i1_prime_formula - i1_prime_direct) > 1.0e-12_real64) then
            write(*, '(A)') "FAIL: I_1 first derivative recurrence relation incorrect"
            return
        end if
        
        call increment_pass("I_n derivative recurrence")
    end subroutine

    subroutine test_complex_higher_derivatives()
        ! Test higher derivatives with complex arguments
        complex(real64) :: z = (1.2_real64, 0.8_real64)
        complex(real64) :: j_second, j_third, i_second, i_third
        
        call increment_test("Complex argument higher derivatives")
        
        j_second = besselj(1, z, 2)
        j_third = besselj(1, z, 3)
        i_second = besseli(1, z, 2)
        i_third = besseli(1, z, 3)
        
        ! These should all be different from each other
        if (abs(j_second - j_third) < 1.0e-14_real64) then
            write(*, '(A)') "FAIL: Complex J_1 second and third derivatives identical (shortcut!)"
            return
        end if
        
        if (abs(i_second - i_third) < 1.0e-14_real64) then
            write(*, '(A)') "FAIL: Complex I_1 second and third derivatives identical (shortcut!)"
            return
        end if
        
        call increment_pass("Complex higher derivatives")
    end subroutine

    subroutine test_derivative_chain_rule()
        ! Test that derivatives follow mathematical chain rule
        real(real64) :: x
        complex(real64) :: z
        complex(real64) :: j0_fourth, j0_first
        
        x = 1.0_real64
        z = cmplx(x, 0.0_real64, real64)
        
        call increment_test("Derivative chain rule consistency")
        
        j0_first = besselj(0, z, 1)
        j0_fourth = besselj(0, z, 4)
        
        ! Fourth derivative should definitely be different from first
        if (abs(j0_first - j0_fourth) < 1.0e-14_real64) then
            write(*, '(A)') "FAIL: J_0 fourth derivative equals first derivative (major shortcut!)"
            return
        end if
        
        call increment_pass("Chain rule consistency")
    end subroutine

    subroutine test_higher_derivative_identities()
        ! Test mathematical identities for higher derivatives
        real(real64) :: x
        complex(real64) :: z
        complex(real64) :: j0_val, j0_second, expected_relation
        
        x = 0.1_real64
        z = cmplx(x, 0.0_real64, real64)
        
        call increment_test("Higher derivative mathematical identities")
        
        j0_val = besselj(0, z, 0)
        j0_second = besselj(0, z, 2)
        
        ! For small x: J_0''(x) ≈ -J_0(x)/4 + O(x^2)
        expected_relation = -j0_val / 4.0_real64
        
        ! This is approximate for small x, but should be closer than returning first derivative
        if (abs(j0_second - expected_relation) > abs(j0_second - besselj(0, z, 1))) then
            write(*, '(A)') "FAIL: J_0 second derivative doesn't follow mathematical behavior"
            return
        end if
        
        call increment_pass("Higher derivative identities")
    end subroutine

    subroutine test_numerical_vs_analytical()
        ! Test numerical differentiation vs analytical (when implemented correctly)
        real(real64) :: x, h
        complex(real64) :: z
        complex(real64) :: j0_second_numerical, j0_second_analytical
        complex(real64) :: j0_plus, j0_minus, j0_center
        
        x = 1.0_real64
        h = 1.0e-6_real64
        z = cmplx(x, 0.0_real64, real64)
        
        call increment_test("Numerical vs analytical derivatives")
        
        ! Numerical second derivative: f''(x) ≈ [f(x+h) - 2f(x) + f(x-h)] / h^2
        j0_plus = besselj(0, cmplx(x + h, 0.0_real64, real64), 0)
        j0_center = besselj(0, z, 0)
        j0_minus = besselj(0, cmplx(x - h, 0.0_real64, real64), 0)
        j0_second_numerical = (j0_plus - 2.0_real64 * j0_center + j0_minus) / h**2
        
        j0_second_analytical = besselj(0, z, 2)
        
        ! Current implementation returns first derivative, so this will be very different from numerical
        if (abs(j0_second_analytical - j0_second_numerical) > 0.1_real64) then
            write(*, '(A)') "FAIL: Analytical second derivative very different from numerical (using wrong formula)"
            return
        end if
        
        call increment_pass("Numerical vs analytical")
    end subroutine

    subroutine increment_test(test_name)
        character(len=*), intent(in) :: test_name
        test_count = test_count + 1
        write(*, '(A, I0, A, A, A)', advance='no') "Test ", test_count, ": ", test_name, " ... "
    end subroutine

    subroutine increment_pass(test_name)
        character(len=*), intent(in) :: test_name
        pass_count = pass_count + 1
        write(*, '(A)') "PASS"
    end subroutine

end program test_bessel_derivatives_correct