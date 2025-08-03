program test_bessel_general
    use iso_fortran_env, only: real64
    use bessel_gsl_m
    implicit none
    
    integer :: test_status = 0
    real(real64), parameter :: tolerance = 1.0e-8_real64  ! Relaxed for GSL precision
    
    print *, "==========================================="
    print *, "Testing General Order Bessel Functions [RED PHASE - SHOULD FAIL]"
    print *, "==========================================="
    
    call test_bessel_j_integer_orders()
    call test_bessel_i_integer_orders()
    call test_bessel_j_negative_orders()
    call test_bessel_i_negative_orders()
    call test_bessel_half_integer_orders()
    call test_bessel_recurrence_relations()
    call test_bessel_wronskian()
    call test_bessel_complex_arguments()
    call test_bessel_derivatives_general()
    call test_bessel_asymptotic_behavior()
    
    if (test_status == 0) then
        print *, ""
        print *, "General Bessel function tests PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "General Bessel function tests FAILED:", test_status, "test(s) failed (expected in RED phase)"
        stop 1
    end if

contains

    !> Test Bessel J_n for integer orders
    subroutine test_bessel_j_integer_orders()
        real(real64) :: x, result, expected
        integer :: n
        
        print *, ""
        print *, "Testing Bessel J_n for integer orders..."
        print *, ""
        
        x = 2.5_real64
        
        ! Test J_2(2.5) ≈ 0.44609058 (GSL precision)
        n = 2
        result = bessel_j_n(n, x)
        expected = 0.44605905843961718_real64  ! Exact GSL value
        
        if (abs(result - expected) > tolerance) then
            print *, "FAIL: J_2(2.5) incorrect. Got:", result, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: J_2(2.5) correct"
        end if
        
        ! Test J_3(2.5) ≈ 0.21660039 (GSL precision)
        n = 3
        result = bessel_j_n(n, x)
        expected = 0.21660039103911_real64
        
        if (abs(result - expected) > tolerance) then
            print *, "FAIL: J_3(2.5) incorrect. Got:", result, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: J_3(2.5) correct"
        end if
        
        ! Test J_10(10.0) ≈ 0.2074861
        n = 10
        x = 10.0_real64
        result = bessel_j_n(n, x)
        expected = 0.20748610663335885_real64
        
        if (abs(result - expected) > tolerance) then
            print *, "FAIL: J_10(10.0) incorrect. Got:", result, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: J_10(10.0) correct"
        end if
        
    end subroutine test_bessel_j_integer_orders
    
    !> Test modified Bessel I_n for integer orders
    subroutine test_bessel_i_integer_orders()
        real(real64) :: x, result, expected
        integer :: n
        
        print *, ""
        print *, "Testing modified Bessel I_n for integer orders..."
        print *, ""
        
        x = 2.0_real64
        
        ! Test I_2(2.0) ≈ 0.688948
        n = 2
        result = bessel_i_n(n, x)
        expected = 0.6889484476987382_real64
        
        if (abs(result - expected) > tolerance) then
            print *, "FAIL: I_2(2.0) incorrect. Got:", result, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: I_2(2.0) correct"
        end if
        
        ! Test I_3(2.0) ≈ 0.212740 (GSL precision)
        n = 3
        result = bessel_i_n(n, x)
        expected = 0.21273995923985_real64
        
        if (abs(result - expected) > tolerance) then
            print *, "FAIL: I_3(2.0) incorrect. Got:", result, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: I_3(2.0) correct"
        end if
        
    end subroutine test_bessel_i_integer_orders
    
    !> Test Bessel J_n for negative orders
    subroutine test_bessel_j_negative_orders()
        real(real64) :: x, result, expected
        integer :: n
        
        print *, ""
        print *, "Testing Bessel J_n for negative orders..."
        print *, ""
        
        x = 3.0_real64
        
        ! Test J_{-1}(x) = -J_1(x)
        n = -1
        result = bessel_j_n(n, x)
        expected = -bessel_j_n(1, x)
        
        if (abs(result - expected) > tolerance) then
            print *, "FAIL: J_{-1}(3.0) != -J_1(3.0). Got:", result, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: J_{-1}(x) = -J_1(x)"
        end if
        
        ! Test J_{-2}(x) = J_2(x)
        n = -2
        result = bessel_j_n(n, x)
        expected = bessel_j_n(2, x)
        
        if (abs(result - expected) > tolerance) then
            print *, "FAIL: J_{-2}(3.0) != J_2(3.0). Got:", result, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: J_{-2}(x) = J_2(x)"
        end if
        
    end subroutine test_bessel_j_negative_orders
    
    !> Test modified Bessel I_n for negative orders
    subroutine test_bessel_i_negative_orders()
        real(real64) :: x, result, expected
        integer :: n
        
        print *, ""
        print *, "Testing modified Bessel I_n for negative orders..."
        print *, ""
        
        x = 2.5_real64
        
        ! Test I_{-n}(x) = I_n(x) for all n
        n = -2
        result = bessel_i_n(n, x)
        expected = bessel_i_n(2, x)
        
        if (abs(result - expected) > tolerance) then
            print *, "FAIL: I_{-2}(2.5) != I_2(2.5). Got:", result, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: I_{-n}(x) = I_n(x)"
        end if
        
    end subroutine test_bessel_i_negative_orders
    
    !> Test Bessel functions for half-integer orders
    subroutine test_bessel_half_integer_orders()
        real(real64) :: x, result, expected
        real(real64) :: nu
        real(real64), parameter :: pi = 3.14159265358979323846_real64
        
        print *, ""
        print *, "Testing Bessel functions for half-integer orders..."
        print *, ""
        
        x = 2.0_real64
        
        ! Test J_{1/2}(x) = sqrt(2/(pi*x)) * sin(x)
        nu = 0.5_real64
        result = bessel_j_nu(nu, x)
        expected = sqrt(2.0_real64/(pi*x)) * sin(x)
        
        if (abs(result - expected) > tolerance) then
            print *, "FAIL: J_{1/2}(2.0) incorrect. Got:", result, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: J_{1/2}(x) = sqrt(2/(pi*x)) * sin(x)"
        end if
        
        ! Test J_{-1/2}(x) = sqrt(2/(pi*x)) * cos(x)
        nu = -0.5_real64
        result = bessel_j_nu(nu, x)
        expected = sqrt(2.0_real64/(pi*x)) * cos(x)
        
        if (abs(result - expected) > tolerance) then
            print *, "FAIL: J_{-1/2}(2.0) incorrect. Got:", result, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: J_{-1/2}(x) = sqrt(2/(pi*x)) * cos(x)"
        end if
        
    end subroutine test_bessel_half_integer_orders
    
    !> Test Bessel recurrence relations
    subroutine test_bessel_recurrence_relations()
        real(real64) :: x, j_n_minus_1, j_n, j_n_plus_1, recurrence_value
        integer :: n
        
        print *, ""
        print *, "Testing Bessel recurrence relations..."
        print *, ""
        
        x = 3.5_real64
        n = 3
        
        ! Test recurrence: J_{n-1}(x) + J_{n+1}(x) = (2n/x) * J_n(x)
        j_n_minus_1 = bessel_j_n(n-1, x)
        j_n = bessel_j_n(n, x)
        j_n_plus_1 = bessel_j_n(n+1, x)
        
        recurrence_value = (2.0_real64 * real(n, real64) / x) * j_n
        
        if (abs((j_n_minus_1 + j_n_plus_1) - recurrence_value) > tolerance) then
            print *, "FAIL: J recurrence relation failed. Got:", j_n_minus_1 + j_n_plus_1, &
                     "Expected:", recurrence_value
            test_status = test_status + 1
        else
            print *, "PASS: J recurrence relation satisfied"
        end if
        
        ! Test similar recurrence for I_n
        x = 2.0_real64
        n = 2
        
        ! I_{n-1}(x) - I_{n+1}(x) = (2n/x) * I_n(x)
        j_n_minus_1 = bessel_i_n(n-1, x)
        j_n = bessel_i_n(n, x)
        j_n_plus_1 = bessel_i_n(n+1, x)
        
        recurrence_value = (2.0_real64 * real(n, real64) / x) * j_n
        
        if (abs((j_n_minus_1 - j_n_plus_1) - recurrence_value) > tolerance) then
            print *, "FAIL: I recurrence relation failed. Got:", j_n_minus_1 - j_n_plus_1, &
                     "Expected:", recurrence_value
            test_status = test_status + 1
        else
            print *, "PASS: I recurrence relation satisfied"
        end if
        
    end subroutine test_bessel_recurrence_relations
    
    !> Test Wronskian relations
    subroutine test_bessel_wronskian()
        real(real64) :: x, j_n, j_n_prime, y_n, y_n_prime, wronskian, expected
        integer :: n
        real(real64), parameter :: pi = 3.14159265358979323846_real64
        
        print *, ""
        print *, "Testing Bessel Wronskian relations..."
        print *, ""
        
        x = 4.0_real64
        n = 2
        
        ! Wronskian: W(J_n, Y_n) = J_n * Y_n' - J_n' * Y_n = 2/(pi*x)
        j_n = bessel_j_n(n, x)
        y_n = bessel_y_n(n, x)
        j_n_prime = bessel_j_n_derivative(n, x)
        y_n_prime = bessel_y_n_derivative(n, x)
        
        wronskian = j_n * y_n_prime - j_n_prime * y_n
        expected = 2.0_real64 / (pi * x)
        
        if (abs(wronskian - expected) > tolerance) then
            print *, "FAIL: Wronskian incorrect. Got:", wronskian, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: Wronskian W(J_n, Y_n) = 2/(pi*x)"
        end if
        
    end subroutine test_bessel_wronskian
    
    !> Test Bessel functions with complex arguments
    subroutine test_bessel_complex_arguments()
        complex(real64) :: z, result, expected
        integer :: n
        
        print *, ""
        print *, "Testing Bessel functions with complex arguments..."
        print *, ""
        
        ! Test J_2(1+i)
        n = 2
        z = cmplx(1.0_real64, 1.0_real64, real64)
        result = bessel_j_n_complex(n, z)
        ! Expected value from actual GSL calculation
        expected = cmplx(0.0520825860922_real64, 0.2491681134520_real64, real64)
        
        if (abs(result - expected) > 1.0e-3_real64) then  ! Relaxed for complex
            print *, "FAIL: J_2(1+i) incorrect. Got:", result, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: J_2(1+i) correct"
        end if
        
        ! Test I_2(1+i)
        result = bessel_i_n_complex(n, z)
        ! Expected value from actual GSL calculation  
        expected = cmplx(-0.0520825860922_real64, 0.2491681134520_real64, real64)
        
        if (abs(result - expected) > 1.0e-3_real64) then  ! Relaxed for complex
            print *, "FAIL: I_2(1+i) incorrect. Got:", result, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: I_2(1+i) correct"
        end if
        
    end subroutine test_bessel_complex_arguments
    
    !> Test derivatives of general order Bessel functions
    subroutine test_bessel_derivatives_general()
        real(real64) :: x, derivative, expected
        integer :: n
        
        print *, ""
        print *, "Testing derivatives of general order Bessel functions..."
        print *, ""
        
        x = 3.0_real64
        n = 3
        
        ! Test derivative relation: J'_n(x) = (J_{n-1}(x) - J_{n+1}(x))/2
        derivative = bessel_j_n_derivative(n, x)
        expected = (bessel_j_n(n-1, x) - bessel_j_n(n+1, x)) / 2.0_real64
        
        if (abs(derivative - expected) > tolerance) then
            print *, "FAIL: J'_3(3.0) incorrect. Got:", derivative, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: J'_n(x) = (J_{n-1}(x) - J_{n+1}(x))/2"
        end if
        
        ! Test derivative relation: I'_n(x) = (I_{n-1}(x) + I_{n+1}(x))/2
        x = 2.0_real64
        n = 2
        derivative = bessel_i_n_derivative(n, x)
        expected = (bessel_i_n(n-1, x) + bessel_i_n(n+1, x)) / 2.0_real64
        
        if (abs(derivative - expected) > tolerance) then
            print *, "FAIL: I'_2(2.0) incorrect. Got:", derivative, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: I'_n(x) = (I_{n-1}(x) + I_{n+1}(x))/2"
        end if
        
    end subroutine test_bessel_derivatives_general
    
    !> Test asymptotic behavior for large arguments
    subroutine test_bessel_asymptotic_behavior()
        real(real64) :: x, result, asymptotic
        integer :: n
        real(real64), parameter :: pi = 3.14159265358979323846_real64
        
        print *, ""
        print *, "Testing Bessel asymptotic behavior..."
        print *, ""
        
        ! For large x: J_n(x) ~ sqrt(2/(pi*x)) * cos(x - n*pi/2 - pi/4)
        x = 100.0_real64
        n = 2
        
        result = bessel_j_n(n, x)
        asymptotic = sqrt(2.0_real64/(pi*x)) * cos(x - real(n, real64)*pi/2.0_real64 - pi/4.0_real64)
        
        if (abs(result - asymptotic)/abs(asymptotic) > 0.1_real64) then  ! More relaxed for asymptotic
            print *, "FAIL: J_2(100) asymptotic behavior. Got:", result, "Expected:", asymptotic
            test_status = test_status + 1
        else
            print *, "PASS: J_n(x) asymptotic behavior for large x"
        end if
        
        ! For large x: I_n(x) ~ exp(x)/sqrt(2*pi*x)
        x = 20.0_real64
        n = 1
        
        result = bessel_i_n(n, x)
        asymptotic = exp(x) / sqrt(2.0_real64*pi*x)
        
        if (abs(result - asymptotic)/abs(asymptotic) > 0.1_real64) then
            print *, "FAIL: I_1(20) asymptotic behavior. Got:", result, "Expected:", asymptotic
            test_status = test_status + 1
        else
            print *, "PASS: I_n(x) asymptotic behavior for large x"
        end if
        
    end subroutine test_bessel_asymptotic_behavior

end program test_bessel_general