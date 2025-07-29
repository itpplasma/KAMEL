program test_kilca_interp_integ
    use iso_fortran_env, only: real64, int32, error_unit
    use kilca_interp_integ_m
    implicit none
    
    ! Test counters
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    
    print *, "======================================================"
    print *, "KiLCA Interpolation/Integration Library - Unit Tests"
    print *, "======================================================"
    print *, ""
    
    ! Test 1: Neville polynomial interpolation
    call test_neville_polynomial_interpolation()
    
    ! Test 2: Lagrange polynomial interpolation  
    call test_lagrange_polynomial_interpolation()
    
    ! Test 3: Binary search algorithm
    call test_binary_search_algorithm()
    
    ! Test 4: GSL integration replacement - trapezoidal rule
    call test_trapezoidal_integration()
    
    ! Test 5: GSL integration replacement - adaptive quadrature
    call test_adaptive_quadrature()
    
    ! Test 6: Cylindrical volume integration
    call test_cylindrical_volume_integration()
    
    ! Test 7: Complex function interpolation
    call test_complex_interpolation()
    
    ! Test 8: Derivative calculation via interpolation
    call test_derivative_interpolation()
    
    ! Final summary
    print *, ""
    print *, "======================================================"
    print *, "TEST SUMMARY"
    print *, "======================================================"
    print *, "Total tests run:    ", total_tests
    print *, "Tests passed:       ", passed_tests
    print *, "Tests failed:       ", failed_tests
    print *, "Success rate:       ", real(passed_tests)/real(total_tests)*100.0, "%"
    print *, ""
    
    if (failed_tests == 0) then
        print *, "*** ALL TESTS PASSED! ***"
    else
        print *, "*** SOME TESTS FAILED! ***"
        stop 1
    end if
    
contains

    !---------------------------------------------------------------------------
    ! Test 1: Neville polynomial interpolation
    !---------------------------------------------------------------------------
    subroutine test_neville_polynomial_interpolation()
        real(real64) :: x(5), y(5), xeval, result(2), expected
        real(real64) :: error
        integer :: i, ind, ierr
        
        call start_test("Neville polynomial interpolation")
        
        ! Create test data - quadratic function y = x^2 + 2x + 1
        do i = 1, 5
            x(i) = real(i - 3, real64)  ! x = [-2, -1, 0, 1, 2]
            y(i) = x(i)**2 + 2.0_real64*x(i) + 1.0_real64
        end do
        
        ! Test interpolation at x = 0.5
        xeval = 0.5_real64
        expected = xeval**2 + 2.0_real64*xeval + 1.0_real64  ! = 2.25
        
        ind = 2  ! Starting index
        call eval_neville_polynom(5, x, y, 4, xeval, 0, 1, ind, result, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            error = abs(result(1) - expected)
            test_passed = (error < 1.0e-12_real64)
        end if
        
        if (.not. test_passed) then
            print *, "  Neville interpolation failed, ierr =", ierr, ", error =", error
        else
            print *, "  Neville interpolation successful, error =", error
        end if
        
        call end_test(test_passed)
    end subroutine test_neville_polynomial_interpolation
    
    !---------------------------------------------------------------------------
    ! Test 2: Lagrange polynomial interpolation
    !---------------------------------------------------------------------------
    subroutine test_lagrange_polynomial_interpolation()
        real(real64) :: x(4), y(4), xeval, result, expected, error
        integer :: i, ind, ierr
        
        call start_test("Lagrange polynomial interpolation")
        
        ! Create test data - cubic function y = x^3 - 2x + 3
        do i = 1, 4
            x(i) = real(i - 2, real64)  ! x = [-1, 0, 1, 2]
            y(i) = x(i)**3 - 2.0_real64*x(i) + 3.0_real64
        end do
        
        ! Test interpolation at x = 0.5
        xeval = 0.5_real64
        expected = xeval**3 - 2.0_real64*xeval + 3.0_real64  ! = 2.125
        
        ind = 1  ! Starting index
        call eval_lagrange_polynom(4, x, y, 3, xeval, 0, 0, ind, result, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            error = abs(result - expected)
            test_passed = (error < 1.0e-12_real64)
        end if
        
        if (.not. test_passed) then
            print *, "  Lagrange interpolation failed, ierr =", ierr, ", error =", error
        else
            print *, "  Lagrange interpolation successful, error =", error
        end if
        
        call end_test(test_passed)
    end subroutine test_lagrange_polynomial_interpolation
    
    !---------------------------------------------------------------------------
    ! Test 3: Binary search algorithm
    !---------------------------------------------------------------------------
    subroutine test_binary_search_algorithm()
        real(real64) :: x(10)
        real(real64) :: target
        integer :: i, result_index, expected_index
        
        call start_test("Binary search algorithm")
        
        ! Create sorted array
        do i = 1, 10
            x(i) = real(i, real64) * 2.0_real64  ! x = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
        end do
        
        ! Test search for value that should be at index 4 (between 8 and 10)
        target = 9.0_real64
        expected_index = 4  ! Should find index 4 (x[4] = 8.0 < target < x[5] = 10.0)
        
        call search_array(target, 10, x, result_index)
        
        test_passed = (result_index == expected_index)
        
        if (.not. test_passed) then
            print *, "  Binary search failed, expected index =", expected_index, &
                     ", got =", result_index
        else
            print *, "  Binary search successful, found correct index =", result_index
        end if
        
        call end_test(test_passed)
    end subroutine test_binary_search_algorithm
    
    !---------------------------------------------------------------------------
    ! Test 4: Trapezoidal integration
    !---------------------------------------------------------------------------
    subroutine test_trapezoidal_integration()
        real(real64) :: x(101), y(101), result, expected, error
        integer :: i, ierr
        
        call start_test("Trapezoidal integration rule")
        
        ! Create test data - integrate x^2 from 0 to 1, expected = 1/3
        do i = 1, 101
            x(i) = real(i - 1, real64) / 100.0_real64
            y(i) = x(i)**2
        end do
        
        expected = 1.0_real64 / 3.0_real64
        
        call integrate_trapezoid(101, x, y, result, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            error = abs(result - expected)
            test_passed = (error < 1.0e-4_real64)  ! Reasonable tolerance for trapezoidal rule
        end if
        
        if (.not. test_passed) then
            print *, "  Trapezoidal integration failed, ierr =", ierr, ", error =", error
        else
            print *, "  Trapezoidal integration successful, error =", error
        end if
        
        call end_test(test_passed)
    end subroutine test_trapezoidal_integration
    
    !---------------------------------------------------------------------------
    ! Test 5: Adaptive quadrature (simplified version)
    !---------------------------------------------------------------------------
    subroutine test_adaptive_quadrature()
        real(real64) :: result, expected, error
        real(real64) :: a, b, epsabs, epsrel
        integer :: ierr
        
        call start_test("Adaptive quadrature integration")
        
        ! Test integrating simple function
        a = 0.0_real64
        b = 1.0_real64
        epsabs = 1.0e-8_real64
        epsrel = 1.0e-8_real64
        expected = 1.0_real64 / 3.0_real64  ! Integral of x^2 from 0 to 1
        
        call integrate_adaptive_quadrature(test_function_x_squared, a, b, &
                                           epsabs, epsrel, result, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            error = abs(result - expected)
            test_passed = (error < 1.0e-6_real64)
        end if
        
        if (.not. test_passed) then
            print *, "  Adaptive quadrature failed, ierr =", ierr, ", error =", error
        else
            print *, "  Adaptive quadrature successful, error =", error
        end if
        
        call end_test(test_passed)
    end subroutine test_adaptive_quadrature
    
    !---------------------------------------------------------------------------
    ! Test 6: Cylindrical volume integration
    !---------------------------------------------------------------------------
    subroutine test_cylindrical_volume_integration()
        real(real64) :: r(51), f(51), vol_fac, result(51), expected, error
        integer :: i, ierr
        
        call start_test("Cylindrical volume integration")
        
        ! Create radial grid and test function f(r) = r (linear function)
        do i = 1, 51
            r(i) = real(i - 1, real64) / 50.0_real64  ! r from 0 to 1
            f(i) = r(i)
        end do
        
        vol_fac = 2.0_real64 * 3.14159265359_real64  ! 2π for cylinder integration
        
        call integrate_over_cylinder(51, r, f, vol_fac, result, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            ! For f(r) = r, integrated over cylinder from 0 to 1: ∫₀¹ r·r dr = ∫₀¹ r² dr = 1/3
            ! With vol_fac = 2π: result should be 2π/3
            expected = vol_fac / 3.0_real64
            error = abs(result(51) - expected)
            test_passed = (error < 1.0e-3_real64)  ! Reasonable tolerance for cylindrical integration
        end if
        
        if (.not. test_passed) then
            print *, "  Cylindrical integration failed, ierr =", ierr, ", error =", error
        else
            print *, "  Cylindrical integration successful, error =", error
        end if
        
        call end_test(test_passed)
    end subroutine test_cylindrical_volume_integration
    
    !---------------------------------------------------------------------------
    ! Test 7: Complex function interpolation
    !---------------------------------------------------------------------------
    subroutine test_complex_interpolation()
        real(real64) :: x(5)
        complex(real64) :: y(5), xeval_cmplx, result, expected, error_cmplx
        real(real64) :: xeval, error
        integer :: i, ind, ierr
        
        call start_test("Complex function interpolation")
        
        ! Create test data - complex function y = (1+i)*x^2
        do i = 1, 5
            x(i) = real(i - 3, real64)  ! x = [-2, -1, 0, 1, 2]
            y(i) = cmplx(1.0_real64, 1.0_real64, real64) * x(i)**2
        end do
        
        ! Test interpolation at x = 0.5
        xeval = 0.5_real64
        expected = cmplx(1.0_real64, 1.0_real64, real64) * xeval**2
        
        ind = 2  ! Starting index
        call eval_complex_interp(5, x, y, 4, xeval, ind, result, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            error_cmplx = result - expected
            error = abs(error_cmplx)
            test_passed = (error < 1.0e-12_real64)
        end if
        
        if (.not. test_passed) then
            print *, "  Complex interpolation failed, ierr =", ierr, ", error =", error
        else
            print *, "  Complex interpolation successful, error =", error
        end if
        
        call end_test(test_passed)
    end subroutine test_complex_interpolation
    
    !---------------------------------------------------------------------------
    ! Test 8: Derivative calculation via interpolation
    !---------------------------------------------------------------------------
    subroutine test_derivative_interpolation()
        real(real64) :: x(5), y(5), xeval, result(2), expected, error
        integer :: i, ind, ierr
        
        call start_test("Derivative calculation via interpolation")
        
        ! Create test data - cubic function y = x^3, dy/dx = 3x^2
        do i = 1, 5
            x(i) = real(i - 3, real64)  ! x = [-2, -1, 0, 1, 2]
            y(i) = x(i)**3
        end do
        
        ! Test derivative at x = 1.5, expected dy/dx = 3*(1.5)^2 = 6.75
        xeval = 1.5_real64
        expected = 3.0_real64 * xeval**2
        
        ind = 3  ! Starting index
        call eval_neville_polynom(5, x, y, 4, xeval, 0, 1, ind, result, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            error = abs(result(2) - expected)  ! result(2) is the derivative
            test_passed = (error < 1.0e-10_real64)
        end if
        
        if (.not. test_passed) then
            print *, "  Derivative interpolation failed, ierr =", ierr, ", error =", error
        else
            print *, "  Derivative interpolation successful, error =", error
        end if
        
        call end_test(test_passed)
    end subroutine test_derivative_interpolation
    
    !---------------------------------------------------------------------------
    ! Test function for quadrature (x^2)
    !---------------------------------------------------------------------------
    function test_function_x_squared(x) result(y)
        real(real64), intent(in) :: x
        real(real64) :: y
        y = x**2
    end function test_function_x_squared
    
    !---------------------------------------------------------------------------
    ! Test utilities
    !---------------------------------------------------------------------------
    subroutine start_test(name)
        character(len=*), intent(in) :: name
        total_tests = total_tests + 1
        write(*,'(A,A)', advance='no') "Testing ", name
        write(*,'(A)', advance='no') " ... "
    end subroutine start_test
    
    subroutine end_test(passed)
        logical, intent(in) :: passed
        if (passed) then
            print *, "PASSED"
            passed_tests = passed_tests + 1
        else
            print *, "FAILED"
            failed_tests = failed_tests + 1
        end if
    end subroutine end_test


end program test_kilca_interp_integ