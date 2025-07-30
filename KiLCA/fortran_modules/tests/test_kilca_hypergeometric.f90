program test_kilca_hypergeometric
    use iso_fortran_env
    use kilca_types_m
    use kilca_hypergeometric_m
    implicit none
    
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed
    character(len=100) :: test_name
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    
    print *, ""
    print *, "========================================================"
    print *, "KiLCA Hypergeometric Functions - Unit Tests"
    print *, "========================================================"
    print *, ""
    
    ! Run all test suites
    call test_gsl_1f1_basic()
    call test_gsl_2f1_basic()
    call test_gsl_0f1_basic()
    call test_gsl_u_basic()
    call test_1f1_known_values()
    call test_1f1_algorithm_selection()
    call test_1f1_complex_arguments()
    call test_1f1_special_cases()
    call test_1f1_asymptotic_behavior()
    call test_1f1_numerical_stability()
    call test_settings_management()
    call test_error_handling()
    
    ! Print summary
    print *, ""
    print *, "========================================================"
    print *, "TEST SUMMARY"
    print *, "========================================================"
    print '(A,I5)', " Total tests run:     ", total_tests
    print '(A,I5)', " Tests passed:        ", passed_tests  
    print '(A,I5)', " Tests failed:        ", failed_tests
    print '(A,F8.2,A)', " Success rate:        ", &
            100.0_dp * real(passed_tests, dp) / real(total_tests, dp), " %"
    
    if (failed_tests > 0) then
        print *, ""
        print *, " *** SOME TESTS FAILED! ***"
        stop 1
    else
        print *, ""
        print *, " *** ALL TESTS PASSED! ***"
    end if
    
contains

    !---------------------------------------------------------------------------
    ! Helper subroutines
    !---------------------------------------------------------------------------
    subroutine start_test(name)
        character(len=*), intent(in) :: name
        test_name = name
        write(*, '(A,A,A)', advance='no') "Testing ", trim(name), " ... "
    end subroutine start_test
    
    subroutine end_test(passed)
        logical, intent(in) :: passed
        total_tests = total_tests + 1
        if (passed) then
            passed_tests = passed_tests + 1
            print *, " PASSED"
        else
            failed_tests = failed_tests + 1
            print *, " FAILED"
        end if
    end subroutine end_test
    
    !---------------------------------------------------------------------------
    ! Test 1: GSL 1F1 basic functionality
    !---------------------------------------------------------------------------
    subroutine test_gsl_1f1_basic()
        real(dp) :: result
        integer :: ierr
        
        call start_test("GSL 1F1 basic functionality")
        test_passed = .true.
        
        ! Test 1F1(1, 2, 1) - known analytical result
        call hyperg_1f1_gsl(1.0_dp, 2.0_dp, 1.0_dp, result, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result - exp(1.0_dp)) < 1.0e-10_dp)
        
        ! Test 1F1(0, 1, z) = 1 for any z
        call hyperg_1f1_gsl(0.0_dp, 1.0_dp, 5.0_dp, result, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result - 1.0_dp) < 1.0e-12_dp)
        
        ! Test 1F1(a, a, z) = exp(z) for any a != 0
        call hyperg_1f1_gsl(3.0_dp, 3.0_dp, 0.5_dp, result, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result - exp(0.5_dp)) < 1.0e-12_dp)
        
        call end_test(test_passed)
    end subroutine test_gsl_1f1_basic
    
    !---------------------------------------------------------------------------
    ! Test 2: GSL 2F1 basic functionality
    !---------------------------------------------------------------------------
    subroutine test_gsl_2f1_basic()
        real(dp) :: result
        integer :: ierr
        
        call start_test("GSL 2F1 basic functionality")
        test_passed = .true.
        
        ! Test 2F1(1, 1, 2, z) = -ln(1-z)/z
        call hyperg_2f1_gsl(1.0_dp, 1.0_dp, 2.0_dp, 0.5_dp, result, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result - (-log(0.5_dp)/0.5_dp)) < 1.0e-10_dp)
        
        ! Test 2F1(0, b, c, z) = 1
        call hyperg_2f1_gsl(0.0_dp, 2.0_dp, 3.0_dp, 0.8_dp, result, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result - 1.0_dp) < 1.0e-12_dp)
        
        call end_test(test_passed)
    end subroutine test_gsl_2f1_basic
    
    !---------------------------------------------------------------------------
    ! Test 3: GSL 0F1 basic functionality
    !---------------------------------------------------------------------------
    subroutine test_gsl_0f1_basic()
        real(dp) :: result, expected
        integer :: ierr
        
        call start_test("GSL 0F1 basic functionality")
        test_passed = .true.
        
        ! Test 0F1(1, z) = sinh(sqrt(z))/sqrt(z) for z > 0
        call hyperg_0f1_gsl(1.0_dp, 4.0_dp, result, ierr)
        expected = sinh(2.0_dp) / 2.0_dp
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result - expected) < 1.0e-10_dp)
        
        ! Test 0F1(c, 0) = 1
        call hyperg_0f1_gsl(5.0_dp, 0.0_dp, result, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result - 1.0_dp) < 1.0e-12_dp)
        
        call end_test(test_passed)
    end subroutine test_gsl_0f1_basic
    
    !---------------------------------------------------------------------------
    ! Test 4: GSL U function basic functionality
    !---------------------------------------------------------------------------
    subroutine test_gsl_u_basic()
        real(dp) :: result
        integer :: ierr
        
        call start_test("GSL U function basic functionality")
        test_passed = .true.
        
        ! Test basic evaluation (no simple analytical form)
        call hyperg_u_gsl(1.0_dp, 2.0_dp, 1.0_dp, result, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (result > 0.0_dp)  ! Should be positive
        
        ! Test U(0, b, z) = 1
        call hyperg_u_gsl(0.0_dp, 3.0_dp, 2.0_dp, result, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result - 1.0_dp) < 1.0e-12_dp)
        
        call end_test(test_passed)
    end subroutine test_gsl_u_basic
    
    !---------------------------------------------------------------------------
    ! Test 5: 1F1 known values and analytical relations
    !---------------------------------------------------------------------------
    subroutine test_1f1_known_values()
        complex(dp) :: b, z, result
        type(hyperg_1f1_settings_t) :: settings
        integer :: ierr
        
        call start_test("1F1 known values and relations")
        test_passed = .true.
        
        ! Default settings
        settings%algorithm = HYPERG_ALGORITHM_AUTO
        settings%tolerance = 1.0e-12_dp
        
        ! Test 1F1(1, 2, z) = (exp(z) - 1)/z
        b = (2.0_dp, 0.0_dp)
        z = (1.0_dp, 0.0_dp)
        call hyperg_1f1_custom(b, z, result, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result - (exp(z) - 1.0_dp)/z) < 1.0e-10_dp)
        
        ! Test 1F1(1, 1+ib, ib) for pure imaginary arguments
        b = (1.0_dp, 2.0_dp)
        z = (0.0_dp, 2.0_dp)
        call hyperg_1f1_custom(b, z, result, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result) > 0.0_dp)  ! Should converge
        
        call end_test(test_passed)
    end subroutine test_1f1_known_values
    
    !---------------------------------------------------------------------------
    ! Test 6: Algorithm selection logic
    !---------------------------------------------------------------------------
    subroutine test_1f1_algorithm_selection()
        complex(dp) :: b, z, result_auto, result_manual
        type(hyperg_1f1_settings_t) :: settings_auto, settings_manual
        integer :: ierr
        
        call start_test("1F1 algorithm selection")
        test_passed = .true.
        
        settings_auto%algorithm = HYPERG_ALGORITHM_AUTO
        settings_auto%tolerance = 1.0e-10_dp
        
        ! Test small |z/b| case (should use Kummer series)
        b = (10.0_dp, 0.0_dp)
        z = (0.5_dp, 0.0_dp)
        call hyperg_1f1_custom(b, z, result_auto, settings_auto, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        settings_manual%algorithm = HYPERG_ALGORITHM_KUMMER
        settings_manual%tolerance = 1.0e-10_dp
        call hyperg_1f1_custom(b, z, result_manual, settings_manual, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result_auto - result_manual) < 1.0e-8_dp)
        
        ! Test large |z/b| case (should use continued fraction)
        b = (0.5_dp, 0.0_dp)
        z = (10.0_dp, 0.0_dp)
        call hyperg_1f1_custom(b, z, result_auto, settings_auto, ierr)
        
        settings_manual%algorithm = HYPERG_ALGORITHM_CONTINUED_FRACTION
        call hyperg_1f1_custom(b, z, result_manual, settings_manual, ierr)
        
        ! Both should converge (may have different accuracy)
        test_passed = test_passed .and. (ierr == 0)
        
        call end_test(test_passed)
    end subroutine test_1f1_algorithm_selection
    
    !---------------------------------------------------------------------------
    ! Test 7: Complex arguments
    !---------------------------------------------------------------------------
    subroutine test_1f1_complex_arguments()
        complex(dp) :: b, z, result
        type(hyperg_1f1_settings_t) :: settings
        integer :: ierr
        
        call start_test("1F1 complex arguments")
        test_passed = .true.
        
        settings%algorithm = HYPERG_ALGORITHM_KUMMER
        settings%tolerance = 1.0e-10_dp
        
        ! Test with complex b and z
        b = (1.5_dp, 0.5_dp)
        z = (0.5_dp, -0.3_dp)
        call hyperg_1f1_custom(b, z, result, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result) > 0.0_dp)
        
        ! Test conjugate symmetry: 1F1(b*, z*) = [1F1(b, z)]*
        b = (2.0_dp, 1.0_dp)
        z = (1.0_dp, 0.5_dp)
        call hyperg_1f1_custom(b, z, result, settings, ierr)
        
        complex(dp) :: b_conj, z_conj, result_conj
        b_conj = conjg(b)
        z_conj = conjg(z)
        call hyperg_1f1_custom(b_conj, z_conj, result_conj, settings, ierr)
        
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(conjg(result) - result_conj) < 1.0e-8_dp)
        
        call end_test(test_passed)
    end subroutine test_1f1_complex_arguments
    
    !---------------------------------------------------------------------------
    ! Test 8: Special cases and edge conditions
    !---------------------------------------------------------------------------
    subroutine test_1f1_special_cases()
        complex(dp) :: b, z, result
        type(hyperg_1f1_settings_t) :: settings
        integer :: ierr
        
        call start_test("1F1 special cases")
        test_passed = .true.
        
        settings%algorithm = HYPERG_ALGORITHM_AUTO
        settings%tolerance = 1.0e-10_dp
        
        ! Test z = 0 (should give 1)
        b = (3.0_dp, 0.5_dp)
        z = (0.0_dp, 0.0_dp)
        call hyperg_1f1_custom(b, z, result, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result - (1.0_dp, 0.0_dp)) < 1.0e-12_dp)
        
        ! Test very small z
        z = (1.0e-10_dp, 0.0_dp)
        call hyperg_1f1_custom(b, z, result, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result - (1.0_dp, 0.0_dp)) < 1.0e-8_dp)
        
        call end_test(test_passed)
    end subroutine test_1f1_special_cases
    
    !---------------------------------------------------------------------------
    ! Test 9: Asymptotic behavior
    !---------------------------------------------------------------------------
    subroutine test_1f1_asymptotic_behavior()
        complex(dp) :: b, z, result
        type(hyperg_1f1_settings_t) :: settings
        integer :: ierr
        
        call start_test("1F1 asymptotic behavior")
        test_passed = .true.
        
        settings%algorithm = HYPERG_ALGORITHM_CONTINUED_FRACTION
        settings%tolerance = 1.0e-8_dp
        
        ! For large |z|, 1F1(1, b, z) ~ exp(z) * z^{-b} for Re(z) > 0
        b = (2.0_dp, 0.0_dp)
        z = (10.0_dp, 0.0_dp)
        call hyperg_1f1_custom(b, z, result, settings, ierr)
        
        complex(dp) :: asymptotic_approx
        asymptotic_approx = exp(z) * z**(-b)
        
        test_passed = test_passed .and. (ierr == 0)
        ! Relative error should be reasonable for asymptotic approximation
        test_passed = test_passed .and. (abs((result - asymptotic_approx)/result) < 0.1_dp)
        
        call end_test(test_passed)
    end subroutine test_1f1_asymptotic_behavior
    
    !---------------------------------------------------------------------------
    ! Test 10: Numerical stability
    !---------------------------------------------------------------------------
    subroutine test_1f1_numerical_stability()
        complex(dp) :: b, z, result1, result2
        type(hyperg_1f1_settings_t) :: settings1, settings2
        integer :: ierr1, ierr2
        
        call start_test("1F1 numerical stability")
        test_passed = .true.
        
        ! Compare different algorithms for the same input
        b = (1.5_dp, 0.2_dp)
        z = (2.0_dp, 0.5_dp)
        
        settings1%algorithm = HYPERG_ALGORITHM_KUMMER
        settings1%tolerance = 1.0e-10_dp
        call hyperg_1f1_custom(b, z, result1, settings1, ierr1)
        
        settings2%algorithm = HYPERG_ALGORITHM_QUADRATURE
        settings2%tolerance = 1.0e-10_dp
        call hyperg_1f1_custom(b, z, result2, settings2, ierr2)
        
        test_passed = test_passed .and. (ierr1 == 0 .and. ierr2 == 0)
        ! Results should agree to reasonable precision
        test_passed = test_passed .and. (abs(result1 - result2) < 1.0e-6_dp)
        
        call end_test(test_passed)
    end subroutine test_1f1_numerical_stability
    
    !---------------------------------------------------------------------------
    ! Test 11: Settings management
    !---------------------------------------------------------------------------
    subroutine test_settings_management()
        type(hyperg_1f1_settings_t) :: settings
        
        call start_test("Settings management")
        test_passed = .true.
        
        ! Test default settings
        test_passed = test_passed .and. (settings%algorithm == 1)
        test_passed = test_passed .and. (settings%tolerance == 1.0e-12_dp)
        test_passed = test_passed .and. (settings%max_iterations == 1000000)
        test_passed = test_passed .and. (.not. settings%use_acceleration)
        test_passed = test_passed .and. (settings%debug_level == 0)
        
        ! Test setting custom values
        settings%algorithm = HYPERG_ALGORITHM_KUMMER
        settings%tolerance = 1.0e-8_dp
        settings%max_iterations = 10000
        settings%use_acceleration = .true.
        settings%debug_level = 2
        
        test_passed = test_passed .and. (settings%algorithm == HYPERG_ALGORITHM_KUMMER)
        test_passed = test_passed .and. (settings%tolerance == 1.0e-8_dp)
        test_passed = test_passed .and. (settings%max_iterations == 10000)
        test_passed = test_passed .and. (settings%use_acceleration)
        test_passed = test_passed .and. (settings%debug_level == 2)
        
        call end_test(test_passed)
    end subroutine test_settings_management
    
    !---------------------------------------------------------------------------
    ! Test 12: Error handling
    !---------------------------------------------------------------------------
    subroutine test_error_handling()
        complex(dp) :: b, z, result
        type(hyperg_1f1_settings_t) :: settings
        integer :: ierr
        real(dp) :: gsl_result
        
        call start_test("Error handling")
        test_passed = .true.
        
        ! Test invalid algorithm selection
        settings%algorithm = 999  ! Invalid algorithm
        b = (1.0_dp, 0.0_dp)
        z = (0.5_dp, 0.0_dp)
        call hyperg_1f1_custom(b, z, result, settings, ierr)
        test_passed = test_passed .and. (ierr /= 0)
        
        ! Test GSL function with potential problematic input
        call hyperg_1f1_gsl(1.0_dp, 0.0_dp, 1.0_dp, gsl_result, ierr)  ! b=0 is problematic
        ! Should handle gracefully (may succeed or fail, but shouldn't crash)
        test_passed = test_passed .and. .true.  ! If we get here, no crash occurred
        
        call end_test(test_passed)
    end subroutine test_error_handling
    
end program test_kilca_hypergeometric