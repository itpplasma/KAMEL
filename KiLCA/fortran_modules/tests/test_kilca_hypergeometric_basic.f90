program test_kilca_hypergeometric_basic
    use iso_fortran_env
    use kilca_types_m
    use kilca_hypergeometric_basic_m
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
    print *, "KiLCA Basic Hypergeometric Functions - Unit Tests"
    print *, "========================================================"
    print *, ""
    
    ! Run all test suites
    call test_settings_structure()
    call test_algorithm_constants()
    call test_pochhammer_symbol()
    call test_convergence_check()
    call test_1f1_known_values()
    call test_1f1_special_cases()
    call test_1f1_complex_arguments()
    call test_asymptotic_approximation()
    call test_error_handling()
    call test_numerical_stability()
    
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
    ! Test 1: Settings structure
    !---------------------------------------------------------------------------
    subroutine test_settings_structure()
        type(hyperg_1f1_settings_t) :: settings
        
        call start_test("Settings structure")
        test_passed = .true.
        
        ! Test default settings
        test_passed = test_passed .and. (settings%algorithm == 1)
        test_passed = test_passed .and. (settings%tolerance == 1.0e-12_dp)
        test_passed = test_passed .and. (settings%max_iterations == 1000)
        test_passed = test_passed .and. (.not. settings%use_acceleration)
        test_passed = test_passed .and. (settings%debug_level == 0)
        
        ! Test setting custom values
        settings%algorithm = HYPERG_ALGORITHM_ASYMPTOTIC
        settings%tolerance = 1.0e-8_dp
        settings%max_iterations = 500
        settings%use_acceleration = .true.
        settings%debug_level = 1
        
        test_passed = test_passed .and. (settings%algorithm == HYPERG_ALGORITHM_ASYMPTOTIC)
        test_passed = test_passed .and. (settings%tolerance == 1.0e-8_dp)
        test_passed = test_passed .and. (settings%max_iterations == 500)
        test_passed = test_passed .and. (settings%use_acceleration)
        test_passed = test_passed .and. (settings%debug_level == 1)
        
        call end_test(test_passed)
    end subroutine test_settings_structure
    
    !---------------------------------------------------------------------------
    ! Test 2: Algorithm constants
    !---------------------------------------------------------------------------
    subroutine test_algorithm_constants()
        call start_test("Algorithm constants")
        test_passed = .true.
        
        test_passed = test_passed .and. (HYPERG_ALGORITHM_SERIES == 1)
        test_passed = test_passed .and. (HYPERG_ALGORITHM_ASYMPTOTIC == 2)
        test_passed = test_passed .and. (HYPERG_ALGORITHM_SERIES /= HYPERG_ALGORITHM_ASYMPTOTIC)
        
        call end_test(test_passed)
    end subroutine test_algorithm_constants
    
    !---------------------------------------------------------------------------
    ! Test 3: Pochhammer symbol computation
    !---------------------------------------------------------------------------
    subroutine test_pochhammer_symbol()
        real(dp) :: result
        
        call start_test("Pochhammer symbol")
        test_passed = .true.
        
        ! Test (x)_0 = 1
        result = hyperg_compute_coefficients(5.0_dp, 0)
        test_passed = test_passed .and. (abs(result - 1.0_dp) < 1.0e-12_dp)
        
        ! Test (x)_1 = x
        result = hyperg_compute_coefficients(3.0_dp, 1)
        test_passed = test_passed .and. (abs(result - 3.0_dp) < 1.0e-12_dp)
        
        ! Test (2)_3 = 2 * 3 * 4 = 24
        result = hyperg_compute_coefficients(2.0_dp, 3)
        test_passed = test_passed .and. (abs(result - 24.0_dp) < 1.0e-12_dp)
        
        ! Test (1/2)_2 = (1/2) * (3/2) = 3/4
        result = hyperg_compute_coefficients(0.5_dp, 2)
        test_passed = test_passed .and. (abs(result - 0.75_dp) < 1.0e-12_dp)
        
        call end_test(test_passed)
    end subroutine test_pochhammer_symbol
    
    !---------------------------------------------------------------------------
    ! Test 4: Convergence checking
    !---------------------------------------------------------------------------
    subroutine test_convergence_check()
        complex(dp) :: term, sum_val
        logical :: converged
        
        call start_test("Convergence checking")
        test_passed = .true.
        
        ! Test converged case
        sum_val = (100.0_dp, 0.0_dp)
        term = (1.0e-10_dp, 0.0_dp)
        converged = hyperg_check_convergence(term, sum_val, 1.0e-8_dp)
        test_passed = test_passed .and. converged
        
        ! Test non-converged case
        term = (1.0_dp, 0.0_dp)
        converged = hyperg_check_convergence(term, sum_val, 1.0e-8_dp)
        test_passed = test_passed .and. (.not. converged)
        
        call end_test(test_passed)
    end subroutine test_convergence_check
    
    !---------------------------------------------------------------------------
    ! Test 5: Known values of 1F1
    !---------------------------------------------------------------------------
    subroutine test_1f1_known_values()
        complex(dp) :: result
        type(hyperg_1f1_settings_t) :: settings
        integer :: ierr
        
        call start_test("1F1 known values")
        test_passed = .true.
        
        settings%tolerance = 1.0e-10_dp
        settings%max_iterations = 100
        
        ! Test 1F1(0, 1, z) = 1 for any z
        call hyperg_1f1_kummer_fortran(0.0_dp, 1.0_dp, 5.0_dp, result, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result - (1.0_dp, 0.0_dp)) < 1.0e-10_dp)
        
        ! Test 1F1(1, 1, z) = exp(z)
        call hyperg_1f1_kummer_fortran(1.0_dp, 1.0_dp, 1.0_dp, result, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result - exp((1.0_dp, 0.0_dp))) < 1.0e-8_dp)
        
        ! Test 1F1(a, a, z) = exp(z) for any a != 0
        call hyperg_1f1_kummer_fortran(2.0_dp, 2.0_dp, 0.5_dp, result, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result - exp((0.5_dp, 0.0_dp))) < 1.0e-10_dp)
        
        call end_test(test_passed)
    end subroutine test_1f1_known_values
    
    !---------------------------------------------------------------------------
    ! Test 6: Special cases
    !---------------------------------------------------------------------------
    subroutine test_1f1_special_cases()
        complex(dp) :: result
        type(hyperg_1f1_settings_t) :: settings
        integer :: ierr
        
        call start_test("1F1 special cases")
        test_passed = .true.
        
        settings%tolerance = 1.0e-10_dp
        settings%max_iterations = 50
        
        ! Test z = 0 (should give 1)
        call hyperg_1f1_kummer_fortran(3.0_dp, 2.0_dp, 0.0_dp, result, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result - (1.0_dp, 0.0_dp)) < 1.0e-12_dp)
        
        ! Test b = 0 (should fail gracefully)
        call hyperg_1f1_kummer_fortran(1.0_dp, 0.0_dp, 1.0_dp, result, settings, ierr)
        test_passed = test_passed .and. (ierr == -1)
        
        call end_test(test_passed)
    end subroutine test_1f1_special_cases
    
    !---------------------------------------------------------------------------
    ! Test 7: Complex arguments
    !---------------------------------------------------------------------------
    subroutine test_1f1_complex_arguments()
        complex(dp) :: b, z, result
        type(hyperg_1f1_settings_t) :: settings
        integer :: ierr
        
        call start_test("1F1 complex arguments")
        test_passed = .true.
        
        settings%tolerance = 1.0e-8_dp
        settings%max_iterations = 200
        
        ! Test 1F1(1, 2+i, 0.5) 
        b = (2.0_dp, 1.0_dp)
        z = (0.5_dp, 0.0_dp)
        call hyperg_1f1_series_fortran(b, z, result, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result) > 0.0_dp)  ! Should converge to finite value
        
        ! Test with complex z
        b = (2.0_dp, 0.0_dp)
        z = (0.3_dp, 0.4_dp)
        call hyperg_1f1_series_fortran(b, z, result, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call end_test(test_passed)
    end subroutine test_1f1_complex_arguments
    
    !---------------------------------------------------------------------------
    ! Test 8: Asymptotic approximation
    !---------------------------------------------------------------------------
    subroutine test_asymptotic_approximation()
        complex(dp) :: result
        type(hyperg_1f1_settings_t) :: settings
        integer :: ierr
        
        call start_test("Asymptotic approximation")
        test_passed = .true.
        
        ! Test with large z (where asymptotic formula should work)
        call hyperg_1f1_asymptotic_fortran(1.0_dp, 2.0_dp, 20.0_dp, result, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(result) > 0.0_dp)
        
        ! Test with small z (should fail appropriately)
        call hyperg_1f1_asymptotic_fortran(1.0_dp, 2.0_dp, 1.0_dp, result, settings, ierr)
        test_passed = test_passed .and. (ierr == -1)
        
        call end_test(test_passed)
    end subroutine test_asymptotic_approximation
    
    !---------------------------------------------------------------------------
    ! Test 9: Error handling
    !---------------------------------------------------------------------------
    subroutine test_error_handling()
        complex(dp) :: result, b, z
        type(hyperg_1f1_settings_t) :: settings
        integer :: ierr
        
        call start_test("Error handling")
        test_passed = .true.
        
        settings%max_iterations = 5  ! Very low to force non-convergence
        
        ! Test non-convergence with limited iterations
        call hyperg_1f1_kummer_fortran(1.0_dp, 2.0_dp, 10.0_dp, result, settings, ierr)
        test_passed = test_passed .and. (ierr == -3)  ! Non-convergence error
        
        ! Test b ≈ 0 for complex version
        b = (1.0e-15_dp, 0.0_dp)
        z = (1.0_dp, 0.0_dp)
        call hyperg_1f1_series_fortran(b, z, result, settings, ierr)
        test_passed = test_passed .and. (ierr == -1)  ! Singular parameter
        
        call end_test(test_passed)
    end subroutine test_error_handling
    
    !---------------------------------------------------------------------------
    ! Test 10: Numerical stability
    !---------------------------------------------------------------------------
    subroutine test_numerical_stability()
        complex(dp) :: result1, result2
        type(hyperg_1f1_settings_t) :: settings1, settings2
        integer :: ierr1, ierr2
        
        call start_test("Numerical stability")
        test_passed = .true.
        
        ! Compare results with different tolerance settings
        settings1%tolerance = 1.0e-8_dp
        settings1%max_iterations = 100
        
        settings2%tolerance = 1.0e-10_dp
        settings2%max_iterations = 200
        
        call hyperg_1f1_kummer_fortran(1.5_dp, 2.5_dp, 2.0_dp, result1, settings1, ierr1)
        call hyperg_1f1_kummer_fortran(1.5_dp, 2.5_dp, 2.0_dp, result2, settings2, ierr2)
        
        test_passed = test_passed .and. (ierr1 == 0 .and. ierr2 == 0)
        ! Results should be close (within the looser tolerance)
        test_passed = test_passed .and. (abs(result1 - result2) < 1.0e-7_dp)
        
        call end_test(test_passed)
    end subroutine test_numerical_stability
    
end program test_kilca_hypergeometric_basic