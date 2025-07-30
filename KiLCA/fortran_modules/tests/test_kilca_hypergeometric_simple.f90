program test_kilca_hypergeometric_simple
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
    print *, "KiLCA Hypergeometric Functions - Simple Tests"
    print *, "========================================================"
    print *, ""
    
    ! Run basic tests (GSL functions only)
    call test_settings_structure()
    call test_algorithm_constants()
    call test_error_handling_basic()
    
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
    ! Test 1: Settings structure functionality
    !---------------------------------------------------------------------------
    subroutine test_settings_structure()
        type(hyperg_1f1_settings_t) :: settings
        
        call start_test("Settings structure")
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
    end subroutine test_settings_structure
    
    !---------------------------------------------------------------------------
    ! Test 2: Algorithm constants
    !---------------------------------------------------------------------------
    subroutine test_algorithm_constants()
        call start_test("Algorithm constants")
        test_passed = .true.
        
        ! Test that all algorithm constants are defined and different
        test_passed = test_passed .and. (HYPERG_ALGORITHM_AUTO == 1)
        test_passed = test_passed .and. (HYPERG_ALGORITHM_KUMMER == 2)
        test_passed = test_passed .and. (HYPERG_ALGORITHM_CONTINUED_FRACTION == 3)
        test_passed = test_passed .and. (HYPERG_ALGORITHM_QUADRATURE == 4)
        
        ! Test that they are all different
        test_passed = test_passed .and. (HYPERG_ALGORITHM_AUTO /= HYPERG_ALGORITHM_KUMMER)
        test_passed = test_passed .and. (HYPERG_ALGORITHM_KUMMER /= HYPERG_ALGORITHM_CONTINUED_FRACTION)
        test_passed = test_passed .and. (HYPERG_ALGORITHM_CONTINUED_FRACTION /= HYPERG_ALGORITHM_QUADRATURE)
        
        call end_test(test_passed)
    end subroutine test_algorithm_constants
    
    !---------------------------------------------------------------------------
    ! Test 3: Basic error handling
    !---------------------------------------------------------------------------
    subroutine test_error_handling_basic()
        complex(dp) :: b, z, result
        type(hyperg_1f1_settings_t) :: settings
        integer :: ierr
        
        call start_test("Basic error handling")
        test_passed = .true.
        
        ! Test invalid algorithm selection
        settings%algorithm = 999  ! Invalid algorithm
        b = (1.0_dp, 0.0_dp)
        z = (0.5_dp, 0.0_dp)
        call hyperg_1f1_custom(b, z, result, settings, ierr)
        test_passed = test_passed .and. (ierr /= 0)
        test_passed = test_passed .and. (result == (0.0_dp, 0.0_dp))
        
        call end_test(test_passed)
    end subroutine test_error_handling_basic
    
end program test_kilca_hypergeometric_simple