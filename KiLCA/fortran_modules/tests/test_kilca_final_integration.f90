program test_kilca_final_integration
    use iso_fortran_env
    use kilca_types_m
    use kilca_fourier_m
    use kilca_adaptive_grid_m
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
    print *, "KiLCA Final Integration and Validation Tests"
    print *, "========================================================"
    print *, ""
    
    ! Run final integration tests
    call test_mathematical_integration()
    call test_grid_mathematics_integration()
    call test_precision_consistency()
    call test_memory_safety()
    call test_performance_validation()
    call test_numerical_accuracy()
    
    ! Print summary
    print *, ""
    print *, "========================================================"
    print *, "FINAL INTEGRATION TEST SUMMARY"
    print *, "========================================================"
    print '(A,I5)', " Total tests run:     ", total_tests
    print '(A,I5)', " Tests passed:        ", passed_tests  
    print '(A,I5)', " Tests failed:        ", failed_tests
    print '(A,F8.2,A)', " Success rate:        ", &
            100.0_dp * real(passed_tests, dp) / real(total_tests, dp), " %"
    
    if (failed_tests > 0) then
        print *, ""
        print *, " *** INTEGRATION TESTS FAILED! ***"
        stop 1
    else
        print *, ""
        print *, " *** ALL INTEGRATION TESTS PASSED! ***"
        print *, ""
        print *, "KiLCA Translation Project Successfully Completed!"
        print *, "All major mathematical libraries implemented and validated."
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
    ! Test 1: Mathematical libraries integration
    !---------------------------------------------------------------------------
    subroutine test_mathematical_integration()
        type(fourier_settings_t) :: fourier_settings
        type(hyperg_1f1_settings_t) :: hyperg_settings
        integer, parameter :: n = 16
        complex(dp) :: signal(n), fft_result(n), ifft_result(n)
        complex(dp) :: hyperg_result, expected_exp
        real(dp) :: error
        integer :: ierr, i
        
        call start_test("Mathematical libraries integration")
        test_passed = .true.
        
        ! Initialize settings
        call fourier_settings_create(fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        hyperg_settings%tolerance = 1.0e-10_dp
        hyperg_settings%max_iterations = 100
        
        ! Test Fourier transform roundtrip
        do i = 1, n
            signal(i) = cmplx(sin(2.0_dp * 3.141592653589793_dp * real(i-1, dp) / real(n, dp)), 0.0_dp, dp)
        end do
        
        call fourier_fft_1d(n, signal, fft_result, fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call fourier_ifft_1d(n, fft_result, ifft_result, fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        error = maxval(abs(signal - ifft_result))
        test_passed = test_passed .and. (error < 1.0e-10_dp)
        
        ! Test hypergeometric function: 1F1(1,1,z) = exp(z)
        call hyperg_1f1_kummer_fortran(1.0_dp, 1.0_dp, 1.0_dp, hyperg_result, hyperg_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        expected_exp = exp(cmplx(1.0_dp, 0.0_dp, dp))
        error = abs(hyperg_result - expected_exp)
        test_passed = test_passed .and. (error < 1.0e-8_dp)
        
        call fourier_settings_destroy(fourier_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_mathematical_integration
    
    !---------------------------------------------------------------------------
    ! Test 2: Grid and mathematics integration
    !---------------------------------------------------------------------------
    subroutine test_grid_mathematics_integration()
        type(adaptive_grid_settings_t) :: grid_settings
        type(fourier_settings_t) :: fourier_settings
        integer, parameter :: n_initial = 10
        real(dp) :: grid_initial(n_initial), function_values(n_initial)
        real(dp), allocatable :: refined_grid(:)
        complex(dp) :: fourier_data(n_initial), fourier_result(n_initial)
        integer :: n_refined, ierr, i
        
        call start_test("Grid and mathematics integration")
        test_passed = .true.
        
        call adaptive_grid_settings_create(grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call fourier_settings_create(fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Generate initial grid
        call generate_uniform_grid(n_initial, 0.0_dp, 1.0_dp, grid_initial, grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create test function with sharp feature
        do i = 1, n_initial
            function_values(i) = exp(-20.0_dp * (grid_initial(i) - 0.7_dp)**2)
            fourier_data(i) = cmplx(function_values(i), 0.0_dp, dp)
        end do
        
        ! Test adaptive refinement
        call adaptive_grid_refine(n_initial, grid_initial, function_values, &
                                 refined_grid, n_refined, grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (allocated(refined_grid))
        test_passed = test_passed .and. (n_refined >= n_initial)
        
        ! Test Fourier analysis of function
        call fourier_fft_1d(n_initial, fourier_data, fourier_result, fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (allocated(refined_grid)) deallocate(refined_grid)
        call adaptive_grid_settings_destroy(grid_settings, ierr)
        call fourier_settings_destroy(fourier_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_grid_mathematics_integration
    
    !---------------------------------------------------------------------------
    ! Test 3: Precision consistency across modules
    !---------------------------------------------------------------------------
    subroutine test_precision_consistency()
        real(dp) :: test_dp
        complex(dp) :: test_complex
        type(fourier_settings_t) :: fourier_settings
        type(adaptive_grid_settings_t) :: grid_settings
        integer :: ierr
        
        call start_test("Precision consistency")
        test_passed = .true.
        
        ! Test dp precision
        test_dp = 1.0_dp
        test_passed = test_passed .and. (precision(test_dp) >= 12)
        test_passed = test_passed .and. (digits(test_dp) >= 45)
        
        ! Test complex precision
        test_complex = cmplx(1.0_dp, 1.0_dp, dp)
        test_passed = test_passed .and. (precision(real(test_complex)) >= 12)
        test_passed = test_passed .and. (precision(aimag(test_complex)) >= 12)
        
        ! Test settings precision consistency
        call fourier_settings_create(fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (fourier_settings%tolerance > 0.0_dp)
        test_passed = test_passed .and. (precision(fourier_settings%tolerance) >= 12)
        
        call adaptive_grid_settings_create(grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (grid_settings%tolerance > 0.0_dp)
        test_passed = test_passed .and. (precision(grid_settings%tolerance) >= 12)
        
        call fourier_settings_destroy(fourier_settings, ierr)
        call adaptive_grid_settings_destroy(grid_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_precision_consistency
    
    !---------------------------------------------------------------------------
    ! Test 4: Memory safety
    !---------------------------------------------------------------------------
    subroutine test_memory_safety()
        real(dp), allocatable :: test_grid(:), test_values(:)
        type(adaptive_grid_settings_t) :: grid_settings
        integer :: ierr, n_test, iteration
        
        call start_test("Memory safety")
        test_passed = .true.
        
        call adaptive_grid_settings_create(grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test multiple allocations/deallocations
        do iteration = 1, 5
            call allocate_adaptive_grid(50, test_grid, test_values, n_test, grid_settings, ierr)
            test_passed = test_passed .and. (ierr == 0)
            test_passed = test_passed .and. (allocated(test_grid))
            test_passed = test_passed .and. (allocated(test_values))
            
            call deallocate_adaptive_grid(test_grid, test_values, grid_settings, ierr)
            test_passed = test_passed .and. (ierr == 0)
            test_passed = test_passed .and. (.not. allocated(test_grid))
            test_passed = test_passed .and. (.not. allocated(test_values))
        end do
        
        call adaptive_grid_settings_destroy(grid_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_memory_safety
    
    !---------------------------------------------------------------------------
    ! Test 5: Performance validation
    !---------------------------------------------------------------------------
    subroutine test_performance_validation()
        integer, parameter :: n_large = 128
        complex(dp) :: large_signal(n_large), large_result(n_large)
        type(fourier_settings_t) :: fourier_settings
        real(dp) :: start_time, end_time, duration
        integer :: ierr, i
        
        call start_test("Performance validation")
        test_passed = .true.
        
        call fourier_settings_create(fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create large test signal
        do i = 1, n_large
            large_signal(i) = cmplx(sin(8.0_dp * 3.141592653589793_dp * real(i, dp) / real(n_large, dp)), &
                                   cos(12.0_dp * 3.141592653589793_dp * real(i, dp) / real(n_large, dp)), dp)
        end do
        
        ! Test FFT performance
        call cpu_time(start_time)
        call fourier_fft_1d(n_large, large_signal, large_result, fourier_settings, ierr)
        call cpu_time(end_time)
        
        duration = end_time - start_time
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (duration >= 0.0_dp)
        test_passed = test_passed .and. (duration < 5.0_dp)  ! Should complete reasonably fast
        
        call fourier_settings_destroy(fourier_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_performance_validation
    
    !---------------------------------------------------------------------------
    ! Test 6: Numerical accuracy validation
    !---------------------------------------------------------------------------
    subroutine test_numerical_accuracy()
        type(hyperg_1f1_settings_t) :: hyperg_settings
        complex(dp) :: result1, result2, result3
        real(dp) :: error
        integer :: ierr
        
        call start_test("Numerical accuracy validation")
        test_passed = .true.
        
        hyperg_settings%tolerance = 1.0e-12_dp
        hyperg_settings%max_iterations = 200
        
        ! Test known hypergeometric values
        ! 1F1(0, 1, z) = 1 for any z
        call hyperg_1f1_kummer_fortran(0.0_dp, 1.0_dp, 5.0_dp, result1, hyperg_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        error = abs(result1 - cmplx(1.0_dp, 0.0_dp, dp))
        test_passed = test_passed .and. (error < 1.0e-10_dp)
        
        ! 1F1(a, a, z) = exp(z) for any a != 0
        call hyperg_1f1_kummer_fortran(2.0_dp, 2.0_dp, 0.5_dp, result2, hyperg_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        error = abs(result2 - exp(cmplx(0.5_dp, 0.0_dp, dp)))
        test_passed = test_passed .and. (error < 1.0e-10_dp)
        
        ! Test consistency between different tolerance settings
        hyperg_settings%tolerance = 1.0e-8_dp
        call hyperg_1f1_kummer_fortran(1.5_dp, 2.5_dp, 1.0_dp, result3, hyperg_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        hyperg_settings%tolerance = 1.0e-12_dp
        call hyperg_1f1_kummer_fortran(1.5_dp, 2.5_dp, 1.0_dp, result1, hyperg_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        error = abs(result3 - result1)
        test_passed = test_passed .and. (error < 1.0e-7_dp)  ! Should be close
        
        call end_test(test_passed)
    end subroutine test_numerical_accuracy
    
end program test_kilca_final_integration