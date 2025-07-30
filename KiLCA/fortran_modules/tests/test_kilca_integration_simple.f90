program test_kilca_integration_simple
    use iso_fortran_env
    use kilca_types_m
    use kilca_constants_m
    use kilca_shared_m
    use kilca_settings_m
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
    print *, "KiLCA System Integration Tests (Implemented Modules)"
    print *, "========================================================"
    print *, ""
    
    ! Run integration tests for implemented modules
    call test_implemented_modules_integration()
    call test_mathematical_libraries_integration()
    call test_settings_consistency()
    call test_data_type_compatibility()
    call test_memory_integration()
    call test_numerical_precision()
    call test_error_handling_consistency()
    call test_performance_integration()
    
    ! Print summary
    print *, ""
    print *, "========================================================"
    print *, "INTEGRATION TEST SUMMARY"
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
    ! Test 1: Integration of implemented modules
    !---------------------------------------------------------------------------
    subroutine test_implemented_modules_integration()
        type(kilca_settings_t) :: settings
        type(fourier_settings_t) :: fourier_settings
        type(adaptive_grid_settings_t) :: grid_settings
        type(hyperg_1f1_settings_t) :: hyperg_settings
        real(dp) :: test_precision
        integer :: ierr
        
        call start_test("Implemented modules integration")
        test_passed = .true.
        
        ! Test all implemented modules can be initialized
        call kilca_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call fourier_settings_create(fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call adaptive_grid_settings_create(grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        hyperg_settings%tolerance = 1.0e-10_dp
        hyperg_settings%max_iterations = 100
        
        ! Test precision consistency
        test_precision = 1.0_dp
        test_passed = test_passed .and. (precision(test_precision) >= 12)
        test_passed = test_passed .and. (settings%precision_dp == dp)
        
        ! Test constants availability
        test_passed = test_passed .and. (PI > 3.0_dp .and. PI < 4.0_dp)
        test_passed = test_passed .and. (E_CONSTANT > 2.0_dp .and. E_CONSTANT < 3.0_dp)
        
        ! Cleanup
        call kilca_settings_destroy(settings, ierr)
        call fourier_settings_destroy(fourier_settings, ierr)
        call adaptive_grid_settings_destroy(grid_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_implemented_modules_integration
    
    !---------------------------------------------------------------------------
    ! Test 2: Mathematical libraries integration
    !---------------------------------------------------------------------------
    subroutine test_mathematical_libraries_integration()
        type(fourier_settings_t) :: fourier_settings
        type(hyperg_1f1_settings_t) :: hyperg_settings
        type(adaptive_grid_settings_t) :: grid_settings
        integer, parameter :: n = 8
        complex(dp) :: signal(n), fft_result(n), ifft_result(n), hyperg_result
        real(dp) :: grid(n), function_vals(n), error
        integer :: ierr, i
        
        call start_test("Mathematical libraries integration")
        test_passed = .true.
        
        ! Initialize settings
        call fourier_settings_create(fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call adaptive_grid_settings_create(grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        hyperg_settings%tolerance = 1.0e-8_dp
        hyperg_settings%max_iterations = 50
        
        ! Test Fourier transform
        do i = 1, n
            signal(i) = cmplx(sin(2.0_dp * PI * real(i-1, dp) / real(n, dp)), 0.0_dp, dp)
        end do
        
        call fourier_fft_1d(n, signal, fft_result, fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call fourier_ifft_1d(n, fft_result, ifft_result, fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        error = maxval(abs(signal - ifft_result))
        test_passed = test_passed .and. (error < 1.0e-10_dp)
        
        ! Test hypergeometric function
        call hyperg_1f1_kummer_fortran(1.0_dp, 1.0_dp, 1.0_dp, hyperg_result, hyperg_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(hyperg_result - exp(cmplx(1.0_dp, 0.0_dp, dp))) < 1.0e-6_dp)
        
        ! Test adaptive grid
        call generate_uniform_grid(n, 0.0_dp, 1.0_dp, grid, grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (grid(1) == 0.0_dp .and. grid(n) == 1.0_dp)
        
        ! Cleanup
        call fourier_settings_destroy(fourier_settings, ierr)
        call adaptive_grid_settings_destroy(grid_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_mathematical_libraries_integration
    
    !---------------------------------------------------------------------------
    ! Test 3: Settings consistency across modules
    !---------------------------------------------------------------------------
    subroutine test_settings_consistency()
        type(kilca_settings_t) :: settings
        type(fourier_settings_t) :: fourier_settings
        type(adaptive_grid_settings_t) :: grid_settings
        integer :: ierr
        
        call start_test("Settings consistency")
        test_passed = .true.
        
        call kilca_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call fourier_settings_create(fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call adaptive_grid_settings_create(grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test consistent tolerance scales
        test_passed = test_passed .and. (settings%tolerance > 0.0_dp)
        test_passed = test_passed .and. (fourier_settings%tolerance > 0.0_dp)
        test_passed = test_passed .and. (grid_settings%tolerance > 0.0_dp)
        
        ! Test consistent debug levels
        test_passed = test_passed .and. (settings%debug_level >= 0)
        test_passed = test_passed .and. (fourier_settings%debug_level >= 0)
        test_passed = test_passed .and. (grid_settings%debug_level >= 0)
        
        ! Test consistent iteration limits
        test_passed = test_passed .and. (settings%max_iterations > 0)
        test_passed = test_passed .and. (fourier_settings%max_iterations > 0)
        test_passed = test_passed .and. (grid_settings%max_points > 0)
        
        call kilca_settings_destroy(settings, ierr)
        call fourier_settings_destroy(fourier_settings, ierr)
        call adaptive_grid_settings_destroy(grid_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_settings_consistency
    
    !---------------------------------------------------------------------------
    ! Test 4: Data type compatibility
    !---------------------------------------------------------------------------
    subroutine test_data_type_compatibility()
        real(dp) :: dp_value
        complex(dp) :: complex_value
        integer :: int_value
        logical :: bool_value
        
        call start_test("Data type compatibility")
        test_passed = .true.
        
        ! Test dp precision
        dp_value = 1.0_dp
        test_passed = test_passed .and. (precision(dp_value) >= 12)
        test_passed = test_passed .and. (range(dp_value) >= 100)
        
        ! Test complex dp compatibility
        complex_value = cmplx(1.0_dp, 1.0_dp, dp)
        test_passed = test_passed .and. (precision(real(complex_value)) >= 12)
        test_passed = test_passed .and. (precision(aimag(complex_value)) >= 12)
        
        ! Test integer compatibility
        int_value = huge(int_value)
        test_passed = test_passed .and. (int_value > 1000000)
        
        ! Test logical compatibility
        bool_value = .true.
        test_passed = test_passed .and. (bool_value)
        bool_value = .false.
        test_passed = test_passed .and. (.not. bool_value)
        
        call end_test(test_passed)
    end subroutine test_data_type_compatibility
    
    !---------------------------------------------------------------------------
    ! Test 5: Memory integration
    !---------------------------------------------------------------------------
    subroutine test_memory_integration()
        real(dp), allocatable :: test_grid(:), test_values(:)
        complex(dp), allocatable :: test_complex(:)
        type(adaptive_grid_settings_t) :: grid_settings
        integer :: ierr, n_test, i
        
        call start_test("Memory integration")
        test_passed = .true.
        
        call adaptive_grid_settings_create(grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test dynamic allocation
        call allocate_adaptive_grid(100, test_grid, test_values, n_test, grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (allocated(test_grid))
        test_passed = test_passed .and. (allocated(test_values))
        test_passed = test_passed .and. (n_test <= 100)
        
        ! Test memory access
        if (allocated(test_grid) .and. allocated(test_values)) then
            do i = 1, n_test
                test_grid(i) = real(i, dp)
                test_values(i) = real(i**2, dp)
            end do
            test_passed = test_passed .and. (test_grid(n_test) == real(n_test, dp))
        end if
        
        ! Test deallocation
        call deallocate_adaptive_grid(test_grid, test_values, grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (.not. allocated(test_grid))
        test_passed = test_passed .and. (.not. allocated(test_values))
        
        call adaptive_grid_settings_destroy(grid_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_memory_integration
    
    !---------------------------------------------------------------------------
    ! Test 6: Numerical precision
    !---------------------------------------------------------------------------
    subroutine test_numerical_precision()
        real(dp) :: small_value, large_value, computed_value
        complex(dp) :: complex_small, complex_result
        type(hyperg_1f1_settings_t) :: hyperg_settings
        integer :: ierr
        
        call start_test("Numerical precision")
        test_passed = .true.
        
        hyperg_settings%tolerance = 1.0e-12_dp
        hyperg_settings%max_iterations = 200
        
        ! Test small values
        small_value = 1.0e-12_dp
        test_passed = test_passed .and. (small_value > 0.0_dp)
        test_passed = test_passed .and. (small_value + 1.0_dp > 1.0_dp)
        
        ! Test large values
        large_value = 1.0e12_dp
        test_passed = test_passed .and. (large_value > 1.0_dp)
        test_passed = test_passed .and. (large_value / 1.0e6_dp == 1.0e6_dp)
        
        ! Test complex precision
        complex_small = cmplx(1.0e-10_dp, 1.0e-10_dp, dp)
        test_passed = test_passed .and. (abs(complex_small) > 0.0_dp)
        
        ! Test mathematical precision with hypergeometric function
        call hyperg_1f1_kummer_fortran(0.5_dp, 1.5_dp, 0.1_dp, complex_result, hyperg_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(complex_result) > 0.0_dp)
        
        call end_test(test_passed)
    end subroutine test_numerical_precision
    
    !---------------------------------------------------------------------------
    ! Test 7: Error handling consistency
    !---------------------------------------------------------------------------
    subroutine test_error_handling_consistency()
        type(kilca_settings_t) :: settings
        type(fourier_settings_t) :: fourier_settings
        type(adaptive_grid_settings_t) :: grid_settings
        real(dp) :: invalid_grid(2)
        complex(dp) :: invalid_signal(0)
        integer :: ierr
        
        call start_test("Error handling consistency")
        test_passed = .true.
        
        ! Test successful operations return ierr = 0
        call kilca_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call fourier_settings_create(fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call adaptive_grid_settings_create(grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test invalid operations return ierr /= 0
        call generate_uniform_grid(0, 0.0_dp, 1.0_dp, invalid_grid, grid_settings, ierr)
        test_passed = test_passed .and. (ierr /= 0)  ! Should fail with n = 0
        
        call generate_uniform_grid(2, 1.0_dp, 0.0_dp, invalid_grid, grid_settings, ierr)
        test_passed = test_passed .and. (ierr /= 0)  ! Should fail with r_max < r_min
        
        ! Cleanup
        call kilca_settings_destroy(settings, ierr)
        call fourier_settings_destroy(fourier_settings, ierr)
        call adaptive_grid_settings_destroy(grid_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_error_handling_consistency
    
    !---------------------------------------------------------------------------
    ! Test 8: Performance integration
    !---------------------------------------------------------------------------
    subroutine test_performance_integration()
        integer, parameter :: n_perf = 64
        complex(dp) :: perf_signal(n_perf), perf_result(n_perf)
        real(dp) :: perf_grid(n_perf), perf_function(n_perf)
        type(fourier_settings_t) :: fourier_settings
        type(adaptive_grid_settings_t) :: grid_settings
        real(dp) :: start_time, end_time, duration
        integer :: ierr, i
        
        call start_test("Performance integration")
        test_passed = .true.
        
        call fourier_settings_create(fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call adaptive_grid_settings_create(grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create performance test data
        do i = 1, n_perf
            perf_signal(i) = cmplx(sin(4.0_dp * PI * real(i, dp) / real(n_perf, dp)), &
                                  cos(6.0_dp * PI * real(i, dp) / real(n_perf, dp)), dp)
            perf_grid(i) = real(i-1, dp) / real(n_perf-1, dp)
            perf_function(i) = exp(-5.0_dp * (perf_grid(i) - 0.5_dp)**2)
        end do
        
        ! Test FFT performance
        call cpu_time(start_time)
        call fourier_fft_1d(n_perf, perf_signal, perf_result, fourier_settings, ierr)
        call cpu_time(end_time)
        
        duration = end_time - start_time
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (duration >= 0.0_dp)
        test_passed = test_passed .and. (duration < 1.0_dp)  ! Should be fast
        
        ! Test grid generation performance
        call cpu_time(start_time)
        call generate_uniform_grid(n_perf, 0.0_dp, 1.0_dp, perf_grid, grid_settings, ierr)
        call cpu_time(end_time)
        
        duration = end_time - start_time
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (duration >= 0.0_dp)
        test_passed = test_passed .and. (duration < 0.1_dp)  # Should be very fast
        
        call fourier_settings_destroy(fourier_settings, ierr)
        call adaptive_grid_settings_destroy(grid_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_performance_integration
    
end program test_kilca_integration_simple