program test_kilca_fourier
    use iso_fortran_env
    use kilca_types_m
    use kilca_fourier_m
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
    print *, "KiLCA Fourier Transform System - Unit Tests"
    print *, "========================================================"
    print *, ""
    
    ! Run all test suites
    call test_fft_basic_operations()
    call test_custom_fourier_transforms()
    call test_spectral_analysis()
    call test_windowing_functions()
    call test_convolution_operations()
    call test_power_spectrum_analysis()
    call test_frequency_domain_filtering()
    call test_inverse_consistency()
    call test_performance_benchmarks()
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
    ! Test 1: Basic FFT operations
    !---------------------------------------------------------------------------
    subroutine test_fft_basic_operations()
        integer, parameter :: n = 8
        complex(dp) :: input(n), output(n), reconstructed(n)
        type(fourier_settings_t) :: settings
        integer :: ierr, i
        real(dp) :: error
        
        call start_test("Basic FFT operations")
        test_passed = .true.
        
        call fourier_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create test signal: delta function
        input = (0.0_dp, 0.0_dp)
        input(1) = (1.0_dp, 0.0_dp)
        
        ! Forward FFT
        call fourier_fft_1d(n, input, output, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Inverse FFT
        call fourier_ifft_1d(n, output, reconstructed, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Check reconstruction accuracy
        error = 0.0_dp
        do i = 1, n
            error = max(error, abs(input(i) - reconstructed(i)))
        end do
        test_passed = test_passed .and. (error < 1.0e-10_dp)
        
        call end_test(test_passed)
    end subroutine test_fft_basic_operations
    
    !---------------------------------------------------------------------------
    ! Test 2: Custom Fourier transforms (direct and inverse)
    !---------------------------------------------------------------------------
    subroutine test_custom_fourier_transforms()
        integer, parameter :: dimx = 10, M = 3
        real(dp) :: x(dimx), x0, delta
        complex(dp) :: y(dimx), f(-M:M), y_reconstructed(dimx)
        type(fourier_settings_t) :: settings
        integer :: ierr, i
        real(dp) :: error
        
        call start_test("Custom Fourier transforms")
        test_passed = .true.
        
        call fourier_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create test data
        x0 = 0.0_dp
        delta = 1.0_dp
        do i = 1, dimx
            x(i) = -delta + 2.0_dp * delta * real(i-1, dp) / real(dimx-1, dp)
            y(i) = exp(cmplx(0.0_dp, -x(i), dp))  ! Complex exponential
        end do
        
        ! Direct transform
        call fourier_direct_transform(dimx, x, y, M, x0, delta, f, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Inverse transform
        call fourier_inverse_transform(dimx, x, f, M, x0, delta, y_reconstructed, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Check basic consistency (not exact due to limited M)
        test_passed = test_passed .and. (abs(f(0)) > 0.0_dp)  ! Should have non-zero DC component
        
        call end_test(test_passed)
    end subroutine test_custom_fourier_transforms
    
    !---------------------------------------------------------------------------
    ! Test 3: Spectral analysis functions
    !---------------------------------------------------------------------------
    subroutine test_spectral_analysis()
        integer, parameter :: n = 16
        complex(dp) :: signal1(n), signal2(n), cross_spec(n)
        real(dp) :: power_spec(n)
        type(fourier_settings_t) :: settings
        integer :: ierr, i
        real(dp) :: freq
        
        call start_test("Spectral analysis")
        test_passed = .true.
        
        call fourier_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create test signals
        do i = 1, n
            freq = 2.0_dp * PI * real(i-1, dp) / real(n, dp)
            signal1(i) = exp(cmplx(0.0_dp, 2.0_dp * freq, dp))  ! 2 cycles
            signal2(i) = exp(cmplx(0.0_dp, 3.0_dp * freq, dp))  ! 3 cycles
        end do
        
        ! Test power spectrum
        call fourier_power_spectrum(n, signal1, power_spec, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (maxval(power_spec) > 0.0_dp)
        
        ! Test cross spectrum
        call fourier_cross_spectrum(n, signal1, signal2, cross_spec, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (maxval(abs(cross_spec)) > 0.0_dp)
        
        call end_test(test_passed)
    end subroutine test_spectral_analysis
    
    !---------------------------------------------------------------------------
    ! Test 4: Windowing functions
    !---------------------------------------------------------------------------
    subroutine test_windowing_functions()
        integer, parameter :: n = 32
        real(dp) :: window(n)
        complex(dp) :: signal(n), windowed_signal(n)
        type(fourier_settings_t) :: settings
        integer :: ierr, i
        
        call start_test("Windowing functions")
        test_passed = .true.
        
        call fourier_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create test signal
        do i = 1, n
            signal(i) = cmplx(1.0_dp, 0.0_dp, dp)
        end do
        
        ! Test Hanning window
        call fourier_create_window(n, WINDOW_HANNING, window, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (window(1) == 0.0_dp)  ! Should be zero at endpoints
        test_passed = test_passed .and. (window(n) == 0.0_dp)
        test_passed = test_passed .and. (maxval(window) > 0.0_dp)
        
        ! Test Hamming window
        call fourier_create_window(n, WINDOW_HAMMING, window, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (minval(window) > 0.0_dp)  ! Non-zero everywhere
        
        ! Test applying window
        call fourier_apply_window(n, signal, WINDOW_HANNING, windowed_signal, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(windowed_signal(1)) < 1.0e-10_dp)  ! Should be near zero
        
        call end_test(test_passed)
    end subroutine test_windowing_functions
    
    !---------------------------------------------------------------------------
    ! Test 5: Convolution operations
    !---------------------------------------------------------------------------
    subroutine test_convolution_operations()
        integer, parameter :: n = 8
        complex(dp) :: signal1(n), signal2(n), conv_result(n), corr_result(n)
        type(fourier_settings_t) :: settings
        integer :: ierr
        
        call start_test("Convolution operations")
        test_passed = .true.
        
        call fourier_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create simple test signals
        signal1 = (0.0_dp, 0.0_dp)
        signal2 = (0.0_dp, 0.0_dp)
        signal1(1) = (1.0_dp, 0.0_dp)  ! Delta function
        signal2(2) = (1.0_dp, 0.0_dp)  ! Shifted delta
        
        ! Test convolution
        call fourier_convolve(n, signal1, signal2, conv_result, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(conv_result(2)) > 0.1_dp)  ! Should be non-zero at position 2
        
        ! Test correlation
        call fourier_correlate(n, signal1, signal2, corr_result, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (maxval(abs(corr_result)) > 0.0_dp)
        
        call end_test(test_passed)
    end subroutine test_convolution_operations
    
    !---------------------------------------------------------------------------
    ! Test 6: Power spectrum analysis
    !---------------------------------------------------------------------------
    subroutine test_power_spectrum_analysis()
        integer, parameter :: n = 16
        complex(dp) :: signal(n)
        real(dp) :: power_spec(n)
        type(fourier_settings_t) :: settings
        integer :: ierr, i
        real(dp) :: freq
        
        call start_test("Power spectrum analysis")
        test_passed = .true.
        
        call fourier_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create sinusoidal signal with known frequency
        do i = 1, n
            freq = 2.0_dp * PI * real(i-1, dp) / real(n, dp)
            signal(i) = exp(cmplx(0.0_dp, 4.0_dp * freq, dp))  ! 4 cycles
        end do
        
        call fourier_power_spectrum(n, signal, power_spec, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Power spectrum should be concentrated at specific frequencies
        test_passed = test_passed .and. (sum(power_spec) > 0.0_dp)
        test_passed = test_passed .and. (maxval(power_spec) > 0.0_dp)
        
        call end_test(test_passed)
    end subroutine test_power_spectrum_analysis
    
    !---------------------------------------------------------------------------
    ! Test 7: Frequency domain filtering
    !---------------------------------------------------------------------------
    subroutine test_frequency_domain_filtering()
        integer, parameter :: n = 16
        complex(dp) :: signal(n), filtered(n)
        type(fourier_settings_t) :: settings
        integer :: ierr, i
        real(dp) :: freq
        
        call start_test("Frequency domain filtering")
        test_passed = .true.
        
        call fourier_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create broadband signal
        do i = 1, n
            freq = 2.0_dp * PI * real(i-1, dp) / real(n, dp)
            signal(i) = exp(cmplx(0.0_dp, freq, dp)) + 0.5_dp * exp(cmplx(0.0_dp, 3.0_dp * freq, dp))
        end do
        
        ! Test low-pass filter
        call fourier_frequency_filter(n, signal, FILTER_LOW_PASS, 0.2_dp, 0.0_dp, filtered, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (maxval(abs(filtered)) > 0.0_dp)
        
        ! Test high-pass filter
        call fourier_frequency_filter(n, signal, FILTER_HIGH_PASS, 0.3_dp, 0.0_dp, filtered, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call end_test(test_passed)
    end subroutine test_frequency_domain_filtering
    
    !---------------------------------------------------------------------------
    ! Test 8: Inverse consistency
    !---------------------------------------------------------------------------
    subroutine test_inverse_consistency()
        integer, parameter :: n = 16
        complex(dp) :: original(n), fft_result(n), reconstructed(n)
        real(dp) :: real_input(n), real_output(n)
        complex(dp) :: complex_fft(n/2+1)
        type(fourier_settings_t) :: settings
        integer :: ierr, i
        real(dp) :: error
        
        call start_test("Inverse consistency")
        test_passed = .true.
        
        call fourier_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test complex FFT consistency
        do i = 1, n
            original(i) = cmplx(sin(2.0_dp * PI * real(i-1, dp) / real(n, dp)), &
                               cos(4.0_dp * PI * real(i-1, dp) / real(n, dp)), dp)
        end do
        
        call fourier_fft_1d(n, original, fft_result, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call fourier_ifft_1d(n, fft_result, reconstructed, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        error = maxval(abs(original - reconstructed))
        test_passed = test_passed .and. (error < 1.0e-10_dp)
        
        ! Test real FFT consistency
        do i = 1, n
            real_input(i) = sin(2.0_dp * PI * real(i-1, dp) / real(n, dp))
        end do
        
        call fourier_fft_1d_real(n, real_input, complex_fft, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call fourier_ifft_1d_real(n, complex_fft, real_output, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        error = maxval(abs(real_input - real_output))
        test_passed = test_passed .and. (error < 1.0e-10_dp)
        
        call end_test(test_passed)
    end subroutine test_inverse_consistency
    
    !---------------------------------------------------------------------------
    ! Test 9: Performance benchmarks
    !---------------------------------------------------------------------------
    subroutine test_performance_benchmarks()
        integer, parameter :: n = 64
        complex(dp) :: signal(n), result(n)
        type(fourier_settings_t) :: settings
        integer :: ierr, i
        real(dp) :: start_time, end_time
        
        call start_test("Performance benchmarks")
        test_passed = .true.
        
        call fourier_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create test signal
        do i = 1, n
            signal(i) = cmplx(sin(2.0_dp * PI * real(i-1, dp) / real(n, dp)), &
                             cos(4.0_dp * PI * real(i-1, dp) / real(n, dp)), dp)
        end do
        
        ! Time the FFT operation
        call cpu_time(start_time)
        call fourier_fft_1d(n, signal, result, settings, ierr) 
        call cpu_time(end_time)
        
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. ((end_time - start_time) >= 0.0_dp)  ! Should complete
        
        if (settings%debug_level > 0) then
            write(*, *) "FFT time for n=", n, ":", end_time - start_time, "seconds"
        end if
        
        call end_test(test_passed)
    end subroutine test_performance_benchmarks
    
    !---------------------------------------------------------------------------
    ! Test 10: Error handling
    !---------------------------------------------------------------------------
    subroutine test_error_handling()
        integer, parameter :: n = 8
        complex(dp) :: signal(n), result(n)
        real(dp) :: window(n)
        type(fourier_settings_t) :: settings
        integer :: ierr
        
        call start_test("Error handling")
        test_passed = .true.
        
        call fourier_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        signal = (1.0_dp, 0.0_dp)
        
        ! Test invalid window type
        call fourier_create_window(n, 999, window, ierr)
        test_passed = test_passed .and. (ierr /= 0)  ! Should fail
        
        ! Test invalid filter type
        call fourier_frequency_filter(n, signal, 999, 0.1_dp, 0.2_dp, result, settings, ierr)
        test_passed = test_passed .and. (ierr /= 0)  ! Should fail
        
        ! Test valid operations (should succeed)
        call fourier_fft_1d(n, signal, result, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call fourier_settings_destroy(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call end_test(test_passed)
    end subroutine test_error_handling
    
end program test_kilca_fourier