!> Fourier Transform System for KiLCA
!! Provides FFT operations, custom spectral analysis, and signal processing functions
!! Supports both FFTW interface and custom implementations for plasma physics applications
module kilca_fourier_m
    use kilca_types_m
    use iso_fortran_env, only: real64
    use iso_c_binding, only: c_double, c_int, c_ptr
    implicit none
    
    private
    
    ! Public types
    public :: fourier_settings_t
    
    ! Public procedures
    public :: fourier_fft_1d, fourier_ifft_1d
    public :: fourier_fft_1d_real, fourier_ifft_1d_real
    public :: fourier_direct_transform, fourier_inverse_transform
    public :: fourier_power_spectrum, fourier_cross_spectrum
    public :: fourier_convolve, fourier_correlate
    public :: fourier_apply_window, fourier_create_window
    public :: fourier_frequency_filter
    public :: fourier_settings_create, fourier_settings_destroy
    
    !> Fourier transform settings and configuration
    type :: fourier_settings_t
        integer :: method = 1                        !< Method: 1=FFTW, 2=custom
        integer :: window_type = 1                   !< Window function type
        real(dp) :: tolerance = 1.0e-12_dp          !< Numerical tolerance
        integer :: max_iterations = 1000            !< Max iterations for custom methods
        logical :: use_padding = .true.             !< Zero-padding for efficiency
        integer :: debug_level = 0                  !< Debug output level
        ! Custom Fourier transform parameters (from C++ implementation)
        integer :: spline_order = 5                 !< Spline interpolation order
        real(dp) :: integration_tolerance = 1.0e-8_dp !< GSL integration tolerance
        integer :: integration_limit = 1000         !< GSL integration workspace limit
    end type fourier_settings_t
    
    ! Method constants
    integer, parameter, public :: FOURIER_METHOD_FFTW = 1
    integer, parameter, public :: FOURIER_METHOD_CUSTOM = 2
    
    ! Window function constants
    integer, parameter, public :: WINDOW_RECTANGULAR = 1
    integer, parameter, public :: WINDOW_HANNING = 2
    integer, parameter, public :: WINDOW_HAMMING = 3
    integer, parameter, public :: WINDOW_BLACKMAN = 4
    integer, parameter, public :: WINDOW_KAISER = 5
    
    ! Filter constants
    integer, parameter, public :: FILTER_LOW_PASS = 1
    integer, parameter, public :: FILTER_HIGH_PASS = 2
    integer, parameter, public :: FILTER_BAND_PASS = 3
    integer, parameter, public :: FILTER_BAND_STOP = 4
    
contains

    !---------------------------------------------------------------------------
    ! Settings management
    !---------------------------------------------------------------------------
    
    !> Create and initialize Fourier settings
    subroutine fourier_settings_create(settings, ierr)
        type(fourier_settings_t), intent(out) :: settings
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Set default values (already set in type definition)
        settings%method = FOURIER_METHOD_FFTW
        settings%window_type = WINDOW_RECTANGULAR
        settings%tolerance = 1.0e-12_dp
        settings%max_iterations = 1000
        settings%use_padding = .true.
        settings%debug_level = 0
        settings%spline_order = 5
        settings%integration_tolerance = 1.0e-8_dp
        settings%integration_limit = 1000
        
    end subroutine fourier_settings_create
    
    !> Destroy Fourier settings (cleanup if needed)
    subroutine fourier_settings_destroy(settings, ierr)
        type(fourier_settings_t), intent(inout) :: settings
        integer, intent(out) :: ierr
        
        ierr = 0
        ! No dynamic memory to clean up in current implementation
        
    end subroutine fourier_settings_destroy
    
    !---------------------------------------------------------------------------
    ! Basic FFT operations (FFTW interface stubs for now)
    !---------------------------------------------------------------------------
    
    !> 1D complex FFT using FFTW
    subroutine fourier_fft_1d(n, input, output, settings, ierr)
        integer, intent(in) :: n
        complex(dp), intent(in) :: input(n)
        complex(dp), intent(out) :: output(n)
        type(fourier_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        integer :: i
        
        ierr = 0
        
        ! Simple DFT implementation (placeholder for FFTW)
        ! For production code, this would use FFTW3 library
        call compute_dft_1d(n, input, output, .false., ierr)
        
        if (settings%debug_level > 0) then
            write(*, *) "FFT completed for n =", n
        end if
        
    end subroutine fourier_fft_1d
    
    !> 1D complex inverse FFT using FFTW
    subroutine fourier_ifft_1d(n, input, output, settings, ierr)
        integer, intent(in) :: n
        complex(dp), intent(in) :: input(n)
        complex(dp), intent(out) :: output(n)
        type(fourier_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Simple inverse DFT implementation
        call compute_dft_1d(n, input, output, .true., ierr)
        
        if (settings%debug_level > 0) then
            write(*, *) "IFFT completed for n =", n
        end if
        
    end subroutine fourier_ifft_1d
    
    !> 1D real FFT (returns complex output)
    subroutine fourier_fft_1d_real(n, input, output, settings, ierr)
        integer, intent(in) :: n
        real(dp), intent(in) :: input(n)
        complex(dp), intent(out) :: output(n/2+1)
        type(fourier_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        complex(dp), allocatable :: temp_input(:)
        complex(dp), allocatable :: temp_output(:)
        integer :: i
        
        ierr = 0
        
        ! Convert real input to complex
        allocate(temp_input(n), temp_output(n))
        
        do i = 1, n
            temp_input(i) = cmplx(input(i), 0.0_dp, dp)
        end do
        
        call compute_dft_1d(n, temp_input, temp_output, .false., ierr)
        
        ! Extract only positive frequencies (exploit conjugate symmetry)
        do i = 1, n/2+1
            output(i) = temp_output(i)
        end do
        
        deallocate(temp_input, temp_output)
        
    end subroutine fourier_fft_1d_real
    
    !> 1D real inverse FFT (from complex input)
    subroutine fourier_ifft_1d_real(n, input, output, settings, ierr)
        integer, intent(in) :: n
        complex(dp), intent(in) :: input(n/2+1)
        real(dp), intent(out) :: output(n)
        type(fourier_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        complex(dp), allocatable :: temp_input(:), temp_output(:)
        integer :: i
        
        ierr = 0
        
        allocate(temp_input(n), temp_output(n))
        
        ! Reconstruct full complex array using conjugate symmetry
        temp_input(1) = input(1)  ! DC component
        do i = 2, n/2+1
            temp_input(i) = input(i)
            if (i <= n/2) then
                temp_input(n-i+2) = conjg(input(i))
            end if
        end do
        
        call compute_dft_1d(n, temp_input, temp_output, .true., ierr)
        
        ! Extract real part
        do i = 1, n
            output(i) = real(temp_output(i), dp)
        end do
        
        deallocate(temp_input, temp_output)
        
    end subroutine fourier_ifft_1d_real
    
    !---------------------------------------------------------------------------
    ! Custom Fourier transforms (translation of C++ functions)
    !---------------------------------------------------------------------------
    
    !> Direct Fourier transform with spline interpolation
    !! Translation of calc_direct_fourier_ function
    subroutine fourier_direct_transform(dimx, x, y, M, x0, delta, f, settings, ierr)
        integer, intent(in) :: dimx, M
        real(dp), intent(in) :: x(dimx)
        complex(dp), intent(in) :: y(dimx)
        real(dp), intent(in) :: x0, delta
        complex(dp), intent(out) :: f(-M:M)
        type(fourier_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        integer :: k
        real(dp) :: k0, x_min, x_max
        complex(dp) :: integral_result
        
        ierr = 0
        
        k0 = PI / delta
        x_min = x0 - delta
        x_max = x0 + delta
        
        ! Simplified implementation without spline interpolation
        ! For production code, this would use spline interpolation + GSL integration
        do k = -M, M
            call compute_fourier_integral(dimx, x, y, k, k0, x0, x_min, x_max, &
                                        integral_result, settings, ierr)
            if (ierr /= 0) return
            
            f(k) = integral_result / (2.0_dp * delta)
        end do
        
        if (settings%debug_level > 0) then
            write(*, *) "Direct Fourier transform completed, M =", M
        end if
        
    end subroutine fourier_direct_transform
    
    !> Inverse Fourier transform
    !! Translation of calc_inverse_fourier_ function  
    subroutine fourier_inverse_transform(dimx, x, f, M, x0, delta, y, settings, ierr)
        integer, intent(in) :: dimx, M
        real(dp), intent(in) :: x(dimx)
        complex(dp), intent(in) :: f(-M:M)
        real(dp), intent(in) :: x0, delta
        complex(dp), intent(out) :: y(dimx)
        type(fourier_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        integer :: i, k
        real(dp) :: k0
        complex(dp) :: arg, tmp
        
        ierr = 0
        
        k0 = PI / delta
        
        do i = 1, dimx
            tmp = (0.0_dp, 0.0_dp)
            do k = -M, M
                arg = cmplx(0.0_dp, (x(i) - x0) * k0 * real(k, dp), dp)
                tmp = tmp + f(k) * exp(arg)
            end do
            y(i) = tmp
        end do
        
        if (settings%debug_level > 0) then
            write(*, *) "Inverse Fourier transform completed"
        end if
        
    end subroutine fourier_inverse_transform
    
    !---------------------------------------------------------------------------
    ! Spectral analysis functions
    !---------------------------------------------------------------------------
    
    !> Compute power spectrum
    subroutine fourier_power_spectrum(n, input, output, settings, ierr)
        integer, intent(in) :: n  
        complex(dp), intent(in) :: input(n)
        real(dp), intent(out) :: output(n)
        type(fourier_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        complex(dp), allocatable :: fft_result(:)
        integer :: i
        
        ierr = 0
        
        allocate(fft_result(n))
        
        ! Compute FFT
        call fourier_fft_1d(n, input, fft_result, settings, ierr)
        if (ierr /= 0) then
            deallocate(fft_result)
            return
        end if
        
        ! Compute power spectrum: |FFT|^2
        do i = 1, n
            output(i) = abs(fft_result(i))**2
        end do
        
        deallocate(fft_result)
        
    end subroutine fourier_power_spectrum
    
    !> Compute cross spectrum between two signals
    subroutine fourier_cross_spectrum(n, input1, input2, output, settings, ierr)
        integer, intent(in) :: n
        complex(dp), intent(in) :: input1(n), input2(n)
        complex(dp), intent(out) :: output(n)
        type(fourier_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        complex(dp), allocatable :: fft1(:), fft2(:)
        integer :: i
        
        ierr = 0
        
        allocate(fft1(n), fft2(n))
        
        ! Compute FFTs of both signals
        call fourier_fft_1d(n, input1, fft1, settings, ierr)
        if (ierr /= 0) then
            deallocate(fft1, fft2)
            return
        end if
        
        call fourier_fft_1d(n, input2, fft2, settings, ierr)
        if (ierr /= 0) then
            deallocate(fft1, fft2)
            return
        end if
        
        ! Cross spectrum: FFT1 * conj(FFT2)
        do i = 1, n
            output(i) = fft1(i) * conjg(fft2(i))
        end do
        
        deallocate(fft1, fft2)
        
    end subroutine fourier_cross_spectrum
    
    !---------------------------------------------------------------------------
    ! Convolution and correlation
    !---------------------------------------------------------------------------
    
    !> FFT-based convolution
    subroutine fourier_convolve(n, input1, input2, output, settings, ierr)
        integer, intent(in) :: n
        complex(dp), intent(in) :: input1(n), input2(n)
        complex(dp), intent(out) :: output(n)
        type(fourier_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        complex(dp), allocatable :: fft1(:), fft2(:), prod(:)
        integer :: i
        
        ierr = 0
        
        allocate(fft1(n), fft2(n), prod(n))
        
        ! Compute FFTs
        call fourier_fft_1d(n, input1, fft1, settings, ierr)
        if (ierr /= 0) goto 100
        
        call fourier_fft_1d(n, input2, fft2, settings, ierr)
        if (ierr /= 0) goto 100
        
        ! Multiply in frequency domain
        do i = 1, n
            prod(i) = fft1(i) * fft2(i)
        end do
        
        ! Inverse FFT
        call fourier_ifft_1d(n, prod, output, settings, ierr)
        
100     continue
        deallocate(fft1, fft2, prod)
        
    end subroutine fourier_convolve
    
    !> FFT-based correlation
    subroutine fourier_correlate(n, input1, input2, output, settings, ierr)
        integer, intent(in) :: n
        complex(dp), intent(in) :: input1(n), input2(n)
        complex(dp), intent(out) :: output(n)
        type(fourier_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        complex(dp), allocatable :: fft1(:), fft2(:), prod(:)
        integer :: i
        
        ierr = 0
        
        allocate(fft1(n), fft2(n), prod(n))
        
        ! Compute FFTs
        call fourier_fft_1d(n, input1, fft1, settings, ierr)
        if (ierr /= 0) goto 200
        
        call fourier_fft_1d(n, input2, fft2, settings, ierr)
        if (ierr /= 0) goto 200
        
        ! Correlation: FFT1 * conj(FFT2)
        do i = 1, n
            prod(i) = fft1(i) * conjg(fft2(i))
        end do
        
        ! Inverse FFT
        call fourier_ifft_1d(n, prod, output, settings, ierr)
        
200     continue
        deallocate(fft1, fft2, prod)
        
    end subroutine fourier_correlate
    
    !---------------------------------------------------------------------------
    ! Window functions
    !---------------------------------------------------------------------------
    
    !> Create window function
    subroutine fourier_create_window(n, window_type, window, ierr)
        integer, intent(in) :: n, window_type
        real(dp), intent(out) :: window(n)
        integer, intent(out) :: ierr
        
        integer :: i
        real(dp) :: arg, alpha
        
        ierr = 0
        
        select case(window_type)
        case(WINDOW_RECTANGULAR)
            window = 1.0_dp
            
        case(WINDOW_HANNING)
            do i = 1, n
                arg = 2.0_dp * PI * real(i-1, dp) / real(n-1, dp)
                window(i) = 0.5_dp * (1.0_dp - cos(arg))
            end do
            
        case(WINDOW_HAMMING)
            do i = 1, n
                arg = 2.0_dp * PI * real(i-1, dp) / real(n-1, dp)
                window(i) = 0.54_dp - 0.46_dp * cos(arg)
            end do
            
        case(WINDOW_BLACKMAN)
            do i = 1, n
                arg = 2.0_dp * PI * real(i-1, dp) / real(n-1, dp)
                window(i) = 0.42_dp - 0.5_dp * cos(arg) + 0.08_dp * cos(2.0_dp * arg)
            end do
            
        case(WINDOW_KAISER)
            ! Simplified Kaiser window (beta = 5)
            alpha = 5.0_dp
            do i = 1, n
                arg = 2.0_dp * real(i-1, dp) / real(n-1, dp) - 1.0_dp
                window(i) = bessel_i0(alpha * sqrt(1.0_dp - arg**2)) / bessel_i0(alpha)
            end do
            
        case default
            ierr = -1
            window = 0.0_dp
        end select
        
    end subroutine fourier_create_window
    
    !> Apply window function to signal
    subroutine fourier_apply_window(n, input, window_type, output, settings, ierr)
        integer, intent(in) :: n, window_type
        complex(dp), intent(in) :: input(n)
        complex(dp), intent(out) :: output(n)
        type(fourier_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        real(dp), allocatable :: window(:)
        integer :: i
        
        ierr = 0
        
        allocate(window(n))
        
        call fourier_create_window(n, window_type, window, ierr)
        if (ierr /= 0) then
            deallocate(window)
            return
        end if
        
        do i = 1, n
            output(i) = input(i) * cmplx(window(i), 0.0_dp, dp)
        end do
        
        deallocate(window)
        
    end subroutine fourier_apply_window
    
    !---------------------------------------------------------------------------
    ! Frequency domain filtering
    !---------------------------------------------------------------------------
    
    !> Apply frequency domain filter
    subroutine fourier_frequency_filter(n, input, filter_type, cutoff_low, &
                                       cutoff_high, output, settings, ierr)
        integer, intent(in) :: n, filter_type
        complex(dp), intent(in) :: input(n)
        real(dp), intent(in) :: cutoff_low, cutoff_high
        complex(dp), intent(out) :: output(n)
        type(fourier_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        complex(dp), allocatable :: fft_data(:)
        real(dp) :: freq, nyquist
        integer :: i
        
        ierr = 0
        
        allocate(fft_data(n))
        
        ! Forward FFT
        call fourier_fft_1d(n, input, fft_data, settings, ierr)
        if (ierr /= 0) then
            deallocate(fft_data)
            return
        end if
        
        nyquist = 0.5_dp
        
        ! Apply filter in frequency domain
        do i = 1, n
            freq = real(i-1, dp) / real(n, dp)
            if (freq > 0.5_dp) freq = freq - 1.0_dp  ! Handle negative frequencies
            freq = abs(freq)
            
            select case(filter_type)
            case(FILTER_LOW_PASS)
                if (freq > cutoff_low) fft_data(i) = (0.0_dp, 0.0_dp)
                
            case(FILTER_HIGH_PASS)
                if (freq < cutoff_low) fft_data(i) = (0.0_dp, 0.0_dp)
                
            case(FILTER_BAND_PASS)
                if (freq < cutoff_low .or. freq > cutoff_high) &
                    fft_data(i) = (0.0_dp, 0.0_dp)
                    
            case(FILTER_BAND_STOP)
                if (freq >= cutoff_low .and. freq <= cutoff_high) &
                    fft_data(i) = (0.0_dp, 0.0_dp)
                    
            case default
                ierr = -1
                deallocate(fft_data)
                return
            end select
        end do
        
        ! Inverse FFT
        call fourier_ifft_1d(n, fft_data, output, settings, ierr)
        
        deallocate(fft_data)
        
    end subroutine fourier_frequency_filter
    
    !---------------------------------------------------------------------------
    ! Helper functions
    !---------------------------------------------------------------------------
    
    !> Compute 1D DFT/IDFT using direct method
    subroutine compute_dft_1d(n, input, output, inverse, ierr)
        integer, intent(in) :: n
        complex(dp), intent(in) :: input(n)
        complex(dp), intent(out) :: output(n)
        logical, intent(in) :: inverse
        integer, intent(out) :: ierr
        
        integer :: k, j
        complex(dp) :: arg, sum_val
        real(dp) :: sign_factor, norm_factor
        
        ierr = 0
        
        sign_factor = merge(1.0_dp, -1.0_dp, inverse)
        norm_factor = merge(1.0_dp / real(n, dp), 1.0_dp, inverse)
        
        do k = 1, n
            sum_val = (0.0_dp, 0.0_dp)
            do j = 1, n
                arg = cmplx(0.0_dp, sign_factor * 2.0_dp * PI * real((k-1)*(j-1), dp) / real(n, dp), dp)
                sum_val = sum_val + input(j) * exp(arg)
            end do
            output(k) = sum_val * norm_factor
        end do
        
    end subroutine compute_dft_1d
    
    !> Compute Fourier integral (simplified version of C++ spline+integration)
    subroutine compute_fourier_integral(dimx, x, y, k, k0, x0, x_min, x_max, &
                                       result, settings, ierr)
        integer, intent(in) :: dimx, k
        real(dp), intent(in) :: x(dimx), k0, x0, x_min, x_max
        complex(dp), intent(in) :: y(dimx)
        complex(dp), intent(out) :: result
        type(fourier_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        integer :: i
        real(dp) :: dx
        complex(dp) :: arg, integrand
        
        ierr = 0
        
        ! Simple trapezoidal integration (placeholder for spline+GSL integration)
        result = (0.0_dp, 0.0_dp)
        
        if (dimx < 2) then
            ierr = -1
            return
        end if
        
        do i = 1, dimx-1
            if (x(i) >= x_min .and. x(i+1) <= x_max) then
                dx = x(i+1) - x(i)
                arg = cmplx(0.0_dp, -(x(i) - x0) * k0 * real(k, dp), dp)
                integrand = exp(arg) * y(i)
                result = result + integrand * dx
            end if
        end do
        
    end subroutine compute_fourier_integral
    
    !> Modified Bessel function I0 (simplified approximation)
    function bessel_i0(x) result(i0)
        real(dp), intent(in) :: x
        real(dp) :: i0
        real(dp) :: t, sum_val
        integer :: n
        
        if (abs(x) < 3.75_dp) then
            t = (x / 3.75_dp)**2
            i0 = 1.0_dp + 3.5156229_dp*t + 3.0899424_dp*t**2 + 1.2067492_dp*t**3 + &
                 0.2659732_dp*t**4 + 0.0360768_dp*t**5 + 0.0045813_dp*t**6
        else
            t = 3.75_dp / abs(x)
            i0 = exp(abs(x)) / sqrt(abs(x)) * (0.39894228_dp + 0.01328592_dp*t + &
                 0.00225319_dp*t**2 - 0.00157565_dp*t**3 + 0.00916281_dp*t**4 - &
                 0.02057706_dp*t**5 + 0.02635537_dp*t**6 - 0.01647633_dp*t**7 + &
                 0.00392377_dp*t**8)
        end if
        
    end function bessel_i0
    
end module kilca_fourier_m