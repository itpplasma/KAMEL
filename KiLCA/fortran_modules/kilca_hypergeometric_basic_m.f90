!> Basic hypergeometric functions module for KiLCA
!! Provides pure Fortran implementations of hypergeometric functions
!! for specialized cases needed in plasma physics applications
module kilca_hypergeometric_basic_m
    use kilca_types_m
    use iso_fortran_env, only: real64
    implicit none
    
    private
    
    ! Public types
    public :: hyperg_1f1_settings_t
    
    ! Public procedures
    public :: hyperg_1f1_kummer_fortran
    public :: hyperg_1f1_series_fortran
    public :: hyperg_1f1_asymptotic_fortran
    public :: hyperg_compute_coefficients
    public :: hyperg_check_convergence
    
    !> Settings for 1F1 calculation algorithms
    type :: hyperg_1f1_settings_t
        integer :: algorithm = 1                     !< Algorithm selection (1=series, 2=asymptotic)
        real(dp) :: tolerance = 1.0e-12_dp          !< Convergence tolerance
        integer :: max_iterations = 1000            !< Maximum iterations
        logical :: use_acceleration = .false.       !< Use series acceleration
        integer :: debug_level = 0                  !< Debug output level
    end type hyperg_1f1_settings_t
    
    ! Algorithm selection constants
    integer, parameter, public :: HYPERG_ALGORITHM_SERIES = 1
    integer, parameter, public :: HYPERG_ALGORITHM_ASYMPTOTIC = 2
    
contains

    !---------------------------------------------------------------------------
    ! Kummer series implementation in pure Fortran
    !---------------------------------------------------------------------------
    
    !> Confluent hypergeometric function 1F1(a,b,z) using Kummer series
    !! 1F1(a,b,z) = sum_{n=0}^∞ (a)_n / (b)_n * z^n / n!
    !! where (x)_n is the Pochhammer symbol (rising factorial)
    subroutine hyperg_1f1_kummer_fortran(a, b, z, result, settings, ierr)
        real(dp), intent(in) :: a, b, z
        complex(dp), intent(out) :: result
        type(hyperg_1f1_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        complex(dp) :: term, sum_val
        real(dp) :: a_n, b_n  ! Pochhammer symbols
        real(dp) :: factorial_n
        integer :: n
        real(dp) :: conv_test
        
        ierr = 0
        
        ! Check for problematic parameters
        if (abs(b) < epsilon(1.0_dp)) then
            ierr = -1  ! b = 0 is singular
            result = (0.0_dp, 0.0_dp)
            return
        end if
        
        ! Initialize series
        sum_val = (1.0_dp, 0.0_dp)
        term = (1.0_dp, 0.0_dp)
        a_n = 1.0_dp
        b_n = 1.0_dp
        factorial_n = 1.0_dp
        
        do n = 1, settings%max_iterations
            ! Update Pochhammer symbols
            a_n = a_n * (a + real(n-1, dp))
            b_n = b_n * (b + real(n-1, dp))
            factorial_n = factorial_n * real(n, dp)
            
            ! Compute term: (a)_n / (b)_n * z^n / n!
            term = cmplx(a_n / b_n / factorial_n, 0.0_dp, dp) * cmplx(z, 0.0_dp, dp)**n
            
            sum_val = sum_val + term
            
            ! Check convergence
            conv_test = abs(term) / max(abs(sum_val), epsilon(1.0_dp))
            if (conv_test < settings%tolerance) then
                result = sum_val
                return
            end if
            
            ! Prevent overflow
            if (abs(term) > huge(1.0_dp) * 0.1_dp) then
                ierr = -2  ! Overflow
                result = (0.0_dp, 0.0_dp)
                return
            end if
        end do
        
        ! Did not converge
        ierr = -3
        result = sum_val
        
    end subroutine hyperg_1f1_kummer_fortran
    
    !---------------------------------------------------------------------------
    ! Series implementation for complex arguments
    !---------------------------------------------------------------------------
    
    !> Series implementation for 1F1(1,b,z) with complex b and z
    subroutine hyperg_1f1_series_fortran(b, z, result, settings, ierr)
        complex(dp), intent(in) :: b, z
        complex(dp), intent(out) :: result
        type(hyperg_1f1_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        complex(dp) :: term, sum_val, z_power
        complex(dp) :: b_n  ! Pochhammer symbol for complex b
        real(dp) :: factorial_n
        integer :: n
        real(dp) :: conv_test
        
        ierr = 0
        
        ! Check for problematic parameters
        if (abs(b) < epsilon(1.0_dp)) then
            ierr = -1
            result = (0.0_dp, 0.0_dp)
            return
        end if
        
        ! Initialize series for 1F1(1, b, z)
        sum_val = (1.0_dp, 0.0_dp)
        term = (1.0_dp, 0.0_dp)
        b_n = (1.0_dp, 0.0_dp)
        z_power = (1.0_dp, 0.0_dp)
        factorial_n = 1.0_dp
        
        do n = 1, settings%max_iterations
            ! Update factors
            b_n = b_n * (b + cmplx(real(n-1, dp), 0.0_dp, dp))
            z_power = z_power * z
            factorial_n = factorial_n * real(n, dp)
            
            ! Compute term: (1)_n / (b)_n * z^n / n! = n / (b + n-1) * z^n / n!
            term = cmplx(real(n, dp), 0.0_dp, dp) / b_n * z_power / cmplx(factorial_n, 0.0_dp, dp)
            
            sum_val = sum_val + term
            
            ! Check convergence
            conv_test = abs(term) / max(abs(sum_val), epsilon(1.0_dp))
            if (conv_test < settings%tolerance) then
                result = sum_val
                return
            end if
            
            ! Prevent overflow
            if (abs(term) > huge(1.0_dp) * 0.1_dp) then
                ierr = -2
                result = (0.0_dp, 0.0_dp)
                return
            end if
        end do
        
        ierr = -3
        result = sum_val
        
    end subroutine hyperg_1f1_series_fortran
    
    !---------------------------------------------------------------------------
    ! Asymptotic expansion for large |z|
    !---------------------------------------------------------------------------
    
    !> Asymptotic expansion of 1F1(a,b,z) for large |z|
    subroutine hyperg_1f1_asymptotic_fortran(a, b, z, result, settings, ierr)
        real(dp), intent(in) :: a, b, z
        complex(dp), intent(out) :: result
        type(hyperg_1f1_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        complex(dp) :: asymptotic_term
        real(dp) :: gamma_ratio
        
        ierr = 0
        
        ! For large z, 1F1(a,b,z) ~ Γ(b)/Γ(a) * exp(z) * z^(a-b)
        ! This is a simplified asymptotic approximation
        
        if (abs(z) < 10.0_dp) then
            ierr = -1  ! Asymptotic expansion not valid for small |z|
            result = (0.0_dp, 0.0_dp)
            return
        end if
        
        ! Approximate gamma function ratio using Stirling's approximation
        ! Γ(b)/Γ(a) ≈ (b/a)^a * exp(a-b) for large a,b
        if (abs(a) > 1.0_dp .and. abs(b) > 1.0_dp) then
            gamma_ratio = (b/a)**a * exp(a - b)
        else
            gamma_ratio = 1.0_dp  ! Simplified
        end if
        
        asymptotic_term = cmplx(gamma_ratio, 0.0_dp, dp) * exp(cmplx(z, 0.0_dp, dp)) * &
                         cmplx(z, 0.0_dp, dp)**(a - b)
        
        result = asymptotic_term
        
    end subroutine hyperg_1f1_asymptotic_fortran
    
    !---------------------------------------------------------------------------
    ! Helper functions
    !---------------------------------------------------------------------------
    
    !> Compute Pochhammer symbol (rising factorial) (x)_n = x(x+1)...(x+n-1)
    function hyperg_compute_coefficients(x, n) result(pochhammer)
        real(dp), intent(in) :: x
        integer, intent(in) :: n
        real(dp) :: pochhammer
        integer :: i
        
        pochhammer = 1.0_dp
        do i = 0, n-1
            pochhammer = pochhammer * (x + real(i, dp))
        end do
        
    end function hyperg_compute_coefficients
    
    !> Check convergence of series
    function hyperg_check_convergence(current_term, sum_val, tolerance) result(converged)
        complex(dp), intent(in) :: current_term, sum_val
        real(dp), intent(in) :: tolerance
        logical :: converged
        
        real(dp) :: relative_error
        
        relative_error = abs(current_term) / max(abs(sum_val), epsilon(1.0_dp))
        converged = (relative_error < tolerance)
        
    end function hyperg_check_convergence
    
end module kilca_hypergeometric_basic_m