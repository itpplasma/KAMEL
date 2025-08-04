module kilca_plasma_physics_m
    use iso_fortran_env, only: real64
    implicit none
    
    private
    
    ! Public procedures
    public :: plasma_z_function
    public :: plasma_z_function_derivative
    public :: plasma_z_function_second_derivative
    public :: velocity_space_integrate
    public :: maxwell_velocity_distribution
    public :: landau_damping_rate
    
contains

    !---------------------------------------------------------------------------
    ! Plasma dispersion Z-function using simplified but accurate algorithm
    ! 
    ! Z(ζ) = i√π w(ζ) where w is the Faddeeva function
    ! 
    ! Uses different algorithms for different regions of the complex plane
    !---------------------------------------------------------------------------
    function plasma_z_function(zeta) result(z_val)
        complex(real64), intent(in) :: zeta
        complex(real64) :: z_val
        
        real(real64), parameter :: sqrt_pi = 1.77245385090551602730_real64
        real(real64), parameter :: pi = 3.14159265358979323846_real64
        complex(real64) :: w
        real(real64) :: x, y
        
        x = real(zeta, real64)
        y = aimag(zeta)
        
        ! Compute Faddeeva function w(z) = exp(-z²) erfc(-iz)
        w = faddeeva_function(zeta)
        
        ! Z(ζ) = i√π w(ζ)
        z_val = cmplx(0.0_real64, sqrt_pi, real64) * w
        
    end function plasma_z_function
    
    !---------------------------------------------------------------------------
    ! Faddeeva function w(z) = exp(-z²) erfc(-iz)
    ! 
    ! Based on Algorithm 916 by Zaghloul & Ali (2011)
    ! ACM Trans. Math. Software 38, 15
    !---------------------------------------------------------------------------
    recursive function faddeeva_function(z) result(w)
        complex(real64), intent(in) :: z
        complex(real64) :: w
        
        real(real64), parameter :: sqrt_pi = 1.77245385090551602730_real64
        real(real64), parameter :: pi = 3.14159265358979323846_real64
        real(real64) :: x, y, x2, y2
        
        x = real(z, real64)
        y = aimag(z)
        x2 = x * x
        y2 = y * y
        
        if (y < 0.0_real64) then
            ! Use symmetry: w(-z*) = 2*exp(-z²) - w(z)*
            w = 2.0_real64 * exp(-z**2) - conjg(faddeeva_function(conjg(-z)))
            return
        end if
        
        ! For y >= 0, use appropriate algorithm based on |z|
        if (x2 + y2 < 1.0_real64) then
            ! Small |z|: Use Taylor series
            w = faddeeva_taylor_series(z)
        else if (x2 + y2 > 100.0_real64) then
            ! Large |z|: Use asymptotic expansion
            w = faddeeva_asymptotic(z)
        else
            ! Intermediate |z|: Use Laplace continued fraction
            w = faddeeva_continued_fraction(z)
        end if
        
    end function faddeeva_function
    
    !---------------------------------------------------------------------------
    ! Taylor series for Faddeeva function (small |z|)
    !---------------------------------------------------------------------------
    function faddeeva_taylor_series(z) result(w)
        complex(real64), intent(in) :: z
        complex(real64) :: w
        
        real(real64), parameter :: sqrt_pi = 1.77245385090551602730_real64
        complex(real64) :: sum, term, z2
        integer :: n
        real(real64) :: factorial
        
        z2 = z * z
        
        ! w(z) = exp(-z²) * [1 + (2i/√π) * sum_{n=0}^∞ z^(2n+1) / (1·3·5···(2n+1))]
        sum = z
        term = z
        factorial = 1.0_real64
        
        do n = 1, 30
            factorial = factorial * real(2*n-1, real64) / real(2*n+1, real64)
            term = term * z2 * factorial
            sum = sum + term
            if (abs(term) < 1.0e-15_real64 * abs(sum)) exit
        end do
        
        w = exp(-z2) * (1.0_real64 + 2.0_real64 * cmplx(0.0_real64, 1.0_real64, real64) * sum / sqrt_pi)
        
    end function faddeeva_taylor_series
    
    !---------------------------------------------------------------------------
    ! Asymptotic expansion for Faddeeva function (large |z|)
    !---------------------------------------------------------------------------
    function faddeeva_asymptotic(z) result(w)
        complex(real64), intent(in) :: z
        complex(real64) :: w
        
        real(real64), parameter :: sqrt_pi = 1.77245385090551602730_real64
        complex(real64) :: sum, term, inv_z2
        integer :: n
        real(real64) :: coeff
        
        inv_z2 = 1.0_real64 / (z * z)
        
        ! w(z) ≈ i/(√π z) * [1 + sum_{n=1}^∞ (1·3·5···(2n-1))/(2z²)^n]
        sum = cmplx(1.0_real64, 0.0_real64, real64)
        term = cmplx(1.0_real64, 0.0_real64, real64)
        
        do n = 1, 15
            coeff = real(2*n-1, real64) / 2.0_real64
            term = term * coeff * inv_z2
            sum = sum + term
            if (abs(term) < 1.0e-15_real64 * abs(sum)) exit
        end do
        
        w = cmplx(0.0_real64, 1.0_real64, real64) * sum / (sqrt_pi * z)
        
    end function faddeeva_asymptotic
    
    !---------------------------------------------------------------------------
    ! Laplace continued fraction for Faddeeva function (intermediate |z|)
    !---------------------------------------------------------------------------
    function faddeeva_continued_fraction(z) result(w)
        complex(real64), intent(in) :: z
        complex(real64) :: w
        
        real(real64), parameter :: sqrt_pi = 1.77245385090551602730_real64
        complex(real64) :: an, bn, pn, qn, pn1, qn1, ratio
        integer :: n
        real(real64) :: eps
        
        eps = 1.0e-12_real64
        
        ! Initialize continued fraction using modified Lentz method
        pn1 = cmplx(1.0_real64, 0.0_real64, real64)
        pn = 1.0_real64 / (sqrt_pi * z)
        qn1 = cmplx(0.0_real64, 0.0_real64, real64)
        qn = cmplx(1.0_real64, 0.0_real64, real64)
        
        do n = 1, 100
            if (mod(n, 2) == 1) then
                an = cmplx(real(n, real64) / 2.0_real64, 0.0_real64, real64)
                bn = z
            else
                an = cmplx(real(n, real64) / 2.0_real64, 0.0_real64, real64)
                bn = z
            end if
            
            ratio = pn
            pn = bn * pn + an * pn1
            pn1 = ratio
            
            ratio = qn
            qn = bn * qn + an * qn1
            qn1 = ratio
            
            if (abs(qn) > 1.0e10_real64) then
                ! Rescale to prevent overflow
                pn = pn / 1.0e10_real64
                pn1 = pn1 / 1.0e10_real64
                qn = qn / 1.0e10_real64
                qn1 = qn1 / 1.0e10_real64
            end if
            
            if (n > 1) then
                ratio = pn / qn
                if (abs(ratio - pn1/qn1) < eps * abs(ratio)) exit
            end if
        end do
        
        w = sqrt_pi * pn / qn
        
    end function faddeeva_continued_fraction
    
    !---------------------------------------------------------------------------
    ! First derivative of plasma dispersion Z-function
    ! Z'(ζ) = -2(1 + ζZ(ζ))
    !---------------------------------------------------------------------------
    function plasma_z_function_derivative(zeta) result(z_prime)
        complex(real64), intent(in) :: zeta
        complex(real64) :: z_prime
        complex(real64) :: z_val
        
        ! Get Z-function value
        z_val = plasma_z_function(zeta)
        
        ! Apply recurrence relation
        z_prime = -2.0_real64 * (1.0_real64 + zeta * z_val)
        
    end function plasma_z_function_derivative
    
    !---------------------------------------------------------------------------
    ! Second derivative of plasma dispersion Z-function
    ! Z''(ζ) = -2(Z(ζ) + ζZ'(ζ))
    !---------------------------------------------------------------------------
    function plasma_z_function_second_derivative(zeta) result(z_double_prime)
        complex(real64), intent(in) :: zeta
        complex(real64) :: z_double_prime
        complex(real64) :: z_val, z_prime
        
        ! Get Z-function and its first derivative
        z_val = plasma_z_function(zeta)
        z_prime = plasma_z_function_derivative(zeta)
        
        ! Apply recurrence relation
        z_double_prime = -2.0_real64 * (z_val + zeta * z_prime)
        
    end function plasma_z_function_second_derivative
    
    !---------------------------------------------------------------------------
    ! Velocity space integration using Gaussian quadrature
    !
    ! Integrates a function over velocity space using Gauss-Legendre quadrature
    ! Can optionally include Maxwellian weight
    !---------------------------------------------------------------------------
    function velocity_space_integrate(func, v_min, v_max, n_points) result(integral)
        interface
            function func(v) result(val)
                use iso_fortran_env, only: real64
                real(real64), intent(in) :: v
                real(real64) :: val
            end function func
        end interface
        real(real64), intent(in) :: v_min, v_max
        integer, intent(in) :: n_points
        real(real64) :: integral
        
        real(real64), dimension(:), allocatable :: nodes, weights
        real(real64) :: v, w, sum_val
        integer :: i
        
        ! Allocate arrays for Gauss-Legendre nodes and weights
        allocate(nodes(n_points))
        allocate(weights(n_points))
        
        ! Get Gauss-Legendre quadrature points and weights
        call gauss_legendre_quadrature(n_points, nodes, weights)
        
        ! Transform from [-1,1] to [v_min, v_max]
        sum_val = 0.0_real64
        do i = 1, n_points
            ! Transform node from [-1,1] to [v_min, v_max]
            v = 0.5_real64 * ((v_max - v_min) * nodes(i) + (v_max + v_min))
            ! Transform weight
            w = 0.5_real64 * (v_max - v_min) * weights(i)
            ! Evaluate function and accumulate
            sum_val = sum_val + w * func(v)
        end do
        
        integral = sum_val
        
        deallocate(nodes)
        deallocate(weights)
        
    end function velocity_space_integrate
    
    !---------------------------------------------------------------------------
    ! Gauss-Legendre quadrature nodes and weights
    !
    ! Computes nodes and weights for n-point Gauss-Legendre quadrature on [-1,1]
    !---------------------------------------------------------------------------
    subroutine gauss_legendre_quadrature(n, nodes, weights)
        integer, intent(in) :: n
        real(real64), dimension(n), intent(out) :: nodes, weights
        
        real(real64), parameter :: pi = 3.14159265358979323846_real64
        real(real64) :: x, x_old, p1, p2, p3, dp
        integer :: i, j, m
        real(real64) :: eps
        
        eps = 1.0e-14_real64
        m = (n + 1) / 2
        
        ! Compute roots of Legendre polynomial and weights
        do i = 1, m
            ! Initial guess for root
            x = cos(pi * (real(i, real64) - 0.25_real64) / (real(n, real64) + 0.5_real64))
            
            ! Newton-Raphson iteration to find root
            do
                x_old = x
                p1 = 1.0_real64
                p2 = 0.0_real64
                
                ! Recurrence relation for Legendre polynomial
                do j = 1, n
                    p3 = p2
                    p2 = p1
                    p1 = ((2.0_real64 * real(j, real64) - 1.0_real64) * x * p2 - &
                          (real(j, real64) - 1.0_real64) * p3) / real(j, real64)
                end do
                
                ! Derivative of Legendre polynomial
                dp = real(n, real64) * (x * p1 - p2) / (x * x - 1.0_real64)
                
                ! Newton step
                x = x_old - p1 / dp
                
                if (abs(x - x_old) < eps) exit
            end do
            
            ! Store nodes (symmetric)
            nodes(i) = -x
            nodes(n + 1 - i) = x
            
            ! Compute weights
            weights(i) = 2.0_real64 / ((1.0_real64 - x * x) * dp * dp)
            weights(n + 1 - i) = weights(i)
        end do
        
    end subroutine gauss_legendre_quadrature
    
    !---------------------------------------------------------------------------
    ! Maxwell velocity distribution function
    !---------------------------------------------------------------------------
    function maxwell_velocity_distribution(v, v_thermal) result(f)
        real(real64), intent(in) :: v         ! Velocity
        real(real64), intent(in) :: v_thermal ! Thermal velocity
        real(real64) :: f
        
        real(real64), parameter :: pi = 3.14159265358979323846_real64
        
        ! Placeholder - simple Maxwellian
        f = (1.0_real64 / (v_thermal * sqrt(pi))) * exp(-(v/v_thermal)**2)
        
    end function maxwell_velocity_distribution
    
    !---------------------------------------------------------------------------
    ! Landau damping rate calculation (PLACEHOLDER)
    !---------------------------------------------------------------------------
    function landau_damping_rate(omega, k, v_thermal) result(gamma)
        real(real64), intent(in) :: omega      ! Wave frequency
        real(real64), intent(in) :: k          ! Wave number
        real(real64), intent(in) :: v_thermal ! Thermal velocity
        real(real64) :: gamma
        
        ! Placeholder - returns simple estimate
        gamma = -omega * exp(-(omega/(k*v_thermal))**2)
        
    end function landau_damping_rate

end module kilca_plasma_physics_m