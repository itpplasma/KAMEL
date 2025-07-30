!> Hypergeometric functions module for KiLCA
!! Provides interfaces to GSL hypergeometric functions and custom implementations
!! for specialized cases needed in plasma physics applications
module kilca_hypergeometric_m
    use kilca_types_m
    use iso_fortran_env, only: real64
    use iso_c_binding, only: c_double, c_int
    implicit none
    
    private
    
    ! Public types
    public :: hyperg_1f1_settings_t
    
    ! Public procedures
    public :: hyperg_1f1_gsl
    public :: hyperg_1f1_custom
    public :: hyperg_1f1_kummer_series
    public :: hyperg_1f1_continued_fraction
    public :: hyperg_1f1_quadrature
    public :: hyperg_1f1_select_algorithm
    public :: hyperg_2f1_gsl
    public :: hyperg_0f1_gsl
    public :: hyperg_u_gsl
    
    !> Settings for 1F1 calculation algorithms
    type :: hyperg_1f1_settings_t
        integer :: algorithm = 1                     !< Algorithm selection (1=auto, 2=kummer, 3=cont_frac, 4=quad)
        real(dp) :: tolerance = 1.0e-12_dp          !< Convergence tolerance
        integer :: max_iterations = 1000000         !< Maximum iterations
        logical :: use_acceleration = .false.       !< Use series acceleration (experimental)
        integer :: debug_level = 0                  !< Debug output level
    end type hyperg_1f1_settings_t
    
    ! Algorithm selection constants
    integer, parameter, public :: HYPERG_ALGORITHM_AUTO = 1
    integer, parameter, public :: HYPERG_ALGORITHM_KUMMER = 2
    integer, parameter, public :: HYPERG_ALGORITHM_CONTINUED_FRACTION = 3
    integer, parameter, public :: HYPERG_ALGORITHM_QUADRATURE = 4
    
    ! GSL interface (external C functions)
    interface
        ! GSL hypergeometric functions
        function gsl_sf_hyperg_1f1_double(a, b, z) bind(c, name='gsl_sf_hyperg_1f1')
            import :: c_double
            real(c_double), value :: a, b, z
            real(c_double) :: gsl_sf_hyperg_1f1_double
        end function gsl_sf_hyperg_1f1_double
        
        function gsl_sf_hyperg_2f1_double(a, b, c, z) bind(c, name='gsl_sf_hyperg_2f1')
            import :: c_double
            real(c_double), value :: a, b, c, z
            real(c_double) :: gsl_sf_hyperg_2f1_double
        end function gsl_sf_hyperg_2f1_double
        
        function gsl_sf_hyperg_0f1_double(c, z) bind(c, name='gsl_sf_hyperg_0f1')
            import :: c_double
            real(c_double), value :: c, z
            real(c_double) :: gsl_sf_hyperg_0f1_double
        end function gsl_sf_hyperg_0f1_double
        
        function gsl_sf_hyperg_u_double(a, b, z) bind(c, name='gsl_sf_hyperg_u')
            import :: c_double
            real(c_double), value :: a, b, z
            real(c_double) :: gsl_sf_hyperg_u_double
        end function gsl_sf_hyperg_u_double
        
        ! Custom C++ implementations from hyper1F1.cpp
        function hypergeometric1f1_kummer_ada(b_re, b_im, z_re, z_im, f_re, f_im) &
                bind(c, name='hypergeometric1f1_kummer_ada_')
            import :: c_double, c_int
            real(c_double), intent(in) :: b_re, b_im, z_re, z_im
            real(c_double), intent(out) :: f_re, f_im
            integer(c_int) :: hypergeometric1f1_kummer_ada
        end function hypergeometric1f1_kummer_ada
        
        function hypergeometric1f1_cont_fract_1_inv_ada(b_re, b_im, z_re, z_im, f_re, f_im) &
                bind(c, name='hypergeometric1f1_cont_fract_1_inv_ada_')
            import :: c_double, c_int
            real(c_double), intent(in) :: b_re, b_im, z_re, z_im
            real(c_double), intent(out) :: f_re, f_im
            integer(c_int) :: hypergeometric1f1_cont_fract_1_inv_ada
        end function hypergeometric1f1_cont_fract_1_inv_ada
        
        function hypergeometric1f1_quad(b_re, b_im, z_re, z_im, f_re, f_im) &
                bind(c, name='hypergeometric1f1_quad_')
            import :: c_double, c_int
            real(c_double), intent(in) :: b_re, b_im, z_re, z_im
            real(c_double), intent(out) :: f_re, f_im
            integer(c_int) :: hypergeometric1f1_quad
        end function hypergeometric1f1_quad
    end interface
    
contains

    !---------------------------------------------------------------------------
    ! GSL hypergeometric function interfaces
    !---------------------------------------------------------------------------
    
    !> Confluent hypergeometric function 1F1(a,b,z) using GSL
    subroutine hyperg_1f1_gsl(a, b, z, result, ierr)
        real(dp), intent(in) :: a, b, z
        real(dp), intent(out) :: result
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! GSL 1F1 function (real arguments only)
        result = gsl_sf_hyperg_1f1_double(real(a, c_double), real(b, c_double), real(z, c_double))
        
        ! Check for NaN or infinity
        if (.not. (result > -huge(result) .and. result < huge(result))) then
            ierr = -1
        end if
        
    end subroutine hyperg_1f1_gsl
    
    !> Gauss hypergeometric function 2F1(a,b,c,z) using GSL
    subroutine hyperg_2f1_gsl(a, b, c, z, result, ierr)
        real(dp), intent(in) :: a, b, c, z
        real(dp), intent(out) :: result
        integer, intent(out) :: ierr
        
        ierr = 0
        
        result = gsl_sf_hyperg_2f1_double(real(a, c_double), real(b, c_double), &
                                          real(c, c_double), real(z, c_double))
        
        if (.not. (result > -huge(result) .and. result < huge(result))) then
            ierr = -1
        end if
        
    end subroutine hyperg_2f1_gsl
    
    !> Hypergeometric function 0F1(c,z) using GSL
    subroutine hyperg_0f1_gsl(c, z, result, ierr)
        real(dp), intent(in) :: c, z
        real(dp), intent(out) :: result
        integer, intent(out) :: ierr
        
        ierr = 0
        
        result = gsl_sf_hyperg_0f1_double(real(c, c_double), real(z, c_double))
        
        if (.not. (result > -huge(result) .and. result < huge(result))) then
            ierr = -1
        end if
        
    end subroutine hyperg_0f1_gsl
    
    !> Confluent hypergeometric function U(a,b,z) using GSL
    subroutine hyperg_u_gsl(a, b, z, result, ierr)
        real(dp), intent(in) :: a, b, z
        real(dp), intent(out) :: result
        integer, intent(out) :: ierr
        
        ierr = 0
        
        result = gsl_sf_hyperg_u_double(real(a, c_double), real(b, c_double), real(z, c_double))
        
        if (.not. (result > -huge(result) .and. result < huge(result))) then
            ierr = -1
        end if
        
    end subroutine hyperg_u_gsl
    
    !---------------------------------------------------------------------------
    ! Custom hypergeometric function implementations
    !---------------------------------------------------------------------------
    
    !> Main interface for 1F1(1,b,z) with complex b and z (custom implementation)
    subroutine hyperg_1f1_custom(b, z, result, settings, ierr)
        complex(dp), intent(in) :: b, z
        complex(dp), intent(out) :: result
        type(hyperg_1f1_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        select case(settings%algorithm)
        case(HYPERG_ALGORITHM_AUTO)
            call hyperg_1f1_select_algorithm(b, z, result, settings, ierr)
        case(HYPERG_ALGORITHM_KUMMER)
            call hyperg_1f1_kummer_series(b, z, result, settings, ierr)
        case(HYPERG_ALGORITHM_CONTINUED_FRACTION)
            call hyperg_1f1_continued_fraction(b, z, result, settings, ierr)
        case(HYPERG_ALGORITHM_QUADRATURE)
            call hyperg_1f1_quadrature(b, z, result, settings, ierr)
        case default
            ierr = -1
            result = (0.0_dp, 0.0_dp)
        end select
        
    end subroutine hyperg_1f1_custom
    
    !> Select optimal algorithm based on parameters
    subroutine hyperg_1f1_select_algorithm(b, z, result, settings, ierr)
        complex(dp), intent(in) :: b, z
        complex(dp), intent(out) :: result
        type(hyperg_1f1_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        real(dp) :: ratio_abs
        
        ierr = 0
        ratio_abs = abs(z/b)
        
        if (ratio_abs < 0.1_dp) then
            ! Use Kummer series for small |z/b|
            call hyperg_1f1_kummer_series(b, z, result, settings, ierr)
        else if (ratio_abs > 10.0_dp) then
            ! Use continued fraction for large |z/b|
            call hyperg_1f1_continued_fraction(b, z, result, settings, ierr)
        else
            ! Use quadrature for intermediate values
            call hyperg_1f1_quadrature(b, z, result, settings, ierr)
        end if
        
    end subroutine hyperg_1f1_select_algorithm
    
    !> Kummer series implementation via C++ interface
    subroutine hyperg_1f1_kummer_series(b, z, result, settings, ierr)
        complex(dp), intent(in) :: b, z
        complex(dp), intent(out) :: result
        type(hyperg_1f1_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        real(c_double) :: b_re, b_im, z_re, z_im, f_re, f_im
        integer(c_int) :: c_ierr
        
        ! Extract real and imaginary parts
        b_re = real(b, c_double)
        b_im = aimag(b)
        z_re = real(z, c_double)
        z_im = aimag(z)
        
        ! Call C++ implementation
        c_ierr = hypergeometric1f1_kummer_ada(b_re, b_im, z_re, z_im, f_re, f_im)
        
        ! Convert result
        result = cmplx(f_re, f_im, dp)
        ierr = int(c_ierr)
        
    end subroutine hyperg_1f1_kummer_series
    
    !> Continued fraction implementation via C++ interface
    subroutine hyperg_1f1_continued_fraction(b, z, result, settings, ierr)
        complex(dp), intent(in) :: b, z
        complex(dp), intent(out) :: result
        type(hyperg_1f1_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        real(c_double) :: b_re, b_im, z_re, z_im, f_re, f_im
        integer(c_int) :: c_ierr
        
        b_re = real(b, c_double)
        b_im = aimag(b)
        z_re = real(z, c_double)
        z_im = aimag(z)
        
        c_ierr = hypergeometric1f1_cont_fract_1_inv_ada(b_re, b_im, z_re, z_im, f_re, f_im)
        
        result = cmplx(f_re, f_im, dp)
        ierr = int(c_ierr)
        
    end subroutine hyperg_1f1_continued_fraction
    
    !> Quadrature implementation via C++ interface
    subroutine hyperg_1f1_quadrature(b, z, result, settings, ierr)
        complex(dp), intent(in) :: b, z
        complex(dp), intent(out) :: result
        type(hyperg_1f1_settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        
        real(c_double) :: b_re, b_im, z_re, z_im, f_re, f_im
        integer(c_int) :: c_ierr
        
        b_re = real(b, c_double)
        b_im = aimag(b)
        z_re = real(z, c_double)
        z_im = aimag(z)
        
        c_ierr = hypergeometric1f1_quad(b_re, b_im, z_re, z_im, f_re, f_im)
        
        result = cmplx(f_re, f_im, dp)
        ierr = int(c_ierr)
        
    end subroutine hyperg_1f1_quadrature
    
end module kilca_hypergeometric_m