!> @file kilca_types_m.f90
!> @brief Basic type definitions for KiLCA Fortran translation
!> @details This module provides all fundamental type definitions, constants,
!>          and parameters used throughout the KiLCA Fortran implementation.
!>          It directly translates C++ type definitions while maintaining
!>          exact numerical equivalence.

module kilca_types_m
    use iso_fortran_env, only: real32, real64, real128, int8, int16, int32, int64
    use iso_c_binding
    implicit none
    private
    
    ! =========================================================================
    ! Basic Type Definitions (from typedefs.h)
    ! =========================================================================
    
    !> @brief Fortran equivalent of C++ unsigned char (uchar)
    type, public :: uchar_t
        integer(int8) :: val = 0_int8
    end type uchar_t
    
    !> @brief Fortran equivalent of C++ signed char (schar)
    type, public :: schar_t
        integer(int8) :: val = 0_int8
    end type schar_t
    
    ! =========================================================================
    ! Numeric Kind Parameters
    ! =========================================================================
    
    !> @brief Single precision real kind (32-bit)
    integer, parameter, public :: sp = real32
    
    !> @brief Double precision real kind (64-bit) - primary computation precision
    integer, parameter, public :: dp = real64
    
    !> @brief Quadruple precision real kind (128-bit) - for high precision needs
    integer, parameter, public :: qp = real128
    
    !> @brief Standard 32-bit integer kind
    integer, parameter, public :: i32 = int32
    
    !> @brief Standard 64-bit integer kind
    integer, parameter, public :: i64 = int64
    
    ! =========================================================================
    ! String Length Parameters
    ! =========================================================================
    
    !> @brief Maximum path length for file names
    integer, parameter, public :: MAX_PATH_LEN = 1024
    
    !> @brief Maximum length for identifier names
    integer, parameter, public :: MAX_NAME_LEN = 128
    
    !> @brief Maximum line length for text input
    integer, parameter, public :: MAX_LINE_LEN = 256
    
    ! =========================================================================
    ! Error Code Definitions
    ! =========================================================================
    
    !> @brief Success return code
    integer, parameter, public :: KILCA_SUCCESS = 0
    
    !> @brief Generic error
    integer, parameter, public :: KILCA_ERROR = -1
    
    !> @brief Memory allocation error
    integer, parameter, public :: KILCA_ERROR_MEMORY = -2
    
    !> @brief File I/O error
    integer, parameter, public :: KILCA_ERROR_FILE = -3
    
    !> @brief Invalid input parameter
    integer, parameter, public :: KILCA_ERROR_INVALID_INPUT = -4
    
    !> @brief Numerical convergence failure
    integer, parameter, public :: KILCA_ERROR_CONVERGENCE = -5
    
    !> @brief Array bounds error
    integer, parameter, public :: KILCA_ERROR_BOUNDS = -6
    
    !> @brief Not implemented feature
    integer, parameter, public :: KILCA_ERROR_NOT_IMPLEMENTED = -7
    
    !> @brief File format or parsing error
    integer, parameter, public :: KILCA_ERROR_FORMAT = -8
    
    !> @brief File not found error
    integer, parameter, public :: KILCA_ERROR_FILE_NOT_FOUND = -9
    
    !> @brief File permission error  
    integer, parameter, public :: KILCA_ERROR_FILE_PERMISSION = -10
    
    !> @brief File open error
    integer, parameter, public :: KILCA_ERROR_FILE_OPEN = -11
    
    !> @brief Namelist read error
    integer, parameter, public :: KILCA_ERROR_NAMELIST_READ = -12
    
    !> @brief Invalid parameter range error
    integer, parameter, public :: KILCA_ERROR_INVALID_PARAMETER = -13
    
    ! Format detection constants
    integer, parameter, public :: NAMELIST_FORMAT = 1
    integer, parameter, public :: LEGACY_FORMAT = 2
    integer, parameter, public :: FORMAT_AUTO_DETECT = 0
    
    ! =========================================================================
    ! Physical Constants (from constants.h)
    ! =========================================================================
    
    !> @brief Speed of light in cm/s
    real(dp), parameter, public :: c_light = 29979245800.0_dp
    
    !> @brief Boltzmann constant in erg/eV
    real(dp), parameter, public :: k_boltz = 1.60216428e-12_dp
    
    !> @brief Proton mass in g
    real(dp), parameter, public :: m_proton = 1.67262158e-24_dp
    
    !> @brief Electron mass in g
    real(dp), parameter, public :: m_electron = 9.10938185917485e-28_dp
    
    !> @brief Elementary charge in esu
    real(dp), parameter, public :: e_charge = 4.8032e-10_dp
    
    !> @brief Adiabatic constant (gamma)
    real(dp), parameter, public :: gamma_adiabatic = 5.0_dp/3.0_dp
    
    ! =========================================================================
    ! Mathematical Constants (from constants.h)
    ! =========================================================================
    
    !> @brief Pi with high precision
    real(dp), parameter, public :: pi = 3.141592653589793238462643383279502884197_dp
    
    !> @brief Euler's constant (gamma)
    real(dp), parameter, public :: euler = 0.5772156649015328606065120900824024310422_dp
    
    !> @brief Square root of 2*pi
    real(dp), parameter, public :: sqrt_2pi = sqrt(2.0_dp * pi)
    
    ! =========================================================================
    ! Complex Number Constants (from constants.h)
    ! =========================================================================
    
    !> @brief Complex zero: O = (0.0, 0.0)
    complex(dp), parameter, public :: cmplx_zero = (0.0_dp, 0.0_dp)
    
    !> @brief Complex one: E = (1.0, 0.0)
    complex(dp), parameter, public :: cmplx_one = (1.0_dp, 0.0_dp)
    
    !> @brief Complex i: I = (0.0, 1.0)
    complex(dp), parameter, public :: cmplx_i = (0.0_dp, 1.0_dp)
    
    ! =========================================================================
    ! Array Dimension Constants
    ! =========================================================================
    
    !> @brief Maximum number of modes
    integer, parameter, public :: MAX_MODES = 100
    
    !> @brief Maximum number of zones
    integer, parameter, public :: MAX_ZONES = 10
    
    !> @brief Maximum number of grid points
    integer, parameter, public :: MAX_GRID_POINTS = 10000
    
    !> @brief Maximum matrix dimension
    integer, parameter, public :: MAX_MATRIX_DIM = 1000
    
    ! =========================================================================
    ! Debug and Control Flags (from code_settings.h)
    ! =========================================================================
    
    !> @brief Debug output flag (0=off, 1=on)
    integer, parameter, public :: DEBUG_FLAG = 0
    
    !> @brief Calculate eigenvalue decomposition flag
    integer, parameter, public :: CALC_EIGEN_DECOMPOSITION = 0
    
    !> @brief Sort dispersion profiles flag
    integer, parameter, public :: SORT_DISPERSION_PROFILES = 0
    
    !> @brief Use Jacobian in ODE solver flag
    integer, parameter, public :: USE_JACOBIAN_IN_ODE_SOLVER = 1
    
    !> @brief Use splines in RHS evaluation flag
    integer, parameter, public :: USE_SPLINES_IN_RHS_EVALUATION = 0
    
    ! =========================================================================
    ! Plasma Model Constants (from zone.h)
    ! =========================================================================
    
    !> @brief Vacuum plasma model identifier
    integer, parameter, public :: PLASMA_MODEL_VACUUM = 0
    
    !> @brief Homogeneous medium plasma model identifier
    integer, parameter, public :: PLASMA_MODEL_MEDIUM = 1
    
    !> @brief Ideal MHD plasma model identifier
    integer, parameter, public :: PLASMA_MODEL_IMHD = 2
    
    !> @brief FLRE plasma model identifier
    integer, parameter, public :: PLASMA_MODEL_FLRE = 3
    
    ! =========================================================================
    ! Solver Method Constants (from method.h)
    ! =========================================================================
    
    !> @brief Normal exact solution method
    integer, parameter, public :: NORMAL_EXACT = 0
    
    !> @brief Eigenvalue expansion exact method
    integer, parameter, public :: EIG_EXP_EXACT = 1
    
    !> @brief Eigenvalue expansion physical method
    integer, parameter, public :: EIG_EXP_PHYS = 2
    
    ! =========================================================================
    ! Utility Functions
    ! =========================================================================
    
    public :: uchar_from_int, schar_from_int
    public :: int_from_uchar, int_from_schar
    
contains
    
    !> @brief Convert integer to unsigned char type
    pure function uchar_from_int(val) result(uc)
        integer, intent(in) :: val
        type(uchar_t) :: uc
        
        ! Ensure value is in valid range [0, 255]
        if (val < 0) then
            uc%val = 0_int8
        else if (val > 255) then
            uc%val = -1_int8  ! This represents 255 in unsigned interpretation
        else if (val > 127) then
            uc%val = int(val - 256, int8)  ! Convert to signed representation
        else
            uc%val = int(val, int8)
        end if
    end function uchar_from_int
    
    !> @brief Convert integer to signed char type
    pure function schar_from_int(val) result(sc)
        integer, intent(in) :: val
        type(schar_t) :: sc
        
        ! Ensure value is in valid range [-128, 127]
        if (val < -128) then
            sc%val = int(z'80', int8)  ! -128 in two's complement
        else if (val > 127) then
            sc%val = 127_int8
        else
            sc%val = int(val, int8)
        end if
    end function schar_from_int
    
    !> @brief Convert unsigned char to integer
    pure function int_from_uchar(uc) result(val)
        type(uchar_t), intent(in) :: uc
        integer :: val
        
        ! Convert to integer, treating as unsigned
        if (uc%val < 0_int8) then
            val = int(uc%val) + 256
        else
            val = int(uc%val)
        end if
    end function int_from_uchar
    
    !> @brief Convert signed char to integer
    pure function int_from_schar(sc) result(val)
        type(schar_t), intent(in) :: sc
        integer :: val
        
        val = int(sc%val)
    end function int_from_schar
    
end module kilca_types_m