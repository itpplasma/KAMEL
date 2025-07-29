module kilca_constants_m
    use iso_fortran_env, only: real64, int32
    implicit none
    private
    
    ! Mathematical constants with exact precision
    real(real64), parameter, public :: pi = 3.141592653589793238462643383279502884197_real64
    real(real64), parameter, public :: eul = 0.5772156649015328606065120900824024310422_real64
    real(real64), parameter, public :: sqrt2pi = sqrt(2.0_real64 * pi)
    
    ! Physical constants with exact precision
    real(real64), parameter, public :: boltz = 1.60216428e-12_real64  ! erg/eV
    real(real64), parameter, public :: c_light = 29979245800.0_real64  ! cm/s (CGS units)
    real(real64), parameter, public :: m_p = 1.67262158e-24_real64  ! proton mass in g
    real(real64), parameter, public :: m_e = 9.10938185917485e-28_real64  ! electron mass in g
    real(real64), parameter, public :: e_charge = 4.8032e-10_real64  ! elementary charge in esu (CGS)
    real(real64), parameter, public :: gamma_adiabatic = 5.0_real64 / 3.0_real64  ! adiabatic constant
    
    ! Complex constants
    complex(real64), parameter, public :: cmplx_zero = cmplx(0.0_real64, 0.0_real64, real64)
    complex(real64), parameter, public :: cmplx_one = cmplx(1.0_real64, 0.0_real64, real64)
    complex(real64), parameter, public :: cmplx_i = cmplx(0.0_real64, 1.0_real64, real64)
    
    ! Code settings constants (from code_settings.h)
    integer(int32), parameter, public :: DEBUG_FLAG = 0
    integer(int32), parameter, public :: CALC_EIGEN_DECOMPOSITION = 0
    integer(int32), parameter, public :: SORT_DISPERSION_PROFILES = 0
    integer(int32), parameter, public :: USE_JACOBIAN_IN_ODE_SOLVER = 1
    integer(int32), parameter, public :: USE_SPLINES_IN_RHS_EVALUATION = 0
    
    ! Solver method constants (from solver/method.h)
    integer(int32), parameter, public :: NORMAL_EXACT = 0
    integer(int32), parameter, public :: EIG_EXP_EXACT = 1
    integer(int32), parameter, public :: EIG_EXP_PHYS = 2
    
    ! Plasma model constants (from mode/zone.h)
    integer(int32), parameter, public :: PLASMA_MODEL_VACUUM = 0
    integer(int32), parameter, public :: PLASMA_MODEL_MEDIUM = 1
    integer(int32), parameter, public :: PLASMA_MODEL_IMHD = 2
    integer(int32), parameter, public :: PLASMA_MODEL_RMHD = 3
    integer(int32), parameter, public :: PLASMA_MODEL_FLRE = 4
    
    ! Boundary condition constants (from mode/zone.h)
    integer(int32), parameter, public :: BOUNDARY_CENTER = 0
    integer(int32), parameter, public :: BOUNDARY_INFINITY = 1
    integer(int32), parameter, public :: BOUNDARY_IDEALWALL = 2
    integer(int32), parameter, public :: BOUNDARY_INTERFACE = 3
    integer(int32), parameter, public :: BOUNDARY_ANTENNA = 4
    
    ! String length constants
    integer(int32), parameter, public :: MAX_PATH_LENGTH = 1024
    integer(int32), parameter, public :: MAX_NAME_LENGTH = 64
    
    ! Numerical precision constants
    real(real64), parameter, public :: DEFAULT_REL_TOL = 1.0e-12_real64
    real(real64), parameter, public :: DEFAULT_ABS_TOL = 1.0e-12_real64
    real(real64), parameter, public :: INTEGRATION_TOL = 1.0e-3_real64
    real(real64), parameter, public :: ROOT_FINDING_TOL = 1.0e-8_real64
    real(real64), parameter, public :: SMALL_NUMBER = 1.0e-16_real64
    real(real64), parameter, public :: COLLISION_FREQ_FACTOR = 1.4e-7_real64
    
    ! Coordinate system constants
    integer(int32), parameter, public :: COORD_CYLINDRICAL = 0
    integer(int32), parameter, public :: COORD_CARTESIAN = 1
    integer(int32), parameter, public :: COORD_SPHERICAL = 2
    
    ! Component naming for cylindrical coordinates
    character(len=3), parameter, public :: CYL_COMPONENTS = "rtz"
    
    ! Alias names for compatibility with C++ code
    ! O -> cmplx_zero
    ! E -> cmplx_one  
    ! I -> cmplx_i
    
end module kilca_constants_m