module kilca_mode_m
    use iso_fortran_env, only: real64, int32
    use kilca_types_m
    use kilca_settings_m
    use kilca_background_m
    implicit none
    
    private
    
    ! RMHD constant missing from kilca_types_m
    ! Note: kilca_types_m has FLRE=3 but C++ has RMHD=3, FLRE=4
    ! For consistency with C++, we define RMHD here
    integer, parameter, public :: PLASMA_MODEL_RMHD = 3
    
    ! Boundary condition constants
    integer, parameter, public :: BOUNDARY_CENTER = 0
    integer, parameter, public :: BOUNDARY_INFINITY = 1
    integer, parameter, public :: BOUNDARY_IDEALWALL = 2
    integer, parameter, public :: BOUNDARY_INTERFACE = 3
    integer, parameter, public :: BOUNDARY_ANTENNA = 4
    
    ! Number of boundary and medium types
    integer, parameter :: Nbc = 5
    integer, parameter :: Nmed = 5
    
    ! String arrays for boundary and medium types
    character(len=64), parameter :: bc_str(Nbc) = [ &
        "center                                                          ", &
        "infinity                                                        ", &
        "idealwall                                                       ", &
        "interface                                                       ", &
        "antenna                                                         " ]
    
    character(len=64), parameter :: med_str(Nmed) = [ &
        "vacuum                                                          ", &
        "medium                                                          ", &
        "imhd                                                            ", &
        "rmhd                                                            ", &
        "flre                                                            " ]
    
    !---------------------------------------------------------------------------
    ! Wave data type - describes properties of a perturbation mode
    !---------------------------------------------------------------------------
    type, public :: wave_data_t
        integer :: m                    ! m harmonic number
        integer :: n                    ! n harmonic number
        complex(real64) :: omega_lab    ! wave frequency in lab frame
        complex(real64) :: omega_mov    ! wave frequency in moving frame
        real(real64) :: r_res           ! resonance location (if present)
        complex(real64) :: det          ! determinant
    end type wave_data_t
    
    !---------------------------------------------------------------------------
    ! Zone type - describes an interval over radius with a plasma model
    !---------------------------------------------------------------------------
    type, public :: zone_t
        real(real64) :: r1              ! left zone boundary
        real(real64) :: r2              ! right zone boundary
        
        integer :: bc1                  ! type of left boundary
        integer :: bc2                  ! type of right boundary
        
        integer :: medium               ! type of plasma model
        integer :: version              ! version of code for model
        integer :: index                ! zone index
        
        ! Pointers to common structures
        type(settings_t), pointer :: sd => null()
        type(background_t), pointer :: bp => null()
        type(wave_data_t), pointer :: wd => null()
        
        character(len=256) :: path      ! path to project
        
        ! Grid and field data
        integer :: dim                  ! radial grid dimension
        real(real64), allocatable :: r(:)        ! radial grid
        real(real64), allocatable :: basis(:)    ! independent basis solutions
        real(real64), allocatable :: EB(:)       ! superposition field
        real(real64), allocatable :: S(:)        ! superposition coefficients
        
        integer :: Nwaves               ! number of waves
        integer :: Ncomps               ! number of (E,B) components
        
        ! ME solution settings
        integer :: max_dim              ! max dimension of radial grid
        real(real64) :: eps_rel         ! relative accuracy
        real(real64) :: eps_abs         ! absolute accuracy
        
        ! ME solution space out settings
        integer :: deg                  ! polynomial degree for spacing
        real(real64) :: reps            ! relative accuracy of sparse solution
        real(real64) :: aeps            ! absolute accuracy of sparse solution
        real(real64) :: step            ! max grid step
        
        integer :: flag_debug           ! debugging flag
    end type zone_t
    
    !---------------------------------------------------------------------------
    ! Mode data type - represents data for a particular perturbation mode
    !---------------------------------------------------------------------------
    type, public :: mode_data_t
        ! Pointers to common structures
        type(settings_t), pointer :: sd => null()
        type(background_t), pointer :: bp => null()
        type(wave_data_t), pointer :: wd => null()
        
        ! Directories for mode
        character(len=256) :: path2linear
        character(len=256) :: path2dispersion
        character(len=256) :: path2poincare
        
        ! Zone data
        integer :: Nzones                           ! number of zones
        type(zone_t), allocatable :: zones(:)       ! zones array
        
        ! Final grid data
        integer :: dim                              ! final radial grid dimension
        real(real64), allocatable :: r(:)           ! final radial grid
        real(real64), allocatable :: EB(:)          ! final EB fields
        real(real64), allocatable :: EB_int(:)      ! EB fields for interpolation
        integer, allocatable :: index(:)            ! zone first point index
        
        ! Stitching system
        integer :: Nc                               ! number of superposition coefficients
        real(real64), allocatable :: A(:)           ! system matrix
        real(real64), allocatable :: B(:)           ! system rhs vector
        real(real64), allocatable :: S(:)           ! superposition coefficients
    end type mode_data_t
    
    ! Public procedures
    public :: wave_data_create, wave_data_destroy
    public :: zone_create, zone_destroy
    public :: mode_data_create, mode_data_destroy
    public :: mode_allocate_zones, mode_set_paths
    public :: ib_index, iEB_index, iFFM_index, iFint_index
    public :: get_boundary_string, get_medium_string
    
contains

    !---------------------------------------------------------------------------
    ! Create wave data structure
    !---------------------------------------------------------------------------
    subroutine wave_data_create(wave, m, n, omega_lab, omega_mov)
        type(wave_data_t), intent(out) :: wave
        integer, intent(in) :: m, n
        complex(real64), intent(in) :: omega_lab, omega_mov
        
        wave%m = m
        wave%n = n
        wave%omega_lab = omega_lab
        wave%omega_mov = omega_mov
        wave%r_res = 0.0_real64
        wave%det = cmplx(0.0_real64, 0.0_real64, real64)
    end subroutine wave_data_create
    
    !---------------------------------------------------------------------------
    ! Destroy wave data structure
    !---------------------------------------------------------------------------
    subroutine wave_data_destroy(wave)
        type(wave_data_t), intent(inout) :: wave
        
        ! Nothing to deallocate for this simple type
        wave%m = 0
        wave%n = 0
    end subroutine wave_data_destroy
    
    !---------------------------------------------------------------------------
    ! Create zone structure
    !---------------------------------------------------------------------------
    subroutine zone_create(zone, r1, r2, bc1, bc2, medium, ierr)
        type(zone_t), intent(out) :: zone
        real(real64), intent(in) :: r1, r2
        integer, intent(in) :: bc1, bc2, medium
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Set boundaries
        zone%r1 = r1
        zone%r2 = r2
        zone%bc1 = bc1
        zone%bc2 = bc2
        zone%medium = medium
        
        ! Initialize other fields
        zone%version = 1
        zone%index = 0
        zone%dim = 0
        zone%Nwaves = 0
        zone%Ncomps = 6  ! Default for E and B components
        
        ! Default ME solution settings
        zone%max_dim = 10000
        zone%eps_rel = 1.0e-12_real64
        zone%eps_abs = 1.0e-12_real64
        
        ! Default space out settings
        zone%deg = 3
        zone%reps = 1.0e-6_real64
        zone%aeps = 1.0e-10_real64
        zone%step = 0.01_real64
        
        zone%flag_debug = 0
        zone%path = ""
        
        ! Validate inputs
        if (r1 >= r2) then
            ierr = -1
            return
        end if
        
        if (bc1 < 0 .or. bc1 >= Nbc) then
            ierr = -2
            return
        end if
        
        if (bc2 < 0 .or. bc2 >= Nbc) then
            ierr = -3
            return
        end if
        
        if (medium < 0 .or. medium >= Nmed) then
            ierr = -4
            return
        end if
        
    end subroutine zone_create
    
    !---------------------------------------------------------------------------
    ! Destroy zone structure
    !---------------------------------------------------------------------------
    subroutine zone_destroy(zone)
        type(zone_t), intent(inout) :: zone
        
        if (allocated(zone%r)) deallocate(zone%r)
        if (allocated(zone%basis)) deallocate(zone%basis)
        if (allocated(zone%EB)) deallocate(zone%EB)
        if (allocated(zone%S)) deallocate(zone%S)
        
        nullify(zone%sd)
        nullify(zone%bp)
        nullify(zone%wd)
    end subroutine zone_destroy
    
    !---------------------------------------------------------------------------
    ! Create mode data structure
    !---------------------------------------------------------------------------
    subroutine mode_data_create(mode, wave, settings, background, ierr)
        type(mode_data_t), intent(out) :: mode
        type(wave_data_t), intent(inout), target :: wave
        type(settings_t), intent(inout), target :: settings
        type(background_t), intent(inout), target :: background
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Link to common structures
        mode%sd => settings
        mode%bp => background
        mode%wd => wave
        
        ! Initialize paths
        mode%path2linear = ""
        mode%path2dispersion = ""
        mode%path2poincare = ""
        
        ! Initialize scalars
        mode%Nzones = 0
        mode%dim = 0
        mode%Nc = 0
        
    end subroutine mode_data_create
    
    !---------------------------------------------------------------------------
    ! Destroy mode data structure
    !---------------------------------------------------------------------------
    subroutine mode_data_destroy(mode)
        type(mode_data_t), intent(inout) :: mode
        integer :: iz
        
        ! Deallocate zones
        if (allocated(mode%zones)) then
            do iz = 1, mode%Nzones
                call zone_destroy(mode%zones(iz))
            end do
            deallocate(mode%zones)
        end if
        
        ! Deallocate arrays
        if (allocated(mode%r)) deallocate(mode%r)
        if (allocated(mode%EB)) deallocate(mode%EB)
        if (allocated(mode%EB_int)) deallocate(mode%EB_int)
        if (allocated(mode%index)) deallocate(mode%index)
        if (allocated(mode%A)) deallocate(mode%A)
        if (allocated(mode%B)) deallocate(mode%B)
        if (allocated(mode%S)) deallocate(mode%S)
        
        ! Nullify pointers
        nullify(mode%sd)
        nullify(mode%bp)
        nullify(mode%wd)
    end subroutine mode_data_destroy
    
    !---------------------------------------------------------------------------
    ! Allocate zones array in mode
    !---------------------------------------------------------------------------
    subroutine mode_allocate_zones(mode, ierr)
        type(mode_data_t), intent(inout) :: mode
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (mode%Nzones <= 0) then
            ierr = -1
            return
        end if
        
        if (allocated(mode%zones)) deallocate(mode%zones)
        allocate(mode%zones(mode%Nzones), stat=ierr)
        
    end subroutine mode_allocate_zones
    
    !---------------------------------------------------------------------------
    ! Set mode paths based on project path and mode parameters
    !---------------------------------------------------------------------------
    subroutine mode_set_paths(mode, path2project, ierr)
        type(mode_data_t), intent(inout) :: mode
        character(len=*), intent(in) :: path2project
        integer, intent(out) :: ierr
        
        character(len=32) :: m_str, n_str, omega_str
        
        ierr = 0
        
        ! Format mode parameters
        write(m_str, '(I0)') mode%wd%m
        write(n_str, '(I0)') mode%wd%n
        write(omega_str, '(ES12.5,"_",ES12.5)') real(mode%wd%omega_lab), aimag(mode%wd%omega_lab)
        
        ! Construct paths
        mode%path2linear = trim(path2project) // "linear/m" // trim(m_str) // &
                          "_n" // trim(n_str) // "_omega" // trim(omega_str) // "/"
        
        mode%path2dispersion = trim(path2project) // "dispersion/m" // trim(m_str) // &
                              "_n" // trim(n_str) // "_omega" // trim(omega_str) // "/"
        
        mode%path2poincare = trim(path2project) // "poincare/m" // trim(m_str) // &
                            "_n" // trim(n_str) // "_omega" // trim(omega_str) // "/"
        
    end subroutine mode_set_paths
    
    !---------------------------------------------------------------------------
    ! Indexing function for basis array
    ! basis(Ncomps, Nwaves, dim) in Fortran, but stored as 1D
    !---------------------------------------------------------------------------
    pure function ib_index(node, sol, comp, part) result(idx)
        integer, intent(in) :: node, sol, comp, part
        integer :: idx
        
        ! Assuming default values for now - would get from zone
        integer, parameter :: Ncomps = 6
        integer, parameter :: Nwaves = 2
        
        ! Fortran 1-based indexing
        ! part: 0=real, 1=imag
        idx = 1 + part + 2*(comp + Ncomps*(sol + Nwaves*node))
    end function ib_index
    
    !---------------------------------------------------------------------------
    ! Indexing function for EB field array
    !---------------------------------------------------------------------------
    pure function iEB_index(node, comp, part) result(idx)
        integer, intent(in) :: node, comp, part
        integer :: idx
        
        integer, parameter :: Ncomps = 6
        
        ! Fortran 1-based indexing
        idx = 1 + part + 2*(comp + Ncomps*node)
    end function iEB_index
    
    !---------------------------------------------------------------------------
    ! Indexing function for final EB array in mode
    !---------------------------------------------------------------------------
    pure function iFFM_index(node, comp, part) result(idx)
        integer, intent(in) :: node, comp, part
        integer :: idx
        
        ! Fortran 1-based indexing
        idx = 1 + part + 2*(comp + 6*node)
    end function iFFM_index
    
    !---------------------------------------------------------------------------
    ! Indexing function for EB_int interpolation array
    !---------------------------------------------------------------------------
    pure function iFint_index(node, comp, part, dim) result(idx)
        integer, intent(in) :: node, comp, part, dim
        integer :: idx
        
        ! Fortran 1-based indexing
        idx = 1 + node + dim*(part + 2*comp)
    end function iFint_index
    
    !---------------------------------------------------------------------------
    ! Get boundary condition string
    !---------------------------------------------------------------------------
    function get_boundary_string(bc) result(str)
        integer, intent(in) :: bc
        character(len=64) :: str
        
        if (bc >= 0 .and. bc < Nbc) then
            str = bc_str(bc + 1)  ! Fortran 1-based
        else
            str = "unknown"
        end if
    end function get_boundary_string
    
    !---------------------------------------------------------------------------
    ! Get medium/plasma model string
    !---------------------------------------------------------------------------
    function get_medium_string(medium) result(str)
        integer, intent(in) :: medium
        character(len=64) :: str
        
        if (medium >= 0 .and. medium < Nmed) then
            str = med_str(medium + 1)  ! Fortran 1-based
        else
            str = "unknown"
        end if
    end function get_medium_string

end module kilca_mode_m