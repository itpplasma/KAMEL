module kilca_zone_m
    use iso_fortran_env, only: real64, int32, error_unit
    use kilca_constants_m  ! Has BOUNDARY_* and PLASMA_MODEL_* constants
    use kilca_settings_m
    use kilca_background_m
    use kilca_mode_m, only: wave_data_t  ! Only import the type, avoid conflicts
    implicit none
    
    private
    
    ! Use constants from kilca_constants_m - no need to redefine
    
    ! Boundary condition names
    integer, parameter :: NBC = 5
    character(len=64), parameter :: BC_STR(NBC) = [ &
        "center    ", "infinity  ", "idealwall ", "interface ", "antenna   " ]
    
    ! Medium names  
    integer, parameter :: NMED = 5
    character(len=64), parameter :: MED_STR(NMED) = [ &
        "vacuum", "medium", "imhd  ", "rmhd  ", "flre  " ]
    
    ! Enhanced zone data structure - extends the basic zone_t from kilca_types_m
    type, public :: zone_extended_t
        ! Zone boundaries
        real(real64) :: r1 = 0.0_real64              ! Left zone boundary
        real(real64) :: r2 = 1.0_real64              ! Right zone boundary
        
        ! Boundary conditions
        integer :: bc1 = BOUNDARY_CENTER              ! Type of left boundary
        integer :: bc2 = BOUNDARY_INTERFACE           ! Type of right boundary
        
        ! Plasma model
        integer :: medium = PLASMA_MODEL_VACUUM       ! Type of plasma model
        integer :: version = 1                        ! Version of code for the model
        
        ! Zone identification
        integer :: index = 0                          ! Zone index
        
        ! Pointers to related data structures
        type(settings_t), pointer :: settings => null()       ! Settings structure
        type(background_t), pointer :: background => null()   ! Background structure  
        type(wave_data_t), pointer :: wave_data => null()     ! Wave data
        
        ! Project path
        character(len=256) :: path = ""               ! Path to project
        
        ! Radial grid and solution data
        integer :: dim = 0                            ! Radial grid dimension
        real(real64), allocatable :: r(:)            ! Radial grid (allocatable, not pointer)
        real(real64), allocatable :: basis(:)        ! Independent basis solutions
        real(real64), allocatable :: EB(:)           ! Superposition field (system vector)
        real(real64), allocatable :: S(:)            ! Superposition coefficients
        
        ! Wave dimensions
        integer :: Nwaves = 0                         ! Number of waves
        integer :: Ncomps = 0                         ! Number of (E, B) components
        
        ! Solution settings
        integer :: max_dim = 1000                     ! Max dimension of radial grid
        real(real64) :: eps_rel = 1.0e-6_real64       ! Relative accuracy
        real(real64) :: eps_abs = 1.0e-12_real64      ! Absolute accuracy
        
        ! Solution spacing settings
        integer :: deg = 3                            ! Polynomial degree for spacing
        real(real64) :: reps = 1.0e-3_real64          ! Relative accuracy of sparse solution
        real(real64) :: aeps = 1.0e-6_real64          ! Absolute accuracy of sparse solution
        real(real64) :: step = 0.01_real64            ! Max grid step in solution
        
        ! Debug flag
        integer :: flag_debug = 0                     ! Debug flag
    end type zone_extended_t
    
    ! Public procedures
    public :: zone_create
    public :: zone_destroy
    public :: zone_set_boundaries
    public :: zone_set_plasma_model
    public :: zone_set_radial_grid
    public :: zone_allocate_basis_fields
    public :: zone_print
    public :: zone_validate
    public :: zone_get_filename
    public :: zone_determine_type
    public :: zone_ib
    public :: zone_iEB
    
contains

    !---------------------------------------------------------------------------
    ! Create and initialize zone structure
    !---------------------------------------------------------------------------
    subroutine zone_create(zone, settings, background, wave_data, path, index, ierr)
        type(zone_extended_t), intent(out) :: zone
        type(settings_t), intent(in), target :: settings
        type(background_t), intent(in), target :: background
        type(wave_data_t), intent(in), target :: wave_data
        character(len=*), intent(in) :: path
        integer, intent(in) :: index
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Set basic parameters
        zone%index = index
        zone%path = trim(path)
        
        ! Set pointers to related structures
        zone%settings => settings
        zone%background => background
        zone%wave_data => wave_data
        
        ! Initialize default values
        zone%r1 = 0.0_real64
        zone%r2 = 1.0_real64
        zone%bc1 = BOUNDARY_CENTER
        zone%bc2 = BOUNDARY_INTERFACE
        zone%medium = PLASMA_MODEL_VACUUM
        zone%version = 1
        
        ! Initialize dimensions
        zone%dim = 0
        zone%Nwaves = 0
        zone%Ncomps = 0
        
        ! Note: allocatable arrays are automatically initialized as unallocated
        
    end subroutine zone_create
    
    !---------------------------------------------------------------------------
    ! Destroy zone structure and free memory
    !---------------------------------------------------------------------------
    subroutine zone_destroy(zone)
        type(zone_extended_t), intent(inout) :: zone
        
        ! Deallocate arrays if allocated
        if (allocated(zone%r)) deallocate(zone%r)
        if (allocated(zone%basis)) deallocate(zone%basis)
        if (allocated(zone%EB)) deallocate(zone%EB)
        if (allocated(zone%S)) deallocate(zone%S)
        
        ! Nullify pointers (only for actual pointers)
        nullify(zone%settings)
        nullify(zone%background)
        nullify(zone%wave_data)
        ! Note: allocatable arrays are automatically deallocated
        
        ! Reset values
        zone%index = 0
        zone%path = ""
        zone%dim = 0
        zone%Nwaves = 0
        zone%Ncomps = 0
        
    end subroutine zone_destroy
    
    !---------------------------------------------------------------------------
    ! Set zone boundaries and boundary conditions
    !---------------------------------------------------------------------------
    subroutine zone_set_boundaries(zone, r1, r2, bc1, bc2, ierr)
        type(zone_extended_t), intent(inout) :: zone
        real(real64), intent(in) :: r1, r2
        integer, intent(in) :: bc1, bc2
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Validate boundaries
        if (r1 >= r2) then
            ierr = -1
            write(error_unit, '(A,2ES15.8)') &
                'Error: zone_set_boundaries: r1 >= r2: ', r1, r2
            return
        end if
        
        ! Validate boundary condition types
        if (bc1 < 0 .or. bc1 >= NBC) then
            ierr = -2
            write(error_unit, '(A,I0)') &
                'Error: zone_set_boundaries: invalid bc1: ', bc1
            return
        end if
        
        if (bc2 < 0 .or. bc2 >= NBC) then
            ierr = -3
            write(error_unit, '(A,I0)') &
                'Error: zone_set_boundaries: invalid bc2: ', bc2
            return
        end if
        
        ! Set boundaries
        zone%r1 = r1
        zone%r2 = r2
        zone%bc1 = bc1
        zone%bc2 = bc2
        
    end subroutine zone_set_boundaries
    
    !---------------------------------------------------------------------------
    ! Set plasma model type and version
    !---------------------------------------------------------------------------
    subroutine zone_set_plasma_model(zone, medium, version, ierr)
        type(zone_extended_t), intent(inout) :: zone
        integer, intent(in) :: medium, version
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Validate plasma model type
        if (medium < 0 .or. medium >= NMED) then
            ierr = -1
            write(error_unit, '(A,I0)') &
                'Error: zone_set_plasma_model: invalid medium: ', medium
            return
        end if
        
        ! Validate version
        if (version < 1) then
            ierr = -2
            write(error_unit, '(A,I0)') &
                'Error: zone_set_plasma_model: invalid version: ', version
            return
        end if
        
        ! Set plasma model
        zone%medium = medium
        zone%version = version
        
    end subroutine zone_set_plasma_model
    
    !---------------------------------------------------------------------------
    ! Set radial grid for zone
    !---------------------------------------------------------------------------
    subroutine zone_set_radial_grid(zone, dim, r_grid, ierr)
        type(zone_extended_t), intent(inout) :: zone
        integer, intent(in) :: dim
        real(real64), intent(in) :: r_grid(dim)
        integer, intent(out) :: ierr
        
        integer :: i
        
        ierr = 0
        
        ! Validate dimension
        if (dim <= 0) then
            ierr = -1
            write(error_unit, '(A,I0)') &
                'Error: zone_set_radial_grid: invalid dim: ', dim
            return
        end if
        
        ! Deallocate existing grid if present
        if (allocated(zone%r)) deallocate(zone%r)
        
        ! Allocate and copy grid
        allocate(zone%r(dim))
        do i = 1, dim
            zone%r(i) = r_grid(i)
        end do
        
        zone%dim = dim
        
        ! Validate grid is monotonic
        do i = 2, dim
            if (zone%r(i) <= zone%r(i-1)) then
                ierr = -2
                write(error_unit, '(A,I0,2ES15.8)') &
                    'Error: zone_set_radial_grid: non-monotonic grid at ', i, &
                    zone%r(i-1), zone%r(i)
                return
            end if
        end do
        
    end subroutine zone_set_radial_grid
    
    !---------------------------------------------------------------------------
    ! Allocate basis field arrays
    !---------------------------------------------------------------------------
    subroutine zone_allocate_basis_fields(zone, ierr)
        type(zone_extended_t), intent(inout) :: zone
        integer, intent(out) :: ierr
        
        integer :: basis_size, EB_size, S_size
        
        ierr = 0
        
        ! Validate dimensions
        if (zone%dim <= 0 .or. zone%Nwaves <= 0 .or. zone%Ncomps <= 0) then
            ierr = -1
            write(error_unit, '(A,3I0)') &
                'Error: zone_allocate_basis_fields: invalid dimensions: ', &
                zone%dim, zone%Nwaves, zone%Ncomps
            return
        end if
        
        ! Calculate array sizes (factor of 2 for complex numbers)
        basis_size = 2 * zone%Ncomps * zone%Nwaves * zone%dim
        EB_size = 2 * zone%Ncomps * zone%dim
        S_size = zone%Nwaves
        
        ! Deallocate existing arrays if present
        if (allocated(zone%basis)) deallocate(zone%basis)
        if (allocated(zone%EB)) deallocate(zone%EB)
        if (allocated(zone%S)) deallocate(zone%S)
        
        ! Allocate arrays
        allocate(zone%basis(basis_size))
        allocate(zone%EB(EB_size))
        allocate(zone%S(S_size))
        
        ! Initialize arrays to zero
        zone%basis = 0.0_real64
        zone%EB = 0.0_real64
        zone%S = 0.0_real64
        
    end subroutine zone_allocate_basis_fields
    
    !---------------------------------------------------------------------------
    ! Indexing function for basis array (translation of C++ ib function)
    !---------------------------------------------------------------------------
    function zone_ib(zone, node, sol, comp, part) result(index)
        type(zone_extended_t), intent(in) :: zone
        integer, intent(in) :: node, sol, comp, part
        integer :: index
        
        ! Fortran 1-based indexing version of C++ formula:
        ! part + 2*(comp + Ncomps*(sol + Nwaves*(node)))
        index = part + 1 + 2*(comp + zone%Ncomps*(sol + zone%Nwaves*node))
        
    end function zone_ib
    
    !---------------------------------------------------------------------------
    ! Indexing function for EB field array (translation of C++ iEB function)
    !---------------------------------------------------------------------------
    function zone_iEB(zone, node, comp, part) result(index)
        type(zone_extended_t), intent(in) :: zone
        integer, intent(in) :: node, comp, part
        integer :: index
        
        ! Fortran 1-based indexing version of C++ formula:
        ! part + 2*(comp + Ncomps*(node))
        index = part + 1 + 2*(comp + zone%Ncomps*node)
        
    end function zone_iEB
    
    !---------------------------------------------------------------------------
    ! Print zone information
    !---------------------------------------------------------------------------
    subroutine zone_print(zone)
        type(zone_extended_t), intent(in) :: zone
        
        write(*, '(A)') "========================================"
        write(*, '(A)') "Zone Information"
        write(*, '(A)') "========================================"
        write(*, '(A,I0)') "Zone index:        ", zone%index
        write(*, '(A,A)') "Project path:      ", trim(zone%path)
        write(*, '(A,2ES15.8)') "Boundaries (r1,r2): ", zone%r1, zone%r2
        
        if (zone%bc1 >= 0 .and. zone%bc1 < NBC) then
            write(*, '(A,A)') "Left boundary:     ", trim(BC_STR(zone%bc1+1))
        else
            write(*, '(A,I0)') "Left boundary:     ", zone%bc1
        end if
        
        if (zone%bc2 >= 0 .and. zone%bc2 < NBC) then
            write(*, '(A,A)') "Right boundary:    ", trim(BC_STR(zone%bc2+1))
        else
            write(*, '(A,I0)') "Right boundary:    ", zone%bc2
        end if
        
        if (zone%medium >= 0 .and. zone%medium < NMED) then
            write(*, '(A,A)') "Plasma model:      ", trim(MED_STR(zone%medium+1))
        else
            write(*, '(A,I0)') "Plasma model:      ", zone%medium
        end if
        
        write(*, '(A,I0)') "Model version:     ", zone%version
        write(*, '(A,I0)') "Grid dimension:    ", zone%dim
        write(*, '(A,I0)') "Number of waves:   ", zone%Nwaves
        write(*, '(A,I0)') "Number of comps:   ", zone%Ncomps
        write(*, '(A,L1)') "Radial grid allocated: ", allocated(zone%r)
        write(*, '(A,L1)') "Basis fields allocated: ", allocated(zone%basis)
        write(*, '(A)') "========================================"
        
    end subroutine zone_print
    
    !---------------------------------------------------------------------------
    ! Validate zone structure
    !---------------------------------------------------------------------------
    subroutine zone_validate(zone, ierr)
        type(zone_extended_t), intent(in) :: zone
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Check basic parameters
        if (zone%index <= 0) then
            ierr = -1
            return
        end if
        
        if (zone%r1 >= zone%r2) then
            ierr = -2
            return
        end if
        
        if (.not. associated(zone%settings)) then
            ierr = -3
            return
        end if
        
        if (.not. associated(zone%background)) then
            ierr = -4
            return
        end if
        
        if (.not. associated(zone%wave_data)) then
            ierr = -5
            return
        end if
        
        if (len_trim(zone%path) == 0) then
            ierr = -6
            return
        end if
        
    end subroutine zone_validate
    
    !---------------------------------------------------------------------------
    ! Generate zone filename (translation of get_zone_file_name)
    !---------------------------------------------------------------------------
    subroutine zone_get_filename(zone_index, filename, ierr)
        integer, intent(in) :: zone_index
        character(len=*), intent(out) :: filename
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Validate zone index
        if (zone_index <= 0) then
            ierr = -1
            filename = ""
            return
        end if
        
        ! Generate filename
        write(filename, '(A,I0,A)') "zone_", zone_index, ".in"
        
    end subroutine zone_get_filename
    
    !---------------------------------------------------------------------------
    ! Determine zone type from file (stub implementation)
    !---------------------------------------------------------------------------
    subroutine zone_determine_type(filename, zone_type, ierr)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: zone_type
        integer, intent(out) :: ierr
        
        ! Stub implementation - always returns vacuum type
        ierr = 0
        zone_type = PLASMA_MODEL_VACUUM
        
        ! In full implementation, would read file and determine type
        ! based on content or filename patterns
        
    end subroutine zone_determine_type

end module kilca_zone_m