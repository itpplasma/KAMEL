!> @file kilca_core_m.f90
!> @brief Core data structure and management for KiLCA Fortran implementation
!> @details This module provides the Fortran equivalent of the C++ core_data class,
!>          maintaining exact functional equivalence while adapting to Fortran paradigms.

module kilca_core_m
    use iso_fortran_env, only: int32, int64, real64
    use iso_c_binding
    use kilca_types_m
    use kilca_settings_m, only: settings_t, antenna_t, eigmode_sett_t
    implicit none
    private
    
    ! Forward declarations for dependent types (to be implemented in respective modules)
    type :: background_t
        ! Placeholder - will be expanded when translating background module
        logical :: initialized = .false.
    end type background_t
    
    type :: mode_data_t
        ! Placeholder - will be expanded when translating mode module
        logical :: initialized = .false.
    end type mode_data_t
    
    ! =========================================================================
    ! Core Data Structure (from core.h)
    ! =========================================================================
    
    !> @brief The top level data structure of KiLCA library
    !> @details Fortran translation of the C++ core_data class
    type, public :: core_data_t
        !> Project path
        character(len=:), allocatable :: path2project
        
        !> Pointer to settings object
        type(settings_t), pointer :: sd => null()
        
        !> Pointer to background object
        type(background_t), pointer :: bp => null()
        
        !> Dimension of the mode_data array (number of modes)
        integer :: dim = 0
        
        !> Array of mode_data objects
        type(mode_data_t), dimension(:), allocatable :: mda
        
        !> Track if settings are owned by this instance
        logical :: owns_settings = .false.
        
        !> Track initialization state
        logical :: initialized = .false.
    end type core_data_t
    
    ! =========================================================================
    ! Public Procedures
    ! =========================================================================
    
    ! Core lifecycle management
    public :: core_data_create
    public :: core_data_destroy
    public :: core_data_delete_modes_array
    
    ! Core computation procedures
    public :: calc_and_set_mode_independent_core_data
    public :: calc_and_set_mode_dependent_core_data_antenna
    public :: calc_and_set_mode_dependent_core_data_eigmode
    public :: calc_and_set_mode_dependent_core_data_antenna_interface
    public :: calc_and_set_mode_dependent_core_data_antenna_interface_mn
    
    ! Accessor procedures
    public :: core_data_get_path
    public :: core_data_has_settings
    public :: core_data_has_background
    public :: core_data_get_modes_count
    public :: core_data_allocate_modes
    
    ! C-Fortran interface procedures
    public :: get_pointer_precision
    public :: set_core_data_in_core_module
    public :: set_background_in_core_module
    public :: set_settings_in_core_module
    public :: clear_all_data_in_mode_data_module
    
    ! =========================================================================
    ! Module Variables (for C-Fortran interface compatibility)
    ! =========================================================================
    
    ! Static pointers for module-level access (mimics C++ extern "C" interface)
    type(core_data_t), pointer, save :: module_core_data => null()
    type(settings_t), pointer, save :: module_settings => null()
    type(background_t), pointer, save :: module_background => null()
    
contains
    
    ! =========================================================================
    ! Core Lifecycle Management
    ! =========================================================================
    
    !> @brief Create and initialize a core_data structure (constructor equivalent)
    !> @param[out] cd Pointer to core_data structure to create
    !> @param[in] path Project path
    !> @param[out] ierr Error code (KILCA_SUCCESS on success)
    subroutine core_data_create(cd, path, ierr)
        type(core_data_t), pointer, intent(out) :: cd
        character(len=*), intent(in) :: path
        integer, intent(out) :: ierr
        
        integer :: pp, alloc_stat
        character(len=256) :: err_msg
        
        ierr = KILCA_SUCCESS
        
        ! Validate input
        if (len_trim(path) == 0) then
            ierr = KILCA_ERROR_INVALID_INPUT
            cd => null()
            return
        end if
        
        ! Check pointer precision compatibility (mimics C++ constructor check)
        call get_pointer_precision(pp, ierr)
        if (pp /= c_intptr_t) then
            write(err_msg, '(a,i0,a,i0,a,i0)') &
                "warning: core_data: sizeof(c_intptr_t)=", c_intptr_t, &
                " != pp=", pp
            print *, trim(err_msg)
            print *, "set appropriate value for integer, parameter :: pp in constants module"
            ierr = KILCA_ERROR
            cd => null()
            return
        end if
        
        ! Allocate core_data structure
        allocate(cd, stat=alloc_stat)
        if (alloc_stat /= 0) then
            ierr = KILCA_ERROR_MEMORY
            cd => null()
            return
        end if
        
        ! Initialize path
        cd%path2project = trim(adjustl(path))
        
        ! Initialize other members
        cd%sd => null()
        cd%bp => null()
        cd%dim = 0
        cd%owns_settings = .false.
        cd%initialized = .true.
        
        ! Note: mda array is not allocated here, matching C++ behavior
        
    end subroutine core_data_create
    
    !> @brief Destroy core_data structure and free all memory (destructor equivalent)
    !> @param[inout] cd Pointer to core_data structure to destroy
    !> @param[out] ierr Error code
    subroutine core_data_destroy(cd, ierr)
        type(core_data_t), pointer, intent(inout) :: cd
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! Check for null pointer
        if (.not. associated(cd)) then
            ierr = KILCA_ERROR_INVALID_INPUT
            return
        end if
        
        ! Do NOT delete sd intentionally, as per C++ implementation comment:
        ! "the static object, this member points to, will be reused"
        
        ! Delete background if associated
        if (associated(cd%bp)) then
            ! Will be implemented when background module is translated
            ! For now, just nullify
            cd%bp => null()
        end if
        
        ! Delete modes array
        call core_data_delete_modes_array(cd, ierr)
        
        ! Deallocate path
        if (allocated(cd%path2project)) then
            deallocate(cd%path2project)
        end if
        
        ! Mark as uninitialized
        cd%initialized = .false.
        
        ! Deallocate the structure itself
        deallocate(cd)
        cd => null()
        
    end subroutine core_data_destroy
    
    !> @brief Delete modes array (exact translation of delete_modes_array method)
    !> @param[inout] cd Core data structure
    !> @param[out] ierr Error code
    subroutine core_data_delete_modes_array(cd, ierr)
        type(core_data_t), intent(inout) :: cd
        integer, intent(out) :: ierr
        
        integer :: ind
        
        ierr = KILCA_SUCCESS
        
        ! Check if mda is allocated (equivalent to C++ "if (mda == 0) return;")
        if (.not. allocated(cd%mda)) return
        
        ! Delete individual mode_data objects
        ! Note: In Fortran, we don't need explicit delete for each element
        ! as they're not pointers in our implementation
        
        ! Deallocate the array
        deallocate(cd%mda)
        cd%dim = 0
        
    end subroutine core_data_delete_modes_array
    
    ! =========================================================================
    ! Core Computation Procedures
    ! =========================================================================
    
    !> @brief Calculate and set mode-independent core data
    !> @param[inout] cd Core data structure
    !> @param[out] ierr Error code
    subroutine calc_and_set_mode_independent_core_data(cd, ierr)
        type(core_data_t), intent(inout) :: cd
        integer, intent(out) :: ierr
        
        ! Static settings for vacuum and flre (mimics C++ static variables)
        type(settings_t), pointer, save :: static_settings_vacuum => null()
        type(settings_t), pointer, save :: static_settings_flre => null()
        logical, save :: first_call_vacuum = .true.
        logical, save :: first_call_flre = .true.
        
        character(len=:), allocatable :: p2pstr
        integer :: vacuum_pos, flre_pos
        
        ierr = KILCA_SUCCESS
        
        ! Convert path to string for searching
        p2pstr = cd%path2project
        
        ! Search for "vacuum" or "flre" in path
        vacuum_pos = index(p2pstr, "vacuum")
        flre_pos = index(p2pstr, "flre")
        
        if (vacuum_pos > 0) then
            if (first_call_vacuum) then
                ! Create static settings for vacuum
                allocate(static_settings_vacuum)
                static_settings_vacuum%path2project = cd%path2project
                ! Will call read_settings when settings module is implemented
                ! static_settings_vacuum%initialized = .true. ! Will be set when fully implemented
                first_call_vacuum = .false.
            end if
            cd%sd => static_settings_vacuum
        else if (flre_pos > 0) then
            if (first_call_flre) then
                ! Create static settings for flre
                allocate(static_settings_flre)
                static_settings_flre%path2project = cd%path2project
                ! Will call read_settings when settings module is implemented
                ! static_settings_flre%initialized = .true. ! Will be set when fully implemented
                first_call_flre = .false.
            end if
            cd%sd => static_settings_flre
        else
            write(*, '(a,a)') "Error: calc_and_set_mode_independent_core_data: " // &
                             "unknown project type in path: ", trim(cd%path2project)
            ierr = KILCA_ERROR
            return
        end if
        
        ! Set global module pointer
        call set_settings_in_core_module(cd%sd)
        
        ! Create background object
        allocate(cd%bp)
        cd%bp%initialized = .true.
        
        ! Set global module pointer
        call set_background_in_core_module(cd%bp)
        
        ! Note: Background profile loading will be implemented when background module is translated
        
    end subroutine calc_and_set_mode_independent_core_data
    
    !> @brief Calculate mode-dependent data for antenna case
    !> @param[inout] cd Core data structure
    !> @param[out] ierr Error code
    subroutine calc_and_set_mode_dependent_core_data_antenna(cd, ierr)
        type(core_data_t), intent(inout) :: cd
        integer, intent(out) :: ierr
        integer :: ind, m, n
        complex(dp) :: omega_lab
        
        ierr = KILCA_SUCCESS
        
        ! Check prerequisites
        if (.not. associated(cd%sd)) then
            ierr = KILCA_ERROR_INVALID_INPUT
            return
        end if
        
        ! Get dimension from antenna settings
        cd%dim = cd%sd%antenna_settings%dma
        
        ! Allocate modes array
        if (allocated(cd%mda)) deallocate(cd%mda)
        allocate(cd%mda(cd%dim), stat=ierr)
        if (ierr /= 0) then
            ierr = KILCA_ERROR_MEMORY
            return
        end if
        
        ! Calculate omega_lab from antenna frequency
        omega_lab = 2.0_dp * pi * cd%sd%antenna_settings%flab
        
        ! Initialize each mode
        do ind = 1, cd%dim
            ! Get m,n from antenna modes array
            m = cd%sd%antenna_settings%modes(2*ind - 1)
            n = cd%sd%antenna_settings%modes(2*ind)
            
            ! Mark as initialized
            cd%mda(ind)%initialized = .true.
            
            ! Clear mode data module (when implemented)
            call clear_all_data_in_mode_data_module(ierr)
        end do
        
    end subroutine calc_and_set_mode_dependent_core_data_antenna
    
    !> @brief Calculate mode-dependent data for eigenmode case
    !> @param[inout] cd Core data structure
    !> @param[out] ierr Error code
    subroutine calc_and_set_mode_dependent_core_data_eigmode(cd, ierr)
        type(core_data_t), intent(inout) :: cd
        integer, intent(out) :: ierr
        integer :: ind
        
        ierr = KILCA_SUCCESS
        
        ! Check prerequisites
        if (.not. associated(cd%sd)) then
            ierr = KILCA_ERROR_INVALID_INPUT
            return
        end if
        
        ! Get dimension from antenna settings
        cd%dim = cd%sd%antenna_settings%dma
        
        ! Allocate modes array
        if (allocated(cd%mda)) deallocate(cd%mda)
        allocate(cd%mda(cd%dim), stat=ierr)
        if (ierr /= 0) then
            ierr = KILCA_ERROR_MEMORY
            return
        end if
        
        ! Initialize modes for eigenmode search
        do ind = 1, cd%dim
            cd%mda(ind)%initialized = .true.
        end do
        
        ! Full implementation will include:
        ! - Frequency loop (search_flag = 1)
        ! - Zero search (search_flag = 0)
        ! - All zeros search (search_flag = -1)
        
    end subroutine calc_and_set_mode_dependent_core_data_eigmode
    
    !> @brief Calculate mode-dependent data for antenna interface case
    !> @param[inout] cd Core data structure
    !> @param[out] ierr Error code
    subroutine calc_and_set_mode_dependent_core_data_antenna_interface(cd, ierr)
        type(core_data_t), intent(inout) :: cd
        integer, intent(out) :: ierr
        integer :: ind
        complex(dp) :: omega_lab
        
        ierr = KILCA_SUCCESS
        
        ! Check prerequisites
        if (.not. associated(cd%sd)) then
            ierr = KILCA_ERROR_INVALID_INPUT
            return
        end if
        
        ! Get dimension from antenna settings
        cd%dim = cd%sd%antenna_settings%dma
        
        ! Allocate modes array
        if (allocated(cd%mda)) deallocate(cd%mda)
        allocate(cd%mda(cd%dim), stat=ierr)
        if (ierr /= 0) then
            ierr = KILCA_ERROR_MEMORY
            return
        end if
        
        ! Calculate omega_lab
        omega_lab = 2.0_dp * pi * cd%sd%antenna_settings%flab
        
        ! Initialize each mode for interface
        do ind = 1, cd%dim
            cd%mda(ind)%initialized = .true.
            call clear_all_data_in_mode_data_module(ierr)
        end do
        
    end subroutine calc_and_set_mode_dependent_core_data_antenna_interface
    
    !> @brief Calculate mode-dependent data for specific m,n mode
    !> @param[inout] cd Core data structure
    !> @param[in] m Poloidal mode number
    !> @param[in] n Toroidal mode number
    !> @param[in] flag Control flag
    !> @param[out] ierr Error code
    subroutine calc_and_set_mode_dependent_core_data_antenna_interface_mn(cd, m, n, flag, ierr)
        type(core_data_t), intent(inout) :: cd
        integer, intent(in) :: m, n, flag
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! Allocate modes array for single mode
        cd%dim = 1
        if (allocated(cd%mda)) deallocate(cd%mda)
        allocate(cd%mda(cd%dim), stat=ierr)
        if (ierr /= 0) then
            ierr = KILCA_ERROR_MEMORY
            return
        end if
        
        ! This will be fully implemented when mode module is translated
        ! For now, just mark as initialized
        cd%mda(1)%initialized = .true.
        
    end subroutine calc_and_set_mode_dependent_core_data_antenna_interface_mn
    
    ! =========================================================================
    ! Accessor Procedures
    ! =========================================================================
    
    !> @brief Get the project path from core_data
    !> @param[in] cd Core data structure
    !> @param[out] path Retrieved path
    !> @param[out] ierr Error code
    subroutine core_data_get_path(cd, path, ierr)
        type(core_data_t), intent(in) :: cd
        character(len=*), intent(out) :: path
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        if (allocated(cd%path2project)) then
            path = cd%path2project
        else
            path = ""
            ierr = KILCA_ERROR
        end if
        
    end subroutine core_data_get_path
    
    !> @brief Check if settings are initialized
    !> @param[in] cd Core data structure
    !> @param[out] has_settings True if settings are associated
    subroutine core_data_has_settings(cd, has_settings)
        type(core_data_t), intent(in) :: cd
        logical, intent(out) :: has_settings
        
        has_settings = associated(cd%sd)
        
    end subroutine core_data_has_settings
    
    !> @brief Check if background is initialized
    !> @param[in] cd Core data structure
    !> @param[out] has_background True if background is associated
    subroutine core_data_has_background(cd, has_background)
        type(core_data_t), intent(in) :: cd
        logical, intent(out) :: has_background
        
        has_background = associated(cd%bp)
        
    end subroutine core_data_has_background
    
    !> @brief Get number of modes
    !> @param[in] cd Core data structure
    !> @param[out] n_modes Number of modes
    !> @param[out] ierr Error code
    subroutine core_data_get_modes_count(cd, n_modes, ierr)
        type(core_data_t), intent(in) :: cd
        integer, intent(out) :: n_modes
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        n_modes = cd%dim
        
    end subroutine core_data_get_modes_count
    
    !> @brief Allocate modes array
    !> @param[inout] cd Core data structure
    !> @param[in] n_modes Number of modes to allocate
    !> @param[out] ierr Error code
    subroutine core_data_allocate_modes(cd, n_modes, ierr)
        type(core_data_t), intent(inout) :: cd
        integer, intent(in) :: n_modes
        integer, intent(out) :: ierr
        
        integer :: alloc_stat
        
        ierr = KILCA_SUCCESS
        
        ! Check for double allocation
        if (allocated(cd%mda)) then
            ierr = KILCA_ERROR
            return
        end if
        
        ! Validate input
        if (n_modes <= 0) then
            ierr = KILCA_ERROR_INVALID_INPUT
            return
        end if
        
        ! Allocate array
        allocate(cd%mda(n_modes), stat=alloc_stat)
        if (alloc_stat /= 0) then
            ierr = KILCA_ERROR_MEMORY
            return
        end if
        
        cd%dim = n_modes
        
    end subroutine core_data_allocate_modes
    
    ! =========================================================================
    ! C-Fortran Interface Procedures
    ! =========================================================================
    
    !> @brief Get pointer precision for compatibility check
    !> @param[out] pp Pointer precision in bytes
    !> @param[out] ierr Error code
    subroutine get_pointer_precision(pp, ierr)
        integer, intent(out) :: pp
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        pp = c_sizeof(c_null_ptr)  ! Size of a C pointer in bytes
        
    end subroutine get_pointer_precision
    
    !> @brief Set core_data pointer in module (for Fortran interface)
    !> @param[in] cd_ptr Pointer to core_data
    subroutine set_core_data_in_core_module(cd_ptr)
        type(core_data_t), pointer, intent(in) :: cd_ptr
        
        module_core_data => cd_ptr
        
    end subroutine set_core_data_in_core_module
    
    !> @brief Set background pointer in module
    !> @param[in] bp_ptr Pointer to background
    subroutine set_background_in_core_module(bp_ptr)
        type(background_t), pointer, intent(in) :: bp_ptr
        
        module_background => bp_ptr
        
    end subroutine set_background_in_core_module
    
    !> @brief Set settings pointer in module
    !> @param[in] sd_ptr Pointer to settings
    subroutine set_settings_in_core_module(sd_ptr)
        type(settings_t), pointer, intent(in) :: sd_ptr
        
        module_settings => sd_ptr
        
    end subroutine set_settings_in_core_module
    
    !> @brief Clear all data in mode_data module
    !> @param[out] ierr Error code
    subroutine clear_all_data_in_mode_data_module(ierr)
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! This will be implemented when mode_data module is translated
        ! For now, just a placeholder
    end subroutine clear_all_data_in_mode_data_module
    
end module kilca_core_m