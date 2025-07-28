!> @file kilca_core_m.f90
!> @brief Core data structure and management for KiLCA Fortran implementation
!> @details This module provides the Fortran equivalent of the C++ core_data class,
!>          maintaining exact functional equivalence while adapting to Fortran paradigms.

module kilca_core_m
    use iso_fortran_env, only: int32, int64, real64
    use iso_c_binding
    use kilca_types_m
    use kilca_settings_m, only: settings_t, antenna_t, eigmode_sett_t, &
                                settings_destroy, settings_deep_copy, settings_validate
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
    public :: core_data_deep_copy
    public :: core_data_copy_assign
    
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
    
    ! Validation procedures
    public :: core_data_validate
    public :: core_data_validate_full
    public :: core_data_check_consistency
    public :: core_data_get_validation_report
    
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
    
    !> @brief Deep copy core_data structure
    !> @param[in] src Source core_data to copy from
    !> @param[out] dst Destination core_data (will be allocated)
    !> @param[out] ierr Error code
    subroutine core_data_deep_copy(src, dst, ierr)
        type(core_data_t), pointer, intent(in) :: src
        type(core_data_t), pointer, intent(out) :: dst
        integer, intent(out) :: ierr
        
        integer :: alloc_stat, i
        
        ierr = KILCA_SUCCESS
        
        ! Check source validity
        if (.not. associated(src)) then
            ierr = KILCA_ERROR_INVALID_INPUT
            dst => null()
            return
        end if
        
        ! Allocate destination
        allocate(dst, stat=alloc_stat)
        if (alloc_stat /= 0) then
            ierr = KILCA_ERROR_MEMORY
            dst => null()
            return
        end if
        
        ! Copy simple members
        dst%path2project = src%path2project
        dst%dim = src%dim
        dst%owns_settings = src%owns_settings
        dst%initialized = src%initialized
        
        ! Deep copy settings if present
        if (associated(src%sd)) then
            if (src%owns_settings) then
                ! Deep copy settings
                call settings_deep_copy(src%sd, dst%sd, ierr)
                if (ierr /= KILCA_SUCCESS) then
                    call core_data_destroy(dst, ierr)
                    dst => null()
                    return
                end if
            else
                ! Just copy pointer (shared settings)
                dst%sd => src%sd
            end if
        end if
        
        ! Deep copy background if present
        if (associated(src%bp)) then
            ! For now, allocate and mark as initialized
            allocate(dst%bp)
            dst%bp%initialized = src%bp%initialized
            ! Full background copy will be implemented when background module is translated
        end if
        
        ! Deep copy mode array if present
        if (allocated(src%mda)) then
            allocate(dst%mda(size(src%mda)), stat=alloc_stat)
            if (alloc_stat /= 0) then
                ierr = KILCA_ERROR_MEMORY
                call core_data_destroy(dst, ierr)
                dst => null()
                return
            end if
            
            ! Copy each mode
            do i = 1, size(src%mda)
                dst%mda(i)%initialized = src%mda(i)%initialized
                ! Full mode copy will be implemented when mode module is translated
            end do
        end if
        
    end subroutine core_data_deep_copy
    
    !> @brief Copy assignment for core_data
    !> @param[inout] dst Destination (existing) core_data
    !> @param[in] src Source core_data to copy from
    !> @param[out] ierr Error code
    subroutine core_data_copy_assign(dst, src, ierr)
        type(core_data_t), pointer, intent(inout) :: dst
        type(core_data_t), pointer, intent(in) :: src
        integer, intent(out) :: ierr
        
        type(core_data_t), pointer :: temp
        
        ierr = KILCA_SUCCESS
        
        ! Check validity
        if (.not. associated(src) .or. .not. associated(dst)) then
            ierr = KILCA_ERROR_INVALID_INPUT
            return
        end if
        
        ! Self-assignment check
        if (associated(dst, src)) then
            return
        end if
        
        ! Create temporary copy
        call core_data_deep_copy(src, temp, ierr)
        if (ierr /= KILCA_SUCCESS) then
            return
        end if
        
        ! Clear existing destination content
        if (associated(dst%sd) .and. dst%owns_settings) then
            call settings_destroy(dst%sd, ierr)
        end if
        if (associated(dst%bp)) then
            deallocate(dst%bp)
        end if
        if (allocated(dst%mda)) then
            deallocate(dst%mda)
        end if
        
        ! Move temp content to dst
        dst%path2project = temp%path2project
        dst%sd => temp%sd
        dst%bp => temp%bp
        dst%dim = temp%dim
        if (allocated(temp%mda)) then
            call move_alloc(temp%mda, dst%mda)
        end if
        dst%owns_settings = temp%owns_settings
        dst%initialized = temp%initialized
        
        ! Clean up temp structure (but not its content which was moved)
        temp%sd => null()
        temp%bp => null()
        deallocate(temp)
        
    end subroutine core_data_copy_assign
    
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
    ! Validation Procedures
    ! =========================================================================
    
    !> @brief Validate core_data structure
    !> @param[in] cd Core data structure to validate
    !> @param[out] is_valid True if structure is valid
    !> @param[out] error_msg Error message if invalid
    !> @param[out] ierr Error code
    subroutine core_data_validate(cd, is_valid, error_msg, ierr)
        type(core_data_t), intent(in) :: cd
        logical, intent(out) :: is_valid
        character(len=*), intent(out) :: error_msg
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        is_valid = .true.
        error_msg = ""
        
        ! Check path validity
        if (.not. allocated(cd%path2project)) then
            is_valid = .false.
            error_msg = "Path not allocated"
            return
        end if
        
        if (len_trim(cd%path2project) == 0) then
            is_valid = .false.
            error_msg = "Empty project path"
            return
        end if
        
        ! Check dimension consistency
        if (cd%dim < 0) then
            is_valid = .false.
            error_msg = "Negative dimension"
            return
        end if
        
        ! Check mode array consistency
        if (cd%dim > 0 .and. .not. allocated(cd%mda)) then
            is_valid = .false.
            error_msg = "Mode array not allocated for dim > 0"
            return
        end if
        
        if (allocated(cd%mda)) then
            if (size(cd%mda) /= cd%dim) then
                is_valid = .false.
                write(error_msg, '(a,i0,a,i0,a)') "Mode array size (", size(cd%mda), &
                                                 ") != dim (", cd%dim, ")"
                return
            end if
        end if
        
        ! Validate settings if present
        if (associated(cd%sd)) then
            ! For now, just check that settings is associated
            ! Full validation would be done separately when settings are populated
            
            ! Check antenna mode consistency only if both are set
            if (cd%dim > 0 .and. cd%sd%antenna_settings%dma > 0) then
                if (cd%sd%antenna_settings%dma /= cd%dim) then
                    is_valid = .false.
                    write(error_msg, '(a,i0,a,i0,a)') "Antenna dma (", cd%sd%antenna_settings%dma, &
                                                     ") != core dim (", cd%dim, ")"
                    return
                end if
            end if
        end if
        
    end subroutine core_data_validate
    
    !> @brief Validate full core_data structure including all subsystems
    !> @param[in] cd Core data structure to validate
    !> @param[out] is_valid True if structure is fully valid
    !> @param[out] error_msg Error message if invalid
    !> @param[out] ierr Error code
    subroutine core_data_validate_full(cd, is_valid, error_msg, ierr)
        type(core_data_t), intent(in) :: cd
        logical, intent(out) :: is_valid
        character(len=*), intent(out) :: error_msg
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! First do basic validation
        call core_data_validate(cd, is_valid, error_msg, ierr)
        if (.not. is_valid) return
        
        ! Check that all required components are present
        if (.not. associated(cd%sd)) then
            is_valid = .false.
            error_msg = "Settings not initialized for full validation"
            return
        end if
        
        if (.not. associated(cd%bp)) then
            is_valid = .false.
            error_msg = "Background not initialized for full validation"
            return
        end if
        
        if (cd%dim > 0 .and. .not. allocated(cd%mda)) then
            is_valid = .false.
            error_msg = "Mode array not allocated for full validation"
            return
        end if
        
        ! Check initialization flag
        if (.not. cd%initialized) then
            is_valid = .false.
            error_msg = "Core data not marked as initialized"
            return
        end if
        
    end subroutine core_data_validate_full
    
    !> @brief Check internal consistency of core_data
    !> @param[in] cd Core data structure to check
    !> @param[out] is_consistent True if internally consistent
    !> @param[out] report Consistency report
    !> @param[out] ierr Error code
    subroutine core_data_check_consistency(cd, is_consistent, report, ierr)
        type(core_data_t), intent(in) :: cd
        logical, intent(out) :: is_consistent
        character(len=*), intent(out) :: report
        integer, intent(out) :: ierr
        
        character(len=1024) :: temp_msg
        integer :: n_issues
        
        ierr = KILCA_SUCCESS
        is_consistent = .true.
        report = "Consistency check report:" // new_line('a')
        n_issues = 0
        
        ! Check ownership consistency
        if (cd%owns_settings .and. .not. associated(cd%sd)) then
            is_consistent = .false.
            n_issues = n_issues + 1
            write(temp_msg, '(i0,a)') n_issues, ". owns_settings=T but sd not associated"
            report = trim(report) // trim(temp_msg) // new_line('a')
        end if
        
        if (.not. cd%owns_settings .and. associated(cd%sd)) then
            ! This is OK - shared settings
            report = trim(report) // "INFO: Using shared settings" // new_line('a')
        end if
        
        ! Check path consistency
        if (associated(cd%sd)) then
            if (allocated(cd%path2project) .and. allocated(cd%sd%path2project)) then
                if (cd%path2project /= cd%sd%path2project) then
                    is_consistent = .false.
                    n_issues = n_issues + 1
                    write(temp_msg, '(i0,a)') n_issues, ". Path mismatch between core and settings"
                    report = trim(report) // trim(temp_msg) // new_line('a')
                end if
            end if
        end if
        
        ! Check static settings consistency
        if (allocated(cd%path2project)) then
            if (index(cd%path2project, "vacuum") > 0 .or. index(cd%path2project, "flre") > 0) then
                if (.not. associated(cd%sd)) then
                    n_issues = n_issues + 1
                    write(temp_msg, '(i0,a)') n_issues, ". vacuum/flre project but no settings"
                    report = trim(report) // trim(temp_msg) // new_line('a')
                    ! Not necessarily inconsistent, just noteworthy
                end if
            end if
        end if
        
        if (n_issues == 0) then
            report = trim(report) // "No consistency issues found" // new_line('a')
        else
            write(temp_msg, '(a,i0,a)') "Found ", n_issues, " consistency issues"
            report = trim(report) // trim(temp_msg) // new_line('a')
        end if
        
    end subroutine core_data_check_consistency
    
    !> @brief Get detailed validation report
    !> @param[in] cd Core data structure
    !> @param[out] report Validation report
    !> @param[out] ierr Error code
    subroutine core_data_get_validation_report(cd, report, ierr)
        type(core_data_t), intent(in) :: cd
        character(len=*), intent(out) :: report
        integer, intent(out) :: ierr
        
        logical :: is_valid, is_consistent
        character(len=1024) :: error_msg, consistency_report
        character(len=256) :: temp_str
        
        ierr = KILCA_SUCCESS
        report = "=== Core Data Validation Report ===" // new_line('a')
        
        ! Basic info
        if (allocated(cd%path2project)) then
            write(temp_str, '(a,a)') "Path: ", trim(cd%path2project)
        else
            temp_str = "Path: <not allocated>"
        end if
        report = trim(report) // trim(temp_str) // new_line('a')
        
        write(temp_str, '(a,i0)') "Dimension: ", cd%dim
        report = trim(report) // trim(temp_str) // new_line('a')
        
        write(temp_str, '(a,l1)') "Initialized: ", cd%initialized
        report = trim(report) // trim(temp_str) // new_line('a')
        
        write(temp_str, '(a,l1)') "Owns settings: ", cd%owns_settings
        report = trim(report) // trim(temp_str) // new_line('a')
        
        ! Component status
        report = trim(report) // new_line('a') // "Components:" // new_line('a')
        
        write(temp_str, '(a,l1)') "  Settings associated: ", associated(cd%sd)
        report = trim(report) // trim(temp_str) // new_line('a')
        
        write(temp_str, '(a,l1)') "  Background associated: ", associated(cd%bp)
        report = trim(report) // trim(temp_str) // new_line('a')
        
        write(temp_str, '(a,l1)') "  Mode array allocated: ", allocated(cd%mda)
        report = trim(report) // trim(temp_str) // new_line('a')
        
        if (allocated(cd%mda)) then
            write(temp_str, '(a,i0)') "  Mode array size: ", size(cd%mda)
            report = trim(report) // trim(temp_str) // new_line('a')
        end if
        
        ! Validation result
        report = trim(report) // new_line('a') // "Validation:" // new_line('a')
        
        call core_data_validate(cd, is_valid, error_msg, ierr)
        write(temp_str, '(a,l1)') "  Basic validation: ", is_valid
        report = trim(report) // trim(temp_str) // new_line('a')
        if (.not. is_valid) then
            write(temp_str, '(a,a)') "  Error: ", trim(error_msg)
            report = trim(report) // trim(temp_str) // new_line('a')
        end if
        
        ! Consistency check
        call core_data_check_consistency(cd, is_consistent, consistency_report, ierr)
        write(temp_str, '(a,l1)') "  Consistency: ", is_consistent
        report = trim(report) // trim(temp_str) // new_line('a')
        
        report = trim(report) // new_line('a') // "=== End of Report ===" // new_line('a')
        
    end subroutine core_data_get_validation_report
    
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