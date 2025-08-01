!> @file kilca_settings_m.f90
!> @brief Settings management module for KiLCA Fortran implementation
!> @details This module provides the Fortran equivalent of the C++ settings class
!>          and all related settings structures (antenna, background, output, eigmode).

module kilca_settings_m
    use iso_fortran_env, only: real64, int32, int64, output_unit
    use iso_c_binding
    use kilca_types_m
    implicit none
    private
    
    ! =========================================================================
    ! Antenna Settings Type (from antenna.h)
    ! =========================================================================
    
    !> @brief Antenna parameters structure
    type, public :: antenna_t
        !> Small radius (cm) of antenna location
        real(dp) :: ra = 0.0_dp
        
        !> Current density layer width
        real(dp) :: wa = 0.0_dp
        
        !> Current in antenna coils (statamp)
        real(dp) :: I0 = 0.0_dp
        
        !> Frequency (Hz) in the laboratory frame
        complex(dp) :: flab = cmplx_zero
        
        !> Dimension of modes array
        integer :: dma = 0
        
        !> Array of modes (m,n pairs)
        integer, dimension(:), allocatable :: modes
        
        !> Debug flag
        integer :: flag_debug = 0
        
        !> Flag for eigenmode search
        integer :: flag_eigmode = 0
    end type antenna_t
    
    ! =========================================================================
    ! Background Settings Type (from back_sett.h)
    ! =========================================================================
    
    !> @brief Background calculation settings
    type, public :: back_sett_t
        !> Machine settings
        real(dp) :: rtor = 0.0_dp         !< Big torus radius (cm) of the machine
        real(dp) :: rp = 0.0_dp           !< Plasma radius (cm)
        real(dp) :: B0 = 0.0_dp           !< Toroidal magnetic field (G) at the center
        
        !> Background field and plasma settings
        character(len=:), allocatable :: path2profiles    !< Path to input background profiles
        integer :: calc_back = 1          !< Sign shows whether profiles from files/interface, value shows how to recalculate
        character(len=:), allocatable :: flag_back        !< Flag for background: normal (full) or homogeneous
        integer :: N = 5                  !< Splines degree: >= NC + 2N+1, where N - order of FLR expansion, must be odd
        real(dp) :: V_gal_sys = 0.0_dp    !< Velocity (cm/c) of a moving frame
        real(dp) :: V_scale = 1.0_dp      !< Scale of the Vz velocity profile: Vz = V_scale*Vz - V_gal_sys
        real(dp) :: m_i = 1.0_dp          !< Ions mass in units of proton mass
        real(dp) :: zele = 1.0_dp         !< Collision coefficient for electrons: = 1.0 for realistic collision frequency
        real(dp) :: zion = 1.0_dp         !< Collision coefficient for ions: = 1.0 for realistic collision frequency
        integer :: flag_debug = 0         !< Flag for debugging mode (additional checks are performed)
        
        !> Particle settings
        real(dp), dimension(:), allocatable :: mass      !< (ions, electrons) masses
        real(dp), dimension(:), allocatable :: charge    !< (ions, electrons) charges
        
        !> Computed profile arrays (derived from input profiles)
        real(dp), dimension(:), allocatable :: density      !< Density profile array
        real(dp), dimension(:), allocatable :: temperature  !< Temperature profile array  
        real(dp), dimension(:), allocatable :: velocity     !< Velocity profile array
        real(dp), dimension(:), allocatable :: profile_grid !< Spatial grid for profiles
        real(dp), dimension(:), allocatable :: q_profile    !< Safety factor profile
        real(dp), dimension(:), allocatable :: pressure     !< Pressure profile
        real(dp), dimension(:), allocatable :: v_thermal    !< Thermal velocity array
        real(dp), dimension(:), allocatable :: nu_collision !< Collision frequency array
        
        !> Other misc parameters
        real(dp) :: huge_factor = 1.0e30_dp    !< Big factor used in special cases
    end type back_sett_t
    
    ! =========================================================================
    ! Output Settings Type (from output_sett.h)
    ! =========================================================================
    
    !> @brief Output control settings
    type, public :: output_sett_t
        !> Output flags
        integer :: flag_background = 0    !< 1 if compute background data, 2 - store
        integer :: flag_emfield = 0       !< 1 if compute em field data, 2 - store
        integer :: flag_additional = 0    !< 1 if compute additional quants, 2 - store
        integer :: flag_dispersion = 0    !< 1 if compute dispersion, 2 - store
        
        !> Quantity flags
        integer :: num_quants = 0         !< Number of flags
        integer, dimension(:), allocatable :: flag_quants  !< Flags for each quantity if compute it
        
        !> Debug flag
        integer :: flag_debug = 0         !< Flag for debugging mode
    end type output_sett_t
    
    ! =========================================================================
    ! Eigenmode Settings Type (from eigmode_sett.h)
    ! =========================================================================
    
    !> @brief Eigenmode search settings
    type, public :: eigmode_sett_t
        !> Output file name
        character(len=:), allocatable :: fname    !< Name of output file
        
        !> Search parameters
        integer :: search_flag = 0        !< Flag specifies the search option
        
        !> Grid parameters: real and imag parts
        integer :: rdim = 100             !< Real dimension
        real(dp) :: rfmin = 0.0_dp        !< Real minimum
        real(dp) :: rfmax = 1.0e9_dp      !< Real maximum
        integer :: idim = 100             !< Imaginary dimension
        real(dp) :: ifmin = -1.0e6_dp     !< Imaginary minimum
        real(dp) :: ifmax = 1.0e6_dp      !< Imaginary maximum
        
        integer :: stop_flag = 0          !< Stop flag
        
        !> Accuracies
        real(dp) :: eps_res = 1.0e-6_dp   !< Residual accuracy
        real(dp) :: eps_abs = 1.0e-8_dp   !< Absolute accuracy
        real(dp) :: eps_rel = 1.0e-6_dp   !< Relative accuracy
        real(dp) :: delta = 1.0e-6_dp     !< Delta for numerical derivative
        
        integer :: test_roots = 0         !< Test roots flag
        integer :: flag_debug = 0         !< Debug flag
        
        !> Starting values for root search
        integer :: Nguess = 0             !< Number of guess values
        integer :: kmin = 1               !< Minimum k value
        integer :: kmax = 10              !< Maximum k value
        complex(dp), dimension(:), allocatable :: fstart  !< Starting frequencies
        
        !> Number of zeros to be found
        integer :: n_zeros = 10           !< Number of zeros to find
        
        !> Winding number flag
        integer :: use_winding = 0        !< Flag for using winding number evaluation
    end type eigmode_sett_t
    
    ! =========================================================================
    ! Main Settings Type (from settings.h)
    ! =========================================================================
    
    !> @brief Top-level settings structure containing all subsettings
    type, public :: settings_t
        !> Project path
        character(len=:), allocatable :: path2project
        
        !> Format detection - which format was used to read settings
        integer :: format_used = FORMAT_AUTO_DETECT
        
        !> Antenna settings
        type(antenna_t) :: antenna_settings
        
        !> Background settings
        type(back_sett_t) :: background_settings
        
        !> Output settings
        type(output_sett_t) :: output_settings
        
        !> Eigenmode settings
        type(eigmode_sett_t) :: eigmode_settings
        
        !> Pointers for C++ compatibility (these point to the embedded structures)
        type(antenna_t), pointer :: as => null()
        type(back_sett_t), pointer :: bs => null()
        type(output_sett_t), pointer :: os => null()
        type(eigmode_sett_t), pointer :: es => null()
    end type settings_t
    
    ! =========================================================================
    ! Public Procedures
    ! =========================================================================
    
    ! Settings lifecycle
    public :: settings_create
    public :: settings_destroy
    public :: settings_deep_copy
    public :: settings_compare
    public :: settings_initialize_defaults
    public :: settings_read_all
    public :: settings_read_all_full
    public :: settings_print_all
    public :: settings_validate
    public :: settings_validate_complete
    public :: settings_validate_consistency
    public :: settings_read_namelist
    public :: read_settings_auto_detect
    
    ! Accessor procedures
    public :: settings_get_antenna
    public :: settings_get_background
    public :: settings_get_output
    public :: settings_get_eigmode
    
    ! Antenna procedures
    public :: antenna_set_parameters
    public :: antenna_read_settings
    public :: antenna_read_settings_full
    public :: antenna_print_settings
    public :: antenna_print_settings_to_unit
    public :: antenna_deep_copy
    public :: antenna_compare
    public :: antenna_initialize_defaults
    public :: antenna_initialize_custom
    public :: antenna_validate
    public :: antenna_settings_set_modes
    public :: antenna_settings_get_modes
    
    ! Background procedures
    public :: back_sett_set_calc_flag
    public :: back_sett_read_settings
    public :: back_sett_read_settings_full
    public :: back_sett_print_settings
    public :: back_sett_print_settings_to_unit
    public :: back_sett_deep_copy
    public :: back_sett_compare
    public :: back_sett_initialize_defaults
    public :: back_sett_validate
    public :: background_settings_compute_derived
    
    ! Output procedures
    public :: output_sett_set_flags
    public :: output_sett_read_settings
    public :: output_sett_read_settings_full
    public :: output_sett_print_settings
    public :: output_sett_print_settings_to_unit
    public :: output_sett_deep_copy
    public :: output_sett_compare
    public :: output_sett_initialize_defaults
    public :: output_sett_validate
    public :: output_settings_set_flag_quants
    public :: output_settings_get_flag_quants
    
    ! Eigenmode procedures
    public :: eigmode_sett_set_search_flag
    public :: eigmode_sett_read_settings
    public :: eigmode_sett_read_settings_full
    public :: eigmode_sett_print_settings
    public :: eigmode_sett_print_settings_to_unit
    public :: eigmode_sett_deep_copy
    public :: eigmode_sett_compare
    public :: eigmode_sett_initialize_defaults
    public :: eigmode_sett_validate
    
    ! C interface procedures
    public :: set_antenna_settings_c
    public :: set_background_settings_c
    public :: set_particles_settings_c
    public :: set_huge_factor_c
    public :: copy_antenna_data_to_antenna_module
    public :: copy_background_data_to_background_module
    
    ! Error handling procedures
    public :: settings_format_error_message
    public :: settings_get_error_name
    public :: settings_validate_with_context
    public :: settings_attempt_recovery
    public :: settings_get_detailed_validation_errors
    public :: antenna_read_settings_with_error_context
    public :: settings_test_allocation_failure
    public :: settings_log_error
    public :: settings_check_error_log_exists
    public :: settings_clear_error_log
    
contains
    
    ! =========================================================================
    ! Settings Lifecycle Management
    ! =========================================================================
    
    !> @brief Create and initialize settings structure
    subroutine settings_create(sd, path, ierr)
        type(settings_t), pointer, intent(out) :: sd
        character(len=*), intent(in) :: path
        integer, intent(out) :: ierr
        
        integer :: alloc_stat
        
        ierr = KILCA_SUCCESS
        
        ! Validate path parameter
        if (len_trim(path) == 0) then
            ierr = KILCA_ERROR_INVALID_INPUT
            sd => null()
            return
        end if
        
        ! Allocate settings structure
        allocate(sd, stat=alloc_stat)
        if (alloc_stat /= 0) then
            ierr = KILCA_ERROR_MEMORY
            sd => null()
            return
        end if
        
        ! Validate and initialize path
        if (len_trim(path) == 0) then
            ierr = KILCA_ERROR_INVALID_INPUT
            deallocate(sd)
            sd => null()
            return
        end if
        
        sd%path2project = trim(adjustl(path))
        
        ! Set up internal pointers to embedded structures
        sd%as => sd%antenna_settings
        sd%bs => sd%background_settings
        sd%os => sd%output_settings
        sd%es => sd%eigmode_settings
        
        ! Initialize with default values (already done by type defaults)
        
    end subroutine settings_create
    
    !> @brief Destroy settings structure and free memory
    subroutine settings_destroy(sd, ierr)
        type(settings_t), pointer, intent(inout) :: sd
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        if (.not. associated(sd)) then
            ierr = KILCA_ERROR_INVALID_INPUT
            return
        end if
        
        ! Deallocate all allocatable arrays
        if (allocated(sd%antenna_settings%modes)) then
            deallocate(sd%antenna_settings%modes)
        end if
        
        ! Deallocate background settings allocatables
        if (allocated(sd%background_settings%mass)) then
            deallocate(sd%background_settings%mass)
        end if
        if (allocated(sd%background_settings%charge)) then
            deallocate(sd%background_settings%charge)
        end if
        
        ! Deallocate output settings allocatables
        if (allocated(sd%output_settings%flag_quants)) then
            deallocate(sd%output_settings%flag_quants)
        end if
        
        ! Deallocate eigenmode settings allocatables
        if (allocated(sd%eigmode_settings%fstart)) then
            deallocate(sd%eigmode_settings%fstart)
        end if
        
        ! Nullify internal pointers
        sd%as => null()
        sd%bs => null()
        sd%os => null()
        sd%es => null()
        
        ! Deallocate the structure
        deallocate(sd)
        sd => null()
        
    end subroutine settings_destroy
    
    !> @brief Deep copy settings structure
    !> @param[in] src Source settings
    !> @param[out] dst Destination settings (will be allocated)
    !> @param[out] ierr Error code
    subroutine settings_deep_copy(src, dst, ierr)
        type(settings_t), pointer, intent(in) :: src
        type(settings_t), pointer, intent(out) :: dst
        integer, intent(out) :: ierr
        
        integer :: alloc_stat
        
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
        
        ! Copy path
        dst%path2project = src%path2project
        
        ! Deep copy antenna settings
        dst%antenna_settings = src%antenna_settings
        if (allocated(src%antenna_settings%modes)) then
            if (allocated(dst%antenna_settings%modes)) deallocate(dst%antenna_settings%modes)
            allocate(dst%antenna_settings%modes(size(src%antenna_settings%modes)))
            dst%antenna_settings%modes = src%antenna_settings%modes
        end if
        
        ! Deep copy background settings
        dst%background_settings = src%background_settings
        if (allocated(src%background_settings%mass)) then
            if (allocated(dst%background_settings%mass)) deallocate(dst%background_settings%mass)
            allocate(dst%background_settings%mass(size(src%background_settings%mass)))
            dst%background_settings%mass = src%background_settings%mass
        end if
        if (allocated(src%background_settings%charge)) then
            if (allocated(dst%background_settings%charge)) deallocate(dst%background_settings%charge)
            allocate(dst%background_settings%charge(size(src%background_settings%charge)))
            dst%background_settings%charge = src%background_settings%charge
        end if
        
        ! Deep copy output settings
        dst%output_settings = src%output_settings
        if (allocated(src%output_settings%flag_quants)) then
            if (allocated(dst%output_settings%flag_quants)) deallocate(dst%output_settings%flag_quants)
            allocate(dst%output_settings%flag_quants(size(src%output_settings%flag_quants)))
            dst%output_settings%flag_quants = src%output_settings%flag_quants
        end if
        
        ! Deep copy eigenmode settings
        dst%eigmode_settings = src%eigmode_settings
        if (allocated(src%eigmode_settings%fstart)) then
            if (allocated(dst%eigmode_settings%fstart)) deallocate(dst%eigmode_settings%fstart)
            allocate(dst%eigmode_settings%fstart(size(src%eigmode_settings%fstart)))
            dst%eigmode_settings%fstart = src%eigmode_settings%fstart
        end if
        
        ! Set up internal pointers
        dst%as => dst%antenna_settings
        dst%bs => dst%background_settings
        dst%os => dst%output_settings
        dst%es => dst%eigmode_settings
        
    end subroutine settings_deep_copy
    
    !> @brief Compare two complete settings for equality
    subroutine settings_compare(sd1, sd2, is_equal, ierr)
        type(settings_t), intent(in) :: sd1, sd2
        logical, intent(out) :: is_equal
        integer, intent(out) :: ierr
        
        logical :: sub_equal
        
        ierr = KILCA_SUCCESS
        is_equal = .false.
        
        ! Compare paths
        if (sd1%path2project /= sd2%path2project) return
        
        ! Compare antenna settings
        call antenna_compare(sd1%antenna_settings, sd2%antenna_settings, sub_equal, ierr)
        if (ierr /= KILCA_SUCCESS) return
        if (.not. sub_equal) return
        
        ! Compare background settings
        call back_sett_compare(sd1%background_settings, sd2%background_settings, sub_equal, ierr)
        if (ierr /= KILCA_SUCCESS) return
        if (.not. sub_equal) return
        
        ! Compare output settings
        call output_sett_compare(sd1%output_settings, sd2%output_settings, sub_equal, ierr)
        if (ierr /= KILCA_SUCCESS) return
        if (.not. sub_equal) return
        
        ! Compare eigenmode settings
        call eigmode_sett_compare(sd1%eigmode_settings, sd2%eigmode_settings, sub_equal, ierr)
        if (ierr /= KILCA_SUCCESS) return
        if (.not. sub_equal) return
        
        is_equal = .true.
        
    end subroutine settings_compare
    
    !> @brief Initialize all settings with defaults and create structure
    subroutine settings_initialize_defaults(sd, path, ierr)
        type(settings_t), pointer, intent(out) :: sd
        character(len=*), intent(in) :: path
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! First create the settings structure
        call settings_create(sd, path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            sd => null()
            return
        end if
        
        ! Initialize all subsystems with defaults
        call antenna_initialize_defaults(sd%antenna_settings, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        call back_sett_initialize_defaults(sd%background_settings, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        call output_sett_initialize_defaults(sd%output_settings, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        call eigmode_sett_initialize_defaults(sd%eigmode_settings, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
    end subroutine settings_initialize_defaults
    
    !> @brief Read all settings from files (mimics C++ read_settings)
    subroutine settings_read_all(sd, ierr)
        type(settings_t), intent(inout) :: sd
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        print *, ">>>>> Reading settings from ", trim(sd%path2project)
        
        ! Read antenna settings
        call antenna_read_settings(sd%antenna_settings, sd%path2project, ierr)
        if (ierr /= KILCA_SUCCESS) return
        call copy_antenna_data_to_antenna_module(sd%as)
        
        ! Read background settings
        call back_sett_read_settings(sd%background_settings, sd%path2project, ierr)
        if (ierr /= KILCA_SUCCESS) return
        call copy_background_data_to_background_module(sd%bs)
        
        ! Read output settings
        call output_sett_read_settings(sd%output_settings, sd%path2project, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Read eigenmode settings
        call eigmode_sett_read_settings(sd%eigmode_settings, sd%path2project, ierr)
        
    end subroutine settings_read_all
    
    !> @brief Print all settings
    subroutine settings_print_all(sd, ierr)
        type(settings_t), intent(in) :: sd
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        print *, "=== KiLCA Settings ==="
        print *, "Project path: ", trim(sd%path2project)
        
        call antenna_print_settings(sd%antenna_settings, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        call back_sett_print_settings(sd%background_settings, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        call output_sett_print_settings(sd%output_settings, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        call eigmode_sett_print_settings(sd%eigmode_settings, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
    end subroutine settings_print_all
    
    !> @brief Validate all settings
    subroutine settings_validate(sd, is_valid, ierr)
        type(settings_t), intent(in) :: sd
        logical, intent(out) :: is_valid
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        is_valid = .true.
        
        ! Validate antenna settings
        if (sd%antenna_settings%ra <= 0.0_dp) then
            is_valid = .false.
            return
        end if
        
        if (sd%antenna_settings%wa <= 0.0_dp) then
            is_valid = .false.
            return
        end if
        
        if (sd%antenna_settings%dma <= 0) then
            is_valid = .false.
            return
        end if
        
        if (.not. allocated(sd%antenna_settings%modes)) then
            is_valid = .false.
            return
        end if
        
        ! Validate background settings
        if (sd%background_settings%calc_back == 0) then
            is_valid = .false.
            return
        end if
        
        if (sd%background_settings%rtor <= 0.0_dp) then
            is_valid = .false.
            return
        end if
        
        if (sd%background_settings%rp <= 0.0_dp) then
            is_valid = .false.
            return
        end if
        
        if (sd%background_settings%B0 <= 0.0_dp) then
            is_valid = .false.
            return
        end if
        
        ! Other validations would go here...
        
    end subroutine settings_validate
    
    ! =========================================================================
    ! Accessor Procedures
    ! =========================================================================
    
    !> @brief Get antenna settings pointer
    subroutine settings_get_antenna(sd, ant, ierr)
        type(settings_t), intent(in), target :: sd
        type(antenna_t), pointer, intent(out) :: ant
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        ant => sd%as
        
    end subroutine settings_get_antenna
    
    !> @brief Get background settings pointer
    subroutine settings_get_background(sd, bs, ierr)
        type(settings_t), intent(in), target :: sd
        type(back_sett_t), pointer, intent(out) :: bs
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        bs => sd%bs
        
    end subroutine settings_get_background
    
    !> @brief Get output settings pointer
    subroutine settings_get_output(sd, os, ierr)
        type(settings_t), intent(in), target :: sd
        type(output_sett_t), pointer, intent(out) :: os
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        os => sd%os
        
    end subroutine settings_get_output
    
    !> @brief Get eigenmode settings pointer
    subroutine settings_get_eigmode(sd, es, ierr)
        type(settings_t), intent(in), target :: sd
        type(eigmode_sett_t), pointer, intent(out) :: es
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        es => sd%es
        
    end subroutine settings_get_eigmode
    
    ! =========================================================================
    ! Antenna Procedures
    ! =========================================================================
    
    !> @brief Set antenna parameters
    subroutine antenna_set_parameters(ant, ra, wa, I0, flab, dma, modes, ierr)
        type(antenna_t), intent(inout) :: ant
        real(dp), intent(in) :: ra, wa, I0
        complex(dp), intent(in) :: flab
        integer, intent(in) :: dma
        integer, dimension(:), intent(in) :: modes
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ant%ra = ra
        ant%wa = wa
        ant%I0 = I0
        ant%flab = flab
        ant%dma = dma
        
        ! Allocate and copy modes array
        if (allocated(ant%modes)) deallocate(ant%modes)
        allocate(ant%modes(size(modes)))
        ant%modes = modes
        
    end subroutine antenna_set_parameters
    
    !> @brief Read antenna settings from file
    subroutine antenna_read_settings(ant, path, ierr)
        type(antenna_t), intent(inout) :: ant  ! NOTE: ant kept for interface compatibility
        character(len=*), intent(in) :: path
        integer, intent(out) :: ierr
        
        character(len=MAX_PATH_LEN) :: filename
        logical :: file_exists
        integer :: unit, iostat
        
        ierr = KILCA_SUCCESS
        
        ! Construct filename
        filename = trim(path) // "antenna.in"
        
        ! Check if file exists
        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            ierr = KILCA_ERROR_FILE
            return
        end if
        
        ! Open and read file (simplified - real implementation would parse properly)
        open(newunit=unit, file=filename, status='old', iostat=iostat)
        if (iostat /= 0) then
            ierr = KILCA_ERROR_FILE
            return
        end if
        
        ! Read parameters (placeholder - actual parsing would go here)
        ! For now, just close the file
        close(unit)
        
    end subroutine antenna_read_settings
    
    !> @brief Print antenna settings to stdout (matches C++ implementation)
    subroutine antenna_print_settings(ant, ierr)
        type(antenna_t), intent(in) :: ant
        integer, intent(out) :: ierr
        
        call antenna_print_settings_to_unit(ant, output_unit, ierr)
        
    end subroutine antenna_print_settings
    
    !> @brief Print antenna settings to specified unit
    subroutine antenna_print_settings_to_unit(ant, unit, ierr)
        type(antenna_t), intent(in) :: ant
        integer, intent(in) :: unit
        integer, intent(out) :: ierr
        integer :: i
        
        ierr = KILCA_SUCCESS
        
        write(unit, '(A)', iostat=ierr) ""
        if (ierr /= 0) return
        write(unit, '(A)', iostat=ierr) "Check for antenna parameters below:"
        if (ierr /= 0) return
        
        write(unit, '(A)', iostat=ierr) ""
        if (ierr /= 0) return
        write(unit, '(A,G0,A)', iostat=ierr) "antenna radius: ", ant%ra, " cm"
        if (ierr /= 0) return
        write(unit, '(A,G0,A)', iostat=ierr) "antenna current layer width: ", ant%wa, " cm"
        if (ierr /= 0) return
        write(unit, '(A,G0,A)', iostat=ierr) "antenna coils current: ", ant%I0, " statamps"
        if (ierr /= 0) return
        write(unit, '(A,G0,A,G0,A)', iostat=ierr) "antenna lab frequency: (", real(ant%flab), ", ", aimag(ant%flab), ") 1/s"
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "dimension of modes array: ", ant%dma
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "flag for debugging mode: ", ant%flag_debug
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "flag for eigmode search: ", ant%flag_eigmode
        if (ierr /= 0) return
        
        write(unit, '(A)', advance='no', iostat=ierr) "array of mode numbers (m, n): "
        if (ierr /= 0) return
        
        if (allocated(ant%modes)) then
            do i = 1, ant%dma
                write(unit, '(A,I0,A,I0,A)', advance='no', iostat=ierr) "(", ant%modes(2*i-1), ", ", ant%modes(2*i), ") "
                if (ierr /= 0) return
            end do
        end if
        
        write(unit, '(A)', iostat=ierr) ""
        if (ierr /= 0) return
        
        ierr = KILCA_SUCCESS
        
    end subroutine antenna_print_settings_to_unit
    
    !> @brief Deep copy antenna settings
    subroutine antenna_deep_copy(src, dst, ierr)
        type(antenna_t), intent(in) :: src
        type(antenna_t), intent(out) :: dst
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! Copy scalar values
        dst%ra = src%ra
        dst%wa = src%wa
        dst%I0 = src%I0
        dst%flab = src%flab
        dst%dma = src%dma
        dst%flag_debug = src%flag_debug
        dst%flag_eigmode = src%flag_eigmode
        
        ! Deep copy allocatable arrays
        if (allocated(src%modes)) then
            if (allocated(dst%modes)) deallocate(dst%modes)
            allocate(dst%modes(size(src%modes)))
            dst%modes = src%modes
        else
            if (allocated(dst%modes)) deallocate(dst%modes)
        end if
        
    end subroutine antenna_deep_copy
    
    !> @brief Compare two antenna settings for equality
    subroutine antenna_compare(ant1, ant2, is_equal, ierr)
        type(antenna_t), intent(in) :: ant1, ant2
        logical, intent(out) :: is_equal
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        is_equal = .false.
        
        ! Compare scalar values with appropriate tolerance for reals
        if (abs(ant1%ra - ant2%ra) > epsilon(ant1%ra)) return
        if (abs(ant1%wa - ant2%wa) > epsilon(ant1%wa)) return
        if (abs(ant1%I0 - ant2%I0) > epsilon(ant1%I0)) return
        if (abs(ant1%flab - ant2%flab) > epsilon(real(ant1%flab))) return
        if (ant1%dma /= ant2%dma) return
        if (ant1%flag_debug /= ant2%flag_debug) return
        if (ant1%flag_eigmode /= ant2%flag_eigmode) return
        
        ! Compare allocatable arrays
        if (allocated(ant1%modes) .neqv. allocated(ant2%modes)) return
        if (allocated(ant1%modes)) then
            if (size(ant1%modes) /= size(ant2%modes)) return
            if (any(ant1%modes /= ant2%modes)) return
        end if
        
        is_equal = .true.
        
    end subroutine antenna_compare
    
    !> @brief Initialize antenna settings with default values
    subroutine antenna_initialize_defaults(ant, ierr)
        type(antenna_t), intent(out) :: ant
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! Set default values (matching type definitions)
        ant%ra = 0.0_dp
        ant%wa = 0.0_dp
        ant%I0 = 0.0_dp
        ant%flab = cmplx_zero
        ant%dma = 0
        ant%flag_debug = 0
        ant%flag_eigmode = 0
        
        ! Ensure arrays are not allocated by default
        if (allocated(ant%modes)) deallocate(ant%modes)
        
    end subroutine antenna_initialize_defaults
    
    !> @brief Initialize antenna settings with custom values
    subroutine antenna_initialize_custom(ant, ra, wa, I0, flab, dma, flag_debug, flag_eigmode, ierr)
        type(antenna_t), intent(out) :: ant
        real(dp), intent(in), optional :: ra, wa, I0
        complex(dp), intent(in), optional :: flab
        integer, intent(in), optional :: dma, flag_debug, flag_eigmode
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! First set defaults
        call antenna_initialize_defaults(ant, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Override with provided values
        if (present(ra)) ant%ra = ra
        if (present(wa)) ant%wa = wa
        if (present(I0)) ant%I0 = I0
        if (present(flab)) ant%flab = flab
        if (present(dma)) ant%dma = dma
        if (present(flag_debug)) ant%flag_debug = flag_debug
        if (present(flag_eigmode)) ant%flag_eigmode = flag_eigmode
        
    end subroutine antenna_initialize_custom
    
    !> @brief Set modes array for antenna settings
    subroutine antenna_settings_set_modes(ant, modes, ierr)
        type(antenna_t), intent(inout) :: ant
        integer, dimension(:), intent(in) :: modes
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! Deallocate existing array if allocated
        if (allocated(ant%modes)) deallocate(ant%modes)
        
        ! Allocate and copy new modes array
        allocate(ant%modes(size(modes)))
        ant%modes = modes
        
        ! Update dma to reflect the number of mode pairs
        ! Each mode pair consists of (m, n), so dma = size(modes) / 2
        if (mod(size(modes), 2) /= 0) then
            ierr = KILCA_ERROR_INVALID_INPUT  ! Must be even number of elements for (m,n) pairs
            return
        end if
        
        ant%dma = size(modes) / 2
        
    end subroutine antenna_settings_set_modes
    
    !> @brief Get modes array from antenna settings  
    subroutine antenna_settings_get_modes(ant, modes_out, ierr)
        type(antenna_t), intent(in) :: ant
        integer, dimension(:), allocatable, intent(out) :: modes_out
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! Check if modes array is allocated
        if (.not. allocated(ant%modes)) then
            ierr = KILCA_ERROR_INVALID_INPUT  ! Modes array not allocated
            return
        end if
        
        ! Allocate output array and copy
        allocate(modes_out(size(ant%modes)))
        modes_out = ant%modes
        
    end subroutine antenna_settings_get_modes
    
    ! =========================================================================
    ! Background Settings Procedures
    ! =========================================================================
    
    !> @brief Set background calculation flag
    subroutine back_sett_set_calc_flag(bs, flag, ierr)
        type(back_sett_t), intent(inout) :: bs
        integer, intent(in) :: flag
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        bs%calc_back = flag
        
    end subroutine back_sett_set_calc_flag
    
    !> @brief Read background settings from file
    subroutine back_sett_read_settings(bs, path, ierr)
        type(back_sett_t), intent(inout) :: bs  ! NOTE: bs kept for interface compatibility
        character(len=*), intent(in) :: path
        integer, intent(out) :: ierr
        
        character(len=MAX_PATH_LEN) :: filename
        
        ierr = KILCA_SUCCESS
        
        ! Construct filename
        filename = trim(path) // "background.in"
        
        ! Placeholder for actual file reading
        
    end subroutine back_sett_read_settings
    
    !> @brief Print background settings to stdout (matches C++ implementation)
    subroutine back_sett_print_settings(bs, ierr)
        type(back_sett_t), intent(in) :: bs
        integer, intent(out) :: ierr
        
        call back_sett_print_settings_to_unit(bs, output_unit, ierr)
        
    end subroutine back_sett_print_settings
    
    !> @brief Print background settings to specified unit
    subroutine back_sett_print_settings_to_unit(bs, unit, ierr)
        type(back_sett_t), intent(in) :: bs
        integer, intent(in) :: unit
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        write(unit, '(A)', iostat=ierr) ""
        if (ierr /= 0) return
        write(unit, '(A)', iostat=ierr) "Check for background parameters below:"
        if (ierr /= 0) return
        
        write(unit, '(A)', iostat=ierr) ""
        if (ierr /= 0) return
        write(unit, '(A,G0,A)', iostat=ierr) "torus big radius: ", bs%rtor, " cm"
        if (ierr /= 0) return
        write(unit, '(A,G0,A)', iostat=ierr) "plasma radius: ", bs%rp, " cm"
        if (ierr /= 0) return
        write(unit, '(A,G0,A)', iostat=ierr) "toroidal magnetic field at the center: ", bs%B0, " G"
        if (ierr /= 0) return
        
        if (allocated(bs%path2profiles)) then
            write(unit, '(A,A)', iostat=ierr) "path to background profiles: ", bs%path2profiles
            if (ierr /= 0) return
        end if
        write(unit, '(A,I0)', iostat=ierr) "flag if recalculate background: ", bs%calc_back
        if (ierr /= 0) return
        if (allocated(bs%flag_back)) then
            write(unit, '(A,A)', iostat=ierr) "flag for background: ", bs%flag_back
            if (ierr /= 0) return
        end if
        write(unit, '(A,I0)', iostat=ierr) "splines degree: ", bs%N
        if (ierr /= 0) return
        write(unit, '(A,G0,A)', iostat=ierr) "velocity of the moving frame: ", bs%V_gal_sys, " cm/s"
        if (ierr /= 0) return
        write(unit, '(A,G0)', iostat=ierr) "scale factor for the Vz velocity profile: ", bs%V_scale
        if (ierr /= 0) return
        write(unit, '(A,G0)', iostat=ierr) "ions mass in units of proton mass: ", bs%m_i
        if (ierr /= 0) return
        write(unit, '(A,G0)', iostat=ierr) "collision coefficient for electrons: ", bs%zele
        if (ierr /= 0) return
        write(unit, '(A,G0)', iostat=ierr) "collision coefficient for ions: ", bs%zion
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "flag for debugging mode: ", bs%flag_debug
        if (ierr /= 0) return
        
        write(unit, '(A)', iostat=ierr) ""
        if (ierr /= 0) return
        
        ierr = KILCA_SUCCESS
        
    end subroutine back_sett_print_settings_to_unit
    
    !> @brief Deep copy background settings
    subroutine back_sett_deep_copy(src, dst, ierr)
        type(back_sett_t), intent(in) :: src
        type(back_sett_t), intent(out) :: dst
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! Copy scalar values
        dst%rtor = src%rtor
        dst%rp = src%rp
        dst%B0 = src%B0
        dst%calc_back = src%calc_back
        dst%N = src%N
        dst%V_gal_sys = src%V_gal_sys
        dst%V_scale = src%V_scale
        dst%m_i = src%m_i
        dst%zele = src%zele
        dst%zion = src%zion
        dst%flag_debug = src%flag_debug
        dst%huge_factor = src%huge_factor
        
        ! Deep copy allocatable strings
        if (allocated(src%path2profiles)) then
            if (allocated(dst%path2profiles)) deallocate(dst%path2profiles)
            dst%path2profiles = src%path2profiles
        else
            if (allocated(dst%path2profiles)) deallocate(dst%path2profiles)
        end if
        
        if (allocated(src%flag_back)) then
            if (allocated(dst%flag_back)) deallocate(dst%flag_back)
            dst%flag_back = src%flag_back
        else
            if (allocated(dst%flag_back)) deallocate(dst%flag_back)
        end if
        
        ! Deep copy allocatable arrays
        if (allocated(src%mass)) then
            if (allocated(dst%mass)) deallocate(dst%mass)
            allocate(dst%mass(size(src%mass)))
            dst%mass = src%mass
        else
            if (allocated(dst%mass)) deallocate(dst%mass)
        end if
        
        if (allocated(src%charge)) then
            if (allocated(dst%charge)) deallocate(dst%charge)
            allocate(dst%charge(size(src%charge)))
            dst%charge = src%charge
        else
            if (allocated(dst%charge)) deallocate(dst%charge)
        end if
        
    end subroutine back_sett_deep_copy
    
    !> @brief Compare two background settings for equality
    subroutine back_sett_compare(bs1, bs2, is_equal, ierr)
        type(back_sett_t), intent(in) :: bs1, bs2
        logical, intent(out) :: is_equal
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        is_equal = .false.
        
        ! Compare scalar values with appropriate tolerance for reals
        if (abs(bs1%rtor - bs2%rtor) > epsilon(bs1%rtor)) return
        if (abs(bs1%rp - bs2%rp) > epsilon(bs1%rp)) return
        if (abs(bs1%B0 - bs2%B0) > epsilon(bs1%B0)) return
        if (bs1%calc_back /= bs2%calc_back) return
        if (bs1%N /= bs2%N) return
        if (abs(bs1%V_gal_sys - bs2%V_gal_sys) > epsilon(bs1%V_gal_sys)) return
        if (abs(bs1%V_scale - bs2%V_scale) > epsilon(bs1%V_scale)) return
        if (abs(bs1%m_i - bs2%m_i) > epsilon(bs1%m_i)) return
        if (abs(bs1%zele - bs2%zele) > epsilon(bs1%zele)) return
        if (abs(bs1%zion - bs2%zion) > epsilon(bs1%zion)) return
        if (bs1%flag_debug /= bs2%flag_debug) return
        if (abs(bs1%huge_factor - bs2%huge_factor) > epsilon(bs1%huge_factor)) return
        
        ! Compare allocatable strings
        if (allocated(bs1%path2profiles) .neqv. allocated(bs2%path2profiles)) return
        if (allocated(bs1%path2profiles)) then
            if (bs1%path2profiles /= bs2%path2profiles) return
        end if
        
        if (allocated(bs1%flag_back) .neqv. allocated(bs2%flag_back)) return
        if (allocated(bs1%flag_back)) then
            if (bs1%flag_back /= bs2%flag_back) return
        end if
        
        ! Compare allocatable arrays
        if (allocated(bs1%mass) .neqv. allocated(bs2%mass)) return
        if (allocated(bs1%mass)) then
            if (size(bs1%mass) /= size(bs2%mass)) return
            if (any(abs(bs1%mass - bs2%mass) > epsilon(bs1%mass))) return
        end if
        
        if (allocated(bs1%charge) .neqv. allocated(bs2%charge)) return
        if (allocated(bs1%charge)) then
            if (size(bs1%charge) /= size(bs2%charge)) return
            if (any(abs(bs1%charge - bs2%charge) > epsilon(bs1%charge))) return
        end if
        
        is_equal = .true.
        
    end subroutine back_sett_compare
    
    !> @brief Initialize background settings with meaningful defaults
    subroutine back_sett_initialize_defaults(bs, ierr)
        type(back_sett_t), intent(out) :: bs
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! Set meaningful default values for a typical tokamak
        bs%rtor = 625.0_dp           ! ASDEX Upgrade major radius (cm)
        bs%rp = 200.0_dp             ! Typical plasma minor radius (cm) 
        bs%B0 = 20000.0_dp           ! Typical toroidal field (G)
        bs%calc_back = 1             ! Calculate background profiles
        bs%N = 5                     ! Spline degree (must be odd)
        bs%V_gal_sys = 0.0_dp        ! No moving frame by default
        bs%V_scale = 1.0_dp          ! Unity velocity scaling
        bs%m_i = 1.0_dp              ! Proton mass units
        bs%zele = 1.0_dp             ! Realistic electron collision frequency
        bs%zion = 1.0_dp             ! Realistic ion collision frequency
        bs%flag_debug = 0            ! No debug by default
        bs%huge_factor = 1.0e30_dp   ! Large factor for special cases
        
        ! Set default strings
        if (allocated(bs%flag_back)) deallocate(bs%flag_back)
        bs%flag_back = "normal"
        
        if (allocated(bs%path2profiles)) deallocate(bs%path2profiles)
        ! Don't set path2profiles by default - should be set explicitly
        
        ! Clean up arrays - should be allocated as needed
        if (allocated(bs%mass)) deallocate(bs%mass)
        if (allocated(bs%charge)) deallocate(bs%charge)
        
    end subroutine back_sett_initialize_defaults
    
    ! =========================================================================
    ! Output Settings Procedures
    ! =========================================================================
    
    !> @brief Set output flags
    subroutine output_sett_set_flags(os, flag_background, flag_emfield, ierr)
        type(output_sett_t), intent(inout) :: os
        integer, intent(in), optional :: flag_background, flag_emfield
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        if (present(flag_background)) os%flag_background = flag_background
        if (present(flag_emfield)) os%flag_emfield = flag_emfield
        
    end subroutine output_sett_set_flags
    
    !> @brief Read output settings from file
    subroutine output_sett_read_settings(os, path, ierr)
        type(output_sett_t), intent(inout) :: os  ! NOTE: os kept for interface compatibility
        character(len=*), intent(in) :: path
        integer, intent(out) :: ierr
        
        character(len=MAX_PATH_LEN) :: filename
        
        ierr = KILCA_SUCCESS
        
        ! Construct filename
        filename = trim(path) // "output.in"
        
        ! Placeholder for actual file reading
        
    end subroutine output_sett_read_settings
    
    !> @brief Print output settings to stdout (matches C++ implementation)
    subroutine output_sett_print_settings(os, ierr)
        type(output_sett_t), intent(in) :: os
        integer, intent(out) :: ierr
        
        call output_sett_print_settings_to_unit(os, output_unit, ierr)
        
    end subroutine output_sett_print_settings
    
    !> @brief Print output settings to specified unit
    subroutine output_sett_print_settings_to_unit(os, unit, ierr)
        type(output_sett_t), intent(in) :: os
        integer, intent(in) :: unit
        integer, intent(out) :: ierr
        integer :: i
        
        ierr = KILCA_SUCCESS
        
        write(unit, '(A)', iostat=ierr) ""
        if (ierr /= 0) return
        write(unit, '(A)', iostat=ierr) "Check for output parameters below:"
        if (ierr /= 0) return
        
        write(unit, '(A)', iostat=ierr) ""
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "if compute background data: ", os%flag_background
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "if compute linear data: ", os%flag_emfield
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "if compute additional quants: ", os%flag_additional
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "if compute dispersion: ", os%flag_dispersion
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "flag for debugging: ", os%flag_debug
        if (ierr /= 0) return
        
        write(unit, '(A,I0)', iostat=ierr) "dimension of flags array for quantities: ", os%num_quants
        if (ierr /= 0) return
        
        write(unit, '(A)', advance='no', iostat=ierr) "array of flags for additional quantities: "
        if (ierr /= 0) return
        
        if (allocated(os%flag_quants)) then
            do i = 1, os%num_quants
                write(unit, '(I0,A)', advance='no', iostat=ierr) os%flag_quants(i), " "
                if (ierr /= 0) return
            end do
        end if
        
        write(unit, '(A)', iostat=ierr) ""
        if (ierr /= 0) return
        
        ierr = KILCA_SUCCESS
        
    end subroutine output_sett_print_settings_to_unit
    
    !> @brief Deep copy output settings
    subroutine output_sett_deep_copy(src, dst, ierr)
        type(output_sett_t), intent(in) :: src
        type(output_sett_t), intent(out) :: dst
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! Copy scalar values
        dst%flag_background = src%flag_background
        dst%flag_emfield = src%flag_emfield
        dst%flag_additional = src%flag_additional
        dst%flag_dispersion = src%flag_dispersion
        dst%num_quants = src%num_quants
        dst%flag_debug = src%flag_debug
        
        ! Deep copy allocatable arrays
        if (allocated(src%flag_quants)) then
            if (allocated(dst%flag_quants)) deallocate(dst%flag_quants)
            allocate(dst%flag_quants(size(src%flag_quants)))
            dst%flag_quants = src%flag_quants
        else
            if (allocated(dst%flag_quants)) deallocate(dst%flag_quants)
        end if
        
    end subroutine output_sett_deep_copy
    
    !> @brief Compare two output settings for equality
    subroutine output_sett_compare(os1, os2, is_equal, ierr)
        type(output_sett_t), intent(in) :: os1, os2
        logical, intent(out) :: is_equal
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        is_equal = .false.
        
        ! Compare scalar values
        if (os1%flag_background /= os2%flag_background) return
        if (os1%flag_emfield /= os2%flag_emfield) return
        if (os1%flag_additional /= os2%flag_additional) return
        if (os1%flag_dispersion /= os2%flag_dispersion) return
        if (os1%num_quants /= os2%num_quants) return
        if (os1%flag_debug /= os2%flag_debug) return
        
        ! Compare allocatable arrays
        if (allocated(os1%flag_quants) .neqv. allocated(os2%flag_quants)) return
        if (allocated(os1%flag_quants)) then
            if (size(os1%flag_quants) /= size(os2%flag_quants)) return
            if (any(os1%flag_quants /= os2%flag_quants)) return
        end if
        
        is_equal = .true.
        
    end subroutine output_sett_compare
    
    !> @brief Initialize output settings with practical defaults
    subroutine output_sett_initialize_defaults(os, ierr)
        type(output_sett_t), intent(out) :: os
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! Set practical default values for typical calculations
        os%flag_background = 1       ! Compute background data by default
        os%flag_emfield = 1          ! Compute EM field data by default
        os%flag_additional = 0       ! Don't compute additional quantities by default
        os%flag_dispersion = 0       ! Don't compute dispersion by default
        os%num_quants = 0            ! No additional quantities by default
        os%flag_debug = 0            ! No debug by default
        
        ! Clean up arrays - should be allocated as needed
        if (allocated(os%flag_quants)) deallocate(os%flag_quants)
        
    end subroutine output_sett_initialize_defaults
    
    ! =========================================================================
    ! Eigenmode Settings Procedures
    ! =========================================================================
    
    !> @brief Set eigenmode search flag
    subroutine eigmode_sett_set_search_flag(es, flag, ierr)
        type(eigmode_sett_t), intent(inout) :: es
        integer, intent(in) :: flag
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        es%search_flag = flag
        
    end subroutine eigmode_sett_set_search_flag
    
    !> @brief Read eigenmode settings from file
    subroutine eigmode_sett_read_settings(es, path, ierr)
        type(eigmode_sett_t), intent(inout) :: es  ! NOTE: es kept for interface compatibility
        character(len=*), intent(in) :: path
        integer, intent(out) :: ierr
        
        character(len=MAX_PATH_LEN) :: filename
        
        ierr = KILCA_SUCCESS
        
        ! Construct filename
        filename = trim(path) // "eigmode.in"
        
        ! Placeholder for actual file reading
        
    end subroutine eigmode_sett_read_settings
    
    !> @brief Print eigenmode settings to stdout (matches C++ implementation)
    subroutine eigmode_sett_print_settings(es, ierr)
        type(eigmode_sett_t), intent(in) :: es
        integer, intent(out) :: ierr
        
        call eigmode_sett_print_settings_to_unit(es, output_unit, ierr)
        
    end subroutine eigmode_sett_print_settings
    
    !> @brief Print eigenmode settings to specified unit
    subroutine eigmode_sett_print_settings_to_unit(es, unit, ierr)
        type(eigmode_sett_t), intent(in) :: es
        integer, intent(in) :: unit
        integer, intent(out) :: ierr
        integer :: k
        
        ierr = KILCA_SUCCESS
        
        write(unit, '(A)', iostat=ierr) ""
        if (ierr /= 0) return
        write(unit, '(A)', iostat=ierr) "Check for eigmode settings below:"
        if (ierr /= 0) return
        
        write(unit, '(A)', iostat=ierr) ""
        if (ierr /= 0) return
        if (allocated(es%fname)) then
            write(unit, '(A,A)', iostat=ierr) "file name: ", es%fname
            if (ierr /= 0) return
        end if
        write(unit, '(A,I0)', iostat=ierr) "search flag: ", es%search_flag
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "real freq mesh dim: ", es%rdim
        if (ierr /= 0) return
        write(unit, '(A,ES15.8)', iostat=ierr) "real freq mesh minimum: ", es%rfmin
        if (ierr /= 0) return
        write(unit, '(A,ES15.8)', iostat=ierr) "real freq mesh maximum: ", es%rfmax
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "imag freq mesh dim: ", es%idim
        if (ierr /= 0) return
        write(unit, '(A,ES15.8)', iostat=ierr) "imag freq mesh minimum: ", es%ifmin
        if (ierr /= 0) return
        write(unit, '(A,ES15.8)', iostat=ierr) "imag freq mesh maximum: ", es%ifmax
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "stopping criteria: ", es%stop_flag
        if (ierr /= 0) return
        write(unit, '(A,ES15.8)', iostat=ierr) "residual error parameter: ", es%eps_res
        if (ierr /= 0) return
        write(unit, '(A,ES15.8)', iostat=ierr) "abs error parameter: ", es%eps_abs
        if (ierr /= 0) return
        write(unit, '(A,ES15.8)', iostat=ierr) "rel error parameter: ", es%eps_rel
        if (ierr /= 0) return
        write(unit, '(A,ES15.8)', iostat=ierr) "delta for derivative: ", es%delta
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "test roots flag: ", es%test_roots
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "flag_debug: ", es%flag_debug
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "n_zeros: ", es%n_zeros
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "use_winding: ", es%use_winding
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "guess array dimension: ", es%Nguess
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "kmin: ", es%kmin
        if (ierr /= 0) return
        write(unit, '(A,I0)', iostat=ierr) "kmax: ", es%kmax
        if (ierr /= 0) return
        
        write(unit, '(A)', iostat=ierr) "guess array:"
        if (ierr /= 0) return
        if (allocated(es%fstart)) then
            do k = 1, es%Nguess
                write(unit, '(A,I0,A,ES15.8,A,ES15.8,A)', iostat=ierr) "k=", k-1, &
                      "\tf=(", real(es%fstart(k)), ", ", aimag(es%fstart(k)), ")"
                if (ierr /= 0) return
            end do
        end if
        
        write(unit, '(A)', iostat=ierr) ""
        if (ierr /= 0) return
        
        ierr = KILCA_SUCCESS
        
    end subroutine eigmode_sett_print_settings_to_unit
    
    !> @brief Deep copy eigenmode settings
    subroutine eigmode_sett_deep_copy(src, dst, ierr)
        type(eigmode_sett_t), intent(in) :: src
        type(eigmode_sett_t), intent(out) :: dst
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! Copy scalar values
        dst%search_flag = src%search_flag
        dst%rdim = src%rdim
        dst%rfmin = src%rfmin
        dst%rfmax = src%rfmax
        dst%idim = src%idim
        dst%ifmin = src%ifmin
        dst%ifmax = src%ifmax
        dst%stop_flag = src%stop_flag
        dst%eps_res = src%eps_res
        dst%eps_abs = src%eps_abs
        dst%eps_rel = src%eps_rel
        dst%delta = src%delta
        dst%test_roots = src%test_roots
        dst%flag_debug = src%flag_debug
        dst%Nguess = src%Nguess
        dst%kmin = src%kmin
        dst%kmax = src%kmax
        dst%n_zeros = src%n_zeros
        dst%use_winding = src%use_winding
        
        ! Deep copy allocatable strings
        if (allocated(src%fname)) then
            if (allocated(dst%fname)) deallocate(dst%fname)
            dst%fname = src%fname
        else
            if (allocated(dst%fname)) deallocate(dst%fname)
        end if
        
        ! Deep copy allocatable arrays
        if (allocated(src%fstart)) then
            if (allocated(dst%fstart)) deallocate(dst%fstart)
            allocate(dst%fstart(size(src%fstart)))
            dst%fstart = src%fstart
        else
            if (allocated(dst%fstart)) deallocate(dst%fstart)
        end if
        
    end subroutine eigmode_sett_deep_copy
    
    !> @brief Compare two eigenmode settings for equality
    subroutine eigmode_sett_compare(es1, es2, is_equal, ierr)
        type(eigmode_sett_t), intent(in) :: es1, es2
        logical, intent(out) :: is_equal
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        is_equal = .false.
        
        ! Compare scalar values with appropriate tolerance for reals
        if (es1%search_flag /= es2%search_flag) return
        if (es1%rdim /= es2%rdim) return
        if (abs(es1%rfmin - es2%rfmin) > epsilon(es1%rfmin)) return
        if (abs(es1%rfmax - es2%rfmax) > epsilon(es1%rfmax)) return
        if (es1%idim /= es2%idim) return
        if (abs(es1%ifmin - es2%ifmin) > epsilon(es1%ifmin)) return
        if (abs(es1%ifmax - es2%ifmax) > epsilon(es1%ifmax)) return
        if (es1%stop_flag /= es2%stop_flag) return
        if (abs(es1%eps_res - es2%eps_res) > epsilon(es1%eps_res)) return
        if (abs(es1%eps_abs - es2%eps_abs) > epsilon(es1%eps_abs)) return
        if (abs(es1%eps_rel - es2%eps_rel) > epsilon(es1%eps_rel)) return
        if (abs(es1%delta - es2%delta) > epsilon(es1%delta)) return
        if (es1%test_roots /= es2%test_roots) return
        if (es1%flag_debug /= es2%flag_debug) return
        if (es1%Nguess /= es2%Nguess) return
        if (es1%kmin /= es2%kmin) return
        if (es1%kmax /= es2%kmax) return
        if (es1%n_zeros /= es2%n_zeros) return
        if (es1%use_winding /= es2%use_winding) return
        
        ! Compare allocatable strings
        if (allocated(es1%fname) .neqv. allocated(es2%fname)) return
        if (allocated(es1%fname)) then
            if (es1%fname /= es2%fname) return
        end if
        
        ! Compare allocatable complex arrays
        if (allocated(es1%fstart) .neqv. allocated(es2%fstart)) return
        if (allocated(es1%fstart)) then
            if (size(es1%fstart) /= size(es2%fstart)) return
            if (any(abs(es1%fstart - es2%fstart) > epsilon(real(es1%fstart)))) return
        end if
        
        is_equal = .true.
        
    end subroutine eigmode_sett_compare
    
    !> @brief Initialize eigenmode settings with practical defaults
    subroutine eigmode_sett_initialize_defaults(es, ierr)
        type(eigmode_sett_t), intent(out) :: es
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! Set practical default values for eigenmode analysis
        es%search_flag = 0           ! No search by default
        es%rdim = 100                ! 100 points for real frequency mesh
        es%rfmin = 0.0_dp            ! Start from zero frequency
        es%rfmax = 1.0e9_dp          ! Up to 1 GHz
        es%idim = 100                ! 100 points for imaginary frequency mesh
        es%ifmin = -1.0e6_dp         ! -1 MHz imaginary
        es%ifmax = 1.0e6_dp          ! +1 MHz imaginary
        es%stop_flag = 0             ! Default stopping criteria
        es%eps_res = 1.0e-6_dp       ! Residual tolerance
        es%eps_abs = 1.0e-8_dp       ! Absolute tolerance
        es%eps_rel = 1.0e-6_dp       ! Relative tolerance
        es%delta = 1.0e-6_dp         ! Delta for derivative calculation
        es%test_roots = 0            ! Don't test roots by default
        es%flag_debug = 0            ! No debug by default
        es%Nguess = 0                ! No initial guesses by default
        es%kmin = 1                  ! Minimum k value
        es%kmax = 10                 ! Maximum k value
        es%n_zeros = 10              ! Number of zeros to search
        es%use_winding = 0           ! Don't use winding number by default
        
        ! Set default filename
        if (allocated(es%fname)) deallocate(es%fname)
        es%fname = "eigenmode_output.dat"
        
        ! Clean up arrays - should be allocated as needed
        if (allocated(es%fstart)) deallocate(es%fstart)
        
    end subroutine eigmode_sett_initialize_defaults
    
    ! =========================================================================
    ! C Interface Procedures
    ! =========================================================================
    
    !> @brief Set antenna settings from C (mimics C++ interface)
    subroutine set_antenna_settings_c(ant_ptr, ra, wa, I0, flab_re, flab_im, dma, flag_debug) &
            bind(C, name="set_antenna_settings_c_")
        type(c_ptr), intent(inout) :: ant_ptr
        real(c_double), intent(in) :: ra, wa, I0, flab_re, flab_im
        integer(c_int), intent(in) :: dma, flag_debug
        
        type(antenna_t), pointer :: ant
        
        ! Convert C pointer to Fortran pointer
        call c_f_pointer(ant_ptr, ant)
        
        ant%ra = ra
        ant%wa = wa
        ant%I0 = I0
        ant%flab = cmplx(flab_re, flab_im, dp)
        ant%dma = dma
        ant%flag_debug = flag_debug
        
    end subroutine set_antenna_settings_c
    
    !> @brief Set background settings from C (mimics C++ interface)
    subroutine set_background_settings_c(bs_ptr, rtor, rp, B0, flag_back, V_gal_sys, &
                                       V_scale, zele, zion, flag_debug) &
            bind(C, name="set_background_settings_c_")
        type(c_ptr), intent(inout) :: bs_ptr
        real(c_double), intent(in) :: rtor, rp, B0, V_gal_sys, V_scale, zele, zion
        character(c_char), dimension(*), intent(in) :: flag_back
        integer(c_int), intent(in) :: flag_debug
        
        type(back_sett_t), pointer :: bs
        integer :: i
        character(len=:), allocatable :: flag_str
        
        ! Convert C pointer to Fortran pointer
        call c_f_pointer(bs_ptr, bs)
        
        ! Set values
        bs%rtor = rtor
        bs%rp = rp
        bs%B0 = B0
        bs%V_gal_sys = V_gal_sys
        bs%V_scale = V_scale
        bs%zele = zele
        bs%zion = zion
        bs%flag_debug = flag_debug
        
        ! Convert C string to Fortran string
        i = 1
        do while (flag_back(i) /= c_null_char .and. i < 256)
            i = i + 1
        end do
        
        allocate(character(len=i-1) :: flag_str)
        do i = 1, len(flag_str)
            flag_str(i:i) = flag_back(i)
        end do
        bs%flag_back = flag_str
        
    end subroutine set_background_settings_c
    
    !> @brief Set particle settings from C
    subroutine set_particles_settings_c(bs_ptr, mass, charge) &
            bind(C, name="set_particles_settings_c_")
        type(c_ptr), intent(inout) :: bs_ptr
        real(c_double), dimension(2), intent(in) :: mass, charge
        
        type(back_sett_t), pointer :: bs
        
        ! Convert C pointer to Fortran pointer
        call c_f_pointer(bs_ptr, bs)
        
        ! Allocate and set particle arrays
        if (allocated(bs%mass)) deallocate(bs%mass)
        if (allocated(bs%charge)) deallocate(bs%charge)
        
        allocate(bs%mass(2))
        allocate(bs%charge(2))
        
        bs%mass = mass
        bs%charge = charge
        
    end subroutine set_particles_settings_c
    
    !> @brief Set huge factor from C
    subroutine set_huge_factor_c(bs_ptr, fac) &
            bind(C, name="set_huge_factor_c_")
        type(c_ptr), intent(inout) :: bs_ptr
        real(c_double), intent(in) :: fac
        
        type(back_sett_t), pointer :: bs
        
        ! Convert C pointer to Fortran pointer
        call c_f_pointer(bs_ptr, bs)
        
        bs%huge_factor = fac
        
    end subroutine set_huge_factor_c
    
    !> @brief Copy antenna data to module (placeholder)
    subroutine copy_antenna_data_to_antenna_module(ant)
        type(antenna_t), pointer, intent(in) :: ant  ! NOTE: ant kept for interface compatibility
        
        ! In full implementation, would copy to a module-level variable
        ! that's accessible from Fortran routines that need it
        
    end subroutine copy_antenna_data_to_antenna_module
    
    !> @brief Copy background data to module (placeholder)
    subroutine copy_background_data_to_background_module(bs)
        type(back_sett_t), pointer, intent(in) :: bs  ! NOTE: bs kept for interface compatibility
        
        ! In full implementation, would copy to a module-level variable
        
    end subroutine copy_background_data_to_background_module
    
    ! =========================================================================
    ! File Parsing Helper Functions
    ! =========================================================================
    
    !> @brief Read a line and extract double value before '#' comment
    subroutine read_line_get_double(unit, value, ierr)
        integer, intent(in) :: unit
        real(dp), intent(out) :: value
        integer, intent(out) :: ierr
        
        character(len=1024) :: line, value_str
        integer :: comment_pos, iostat
        
        ierr = KILCA_SUCCESS
        
        read(unit, '(a)', iostat=iostat) line
        if (iostat /= 0) then
            ierr = KILCA_ERROR_FILE
            return
        end if
        
        ! Find comment position
        comment_pos = index(line, '#')
        if (comment_pos > 0) then
            value_str = line(1:comment_pos-1)
        else
            value_str = line
        end if
        
        ! Convert to double
        read(value_str, *, iostat=iostat) value
        if (iostat /= 0) then
            ierr = KILCA_ERROR_FORMAT
        end if
        
    end subroutine read_line_get_double
    
    !> @brief Read a line and extract complex value before '#' comment
    subroutine read_line_get_complex(unit, value, ierr)
        integer, intent(in) :: unit
        complex(dp), intent(out) :: value
        integer, intent(out) :: ierr
        
        character(len=1024) :: line, value_str
        integer :: comment_pos, iostat, paren1, paren2, comma
        real(dp) :: re_part, im_part
        
        ierr = KILCA_SUCCESS
        
        read(unit, '(a)', iostat=iostat) line
        if (iostat /= 0) then
            ierr = KILCA_ERROR_FILE
            return
        end if
        
        ! Find comment position
        comment_pos = index(line, '#')
        if (comment_pos > 0) then
            value_str = line(1:comment_pos-1)
        else
            value_str = line
        end if
        
        ! Parse complex number in format (real, imag)
        paren1 = index(value_str, '(')
        paren2 = index(value_str, ')')
        comma = index(value_str, ',')
        
        if (paren1 > 0 .and. paren2 > paren1 .and. comma > paren1 .and. comma < paren2) then
            read(value_str(paren1+1:comma-1), *, iostat=iostat) re_part
            if (iostat /= 0) then
                ierr = KILCA_ERROR_FORMAT
                return
            end if
            
            read(value_str(comma+1:paren2-1), *, iostat=iostat) im_part
            if (iostat /= 0) then
                ierr = KILCA_ERROR_FORMAT
                return
            end if
            
            value = cmplx(re_part, im_part, dp)
        else
            ierr = KILCA_ERROR_FORMAT
        end if
        
    end subroutine read_line_get_complex
    
    !> @brief Read a line and extract integer value before '#' comment
    subroutine read_line_get_int(unit, value, ierr)
        integer, intent(in) :: unit
        integer, intent(out) :: value
        integer, intent(out) :: ierr
        
        character(len=1024) :: line, value_str
        integer :: comment_pos, iostat
        
        ierr = KILCA_SUCCESS
        
        read(unit, '(a)', iostat=iostat) line
        if (iostat /= 0) then
            ierr = KILCA_ERROR_FILE
            return
        end if
        
        ! Find comment position
        comment_pos = index(line, '#')
        if (comment_pos > 0) then
            value_str = line(1:comment_pos-1)
        else
            value_str = line
        end if
        
        ! Convert to integer
        read(value_str, *, iostat=iostat) value
        if (iostat /= 0) then
            ierr = KILCA_ERROR_FORMAT
        end if
        
    end subroutine read_line_get_int
    
    !> @brief Read a line and extract string value before '#' comment
    subroutine read_line_get_string(unit, value, ierr)
        integer, intent(in) :: unit
        character(len=:), allocatable, intent(out) :: value
        integer, intent(out) :: ierr
        
        character(len=1024) :: line, value_str
        integer :: comment_pos, iostat
        
        ierr = KILCA_SUCCESS
        
        read(unit, '(a)', iostat=iostat) line
        if (iostat /= 0) then
            ierr = KILCA_ERROR_FILE
            return
        end if
        
        ! Find comment position
        comment_pos = index(line, '#')
        if (comment_pos > 0) then
            value_str = line(1:comment_pos-1)
        else
            value_str = line
        end if
        
        ! Trim whitespace and allocate
        value_str = trim(adjustl(value_str))
        value = value_str
        
    end subroutine read_line_get_string
    
    !> @brief Skip a line (read and discard)
    subroutine read_line_skip(unit, ierr)
        integer, intent(in) :: unit
        integer, intent(out) :: ierr
        
        character(len=1024) :: line
        integer :: iostat
        
        ierr = KILCA_SUCCESS
        
        read(unit, '(a)', iostat=iostat) line
        if (iostat /= 0) then
            ierr = KILCA_ERROR_FILE
        end if
        
    end subroutine read_line_skip
    
    ! =========================================================================
    ! Full File Parsing Implementations
    ! =========================================================================
    
    !> @brief Read antenna settings from file with full parsing
    subroutine antenna_read_settings_full(ant, path, ierr)
        type(antenna_t), intent(inout) :: ant
        character(len=*), intent(in) :: path
        integer, intent(out) :: ierr
        
        character(len=MAX_PATH_LEN) :: filename, modes_filename
        character(len=256) :: line
        logical :: file_exists
        integer :: unit, iostat, i, paren1, paren2, comma
        
        ierr = KILCA_SUCCESS
        
        ! Construct antenna.in filename
        filename = trim(path) // "antenna.in"
        
        ! Check if file exists
        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            ierr = KILCA_ERROR_FILE
            return
        end if
        
        ! Open and read antenna.in file
        open(newunit=unit, file=filename, status='old', iostat=iostat)
        if (iostat /= 0) then
            ierr = KILCA_ERROR_FILE
            return
        end if
        
        ! Read antenna settings following C++ format
        call read_line_skip(unit, ierr)                    ! Skip comment line
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, ant%ra, ierr)      ! ra
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, ant%wa, ierr)      ! wa
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, ant%I0, ierr)      ! I0
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_complex(unit, ant%flab, ierr)   ! flab
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, ant%dma, ierr)        ! dma
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, ant%flag_debug, ierr) ! flag_debug
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, ant%flag_eigmode, ierr) ! flag_eigmode
        if (ierr /= KILCA_SUCCESS) goto 100
        
        close(unit)
        
        ! Now read modes.in file
        modes_filename = trim(path) // "modes.in"
        
        inquire(file=modes_filename, exist=file_exists)
        if (.not. file_exists) then
            ierr = KILCA_ERROR_FILE
            return
        end if
        
        open(newunit=unit, file=modes_filename, status='old', iostat=iostat)
        if (iostat /= 0) then
            ierr = KILCA_ERROR_FILE
            return
        end if
        
        ! Allocate modes array
        if (allocated(ant%modes)) deallocate(ant%modes)
        allocate(ant%modes(2*ant%dma))
        
        ! Read mode pairs (m,n) in format "(m, n)"
        do i = 1, ant%dma
            read(unit, '(a)', iostat=iostat) line
            if (iostat /= 0) then
                ierr = KILCA_ERROR_FORMAT
                goto 200
            end if
            
            ! Parse line in format "(m, n)"
            paren1 = index(line, '(')
            paren2 = index(line, ')')
            comma = index(line, ',')
            
            if (paren1 > 0 .and. paren2 > paren1 .and. comma > paren1 .and. comma < paren2) then
                read(line(paren1+1:comma-1), *, iostat=iostat) ant%modes(2*i-1)
                if (iostat /= 0) then
                    ierr = KILCA_ERROR_FORMAT
                    goto 200
                end if
                
                read(line(comma+1:paren2-1), *, iostat=iostat) ant%modes(2*i)
                if (iostat /= 0) then
                    ierr = KILCA_ERROR_FORMAT
                    goto 200
                end if
            else
                ierr = KILCA_ERROR_FORMAT
                goto 200
            end if
        end do
        
        close(unit)
        return
        
100     close(unit)
        return
        
200     close(unit)
        if (allocated(ant%modes)) deallocate(ant%modes)
        return
        
    end subroutine antenna_read_settings_full
    
    !> @brief Read background settings from file with full parsing
    subroutine back_sett_read_settings_full(bs, path, ierr)
        type(back_sett_t), intent(inout) :: bs
        character(len=*), intent(in) :: path
        integer, intent(out) :: ierr
        
        character(len=MAX_PATH_LEN) :: filename
        logical :: file_exists
        integer :: unit, iostat
        
        ierr = KILCA_SUCCESS
        
        ! Construct filename
        filename = trim(path) // "background.in"
        
        ! Check if file exists
        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            ierr = KILCA_ERROR_FILE
            return
        end if
        
        ! Open and read file
        open(newunit=unit, file=filename, status='old', iostat=iostat)
        if (iostat /= 0) then
            ierr = KILCA_ERROR_FILE
            return
        end if
        
        ! Read background settings following C++ format
        call read_line_skip(unit, ierr)                    ! Machine settings comment
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, bs%rtor, ierr)     ! rtor
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, bs%rp, ierr)       ! rp  
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, bs%B0, ierr)       ! B0
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_skip(unit, ierr)                    ! Background settings comment
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_string(unit, bs%path2profiles, ierr) ! path2profiles
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, bs%calc_back, ierr)   ! calc_back
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_string(unit, bs%flag_back, ierr) ! flag_back
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, bs%N, ierr)           ! N
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, bs%V_gal_sys, ierr) ! V_gal_sys
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, bs%V_scale, ierr)  ! V_scale
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, bs%m_i, ierr)      ! m_i
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, bs%zele, ierr)     ! zele
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, bs%zion, ierr)     ! zion
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_skip(unit, ierr)                    ! Checking settings comment
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, bs%flag_debug, ierr)  ! flag_debug
        if (ierr /= KILCA_SUCCESS) goto 100
        
        close(unit)
        return
        
100     close(unit)
        return
        
    end subroutine back_sett_read_settings_full
    
    !> @brief Read output settings from file with full parsing
    subroutine output_sett_read_settings_full(os, path, ierr)
        type(output_sett_t), intent(inout) :: os
        character(len=*), intent(in) :: path
        integer, intent(out) :: ierr
        
        character(len=MAX_PATH_LEN) :: filename
        logical :: file_exists
        integer :: unit, iostat, i
        
        ierr = KILCA_SUCCESS
        
        ! Construct filename
        filename = trim(path) // "output.in"
        
        ! Check if file exists
        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            ierr = KILCA_ERROR_FILE
            return
        end if
        
        ! Open and read file
        open(newunit=unit, file=filename, status='old', iostat=iostat)
        if (iostat /= 0) then
            ierr = KILCA_ERROR_FILE
            return
        end if
        
        ! Read output settings
        call read_line_skip(unit, ierr)                        ! Comment
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, os%flag_background, ierr) ! flag_background
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, os%flag_emfield, ierr)    ! flag_emfield
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, os%flag_additional, ierr) ! flag_additional
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, os%flag_dispersion, ierr) ! flag_dispersion
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, os%num_quants, ierr)      ! num_quants
        if (ierr /= KILCA_SUCCESS) goto 100
        
        ! Allocate and read flag_quants array
        if (os%num_quants > 0) then
            if (allocated(os%flag_quants)) deallocate(os%flag_quants)
            allocate(os%flag_quants(os%num_quants))
            
            do i = 1, os%num_quants
                call read_line_get_int(unit, os%flag_quants(i), ierr)
                if (ierr /= KILCA_SUCCESS) goto 100
            end do
        end if
        
        call read_line_get_int(unit, os%flag_debug, ierr)      ! flag_debug
        if (ierr /= KILCA_SUCCESS) goto 100
        
        close(unit)
        return
        
100     close(unit)
        if (allocated(os%flag_quants)) deallocate(os%flag_quants)
        return
        
    end subroutine output_sett_read_settings_full
    
    !> @brief Read eigenmode settings from file with full parsing
    subroutine eigmode_sett_read_settings_full(es, path, ierr)
        type(eigmode_sett_t), intent(inout) :: es
        character(len=*), intent(in) :: path
        integer, intent(out) :: ierr
        
        character(len=MAX_PATH_LEN) :: filename
        logical :: file_exists
        integer :: unit, iostat, i
        
        ierr = KILCA_SUCCESS
        
        ! Construct filename
        filename = trim(path) // "eigmode.in"
        
        ! Check if file exists
        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            ierr = KILCA_ERROR_FILE
            return
        end if
        
        ! Open and read file
        open(newunit=unit, file=filename, status='old', iostat=iostat)
        if (iostat /= 0) then
            ierr = KILCA_ERROR_FILE
            return
        end if
        
        ! Read eigenmode settings
        call read_line_skip(unit, ierr)                      ! Comment
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_string(unit, es%fname, ierr)      ! fname
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, es%search_flag, ierr)   ! search_flag
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, es%rdim, ierr)          ! rdim
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, es%rfmin, ierr)      ! rfmin
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, es%rfmax, ierr)      ! rfmax
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, es%idim, ierr)          ! idim
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, es%ifmin, ierr)      ! ifmin
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, es%ifmax, ierr)      ! ifmax
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, es%stop_flag, ierr)     ! stop_flag
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, es%eps_res, ierr)    ! eps_res
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, es%eps_abs, ierr)    ! eps_abs
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, es%eps_rel, ierr)    ! eps_rel
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_double(unit, es%delta, ierr)      ! delta
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, es%test_roots, ierr)    ! test_roots
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, es%flag_debug, ierr)    ! flag_debug
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, es%Nguess, ierr)        ! Nguess
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, es%kmin, ierr)          ! kmin
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, es%kmax, ierr)          ! kmax
        if (ierr /= KILCA_SUCCESS) goto 100
        
        ! Allocate and read fstart array
        if (es%Nguess > 0) then
            if (allocated(es%fstart)) deallocate(es%fstart)
            allocate(es%fstart(es%Nguess))
            
            do i = 1, es%Nguess
                call read_line_get_complex(unit, es%fstart(i), ierr)
                if (ierr /= KILCA_SUCCESS) goto 100
            end do
        end if
        
        call read_line_get_int(unit, es%n_zeros, ierr)       ! n_zeros
        if (ierr /= KILCA_SUCCESS) goto 100
        
        call read_line_get_int(unit, es%use_winding, ierr)   ! use_winding
        if (ierr /= KILCA_SUCCESS) goto 100
        
        close(unit)
        return
        
100     close(unit)
        if (allocated(es%fstart)) deallocate(es%fstart)
        return
        
    end subroutine eigmode_sett_read_settings_full
    
    !> @brief Read all settings from files with full parsing
    subroutine settings_read_all_full(sd, ierr)
        type(settings_t), intent(inout) :: sd
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        print *, ">>>>> Reading settings from ", trim(sd%path2project)
        
        ! Read antenna settings
        call antenna_read_settings_full(sd%antenna_settings, sd%path2project, ierr)
        if (ierr /= KILCA_SUCCESS) return
        call copy_antenna_data_to_antenna_module(sd%as)
        
        ! Read background settings
        call back_sett_read_settings_full(sd%background_settings, sd%path2project, ierr)
        if (ierr /= KILCA_SUCCESS) return
        call copy_background_data_to_background_module(sd%bs)
        
        ! Read output settings
        call output_sett_read_settings_full(sd%output_settings, sd%path2project, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Read eigenmode settings
        call eigmode_sett_read_settings_full(sd%eigmode_settings, sd%path2project, ierr)
        
    end subroutine settings_read_all_full
    
    ! =========================================================================
    ! Individual Settings Validation Procedures
    ! =========================================================================
    
    !> @brief Validate antenna settings
    subroutine antenna_validate(ant, is_valid, error_msg, ierr)
        type(antenna_t), intent(in) :: ant
        logical, intent(out) :: is_valid
        character(len=*), intent(out) :: error_msg
        integer, intent(out) :: ierr
        integer :: i
        
        ierr = KILCA_SUCCESS
        is_valid = .true.
        error_msg = ""
        
        ! Validate antenna radius
        if (ant%ra <= 0.0_dp) then
            is_valid = .false.
            error_msg = "Antenna radius (ra) must be positive, got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') ant%ra
            return
        end if
        
        ! Validate current layer width
        if (ant%wa <= 0.0_dp) then
            is_valid = .false.
            error_msg = "Antenna width (wa) must be positive, got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') ant%wa
            return
        end if
        
        ! Validate antenna current
        if (ant%I0 < 0.0_dp) then
            is_valid = .false.
            error_msg = "Antenna current (I0) cannot be negative, got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') ant%I0
            return
        end if
        
        ! Validate dimension of modes array
        if (ant%dma <= 0) then
            is_valid = .false.
            error_msg = "Modes array dimension (dma) must be positive, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') ant%dma
            return
        end if
        
        ! Validate modes array allocation
        if (.not. allocated(ant%modes)) then
            is_valid = .false.
            error_msg = "Modes array is not allocated"
            return
        end if
        
        ! Validate modes array size
        if (size(ant%modes) /= 2 * ant%dma) then
            is_valid = .false.
            error_msg = "Modes array size mismatch: expected "
            write(error_msg(len_trim(error_msg)+1:), '(i0,a,i0)') 2*ant%dma, ", got ", size(ant%modes)
            return
        end if
        
        ! Validate frequency (check for NaN/Inf)
        if (.not. (real(ant%flab) == real(ant%flab))) then  ! NaN check
            is_valid = .false.
            error_msg = "Antenna frequency real part is NaN"
            return
        end if
        
        if (.not. (aimag(ant%flab) == aimag(ant%flab))) then  ! NaN check
            is_valid = .false.
            error_msg = "Antenna frequency imaginary part is NaN"
            return
        end if
        
        ! Validate flag values
        if (ant%flag_debug < 0 .or. ant%flag_debug > 2) then
            is_valid = .false.
            error_msg = "Debug flag must be 0, 1, or 2, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') ant%flag_debug
            return
        end if
        
        if (ant%flag_eigmode < 0 .or. ant%flag_eigmode > 1) then
            is_valid = .false.
            error_msg = "Eigenmode flag must be 0 or 1, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') ant%flag_eigmode
            return
        end if
        
        ! Enhanced physics-based validation
        ! Validate reasonable antenna radius range (1 cm to 1000 cm)
        if (ant%ra < 1.0_dp .or. ant%ra > 1000.0_dp) then
            is_valid = .false.
            error_msg = "Antenna radius (ra) should be between 1 and 1000 cm for typical tokamaks, got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') ant%ra
            return
        end if
        
        ! Validate reasonable current layer width (must be smaller than radius)
        if (ant%wa >= ant%ra) then
            is_valid = .false.
            error_msg = "Antenna width (wa) must be smaller than radius (ra): wa="
            write(error_msg(len_trim(error_msg)+1:), '(g0,a,g0)') ant%wa, " >= ra=", ant%ra
            return
        end if
        
        ! Validate reasonable current amplitude (0 to 1e15 statamps)
        if (ant%I0 > 1.0e15_dp) then
            is_valid = .false.
            error_msg = "Antenna current (I0) seems unreasonably large for typical applications, got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') ant%I0
            return
        end if
        
        ! Validate reasonable frequency range (1 MHz to 1000 MHz for typical ICRF)
        if (real(ant%flab) > 0.0_dp) then  ! Only check if frequency is positive
            if (real(ant%flab) < 1.0e6_dp .or. real(ant%flab) > 1.0e9_dp) then
                is_valid = .false.
                error_msg = "Antenna frequency should typically be between 1 MHz and 1 GHz for ICRF, got: "
                write(error_msg(len_trim(error_msg)+1:), '(g0,a)') real(ant%flab), " Hz"
                return
            end if
        end if
        
        ! Validate modes array values (m and n should be reasonable integers)
        if (allocated(ant%modes)) then
            do i = 1, ant%dma
                ! m values (poloidal mode numbers) typically -20 to +20
                if (abs(ant%modes(2*i-1)) > 20) then
                    is_valid = .false.
                    error_msg = "Poloidal mode number (m) seems unreasonably large, got: "
                    write(error_msg(len_trim(error_msg)+1:), '(i0)') ant%modes(2*i-1)
                    return
                end if
                
                ! n values (toroidal mode numbers) typically -100 to +100
                if (abs(ant%modes(2*i)) > 100) then
                    is_valid = .false.
                    error_msg = "Toroidal mode number (n) seems unreasonably large, got: "
                    write(error_msg(len_trim(error_msg)+1:), '(i0)') ant%modes(2*i)
                    return
                end if
            end do
        end if
        
    end subroutine antenna_validate
    
    !> @brief Validate background settings
    subroutine back_sett_validate(bs, is_valid, error_msg, ierr)
        type(back_sett_t), intent(in) :: bs
        logical, intent(out) :: is_valid
        character(len=*), intent(out) :: error_msg
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        is_valid = .true.
        error_msg = ""
        
        ! Validate torus radius
        if (bs%rtor <= 0.0_dp) then
            is_valid = .false.
            error_msg = "Torus radius (rtor) must be positive, got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') bs%rtor
            return
        end if
        
        ! Physics-based range validation for torus radius (typical tokamak range)
        if (bs%rtor < 50.0_dp .or. bs%rtor > 1000.0_dp) then
            is_valid = .false.
            error_msg = "Torus radius (rtor) outside typical tokamak range [50-1000 cm], got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') bs%rtor
            return
        end if
        
        ! Validate plasma radius
        if (bs%rp <= 0.0_dp) then
            is_valid = .false.
            error_msg = "Plasma radius (rp) must be positive, got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') bs%rp
            return
        end if
        
        ! Physics-based range validation for plasma radius
        if (bs%rp > 200.0_dp) then
            is_valid = .false.
            error_msg = "Plasma radius (rp) unreasonably large (>200 cm), got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') bs%rp
            return
        end if
        
        ! Physical consistency: plasma radius must be less than torus radius
        if (bs%rp >= bs%rtor) then
            is_valid = .false.
            error_msg = "Plasma radius (rp) must be less than torus radius (rtor), got rp="
            write(error_msg(len_trim(error_msg)+1:), '(g0)') bs%rp
            error_msg = trim(error_msg) // ", rtor="
            write(error_msg(len_trim(error_msg)+1:), '(g0)') bs%rtor
            return
        end if
        
        ! Validate magnetic field
        if (bs%B0 <= 0.0_dp) then
            is_valid = .false.
            error_msg = "Magnetic field (B0) must be positive, got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') bs%B0
            return
        end if
        
        ! Physics-based range validation for magnetic field (typical tokamak range)
        if (bs%B0 < 1000.0_dp .or. bs%B0 > 100000.0_dp) then
            is_valid = .false.
            error_msg = "Magnetic field (B0) outside typical tokamak range [0.1-10 Tesla], got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') bs%B0
            return
        end if
        
        ! Validate calc_back flag
        if (bs%calc_back == 0) then
            is_valid = .false.
            error_msg = "Background calculation flag (calc_back) cannot be zero"
            return
        end if
        
        ! Validate spline degree (must be odd)
        if (bs%N <= 0) then
            is_valid = .false.
            error_msg = "Spline degree (N) must be positive, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') bs%N
            return
        end if
        
        if (mod(bs%N, 2) == 0) then
            is_valid = .false.
            error_msg = "Spline degree (N) must be odd, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') bs%N
            return
        end if
        
        ! Validate ion mass
        if (bs%m_i <= 0.0_dp) then
            is_valid = .false.
            error_msg = "Ion mass (m_i) must be positive, got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') bs%m_i
            return
        end if
        
        ! Physics-based range validation for ion mass (realistic ion range)
        if (bs%m_i < 0.5_dp .or. bs%m_i > 50.0_dp) then
            is_valid = .false.
            error_msg = "Ion mass (m_i) outside realistic range [0.5-50 amu], got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') bs%m_i
            return
        end if
        
        ! Validate collision coefficients
        if (bs%zele < 0.0_dp) then
            is_valid = .false.
            error_msg = "Electron collision coefficient (zele) cannot be negative, got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') bs%zele
            return
        end if
        
        if (bs%zion < 0.0_dp) then
            is_valid = .false.
            error_msg = "Ion collision coefficient (zion) cannot be negative, got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') bs%zion
            return
        end if
        
        ! Physics-based validation for collision coefficients (typical plasma range)
        if (bs%zele > 100.0_dp) then
            is_valid = .false.
            error_msg = "Electron collision coefficient (zele) unreasonably large (>100), got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') bs%zele
            return
        end if
        
        if (bs%zion > 100.0_dp) then
            is_valid = .false.
            error_msg = "Ion collision coefficient (zion) unreasonably large (>100), got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') bs%zion
            return
        end if
        
        ! Validate velocity scale
        if (bs%V_scale == 0.0_dp) then
            is_valid = .false.
            error_msg = "Velocity scale (V_scale) cannot be zero"
            return
        end if
        
        ! Validate galaxy system voltage (physics-based range)
        if (abs(bs%V_gal_sys) > 1.0e12_dp) then
            is_valid = .false.
            error_msg = "Galaxy system voltage (V_gal_sys) unreasonably large (>1e12), got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') bs%V_gal_sys
            return
        end if
        
        ! Validate huge factor (should be positive and reasonable)
        if (bs%huge_factor <= 0.0_dp .or. bs%huge_factor > 1.0e50_dp) then
            is_valid = .false.
            error_msg = "Huge factor outside reasonable range (0, 1e50], got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') bs%huge_factor
            return
        end if
        
        ! Validate file path (must be allocated and non-empty if used)
        if (allocated(bs%path2profiles)) then
            if (len_trim(bs%path2profiles) == 0) then
                is_valid = .false.
                error_msg = "Profile path (path2profiles) cannot be empty string"
                return
            end if
        end if
        
        ! Validate flag values
        if (bs%flag_debug < 0 .or. bs%flag_debug > 1) then
            is_valid = .false.
            error_msg = "Debug flag must be 0 or 1, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') bs%flag_debug
            return
        end if
        
        ! Validate allocated arrays if present
        if (allocated(bs%mass)) then
            if (size(bs%mass) /= 2) then
                is_valid = .false.
                error_msg = "Mass array must have size 2, got: "
                write(error_msg(len_trim(error_msg)+1:), '(i0)') size(bs%mass)
                return
            end if
            
            if (any(bs%mass <= 0.0_dp)) then
                is_valid = .false.
                error_msg = "All masses must be positive"
                return
            end if
        end if
        
        if (allocated(bs%charge)) then
            if (size(bs%charge) /= 2) then
                is_valid = .false.
                error_msg = "Charge array must have size 2, got: "
                write(error_msg(len_trim(error_msg)+1:), '(i0)') size(bs%charge)
                return
            end if
        end if
        
    end subroutine back_sett_validate
    
    !> @brief Compute derived values for background settings (mass and charge arrays)
    subroutine background_settings_compute_derived(bs, ierr)
        use kilca_types_m, only: m_proton, m_electron, e_charge
        type(back_sett_t), intent(inout) :: bs
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! Allocate arrays if not allocated or wrong size
        if (allocated(bs%mass)) then
            if (size(bs%mass) /= 2) then
                deallocate(bs%mass)
                allocate(bs%mass(2))
            end if
        else
            allocate(bs%mass(2))
        end if
        
        if (allocated(bs%charge)) then
            if (size(bs%charge) /= 2) then
                deallocate(bs%charge)
                allocate(bs%charge(2))
            end if
        else
            allocate(bs%charge(2))
        end if
        
        ! Compute masses: ions and electrons
        bs%mass(1) = bs%m_i * m_proton    ! Ion mass = m_i * proton mass
        bs%mass(2) = m_electron           ! Electron mass
        
        ! Compute charges: ions and electrons
        bs%charge(1) = e_charge           ! Ion charge = +e
        bs%charge(2) = -e_charge          ! Electron charge = -e
        
        ! Compute additional derived arrays for plasma physics
        call compute_profile_arrays(bs, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        call compute_thermal_parameters(bs, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
    end subroutine background_settings_compute_derived
    
    !> @brief Compute profile arrays from background parameters
    subroutine compute_profile_arrays(bs, ierr)
        type(back_sett_t), intent(inout) :: bs
        integer, intent(out) :: ierr
        
        integer, parameter :: n_profile_points = 100
        integer :: i
        real(dp) :: rho, rho_norm
        
        ierr = KILCA_SUCCESS
        
        ! Allocate profile arrays
        if (allocated(bs%profile_grid)) deallocate(bs%profile_grid)
        if (allocated(bs%density)) deallocate(bs%density)
        if (allocated(bs%temperature)) deallocate(bs%temperature)
        if (allocated(bs%velocity)) deallocate(bs%velocity)
        if (allocated(bs%q_profile)) deallocate(bs%q_profile)
        if (allocated(bs%pressure)) deallocate(bs%pressure)
        
        allocate(bs%profile_grid(n_profile_points))
        allocate(bs%density(n_profile_points))
        allocate(bs%temperature(n_profile_points))
        allocate(bs%velocity(n_profile_points))
        allocate(bs%q_profile(n_profile_points))
        allocate(bs%pressure(n_profile_points))
        
        ! Create simple radial grid (normalized radius)
        do i = 1, n_profile_points
            rho_norm = real(i-1, dp) / real(n_profile_points-1, dp)
            bs%profile_grid(i) = rho_norm * bs%rp  ! Convert to physical radius
            rho = rho_norm
            
            ! Simple parabolic density profile: n(rho) = n0 * (1 - rho^2)
            bs%density(i) = 1.0e19_dp * (1.0_dp - rho**2)  ! particles/m³
            
            ! Simple parabolic temperature profile: T(rho) = T0 * (1 - rho^2) + T_edge
            bs%temperature(i) = 2000.0_dp * (1.0_dp - rho**2) + 100.0_dp  ! eV
            
            ! Simple velocity profile (mostly zero with small rotation)
            bs%velocity(i) = bs%V_gal_sys * rho * 1.0e-2_dp  ! cm/s
            
            ! Safety factor profile: q(rho) = q0 * (1 + rho^2)
            bs%q_profile(i) = 1.0_dp * (1.0_dp + 2.0_dp * rho**2)
            
            ! Pressure = n * T (in SI: Pa = m^-3 * eV * e_charge)
            bs%pressure(i) = bs%density(i) * bs%temperature(i) * 1.602e-19_dp  ! Pa
        end do
        
    end subroutine compute_profile_arrays
    
    !> @brief Compute thermal parameters from profiles
    subroutine compute_thermal_parameters(bs, ierr)
        use kilca_types_m, only: k_boltz, m_proton, m_electron
        type(back_sett_t), intent(inout) :: bs
        integer, intent(out) :: ierr
        
        integer :: i, n_points
        real(dp) :: T_eV, m_ion, m_elec
        real(dp) :: nu_ion_ion, nu_elec_elec
        
        ierr = KILCA_SUCCESS
        
        if (.not. allocated(bs%density) .or. .not. allocated(bs%temperature)) then
            ierr = KILCA_ERROR_INVALID_INPUT
            return
        end if
        
        n_points = size(bs%density)
        
        ! Allocate thermal parameter arrays
        if (allocated(bs%v_thermal)) deallocate(bs%v_thermal)
        if (allocated(bs%nu_collision)) deallocate(bs%nu_collision)
        
        allocate(bs%v_thermal(n_points))
        allocate(bs%nu_collision(n_points))
        
        m_ion = bs%m_i * m_proton
        m_elec = m_electron
        
        do i = 1, n_points
            T_eV = bs%temperature(i)
            
            ! Thermal velocity: v_th = sqrt(2*k*T/m)
            ! For electrons (lighter, faster)
            bs%v_thermal(i) = sqrt(2.0_dp * k_boltz * T_eV * 1.602e-19_dp / m_elec)  ! m/s
            
            ! Collision frequency estimation (simplified)
            ! nu ~ n * sigma * v_rel, where sigma ~ e^4/(4*pi*eps0*T)^2
            nu_elec_elec = bs%zele * 2.9e-12_dp * bs%density(i) * 1.0e-6_dp / (T_eV**1.5_dp)  ! s^-1
            nu_ion_ion = bs%zion * 4.8e-14_dp * bs%density(i) * 1.0e-6_dp / (T_eV**1.5_dp)    ! s^-1
            
            ! Use electron collision frequency (typically higher)
            bs%nu_collision(i) = nu_elec_elec
        end do
        
    end subroutine compute_thermal_parameters
    
    !> @brief Validate output settings
    subroutine output_sett_validate(os, is_valid, error_msg, ierr)
        type(output_sett_t), intent(in) :: os
        logical, intent(out) :: is_valid
        character(len=*), intent(out) :: error_msg
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        is_valid = .true.
        error_msg = ""
        
        ! Validate flag values (should be 0, 1, or 2)
        if (os%flag_background < 0 .or. os%flag_background > 2) then
            is_valid = .false.
            error_msg = "Background flag must be 0, 1, or 2, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') os%flag_background
            return
        end if
        
        if (os%flag_emfield < 0 .or. os%flag_emfield > 2) then
            is_valid = .false.
            error_msg = "EM field flag must be 0, 1, or 2, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') os%flag_emfield
            return
        end if
        
        if (os%flag_additional < 0 .or. os%flag_additional > 2) then
            is_valid = .false.
            error_msg = "Additional flag must be 0, 1, or 2, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') os%flag_additional
            return
        end if
        
        if (os%flag_dispersion < 0 .or. os%flag_dispersion > 2) then
            is_valid = .false.
            error_msg = "Dispersion flag must be 0, 1, or 2, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') os%flag_dispersion
            return
        end if
        
        ! Validate number of quantities
        if (os%num_quants < 0) then
            is_valid = .false.
            error_msg = "Number of quantities (num_quants) cannot be negative, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') os%num_quants
            return
        end if
        
        ! Validate flag_quants array if num_quants > 0
        if (os%num_quants > 0) then
            if (.not. allocated(os%flag_quants)) then
                is_valid = .false.
                error_msg = "Quantities flags array not allocated but num_quants > 0"
                return
            end if
            
            if (size(os%flag_quants) /= os%num_quants) then
                is_valid = .false.
                error_msg = "Quantities flags array size mismatch: expected "
                write(error_msg(len_trim(error_msg)+1:), '(i0,a,i0)') os%num_quants, ", got ", size(os%flag_quants)
                return
            end if
            
            ! Validate individual quantity flags
            if (any(os%flag_quants < 0) .or. any(os%flag_quants > 2)) then
                is_valid = .false.
                error_msg = "All quantity flags must be 0, 1, or 2"
                return
            end if
        end if
        
        ! Validate debug flag
        if (os%flag_debug < 0 .or. os%flag_debug > 1) then
            is_valid = .false.
            error_msg = "Debug flag must be 0 or 1, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') os%flag_debug
            return
        end if
        
    end subroutine output_sett_validate
    
    !> @brief Set flag_quants array with validation
    subroutine output_settings_set_flag_quants(os, flag_quants, ierr)
        type(output_sett_t), intent(inout) :: os
        integer, dimension(:), intent(in) :: flag_quants
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! Clean up existing array
        if (allocated(os%flag_quants)) deallocate(os%flag_quants)
        
        ! Allocate and copy
        allocate(os%flag_quants(size(flag_quants)))
        os%flag_quants = flag_quants
        
        ! Update num_quants to be consistent
        os%num_quants = size(flag_quants)
        
        ! Validate array values (should be 0, 1, or 2)
        if (any(flag_quants < 0) .or. any(flag_quants > 2)) then
            ierr = KILCA_ERROR_INVALID_INPUT
            return
        end if
        
    end subroutine output_settings_set_flag_quants
    
    !> @brief Get flag_quants array
    subroutine output_settings_get_flag_quants(os, flag_quants, ierr)
        type(output_sett_t), intent(in) :: os
        integer, dimension(:), allocatable, intent(out) :: flag_quants
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        if (allocated(flag_quants)) deallocate(flag_quants)
        
        if (allocated(os%flag_quants)) then
            allocate(flag_quants(size(os%flag_quants)))
            flag_quants = os%flag_quants
        else
            ! Return empty array if not allocated
            allocate(flag_quants(0))
        end if
        
    end subroutine output_settings_get_flag_quants
    
    !> @brief Validate eigenmode settings
    subroutine eigmode_sett_validate(es, is_valid, error_msg, ierr)
        type(eigmode_sett_t), intent(in) :: es
        logical, intent(out) :: is_valid
        character(len=*), intent(out) :: error_msg
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        is_valid = .true.
        error_msg = ""
        
        ! Validate search flag
        if (es%search_flag < 0) then
            is_valid = .false.
            error_msg = "Search flag cannot be negative, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') es%search_flag
            return
        end if
        
        ! Validate grid dimensions
        if (es%rdim <= 0) then
            is_valid = .false.
            error_msg = "Real grid dimension (rdim) must be positive, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') es%rdim
            return
        end if
        
        if (es%idim <= 0) then
            is_valid = .false.
            error_msg = "Imaginary grid dimension (idim) must be positive, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') es%idim
            return
        end if
        
        ! Validate frequency ranges
        if (es%rfmin >= es%rfmax) then
            is_valid = .false.
            error_msg = "Real frequency range invalid: rfmin >= rfmax ("
            write(error_msg(len_trim(error_msg)+1:), '(g0,a,g0,a)') es%rfmin, " >= ", es%rfmax, ")"
            return
        end if
        
        if (es%ifmin >= es%ifmax) then
            is_valid = .false.
            error_msg = "Imaginary frequency range invalid: ifmin >= ifmax ("
            write(error_msg(len_trim(error_msg)+1:), '(g0,a,g0,a)') es%ifmin, " >= ", es%ifmax, ")"
            return
        end if
        
        ! Validate tolerances
        if (es%eps_res <= 0.0_dp) then
            is_valid = .false.
            error_msg = "Residual tolerance (eps_res) must be positive, got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') es%eps_res
            return
        end if
        
        if (es%eps_abs <= 0.0_dp) then
            is_valid = .false.
            error_msg = "Absolute tolerance (eps_abs) must be positive, got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') es%eps_abs
            return
        end if
        
        if (es%eps_rel <= 0.0_dp) then
            is_valid = .false.
            error_msg = "Relative tolerance (eps_rel) must be positive, got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') es%eps_rel
            return
        end if
        
        if (es%delta <= 0.0_dp) then
            is_valid = .false.
            error_msg = "Delta for derivatives must be positive, got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') es%delta
            return
        end if
        
        ! Validate k values
        if (es%kmin < 1) then
            is_valid = .false.
            error_msg = "Minimum k value must be >= 1, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') es%kmin
            return
        end if
        
        if (es%kmax < es%kmin) then
            is_valid = .false.
            error_msg = "Maximum k value must be >= kmin, got kmax="
            write(error_msg(len_trim(error_msg)+1:), '(i0,a,i0)') es%kmax, ", kmin=", es%kmin
            return
        end if
        
        ! Validate Nguess and fstart array
        if (es%Nguess < 0) then
            is_valid = .false.
            error_msg = "Number of guess values (Nguess) cannot be negative, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') es%Nguess
            return
        end if
        
        if (es%Nguess > 0) then
            if (.not. allocated(es%fstart)) then
                is_valid = .false.
                error_msg = "Start frequencies array not allocated but Nguess > 0"
                return
            end if
            
            if (size(es%fstart) /= es%Nguess) then
                is_valid = .false.
                error_msg = "Start frequencies array size mismatch: expected "
                write(error_msg(len_trim(error_msg)+1:), '(i0,a,i0)') es%Nguess, ", got ", size(es%fstart)
                return
            end if
        end if
        
        ! Validate n_zeros
        if (es%n_zeros <= 0) then
            is_valid = .false.
            error_msg = "Number of zeros to find must be positive, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') es%n_zeros
            return
        end if
        
        ! Validate flag values
        if (es%test_roots < 0 .or. es%test_roots > 1) then
            is_valid = .false.
            error_msg = "Test roots flag must be 0 or 1, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') es%test_roots
            return
        end if
        
        if (es%flag_debug < 0 .or. es%flag_debug > 1) then
            is_valid = .false.
            error_msg = "Debug flag must be 0 or 1, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') es%flag_debug
            return
        end if
        
        if (es%use_winding < 0 .or. es%use_winding > 1) then
            is_valid = .false.
            error_msg = "Use winding flag must be 0 or 1, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') es%use_winding
            return
        end if
        
    end subroutine eigmode_sett_validate
    
    !> @brief Validate complete settings structure
    subroutine settings_validate_complete(sd, is_valid, error_msg, ierr)
        type(settings_t), intent(in) :: sd
        logical, intent(out) :: is_valid
        character(len=*), intent(out) :: error_msg
        integer, intent(out) :: ierr
        
        character(len=1024) :: local_error
        
        ierr = KILCA_SUCCESS
        is_valid = .true.
        error_msg = ""
        
        ! Validate path existence
        if (.not. allocated(sd%path2project)) then
            is_valid = .false.
            error_msg = "Project path not allocated"
            return
        end if
        
        if (len_trim(sd%path2project) == 0) then
            is_valid = .false.
            error_msg = "Project path is empty"
            return
        end if
        
        ! Validate antenna settings
        call antenna_validate(sd%antenna_settings, is_valid, local_error, ierr)
        if (ierr /= KILCA_SUCCESS) return
        if (.not. is_valid) then
            error_msg = "Antenna validation failed: " // trim(local_error)
            return
        end if
        
        ! Validate background settings
        call back_sett_validate(sd%background_settings, is_valid, local_error, ierr)
        if (ierr /= KILCA_SUCCESS) return
        if (.not. is_valid) then
            error_msg = "Background validation failed: " // trim(local_error)
            return
        end if
        
        ! Validate output settings
        call output_sett_validate(sd%output_settings, is_valid, local_error, ierr)
        if (ierr /= KILCA_SUCCESS) return
        if (.not. is_valid) then
            error_msg = "Output validation failed: " // trim(local_error)
            return
        end if
        
        ! Validate eigenmode settings
        call eigmode_sett_validate(sd%eigmode_settings, is_valid, local_error, ierr)
        if (ierr /= KILCA_SUCCESS) return
        if (.not. is_valid) then
            error_msg = "Eigenmode validation failed: " // trim(local_error)
            return
        end if
        
        ! Check pointer consistency (pointers should be associated)
        if (.not. associated(sd%as)) then
            is_valid = .false.
            error_msg = "Antenna settings pointer not associated"
            return
        end if
        
        if (.not. associated(sd%bs)) then
            is_valid = .false.
            error_msg = "Background settings pointer not associated"
            return
        end if
        
        if (.not. associated(sd%os)) then
            is_valid = .false.
            error_msg = "Output settings pointer not associated"
            return
        end if
        
        if (.not. associated(sd%es)) then
            is_valid = .false.
            error_msg = "Eigenmode settings pointer not associated"
            return
        end if
        
    end subroutine settings_validate_complete
    
    !> @brief Validate consistency between different settings
    subroutine settings_validate_consistency(sd, is_valid, error_msg, ierr)
        type(settings_t), intent(in) :: sd
        logical, intent(out) :: is_valid
        character(len=*), intent(out) :: error_msg
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        is_valid = .true.
        error_msg = ""
        
        ! Check antenna-output consistency
        if (sd%antenna_settings%dma > 0 .and. allocated(sd%antenna_settings%modes)) then
            if (size(sd%antenna_settings%modes) /= 2 * sd%antenna_settings%dma) then
                is_valid = .false.
                error_msg = "Antenna modes array size inconsistent with dma"
                return
            end if
        end if
        
        ! Check output quantities consistency
        if (sd%output_settings%num_quants > 0) then
            if (.not. allocated(sd%output_settings%flag_quants)) then
                is_valid = .false.
                error_msg = "Output quantities flags not allocated but num_quants > 0"
                return
            end if
            
            if (size(sd%output_settings%flag_quants) /= sd%output_settings%num_quants) then
                is_valid = .false.
                error_msg = "Output quantities array size inconsistent with num_quants"
                return
            end if
        end if
        
        ! Check eigenmode start frequencies consistency
        if (sd%eigmode_settings%Nguess > 0) then
            if (.not. allocated(sd%eigmode_settings%fstart)) then
                is_valid = .false.
                error_msg = "Eigenmode start frequencies not allocated but Nguess > 0"
                return
            end if
            
            if (size(sd%eigmode_settings%fstart) /= sd%eigmode_settings%Nguess) then
                is_valid = .false.
                error_msg = "Eigenmode start frequencies array size inconsistent with Nguess"
                return
            end if
        end if
        
        ! Check background particle arrays consistency
        if (allocated(sd%background_settings%mass) .and. allocated(sd%background_settings%charge)) then
            if (size(sd%background_settings%mass) /= size(sd%background_settings%charge)) then
                is_valid = .false.
                error_msg = "Background mass and charge arrays have different sizes"
                return
            end if
        end if
        
        ! Physical consistency checks
        if (sd%background_settings%rp >= sd%background_settings%rtor) then
            is_valid = .false.
            error_msg = "Plasma radius must be less than torus radius"
            return
        end if
        
        ! Check antenna is inside plasma
        if (sd%antenna_settings%ra >= sd%background_settings%rp) then
            is_valid = .false.
            error_msg = "Antenna radius must be less than plasma radius"
            return
        end if
        
    end subroutine settings_validate_consistency
    
    ! =========================================================================
    ! Error Handling Procedures
    ! =========================================================================
    
    !> @brief Format error message with context and specific parameter
    subroutine settings_format_error_message(error_code, context, parameter, error_msg)
        integer, intent(in) :: error_code
        character(len=*), intent(in) :: context
        character(len=*), intent(in) :: parameter
        character(len=*), intent(out) :: error_msg
        
        character(len=256) :: error_name
        integer :: ierr_dummy
        
        ! Get error name
        call settings_get_error_name(error_code, error_name, ierr_dummy)
        
        ! Format complete error message
        write(error_msg, '(A,": ",A," error in ",A," - ",A)') &
            trim(error_name), &
            trim(get_error_type_description(error_code)), &
            trim(context), &
            trim(parameter)
    end subroutine settings_format_error_message
    
    !> @brief Get error name from error code
    subroutine settings_get_error_name(error_code, error_name, ierr)
        integer, intent(in) :: error_code
        character(len=*), intent(out) :: error_name
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        select case (error_code)
        case (KILCA_SUCCESS)
            error_name = "KILCA_SUCCESS"
        case (KILCA_ERROR_INVALID_INPUT)
            error_name = "KILCA_ERROR_INVALID_INPUT"
        case (KILCA_ERROR_MEMORY)
            error_name = "KILCA_ERROR_MEMORY"
        case (KILCA_ERROR_FILE)
            error_name = "KILCA_ERROR_FILE"
        case (KILCA_ERROR_CONVERGENCE)
            error_name = "KILCA_ERROR_CONVERGENCE"
        case (KILCA_ERROR_BOUNDS)
            error_name = "KILCA_ERROR_BOUNDS"
        case (KILCA_ERROR_NOT_IMPLEMENTED)
            error_name = "KILCA_ERROR_NOT_IMPLEMENTED"
        case (KILCA_ERROR_FORMAT)
            error_name = "KILCA_ERROR_FORMAT"
        case (KILCA_ERROR)
            error_name = "KILCA_ERROR"
        case default
            error_name = "UNKNOWN_ERROR_CODE"
            ierr = KILCA_ERROR_INVALID_INPUT
        end select
    end subroutine settings_get_error_name
    
    !> @brief Validate settings with context information
    subroutine settings_validate_with_context(sd, context, error_msg, ierr)
        type(settings_t), intent(in) :: sd
        character(len=*), intent(in) :: context
        character(len=*), intent(out) :: error_msg
        integer, intent(out) :: ierr
        
        logical :: is_valid
        character(len=1024) :: validation_error
        
        ierr = KILCA_SUCCESS
        
        ! Validate complete settings
        call settings_validate_complete(sd, is_valid, validation_error, ierr)
        if (ierr /= KILCA_SUCCESS) then
            error_msg = "Validation failed in context: " // trim(context)
            return
        end if
        
        if (.not. is_valid) then
            ierr = KILCA_ERROR_INVALID_INPUT
            error_msg = "Context: " // trim(context) // " - " // trim(validation_error)
            return
        end if
        
        error_msg = ""
    end subroutine settings_validate_with_context
    
    !> @brief Attempt to recover from invalid antenna state
    subroutine settings_attempt_recovery(ant, recovered, ierr)
        type(antenna_t), intent(inout) :: ant
        logical, intent(out) :: recovered
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        recovered = .false.
        
        ! Try to fix invalid ra
        if (ant%ra < 0.0_dp) then
            ant%ra = 0.0_dp
            recovered = .true.
        end if
        
        ! Try to fix modes array allocation
        if (ant%dma > 0 .and. .not. allocated(ant%modes)) then
            allocate(ant%modes(2 * ant%dma))
            ant%modes = 0
            recovered = .true.
        end if
        
        ! Try to fix inconsistent array size
        if (ant%dma > 0 .and. allocated(ant%modes)) then
            if (size(ant%modes) /= 2 * ant%dma) then
                deallocate(ant%modes)
                allocate(ant%modes(2 * ant%dma))
                ant%modes = 0
                recovered = .true.
            end if
        end if
        
        ! Try to fix zero dma with allocated array
        if (ant%dma == 0 .and. allocated(ant%modes)) then
            deallocate(ant%modes)
            recovered = .true.
        end if
    end subroutine settings_attempt_recovery
    
    !> @brief Get detailed validation errors for background settings
    subroutine settings_get_detailed_validation_errors(bs, detailed_error, ierr)
        type(back_sett_t), intent(in) :: bs
        character(len=*), intent(out) :: detailed_error
        integer, intent(out) :: ierr
        
        integer :: error_count
        ! Variables is_valid, temp_error removed - were unused
        
        ierr = KILCA_SUCCESS
        detailed_error = ""
        error_count = 0
        
        ! Check rtor
        if (bs%rtor <= 0.0_dp) then
            error_count = error_count + 1
            if (len_trim(detailed_error) > 0) detailed_error = trim(detailed_error) // "; "
            detailed_error = trim(detailed_error) // "rtor must be positive"
        end if
        
        ! Check rp vs rtor
        if (bs%rp >= bs%rtor) then
            error_count = error_count + 1
            if (len_trim(detailed_error) > 0) detailed_error = trim(detailed_error) // "; "
            detailed_error = trim(detailed_error) // "rp must be less than rtor"
        end if
        
        ! Check B0
        if (bs%B0 <= 0.0_dp) then
            error_count = error_count + 1
            if (len_trim(detailed_error) > 0) detailed_error = trim(detailed_error) // "; "
            detailed_error = trim(detailed_error) // "B0 must be positive"
        end if
        
        ! Check N (must be odd)
        if (mod(bs%N, 2) == 0) then
            error_count = error_count + 1
            if (len_trim(detailed_error) > 0) detailed_error = trim(detailed_error) // "; "
            detailed_error = trim(detailed_error) // "N must be odd"
        end if
        
        if (error_count == 0) then
            detailed_error = "No validation errors found"
        end if
    end subroutine settings_get_detailed_validation_errors
    
    !> @brief Read antenna settings with error context
    subroutine antenna_read_settings_with_error_context(ant, path, error_msg, ierr)
        type(antenna_t), intent(out) :: ant
        character(len=*), intent(in) :: path
        character(len=*), intent(out) :: error_msg
        integer, intent(out) :: ierr
        
        logical :: file_exists
        
        ! Check if path exists
        inquire(file=trim(path), exist=file_exists)
        if (.not. file_exists) then
            ierr = KILCA_ERROR_FILE
            error_msg = "File not found: " // trim(path)
            return
        end if
        
        ! Try to read (simplified for testing)
        call antenna_read_settings(ant, path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            error_msg = "Failed to read antenna settings from: " // trim(path)
        else
            error_msg = ""
        end if
    end subroutine antenna_read_settings_with_error_context
    
    !> @brief Test allocation failure simulation
    subroutine settings_test_allocation_failure(sd, error_msg, ierr)
        type(settings_t), pointer, intent(out) :: sd
        character(len=*), intent(out) :: error_msg
        integer, intent(out) :: ierr
        
        ! Simulate allocation failure
        ierr = KILCA_ERROR_MEMORY
        error_msg = "Memory allocation failed during settings creation"
        sd => null()
    end subroutine settings_test_allocation_failure
    
    !> @brief Read settings from namelist file using Fortran intrinsic namelist I/O
    subroutine settings_read_namelist(filename, settings, ierr)
        character(len=*), intent(in) :: filename
        type(settings_t), intent(inout) :: settings
        integer, intent(out) :: ierr
        
        integer :: unit
        logical :: file_exists
        
        ! Local variables matching namelist exactly - antenna
        real(dp) :: ra, wa, I0
        complex(dp) :: flab
        integer :: dma, flag_debug_ant, flag_eigmode
        integer, dimension(50) :: modes           ! Increased size for advanced tests
        
        ! Local variables matching namelist exactly - background
        real(dp) :: rtor, rp, B0, V_gal_sys, m_i, V_scale, zele, zion
        integer :: calc_back, flag_debug_bg, N
        character(len=512) :: path2profiles, flag_back
        real(dp), dimension(20) :: mass, charge  ! Increased size for advanced tests
        
        ! Local variables matching namelist exactly - output
        integer :: flag_background, flag_emfield, flag_additional, flag_dispersion
        integer :: num_quants, flag_debug_out
        integer, dimension(30) :: flag_quants     ! Increased size for advanced tests
        
        ! Local variables matching namelist exactly - eigenmode
        character(len=512) :: fname
        integer :: search_flag, rdim, idim, stop_flag, test_roots, flag_debug_eig
        integer :: Nguess, kmin, kmax, n_zeros, use_winding
        real(dp) :: rfmin, rfmax, ifmin, ifmax, eps_res, eps_abs, eps_rel, delta
        complex(dp), dimension(20) :: fstart      ! Increased size for advanced complex tests
        
        ! Direct namelist declarations - simple and clean
        namelist /antenna/ ra, wa, I0, flab, dma, flag_debug_ant, flag_eigmode, modes
        namelist /background/ rtor, rp, B0, V_gal_sys, m_i, calc_back, flag_debug_bg, mass, charge, &
                              path2profiles, flag_back, N, V_scale, zele, zion
        namelist /output/ flag_background, flag_emfield, flag_additional, flag_dispersion, num_quants, &
                          flag_quants, flag_debug_out
        namelist /eigenmode/ fname, search_flag, rdim, idim, rfmin, rfmax, ifmin, ifmax, stop_flag, &
                             eps_res, eps_abs, eps_rel, delta, test_roots, flag_debug_eig, &
                             Nguess, kmin, kmax, n_zeros, use_winding, fstart
        
        ierr = KILCA_SUCCESS
        
        ! Check file existence first
        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            ierr = KILCA_ERROR_FILE_NOT_FOUND
            return
        end if
        
        ! Read namelists with detailed error checking
        open(newunit=unit, file=filename, status='old', action='read', iostat=ierr)
        if (ierr /= 0) then
            ierr = KILCA_ERROR_FILE_OPEN
            return
        end if
        
        ! Read antenna namelist
        read(unit, nml=antenna, iostat=ierr)
        if (ierr /= 0) then
            close(unit)
            if (ierr > 0) then
                ierr = KILCA_ERROR_FORMAT  ! Antenna section format error
            else if (ierr < 0) then
                ierr = KILCA_ERROR_FILE_NOT_FOUND  ! Antenna section missing
            else
                ierr = KILCA_ERROR_NAMELIST_READ + 1  ! Antenna section read error
            end if
            return
        end if
        
        ! Read background namelist  
        read(unit, nml=background, iostat=ierr)
        if (ierr /= 0) then
            close(unit)
            if (ierr > 0) then
                ierr = KILCA_ERROR_FORMAT  ! Background section format error
            else if (ierr < 0) then
                ierr = KILCA_ERROR_FILE_NOT_FOUND  ! Background section missing
            else
                ierr = KILCA_ERROR_NAMELIST_READ + 2  ! Background section read error
            end if
            return
        end if
        
        ! Read output namelist
        read(unit, nml=output, iostat=ierr)
        if (ierr /= 0) then
            close(unit)
            if (ierr > 0) then
                ierr = KILCA_ERROR_FORMAT  ! Output section format error
            else if (ierr < 0) then
                ierr = KILCA_ERROR_FILE_NOT_FOUND  ! Output section missing
            else
                ierr = KILCA_ERROR_NAMELIST_READ + 3  ! Output section read error
            end if
            return
        end if
        
        ! Read eigenmode namelist
        read(unit, nml=eigenmode, iostat=ierr)
        if (ierr /= 0) then
            close(unit)
            if (ierr > 0) then
                ierr = KILCA_ERROR_FORMAT  ! Eigenmode section format error
            else if (ierr < 0) then
                ierr = KILCA_ERROR_FILE_NOT_FOUND  ! Eigenmode section missing
            else
                ierr = KILCA_ERROR_NAMELIST_READ + 4  ! Eigenmode section read error
            end if
            return
        end if
        
        ! Direct assignment to settings structure
        ! Antenna settings
        settings%antenna_settings%ra = ra
        settings%antenna_settings%wa = wa
        settings%antenna_settings%I0 = I0
        settings%antenna_settings%flab = flab
        settings%antenna_settings%dma = dma
        settings%antenna_settings%flag_debug = flag_debug_ant
        settings%antenna_settings%flag_eigmode = flag_eigmode
        
        ! Antenna modes array - efficient allocation and copy
        if (allocated(settings%antenna_settings%modes)) then
            if (size(settings%antenna_settings%modes) /= dma * 2) then
                deallocate(settings%antenna_settings%modes)
                allocate(settings%antenna_settings%modes(dma * 2))  ! dma pairs of (m,n)
            end if
        else
            allocate(settings%antenna_settings%modes(dma * 2))  ! dma pairs of (m,n)
        end if
        settings%antenna_settings%modes(1:dma*2) = modes(1:dma*2)
        
        ! Background settings
        settings%background_settings%rtor = rtor
        settings%background_settings%rp = rp
        settings%background_settings%B0 = B0
        settings%background_settings%V_gal_sys = V_gal_sys
        settings%background_settings%m_i = m_i
        settings%background_settings%calc_back = calc_back
        settings%background_settings%flag_debug = flag_debug_bg
        settings%background_settings%N = N
        settings%background_settings%V_scale = V_scale
        settings%background_settings%zele = zele
        settings%background_settings%zion = zion
        
        ! Background string parameters - allocate and copy
        if (allocated(settings%background_settings%path2profiles)) then
            deallocate(settings%background_settings%path2profiles)
        end if
        settings%background_settings%path2profiles = trim(adjustl(path2profiles))
        
        if (allocated(settings%background_settings%flag_back)) then
            deallocate(settings%background_settings%flag_back)
        end if
        settings%background_settings%flag_back = trim(adjustl(flag_back))
        
        ! Background arrays - allocate and copy from namelist arrays
        if (allocated(settings%background_settings%mass)) then
            deallocate(settings%background_settings%mass)
        end if
        allocate(settings%background_settings%mass(3))
        settings%background_settings%mass(1:3) = mass(1:3)
        
        if (allocated(settings%background_settings%charge)) then
            deallocate(settings%background_settings%charge)
        end if
        allocate(settings%background_settings%charge(3))
        settings%background_settings%charge(1:3) = charge(1:3)
        
        ! Output settings
        settings%output_settings%flag_background = flag_background
        settings%output_settings%flag_emfield = flag_emfield
        settings%output_settings%flag_additional = flag_additional
        settings%output_settings%flag_dispersion = flag_dispersion
        settings%output_settings%num_quants = num_quants
        settings%output_settings%flag_debug = flag_debug_out
        
        ! Output arrays - efficient allocation and copy from namelist arrays
        if (allocated(settings%output_settings%flag_quants)) then
            if (size(settings%output_settings%flag_quants) /= num_quants) then
                deallocate(settings%output_settings%flag_quants)
                allocate(settings%output_settings%flag_quants(num_quants))
            end if
        else
            allocate(settings%output_settings%flag_quants(num_quants))
        end if
        settings%output_settings%flag_quants(1:num_quants) = flag_quants(1:num_quants)
        
        ! Eigenmode settings
        settings%eigmode_settings%search_flag = search_flag
        settings%eigmode_settings%rdim = rdim
        settings%eigmode_settings%idim = idim
        settings%eigmode_settings%rfmin = rfmin
        settings%eigmode_settings%rfmax = rfmax
        settings%eigmode_settings%ifmin = ifmin
        settings%eigmode_settings%ifmax = ifmax
        settings%eigmode_settings%stop_flag = stop_flag
        settings%eigmode_settings%eps_res = eps_res
        settings%eigmode_settings%eps_abs = eps_abs
        settings%eigmode_settings%eps_rel = eps_rel
        settings%eigmode_settings%delta = delta
        settings%eigmode_settings%test_roots = test_roots
        settings%eigmode_settings%flag_debug = flag_debug_eig
        settings%eigmode_settings%Nguess = Nguess
        settings%eigmode_settings%kmin = kmin
        settings%eigmode_settings%kmax = kmax
        settings%eigmode_settings%n_zeros = n_zeros
        settings%eigmode_settings%use_winding = use_winding
        
        ! Eigenmode string parameter - allocate and copy
        if (allocated(settings%eigmode_settings%fname)) then
            deallocate(settings%eigmode_settings%fname)
        end if
        settings%eigmode_settings%fname = trim(adjustl(fname))
        
        ! Eigenmode complex array - efficient allocation and copy
        if (allocated(settings%eigmode_settings%fstart)) then
            if (size(settings%eigmode_settings%fstart) /= Nguess) then
                deallocate(settings%eigmode_settings%fstart)
                allocate(settings%eigmode_settings%fstart(Nguess))
            end if
        else
            allocate(settings%eigmode_settings%fstart(Nguess))
        end if
        settings%eigmode_settings%fstart(1:Nguess) = fstart(1:Nguess)
        
        ! Comprehensive parameter range validation
        call validate_namelist_parameters(settings, ierr)
        if (ierr /= KILCA_SUCCESS) then
            close(unit)
            ! Add specific error context
            select case (ierr)
            case (KILCA_ERROR_INVALID_PARAMETER)
                ierr = KILCA_ERROR_INVALID_PARAMETER  ! Parameter out of physics range
            case default
                ierr = KILCA_ERROR_INVALID_PARAMETER  ! Generic parameter validation error
            end select
            return
        end if
        
        close(unit)
        
    end subroutine settings_read_namelist
    
    !> @brief Validate namelist parameters against physics ranges
    subroutine validate_namelist_parameters(settings, ierr)
        type(settings_t), intent(in) :: settings
        integer, intent(out) :: ierr
        integer :: i  ! Loop counter for array validation
        
        ierr = KILCA_SUCCESS
        
        ! Antenna validation
        if (settings%antenna_settings%ra <= 0.0_dp .or. settings%antenna_settings%ra > 500.0_dp) then
            ierr = KILCA_ERROR_INVALID_PARAMETER
            return
        end if
        
        if (settings%antenna_settings%I0 <= 0.0_dp .or. settings%antenna_settings%I0 > 1.0e15_dp) then
            ierr = KILCA_ERROR_INVALID_PARAMETER
            return
        end if
        
        if (settings%antenna_settings%dma < 0 .or. settings%antenna_settings%dma > 50) then
            ierr = KILCA_ERROR_INVALID_PARAMETER
            return
        end if
        
        ! Background validation  
        if (settings%background_settings%rtor <= 0.0_dp .or. settings%background_settings%rtor > 1000.0_dp) then
            ierr = KILCA_ERROR_INVALID_PARAMETER
            return
        end if
        
        if (settings%background_settings%rp <= 0.0_dp .or. settings%background_settings%rp > 500.0_dp) then
            ierr = KILCA_ERROR_INVALID_PARAMETER
            return
        end if
        
        if (settings%background_settings%rp >= settings%background_settings%rtor) then
            ierr = KILCA_ERROR_INVALID_PARAMETER  ! Plasma radius must be less than torus radius
            return
        end if
        
        if (settings%background_settings%B0 <= 0.0_dp .or. settings%background_settings%B0 > 100000.0_dp) then
            ierr = KILCA_ERROR_INVALID_PARAMETER
            return
        end if
        
        if (settings%background_settings%m_i <= 0.0_dp .or. settings%background_settings%m_i > 50.0_dp) then
            ierr = KILCA_ERROR_INVALID_PARAMETER
            return
        end if
        
        ! Enhanced validation - Array consistency checks
        if (allocated(settings%antenna_settings%modes)) then
            if (size(settings%antenna_settings%modes) /= settings%antenna_settings%dma * 2) then
                ierr = KILCA_ERROR_INVALID_PARAMETER  ! Array size inconsistent with dma
                return
            end if
        end if
        
        if (allocated(settings%output_settings%flag_quants)) then
            if (size(settings%output_settings%flag_quants) /= settings%output_settings%num_quants) then
                ierr = KILCA_ERROR_INVALID_PARAMETER  ! Array size inconsistent with num_quants
                return
            end if
        end if
        
        if (allocated(settings%eigmode_settings%fstart)) then
            if (size(settings%eigmode_settings%fstart) /= settings%eigmode_settings%Nguess) then
                ierr = KILCA_ERROR_INVALID_PARAMETER  ! Array size inconsistent with Nguess
                return
            end if
            
            ! Validate complex numbers are finite and not NaN
            do i = 1, size(settings%eigmode_settings%fstart)
                if (.not. (abs(real(settings%eigmode_settings%fstart(i))) < huge(1.0_dp) .and. &
                          abs(aimag(settings%eigmode_settings%fstart(i))) < huge(1.0_dp))) then
                    ierr = KILCA_ERROR_INVALID_PARAMETER  ! Complex number out of range
                    return
                end if
            end do
        end if
        
        ! String parameter validation
        if (allocated(settings%background_settings%path2profiles)) then
            if (len_trim(settings%background_settings%path2profiles) == 0) then
                ierr = KILCA_ERROR_INVALID_PARAMETER  ! Empty path string
                return
            end if
        end if
        
        ! Output validation
        if (settings%output_settings%num_quants < 0 .or. settings%output_settings%num_quants > 50) then
            ierr = KILCA_ERROR_INVALID_PARAMETER
            return
        end if
        
        ! Eigenmode validation
        if (settings%eigmode_settings%rdim <= 0 .or. settings%eigmode_settings%rdim > 10000) then
            ierr = KILCA_ERROR_INVALID_PARAMETER
            return
        end if
        
        if (settings%eigmode_settings%idim <= 0 .or. settings%eigmode_settings%idim > 10000) then
            ierr = KILCA_ERROR_INVALID_PARAMETER
            return
        end if
        
        if (settings%eigmode_settings%rfmin >= settings%eigmode_settings%rfmax) then
            ierr = KILCA_ERROR_INVALID_PARAMETER
            return
        end if
        
        if (settings%eigmode_settings%ifmin >= settings%eigmode_settings%ifmax) then
            ierr = KILCA_ERROR_INVALID_PARAMETER
            return
        end if
        
    end subroutine validate_namelist_parameters
    
    !> @brief Log error to file
    subroutine settings_log_error(error_code, context, parameter, log_file, ierr)
        integer, intent(in) :: error_code
        character(len=*), intent(in) :: context
        character(len=*), intent(in) :: parameter
        character(len=*), intent(in) :: log_file
        integer, intent(out) :: ierr
        
        integer :: unit_num
        character(len=1024) :: error_msg
        
        ierr = KILCA_SUCCESS
        
        ! Format error message
        call settings_format_error_message(error_code, context, parameter, error_msg)
        
        ! Open log file for writing
        open(newunit=unit_num, file=trim(log_file), status='unknown', &
             position='append', iostat=ierr)
        if (ierr /= 0) then
            ierr = KILCA_ERROR_FILE
            return
        end if
        
        ! Write error to log
        write(unit_num, '(A)') trim(error_msg)
        close(unit_num)
        
        ierr = KILCA_SUCCESS
    end subroutine settings_log_error
    
    !> @brief Check if error log exists
    subroutine settings_check_error_log_exists(log_file, logged, ierr)
        character(len=*), intent(in) :: log_file
        logical, intent(out) :: logged
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! Check if log file exists
        inquire(file=trim(log_file), exist=logged)
    end subroutine settings_check_error_log_exists
    
    !> @brief Clear error log file
    subroutine settings_clear_error_log(log_file, ierr)
        character(len=*), intent(in) :: log_file
        integer, intent(out) :: ierr
        
        integer :: unit_num
        logical :: file_exists
        
        ierr = KILCA_SUCCESS
        
        ! Check if file exists
        inquire(file=trim(log_file), exist=file_exists)
        if (.not. file_exists) return
        
        ! Try to delete the file by opening and closing with status='delete'
        open(newunit=unit_num, file=trim(log_file), status='old', iostat=ierr)
        if (ierr == 0) then
            close(unit_num, status='delete', iostat=ierr)
            ! If deletion fails, just continue - not critical for test
            if (ierr /= 0) ierr = KILCA_SUCCESS
        else
            ierr = KILCA_SUCCESS  ! File doesn't exist, no need to delete
        end if
    end subroutine settings_clear_error_log
    
    ! =========================================================================
    ! Private Helper Functions
    ! =========================================================================
    
    !> @brief Get error type description
    function get_error_type_description(error_code) result(description)
        integer, intent(in) :: error_code
        character(len=256) :: description
        
        select case (error_code)
        case (KILCA_SUCCESS)
            description = "Success"
        case (KILCA_ERROR_INVALID_INPUT)
            description = "Invalid input"
        case (KILCA_ERROR_MEMORY)
            description = "Memory"
        case (KILCA_ERROR_FILE)
            description = "File"
        case (KILCA_ERROR_CONVERGENCE)
            description = "Convergence"
        case (KILCA_ERROR_BOUNDS)
            description = "Bounds"
        case (KILCA_ERROR_NOT_IMPLEMENTED)
            description = "Not implemented"
        case (KILCA_ERROR_FORMAT)
            description = "Format"
        case (KILCA_ERROR)
            description = "General error"
        case default
            description = "Unrecognized"
        end select
    end function get_error_type_description
    
    !> @brief Read settings with automatic format detection
    !> @details Automatically detects whether to use namelist or legacy format
    !>          based on file presence. Priority: namelist > legacy > error
    !>          Enhanced with comprehensive error handling and validation.
    subroutine read_settings_auto_detect(path, settings_out, ierr)
        character(len=*), intent(in) :: path
        type(settings_t), intent(out) :: settings_out
        integer, intent(out) :: ierr
        
        logical :: path_exists, file_readable
        integer :: format_result
        character(len=1024) :: error_context
        
        ierr = KILCA_SUCCESS
        error_context = ""
        
        ! Validate input path
        call validate_settings_path(path, ierr, error_context)
        if (ierr /= KILCA_SUCCESS) then
            settings_out%format_used = FORMAT_AUTO_DETECT
            return
        end if
        
        ! Initialize settings structure
        call initialize_settings_structure(settings_out, path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            error_context = "Failed to initialize settings structure"
            return
        end if
        
        ! Attempt format detection with comprehensive error handling
        call detect_and_read_format(path, settings_out, format_result, ierr, error_context)
        
        ! Handle detection results
        select case (format_result)
        case (NAMELIST_FORMAT)
            if (ierr == KILCA_SUCCESS) then
                print *, "Auto-detect: Successfully loaded namelist format"
            else
                print *, "Auto-detect: Namelist format detected but reading failed:", trim(error_context)
            end if
        case (LEGACY_FORMAT)  
            if (ierr == KILCA_SUCCESS) then
                print *, "Auto-detect: Successfully loaded legacy format"
            else
                print *, "Auto-detect: Legacy format detected but reading failed:", trim(error_context)
            end if
        case default
            print *, "Auto-detect: No valid settings files found in ", trim(path)
            ierr = KILCA_ERROR_FILE_NOT_FOUND
            settings_out%format_used = FORMAT_AUTO_DETECT
        end select
        
    end subroutine read_settings_auto_detect
    
    !> @brief Validate settings path for format detection
    subroutine validate_settings_path(path, ierr, error_msg)
        character(len=*), intent(in) :: path
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: error_msg
        
        logical :: path_exists
        
        ierr = KILCA_SUCCESS
        error_msg = ""
        
        ! Check path validity
        if (len_trim(path) == 0) then
            ierr = KILCA_ERROR_INVALID_INPUT
            error_msg = "Empty path provided"
            return
        end if
        
        ! Check if path exists and is accessible
        inquire(file=trim(path), exist=path_exists)
        if (.not. path_exists) then
            ierr = KILCA_ERROR_FILE_NOT_FOUND
            error_msg = "Settings path does not exist: " // trim(path)
            return
        end if
        
    end subroutine validate_settings_path
    
    !> @brief Initialize settings structure for format detection
    subroutine initialize_settings_structure(settings, path, ierr)
        type(settings_t), intent(out) :: settings
        character(len=*), intent(in) :: path
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! Initialize with default values
        settings%format_used = FORMAT_AUTO_DETECT
        
        ! Safely allocate path
        if (allocated(settings%path2project)) deallocate(settings%path2project)
        allocate(character(len=len_trim(path)) :: settings%path2project)
        settings%path2project = trim(path)
        
    end subroutine initialize_settings_structure
    
    !> @brief Detect format and attempt to read settings
    subroutine detect_and_read_format(path, settings, format_result, ierr, error_msg)
        character(len=*), intent(in) :: path
        type(settings_t), intent(inout) :: settings
        integer, intent(out) :: format_result
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: error_msg
        
        logical :: file_exists, file_readable
        character(len=1024) :: test_file
        
        ierr = KILCA_SUCCESS
        error_msg = ""
        format_result = FORMAT_AUTO_DETECT
        
        ! Test for namelist format (highest priority)
        test_file = trim(path) // "/settings.conf"
        call validate_settings_file(test_file, file_exists, file_readable, ierr)
        
        if (file_exists .and. file_readable) then
            format_result = NAMELIST_FORMAT
            settings%format_used = NAMELIST_FORMAT
            call settings_read_namelist(test_file, settings, ierr)
            if (ierr /= KILCA_SUCCESS) then
                error_msg = "Namelist file exists but reading failed"
            end if
            return
        end if
        
        ! Test for legacy format (second priority)
        call check_legacy_format_availability(path, file_exists, ierr)
        
        if (file_exists) then
            format_result = LEGACY_FORMAT
            settings%format_used = LEGACY_FORMAT
            call read_settings_legacy_format(path, settings, ierr)
            if (ierr /= KILCA_SUCCESS) then
                error_msg = "Legacy files exist but reading failed"
            end if
            return
        end if
        
        ! No valid format found
        format_result = FORMAT_AUTO_DETECT
        ierr = KILCA_ERROR_FILE_NOT_FOUND
        error_msg = "No valid settings files found"
        
    end subroutine detect_and_read_format
    
    !> @brief Validate individual settings file
    subroutine validate_settings_file(filename, exists, readable, ierr)
        character(len=*), intent(in) :: filename
        logical, intent(out) :: exists, readable
        integer, intent(out) :: ierr
        
        integer :: unit, iostat
        
        ierr = KILCA_SUCCESS
        exists = .false.
        readable = .false.
        
        ! Check file existence
        inquire(file=filename, exist=exists)
        
        if (exists) then
            ! Test readability by attempting to open
            open(newunit=unit, file=filename, status='old', action='read', iostat=iostat)
            if (iostat == 0) then
                readable = .true.
                close(unit)
            else
                readable = .false.
                ierr = KILCA_ERROR_FILE_PERMISSION
            end if
        end if
        
    end subroutine validate_settings_file
    
    !> @brief Check if legacy format files are available
    subroutine check_legacy_format_availability(path, available, ierr)
        character(len=*), intent(in) :: path
        logical, intent(out) :: available
        integer, intent(out) :: ierr
        
        logical :: antenna_exists, background_exists
        character(len=1024) :: antenna_file, background_file
        
        ierr = KILCA_SUCCESS
        available = .false.
        
        ! Check for essential legacy files
        antenna_file = trim(path) // "/antenna.in"
        background_file = trim(path) // "/background.in"
        
        inquire(file=antenna_file, exist=antenna_exists)
        inquire(file=background_file, exist=background_exists)
        
        ! Both antenna and background files must exist for valid legacy format
        available = antenna_exists .and. background_exists
        
        if (.not. available .and. (antenna_exists .or. background_exists)) then
            ! Incomplete legacy format detected
            ierr = KILCA_ERROR_FORMAT
        end if
        
    end subroutine check_legacy_format_availability
    
    !> @brief Read settings from legacy format files
    !> @details Reads settings from individual .in files in legacy format
    subroutine read_settings_legacy_format(path, settings_out, ierr)
        character(len=*), intent(in) :: path
        type(settings_t), intent(out) :: settings_out
        integer, intent(out) :: ierr
        
        character(len=1024) :: error_context
        
        ierr = KILCA_SUCCESS
        error_context = ""
        
        ! Initialize settings with defaults
        call antenna_initialize_defaults(settings_out%antenna_settings, ierr)
        call back_sett_initialize_defaults(settings_out%background_settings, ierr)
        call output_sett_initialize_defaults(settings_out%output_settings, ierr)
        call eigmode_sett_initialize_defaults(settings_out%eigmode_settings, ierr)
        
        ! Read antenna settings
        call read_legacy_antenna(path, settings_out%antenna_settings, ierr, error_context)
        if (ierr /= KILCA_SUCCESS) then
            print *, "Legacy format: Failed to read antenna.in -", trim(error_context)
            return
        end if
        
        ! Read background settings
        call read_legacy_background(path, settings_out%background_settings, ierr, error_context)
        if (ierr /= KILCA_SUCCESS) then
            print *, "Legacy format: Failed to read background.in -", trim(error_context)
            return
        end if
        
        ! Read output settings
        call read_legacy_output(path, settings_out%output_settings, ierr, error_context)
        if (ierr /= KILCA_SUCCESS) then
            print *, "Legacy format: Failed to read output.in -", trim(error_context)
            return
        end if
        
        ! Read eigenmode settings
        call read_legacy_eigenmode(path, settings_out%eigmode_settings, ierr, error_context)
        if (ierr /= KILCA_SUCCESS) then
            print *, "Legacy format: Failed to read eigenmode.in -", trim(error_context)
            return
        end if
        
        print *, "Legacy format: Successfully loaded all settings"
        
    end subroutine read_settings_legacy_format
    
    !> @brief Read antenna settings from legacy format
    subroutine read_legacy_antenna(path, antenna, ierr, error_msg)
        character(len=*), intent(in) :: path
        type(antenna_t), intent(inout) :: antenna
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: error_msg
        
        character(len=1024) :: filename
        integer :: unit, iostat, i
        real(dp) :: flab_real, flab_imag
        integer, dimension(20) :: modes_temp
        
        ierr = KILCA_SUCCESS
        error_msg = ""
        
        filename = trim(path) // "/antenna.in"
        
        ! Open file
        open(newunit=unit, file=filename, status='old', action='read', iostat=iostat)
        if (iostat /= 0) then
            ierr = KILCA_ERROR_FILE_NOT_FOUND
            error_msg = "Cannot open antenna.in"
            return
        end if
        
        ! Read ra, wa, I0
        read(unit, *, iostat=iostat) antenna%ra, antenna%wa, antenna%I0
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading ra, wa, I0"
            return
        end if
        
        ! Read flab (complex as two reals)
        read(unit, *, iostat=iostat) flab_real, flab_imag
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading flab"
            return
        end if
        antenna%flab = cmplx(flab_real, flab_imag, dp)
        
        ! Read dma
        read(unit, *, iostat=iostat) antenna%dma
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading dma"
            return
        end if
        
        ! Read flag_debug and flag_eigmode
        read(unit, *, iostat=iostat) antenna%flag_debug, antenna%flag_eigmode
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading flags"
            return
        end if
        
        ! Read modes array
        if (antenna%dma > 0) then
            read(unit, *, iostat=iostat) (modes_temp(i), i=1, antenna%dma*2)
            if (iostat /= 0) then
                close(unit)
                ierr = KILCA_ERROR_FORMAT
                error_msg = "Error reading modes"
                return
            end if
            
            ! Allocate and copy modes
            if (allocated(antenna%modes)) deallocate(antenna%modes)
            allocate(antenna%modes(antenna%dma * 2))
            antenna%modes(1:antenna%dma*2) = modes_temp(1:antenna%dma*2)
        end if
        
        close(unit)
        
    end subroutine read_legacy_antenna
    
    !> @brief Read background settings from legacy format
    subroutine read_legacy_background(path, background, ierr, error_msg)
        character(len=*), intent(in) :: path
        type(back_sett_t), intent(inout) :: background
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: error_msg
        
        character(len=1024) :: filename
        character(len=256) :: path_temp, flag_temp
        integer :: unit, iostat
        real(dp), dimension(10) :: mass_temp, charge_temp
        
        ierr = KILCA_SUCCESS
        error_msg = ""
        
        filename = trim(path) // "/background.in"
        
        ! Open file
        open(newunit=unit, file=filename, status='old', action='read', iostat=iostat)
        if (iostat /= 0) then
            ierr = KILCA_ERROR_FILE_NOT_FOUND
            error_msg = "Cannot open background.in"
            return
        end if
        
        ! Read rtor, rp, B0
        read(unit, *, iostat=iostat) background%rtor, background%rp, background%B0
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading rtor, rp, B0"
            return
        end if
        
        ! Read path2profiles
        read(unit, '(a)', iostat=iostat) path_temp
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading path2profiles"
            return
        end if
        background%path2profiles = trim(adjustl(path_temp))
        
        ! Read calc_back
        read(unit, *, iostat=iostat) background%calc_back
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading calc_back"
            return
        end if
        
        ! Read flag_back
        read(unit, '(a)', iostat=iostat) flag_temp
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading flag_back"
            return
        end if
        background%flag_back = trim(adjustl(flag_temp))
        
        ! Read N
        read(unit, *, iostat=iostat) background%N
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading N"
            return
        end if
        
        ! Read V_gal_sys, V_scale, m_i, zele, zion
        read(unit, *, iostat=iostat) background%V_gal_sys, background%V_scale, &
                                     background%m_i, background%zele, background%zion
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading V_gal_sys, V_scale, m_i, zele, zion"
            return
        end if
        
        ! Read flag_debug
        read(unit, *, iostat=iostat) background%flag_debug
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading flag_debug"
            return
        end if
        
        ! Read mass array
        read(unit, *, iostat=iostat) mass_temp(1:3)
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading mass array"
            return
        end if
        
        ! Read charge array
        read(unit, *, iostat=iostat) charge_temp(1:3)
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading charge array"
            return
        end if
        
        ! Allocate and copy arrays
        if (allocated(background%mass)) deallocate(background%mass)
        allocate(background%mass(3))
        background%mass = mass_temp(1:3)
        
        if (allocated(background%charge)) deallocate(background%charge)
        allocate(background%charge(3))
        background%charge = charge_temp(1:3)
        
        close(unit)
        
    end subroutine read_legacy_background
    
    !> @brief Read output settings from legacy format
    subroutine read_legacy_output(path, output, ierr, error_msg)
        character(len=*), intent(in) :: path
        type(output_sett_t), intent(inout) :: output
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: error_msg
        
        character(len=1024) :: filename
        integer :: unit, iostat, i
        integer, dimension(20) :: flag_quants_temp
        
        ierr = KILCA_SUCCESS
        error_msg = ""
        
        filename = trim(path) // "/output.in"
        
        ! Open file
        open(newunit=unit, file=filename, status='old', action='read', iostat=iostat)
        if (iostat /= 0) then
            ierr = KILCA_ERROR_FILE_NOT_FOUND
            error_msg = "Cannot open output.in"
            return
        end if
        
        ! Read flags
        read(unit, *, iostat=iostat) output%flag_background, output%flag_emfield, &
                                     output%flag_additional, output%flag_dispersion
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading output flags"
            return
        end if
        
        ! Read flag_debug
        read(unit, *, iostat=iostat) output%flag_debug
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading flag_debug"
            return
        end if
        
        ! Read num_quants
        read(unit, *, iostat=iostat) output%num_quants
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading num_quants"
            return
        end if
        
        ! Read flag_quants array
        if (output%num_quants > 0) then
            read(unit, *, iostat=iostat) (flag_quants_temp(i), i=1, output%num_quants)
            if (iostat /= 0) then
                close(unit)
                ierr = KILCA_ERROR_FORMAT
                error_msg = "Error reading flag_quants"
                return
            end if
            
            ! Allocate and copy flag_quants
            if (allocated(output%flag_quants)) deallocate(output%flag_quants)
            allocate(output%flag_quants(output%num_quants))
            output%flag_quants = flag_quants_temp(1:output%num_quants)
        end if
        
        close(unit)
        
    end subroutine read_legacy_output
    
    !> @brief Read eigenmode settings from legacy format
    subroutine read_legacy_eigenmode(path, eigmode, ierr, error_msg)
        character(len=*), intent(in) :: path
        type(eigmode_sett_t), intent(inout) :: eigmode
        integer, intent(out) :: ierr
        character(len=*), intent(out) :: error_msg
        
        character(len=1024) :: filename
        character(len=256) :: fname_temp
        integer :: unit, iostat, i
        real(dp) :: fstart_real, fstart_imag
        complex(dp), dimension(20) :: fstart_temp
        
        ierr = KILCA_SUCCESS
        error_msg = ""
        
        filename = trim(path) // "/eigenmode.in"
        
        ! Open file
        open(newunit=unit, file=filename, status='old', action='read', iostat=iostat)
        if (iostat /= 0) then
            ierr = KILCA_ERROR_FILE_NOT_FOUND
            error_msg = "Cannot open eigenmode.in"
            return
        end if
        
        ! Read fname
        read(unit, '(a)', iostat=iostat) fname_temp
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading fname"
            return
        end if
        eigmode%fname = trim(adjustl(fname_temp))
        
        ! Read search_flag
        read(unit, *, iostat=iostat) eigmode%search_flag
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading search_flag"
            return
        end if
        
        ! Read rdim, idim
        read(unit, *, iostat=iostat) eigmode%rdim, eigmode%idim
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading rdim, idim"
            return
        end if
        
        ! Read rfmin, rfmax
        read(unit, *, iostat=iostat) eigmode%rfmin, eigmode%rfmax
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading rfmin, rfmax"
            return
        end if
        
        ! Read ifmin, ifmax
        read(unit, *, iostat=iostat) eigmode%ifmin, eigmode%ifmax
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading ifmin, ifmax"
            return
        end if
        
        ! Read stop_flag
        read(unit, *, iostat=iostat) eigmode%stop_flag
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading stop_flag"
            return
        end if
        
        ! Read eps_res, eps_abs, eps_rel, delta
        read(unit, *, iostat=iostat) eigmode%eps_res, eigmode%eps_abs, &
                                     eigmode%eps_rel, eigmode%delta
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading tolerances"
            return
        end if
        
        ! Read test_roots, flag_debug
        read(unit, *, iostat=iostat) eigmode%test_roots, eigmode%flag_debug
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading test_roots, flag_debug"
            return
        end if
        
        ! Read Nguess, kmin, kmax, n_zeros
        read(unit, *, iostat=iostat) eigmode%Nguess, eigmode%kmin, &
                                     eigmode%kmax, eigmode%n_zeros
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading Nguess, kmin, kmax, n_zeros"
            return
        end if
        
        ! Read use_winding
        read(unit, *, iostat=iostat) eigmode%use_winding
        if (iostat /= 0) then
            close(unit)
            ierr = KILCA_ERROR_FORMAT
            error_msg = "Error reading use_winding"
            return
        end if
        
        ! Read fstart array (complex numbers as pairs of reals)
        if (eigmode%Nguess > 0) then
            do i = 1, eigmode%Nguess
                read(unit, *, iostat=iostat) fstart_real, fstart_imag
                if (iostat /= 0) then
                    close(unit)
                    ierr = KILCA_ERROR_FORMAT
                    error_msg = "Error reading fstart array"
                    return
                end if
                fstart_temp(i) = cmplx(fstart_real, fstart_imag, dp)
            end do
            
            ! Allocate and copy fstart
            if (allocated(eigmode%fstart)) deallocate(eigmode%fstart)
            allocate(eigmode%fstart(eigmode%Nguess))
            eigmode%fstart = fstart_temp(1:eigmode%Nguess)
        end if
        
        close(unit)
        
    end subroutine read_legacy_eigenmode
    
end module kilca_settings_m