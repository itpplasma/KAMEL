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
    public :: settings_read_all
    public :: settings_read_all_full
    public :: settings_print_all
    public :: settings_validate
    public :: settings_validate_complete
    public :: settings_validate_consistency
    
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
    public :: antenna_validate
    
    ! Background procedures
    public :: back_sett_set_calc_flag
    public :: back_sett_read_settings
    public :: back_sett_read_settings_full
    public :: back_sett_print_settings
    public :: back_sett_print_settings_to_unit
    public :: back_sett_validate
    
    ! Output procedures
    public :: output_sett_set_flags
    public :: output_sett_read_settings
    public :: output_sett_read_settings_full
    public :: output_sett_print_settings
    public :: output_sett_print_settings_to_unit
    public :: output_sett_validate
    
    ! Eigenmode procedures
    public :: eigmode_sett_set_search_flag
    public :: eigmode_sett_read_settings
    public :: eigmode_sett_read_settings_full
    public :: eigmode_sett_print_settings
    public :: eigmode_sett_print_settings_to_unit
    public :: eigmode_sett_validate
    
    ! C interface procedures
    public :: set_antenna_settings_c
    public :: set_background_settings_c
    public :: set_particles_settings_c
    public :: set_huge_factor_c
    public :: copy_antenna_data_to_antenna_module
    public :: copy_background_data_to_background_module
    
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
        
        ! Allocate settings structure
        allocate(sd, stat=alloc_stat)
        if (alloc_stat /= 0) then
            ierr = KILCA_ERROR_MEMORY
            sd => null()
            return
        end if
        
        ! Initialize path
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
        type(antenna_t), intent(inout) :: ant
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
        type(back_sett_t), intent(inout) :: bs
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
        type(output_sett_t), intent(inout) :: os
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
        type(eigmode_sett_t), intent(inout) :: es
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
                write(unit, '(A,I0,A,ES15.8,A,ES15.8,A)', iostat=ierr) "k=", k-1, "\tf=(", real(es%fstart(k)), ", ", aimag(es%fstart(k)), ")"
                if (ierr /= 0) return
            end do
        end if
        
        write(unit, '(A)', iostat=ierr) ""
        if (ierr /= 0) return
        
        ierr = KILCA_SUCCESS
        
    end subroutine eigmode_sett_print_settings_to_unit
    
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
        type(antenna_t), pointer, intent(in) :: ant
        
        ! In full implementation, would copy to a module-level variable
        ! that's accessible from Fortran routines that need it
        
    end subroutine copy_antenna_data_to_antenna_module
    
    !> @brief Copy background data to module (placeholder)
    subroutine copy_background_data_to_background_module(bs)
        type(back_sett_t), pointer, intent(in) :: bs
        
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
        if (ant%flag_debug < 0 .or. ant%flag_debug > 1) then
            is_valid = .false.
            error_msg = "Debug flag must be 0 or 1, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') ant%flag_debug
            return
        end if
        
        if (ant%flag_eigmode < 0 .or. ant%flag_eigmode > 1) then
            is_valid = .false.
            error_msg = "Eigenmode flag must be 0 or 1, got: "
            write(error_msg(len_trim(error_msg)+1:), '(i0)') ant%flag_eigmode
            return
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
        
        ! Validate plasma radius
        if (bs%rp <= 0.0_dp) then
            is_valid = .false.
            error_msg = "Plasma radius (rp) must be positive, got: "
            write(error_msg(len_trim(error_msg)+1:), '(g0)') bs%rp
            return
        end if
        
        ! Validate magnetic field
        if (bs%B0 <= 0.0_dp) then
            is_valid = .false.
            error_msg = "Magnetic field (B0) must be positive, got: "
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
        
        ! Validate velocity scale
        if (bs%V_scale == 0.0_dp) then
            is_valid = .false.
            error_msg = "Velocity scale (V_scale) cannot be zero"
            return
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
    
end module kilca_settings_m