!> @file kilca_settings_m.f90
!> @brief Settings management module for KiLCA Fortran implementation
!> @details This module provides the Fortran equivalent of the C++ settings class
!>          and all related settings structures (antenna, background, output, eigmode).

module kilca_settings_m
    use iso_fortran_env, only: real64, int32, int64
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
    
    ! Background procedures
    public :: back_sett_set_calc_flag
    public :: back_sett_read_settings
    public :: back_sett_read_settings_full
    public :: back_sett_print_settings
    
    ! Output procedures
    public :: output_sett_set_flags
    public :: output_sett_read_settings
    public :: output_sett_read_settings_full
    public :: output_sett_print_settings
    
    ! Eigenmode procedures
    public :: eigmode_sett_set_search_flag
    public :: eigmode_sett_read_settings
    public :: eigmode_sett_read_settings_full
    public :: eigmode_sett_print_settings
    
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
    subroutine settings_print_all(sd)
        type(settings_t), intent(in) :: sd
        
        print *, "=== KiLCA Settings ==="
        print *, "Project path: ", trim(sd%path2project)
        
        call antenna_print_settings(sd%antenna_settings)
        call back_sett_print_settings(sd%background_settings)
        call output_sett_print_settings(sd%output_settings)
        call eigmode_sett_print_settings(sd%eigmode_settings)
        
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
    
    !> @brief Print antenna settings
    subroutine antenna_print_settings(ant)
        type(antenna_t), intent(in) :: ant
        
        print *, "--- Antenna Settings ---"
        print *, "  Radius (ra):     ", ant%ra, " cm"
        print *, "  Width (wa):      ", ant%wa, " cm"
        print *, "  Current (I0):    ", ant%I0, " statamp"
        print *, "  Frequency (flab):", ant%flab, " Hz"
        print *, "  Number of modes: ", ant%dma
        if (allocated(ant%modes)) then
            print *, "  Modes (m,n):    ", ant%modes
        end if
        print *, "  Debug flag:      ", ant%flag_debug
        print *, "  Eigenmode flag:  ", ant%flag_eigmode
        
    end subroutine antenna_print_settings
    
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
    
    !> @brief Print background settings
    subroutine back_sett_print_settings(bs)
        type(back_sett_t), intent(in) :: bs
        
        print *, "--- Background Settings ---"
        print *, "  Machine rtor:    ", bs%rtor
        print *, "  Plasma radius:   ", bs%rp
        print *, "  B0:              ", bs%B0
        print *, "  Calc flag:       ", bs%calc_back
        if (allocated(bs%flag_back)) then
            print *, "  Flag back:       ", trim(bs%flag_back)
        end if
        print *, "  Spline degree N: ", bs%N
        print *, "  V_gal_sys:       ", bs%V_gal_sys
        print *, "  V_scale:         ", bs%V_scale
        print *, "  Ion mass m_i:    ", bs%m_i
        print *, "  Collision coeffs:", bs%zele, bs%zion
        print *, "  Huge factor:     ", bs%huge_factor
        
    end subroutine back_sett_print_settings
    
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
    
    !> @brief Print output settings
    subroutine output_sett_print_settings(os)
        type(output_sett_t), intent(in) :: os
        
        print *, "--- Output Settings ---"
        print *, "  Flag background: ", os%flag_background
        print *, "  Flag EM field:   ", os%flag_emfield
        print *, "  Flag additional: ", os%flag_additional
        print *, "  Flag dispersion: ", os%flag_dispersion
        print *, "  Flag debug:      ", os%flag_debug
        print *, "  Num quants:      ", os%num_quants
        if (allocated(os%flag_quants)) then
            print *, "  Flag quants:     ", os%flag_quants
        end if
        
    end subroutine output_sett_print_settings
    
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
    
    !> @brief Print eigenmode settings
    subroutine eigmode_sett_print_settings(es)
        type(eigmode_sett_t), intent(in) :: es
        
        print *, "--- Eigenmode Settings ---"
        print *, "  Search flag:     ", es%search_flag
        if (allocated(es%fname)) then
            print *, "  Output file:     ", trim(es%fname)
        end if
        print *, "  Real grid:       ", es%rdim, es%rfmin, es%rfmax
        print *, "  Imag grid:       ", es%idim, es%ifmin, es%ifmax
        print *, "  Accuracies:      ", es%eps_res, es%eps_abs, es%eps_rel
        print *, "  Delta:           ", es%delta
        print *, "  Nguess:          ", es%Nguess
        print *, "  N zeros:         ", es%n_zeros
        print *, "  Use winding:     ", es%use_winding
        
    end subroutine eigmode_sett_print_settings
    
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
    
end module kilca_settings_m