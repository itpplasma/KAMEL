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
        !> Background calculation flag (>0: from files, <0: from interface, 0: error)
        integer :: calc_back = 1
        
        !> Number of radial grid points
        integer :: n_grid = 100
        
        !> Inner radius
        real(dp) :: r_min = 0.0_dp
        
        !> Outer radius
        real(dp) :: r_max = 100.0_dp
        
        !> Grid type (0=uniform, 1=adaptive)
        integer :: grid_type = 0
        
        !> Profile file names
        character(len=MAX_PATH_LEN) :: density_file = ""
        character(len=MAX_PATH_LEN) :: temperature_file = ""
        character(len=MAX_PATH_LEN) :: pressure_file = ""
        character(len=MAX_PATH_LEN) :: magnetic_file = ""
    end type back_sett_t
    
    ! =========================================================================
    ! Output Settings Type (from output_sett.h)
    ! =========================================================================
    
    !> @brief Output control settings
    type, public :: output_sett_t
        !> Save background profiles flag
        logical :: save_profiles = .true.
        
        !> Save electromagnetic fields flag
        logical :: save_fields = .true.
        
        !> Save matrix elements flag
        logical :: save_matrix = .false.
        
        !> Save eigenvalues flag
        logical :: save_eigenvalues = .true.
        
        !> Output format (0=text, 1=binary, 2=HDF5)
        integer :: output_format = 0
        
        !> Output directory
        character(len=MAX_PATH_LEN) :: output_dir = "output/"
        
        !> Verbosity level (0=quiet, 1=normal, 2=verbose, 3=debug)
        integer :: verbosity = 1
    end type output_sett_t
    
    ! =========================================================================
    ! Eigenmode Settings Type (from eigmode_sett.h)
    ! =========================================================================
    
    !> @brief Eigenmode search settings
    type, public :: eigmode_sett_t
        !> Search flag (1=frequency scan, 0=zero search, -1=all zeros)
        integer :: search_flag = 0
        
        !> Minimum frequency for scan
        real(dp) :: freq_min = 0.0_dp
        
        !> Maximum frequency for scan
        real(dp) :: freq_max = 1.0e9_dp
        
        !> Number of frequency points
        integer :: n_freq = 100
        
        !> Convergence tolerance
        real(dp) :: tolerance = 1.0e-6_dp
        
        !> Maximum iterations
        integer :: max_iter = 100
        
        !> Initial guess for complex frequency
        complex(dp) :: omega_guess = cmplx_zero
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
    public :: settings_read_all
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
    public :: antenna_print_settings
    
    ! Background procedures
    public :: back_sett_set_calc_flag
    public :: back_sett_read_settings
    public :: back_sett_print_settings
    
    ! Output procedures
    public :: output_sett_set_flags
    public :: output_sett_read_settings
    public :: output_sett_print_settings
    
    ! Eigenmode procedures
    public :: eigmode_sett_set_search_flag
    public :: eigmode_sett_read_settings
    public :: eigmode_sett_print_settings
    
    ! C interface procedures
    public :: set_antenna_settings_c
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
        
        ! Deallocate antenna modes array if allocated
        if (allocated(sd%antenna_settings%modes)) then
            deallocate(sd%antenna_settings%modes)
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
        
        if (sd%background_settings%r_max <= sd%background_settings%r_min) then
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
        print *, "  Calc flag:       ", bs%calc_back
        print *, "  Grid points:     ", bs%n_grid
        print *, "  R min/max:       ", bs%r_min, bs%r_max
        print *, "  Grid type:       ", bs%grid_type
        
    end subroutine back_sett_print_settings
    
    ! =========================================================================
    ! Output Settings Procedures
    ! =========================================================================
    
    !> @brief Set output flags
    subroutine output_sett_set_flags(os, save_profiles, save_fields, ierr)
        type(output_sett_t), intent(inout) :: os
        logical, intent(in), optional :: save_profiles, save_fields
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        if (present(save_profiles)) os%save_profiles = save_profiles
        if (present(save_fields)) os%save_fields = save_fields
        
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
        print *, "  Save profiles:   ", os%save_profiles
        print *, "  Save fields:     ", os%save_fields
        print *, "  Save matrix:     ", os%save_matrix
        print *, "  Output format:   ", os%output_format
        print *, "  Output dir:      ", trim(os%output_dir)
        print *, "  Verbosity:       ", os%verbosity
        
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
        print *, "  Freq min/max:    ", es%freq_min, es%freq_max
        print *, "  Freq points:     ", es%n_freq
        print *, "  Tolerance:       ", es%tolerance
        print *, "  Max iterations:  ", es%max_iter
        
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
    
end module kilca_settings_m