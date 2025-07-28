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
    
end module kilca_settings_m