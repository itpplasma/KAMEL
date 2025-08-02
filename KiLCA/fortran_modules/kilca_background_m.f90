module kilca_background_m
    use iso_fortran_env, only: real64, int32
    use kilca_types_m
    use kilca_settings_m
    use kilca_spline_m
    implicit none
    
    private
    
    ! Background data structure
    type, public :: background_t
        type(settings_t), pointer :: sd => null()  ! Pointer to settings data
        
        ! Profiles
        integer :: Nprofiles                        ! Number of basic input profiles
        character(len=256), allocatable :: profile_names(:)  ! Names of profiles
        character(len=256), pointer :: path2background => null()  ! Output path
        
        ! Grid data
        integer :: dimx                             ! Dimension of profile x grid
        real(real64), pointer :: x(:) => null()     ! x grid for profiles
        
        integer :: dimy                             ! Number of all background profiles
        real(real64), pointer :: y(:) => null()     ! y grids (1D array)
        
        ! Spline data
        integer :: ind                              ! Search index
        real(real64), pointer :: C(:) => null()     ! Spline coefficients
        real(real64), pointer :: R(:) => null()     ! Spline values at point
        type(spline_data_t) :: spline             ! Spline structure
        
        ! Profile indices for accessing y, C, R arrays
        ! Basic profiles
        integer :: i_q                              ! Safety factor q
        integer :: i_n                              ! Density
        integer, allocatable :: i_T(:)              ! Temperature (ions, electrons)
        integer, allocatable :: i_Vth(:)            ! Thermal velocity (i,e,total)
        integer, allocatable :: i_Vz(:)             ! Axial velocity (i,e,total)
        integer :: i_Er                             ! Radial electric field
        
        ! From equilibrium
        integer :: i_hth                            ! h_theta metric
        integer :: i_hz                             ! h_z metric
        integer :: i_Bth                            ! B_theta field
        integer :: i_Bz                             ! B_z field
        integer :: i_B                              ! Total B field
        
        integer, allocatable :: i_J0th(:)           ! Current density theta (i,e,total)
        integer, allocatable :: i_J0z(:)            ! Current density z (i,e,total)
        
        ! From f0 parameter search
        integer :: i_dPhi0                          ! Potential correction
        integer, allocatable :: i_n_p(:)            ! Density parameter (i,e)
        integer, allocatable :: i_Vp_p(:)           ! Poloidal velocity parameter (i,e)
        integer, allocatable :: i_Vt_p(:)           ! Toroidal velocity parameter (i,e)
        
        ! Additional
        integer, allocatable :: i_nu(:)             ! Collision frequency (i,e)
        
        integer :: flag_dPhi0_calc                  ! Flag to recalculate dPhi0
        
    end type background_t
    
    ! Public procedures
    public :: background_create
    public :: background_destroy
    public :: background_set_profiles_indices
    public :: background_load_input_profiles
    public :: background_spline_profiles
    public :: background_find_f0_parameters
    public :: background_calculate_equilibrium
    public :: background_save
    public :: background_interp_basic_profiles
    
contains

    !---------------------------------------------------------------------------
    ! Create and initialize background structure
    !---------------------------------------------------------------------------
    subroutine background_create(bg, settings, ierr)
        type(background_t), intent(out) :: bg
        type(settings_t), intent(in), target :: settings
        integer, intent(out) :: ierr
        
        integer :: nion
        
        ierr = 0
        
        ! Link to settings
        bg%sd => settings
        
        ! Get number of ion species (default to 1 for now)
        nion = 1  ! TODO: Get from settings when particle info is added
        
        ! Initialize scalars
        bg%Nprofiles = 7  ! TODO: Get from settings when profile info is added
        bg%dimx = 0
        bg%dimy = 0
        bg%ind = 1
        bg%flag_dPhi0_calc = 0
        
        ! Initialize indices to -1 (unset)
        bg%i_q = -1
        bg%i_n = -1
        bg%i_Er = -1
        bg%i_hth = -1
        bg%i_hz = -1
        bg%i_Bth = -1
        bg%i_Bz = -1
        bg%i_B = -1
        bg%i_dPhi0 = -1
        
        ! Allocate profile name array
        if (bg%Nprofiles > 0) then
            allocate(bg%profile_names(bg%Nprofiles))
            bg%profile_names = ""
        end if
        
        ! Allocate index arrays based on number of species
        ! Temperature: ions + electrons = nion + 1
        allocate(bg%i_T(0:nion))
        bg%i_T = -1
        
        ! Velocities: ions + electrons + total = nion + 2
        allocate(bg%i_Vth(0:nion+1))
        allocate(bg%i_Vz(0:nion+1))
        allocate(bg%i_J0th(0:nion+1))
        allocate(bg%i_J0z(0:nion+1))
        bg%i_Vth = -1
        bg%i_Vz = -1
        bg%i_J0th = -1
        bg%i_J0z = -1
        
        ! Parameters: ions + electrons = nion + 1
        allocate(bg%i_n_p(0:nion))
        allocate(bg%i_Vp_p(0:nion))
        allocate(bg%i_Vt_p(0:nion))
        bg%i_n_p = -1
        bg%i_Vp_p = -1
        bg%i_Vt_p = -1
        
        ! Collision frequencies: ions + electrons = nion + 1
        allocate(bg%i_nu(0:nion))
        bg%i_nu = -1
        
        ! Set output path (placeholder)
        ! bg%path2background => settings%output%path2output  ! TODO: Add when output path is available
        
    end subroutine background_create
    
    !---------------------------------------------------------------------------
    ! Destroy background structure and free memory
    !---------------------------------------------------------------------------
    subroutine background_destroy(bg)
        type(background_t), intent(inout) :: bg
        
        ! Deallocate arrays
        if (allocated(bg%profile_names)) deallocate(bg%profile_names)
        if (allocated(bg%i_T)) deallocate(bg%i_T)
        if (allocated(bg%i_Vth)) deallocate(bg%i_Vth)
        if (allocated(bg%i_Vz)) deallocate(bg%i_Vz)
        if (allocated(bg%i_J0th)) deallocate(bg%i_J0th)
        if (allocated(bg%i_J0z)) deallocate(bg%i_J0z)
        if (allocated(bg%i_n_p)) deallocate(bg%i_n_p)
        if (allocated(bg%i_Vp_p)) deallocate(bg%i_Vp_p)
        if (allocated(bg%i_Vt_p)) deallocate(bg%i_Vt_p)
        if (allocated(bg%i_nu)) deallocate(bg%i_nu)
        
        ! Deallocate pointers
        if (associated(bg%x)) deallocate(bg%x)
        if (associated(bg%y)) deallocate(bg%y)
        if (associated(bg%C)) deallocate(bg%C)
        if (associated(bg%R)) deallocate(bg%R)
        
        ! Destroy spline
        call spline_destroy(bg%spline)
        
        ! Nullify pointers
        nullify(bg%sd)
        nullify(bg%path2background)
        nullify(bg%x)
        nullify(bg%y)
        nullify(bg%C)
        nullify(bg%R)
        
    end subroutine background_destroy
    
    !---------------------------------------------------------------------------
    ! Set profile indices for array access
    !---------------------------------------------------------------------------
    subroutine background_set_profiles_indices(bg)
        type(background_t), intent(inout) :: bg
        
        integer :: idx, i, nion
        
        nion = 1  ! TODO: Get from settings when particle info is added
        idx = 0
        
        ! Basic profiles
        bg%i_q = idx
        idx = idx + 1
        
        bg%i_n = idx
        idx = idx + 1
        
        ! Temperature profiles (ions + electrons)
        do i = 0, nion
            bg%i_T(i) = idx
            idx = idx + 1
        end do
        
        ! Velocity profiles
        do i = 0, nion + 1
            bg%i_Vth(i) = idx
            idx = idx + 1
        end do
        
        do i = 0, nion + 1
            bg%i_Vz(i) = idx
            idx = idx + 1
        end do
        
        bg%i_Er = idx
        idx = idx + 1
        
        ! Equilibrium fields
        bg%i_hth = idx
        idx = idx + 1
        
        bg%i_hz = idx
        idx = idx + 1
        
        bg%i_Bth = idx
        idx = idx + 1
        
        bg%i_Bz = idx
        idx = idx + 1
        
        bg%i_B = idx
        idx = idx + 1
        
        ! Current densities
        do i = 0, nion + 1
            bg%i_J0th(i) = idx
            idx = idx + 1
        end do
        
        do i = 0, nion + 1
            bg%i_J0z(i) = idx
            idx = idx + 1
        end do
        
        ! F0 parameters
        bg%i_dPhi0 = idx
        idx = idx + 1
        
        do i = 0, nion
            bg%i_n_p(i) = idx
            idx = idx + 1
        end do
        
        do i = 0, nion
            bg%i_Vp_p(i) = idx
            idx = idx + 1
        end do
        
        do i = 0, nion
            bg%i_Vt_p(i) = idx
            idx = idx + 1
        end do
        
        ! Collision frequencies
        do i = 0, nion
            bg%i_nu(i) = idx
            idx = idx + 1
        end do
        
        ! Total number of profiles
        bg%dimy = idx
        
    end subroutine background_set_profiles_indices
    
    !---------------------------------------------------------------------------
    ! Load input background profiles (placeholder)
    !---------------------------------------------------------------------------
    subroutine background_load_input_profiles(bg, ierr)
        type(background_t), intent(inout) :: bg
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! This would load profiles from files
        ! For now, just set up the structure
        
        if (bg%dimx == 0) bg%dimx = 100  ! Default grid size
        
        ! Allocate arrays
        if (.not. associated(bg%x)) allocate(bg%x(bg%dimx))
        if (.not. associated(bg%y)) allocate(bg%y(bg%dimx * bg%dimy))
        
        ! Initialize to zero
        bg%x = 0.0_real64
        bg%y = 0.0_real64
        
    end subroutine background_load_input_profiles
    
    !---------------------------------------------------------------------------
    ! Create splines for profiles
    !---------------------------------------------------------------------------
    subroutine background_spline_profiles(bg, ierr)
        type(background_t), intent(inout) :: bg
        integer, intent(out) :: ierr
        
        integer :: i, j, k, offset
        
        ierr = 0
        
        ! Allocate coefficient array
        if (.not. associated(bg%C)) then
            allocate(bg%C(bg%dimx * bg%dimy * 4))  ! Cubic spline needs 4 coefficients
        end if
        
        ! Allocate result array
        if (.not. associated(bg%R)) then
            allocate(bg%R(bg%dimy * 2))  ! Value and first derivative
        end if
        
        ! Create spline structure
        call spline_create(bg%spline, 3, SPLINE_NATURAL, bg%dimx, bg%x, ierr)
        if (ierr /= 0) return
        
        ! Calculate spline coefficients for each profile
        do i = 1, bg%dimy
            offset = (i-1) * bg%dimx + 1
            call spline_calc_coefficients(bg%spline, bg%y(offset:offset+bg%dimx-1), ierr)
            if (ierr /= 0) return
            
            ! Store coefficients from spline structure
            ! For cubic splines (N=3), we have 4 coefficients per interval
            ! The spline%C matrix contains these coefficients
            if (associated(bg%spline%C)) then
                ! Copy coefficients from spline structure
                ! spline%C is (dimx, 4) for cubic splines
                do j = 1, bg%dimx
                    do k = 1, 4
                        bg%C((i-1)*bg%dimx*4 + (j-1)*4 + k) = bg%spline%C(j, k)
                    end do
                end do
            else
                ! Fallback if spline coefficients not available
                bg%C((i-1)*bg%dimx*4+1:i*bg%dimx*4) = 0.0_real64
                ierr = -1
                print *, "Warning: Spline coefficients not calculated properly for profile", i
            end if
        end do
        
    end subroutine background_spline_profiles
    
    !---------------------------------------------------------------------------
    ! Find f0 parameters (placeholder)
    !---------------------------------------------------------------------------
    subroutine background_find_f0_parameters(bg, ierr)
        type(background_t), intent(inout) :: bg
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! This would implement the f0 parameter search
        ! For now, just set default values
        
        if (bg%i_dPhi0 >= 0 .and. bg%i_dPhi0 < bg%dimy) then
            bg%y(bg%i_dPhi0 * bg%dimx + 1:(bg%i_dPhi0 + 1) * bg%dimx) = 0.0_real64
        end if
        
    end subroutine background_find_f0_parameters
    
    !---------------------------------------------------------------------------
    ! Calculate equilibrium fields (placeholder)
    !---------------------------------------------------------------------------
    subroutine background_calculate_equilibrium(bg, ierr)
        type(background_t), intent(inout) :: bg
        integer, intent(out) :: ierr
        
        integer :: i
        real(real64) :: B0, rtor
        
        ierr = 0
        
        ! Get machine parameters from background settings
        B0 = bg%sd%background_settings%B0
        rtor = bg%sd%background_settings%rtor
        
        ! Simple cylindrical equilibrium
        do i = 1, bg%dimx
            ! Metric coefficients (cylindrical)
            if (bg%i_hth >= 0) bg%y(bg%i_hth * bg%dimx + i) = bg%x(i)  ! h_theta = r
            if (bg%i_hz >= 0) bg%y(bg%i_hz * bg%dimx + i) = 1.0_real64  ! h_z = 1
            
            ! Magnetic field (simple model)
            if (bg%i_Bth >= 0) bg%y(bg%i_Bth * bg%dimx + i) = 0.0_real64  ! No poloidal field
            if (bg%i_Bz >= 0) bg%y(bg%i_Bz * bg%dimx + i) = B0  ! Constant toroidal field
            if (bg%i_B >= 0) bg%y(bg%i_B * bg%dimx + i) = B0  ! Total field
        end do
        
    end subroutine background_calculate_equilibrium
    
    !---------------------------------------------------------------------------
    ! Save background data (placeholder)
    !---------------------------------------------------------------------------
    subroutine background_save(bg, ierr)
        type(background_t), intent(in) :: bg  ! NOTE: bg kept for interface compatibility
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! This would save background data to files
        ! For now, just return success
        
    end subroutine background_save
    
    !---------------------------------------------------------------------------
    ! Interpolate basic background profiles at given radial position
    !---------------------------------------------------------------------------
    subroutine background_interp_basic_profiles(bg, r, q, n, Ti, Te, Vth, Vz, dPhi0)
        type(background_t), intent(in) :: bg  ! NOTE: bg kept for interface compatibility
        real(real64), intent(in) :: r  ! NOTE: r kept for interface compatibility
        real(real64), intent(out) :: q, n, Ti, Te, Vth, Vz, dPhi0
        
        ! For now, return dummy values
        ! In full implementation, would use spline interpolation
        
        q = 2.0_real64
        n = 1.0e19_real64
        Ti = 1000.0_real64
        Te = 1000.0_real64
        Vth = 0.0_real64
        Vz = 0.0_real64
        dPhi0 = 0.0_real64
        
    end subroutine background_interp_basic_profiles

end module kilca_background_m