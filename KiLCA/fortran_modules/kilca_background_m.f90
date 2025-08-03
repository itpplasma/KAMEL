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
    ! Initialize default realistic plasma profiles
    !---------------------------------------------------------------------------
    subroutine initialize_default_profiles(bg)
        type(background_t), intent(inout) :: bg
        
        integer :: i
        real(real64) :: rho, q_center, q_edge, n_center, n_edge
        real(real64) :: T_center, Te_center, Te_edge, Ti_center, Ti_edge
        
        ! Default parameters for typical tokamak profiles
        q_center = 1.0_real64       ! Safety factor at center
        q_edge = 3.0_real64         ! Safety factor at edge
        n_center = 2.0e19_real64    ! Density at center [m^-3]
        n_edge = 0.1e19_real64      ! Density at edge [m^-3]
        Te_center = 3000.0_real64   ! Electron temperature at center [eV]
        Te_edge = 100.0_real64      ! Electron temperature at edge [eV]
        Ti_center = 2000.0_real64   ! Ion temperature at center [eV]
        Ti_edge = 100.0_real64      ! Ion temperature at edge [eV]
        
        ! Initialize all profiles to zero first
        bg%y = 0.0_real64
        
        ! Validate grid is properly set up
        if (.not. associated(bg%x) .or. bg%dimx <= 0) then
            print *, "Warning: Radial grid not properly initialized"
            return
        end if
        
        ! Generate realistic radial profiles with optimized loop
        do i = 1, bg%dimx
            rho = bg%x(i)  ! Normalized radius 0-1
            
            ! Safety factor profile (monotonic increasing)
            if (bg%i_q >= 0) then
                bg%y(bg%i_q * bg%dimx + i) = q_center + (q_edge - q_center) * rho**2
            end if
            
            ! Density profile (peaked)
            if (bg%i_n >= 0) then
                bg%y(bg%i_n * bg%dimx + i) = n_center * (1.0_real64 - rho**2)**2 + n_edge
            end if
            
            ! Electron temperature profile (peaked)
            if (bg%i_T(0) >= 0) then  ! Electrons are species 0
                bg%y(bg%i_T(0) * bg%dimx + i) = Te_center * (1.0_real64 - rho**2)**1.5_real64 + Te_edge
            end if
            
            ! Ion temperature profile (peaked)
            if (bg%i_T(1) >= 0) then  ! Ions are species 1
                bg%y(bg%i_T(1) * bg%dimx + i) = Ti_center * (1.0_real64 - rho**2)**1.5_real64 + Ti_edge
            end if
            
            ! Thermal velocities (calculated from temperature)
            if (bg%i_Vth(0) >= 0) then  ! Electron thermal velocity
                T_center = bg%y(bg%i_T(0) * bg%dimx + i) * 1.602e-19_real64  ! Convert eV to J
                bg%y(bg%i_Vth(0) * bg%dimx + i) = sqrt(2.0_real64 * T_center / 9.109e-31_real64)  ! v_th = sqrt(2T/m)
            end if
            
            if (bg%i_Vth(1) >= 0) then  ! Ion thermal velocity
                T_center = bg%y(bg%i_T(1) * bg%dimx + i) * 1.602e-19_real64  ! Convert eV to J
                bg%y(bg%i_Vth(1) * bg%dimx + i) = sqrt(2.0_real64 * T_center / 1.673e-27_real64)  ! v_th = sqrt(2T/m)
            end if
            
            ! Axial velocities (small rotation)
            if (bg%i_Vz(0) >= 0) then  ! Electron axial velocity
                bg%y(bg%i_Vz(0) * bg%dimx + i) = 1.0e4_real64 * (1.0_real64 - rho)  ! Simple rotation profile
            end if
            
            if (bg%i_Vz(1) >= 0) then  ! Ion axial velocity
                bg%y(bg%i_Vz(1) * bg%dimx + i) = 1.0e4_real64 * (1.0_real64 - rho)  ! Simple rotation profile
            end if
            
            ! Electric field (small radial field)
            if (bg%i_Er >= 0) then
                bg%y(bg%i_Er * bg%dimx + i) = -1000.0_real64 * rho  ! Linear radial electric field [V/m]
            end if
        end do
        
    end subroutine initialize_default_profiles
    
    !---------------------------------------------------------------------------
    ! Load input background profiles (placeholder)
    !---------------------------------------------------------------------------
    subroutine background_load_input_profiles(bg, ierr)
        type(background_t), intent(inout) :: bg
        integer, intent(out) :: ierr
        
        integer :: i
        real(real64) :: rho_max, rho_step
        logical :: grid_already_exists
        
        ierr = 0
        
        ! This would load profiles from files
        ! For now, generate a normalized radial grid from 0 to 1
        
        ! Check if grid is already set up (e.g., by tests)
        grid_already_exists = (bg%dimx > 0) .and. associated(bg%x) .and. associated(bg%y)
        
        if (.not. grid_already_exists) then
            ! Set up new grid
            if (bg%dimx == 0) bg%dimx = 100  ! Default grid size
            
            ! Allocate arrays
            if (.not. associated(bg%x)) allocate(bg%x(bg%dimx))
            if (.not. associated(bg%y)) allocate(bg%y(bg%dimx * bg%dimy))
            
            ! Generate normalized radial grid: rho = r/a, from 0 to 1
            rho_max = 1.0_real64
            rho_step = rho_max / real(bg%dimx - 1, real64)
            
            do i = 1, bg%dimx
                bg%x(i) = real(i - 1, real64) * rho_step
            end do
        end if
        
        ! Initialize profile data with realistic plasma values
        call initialize_default_profiles(bg)
        
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
        
        integer :: i
        real(real64) :: phi0_center, phi0_edge, rho
        
        ierr = 0
        
        ! This would implement the f0 parameter search
        ! For now, set up reasonable default values for F0 parameters
        
        phi0_center = 0.0_real64      ! Potential correction at center [V]
        phi0_edge = -100.0_real64     ! Potential correction at edge [V]
        
        ! Set dPhi0 profile
        if (bg%i_dPhi0 >= 0 .and. bg%i_dPhi0 < bg%dimy) then
            do i = 1, bg%dimx
                rho = bg%x(i)
                bg%y(bg%i_dPhi0 * bg%dimx + i) = phi0_center + (phi0_edge - phi0_center) * rho**2
            end do
        end if
        
        ! Set density parameters (default to 1.0 for now)
        if (allocated(bg%i_n_p)) then
            do i = 0, size(bg%i_n_p) - 1
                if (bg%i_n_p(i) >= 0) then
                    bg%y(bg%i_n_p(i) * bg%dimx + 1:(bg%i_n_p(i) + 1) * bg%dimx) = 1.0_real64
                end if
            end do
        end if
        
        ! Set velocity parameters (default to small values)
        if (allocated(bg%i_Vp_p)) then
            do i = 0, size(bg%i_Vp_p) - 1
                if (bg%i_Vp_p(i) >= 0) then
                    bg%y(bg%i_Vp_p(i) * bg%dimx + 1:(bg%i_Vp_p(i) + 1) * bg%dimx) = 0.1_real64
                end if
            end do
        end if
        
        if (allocated(bg%i_Vt_p)) then
            do i = 0, size(bg%i_Vt_p) - 1
                if (bg%i_Vt_p(i) >= 0) then
                    bg%y(bg%i_Vt_p(i) * bg%dimx + 1:(bg%i_Vt_p(i) + 1) * bg%dimx) = 0.1_real64
                end if
            end do
        end if
        
        ! Set flag to indicate F0 parameters have been calculated
        bg%flag_dPhi0_calc = 1
        
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
    ! Find interpolation indices and weight for given radius
    ! Uses optimized linear search for monotonic grid
    !---------------------------------------------------------------------------
    subroutine find_interpolation_indices(bg, r, i_low, i_high, weight)
        type(background_t), intent(in) :: bg
        real(real64), intent(in) :: r
        integer, intent(out) :: i_low, i_high
        real(real64), intent(out) :: weight
        
        integer :: i
        real(real64), parameter :: eps = 1.0e-14_real64
        
        ! Handle edge cases first for efficiency
        if (r <= bg%x(1)) then
            i_low = 1; i_high = 1; weight = 0.0_real64
            return
        else if (r >= bg%x(bg%dimx)) then
            i_low = bg%dimx; i_high = bg%dimx; weight = 0.0_real64
            return
        end if
        
        ! Linear search for interval (could use binary search for large grids)
        i_low = 1
        do i = 1, bg%dimx - 1
            if (r >= bg%x(i) .and. r <= bg%x(i + 1)) then
                i_low = i
                i_high = i + 1
                exit
            end if
        end do
        
        ! Calculate interpolation weight with numerical stability
        if (abs(bg%x(i_high) - bg%x(i_low)) > eps) then
            weight = (r - bg%x(i_low)) / (bg%x(i_high) - bg%x(i_low))
        else
            weight = 0.0_real64  ! Avoid division by zero
        end if
        
    end subroutine find_interpolation_indices
    
    !---------------------------------------------------------------------------
    ! Perform linear interpolation between two values
    !---------------------------------------------------------------------------
    function linear_interpolate(val_low, val_high, weight) result(val)
        real(real64), intent(in) :: val_low, val_high, weight
        real(real64) :: val
        
        val = val_low + weight * (val_high - val_low)
        
    end function linear_interpolate
    
    !---------------------------------------------------------------------------
    ! Interpolate basic background profiles at given radial position
    ! 
    ! Uses linear interpolation between grid points for plasma parameters.
    ! Input radius is clamped to [0,1] range. Returns physical defaults
    ! if background is not properly initialized.
    !
    ! @param bg      Background data structure
    ! @param r       Normalized radial coordinate [0,1]
    ! @param q       Safety factor (output)
    ! @param n       Density [m^-3] (output)
    ! @param Ti      Ion temperature [eV] (output)
    ! @param Te      Electron temperature [eV] (output)
    ! @param Vth     Thermal velocity [m/s] (output)
    ! @param Vz      Axial velocity [m/s] (output)
    ! @param dPhi0   Potential correction [V] (output)
    !---------------------------------------------------------------------------
    subroutine background_interp_basic_profiles(bg, r, q, n, Ti, Te, Vth, Vz, dPhi0)
        type(background_t), intent(in) :: bg
        real(real64), intent(in) :: r
        real(real64), intent(out) :: q, n, Ti, Te, Vth, Vz, dPhi0
        
        integer :: i_low, i_high
        real(real64) :: r_clamped, weight
        
        ! Validate input
        if (.not. associated(bg%x) .or. .not. associated(bg%y) .or. bg%dimx <= 0) then
            ! Return default values if background not properly initialized
            q = 2.0_real64; n = 1.0e19_real64; Ti = 1000.0_real64; Te = 1000.0_real64
            Vth = 0.0_real64; Vz = 0.0_real64; dPhi0 = 0.0_real64
            return
        end if
        
        ! Clamp radius to [0, 1] range for robustness
        r_clamped = max(0.0_real64, min(1.0_real64, r))
        
        ! Find interpolation indices
        call find_interpolation_indices(bg, r_clamped, i_low, i_high, weight)
        
        ! Interpolate safety factor q
        if (bg%i_q >= 0) then
            q = linear_interpolate(bg%y(bg%i_q * bg%dimx + i_low), &
                                  bg%y(bg%i_q * bg%dimx + i_high), weight)
        else
            q = 2.0_real64  ! Default fallback
        end if
        
        ! Interpolate density n
        if (bg%i_n >= 0) then
            n = linear_interpolate(bg%y(bg%i_n * bg%dimx + i_low), &
                                  bg%y(bg%i_n * bg%dimx + i_high), weight)
        else
            n = 1.0e19_real64  ! Default fallback
        end if
        
        ! Interpolate ion temperature Ti (species 1)
        if (allocated(bg%i_T) .and. size(bg%i_T) > 1 .and. bg%i_T(1) >= 0) then
            Ti = linear_interpolate(bg%y(bg%i_T(1) * bg%dimx + i_low), &
                                   bg%y(bg%i_T(1) * bg%dimx + i_high), weight)
        else
            Ti = 1000.0_real64  ! Default fallback
        end if
        
        ! Interpolate electron temperature Te (species 0)
        if (allocated(bg%i_T) .and. size(bg%i_T) > 0 .and. bg%i_T(0) >= 0) then
            Te = linear_interpolate(bg%y(bg%i_T(0) * bg%dimx + i_low), &
                                   bg%y(bg%i_T(0) * bg%dimx + i_high), weight)
        else
            Te = 1000.0_real64  ! Default fallback
        end if
        
        ! Interpolate thermal velocity Vth (average of species)
        if (allocated(bg%i_Vth) .and. size(bg%i_Vth) > 1 .and. bg%i_Vth(1) >= 0) then
            Vth = linear_interpolate(bg%y(bg%i_Vth(1) * bg%dimx + i_low), &
                                    bg%y(bg%i_Vth(1) * bg%dimx + i_high), weight)
        else
            Vth = 0.0_real64  ! Default fallback
        end if
        
        ! Interpolate axial velocity Vz (average of species)
        if (allocated(bg%i_Vz) .and. size(bg%i_Vz) > 1 .and. bg%i_Vz(1) >= 0) then
            Vz = linear_interpolate(bg%y(bg%i_Vz(1) * bg%dimx + i_low), &
                                   bg%y(bg%i_Vz(1) * bg%dimx + i_high), weight)
        else
            Vz = 0.0_real64  ! Default fallback
        end if
        
        ! Interpolate potential correction dPhi0
        if (bg%i_dPhi0 >= 0) then
            dPhi0 = linear_interpolate(bg%y(bg%i_dPhi0 * bg%dimx + i_low), &
                                      bg%y(bg%i_dPhi0 * bg%dimx + i_high), weight)
        else
            dPhi0 = 0.0_real64  ! Default fallback
        end if
        
    end subroutine background_interp_basic_profiles

end module kilca_background_m