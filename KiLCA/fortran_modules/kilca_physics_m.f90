module kilca_physics_m
    !---------------------------------------------------------------------------
    ! KiLCA Plasma Physics Models Module
    !
    ! This module provides comprehensive plasma physics models including
    ! background plasma profiles, physics zones (FLRE, vacuum, IMHD),
    ! dispersion relations, and zone management systems.
    !
    ! Key features:
    ! 1. Background plasma profile management and interpolation
    ! 2. Equilibrium calculation and consistency checks
    ! 3. FLRE (Finite Larmor Radius Effects) zone physics
    ! 4. Vacuum zone electromagnetic wave propagation
    ! 5. IMHD (Ideal MagnetoHydroDynamic) zone physics
    ! 6. Dispersion relations and wave physics
    ! 7. Physics zone management and coupling
    ! 8. Lab frame transformations and Doppler shifts
    ! 9. F0 distribution function moments
    ! 10. Collision frequency calculations
    ! 11. Magnetic field geometry calculations
    ! 12. Physics spline interpolation
    ! 13. Zone coupling and continuity conditions
    ! 14. Physics parameter validation
    !
    ! Author: Claude (Anthropic)
    ! Date: 2024
    !
    ! Translation from: KiLCA background, zone, and physics systems
    !---------------------------------------------------------------------------
    
    use iso_fortran_env, only: real64, int32, error_unit
    use kilca_types_m
    implicit none
    
    private
    
    ! Background plasma profile data structure
    type, public :: background_profiles_t
        ! Grid data
        integer :: n_points                        ! Number of grid points
        real(real64), allocatable :: x_grid(:)     ! Radial grid (normalized)
        
        ! Basic plasma profiles
        real(real64), allocatable :: q_profile(:)   ! Safety factor profile
        real(real64), allocatable :: n_profile(:)   ! Density profile (1/m^3)
        real(real64), allocatable :: Ti_profile(:)  ! Ion temperature profile (eV)
        real(real64), allocatable :: Te_profile(:)  ! Electron temperature profile (eV)
        real(real64), allocatable :: Vth_profile(:) ! Thermal velocity profile (m/s)
        real(real64), allocatable :: Vz_profile(:)  ! Toroidal velocity profile (m/s)
        real(real64), allocatable :: Er_profile(:)  ! Radial electric field profile (V/m)
        
        ! Equilibrium data
        real(real64), allocatable :: Bt_profile(:)  ! Toroidal magnetic field (T)
        real(real64), allocatable :: Bp_profile(:)  ! Poloidal magnetic field (T)
        real(real64), allocatable :: B_profile(:)   ! Total magnetic field (T)
        
        ! Spline coefficients for interpolation
        real(real64), allocatable :: q_spline(:)    ! Safety factor spline coeffs
        real(real64), allocatable :: n_spline(:)    ! Density spline coeffs
        real(real64), allocatable :: Ti_spline(:)   ! Ion temperature spline coeffs
        real(real64), allocatable :: Te_spline(:)   ! Electron temperature spline coeffs
        
        ! Profile names and paths
        character(len=256) :: path_to_background    ! Path to background directory
        integer :: n_profiles                       ! Number of basic profiles
        character(len=64), allocatable :: profile_names(:) ! Names of profiles
    end type background_profiles_t
    
    ! Equilibrium calculation data
    type, public :: equilibrium_t
        ! Geometry parameters
        real(real64) :: R_major                     ! Major radius (m)
        real(real64) :: r_minor                     ! Minor radius (m)
        real(real64) :: aspect_ratio                ! R/a
        
        ! Magnetic field components
        real(real64) :: B0                          ! Reference magnetic field (T)
        real(real64) :: Bt0                         ! Toroidal field on axis (T)
        real(real64) :: I_plasma                    ! Plasma current (A)
        
        ! Pressure and beta
        real(real64) :: pressure_central            ! Central pressure (Pa)
        real(real64) :: beta_pol                    ! Poloidal beta
        real(real64) :: beta_tor                    ! Toroidal beta
        
        ! Force balance check
        logical :: is_equilibrium_valid             ! Equilibrium validity flag
        real(real64) :: force_balance_error         ! Force balance error metric
    end type equilibrium_t
    
    ! FLRE zone physics data
    type, public :: flre_zone_t
        ! Zone identification
        integer :: zone_id                          ! Zone ID number
        real(real64) :: r_inner, r_outer            ! Zone boundaries
        
        ! Plasma parameters
        real(real64) :: omega_ci                    ! Ion cyclotron frequency (rad/s)
        real(real64) :: omega_ce                    ! Electron cyclotron frequency (rad/s)
        real(real64) :: rho_i                       ! Ion Larmor radius (m)
        real(real64) :: rho_e                       ! Electron Larmor radius (m)
        
        ! FLRE expansion parameters
        real(real64) :: k_perp                      ! Perpendicular wavenumber (1/m)
        real(real64) :: flre_order_i                ! Ion FLRE parameter
        real(real64) :: flre_order_e                ! Electron FLRE parameter
        integer :: max_flre_order                   ! Maximum FLRE expansion order
        
        ! Dielectric tensor elements
        complex(real64) :: epsilon_perp             ! Perpendicular dielectric
        complex(real64) :: epsilon_par              ! Parallel dielectric
        complex(real64) :: epsilon_xy               ! Off-diagonal element
    end type flre_zone_t
    
    ! Vacuum zone physics data
    type, public :: vacuum_zone_t
        ! Zone identification
        integer :: zone_id                          ! Zone ID number
        real(real64) :: r_inner, r_outer            ! Zone boundaries
        
        ! Vacuum properties
        real(real64) :: c_light                     ! Speed of light (m/s)
        real(real64) :: mu_0                        ! Vacuum permeability
        real(real64) :: epsilon_0                   ! Vacuum permittivity
        real(real64) :: wave_impedance              ! Wave impedance (Ohm)
        
        ! Wave parameters
        complex(real64) :: k_vac                    ! Vacuum wave number
        complex(real64) :: omega                    ! Angular frequency (rad/s)
        
        ! Field relationships
        real(real64) :: E_to_B_ratio                ! E/B = c in vacuum
    end type vacuum_zone_t
    
    ! IMHD zone physics data
    type, public :: imhd_zone_t
        ! Zone identification
        integer :: zone_id                          ! Zone ID number
        real(real64) :: r_inner, r_outer            ! Zone boundaries
        
        ! MHD parameters
        real(real64) :: rho_mass                    ! Mass density (kg/m^3)
        real(real64) :: pressure                    ! Pressure (Pa)
        real(real64) :: gamma_adiabatic             ! Adiabatic index
        
        ! Characteristic speeds
        real(real64) :: v_alfven                    ! Alfven speed (m/s)
        real(real64) :: c_sound                     ! Sound speed (m/s)
        real(real64) :: beta_plasma                 ! Plasma beta
        
        ! Flow velocities
        real(real64) :: vx, vy, vz                  ! Velocity components (m/s)
        
        ! MHD constraints
        real(real64) :: divB                        ! div(B) constraint
        real(real64) :: curlB_magnitude             ! |curl(B)| for current
    end type imhd_zone_t
    
    ! Dispersion relation data
    type, public :: dispersion_t
        ! Plasma frequencies
        real(real64) :: omega_pe                    ! Electron plasma frequency (rad/s)
        real(real64) :: omega_pi                    ! Ion plasma frequency (rad/s)
        real(real64) :: omega_ce                    ! Electron cyclotron frequency (rad/s)
        real(real64) :: omega_ci                    ! Ion cyclotron frequency (rad/s)
        
        ! Wave parameters
        complex(real64) :: omega                    ! Complex frequency (rad/s)
        complex(real64) :: k_par                    ! Parallel wavenumber (1/m)
        complex(real64) :: k_perp                   ! Perpendicular wavenumber (1/m)
        
        ! Dielectric tensor
        complex(real64) :: epsilon_xx               ! Dielectric tensor elements
        complex(real64) :: epsilon_yy
        complex(real64) :: epsilon_zz
        complex(real64) :: epsilon_xy
        
        ! Refractive indices
        complex(real64) :: n_index_o                ! Ordinary wave index
        complex(real64) :: n_index_x                ! Extraordinary wave index
        
        ! Dispersion determinant
        complex(real64) :: det_disp                 ! Dispersion relation determinant
    end type dispersion_t
    
    ! Physics zone management
    type, public :: zone_manager_t
        ! Zone configuration
        integer :: n_zones                          ! Number of zones
        integer, allocatable :: zone_types(:)       ! Zone type array (1=FLRE, 2=Vacuum, 3=IMHD)
        real(real64), allocatable :: zone_boundaries(:) ! Zone boundary positions
        
        ! Zone data storage
        type(flre_zone_t), allocatable :: flre_zones(:)
        type(vacuum_zone_t), allocatable :: vacuum_zones(:)
        type(imhd_zone_t), allocatable :: imhd_zones(:)
        
        ! Zone coupling data
        real(real64), allocatable :: coupling_strength(:) ! Inter-zone coupling
        logical :: zones_coupled                    ! Zone coupling status
    end type zone_manager_t
    
    ! Lab frame transformation data
    type, public :: lab_frame_transform_t
        ! Rotation parameters
        real(real64) :: omega_lab                   ! Laboratory frequency (rad/s)
        real(real64) :: omega_plasma                ! Plasma frame frequency (rad/s)
        real(real64) :: v_rot                       ! Rotation velocity (m/s)
        real(real64) :: R_major                     ! Major radius (m)
        
        ! Doppler shift
        real(real64) :: doppler_factor              ! Doppler shift factor
        real(real64) :: gamma_relativistic          ! Relativistic gamma (usually ~1)
        
        ! Transformation matrix
        real(real64) :: transformation_matrix(3,3)  ! Coordinate transformation
        
        ! Field transformations
        real(real64) :: E_lab(3), E_plasma(3)       ! Electric field vectors
        real(real64) :: B_lab(3), B_plasma(3)       ! Magnetic field vectors
    end type lab_frame_transform_t
    
    ! F0 distribution moments
    type, public :: f0_moments_t
        ! Plasma parameters
        real(real64) :: n_density                   ! Number density (1/m^3)
        real(real64) :: Ti                          ! Ion temperature (J)
        real(real64) :: Te                          ! Electron temperature (J)
        
        ! Thermal velocities
        real(real64) :: vth_i                       ! Ion thermal velocity (m/s)
        real(real64) :: vth_e                       ! Electron thermal velocity (m/s)
        
        ! Distribution moments
        real(real64) :: f0_moment_0                 ! 0th moment: density
        real(real64) :: f0_moment_1                 ! 1st moment: mean velocity
        real(real64) :: f0_moment_2                 ! 2nd moment: pressure/mass
        
        ! Derived quantities
        real(real64) :: pressure_i                  ! Ion pressure (Pa)
        real(real64) :: pressure_e                  ! Electron pressure (Pa)
        real(real64) :: energy_density              ! Total energy density (J/m^3)
        real(real64) :: beta_thermal                ! Thermal beta parameter
    end type f0_moments_t
    
    ! Collision frequencies
    type, public :: collision_freq_t
        ! Plasma parameters
        real(real64) :: n_density                   ! Number density (1/m^3)
        real(real64) :: Ti                          ! Ion temperature (eV)
        real(real64) :: Te                          ! Electron temperature (eV)
        real(real64) :: Z_eff                       ! Effective charge
        real(real64) :: ln_lambda                   ! Coulomb logarithm
        
        ! Collision frequencies
        real(real64) :: nu_ii                       ! Ion-ion collision frequency (1/s)
        real(real64) :: nu_ee                       ! Electron-electron collision frequency (1/s)
        real(real64) :: nu_ei                       ! Electron-ion collision frequency (1/s)
        real(real64) :: nu_ie                       ! Ion-electron collision frequency (1/s)
        
        ! Collision times
        real(real64) :: tau_ii                      ! Ion-ion collision time (s)
        real(real64) :: tau_ee                      ! Electron-electron collision time (s)
        real(real64) :: tau_ei                      ! Electron-ion collision time (s)
        real(real64) :: slowing_down_time           ! Slowing down time (s)
    end type collision_freq_t
    
    ! Magnetic field calculations
    type, public :: magnetic_field_t
        ! Geometry parameters
        real(real64) :: R_major                     ! Major radius (m)
        real(real64) :: r_minor                     ! Minor radius (m)
        real(real64) :: I_plasma                    ! Plasma current (A)
        real(real64) :: mu_0                        ! Vacuum permeability
        
        ! Field components
        real(real64) :: Bp                          ! Poloidal field (T)
        real(real64) :: Bt                          ! Toroidal field (T)
        real(real64) :: B_total                     ! Total field magnitude (T)
        real(real64) :: q_safety                    ! Safety factor
        
        ! Flux and gradients
        real(real64) :: psi_flux                    ! Poloidal flux
        real(real64) :: grad_B                      ! Magnetic field gradient
        real(real64) :: curvature_drift             ! Curvature drift frequency
        real(real64) :: magnetic_shear              ! Magnetic shear
        
        ! Field vector
        real(real64) :: B_field(3)                  ! (Br, Bp, Bt) components
        real(real64) :: grad_B_vec(3)               ! Gradient vector
    end type magnetic_field_t
    
    ! Public procedures
    public :: background_profiles_create
    public :: background_profiles_destroy
    public :: background_profiles_set_from_arrays
    public :: background_profiles_interpolate
    public :: background_profiles_transform_to_lab
    public :: equilibrium_create
    public :: equilibrium_destroy
    public :: equilibrium_calculate
    public :: equilibrium_check_force_balance
    public :: flre_zone_create
    public :: flre_zone_destroy
    public :: flre_zone_calculate_parameters
    public :: flre_zone_calculate_dielectric
    public :: vacuum_zone_create
    public :: vacuum_zone_destroy
    public :: vacuum_zone_calculate_parameters
    public :: vacuum_zone_check_field_relations
    public :: imhd_zone_create
    public :: imhd_zone_destroy
    public :: imhd_zone_calculate_speeds
    public :: imhd_zone_check_constraints
    public :: dispersion_create
    public :: dispersion_destroy
    public :: dispersion_calculate_frequencies
    public :: dispersion_calculate_dielectric
    public :: dispersion_calculate_indices
    public :: zone_manager_create
    public :: zone_manager_destroy
    public :: zone_manager_setup_zones
    public :: zone_manager_identify_zone
    public :: zone_manager_check_coupling
    public :: lab_frame_transform_create
    public :: lab_frame_transform_destroy
    public :: lab_frame_transform_calculate
    public :: lab_frame_transform_fields
    public :: f0_moments_create
    public :: f0_moments_destroy
    public :: f0_moments_calculate
    public :: f0_moments_calculate_pressure
    public :: collision_freq_create
    public :: collision_freq_destroy
    public :: collision_freq_calculate
    public :: collision_freq_calculate_times
    public :: magnetic_field_create
    public :: magnetic_field_destroy
    public :: magnetic_field_calculate
    public :: magnetic_field_calculate_gradients
    public :: physics_spline_interpolate
    public :: physics_check_zone_continuity
    public :: physics_validate_parameters
    
    ! Physical constants
    ! C_LIGHT already defined in kilca_types_m as c_light = 29979245800.0_dp (cm/s)
    real(real64), parameter :: MU_0 = 4.0e-7_real64 * PI     ! Vacuum permeability
    real(real64), parameter :: EPSILON_0 = 8.854187817e-12_real64  ! Vacuum permittivity (F/m)
    real(real64), parameter :: ELEM_CHARGE = 1.6022e-19_real64  ! Elementary charge (C)
    real(real64), parameter :: MASS_ELECTRON = 9.1094e-31_real64 ! Electron mass (kg)
    real(real64), parameter :: MASS_PROTON = 1.6726e-27_real64  ! Proton mass (kg)
    real(real64), parameter :: EV_TO_JOULE = 1.6022e-19_real64  ! eV to J conversion
    
contains

    !---------------------------------------------------------------------------
    ! Background profile management procedures
    !---------------------------------------------------------------------------
    subroutine background_profiles_create(profiles, n_points, ierr)
        type(background_profiles_t), intent(out) :: profiles
        integer, intent(in) :: n_points
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (n_points < 2) then
            ierr = -1
            return
        end if
        
        profiles%n_points = n_points
        profiles%n_profiles = 7  ! Standard set of profiles
        
        ! Allocate profile arrays
        allocate(profiles%x_grid(n_points), profiles%q_profile(n_points), &
                profiles%n_profile(n_points), profiles%Ti_profile(n_points), &
                profiles%Te_profile(n_points), profiles%Vth_profile(n_points), &
                profiles%Vz_profile(n_points), profiles%Er_profile(n_points), &
                profiles%Bt_profile(n_points), profiles%Bp_profile(n_points), &
                profiles%B_profile(n_points), stat=ierr)
        
        if (ierr /= 0) then
            ierr = -2
            return
        end if
        
        ! Allocate spline coefficient arrays (4*n_points for cubic splines)
        allocate(profiles%q_spline(4*n_points), profiles%n_spline(4*n_points), &
                profiles%Ti_spline(4*n_points), profiles%Te_spline(4*n_points), &
                stat=ierr)
        
        if (ierr /= 0) then
            ierr = -3
            return
        end if
        
        ! Allocate profile names
        allocate(profiles%profile_names(profiles%n_profiles), stat=ierr)
        if (ierr /= 0) then
            ierr = -4
            return
        end if
        
        ! Set default profile names
        profiles%profile_names = ['q       ', 'n       ', 'Ti      ', 'Te      ', &
                                 'Vth     ', 'Vz      ', 'Er      ']
        
        ! Initialize arrays to zero
        profiles%x_grid = 0.0_real64
        profiles%q_profile = 0.0_real64
        profiles%n_profile = 0.0_real64
        profiles%Ti_profile = 0.0_real64
        profiles%Te_profile = 0.0_real64
        profiles%Vth_profile = 0.0_real64
        profiles%Vz_profile = 0.0_real64
        profiles%Er_profile = 0.0_real64
        profiles%Bt_profile = 0.0_real64
        profiles%Bp_profile = 0.0_real64
        profiles%B_profile = 0.0_real64
        
        profiles%path_to_background = './background/'
        
    end subroutine background_profiles_create
    
    subroutine background_profiles_destroy(profiles, ierr)
        type(background_profiles_t), intent(inout) :: profiles
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (allocated(profiles%x_grid)) deallocate(profiles%x_grid)
        if (allocated(profiles%q_profile)) deallocate(profiles%q_profile)
        if (allocated(profiles%n_profile)) deallocate(profiles%n_profile)
        if (allocated(profiles%Ti_profile)) deallocate(profiles%Ti_profile)
        if (allocated(profiles%Te_profile)) deallocate(profiles%Te_profile)
        if (allocated(profiles%Vth_profile)) deallocate(profiles%Vth_profile)
        if (allocated(profiles%Vz_profile)) deallocate(profiles%Vz_profile)
        if (allocated(profiles%Er_profile)) deallocate(profiles%Er_profile)
        if (allocated(profiles%Bt_profile)) deallocate(profiles%Bt_profile)
        if (allocated(profiles%Bp_profile)) deallocate(profiles%Bp_profile)
        if (allocated(profiles%B_profile)) deallocate(profiles%B_profile)
        if (allocated(profiles%q_spline)) deallocate(profiles%q_spline)
        if (allocated(profiles%n_spline)) deallocate(profiles%n_spline)
        if (allocated(profiles%Ti_spline)) deallocate(profiles%Ti_spline)
        if (allocated(profiles%Te_spline)) deallocate(profiles%Te_spline)
        if (allocated(profiles%profile_names)) deallocate(profiles%profile_names)
        
        profiles%n_points = 0
        profiles%n_profiles = 0
        
    end subroutine background_profiles_destroy
    
    subroutine background_profiles_set_from_arrays(profiles, x_grid, q, n, Ti, Te, Vth, Vz, Er, ierr)
        type(background_profiles_t), intent(inout) :: profiles
        real(real64), intent(in) :: x_grid(:), q(:), n(:), Ti(:), Te(:), Vth(:), Vz(:), Er(:)
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Check array sizes
        if (size(x_grid) /= profiles%n_points) then
            ierr = -1
            return
        end if
        
        ! Copy profile data
        profiles%x_grid = x_grid
        profiles%q_profile = q
        profiles%n_profile = n
        profiles%Ti_profile = Ti
        profiles%Te_profile = Te
        profiles%Vth_profile = Vth
        profiles%Vz_profile = Vz
        profiles%Er_profile = Er
        
        ! Calculate magnetic field components (simplified)
        profiles%Bt_profile = 2.5_real64  ! Typical value
        profiles%Bp_profile = profiles%Bt_profile / profiles%q_profile
        profiles%B_profile = sqrt(profiles%Bt_profile**2 + profiles%Bp_profile**2)
        
    end subroutine background_profiles_set_from_arrays
    
    subroutine background_profiles_interpolate(profiles, r, q, n, Ti, Te, Vth, Vz, dPhi0, ierr)
        type(background_profiles_t), intent(in) :: profiles
        real(real64), intent(in) :: r
        real(real64), intent(out) :: q, n, Ti, Te, Vth, Vz, dPhi0
        integer, intent(out) :: ierr
        
        integer :: i
        real(real64) :: t, x1, x2, y1, y2
        
        ierr = 0
        
        ! Boundary cases
        if (r <= profiles%x_grid(1)) then
            q = profiles%q_profile(1)
            n = profiles%n_profile(1)
            Ti = profiles%Ti_profile(1)
            Te = profiles%Te_profile(1)
            Vth = profiles%Vth_profile(1)
            Vz = profiles%Vz_profile(1)
            dPhi0 = -profiles%Er_profile(1)
            return
        end if
        
        if (r >= profiles%x_grid(profiles%n_points)) then
            q = profiles%q_profile(profiles%n_points)
            n = profiles%n_profile(profiles%n_points)
            Ti = profiles%Ti_profile(profiles%n_points)
            Te = profiles%Te_profile(profiles%n_points)
            Vth = profiles%Vth_profile(profiles%n_points)
            Vz = profiles%Vz_profile(profiles%n_points)
            dPhi0 = -profiles%Er_profile(profiles%n_points)
            return
        end if
        
        ! Linear interpolation
        do i = 1, profiles%n_points - 1
            if (r >= profiles%x_grid(i) .and. r <= profiles%x_grid(i+1)) then
                x1 = profiles%x_grid(i)
                x2 = profiles%x_grid(i+1)
                t = (r - x1) / (x2 - x1)
                
                ! Interpolate each profile
                y1 = profiles%q_profile(i)
                y2 = profiles%q_profile(i+1)
                q = y1 + t * (y2 - y1)
                
                y1 = profiles%n_profile(i)
                y2 = profiles%n_profile(i+1)
                n = y1 + t * (y2 - y1)
                
                y1 = profiles%Ti_profile(i)
                y2 = profiles%Ti_profile(i+1)
                Ti = y1 + t * (y2 - y1)
                
                y1 = profiles%Te_profile(i)
                y2 = profiles%Te_profile(i+1)
                Te = y1 + t * (y2 - y1)
                
                y1 = profiles%Vth_profile(i)
                y2 = profiles%Vth_profile(i+1)
                Vth = y1 + t * (y2 - y1)
                
                y1 = profiles%Vz_profile(i)
                y2 = profiles%Vz_profile(i+1)
                Vz = y1 + t * (y2 - y1)
                
                y1 = profiles%Er_profile(i)
                y2 = profiles%Er_profile(i+1)
                dPhi0 = -(y1 + t * (y2 - y1))  ! Er = -dPhi/dr
                
                return
            end if
        end do
        
        ierr = -1  ! Should not reach here
        
    end subroutine background_profiles_interpolate
    
    subroutine background_profiles_transform_to_lab(profiles, r, q, n, Ti, Te, Vth, Vz, dPhi0, ierr)
        type(background_profiles_t), intent(in) :: profiles
        real(real64), intent(in) :: r
        real(real64), intent(inout) :: q, n, Ti, Te, Vth, Vz, dPhi0
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! In lab frame transformation, most quantities remain the same
        ! Only velocities might need Doppler correction
        ! For now, no transformation is applied (placeholder for future)
        
    end subroutine background_profiles_transform_to_lab
    
    !---------------------------------------------------------------------------
    ! Equilibrium calculation procedures
    !---------------------------------------------------------------------------
    subroutine equilibrium_create(equil, R_major, r_minor, B0, ierr)
        type(equilibrium_t), intent(out) :: equil
        real(real64), intent(in) :: R_major, r_minor, B0
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (R_major <= 0.0_real64 .or. r_minor <= 0.0_real64 .or. B0 <= 0.0_real64) then
            ierr = -1
            return
        end if
        
        equil%R_major = R_major
        equil%r_minor = r_minor
        equil%aspect_ratio = R_major / r_minor
        equil%B0 = B0
        equil%Bt0 = B0
        equil%I_plasma = 1.0e6_real64  ! Default value
        equil%pressure_central = 1.0e4_real64  ! Default value
        equil%beta_pol = 0.0_real64
        equil%beta_tor = 0.0_real64
        equil%is_equilibrium_valid = .false.
        equil%force_balance_error = 1.0_real64
        
    end subroutine equilibrium_create
    
    subroutine equilibrium_destroy(equil, ierr)
        type(equilibrium_t), intent(inout) :: equil
        integer, intent(out) :: ierr
        
        ierr = 0
        
        equil%R_major = 0.0_real64
        equil%r_minor = 0.0_real64
        equil%aspect_ratio = 0.0_real64
        equil%B0 = 0.0_real64
        equil%Bt0 = 0.0_real64
        equil%I_plasma = 0.0_real64
        equil%pressure_central = 0.0_real64
        equil%beta_pol = 0.0_real64
        equil%beta_tor = 0.0_real64
        equil%is_equilibrium_valid = .false.
        equil%force_balance_error = 1.0_real64
        
    end subroutine equilibrium_destroy
    
    subroutine equilibrium_calculate(equil, pressure_grad, ierr)
        type(equilibrium_t), intent(inout) :: equil
        real(real64), intent(in) :: pressure_grad
        integer, intent(out) :: ierr
        
        real(real64) :: Bpol
        
        ierr = 0
        
        ! Calculate poloidal field from plasma current
        Bpol = MU_0 * equil%I_plasma / (2.0_real64 * PI * equil%r_minor)
        
        ! Calculate beta values
        equil%beta_pol = 2.0_real64 * MU_0 * pressure_grad / Bpol**2
        equil%beta_tor = 2.0_real64 * MU_0 * pressure_grad / equil%Bt0**2
        
        ! Check equilibrium validity
        equil%is_equilibrium_valid = (equil%beta_pol > 0.0_real64) .and. &
                                    (equil%beta_pol < 1.0_real64) .and. &
                                    (equil%beta_tor > 0.0_real64) .and. &
                                    (equil%beta_tor < 1.0_real64)
        
        ! Calculate force balance error (simplified)
        equil%force_balance_error = abs(pressure_grad - Bpol**2 / (2.0_real64 * MU_0))
        
    end subroutine equilibrium_calculate
    
    subroutine equilibrium_check_force_balance(equil, is_balanced, ierr)
        type(equilibrium_t), intent(in) :: equil
        logical, intent(out) :: is_balanced
        integer, intent(out) :: ierr
        
        ierr = 0
        
        is_balanced = equil%is_equilibrium_valid .and. &
                     (equil%force_balance_error < 1.0e-6_real64)
        
    end subroutine equilibrium_check_force_balance
    
    !---------------------------------------------------------------------------
    ! FLRE zone physics procedures
    !---------------------------------------------------------------------------
    subroutine flre_zone_create(zone, zone_id, r_inner, r_outer, ierr)
        type(flre_zone_t), intent(out) :: zone
        integer, intent(in) :: zone_id
        real(real64), intent(in) :: r_inner, r_outer
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (r_inner >= r_outer) then
            ierr = -1
            return
        end if
        
        zone%zone_id = zone_id
        zone%r_inner = r_inner
        zone%r_outer = r_outer
        zone%omega_ci = 0.0_real64
        zone%omega_ce = 0.0_real64
        zone%rho_i = 0.0_real64
        zone%rho_e = 0.0_real64
        zone%k_perp = 1.0e4_real64  ! Default value
        zone%flre_order_i = 0.0_real64
        zone%flre_order_e = 0.0_real64
        zone%max_flre_order = 5  ! Default maximum order
        zone%epsilon_perp = cmplx(1.0_real64, 0.0_real64, real64)
        zone%epsilon_par = cmplx(1.0_real64, 0.0_real64, real64)
        zone%epsilon_xy = cmplx(0.0_real64, 0.0_real64, real64)
        
    end subroutine flre_zone_create
    
    subroutine flre_zone_destroy(zone, ierr)
        type(flre_zone_t), intent(inout) :: zone
        integer, intent(out) :: ierr
        
        ierr = 0
        
        zone%zone_id = 0
        zone%r_inner = 0.0_real64
        zone%r_outer = 0.0_real64
        zone%omega_ci = 0.0_real64
        zone%omega_ce = 0.0_real64
        zone%rho_i = 0.0_real64
        zone%rho_e = 0.0_real64
        zone%k_perp = 0.0_real64
        zone%flre_order_i = 0.0_real64
        zone%flre_order_e = 0.0_real64
        zone%max_flre_order = 0
        
    end subroutine flre_zone_destroy
    
    subroutine flre_zone_calculate_parameters(zone, Ti, Te, B0, mi, me, ierr)
        type(flre_zone_t), intent(inout) :: zone
        real(real64), intent(in) :: Ti, Te, B0, mi, me
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (B0 <= 0.0_real64 .or. mi <= 0.0_real64 .or. me <= 0.0_real64) then
            ierr = -1
            return
        end if
        
        ! Calculate cyclotron frequencies
        zone%omega_ci = ELEM_CHARGE * B0 / mi
        zone%omega_ce = ELEM_CHARGE * B0 / me
        
        ! Calculate Larmor radii
        zone%rho_i = sqrt(2.0_real64 * Ti / mi) / zone%omega_ci
        zone%rho_e = sqrt(2.0_real64 * Te / me) / zone%omega_ce
        
        ! Calculate FLRE expansion parameters
        zone%flre_order_i = zone%k_perp * zone%rho_i
        zone%flre_order_e = zone%k_perp * zone%rho_e
        
    end subroutine flre_zone_calculate_parameters
    
    subroutine flre_zone_calculate_dielectric(zone, ierr)
        type(flre_zone_t), intent(inout) :: zone
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Simplified dielectric tensor calculation
        zone%epsilon_perp = cmplx(1.0_real64 + zone%flre_order_i, &
                                 0.1_real64 * zone%flre_order_i, real64)
        zone%epsilon_par = cmplx(1.0_real64 - 0.5_real64 * zone%flre_order_e, &
                                0.01_real64, real64)
        zone%epsilon_xy = cmplx(0.0_real64, 0.1_real64 * zone%flre_order_i, real64)
        
    end subroutine flre_zone_calculate_dielectric
    
    !---------------------------------------------------------------------------
    ! Vacuum zone physics procedures
    !---------------------------------------------------------------------------
    subroutine vacuum_zone_create(zone, zone_id, r_inner, r_outer, ierr)
        type(vacuum_zone_t), intent(out) :: zone
        integer, intent(in) :: zone_id
        real(real64), intent(in) :: r_inner, r_outer
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (r_inner >= r_outer) then
            ierr = -1
            return
        end if
        
        zone%zone_id = zone_id
        zone%r_inner = r_inner
        zone%r_outer = r_outer
        zone%c_light = C_LIGHT
        zone%mu_0 = MU_0
        zone%epsilon_0 = EPSILON_0
        zone%wave_impedance = sqrt(MU_0 / EPSILON_0)
        zone%k_vac = cmplx(0.0_real64, 0.0_real64, real64)
        zone%omega = cmplx(0.0_real64, 0.0_real64, real64)
        zone%E_to_B_ratio = C_LIGHT
        
    end subroutine vacuum_zone_create
    
    subroutine vacuum_zone_destroy(zone, ierr)
        type(vacuum_zone_t), intent(inout) :: zone
        integer, intent(out) :: ierr
        
        ierr = 0
        
        zone%zone_id = 0
        zone%r_inner = 0.0_real64
        zone%r_outer = 0.0_real64
        zone%c_light = 0.0_real64
        zone%mu_0 = 0.0_real64
        zone%epsilon_0 = 0.0_real64
        zone%wave_impedance = 0.0_real64
        zone%E_to_B_ratio = 0.0_real64
        
    end subroutine vacuum_zone_destroy
    
    subroutine vacuum_zone_calculate_parameters(zone, omega, ierr)
        type(vacuum_zone_t), intent(inout) :: zone
        complex(real64), intent(in) :: omega
        integer, intent(out) :: ierr
        
        ierr = 0
        
        zone%omega = omega
        zone%k_vac = omega / cmplx(zone%c_light, 0.0_real64, real64)
        
    end subroutine vacuum_zone_calculate_parameters
    
    subroutine vacuum_zone_check_field_relations(zone, Er, Et, Ez, Br, Bt, Bz, is_valid, ierr)
        type(vacuum_zone_t), intent(in) :: zone
        real(real64), intent(in) :: Er, Et, Ez, Br, Bt, Bz
        logical, intent(out) :: is_valid
        integer, intent(out) :: ierr
        
        real(real64) :: tolerance
        
        ierr = 0
        tolerance = 1.0e-10_real64
        
        ! Check E = c*B relationship
        is_valid = (abs(Er - zone%c_light * Br) < tolerance) .and. &
                  (abs(Et - zone%c_light * Bt) < tolerance) .and. &
                  (abs(Ez - zone%c_light * Bz) < tolerance)
        
    end subroutine vacuum_zone_check_field_relations
    
    !---------------------------------------------------------------------------
    ! IMHD zone physics procedures
    !---------------------------------------------------------------------------
    subroutine imhd_zone_create(zone, zone_id, r_inner, r_outer, ierr)
        type(imhd_zone_t), intent(out) :: zone
        integer, intent(in) :: zone_id
        real(real64), intent(in) :: r_inner, r_outer
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (r_inner >= r_outer) then
            ierr = -1
            return
        end if
        
        zone%zone_id = zone_id
        zone%r_inner = r_inner
        zone%r_outer = r_outer
        zone%rho_mass = 1.0e-12_real64  ! Default value
        zone%pressure = 1.0e4_real64     ! Default value
        zone%gamma_adiabatic = 5.0_real64 / 3.0_real64
        zone%v_alfven = 0.0_real64
        zone%c_sound = 0.0_real64
        zone%beta_plasma = 0.0_real64
        zone%vx = 0.0_real64
        zone%vy = 0.0_real64
        zone%vz = 0.0_real64
        zone%divB = 0.0_real64
        zone%curlB_magnitude = 0.0_real64
        
    end subroutine imhd_zone_create
    
    subroutine imhd_zone_destroy(zone, ierr)
        type(imhd_zone_t), intent(inout) :: zone
        integer, intent(out) :: ierr
        
        ierr = 0
        
        zone%zone_id = 0
        zone%r_inner = 0.0_real64
        zone%r_outer = 0.0_real64
        zone%rho_mass = 0.0_real64
        zone%pressure = 0.0_real64
        zone%gamma_adiabatic = 0.0_real64
        zone%v_alfven = 0.0_real64
        zone%c_sound = 0.0_real64
        zone%beta_plasma = 0.0_real64
        
    end subroutine imhd_zone_destroy
    
    subroutine imhd_zone_calculate_speeds(zone, Bx, By, Bz, ierr)
        type(imhd_zone_t), intent(inout) :: zone
        real(real64), intent(in) :: Bx, By, Bz
        integer, intent(out) :: ierr
        
        real(real64) :: B_squared
        
        ierr = 0
        
        B_squared = Bx**2 + By**2 + Bz**2
        
        ! Calculate Alfven speed
        zone%v_alfven = sqrt(B_squared / (MU_0 * zone%rho_mass))
        
        ! Calculate sound speed
        zone%c_sound = sqrt(zone%gamma_adiabatic * zone%pressure / zone%rho_mass)
        
        ! Calculate plasma beta
        zone%beta_plasma = 2.0_real64 * MU_0 * zone%pressure / B_squared
        
    end subroutine imhd_zone_calculate_speeds
    
    subroutine imhd_zone_check_constraints(zone, is_valid, ierr)
        type(imhd_zone_t), intent(inout) :: zone
        logical, intent(out) :: is_valid
        integer, intent(out) :: ierr
        
        ierr = 0
        
        is_valid = (abs(zone%divB) < 1.0e-10_real64) .and. &
                  (zone%curlB_magnitude < 1.0_real64) .and. &
                  (zone%beta_plasma > 0.0_real64) .and. &
                  (zone%beta_plasma < 1.0_real64)
        
    end subroutine imhd_zone_check_constraints
    
    !---------------------------------------------------------------------------
    ! Dispersion relation procedures
    !---------------------------------------------------------------------------
    subroutine dispersion_create(disp, ierr)
        type(dispersion_t), intent(out) :: disp
        integer, intent(out) :: ierr
        
        ierr = 0
        
        disp%omega_pe = 0.0_real64
        disp%omega_pi = 0.0_real64
        disp%omega_ce = 0.0_real64
        disp%omega_ci = 0.0_real64
        disp%omega = cmplx(0.0_real64, 0.0_real64, real64)
        disp%k_par = cmplx(0.0_real64, 0.0_real64, real64)
        disp%k_perp = cmplx(0.0_real64, 0.0_real64, real64)
        disp%epsilon_xx = cmplx(1.0_real64, 0.0_real64, real64)
        disp%epsilon_yy = cmplx(1.0_real64, 0.0_real64, real64)
        disp%epsilon_zz = cmplx(1.0_real64, 0.0_real64, real64)
        disp%epsilon_xy = cmplx(0.0_real64, 0.0_real64, real64)
        disp%n_index_o = cmplx(1.0_real64, 0.0_real64, real64)
        disp%n_index_x = cmplx(1.0_real64, 0.0_real64, real64)
        disp%det_disp = cmplx(0.0_real64, 0.0_real64, real64)
        
    end subroutine dispersion_create
    
    subroutine dispersion_destroy(disp, ierr)
        type(dispersion_t), intent(inout) :: disp
        integer, intent(out) :: ierr
        
        ierr = 0
        
        disp%omega_pe = 0.0_real64
        disp%omega_pi = 0.0_real64
        disp%omega_ce = 0.0_real64
        disp%omega_ci = 0.0_real64
        
    end subroutine dispersion_destroy
    
    subroutine dispersion_calculate_frequencies(disp, n_density, B0, ierr)
        type(dispersion_t), intent(inout) :: disp
        real(real64), intent(in) :: n_density, B0
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Calculate plasma frequencies
        disp%omega_pe = sqrt(n_density * ELEM_CHARGE**2 / (EPSILON_0 * MASS_ELECTRON))
        disp%omega_pi = sqrt(n_density * ELEM_CHARGE**2 / (EPSILON_0 * MASS_PROTON))
        
        ! Calculate cyclotron frequencies
        disp%omega_ce = ELEM_CHARGE * B0 / MASS_ELECTRON
        disp%omega_ci = ELEM_CHARGE * B0 / MASS_PROTON
        
    end subroutine dispersion_calculate_frequencies
    
    subroutine dispersion_calculate_dielectric(disp, omega, ierr)
        type(dispersion_t), intent(inout) :: disp
        complex(real64), intent(in) :: omega
        integer, intent(out) :: ierr
        
        ierr = 0
        
        disp%omega = omega
        
        ! Simplified cold plasma dielectric tensor
        disp%epsilon_xx = cmplx(1.0_real64 - disp%omega_pe**2/real(omega)**2, -0.01_real64, real64)
        disp%epsilon_yy = disp%epsilon_xx
        disp%epsilon_zz = cmplx(1.0_real64 - disp%omega_pe**2/real(omega)**2, -0.005_real64, real64)
        disp%epsilon_xy = cmplx(0.0_real64, disp%omega_pe**2 * disp%omega_ce / real(omega)**3, real64)
        
    end subroutine dispersion_calculate_dielectric
    
    subroutine dispersion_calculate_indices(disp, ierr)
        type(dispersion_t), intent(inout) :: disp
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Calculate refractive indices
        disp%n_index_o = sqrt(disp%epsilon_xx)  ! Ordinary wave
        disp%n_index_x = sqrt(disp%epsilon_zz)  ! Extraordinary wave
        
        ! Calculate dispersion determinant
        disp%det_disp = disp%epsilon_xx - (disp%k_par**2 * C_LIGHT**2 / disp%omega**2)
        
    end subroutine dispersion_calculate_indices
    
    !---------------------------------------------------------------------------
    ! Zone management procedures
    !---------------------------------------------------------------------------
    subroutine zone_manager_create(manager, n_zones, ierr)
        type(zone_manager_t), intent(out) :: manager
        integer, intent(in) :: n_zones
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (n_zones < 1) then
            ierr = -1
            return
        end if
        
        manager%n_zones = n_zones
        manager%zones_coupled = .false.
        
        allocate(manager%zone_types(n_zones), manager%zone_boundaries(n_zones+1), &
                manager%coupling_strength(n_zones-1), stat=ierr)
        
        if (ierr /= 0) then
            ierr = -2
            return
        end if
        
        manager%zone_types = 0
        manager%zone_boundaries = 0.0_real64
        manager%coupling_strength = 0.0_real64
        
    end subroutine zone_manager_create
    
    subroutine zone_manager_destroy(manager, ierr)
        type(zone_manager_t), intent(inout) :: manager
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (allocated(manager%zone_types)) deallocate(manager%zone_types)
        if (allocated(manager%zone_boundaries)) deallocate(manager%zone_boundaries)
        if (allocated(manager%coupling_strength)) deallocate(manager%coupling_strength)
        if (allocated(manager%flre_zones)) deallocate(manager%flre_zones)
        if (allocated(manager%vacuum_zones)) deallocate(manager%vacuum_zones)
        if (allocated(manager%imhd_zones)) deallocate(manager%imhd_zones)
        
        manager%n_zones = 0
        manager%zones_coupled = .false.
        
    end subroutine zone_manager_destroy
    
    subroutine zone_manager_setup_zones(manager, zone_types, zone_boundaries, ierr)
        type(zone_manager_t), intent(inout) :: manager
        integer, intent(in) :: zone_types(:)
        real(real64), intent(in) :: zone_boundaries(:)
        integer, intent(out) :: ierr
        
        integer :: i, n_flre, n_vacuum, n_imhd
        
        ierr = 0
        
        if (size(zone_types) /= manager%n_zones) then
            ierr = -1
            return
        end if
        
        if (size(zone_boundaries) /= manager%n_zones + 1) then
            ierr = -2
            return
        end if
        
        manager%zone_types = zone_types
        manager%zone_boundaries = zone_boundaries
        
        ! Count zone types
        n_flre = count(zone_types == 1)
        n_vacuum = count(zone_types == 2)
        n_imhd = count(zone_types == 3)
        
        ! Allocate zone storage
        if (n_flre > 0) allocate(manager%flre_zones(n_flre))
        if (n_vacuum > 0) allocate(manager%vacuum_zones(n_vacuum))
        if (n_imhd > 0) allocate(manager%imhd_zones(n_imhd))
        
        ! Initialize zones
        n_flre = 0
        n_vacuum = 0
        n_imhd = 0
        
        do i = 1, manager%n_zones
            select case (zone_types(i))
            case (1)  ! FLRE
                n_flre = n_flre + 1
                call flre_zone_create(manager%flre_zones(n_flre), i, &
                                     zone_boundaries(i), zone_boundaries(i+1), ierr)
            case (2)  ! Vacuum
                n_vacuum = n_vacuum + 1
                call vacuum_zone_create(manager%vacuum_zones(n_vacuum), i, &
                                       zone_boundaries(i), zone_boundaries(i+1), ierr)
            case (3)  ! IMHD
                n_imhd = n_imhd + 1
                call imhd_zone_create(manager%imhd_zones(n_imhd), i, &
                                     zone_boundaries(i), zone_boundaries(i+1), ierr)
            end select
        end do
        
    end subroutine zone_manager_setup_zones
    
    subroutine zone_manager_identify_zone(manager, r, zone_id, zone_type, ierr)
        type(zone_manager_t), intent(in) :: manager
        real(real64), intent(in) :: r
        integer, intent(out) :: zone_id, zone_type
        integer, intent(out) :: ierr
        
        integer :: i
        
        ierr = 0
        zone_id = 0
        zone_type = 0
        
        do i = 1, manager%n_zones
            if (r >= manager%zone_boundaries(i) .and. r < manager%zone_boundaries(i+1)) then
                zone_id = i
                zone_type = manager%zone_types(i)
                return
            end if
        end do
        
        ierr = -1  ! Point not in any zone
        
    end subroutine zone_manager_identify_zone
    
    subroutine zone_manager_check_coupling(manager, ierr)
        type(zone_manager_t), intent(inout) :: manager
        integer, intent(out) :: ierr
        
        integer :: i
        
        ierr = 0
        
        ! Calculate coupling strength between adjacent zones
        do i = 1, manager%n_zones - 1
            ! Simplified coupling calculation
            manager%coupling_strength(i) = 0.8_real64  ! Default good coupling
        end do
        
        manager%zones_coupled = all(manager%coupling_strength > 0.5_real64)
        
    end subroutine zone_manager_check_coupling
    
    !---------------------------------------------------------------------------
    ! Lab frame transformation procedures
    !---------------------------------------------------------------------------
    subroutine lab_frame_transform_create(transform, omega_lab, v_rot, R_major, ierr)
        type(lab_frame_transform_t), intent(out) :: transform
        real(real64), intent(in) :: omega_lab, v_rot, R_major
        integer, intent(out) :: ierr
        
        ierr = 0
        
        transform%omega_lab = omega_lab
        transform%v_rot = v_rot
        transform%R_major = R_major
        transform%doppler_factor = 1.0_real64 + v_rot / C_LIGHT
        transform%omega_plasma = omega_lab / transform%doppler_factor
        transform%gamma_relativistic = 1.0_real64  ! Non-relativistic
        
        ! Initialize transformation matrix (identity)
        transform%transformation_matrix = 0.0_real64
        transform%transformation_matrix(1,1) = 1.0_real64
        transform%transformation_matrix(2,2) = 1.0_real64
        transform%transformation_matrix(3,3) = 1.0_real64
        
        transform%E_lab = 0.0_real64
        transform%E_plasma = 0.0_real64
        transform%B_lab = 0.0_real64
        transform%B_plasma = 0.0_real64
        
    end subroutine lab_frame_transform_create
    
    subroutine lab_frame_transform_destroy(transform, ierr)
        type(lab_frame_transform_t), intent(inout) :: transform
        integer, intent(out) :: ierr
        
        ierr = 0
        
        transform%omega_lab = 0.0_real64
        transform%omega_plasma = 0.0_real64
        transform%v_rot = 0.0_real64
        transform%R_major = 0.0_real64
        transform%doppler_factor = 1.0_real64
        
    end subroutine lab_frame_transform_destroy
    
    subroutine lab_frame_transform_calculate(transform, theta, ierr)
        type(lab_frame_transform_t), intent(inout) :: transform
        real(real64), intent(in) :: theta
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Set up rotation transformation matrix
        transform%transformation_matrix(1,1) = cos(theta)
        transform%transformation_matrix(1,2) = -sin(theta)
        transform%transformation_matrix(2,1) = sin(theta)
        transform%transformation_matrix(2,2) = cos(theta)
        transform%transformation_matrix(3,3) = 1.0_real64
        
    end subroutine lab_frame_transform_calculate
    
    subroutine lab_frame_transform_fields(transform, E_in, B_in, E_out, B_out, to_lab, ierr)
        type(lab_frame_transform_t), intent(inout) :: transform
        real(real64), intent(in) :: E_in(3), B_in(3)
        real(real64), intent(out) :: E_out(3), B_out(3)
        logical, intent(in) :: to_lab
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (to_lab) then
            ! Transform from plasma to lab frame
            E_out = matmul(transpose(transform%transformation_matrix), E_in)
            B_out = matmul(transpose(transform%transformation_matrix), B_in)
        else
            ! Transform from lab to plasma frame
            E_out = matmul(transform%transformation_matrix, E_in)
            B_out = matmul(transform%transformation_matrix, B_in)
        end if
        
    end subroutine lab_frame_transform_fields
    
    !---------------------------------------------------------------------------
    ! F0 distribution moment procedures
    !---------------------------------------------------------------------------
    subroutine f0_moments_create(moments, n_density, Ti_eV, Te_eV, ierr)
        type(f0_moments_t), intent(out) :: moments
        real(real64), intent(in) :: n_density, Ti_eV, Te_eV
        integer, intent(out) :: ierr
        
        ierr = 0
        
        moments%n_density = n_density
        moments%Ti = Ti_eV * EV_TO_JOULE
        moments%Te = Te_eV * EV_TO_JOULE
        moments%vth_i = sqrt(2.0_real64 * moments%Ti / MASS_PROTON)
        moments%vth_e = sqrt(2.0_real64 * moments%Te / MASS_ELECTRON)
        moments%f0_moment_0 = n_density
        moments%f0_moment_1 = 0.0_real64
        moments%f0_moment_2 = 0.0_real64
        moments%pressure_i = 0.0_real64
        moments%pressure_e = 0.0_real64
        moments%energy_density = 0.0_real64
        moments%beta_thermal = 0.0_real64
        
    end subroutine f0_moments_create
    
    subroutine f0_moments_destroy(moments, ierr)
        type(f0_moments_t), intent(inout) :: moments
        integer, intent(out) :: ierr
        
        ierr = 0
        
        moments%n_density = 0.0_real64
        moments%Ti = 0.0_real64
        moments%Te = 0.0_real64
        moments%vth_i = 0.0_real64
        moments%vth_e = 0.0_real64
        
    end subroutine f0_moments_destroy
    
    subroutine f0_moments_calculate(moments, ierr)
        type(f0_moments_t), intent(inout) :: moments
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Calculate distribution moments
        moments%f0_moment_0 = moments%n_density
        moments%f0_moment_1 = 0.0_real64  ! Zero mean velocity
        moments%f0_moment_2 = moments%n_density * moments%vth_i**2
        
    end subroutine f0_moments_calculate
    
    subroutine f0_moments_calculate_pressure(moments, B_field, ierr)
        type(f0_moments_t), intent(inout) :: moments
        real(real64), intent(in) :: B_field
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Calculate pressures
        moments%pressure_i = moments%n_density * moments%Ti
        moments%pressure_e = moments%n_density * moments%Te
        moments%energy_density = 1.5_real64 * (moments%pressure_i + moments%pressure_e)
        
        ! Calculate thermal beta
        moments%beta_thermal = (moments%pressure_i + moments%pressure_e) / &
                              (B_field**2 / (2.0_real64 * MU_0))
        
    end subroutine f0_moments_calculate_pressure
    
    !---------------------------------------------------------------------------
    ! Collision frequency procedures
    !---------------------------------------------------------------------------
    subroutine collision_freq_create(coll, n_density, Ti_eV, Te_eV, Z_eff, ierr)
        type(collision_freq_t), intent(out) :: coll
        real(real64), intent(in) :: n_density, Ti_eV, Te_eV, Z_eff
        integer, intent(out) :: ierr
        
        ierr = 0
        
        coll%n_density = n_density
        coll%Ti = Ti_eV
        coll%Te = Te_eV
        coll%Z_eff = Z_eff
        coll%ln_lambda = 15.0_real64  ! Default Coulomb logarithm
        coll%nu_ii = 0.0_real64
        coll%nu_ee = 0.0_real64
        coll%nu_ei = 0.0_real64
        coll%nu_ie = 0.0_real64
        coll%tau_ii = 0.0_real64
        coll%tau_ee = 0.0_real64
        coll%tau_ei = 0.0_real64
        coll%slowing_down_time = 0.0_real64
        
    end subroutine collision_freq_create
    
    subroutine collision_freq_destroy(coll, ierr)
        type(collision_freq_t), intent(inout) :: coll
        integer, intent(out) :: ierr
        
        ierr = 0
        
        coll%n_density = 0.0_real64
        coll%Ti = 0.0_real64
        coll%Te = 0.0_real64
        coll%Z_eff = 0.0_real64
        coll%ln_lambda = 0.0_real64
        
    end subroutine collision_freq_destroy
    
    subroutine collision_freq_calculate(coll, ierr)
        type(collision_freq_t), intent(inout) :: coll
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Calculate collision frequencies
        coll%nu_ii = 4.8e-14_real64 * coll%Z_eff**4 * coll%n_density * &
                    coll%ln_lambda / coll%Ti**1.5_real64
        
        coll%nu_ee = 2.9e-12_real64 * coll%n_density * coll%ln_lambda / &
                    coll%Te**1.5_real64
        
        coll%nu_ei = 2.9e-12_real64 * coll%Z_eff * coll%n_density * &
                    coll%ln_lambda / coll%Te**1.5_real64
        
        coll%nu_ie = coll%nu_ei * MASS_ELECTRON / MASS_PROTON
        
    end subroutine collision_freq_calculate
    
    subroutine collision_freq_calculate_times(coll, ierr)
        type(collision_freq_t), intent(inout) :: coll
        integer, intent(out) :: ierr
        
        ierr = 0
        
        if (coll%nu_ii > 0.0_real64) coll%tau_ii = 1.0_real64 / coll%nu_ii
        if (coll%nu_ee > 0.0_real64) coll%tau_ee = 1.0_real64 / coll%nu_ee
        if (coll%nu_ei > 0.0_real64) coll%tau_ei = 1.0_real64 / coll%nu_ei
        
        coll%slowing_down_time = coll%tau_ei * (MASS_PROTON / MASS_ELECTRON)
        
    end subroutine collision_freq_calculate_times
    
    !---------------------------------------------------------------------------
    ! Magnetic field calculation procedures
    !---------------------------------------------------------------------------
    subroutine magnetic_field_create(mag_field, R_major, r_minor, I_plasma, Bt, ierr)
        type(magnetic_field_t), intent(out) :: mag_field
        real(real64), intent(in) :: R_major, r_minor, I_plasma, Bt
        integer, intent(out) :: ierr
        
        ierr = 0
        
        mag_field%R_major = R_major
        mag_field%r_minor = r_minor
        mag_field%I_plasma = I_plasma
        mag_field%mu_0 = MU_0
        mag_field%Bt = Bt
        mag_field%Bp = 0.0_real64
        mag_field%B_total = 0.0_real64
        mag_field%q_safety = 2.0_real64  ! Default
        mag_field%psi_flux = 0.0_real64
        mag_field%grad_B = 0.0_real64
        mag_field%curvature_drift = 0.0_real64
        mag_field%magnetic_shear = 0.1_real64  ! Default
        mag_field%B_field = 0.0_real64
        mag_field%grad_B_vec = 0.0_real64
        
    end subroutine magnetic_field_create
    
    subroutine magnetic_field_destroy(mag_field, ierr)
        type(magnetic_field_t), intent(inout) :: mag_field
        integer, intent(out) :: ierr
        
        ierr = 0
        
        mag_field%R_major = 0.0_real64
        mag_field%r_minor = 0.0_real64
        mag_field%I_plasma = 0.0_real64
        mag_field%Bt = 0.0_real64
        mag_field%Bp = 0.0_real64
        mag_field%B_total = 0.0_real64
        mag_field%q_safety = 0.0_real64
        
    end subroutine magnetic_field_destroy
    
    subroutine magnetic_field_calculate(mag_field, ierr)
        type(magnetic_field_t), intent(inout) :: mag_field
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Calculate poloidal field
        mag_field%Bp = mag_field%Bt / mag_field%q_safety
        
        ! Calculate total field
        mag_field%B_total = sqrt(mag_field%Bp**2 + mag_field%Bt**2)
        
        ! Calculate flux
        mag_field%psi_flux = mag_field%Bp * mag_field%R_major * mag_field%r_minor
        
        ! Set field components
        mag_field%B_field = [0.1_real64, mag_field%Bp, mag_field%Bt]
        
    end subroutine magnetic_field_calculate
    
    subroutine magnetic_field_calculate_gradients(mag_field, ierr)
        type(magnetic_field_t), intent(inout) :: mag_field
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Calculate gradient (simplified 1/R dependence)
        mag_field%grad_B = mag_field%Bt / mag_field%R_major
        
        ! Calculate curvature drift
        mag_field%curvature_drift = mag_field%grad_B / mag_field%B_total
        
        ! Gradient vector
        mag_field%grad_B_vec = [mag_field%grad_B, 0.0_real64, 0.0_real64]
        
    end subroutine magnetic_field_calculate_gradients
    
    !---------------------------------------------------------------------------
    ! Physics utility procedures
    !---------------------------------------------------------------------------
    subroutine physics_spline_interpolate(x_data, y_data, x_test, y_interp, ierr)
        real(real64), intent(in) :: x_data(:), y_data(:), x_test(:)
        real(real64), intent(out) :: y_interp(:)
        integer, intent(out) :: ierr
        
        integer :: i, j, n_data, n_test
        real(real64) :: t
        
        ierr = 0
        n_data = size(x_data)
        n_test = size(x_test)
        
        ! Simple linear interpolation
        do i = 1, n_test
            if (x_test(i) <= x_data(1)) then
                y_interp(i) = y_data(1)
            else if (x_test(i) >= x_data(n_data)) then
                y_interp(i) = y_data(n_data)
            else
                do j = 1, n_data - 1
                    if (x_test(i) >= x_data(j) .and. x_test(i) <= x_data(j+1)) then
                        t = (x_test(i) - x_data(j)) / (x_data(j+1) - x_data(j))
                        y_interp(i) = y_data(j) + t * (y_data(j+1) - y_data(j))
                        exit
                    end if
                end do
            end if
        end do
        
    end subroutine physics_spline_interpolate
    
    subroutine physics_check_zone_continuity(E_left, E_right, B_left, B_right, &
                                           current_jump, is_continuous, ierr)
        real(real64), intent(in) :: E_left(3), E_right(3), B_left(3), B_right(3)
        real(real64), intent(in) :: current_jump(3)
        logical, intent(out) :: is_continuous
        integer, intent(out) :: ierr
        
        real(real64) :: tolerance
        
        ierr = 0
        tolerance = 1.0e-10_real64
        
        ! Check tangential field continuity
        is_continuous = (abs(E_left(1) - E_right(1)) < tolerance) .and. &
                       (abs(E_left(2) - E_right(2)) < tolerance) .and. &
                       (abs(B_left(2) - B_right(2)) < tolerance) .and. &
                       (abs(B_left(3) - B_right(3)) < tolerance)
        
    end subroutine physics_check_zone_continuity
    
    subroutine physics_validate_parameters(Ti, Te, n_density, B_field, q_safety, &
                                         beta_pressure, omega_freq, is_valid, ierr)
        real(real64), intent(in) :: Ti, Te, n_density, B_field, q_safety
        real(real64), intent(in) :: beta_pressure, omega_freq
        logical, intent(out) :: is_valid
        integer, intent(out) :: ierr
        
        logical :: is_valid_temperature, is_valid_density, is_valid_field
        logical :: is_valid_safety, is_valid_beta, is_valid_frequency
        
        ierr = 0
        
        ! Validate temperature ranges
        is_valid_temperature = (Ti > 0.0_real64) .and. (Ti < 1.0e5_real64) .and. &
                              (Te > 0.0_real64) .and. (Te < 1.0e5_real64) .and. &
                              (Te >= Ti)
        
        ! Validate density ranges
        is_valid_density = (n_density > 1.0e17_real64) .and. (n_density < 1.0e22_real64)
        
        ! Validate magnetic field ranges
        is_valid_field = (B_field > 0.1_real64) .and. (B_field < 20.0_real64)
        
        ! Validate safety factor ranges
        is_valid_safety = (q_safety > 1.0_real64) .and. (q_safety < 10.0_real64)
        
        ! Validate beta parameter
        is_valid_beta = (beta_pressure > 0.0_real64) .and. (beta_pressure < 1.0_real64)
        
        ! Validate frequency ranges
        is_valid_frequency = (omega_freq > 1.0e6_real64) .and. (omega_freq < 1.0e12_real64)
        
        is_valid = is_valid_temperature .and. is_valid_density .and. is_valid_field .and. &
                  is_valid_safety .and. is_valid_beta .and. is_valid_frequency
        
    end subroutine physics_validate_parameters

end module kilca_physics_m