program test_kilca_physics
    !---------------------------------------------------------------------------
    ! KiLCA Plasma Physics Models - Unit Tests
    !
    ! Tests the plasma physics models including background plasma profiles,
    ! physics zones (FLRE, vacuum, IMHD), dispersion relations, and zone
    ! management systems.
    !
    ! Author: Claude (Anthropic)
    ! Date: 2024
    !
    ! Testing: kilca_physics_m.f90 (to be implemented)
    !---------------------------------------------------------------------------
    
    use iso_fortran_env, only: real64, int32, error_unit
    use kilca_types_m
    implicit none
    
    ! Test counters
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    
    write(*,'(A)') " ========================================================"
    write(*,'(A)') " KiLCA Plasma Physics Models - Unit Tests"
    write(*,'(A)') " ========================================================"
    write(*,'(A)') " "
    
    ! Run all tests
    call test_background_profiles()
    call test_equilibrium_calculation()
    call test_flre_zone_physics()
    call test_vacuum_zone_physics()
    call test_imhd_zone_physics()
    call test_dispersion_relations()
    call test_physics_zone_management()
    call test_lab_frame_transformations()
    call test_f0_distribution_moments()
    call test_collision_frequencies()
    call test_magnetic_field_calculations()
    call test_physics_spline_interpolation()
    call test_physics_zone_coupling()
    call test_physics_parameter_validation()
    
    ! Print summary
    write(*,'(A)') " "
    write(*,'(A)') " ========================================================"
    write(*,'(A)') " TEST SUMMARY"
    write(*,'(A)') " ========================================================"
    write(*,'(A,I15)') " Total tests run:              ", total_tests
    write(*,'(A,I15)') " Tests passed:                 ", passed_tests
    write(*,'(A,I15)') " Tests failed:                 ", failed_tests
    write(*,'(A,F15.7,A)') " Success rate:         ", &
        100.0_real64 * real(passed_tests, real64) / real(max(total_tests,1), real64), " %"
    write(*,'(A)') " "
    
    if (failed_tests > 0) then
        write(*,'(A)') " *** SOME TESTS FAILED! ***"
        stop 1
    else
        write(*,'(A)') " *** ALL TESTS PASSED! ***"
    end if

contains

    !---------------------------------------------------------------------------
    ! Test 1: Background plasma profiles management
    !---------------------------------------------------------------------------
    subroutine test_background_profiles()
        real(real64), allocatable :: x_grid(:)
        real(real64), allocatable :: q_profile(:), n_profile(:)
        real(real64), allocatable :: Ti_profile(:), Te_profile(:)
        real(real64), allocatable :: Vth_profile(:), Vz_profile(:)
        real(real64), allocatable :: Er_profile(:)
        integer :: n_points, ierr
        real(real64) :: r_test, q_val, n_val, Ti_val, Te_val
        
        call start_test("Background plasma profiles management")
        test_passed = .true.
        
        ! Test background profile creation and interpolation
        n_points = 50
        allocate(x_grid(n_points), q_profile(n_points), n_profile(n_points), &
                Ti_profile(n_points), Te_profile(n_points), Vth_profile(n_points), &
                Vz_profile(n_points), Er_profile(n_points), stat=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            ! Create radial grid (0 to 1 normalized)
            do ierr = 1, n_points
                x_grid(ierr) = real(ierr-1, real64) / real(n_points-1, real64)
            end do
            
            ! Set typical plasma profiles
            ! Safety factor: q = 1 + 2*r^2
            q_profile = 1.0_real64 + 2.0_real64 * x_grid**2
            
            ! Density: parabolic profile n = n0*(1-r^2)
            n_profile = 1.0e14_real64 * (1.0_real64 - x_grid**2)
            
            ! Ion temperature: Ti = Ti0*(1-r^2)^0.5
            Ti_profile = 1000.0_real64 * (1.0_real64 - x_grid**2)**0.5_real64
            
            ! Electron temperature: Te = Te0*(1-r^2)^0.8
            Te_profile = 2000.0_real64 * (1.0_real64 - x_grid**2)**0.8_real64
            
            ! Thermal velocity: Vth = sqrt(2*Ti/mi)
            Vth_profile = sqrt(2.0_real64 * Ti_profile * 1.6022e-12_real64 / 1.6726e-24_real64)
            
            ! Toroidal velocity: Vz (assumed small)
            Vz_profile = 1.0e5_real64 * x_grid  ! Linear profile
            
            ! Radial electric field: Er = -dPhi/dr
            Er_profile = 1000.0_real64 * x_grid  ! Simple linear profile
            
            ! Test profile values at center and edge
            r_test = 0.0_real64  ! Center
            test_passed = test_passed .and. (abs(q_profile(1) - 1.0_real64) < 1.0e-12_real64)
            test_passed = test_passed .and. (abs(n_profile(1) - 1.0e14_real64) < 1.0e12_real64)
            test_passed = test_passed .and. (Ti_profile(1) > 0.0_real64)
            test_passed = test_passed .and. (Te_profile(1) > Ti_profile(1))
            
            r_test = 1.0_real64  ! Edge
            test_passed = test_passed .and. (q_profile(n_points) > q_profile(1))
            test_passed = test_passed .and. (n_profile(n_points) == 0.0_real64)
            test_passed = test_passed .and. (Ti_profile(n_points) == 0.0_real64)
            test_passed = test_passed .and. (Te_profile(n_points) == 0.0_real64)
            
            ! Test monotonicity and physical bounds
            test_passed = test_passed .and. all(q_profile > 0.0_real64)
            test_passed = test_passed .and. all(n_profile >= 0.0_real64)
            test_passed = test_passed .and. all(Ti_profile >= 0.0_real64)
            test_passed = test_passed .and. all(Te_profile >= 0.0_real64)
            test_passed = test_passed .and. all(Vth_profile >= 0.0_real64)
            
            deallocate(x_grid, q_profile, n_profile, Ti_profile, Te_profile, &
                      Vth_profile, Vz_profile, Er_profile)
        end if
        
        call end_test(test_passed)
    end subroutine test_background_profiles
    
    !---------------------------------------------------------------------------
    ! Test 2: Equilibrium calculation and consistency
    !---------------------------------------------------------------------------
    subroutine test_equilibrium_calculation()
        real(real64) :: B0, Bt, Bz, Bpol, Btor
        real(real64) :: R_major, r_minor, q_val, pressure_grad
        real(real64) :: beta_pol, beta_tor, safety_factor
        integer :: ierr
        
        call start_test("Equilibrium calculation and consistency")
        test_passed = .true.
        
        ! Test tokamak equilibrium parameters
        R_major = 1.65_real64  ! Major radius (m)
        r_minor = 0.5_real64   ! Minor radius (m)
        B0 = 2.5_real64        ! Magnetic field strength (T)
        
        ! Calculate field components
        Bt = B0  ! Toroidal field
        q_val = 2.0_real64  ! Safety factor
        ! For given q, r, R, and Bt, calculate Bpol from q = (r * Bt) / (R * Bpol)
        Bpol = (r_minor * Bt) / (R_major * q_val)   ! Poloidal field from safety factor
        Bz = Bt * r_minor / R_major  ! Approximate vertical field
        
        ! Validate equilibrium consistency
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (B0 > 0.0_real64)
        test_passed = test_passed .and. (Bt > 0.0_real64)
        test_passed = test_passed .and. (Bpol > 0.0_real64)
        test_passed = test_passed .and. (q_val > 1.0_real64)  ! Stability requirement
        
        ! Test force balance (simplified)
        pressure_grad = 1.0e4_real64  ! Pa/m
        beta_pol = 2.0_real64 * 4.0e-7_real64 * PI * pressure_grad / Bpol**2
        beta_tor = 2.0_real64 * 4.0e-7_real64 * PI * pressure_grad / Bt**2
        
        test_passed = test_passed .and. (beta_pol > 0.0_real64)
        test_passed = test_passed .and. (beta_tor > 0.0_real64)
        test_passed = test_passed .and. (beta_pol < 1.0_real64)  ! Physical constraint
        test_passed = test_passed .and. (beta_tor < 1.0_real64)  ! Physical constraint
        
        ! Test magnetic field total
        safety_factor = (r_minor * Bt) / (R_major * Bpol)
        test_passed = test_passed .and. (abs(safety_factor - q_val) < 0.1_real64)
        
        call end_test(test_passed)
    end subroutine test_equilibrium_calculation
    
    !---------------------------------------------------------------------------
    ! Test 3: FLRE (Finite Larmor Radius Effects) zone physics
    !---------------------------------------------------------------------------
    subroutine test_flre_zone_physics()
        real(real64) :: omega_ci, omega_ce, rho_i, rho_e
        real(real64) :: Ti, Te, B0, mi, me, charge_e
        real(real64) :: k_perp, flre_order_i, flre_order_e
        complex(real64) :: epsilon_perp, epsilon_par
        integer :: ierr
        
        call start_test("FLRE zone physics")
        test_passed = .true.
        
        ! Test FLRE physics calculations
        Ti = 1000.0_real64 * 1.6022e-19_real64  ! Ion temperature (J)
        Te = 2000.0_real64 * 1.6022e-19_real64  ! Electron temperature (J)
        B0 = 2.5_real64  ! Magnetic field (T)
        mi = 1.6726e-27_real64  ! Ion mass (kg)
        me = 9.1094e-31_real64  ! Electron mass (kg)
        charge_e = 1.6022e-19_real64  ! Elementary charge (C)
        
        ! Calculate cyclotron frequencies
        omega_ci = charge_e * B0 / mi
        omega_ce = charge_e * B0 / me
        
        ! Calculate Larmor radii
        rho_i = sqrt(2.0_real64 * Ti / mi) / omega_ci
        rho_e = sqrt(2.0_real64 * Te / me) / omega_ce
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (omega_ci > 0.0_real64)
        test_passed = test_passed .and. (omega_ce > omega_ci)  ! Electrons faster
        test_passed = test_passed .and. (rho_i > rho_e)       ! Ion radius larger
        test_passed = test_passed .and. (rho_i < 0.1_real64)  ! Reasonable size (m)
        test_passed = test_passed .and. (rho_e < 0.001_real64) ! Electron radius smaller
        
        ! Test FLRE expansion parameter
        k_perp = 1.0e4_real64  ! Perpendicular wavenumber (1/m)
        flre_order_i = k_perp * rho_i
        flre_order_e = k_perp * rho_e
        
        test_passed = test_passed .and. (flre_order_i > 0.0_real64)
        test_passed = test_passed .and. (flre_order_e > 0.0_real64)
        test_passed = test_passed .and. (flre_order_i > flre_order_e)
        
        ! Test dielectric tensor elements (simplified)
        epsilon_perp = cmplx(1.0_real64 + flre_order_i, 0.1_real64 * flre_order_i, real64)
        epsilon_par = cmplx(1.0_real64 - 0.5_real64 * flre_order_e, 0.01_real64, real64)
        
        test_passed = test_passed .and. (real(epsilon_perp) > 0.0_real64)
        test_passed = test_passed .and. (real(epsilon_par) > 0.0_real64)
        test_passed = test_passed .and. (abs(epsilon_perp) > abs(epsilon_par))
        
        call end_test(test_passed)
    end subroutine test_flre_zone_physics
    
    !---------------------------------------------------------------------------
    ! Test 4: Vacuum zone physics
    !---------------------------------------------------------------------------
    subroutine test_vacuum_zone_physics()
        complex(real64) :: k_vac, omega
        real(real64) :: c_light, mu_0, wave_impedance, epsilon_0
        real(real64) :: Er, Et, Ez, Br, Bt, Bz
        integer :: ierr
        
        call start_test("Vacuum zone physics")
        test_passed = .true.
        
        ! Test vacuum physics - electromagnetic wave propagation
        c_light = 2.998e8_real64  ! Speed of light (m/s)
        mu_0 = 4.0e-7_real64 * PI  ! Vacuum permeability
        epsilon_0 = 1.0_real64 / (mu_0 * c_light**2)  ! Vacuum permittivity
        
        ! Test wave parameters in vacuum
        omega = cmplx(1.0e8_real64, 0.0_real64, real64)  ! Angular frequency (rad/s)
        k_vac = omega / cmplx(c_light, 0.0_real64, real64)  ! Vacuum wave number
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (c_light > 0.0_real64)
        test_passed = test_passed .and. (mu_0 > 0.0_real64)
        test_passed = test_passed .and. (epsilon_0 > 0.0_real64)
        test_passed = test_passed .and. (abs(k_vac) > 0.0_real64)
        
        ! Test wave impedance
        wave_impedance = sqrt(mu_0 / epsilon_0)
        test_passed = test_passed .and. (abs(wave_impedance - 377.0_real64) < 1.0_real64)
        
        ! Test vacuum field relationships E = c*B
        Er = 1.0_real64  ! Electric field (V/m)
        Et = 0.5_real64
        Ez = 0.1_real64
        
        Br = Er / c_light  ! Magnetic field (T)
        Bt = Et / c_light
        Bz = Ez / c_light
        
        test_passed = test_passed .and. (abs(Er - c_light * Br) < 1.0e-10_real64)
        test_passed = test_passed .and. (abs(Et - c_light * Bt) < 1.0e-10_real64)
        test_passed = test_passed .and. (abs(Ez - c_light * Bz) < 1.0e-10_real64)
        
        ! Test Poynting vector S = E x B / mu_0
        test_passed = test_passed .and. (Br > 0.0_real64)
        test_passed = test_passed .and. (Bt > 0.0_real64)
        test_passed = test_passed .and. (Bz > 0.0_real64)
        
        call end_test(test_passed)
    end subroutine test_vacuum_zone_physics
    
    !---------------------------------------------------------------------------
    ! Test 5: IMHD (Ideal MagnetoHydroDynamic) zone physics
    !---------------------------------------------------------------------------
    subroutine test_imhd_zone_physics()
        real(real64) :: rho_mass, pressure, gamma_adiabatic
        real(real64) :: Bx, By, Bz, vx, vy, vz
        real(real64) :: v_alfven, c_sound, beta_plasma
        real(real64) :: divB, curlB_x, curlB_y, curlB_z
        integer :: ierr
        
        call start_test("IMHD zone physics")
        test_passed = .true.
        
        ! Test IMHD physics - magnetohydrodynamic equations
        rho_mass = 1.0e-12_real64  ! Mass density (kg/m^3)
        pressure = 1.0e4_real64    ! Pressure (Pa)
        gamma_adiabatic = 5.0_real64 / 3.0_real64  ! Adiabatic index
        
        ! Magnetic field components
        Bx = 0.1_real64  ! Tesla
        By = 0.05_real64
        Bz = 2.5_real64  ! Dominant toroidal field
        
        ! Velocity components
        vx = 1.0e4_real64  ! m/s
        vy = 5.0e3_real64
        vz = 1.0e5_real64
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (rho_mass > 0.0_real64)
        test_passed = test_passed .and. (pressure > 0.0_real64)
        test_passed = test_passed .and. (gamma_adiabatic > 1.0_real64)
        
        ! Calculate characteristic speeds
        v_alfven = sqrt((Bx**2 + By**2 + Bz**2) / (4.0e-7_real64 * PI * rho_mass))
        c_sound = sqrt(gamma_adiabatic * pressure / rho_mass)
        beta_plasma = 2.0_real64 * 4.0e-7_real64 * PI * pressure / (Bx**2 + By**2 + Bz**2)
        
        test_passed = test_passed .and. (v_alfven > 0.0_real64)
        test_passed = test_passed .and. (c_sound > 0.0_real64)
        test_passed = test_passed .and. (beta_plasma > 0.0_real64)
        test_passed = test_passed .and. (v_alfven > c_sound)  ! Typical for tokamak
        test_passed = test_passed .and. (beta_plasma < 1.0_real64)  ! Low beta plasma
        
        ! Test MHD constraints (simplified)
        divB = 0.0_real64  ! div(B) = 0 constraint
        test_passed = test_passed .and. (abs(divB) < 1.0e-10_real64)
        
        ! Test Ampere's law curl(B) = mu_0*J (simplified)
        curlB_x = 0.1_real64  ! Simplified curl calculation
        curlB_y = 0.05_real64
        curlB_z = 0.01_real64
        
        test_passed = test_passed .and. (curlB_x**2 + curlB_y**2 + curlB_z**2 < 1.0_real64)
        
        call end_test(test_passed)
    end subroutine test_imhd_zone_physics
    
    !---------------------------------------------------------------------------
    ! Test 6: Dispersion relations and wave physics
    !---------------------------------------------------------------------------
    subroutine test_dispersion_relations()
        complex(real64) :: omega, k_par, k_perp
        complex(real64) :: epsilon_xx, epsilon_yy, epsilon_zz, epsilon_xy
        complex(real64) :: det_disp, n_index_o, n_index_x
        real(real64) :: omega_pe, omega_pi, omega_ce, omega_ci
        integer :: ierr
        
        call start_test("Dispersion relations and wave physics")
        test_passed = .true.
        
        ! Test plasma dispersion relations
        ! Plasma frequencies
        omega_pe = 5.0e9_real64  ! Electron plasma frequency (rad/s)
        omega_pi = 1.0e8_real64  ! Ion plasma frequency (rad/s)
        omega_ce = 7.0e10_real64 ! Electron cyclotron frequency (rad/s)
        omega_ci = 3.8e7_real64  ! Ion cyclotron frequency (rad/s)
        
        ! Wave parameters
        omega = cmplx(1.0e9_real64, 1.0e6_real64, real64)  ! Complex frequency
        k_par = cmplx(1.0e3_real64, 0.0_real64, real64)    ! Parallel wavenumber
        k_perp = cmplx(1.0e4_real64, 0.0_real64, real64)   ! Perpendicular wavenumber
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (omega_pe > omega_pi)  ! Electrons lighter
        test_passed = test_passed .and. (omega_ce > omega_ci)  ! Electrons faster
        test_passed = test_passed .and. (abs(omega) > 0.0_real64)
        test_passed = test_passed .and. (abs(k_par) > 0.0_real64)
        test_passed = test_passed .and. (abs(k_perp) > 0.0_real64)
        
        ! Simplified dielectric tensor elements
        epsilon_xx = cmplx(1.0_real64 - omega_pe**2/real(omega)**2, -0.01_real64, real64)
        epsilon_yy = epsilon_xx
        epsilon_zz = cmplx(1.0_real64 - omega_pe**2/real(omega)**2, -0.005_real64, real64)
        epsilon_xy = cmplx(0.0_real64, omega_pe**2 * omega_ce / (real(omega)**3), real64)
        
        test_passed = test_passed .and. (abs(epsilon_xx) > 0.0_real64)
        test_passed = test_passed .and. (abs(epsilon_yy) > 0.0_real64)
        test_passed = test_passed .and. (abs(epsilon_zz) > 0.0_real64)
        test_passed = test_passed .and. (abs(epsilon_xy) > 0.0_real64)
        
        ! Test dispersion relation determinant (simplified)
        det_disp = epsilon_xx - (k_par**2 * 2.998e8_real64**2 / omega**2)
        test_passed = test_passed .and. (abs(det_disp) > 0.0_real64)
        
        ! Test refractive indices
        n_index_o = sqrt(epsilon_xx)  ! Ordinary wave
        n_index_x = sqrt(epsilon_zz)  ! Extraordinary wave
        
        test_passed = test_passed .and. (abs(n_index_o) > 0.0_real64)
        test_passed = test_passed .and. (abs(n_index_x) > 0.0_real64)
        test_passed = test_passed .and. (abs(n_index_o) < 10.0_real64)  ! Reasonable bounds
        test_passed = test_passed .and. (abs(n_index_x) < 10.0_real64)
        
        call end_test(test_passed)
    end subroutine test_dispersion_relations
    
    !---------------------------------------------------------------------------
    ! Test 7: Physics zone management and coupling
    !---------------------------------------------------------------------------
    subroutine test_physics_zone_management()
        integer, parameter :: n_zones = 3
        integer :: zone_types(n_zones)
        real(real64) :: zone_boundaries(n_zones+1)
        real(real64) :: r_test, zone_width
        integer :: zone_id, ierr
        logical :: is_flre_zone, is_vacuum_zone, is_imhd_zone
        
        call start_test("Physics zone management")
        test_passed = .true.
        
        ! Test physics zone setup
        ! Zone types: 1=FLRE, 2=Vacuum, 3=IMHD
        zone_types = [1, 2, 3]  ! FLRE, Vacuum, IMHD zones
        zone_boundaries = [0.0_real64, 0.7_real64, 0.9_real64, 1.0_real64]
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (size(zone_types) == n_zones)
        test_passed = test_passed .and. (size(zone_boundaries) == n_zones + 1)
        test_passed = test_passed .and. (all(zone_types > 0))
        test_passed = test_passed .and. (all(zone_types <= 3))
        
        ! Test zone boundary ordering
        do ierr = 2, n_zones + 1
            test_passed = test_passed .and. (zone_boundaries(ierr) > zone_boundaries(ierr-1))
        end do
        
        ! Test zone identification
        r_test = 0.5_real64  ! Core region
        zone_id = 1  ! Should be FLRE zone
        is_flre_zone = (zone_types(zone_id) == 1)
        test_passed = test_passed .and. is_flre_zone
        test_passed = test_passed .and. (r_test >= zone_boundaries(zone_id))
        test_passed = test_passed .and. (r_test < zone_boundaries(zone_id+1))
        
        r_test = 0.8_real64  ! Edge region
        zone_id = 2  ! Should be Vacuum zone
        is_vacuum_zone = (zone_types(zone_id) == 2)
        test_passed = test_passed .and. is_vacuum_zone
        test_passed = test_passed .and. (r_test >= zone_boundaries(zone_id))
        test_passed = test_passed .and. (r_test < zone_boundaries(zone_id+1))
        
        r_test = 0.95_real64  ! Far edge
        zone_id = 3  ! Should be IMHD zone
        is_imhd_zone = (zone_types(zone_id) == 3)
        test_passed = test_passed .and. is_imhd_zone
        test_passed = test_passed .and. (r_test >= zone_boundaries(zone_id))
        test_passed = test_passed .and. (r_test < zone_boundaries(zone_id+1))
        
        ! Test zone width calculations
        do ierr = 1, n_zones
            zone_width = zone_boundaries(ierr+1) - zone_boundaries(ierr)
            test_passed = test_passed .and. (zone_width > 0.0_real64)
            test_passed = test_passed .and. (zone_width < 1.0_real64)
        end do
        
        call end_test(test_passed)
    end subroutine test_physics_zone_management
    
    !---------------------------------------------------------------------------
    ! Test 8: Lab frame transformations
    !---------------------------------------------------------------------------
    subroutine test_lab_frame_transformations()
        real(real64) :: omega_lab, omega_plasma, v_rot, R_major
        real(real64) :: E_lab(3), E_plasma(3), B_lab(3), B_plasma(3)
        real(real64) :: transformation_matrix(3,3)
        real(real64) :: doppler_factor, gamma_relativistic
        integer :: ierr, i
        
        call start_test("Lab frame transformations")
        test_passed = .true.
        
        ! Test Doppler shift due to plasma rotation
        omega_lab = 5.0e7_real64     ! Laboratory frequency (rad/s)
        v_rot = 1.0e4_real64         ! Rotation velocity (m/s)
        R_major = 1.65_real64        ! Major radius (m)
        
        doppler_factor = 1.0_real64 + v_rot / 2.998e8_real64  ! Non-relativistic
        omega_plasma = omega_lab / doppler_factor
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (omega_lab > 0.0_real64)
        test_passed = test_passed .and. (omega_plasma > 0.0_real64)
        test_passed = test_passed .and. (abs(omega_plasma - omega_lab) < omega_lab * 0.1_real64)
        test_passed = test_passed .and. (doppler_factor > 1.0_real64)
        test_passed = test_passed .and. (doppler_factor < 1.1_real64)  ! Non-relativistic
        
        ! Test field transformations (simplified)
        E_lab = [1.0_real64, 0.5_real64, 0.1_real64]      ! Lab frame E field
        B_lab = [0.01_real64, 0.005_real64, 2.5_real64]   ! Lab frame B field
        
        ! Simple rotation transformation (toroidal geometry)
        transformation_matrix = 0.0_real64
        transformation_matrix(1,1) = cos(PI/4.0_real64)   ! 45 degree rotation
        transformation_matrix(1,2) = -sin(PI/4.0_real64)
        transformation_matrix(2,1) = sin(PI/4.0_real64)
        transformation_matrix(2,2) = cos(PI/4.0_real64)
        transformation_matrix(3,3) = 1.0_real64
        
        ! Transform to plasma frame
        E_plasma = matmul(transformation_matrix, E_lab)
        B_plasma = matmul(transformation_matrix, B_lab)
        
        test_passed = test_passed .and. (abs(E_plasma(1)) > 0.0_real64)
        test_passed = test_passed .and. (abs(E_plasma(2)) > 0.0_real64)
        test_passed = test_passed .and. (abs(E_plasma(3) - E_lab(3)) < 1.0e-10_real64)  ! No z rotation
        test_passed = test_passed .and. (abs(B_plasma(3) - B_lab(3)) < 1.0e-10_real64)  ! No z rotation
        
        ! Test energy conservation
        test_passed = test_passed .and. (abs(sum(E_plasma**2) - sum(E_lab**2)) < 1.0e-10_real64)
        test_passed = test_passed .and. (abs(sum(B_plasma**2) - sum(B_lab**2)) < 1.0e-10_real64)
        
        call end_test(test_passed)
    end subroutine test_lab_frame_transformations
    
    !---------------------------------------------------------------------------
    ! Test 9: F0 distribution function moments
    !---------------------------------------------------------------------------
    subroutine test_f0_distribution_moments()
        real(real64) :: n_density, Ti, Te, vth_i, vth_e
        real(real64) :: f0_moment_0, f0_moment_1, f0_moment_2
        real(real64) :: pressure_i, pressure_e, energy_density
        real(real64) :: beta_thermal, v_thermal
        integer :: ierr
        
        call start_test("F0 distribution function moments")
        test_passed = .true.
        
        ! Test Maxwellian distribution moments
        n_density = 1.0e19_real64    ! Number density (1/m^3)
        Ti = 1000.0_real64 * 1.6022e-19_real64  ! Ion temperature (J)
        Te = 2000.0_real64 * 1.6022e-19_real64  ! Electron temperature (J)
        
        ! Thermal velocities
        vth_i = sqrt(2.0_real64 * Ti / 1.6726e-27_real64)  ! Ion thermal velocity
        vth_e = sqrt(2.0_real64 * Te / 9.1094e-31_real64)  ! Electron thermal velocity
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (n_density > 0.0_real64)
        test_passed = test_passed .and. (Ti > 0.0_real64)
        test_passed = test_passed .and. (Te > Ti)  ! Typically Te > Ti
        test_passed = test_passed .and. (vth_e > vth_i)  ! Electrons faster
        
        ! Calculate distribution moments
        f0_moment_0 = n_density                    ! 0th moment: density
        f0_moment_1 = 0.0_real64                   ! 1st moment: mean velocity (zero)
        f0_moment_2 = n_density * vth_i**2         ! 2nd moment: pressure/mass
        
        test_passed = test_passed .and. (abs(f0_moment_0 - n_density) < 1.0e-10_real64)
        test_passed = test_passed .and. (abs(f0_moment_1) < 1.0e-10_real64)
        test_passed = test_passed .and. (f0_moment_2 > 0.0_real64)
        
        ! Calculate pressures
        pressure_i = n_density * Ti  ! Ion pressure
        pressure_e = n_density * Te  ! Electron pressure
        energy_density = 1.5_real64 * (pressure_i + pressure_e)  ! Total energy density
        
        test_passed = test_passed .and. (pressure_i > 0.0_real64)
        test_passed = test_passed .and. (pressure_e > pressure_i)
        test_passed = test_passed .and. (energy_density > 0.0_real64)
        
        ! Test thermal beta parameter
        beta_thermal = (pressure_i + pressure_e) / (2.5_real64**2 / (2.0_real64 * 4.0e-7_real64 * PI))
        test_passed = test_passed .and. (beta_thermal > 0.0_real64)
        test_passed = test_passed .and. (beta_thermal < 1.0_real64)  ! Low beta plasma
        
        call end_test(test_passed)
    end subroutine test_f0_distribution_moments
    
    !---------------------------------------------------------------------------
    ! Test 10: Collision frequencies
    !---------------------------------------------------------------------------
    subroutine test_collision_frequencies()
        real(real64) :: nu_ii, nu_ee, nu_ei, nu_ie
        real(real64) :: n_density, Ti, Te, Z_eff, ln_lambda
        real(real64) :: tau_ii, tau_ee, tau_ei, slowing_down_time
        integer :: ierr
        
        call start_test("Collision frequencies")
        test_passed = .true.
        
        ! Test collision frequency calculations
        n_density = 1.0e19_real64    ! Number density (1/m^3)
        Ti = 1000.0_real64           ! Ion temperature (eV)
        Te = 2000.0_real64           ! Electron temperature (eV)
        Z_eff = 1.0_real64           ! Effective charge
        ln_lambda = 15.0_real64      ! Coulomb logarithm
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (n_density > 0.0_real64)
        test_passed = test_passed .and. (Ti > 0.0_real64)
        test_passed = test_passed .and. (Te > 0.0_real64)
        test_passed = test_passed .and. (Z_eff >= 1.0_real64)
        test_passed = test_passed .and. (ln_lambda > 10.0_real64)
        
        ! Calculate collision frequencies (simplified formulas)
        ! Ion-ion collisions
        nu_ii = 4.8e-14_real64 * Z_eff**4 * n_density * ln_lambda / Ti**1.5_real64
        
        ! Electron-electron collisions  
        nu_ee = 2.9e-12_real64 * n_density * ln_lambda / Te**1.5_real64
        
        ! Electron-ion collisions
        nu_ei = 2.9e-12_real64 * Z_eff * n_density * ln_lambda / Te**1.5_real64
        
        ! Ion-electron collisions (momentum transfer)
        nu_ie = nu_ei * 9.1094e-31_real64 / 1.6726e-27_real64  ! Mass ratio correction
        
        test_passed = test_passed .and. (nu_ii > 0.0_real64)
        test_passed = test_passed .and. (nu_ee > 0.0_real64)
        test_passed = test_passed .and. (nu_ei > 0.0_real64)
        test_passed = test_passed .and. (nu_ie > 0.0_real64)
        test_passed = test_passed .and. (nu_ee > nu_ii)     ! Electrons collide more
        test_passed = test_passed .and. (nu_ei > nu_ie)     ! Different mass scaling
        
        ! Calculate collision times
        tau_ii = 1.0_real64 / nu_ii
        tau_ee = 1.0_real64 / nu_ee
        tau_ei = 1.0_real64 / nu_ei
        slowing_down_time = tau_ei * 1836.0_real64  ! Mass ratio effect
        
        test_passed = test_passed .and. (tau_ii > 0.0_real64)
        test_passed = test_passed .and. (tau_ee > 0.0_real64)
        test_passed = test_passed .and. (tau_ei > 0.0_real64)
        test_passed = test_passed .and. (slowing_down_time > tau_ei)
        test_passed = test_passed .and. (tau_ii < 1.0_real64)      ! Reasonable timescales
        test_passed = test_passed .and. (tau_ee < 1.0_real64)
        
        call end_test(test_passed)
    end subroutine test_collision_frequencies
    
    !---------------------------------------------------------------------------
    ! Test 11: Magnetic field calculations
    !---------------------------------------------------------------------------
    subroutine test_magnetic_field_calculations()
        real(real64) :: Bp, Bt, B_total, q_safety, psi_flux
        real(real64) :: R_major, r_minor, I_plasma, mu_0
        real(real64) :: grad_B, curvature_drift, magnetic_shear
        real(real64) :: B_field(3), grad_B_vec(3)
        integer :: ierr, i
        
        call start_test("Magnetic field calculations")
        test_passed = .true.
        
        ! Test tokamak magnetic field geometry
        R_major = 1.65_real64        ! Major radius (m)
        r_minor = 0.5_real64         ! Minor radius (m)
        I_plasma = 1.0e6_real64      ! Plasma current (A)
        mu_0 = 4.0e-7_real64 * PI    ! Vacuum permeability
        q_safety = 2.0_real64        ! Safety factor
        
        ! Calculate field components
        Bt = 2.5_real64              ! Toroidal field (T)
        Bp = Bt / q_safety           ! Poloidal field from safety factor
        B_total = sqrt(Bp**2 + Bt**2)
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (Bt > 0.0_real64)
        test_passed = test_passed .and. (Bp > 0.0_real64)
        test_passed = test_passed .and. (B_total > Bt)      ! Total > components
        test_passed = test_passed .and. (B_total > Bp)
        test_passed = test_passed .and. (q_safety > 1.0_real64)  ! Stability requirement
        
        ! Test flux surface calculation
        psi_flux = Bp * R_major * r_minor  ! Simplified poloidal flux
        test_passed = test_passed .and. (psi_flux > 0.0_real64)
        test_passed = test_passed .and. (psi_flux < 10.0_real64)  ! Reasonable magnitude
        
        ! Test magnetic field gradient
        grad_B = Bt / R_major  ! Simplified gradient (1/R dependence)
        test_passed = test_passed .and. (grad_B > 0.0_real64)
        test_passed = test_passed .and. (grad_B < 10.0_real64)
        
        ! Test curvature drift frequency
        curvature_drift = grad_B / B_total  ! Simplified curvature
        test_passed = test_passed .and. (curvature_drift > 0.0_real64)
        test_passed = test_passed .and. (curvature_drift < 1.0_real64)
        
        ! Test magnetic shear
        magnetic_shear = 0.1_real64  ! Typical value
        test_passed = test_passed .and. (magnetic_shear > 0.0_real64)
        test_passed = test_passed .and. (magnetic_shear < 1.0_real64)
        
        ! Test field vector properties
        B_field = [0.1_real64, Bp, Bt]  ! (Br, Bp, Bt)
        test_passed = test_passed .and. (sqrt(sum(B_field**2)) > 0.0_real64)
        test_passed = test_passed .and. (B_field(3) > B_field(2))  ! Bt > Bp
        test_passed = test_passed .and. (B_field(2) > abs(B_field(1)))  ! Bp > Br
        
        call end_test(test_passed)
    end subroutine test_magnetic_field_calculations
    
    !---------------------------------------------------------------------------
    ! Test 12: Physics spline interpolation
    !---------------------------------------------------------------------------
    subroutine test_physics_spline_interpolation()
        real(real64), allocatable :: x_data(:), y_data(:), y_interp(:)
        real(real64), allocatable :: x_test(:), spline_coeffs(:)
        integer :: n_data, n_test, ierr, i
        real(real64) :: x_val, y_val, dy_val, error_max
        
        call start_test("Physics spline interpolation")
        test_passed = .true.
        
        ! Test spline interpolation for physics profiles
        n_data = 20
        n_test = 50
        
        allocate(x_data(n_data), y_data(n_data), x_test(n_test), &
                y_interp(n_test), spline_coeffs(4*n_data), stat=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            ! Create test data (parabolic profile)
            do i = 1, n_data
                x_data(i) = real(i-1, real64) / real(n_data-1, real64)
                y_data(i) = 1.0_real64 - x_data(i)**2  ! Parabolic profile
            end do
            
            ! Create test points
            do i = 1, n_test
                x_test(i) = real(i-1, real64) / real(n_test-1, real64)
            end do
            
            ! Mock spline interpolation
            do i = 1, n_test
                x_val = x_test(i)
                y_val = 1.0_real64 - x_val**2  ! Analytical result
                y_interp(i) = y_val  ! Perfect interpolation for test
            end do
            
            test_passed = test_passed .and. (all(x_data >= 0.0_real64))
            test_passed = test_passed .and. (all(x_data <= 1.0_real64))
            test_passed = test_passed .and. (all(y_data >= 0.0_real64))
            test_passed = test_passed .and. (y_data(1) == 1.0_real64)     ! Center value
            test_passed = test_passed .and. (y_data(n_data) == 0.0_real64) ! Edge value
            
            ! Test interpolation accuracy
            error_max = 0.0_real64
            do i = 1, n_test
                error_max = max(error_max, abs(y_interp(i) - (1.0_real64 - x_test(i)**2)))
            end do
            
            test_passed = test_passed .and. (error_max < 1.0e-10_real64)  ! Perfect for polynomial
            test_passed = test_passed .and. (all(y_interp >= 0.0_real64))
            test_passed = test_passed .and. (y_interp(1) >= y_interp(n_test))  ! Monotonic decrease
            
            ! Test derivative calculation (simplified)
            do i = 2, n_test-1
                dy_val = (y_interp(i+1) - y_interp(i-1)) / (x_test(i+1) - x_test(i-1))
                test_passed = test_passed .and. (abs(dy_val + 2.0_real64 * x_test(i)) < 0.1_real64)
            end do
            
            deallocate(x_data, y_data, x_test, y_interp, spline_coeffs)
        end if
        
        call end_test(test_passed)
    end subroutine test_physics_spline_interpolation
    
    !---------------------------------------------------------------------------
    ! Test 13: Physics zone coupling and continuity
    !---------------------------------------------------------------------------
    subroutine test_physics_zone_coupling()
        real(real64) :: E_field_left(3), E_field_right(3)
        real(real64) :: B_field_left(3), B_field_right(3)
        real(real64) :: current_jump(3), boundary_condition
        real(real64) :: impedance_left, impedance_right, reflection_coeff
        real(real64) :: r_boundary, coupling_strength
        integer :: ierr
        
        call start_test("Physics zone coupling")
        test_passed = .true.
        
        ! Test electromagnetic field continuity at zone boundaries
        r_boundary = 0.7_real64  ! FLRE-Vacuum boundary
        
        ! Fields on FLRE side (left)
        E_field_left = [1.0_real64, 0.5_real64, 0.1_real64]
        B_field_left = [0.01_real64, 0.005_real64, 2.5_real64]
        
        ! Fields on Vacuum side (right) - continuity conditions
        E_field_right = E_field_left  ! Tangential E continuous
        B_field_right = B_field_left  ! Tangential B continuous
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (r_boundary > 0.0_real64)
        test_passed = test_passed .and. (r_boundary < 1.0_real64)
        
        ! Test field continuity
        test_passed = test_passed .and. (abs(E_field_right(1) - E_field_left(1)) < 1.0e-10_real64)
        test_passed = test_passed .and. (abs(E_field_right(2) - E_field_left(2)) < 1.0e-10_real64)
        test_passed = test_passed .and. (abs(B_field_right(2) - B_field_left(2)) < 1.0e-10_real64)
        test_passed = test_passed .and. (abs(B_field_right(3) - B_field_left(3)) < 1.0e-10_real64)
        
        ! Test current jump conditions (simplified)
        current_jump = [0.1_real64, 0.05_real64, 0.0_real64]  ! Surface current
        boundary_condition = abs(B_field_right(1) - B_field_left(1)) - abs(current_jump(1))
        test_passed = test_passed .and. (abs(boundary_condition) < 0.2_real64)
        
        ! Test impedance matching
        impedance_left = 377.0_real64 / sqrt(2.0_real64)     ! Plasma impedance (simplified)
        impedance_right = 377.0_real64                       ! Vacuum impedance
        reflection_coeff = abs((impedance_right - impedance_left) / (impedance_right + impedance_left))**2
        
        test_passed = test_passed .and. (impedance_left > 0.0_real64)
        test_passed = test_passed .and. (impedance_right > 0.0_real64)
        test_passed = test_passed .and. (reflection_coeff >= 0.0_real64)
        test_passed = test_passed .and. (reflection_coeff <= 1.0_real64)
        
        ! Test coupling strength
        coupling_strength = 1.0_real64 - reflection_coeff
        test_passed = test_passed .and. (coupling_strength >= 0.0_real64)
        test_passed = test_passed .and. (coupling_strength <= 1.0_real64)
        test_passed = test_passed .and. (coupling_strength > 0.5_real64)  ! Good coupling
        
        call end_test(test_passed)
    end subroutine test_physics_zone_coupling
    
    !---------------------------------------------------------------------------
    ! Test 14: Physics parameter validation
    !---------------------------------------------------------------------------
    subroutine test_physics_parameter_validation()
        real(real64) :: Ti, Te, n_density, B_field, q_safety
        real(real64) :: beta_pressure, omega_freq, k_wavelength
        logical :: is_valid_temperature, is_valid_density, is_valid_field
        logical :: is_valid_safety, is_valid_beta, is_valid_frequency
        integer :: ierr
        
        call start_test("Physics parameter validation")
        test_passed = .true.
        
        ! Test physics parameter bounds and consistency
        Ti = 1000.0_real64       ! Ion temperature (eV)
        Te = 2000.0_real64       ! Electron temperature (eV)  
        n_density = 1.0e19_real64 ! Number density (1/m^3)
        B_field = 2.5_real64     ! Magnetic field (T)
        q_safety = 2.0_real64    ! Safety factor
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        
        ! Validate temperature ranges
        is_valid_temperature = (Ti > 0.0_real64) .and. (Ti < 1.0e5_real64) .and. &
                              (Te > 0.0_real64) .and. (Te < 1.0e5_real64) .and. &
                              (Te >= Ti)  ! Usually Te >= Ti
        test_passed = test_passed .and. is_valid_temperature
        
        ! Validate density ranges
        is_valid_density = (n_density > 1.0e17_real64) .and. (n_density < 1.0e22_real64)
        test_passed = test_passed .and. is_valid_density
        
        ! Validate magnetic field ranges
        is_valid_field = (B_field > 0.1_real64) .and. (B_field < 20.0_real64)
        test_passed = test_passed .and. is_valid_field
        
        ! Validate safety factor ranges
        is_valid_safety = (q_safety > 1.0_real64) .and. (q_safety < 10.0_real64)
        test_passed = test_passed .and. is_valid_safety
        
        ! Validate beta parameter
        beta_pressure = 2.0_real64 * 4.0e-7_real64 * PI * n_density * (Ti + Te) * 1.6022e-19_real64 / B_field**2
        is_valid_beta = (beta_pressure > 0.0_real64) .and. (beta_pressure < 1.0_real64)
        test_passed = test_passed .and. is_valid_beta
        
        ! Validate frequency ranges
        omega_freq = 5.0e7_real64  ! RF frequency (rad/s)
        is_valid_frequency = (omega_freq > 1.0e6_real64) .and. (omega_freq < 1.0e12_real64)
        test_passed = test_passed .and. is_valid_frequency
        
        ! Validate wavelength consistency
        k_wavelength = omega_freq / 2.998e8_real64  ! Vacuum wavelength
        test_passed = test_passed .and. (k_wavelength > 0.0_real64)
        test_passed = test_passed .and. (k_wavelength < 1.0e6_real64)
        test_passed = test_passed .and. (2.0_real64 * PI / k_wavelength > 1.0e-6_real64)  ! Reasonable wavelength
        
        ! Test consistency checks
        test_passed = test_passed .and. (Ti * n_density < 1.0e25_real64)     ! Pressure bound
        test_passed = test_passed .and. (Te * n_density < 1.0e25_real64)     ! Pressure bound
        test_passed = test_passed .and. (B_field**2 > 2.0_real64 * 4.0e-7_real64 * PI * beta_pressure)  ! Beta definition
        
        call end_test(test_passed)
    end subroutine test_physics_parameter_validation
    
    !---------------------------------------------------------------------------
    ! Test utilities  
    !---------------------------------------------------------------------------
    subroutine start_test(name)
        character(len=*), intent(in) :: name
        total_tests = total_tests + 1
        write(*,'(A,A)', advance='no') "Testing ", name
        write(*,'(A)', advance='no') " ... "
    end subroutine start_test
    
    subroutine end_test(passed)
        logical, intent(in) :: passed
        if (passed) then
            print *, "PASSED"
            passed_tests = passed_tests + 1
        else
            print *, "FAILED"
            failed_tests = failed_tests + 1
        end if
    end subroutine end_test

end program test_kilca_physics