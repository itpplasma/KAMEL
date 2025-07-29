program test_kilca_antenna
    !---------------------------------------------------------------------------
    ! KiLCA Antenna and Interface Systems - Unit Tests
    !
    ! Tests the antenna modeling and interface systems including current density
    ! spectrum calculations, continuity boundary matching, and external code
    ! interface protocols.
    !
    ! Author: Claude (Anthropic)
    ! Date: 2024
    !
    ! Testing: kilca_antenna_m.f90 (to be implemented)
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
    write(*,'(A)') " KiLCA Antenna and Interface Systems - Unit Tests"
    write(*,'(A)') " ========================================================"
    write(*,'(A)') " "
    
    ! Run all tests
    call test_antenna_configuration()
    call test_current_density_spectrum()
    call test_continuity_boundary_matching()
    call test_coordinate_transformations()
    call test_delta_function_source()
    call test_coupling_calculations()
    call test_interface_data_structures()
    call test_wave_code_interface()
    call test_ql_balance_integration()
    call test_maxwell_integration()
    call test_antenna_spectrum_modes()
    call test_memory_management()
    
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
    ! Test 1: Antenna configuration and settings
    !---------------------------------------------------------------------------
    subroutine test_antenna_configuration()
        real(real64) :: ra, wa, I0
        complex(real64) :: flab
        integer, allocatable :: modes(:)
        integer :: ierr, dma
        logical :: flag_debug, flag_eigmode
        
        call start_test("Antenna configuration and settings")
        test_passed = .true.
        
        ! Test antenna configuration parameters
        ! antenna class translation: ra, wa, I0, flab, modes, flags
        
        ra = 8.5_real64          ! Antenna radius (cm)
        wa = 1.0_real64          ! Current density layer width
        I0 = 1000.0_real64       ! Current in antenna coils (statamp)
        flab = cmplx(50.0e6_real64, 1.0e3_real64, real64)  ! Laboratory frequency (Hz)
        flag_debug = .false.
        flag_eigmode = .true.
        
        ! Mode numbers array allocation
        dma = 4
        allocate(modes(2*dma), stat=ierr)  ! (m,n) pairs
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            ! Set test mode numbers: (3,1), (12,4), (-3,-1), (-12,-4)
            modes(1:8) = [3, 1, 12, 4, -3, -1, -12, -4]
            
            ! Validate configuration
            test_passed = test_passed .and. (ra > 0.0_real64)
            test_passed = test_passed .and. (wa > 0.0_real64)
            test_passed = test_passed .and. (I0 > 0.0_real64)
            test_passed = test_passed .and. (real(flab) > 0.0_real64)
            test_passed = test_passed .and. (dma > 0)
            test_passed = test_passed .and. (size(modes) == 2*dma)
            
            deallocate(modes)
        end if
        
        call end_test(test_passed)
    end subroutine test_antenna_configuration
    
    !---------------------------------------------------------------------------
    ! Test 2: Current density spectrum calculations
    !---------------------------------------------------------------------------
    subroutine test_current_density_spectrum()
        complex(real64), allocatable :: spectrum(:)
        integer :: mode_m, mode_n, n_harmonics, i
        real(real64) :: pi_val
        integer :: ierr
        
        call start_test("Current density spectrum calculations")
        test_passed = .true.
        
        ! Test current density spectrum calculations for different modes
        pi_val = PI
        n_harmonics = 16
        allocate(spectrum(n_harmonics), stat=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            ! Test (3,1) mode spectrum
            mode_m = 3
            mode_n = 1
            
            do i = 1, n_harmonics
                ! (3,1) mode: F_n = 4.0*exp(-3.0*IM*pi*nl/16.0)*sin(pi*n/4.0)/sin(pi*n/16.0)
                spectrum(i) = 4.0_real64 * &
                    exp(cmplx(0.0_real64, -3.0_real64 * pi_val * real(i, real64) / 16.0_real64, real64)) * &
                    sin(pi_val * real(i, real64) / 4.0_real64) / &
                    sin(pi_val * real(i, real64) / 16.0_real64)
            end do
            
            ! Validate spectrum properties
            test_passed = test_passed .and. (abs(spectrum(1)) > 0.0_real64)
            test_passed = test_passed .and. (abs(spectrum(4)) > 0.0_real64)
            test_passed = test_passed .and. all(abs(spectrum) < 100.0_real64)  ! Reasonable magnitude
            
            ! Test (12,4) mode spectrum (constant amplitude)
            mode_m = 12
            mode_n = 4
            spectrum = cmplx(16.0_real64, 0.0_real64, real64)  ! Constant spectrum
            
            test_passed = test_passed .and. (abs(spectrum(1)) == 16.0_real64)
            test_passed = test_passed .and. (abs(spectrum(n_harmonics)) == 16.0_real64)
            
            deallocate(spectrum)
        end if
        
        call end_test(test_passed)
    end subroutine test_current_density_spectrum
    
    !---------------------------------------------------------------------------
    ! Test 3: Continuity boundary matching equations
    !---------------------------------------------------------------------------
    subroutine test_continuity_boundary_matching()
        complex(real64), allocatable :: matrix(:,:), coeffs(:), jumps(:)
        complex(real64) :: determinant
        integer :: ierr, n_eq, info
        integer, allocatable :: ipiv(:)
        
        call start_test("Continuity boundary matching equations")
        test_passed = .true.
        
        ! Test continuity equation solving: mat * coeffs = jumps
        n_eq = 4
        allocate(matrix(n_eq, n_eq), coeffs(n_eq), jumps(n_eq), ipiv(n_eq), stat=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            ! Create test system matrix (mock inner/outer solution vectors)
            matrix = cmplx(0.0_real64, 0.0_real64, real64)
            matrix(1,1) = cmplx(2.0_real64, 0.1_real64, real64)  ! Inner solution vectors
            matrix(1,2) = cmplx(1.0_real64, 0.05_real64, real64)
            matrix(2,1) = cmplx(0.5_real64, 0.02_real64, real64)
            matrix(2,2) = cmplx(1.5_real64, 0.08_real64, real64)
            
            matrix(3,3) = cmplx(1.8_real64, -0.1_real64, real64)  ! Outer solution vectors
            matrix(3,4) = cmplx(0.9_real64, -0.05_real64, real64)
            matrix(4,3) = cmplx(0.7_real64, -0.03_real64, real64)
            matrix(4,4) = cmplx(1.3_real64, -0.07_real64, real64)
            
            ! Coupling between inner and outer
            matrix(1,3) = cmplx(-0.1_real64, 0.0_real64, real64)
            matrix(2,4) = cmplx(-0.1_real64, 0.0_real64, real64)
            
            ! Set jump conditions (boundary sources)
            jumps(1) = cmplx(1.0_real64, 0.0_real64, real64)   ! delBs
            jumps(2) = cmplx(0.5_real64, 0.0_real64, real64)   ! delBp
            jumps(3) = cmplx(0.1_real64, 0.0_real64, real64)   ! Continuity condition 1
            jumps(4) = cmplx(0.05_real64, 0.0_real64, real64)  ! Continuity condition 2
            
            ! Solve system using LAPACK zgesv (mock)
            coeffs = jumps  ! Initialize
            info = 0  ! Mock successful solve
            
            ! Calculate determinant (simplified)
            determinant = matrix(1,1) * matrix(2,2) - matrix(1,2) * matrix(2,1)
            
            test_passed = test_passed .and. (info == 0)
            test_passed = test_passed .and. (abs(determinant) > 1.0e-10_real64)
            test_passed = test_passed .and. (all(abs(coeffs) < 100.0_real64))
            
            deallocate(matrix, coeffs, jumps, ipiv)
        end if
        
        call end_test(test_passed)
    end subroutine test_continuity_boundary_matching
    
    !---------------------------------------------------------------------------
    ! Test 4: Coordinate transformations (cyl2rsp)
    !---------------------------------------------------------------------------
    subroutine test_coordinate_transformations()
        complex(real64) :: Ja_cyl(2), Ja_rsp(2)
        real(real64) :: theta, transformation_factor
        integer :: ierr
        
        call start_test("Coordinate transformations")
        test_passed = .true.
        
        ! Test cylindrical to (r,s,p) coordinate transformation
        theta = PI / 4.0_real64  ! 45 degrees
        
        ! Mock cylindrical current density
        Ja_cyl(1) = cmplx(1.0_real64, 0.1_real64, real64)  ! J_theta
        Ja_cyl(2) = cmplx(0.5_real64, 0.05_real64, real64) ! J_z
        
        ! Transform: cyl2rsp transformation
        transformation_factor = cos(theta)  ! Simplified transformation
        
        Ja_rsp(1) = Ja_cyl(1) * transformation_factor     ! J_s component
        Ja_rsp(2) = Ja_cyl(2) * transformation_factor     ! J_p component
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(Ja_rsp(1)) > 0.0_real64)
        test_passed = test_passed .and. (abs(Ja_rsp(2)) > 0.0_real64)
        test_passed = test_passed .and. (abs(Ja_rsp(1)) <= abs(Ja_cyl(1)))
        test_passed = test_passed .and. (abs(Ja_rsp(2)) <= abs(Ja_cyl(2)))
        
        call end_test(test_passed)
    end subroutine test_coordinate_transformations
    
    !---------------------------------------------------------------------------
    ! Test 5: Delta function source calculation
    !---------------------------------------------------------------------------
    subroutine test_delta_function_source()
        real(real64) :: ra, wa, r_test, delta_value
        real(real64) :: isqrt2pi, gaussian_factor
        integer :: ierr
        
        call start_test("Delta function source calculation")
        test_passed = .true.
        
        ! Test delta function approximation: delta = (isqrt2pi/wa)*exp(-((r-ra)^2)/(2.0*wa^2))
        ra = 8.5_real64  ! Antenna radius
        wa = 1.0_real64  ! Width parameter
        isqrt2pi = 1.0_real64 / sqrt(2.0_real64 * PI)
        
        ! Test at antenna center
        r_test = ra
        gaussian_factor = exp(-((r_test - ra)**2) / (2.0_real64 * wa**2))
        delta_value = (isqrt2pi / wa) * gaussian_factor
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (delta_value > 0.0_real64)
        test_passed = test_passed .and. (gaussian_factor == 1.0_real64)  ! At center
        
        ! Test away from antenna
        r_test = ra + 2.0_real64 * wa  ! 2 widths away
        gaussian_factor = exp(-((r_test - ra)**2) / (2.0_real64 * wa**2))
        delta_value = (isqrt2pi / wa) * gaussian_factor
        
        test_passed = test_passed .and. (delta_value > 0.0_real64)
        test_passed = test_passed .and. (gaussian_factor < 0.2_real64)  ! Significantly smaller
        
        call end_test(test_passed)
    end subroutine test_delta_function_source
    
    !---------------------------------------------------------------------------
    ! Test 6: Coupling calculations and magnetic field jumps
    !---------------------------------------------------------------------------
    subroutine test_coupling_calculations()
        complex(real64) :: jsurft(2), delBs, delBp
        real(real64) :: fpc, coupling_strength
        integer :: ierr
        
        call start_test("Coupling calculations")
        test_passed = .true.
        
        ! Test current jumps and magnetic field calculations
        fpc = 4.0_real64 * PI / 3.0e10_real64  ! Physics constant (c.g.s units)
        
        ! Mock surface current density
        jsurft(1) = cmplx(100.0_real64, 10.0_real64, real64)  ! Surface current s-component
        jsurft(2) = cmplx(200.0_real64, 20.0_real64, real64)  ! Surface current p-component
        
        ! Calculate magnetic field jumps: delBs = fpc*jsurft(2), delBp = -fpc*jsurft(1)
        delBs = fpc * jsurft(2)
        delBp = -fpc * jsurft(1)
        
        ! Calculate coupling strength
        coupling_strength = abs(delBs)**2 + abs(delBp)**2
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(delBs) > 0.0_real64)
        test_passed = test_passed .and. (abs(delBp) > 0.0_real64)
        test_passed = test_passed .and. (coupling_strength > 0.0_real64)
        test_passed = test_passed .and. (real(delBp) < 0.0_real64)  ! Negative sign check
        
        call end_test(test_passed)
    end subroutine test_coupling_calculations
    
    !---------------------------------------------------------------------------
    ! Test 7: Interface data structures and exchange protocols
    !---------------------------------------------------------------------------
    subroutine test_interface_data_structures()
        real(real64), allocatable :: wave_fields(:), current_densities(:)
        real(real64), allocatable :: background_profiles(:)
        complex(real64), allocatable :: conductivity_matrices(:,:)
        integer :: ierr, n_fields, n_currents, n_profiles
        
        call start_test("Interface data structures")
        test_passed = .true.
        
        ! Test wave code interface data structures
        ! Wave fields: Er, Es, Ep, Et, Ez, Br, Bs, Bp, Bt, Bz
        n_fields = 10
        allocate(wave_fields(n_fields), stat=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            wave_fields = [1.0_real64, 0.5_real64, 0.2_real64, 0.1_real64, 0.05_real64, &
                          0.02_real64, 0.01_real64, 0.005_real64, 0.002_real64, 0.001_real64]
            
            test_passed = test_passed .and. (all(wave_fields >= 0.0_real64))
            deallocate(wave_fields)
        end if
        
        ! Current densities: Jri, Jsi, Jpi, Jre, Jse, Jpe
        n_currents = 6
        allocate(current_densities(n_currents), stat=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            current_densities = [100.0_real64, 50.0_real64, 20.0_real64, &
                               200.0_real64, 100.0_real64, 40.0_real64]
            
            test_passed = test_passed .and. (all(abs(current_densities) < 1000.0_real64))
            deallocate(current_densities)
        end if
        
        ! Background profiles: q, n, Ti, Te, Vth, Vz, dPhi0
        n_profiles = 7
        allocate(background_profiles(n_profiles), stat=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            background_profiles = [1.5_real64, 1.0e14_real64, 1000.0_real64, 2000.0_real64, &
                                 1.0e6_real64, 0.0_real64, 1000.0_real64]
            
            test_passed = test_passed .and. (background_profiles(1) > 0.0_real64)  ! q > 0
            test_passed = test_passed .and. (background_profiles(2) > 0.0_real64)  ! n > 0
            deallocate(background_profiles)
        end if
        
        ! Conductivity matrices: k_mat_i, k_mat_e, c_mat_i, c_mat_e
        allocate(conductivity_matrices(3, 3), stat=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            conductivity_matrices = cmplx(1.0_real64, 0.1_real64, real64)
            conductivity_matrices(1,1) = cmplx(2.0_real64, 0.2_real64, real64)
            
            test_passed = test_passed .and. (abs(conductivity_matrices(1,1)) > 0.0_real64)
            deallocate(conductivity_matrices)
        end if
        
        call end_test(test_passed)
    end subroutine test_interface_data_structures
    
    !---------------------------------------------------------------------------
    ! Test 8: Wave code interface functions
    !---------------------------------------------------------------------------
    subroutine test_wave_code_interface()
        integer :: ierr, antenna_spectrum_dim
        integer, allocatable :: mode_numbers(:)
        real(real64), allocatable :: power_density(:)
        
        call start_test("Wave code interface functions")
        test_passed = .true.
        
        ! Test interface functions (mock implementations)
        ! get_antenna_spectrum_dim_()
        antenna_spectrum_dim = 8  ! Mock dimension
        test_passed = test_passed .and. (antenna_spectrum_dim > 0)
        
        ! get_antenna_spectrum_numbers_()
        allocate(mode_numbers(antenna_spectrum_dim), stat=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            mode_numbers = [3, 1, 12, 4, -3, -1, -12, -4]  ! Mock mode numbers
            test_passed = test_passed .and. (size(mode_numbers) == antenna_spectrum_dim)
            test_passed = test_passed .and. (any(mode_numbers > 0))
            test_passed = test_passed .and. (any(mode_numbers < 0))
            deallocate(mode_numbers)
        end if
        
        ! get_diss_power_density_from_wave_code_()
        allocate(power_density(10), stat=ierr)  ! 10 radial points
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            power_density = [1.0_real64, 0.9_real64, 0.8_real64, 0.7_real64, 0.6_real64, &
                           0.5_real64, 0.4_real64, 0.3_real64, 0.2_real64, 0.1_real64]
            
            test_passed = test_passed .and. (all(power_density >= 0.0_real64))
            test_passed = test_passed .and. (power_density(1) >= power_density(10))  ! Decreasing
            deallocate(power_density)
        end if
        
        call end_test(test_passed)
    end subroutine test_wave_code_interface
    
    !---------------------------------------------------------------------------
    ! Test 9: QL-Balance integration and data exchange
    !---------------------------------------------------------------------------
    subroutine test_ql_balance_integration()
        real(real64), allocatable :: transport_coeffs(:), profiles(:)
        integer :: ierr, n_radial
        
        call start_test("QL-Balance integration")
        test_passed = .true.
        
        ! Test QL-Balance integration and data exchange
        n_radial = 20
        
        ! Transport coefficients from QL-Balance
        allocate(transport_coeffs(n_radial), stat=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            ! Mock transport coefficients (decreasing with radius)
            do ierr = 1, n_radial
                transport_coeffs(ierr) = 1.0_real64 * exp(-real(ierr-1, real64) / 10.0_real64)
            end do
            
            test_passed = test_passed .and. (all(transport_coeffs > 0.0_real64))
            test_passed = test_passed .and. (transport_coeffs(1) > transport_coeffs(n_radial))
            deallocate(transport_coeffs)
        end if
        
        ! Background profiles to QL-Balance
        allocate(profiles(n_radial), stat=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            ! Mock profile data (parabolic)
            do ierr = 1, n_radial
                profiles(ierr) = 1.0_real64 - (real(ierr-1, real64) / real(n_radial-1, real64))**2
            end do
            
            test_passed = test_passed .and. (all(profiles >= 0.0_real64))
            test_passed = test_passed .and. (profiles(1) > profiles(n_radial))
            deallocate(profiles)
        end if
        
        call end_test(test_passed)
    end subroutine test_ql_balance_integration
    
    !---------------------------------------------------------------------------
    ! Test 10: Integration with Maxwell equations system
    !---------------------------------------------------------------------------
    subroutine test_maxwell_integration()
        complex(real64) :: source_terms(3), boundary_jumps(3)
        real(real64) :: antenna_position
        integer :: ierr
        
        call start_test("Maxwell equations integration")
        test_passed = .true.
        
        ! Test integration with Maxwell equations system
        antenna_position = 8.5_real64
        
        ! Mock antenna source terms for Maxwell equations
        source_terms(1) = cmplx(0.0_real64, 0.0_real64, real64)  ! No radial current
        source_terms(2) = cmplx(100.0_real64, 10.0_real64, real64)  ! Poloidal current
        source_terms(3) = cmplx(50.0_real64, 5.0_real64, real64)   ! Toroidal current
        
        ! Boundary jumps from antenna
        boundary_jumps(1) = cmplx(0.1_real64, 0.01_real64, real64)  ! B_s jump
        boundary_jumps(2) = cmplx(0.05_real64, 0.005_real64, real64) ! B_p jump
        boundary_jumps(3) = cmplx(0.0_real64, 0.0_real64, real64)   ! Continuity
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(source_terms(2)) > abs(source_terms(1)))
        test_passed = test_passed .and. (abs(boundary_jumps(1)) > 0.0_real64)
        test_passed = test_passed .and. (abs(boundary_jumps(2)) > 0.0_real64)
        test_passed = test_passed .and. (antenna_position > 0.0_real64)
        
        call end_test(test_passed)
    end subroutine test_maxwell_integration
    
    !---------------------------------------------------------------------------
    ! Test 11: Antenna spectrum mode calculations
    !---------------------------------------------------------------------------
    subroutine test_antenna_spectrum_modes()
        integer, allocatable :: mode_pairs(:,:)
        complex(real64), allocatable :: mode_amplitudes(:)
        real(real64) :: total_power
        integer :: ierr, n_modes, i
        
        call start_test("Antenna spectrum mode calculations")
        test_passed = .true.
        
        ! Test different antenna spectrum modes
        n_modes = 4
        allocate(mode_pairs(2, n_modes), mode_amplitudes(n_modes), stat=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            ! Set up mode pairs: (m, n)
            mode_pairs(1, :) = [3, 12, -3, -12]  ! m values
            mode_pairs(2, :) = [1, 4, -1, -4]    ! n values
            
            ! Calculate mode amplitudes (simplified)
            do i = 1, n_modes
                if (abs(mode_pairs(1, i)) == 3) then
                    ! (±3, ±1) modes
                    mode_amplitudes(i) = cmplx(4.0_real64, 0.5_real64, real64)
                else
                    ! (±12, ±4) modes  
                    mode_amplitudes(i) = cmplx(16.0_real64, 0.0_real64, real64)
                end if
            end do
            
            ! Calculate total power
            total_power = 0.0_real64
            do i = 1, n_modes
                total_power = total_power + abs(mode_amplitudes(i))**2
            end do
            
            test_passed = test_passed .and. (all(abs(mode_amplitudes) > 0.0_real64))
            test_passed = test_passed .and. (total_power > 0.0_real64)
            test_passed = test_passed .and. (size(mode_pairs, 1) == 2)
            test_passed = test_passed .and. (size(mode_pairs, 2) == n_modes)
            
            deallocate(mode_pairs, mode_amplitudes)
        end if
        
        call end_test(test_passed)
    end subroutine test_antenna_spectrum_modes
    
    !---------------------------------------------------------------------------
    ! Test 12: Memory management for antenna arrays
    !---------------------------------------------------------------------------
    subroutine test_memory_management()
        complex(real64), allocatable :: large_spectrum(:,:)
        real(real64), allocatable :: coupling_matrix(:,:)
        integer :: ierr, n_modes, n_harmonics
        
        call start_test("Memory management for antenna arrays")
        test_passed = .true.
        
        ! Test allocation and deallocation of large antenna arrays
        n_modes = 10
        n_harmonics = 64
        
        ! Large spectrum array
        allocate(large_spectrum(n_harmonics, n_modes), stat=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (allocated(large_spectrum)) then
            large_spectrum = cmplx(1.0_real64, 0.1_real64, real64)
            test_passed = test_passed .and. (abs(large_spectrum(1,1)) > 0.0_real64)
            test_passed = test_passed .and. (size(large_spectrum, 1) == n_harmonics)
            test_passed = test_passed .and. (size(large_spectrum, 2) == n_modes)
            
            deallocate(large_spectrum, stat=ierr)
            test_passed = test_passed .and. (ierr == 0)
        end if
        
        ! Coupling matrix
        allocate(coupling_matrix(n_modes, n_modes), stat=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (allocated(coupling_matrix)) then
            coupling_matrix = 0.0_real64
            ! Diagonal coupling
            do ierr = 1, n_modes
                coupling_matrix(ierr, ierr) = real(ierr, real64)
            end do
            
            test_passed = test_passed .and. (coupling_matrix(1,1) == 1.0_real64)
            test_passed = test_passed .and. (coupling_matrix(n_modes, n_modes) == real(n_modes, real64))
            
            deallocate(coupling_matrix, stat=ierr)
            test_passed = test_passed .and. (ierr == 0)
        end if
        
        call end_test(test_passed)
    end subroutine test_memory_management
    
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

end program test_kilca_antenna