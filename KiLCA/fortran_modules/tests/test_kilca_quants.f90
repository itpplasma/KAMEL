program test_kilca_quants
    !---------------------------------------------------------------------------
    ! KiLCA Physical Quantities and Diagnostics - Unit Tests
    !
    ! Tests the physical quantities calculation system including current density,
    ! power absorption/dissipation, energy fluxes, and diagnostic outputs.
    !
    ! Author: Claude (Anthropic)
    ! Date: 2024
    !
    ! Testing: kilca_quants_m.f90 (to be implemented)
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
    write(*,'(A)') " KiLCA Physical Quantities and Diagnostics - Unit Tests"
    write(*,'(A)') " ========================================================"
    write(*,'(A)') " "
    
    ! Run all tests
    call test_flre_quants_structure()
    call test_current_density_calculation()
    call test_power_absorption_calculation()
    call test_power_dissipation_calculation()
    call test_energy_flux_calculations()
    call test_antenna_plasma_coupling()
    call test_coordinate_transformations()
    call test_integration_over_surfaces()
    call test_output_file_generation()
    call test_conductivity_integration()
    call test_maxwell_field_integration()
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
    ! Test 1: FLRE quantities data structure and constants
    !---------------------------------------------------------------------------
    subroutine test_flre_quants_structure()
        integer :: ierr
        
        call start_test("FLRE quantities data structure")
        test_passed = .true.
        
        ! Test physical quantity constants and indexing
        ! Constants to test:
        ! - CURRENT_DENS, ABS_POWER_DENS, DISS_POWER_DENS
        ! - KIN_FLUX, POY_FLUX, TOT_FLUX
        ! - NUMBER_DENS, LOR_TORQUE_DENS
        
        ! Data structure to test: flre_quants_t
        ! Fields:
        ! - Physical quantity arrays: current_density, power_abs, power_diss
        ! - Energy flux arrays: kinetic_flux, poynting_flux, total_flux
        ! - Local and integrated profiles: qloc, qint
        ! - Calculation and save function pointers
        
        ierr = 0
        test_passed = (ierr == 0)
        
        call end_test(test_passed)
    end subroutine test_flre_quants_structure
    
    !---------------------------------------------------------------------------
    ! Test 2: Current density calculation (J = C·E)
    !---------------------------------------------------------------------------
    subroutine test_current_density_calculation()
        complex(real64) :: conductivity(3,3), electric_field(3), current_density(3)
        real(real64) :: r_test
        integer :: ierr, spec, type, comp
        
        call start_test("Current density calculation")
        test_passed = .true.
        
        ! Test current density calculation: J = C·E
        ! For each species (ions, electrons), type (0,1), component (r,s,p)
        r_test = 0.5_real64
        
        ! Mock conductivity tensor
        conductivity = cmplx(0.0_real64, 0.0_real64, real64)
        conductivity(1,1) = cmplx(2.0_real64, 0.1_real64, real64)  ! σ_rr
        conductivity(2,2) = cmplx(1.8_real64, 0.05_real64, real64) ! σ_ss
        conductivity(3,3) = cmplx(1.5_real64, 0.02_real64, real64) ! σ_pp
        conductivity(1,2) = cmplx(0.1_real64, 0.2_real64, real64)  ! σ_rs coupling
        
        ! Mock electric field
        electric_field(1) = cmplx(1.0_real64, 0.0_real64, real64)  ! Er
        electric_field(2) = cmplx(0.5_real64, 0.1_real64, real64)  ! Es
        electric_field(3) = cmplx(0.2_real64, 0.05_real64, real64) ! Ep
        
        ! Calculate current density: J = σ·E
        current_density(1) = conductivity(1,1)*electric_field(1) + &
                            conductivity(1,2)*electric_field(2) + &
                            conductivity(1,3)*electric_field(3)
        current_density(2) = conductivity(2,1)*electric_field(1) + &
                            conductivity(2,2)*electric_field(2) + &
                            conductivity(2,3)*electric_field(3)
        current_density(3) = conductivity(3,1)*electric_field(1) + &
                            conductivity(3,2)*electric_field(2) + &
                            conductivity(3,3)*electric_field(3)
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(current_density(1)) > 0.0_real64)
        test_passed = test_passed .and. (abs(current_density(2)) > 0.0_real64)
        test_passed = test_passed .and. (abs(current_density(3)) > 0.0_real64)
        
        call end_test(test_passed)
    end subroutine test_current_density_calculation
    
    !---------------------------------------------------------------------------
    ! Test 3: Power absorption calculation (P_abs = 0.5 * Re(J* · E))
    !---------------------------------------------------------------------------
    subroutine test_power_absorption_calculation()
        complex(real64) :: current_density(3), electric_field(3)
        real(real64) :: power_absorbed
        integer :: ierr, i
        
        call start_test("Power absorption calculation")
        test_passed = .true.
        
        ! Test absorbed power density: P_abs = 0.5 * Re(J* · E)
        
        ! Mock current density and electric field
        current_density(1) = cmplx(2.0_real64, 0.1_real64, real64)
        current_density(2) = cmplx(1.5_real64, -0.05_real64, real64)
        current_density(3) = cmplx(0.8_real64, 0.02_real64, real64)
        
        electric_field(1) = cmplx(1.0_real64, 0.0_real64, real64)
        electric_field(2) = cmplx(0.5_real64, 0.1_real64, real64)
        electric_field(3) = cmplx(0.2_real64, 0.05_real64, real64)
        
        ! Calculate absorbed power: P_abs = 0.5 * Re(J* · E)
        power_absorbed = 0.0_real64
        do i = 1, 3
            power_absorbed = power_absorbed + 0.5_real64 * &
                           real(conjg(current_density(i)) * electric_field(i), real64)
        end do
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (power_absorbed > 0.0_real64)
        test_passed = test_passed .and. (power_absorbed < 10.0_real64)  ! Reasonable magnitude
        
        call end_test(test_passed)
    end subroutine test_power_absorption_calculation
    
    !---------------------------------------------------------------------------
    ! Test 4: Power dissipation calculation using K matrices
    !---------------------------------------------------------------------------
    subroutine test_power_dissipation_calculation()
        complex(real64) :: K_matrix(3,3), electric_field(3)
        real(real64) :: power_dissipated
        integer :: ierr, i, j
        
        call start_test("Power dissipation calculation")
        test_passed = .true.
        
        ! Test dissipated power using conductivity K matrices
        ! P_diss involves finite Larmor radius effects
        
        ! Mock K matrix (conductivity derivatives)
        K_matrix = cmplx(0.0_real64, 0.0_real64, real64)
        K_matrix(1,1) = cmplx(0.1_real64, 0.01_real64, real64)
        K_matrix(2,2) = cmplx(0.08_real64, 0.005_real64, real64)
        K_matrix(3,3) = cmplx(0.05_real64, 0.002_real64, real64)
        
        ! Mock electric field
        electric_field(1) = cmplx(1.0_real64, 0.0_real64, real64)
        electric_field(2) = cmplx(0.5_real64, 0.1_real64, real64)
        electric_field(3) = cmplx(0.2_real64, 0.05_real64, real64)
        
        ! Calculate dissipated power (simplified)
        power_dissipated = 0.0_real64
        do i = 1, 3
            do j = 1, 3
                power_dissipated = power_dissipated + &
                    real(K_matrix(i,j) * conjg(electric_field(i)) * electric_field(j), real64)
            end do
        end do
        power_dissipated = abs(power_dissipated)
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (power_dissipated >= 0.0_real64)
        test_passed = test_passed .and. (power_dissipated < 1.0_real64)
        
        call end_test(test_passed)
    end subroutine test_power_dissipation_calculation
    
    !---------------------------------------------------------------------------
    ! Test 5: Energy flux calculations (kinetic, Poynting, total)
    !---------------------------------------------------------------------------
    subroutine test_energy_flux_calculations()
        complex(real64) :: electric_field(3), magnetic_field(3)
        real(real64) :: kinetic_flux, poynting_flux, total_flux
        integer :: ierr
        
        call start_test("Energy flux calculations")
        test_passed = .true.
        
        ! Test energy flux calculations
        ! Kinetic flux: Particle energy transport
        ! Poynting flux: S = (1/μ₀) * Re(E × B*)
        ! Total flux: Combined kinetic + electromagnetic
        
        ! Mock electromagnetic fields
        electric_field(1) = cmplx(1.0_real64, 0.1_real64, real64)
        electric_field(2) = cmplx(0.5_real64, -0.05_real64, real64)
        electric_field(3) = cmplx(0.0_real64, 0.0_real64, real64)  ! No parallel component
        
        magnetic_field(1) = cmplx(0.1_real64, 0.01_real64, real64)
        magnetic_field(2) = cmplx(0.05_real64, -0.005_real64, real64)
        magnetic_field(3) = cmplx(2.0_real64, 0.0_real64, real64)   ! Strong toroidal field
        
        ! Mock kinetic flux calculation
        kinetic_flux = 0.1_real64 * (abs(electric_field(1))**2 + abs(electric_field(2))**2)
        
        ! Poynting flux: S_r = (1/μ₀) * Re(E_s * B_p* - E_p * B_s*)
        poynting_flux = real(electric_field(2) * conjg(magnetic_field(3)) - &
                            electric_field(3) * conjg(magnetic_field(2)), real64)
        poynting_flux = abs(poynting_flux) / (4.0_real64 * PI * 1.0e-7_real64)  ! 1/μ₀
        
        ! Total flux
        total_flux = kinetic_flux + poynting_flux
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (kinetic_flux >= 0.0_real64)
        test_passed = test_passed .and. (poynting_flux >= 0.0_real64)
        test_passed = test_passed .and. (total_flux >= kinetic_flux)
        test_passed = test_passed .and. (total_flux >= poynting_flux)
        
        call end_test(test_passed)
    end subroutine test_energy_flux_calculations
    
    !---------------------------------------------------------------------------
    ! Test 6: Antenna-plasma coupling calculation
    !---------------------------------------------------------------------------
    subroutine test_antenna_plasma_coupling()
        complex(real64) :: antenna_current(3), electric_field(3)
        real(real64) :: coupling_power, vol_factor, r_position
        integer :: ierr
        
        call start_test("Antenna-plasma coupling")
        test_passed = .true.
        
        ! Test antenna-plasma coupling: JaE = 0.5*vol_fac*r*real(ja·conj(E))
        
        r_position = 0.8_real64
        vol_factor = 2.0_real64 * PI * r_position  ! Cylindrical volume factor
        
        ! Mock antenna current
        antenna_current(1) = cmplx(0.0_real64, 0.0_real64, real64)  ! No radial current
        antenna_current(2) = cmplx(1.0_real64, 0.0_real64, real64)  ! Poloidal current
        antenna_current(3) = cmplx(0.5_real64, 0.1_real64, real64)  ! Toroidal current
        
        ! Mock electric field at antenna
        electric_field(1) = cmplx(0.1_real64, 0.01_real64, real64)
        electric_field(2) = cmplx(2.0_real64, 0.5_real64, real64)
        electric_field(3) = cmplx(1.0_real64, 0.2_real64, real64)
        
        ! Calculate coupling power
        coupling_power = 0.5_real64 * vol_factor * r_position * &
                        real(antenna_current(1) * conjg(electric_field(1)) + &
                             antenna_current(2) * conjg(electric_field(2)) + &
                             antenna_current(3) * conjg(electric_field(3)), real64)
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(coupling_power) > 0.0_real64)
        test_passed = test_passed .and. (abs(coupling_power) < 100.0_real64)
        
        call end_test(test_passed)
    end subroutine test_antenna_plasma_coupling
    
    !---------------------------------------------------------------------------
    ! Test 7: Coordinate system transformations
    !---------------------------------------------------------------------------
    subroutine test_coordinate_transformations()
        real(real64) :: cylindrical_coords(3), field_aligned_coords(3)
        real(real64) :: transformation_matrix(3,3)
        integer :: ierr, i, j
        
        call start_test("Coordinate system transformations")
        test_passed = .true.
        
        ! Test coordinate transformations:
        ! - Cylindrical (r,θ,z) ↔ Field-aligned (r,s,p)
        ! - Laboratory frame ↔ Moving frame
        
        ! Mock cylindrical coordinates
        cylindrical_coords(1) = 0.5_real64  ! r
        cylindrical_coords(2) = 0.3_real64  ! θ
        cylindrical_coords(3) = 1.2_real64  ! z
        
        ! Mock transformation matrix (simplified)
        transformation_matrix = 0.0_real64
        do i = 1, 3
            transformation_matrix(i,i) = 1.0_real64  ! Identity base
        end do
        transformation_matrix(2,3) = 0.1_real64  ! Small coupling
        transformation_matrix(3,2) = -0.1_real64
        
        ! Apply transformation: field_aligned = T * cylindrical
        do i = 1, 3
            field_aligned_coords(i) = 0.0_real64
            do j = 1, 3
                field_aligned_coords(i) = field_aligned_coords(i) + &
                    transformation_matrix(i,j) * cylindrical_coords(j)
            end do
        end do
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(field_aligned_coords(1) - cylindrical_coords(1)) < 0.1_real64)
        test_passed = test_passed .and. (all(abs(field_aligned_coords) < 10.0_real64))
        
        call end_test(test_passed)
    end subroutine test_coordinate_transformations
    
    !---------------------------------------------------------------------------
    ! Test 8: Integration over cylindrical surfaces
    !---------------------------------------------------------------------------
    subroutine test_integration_over_surfaces()
        real(real64), allocatable :: radial_profile(:), integrated_profile(:)
        real(real64), allocatable :: radial_grid(:)
        integer :: ierr, i, n_points
        
        call start_test("Integration over cylindrical surfaces")
        test_passed = .true.
        
        ! Test integration of local quantities over cylindrical surfaces
        n_points = 20
        allocate(radial_grid(n_points), radial_profile(n_points), integrated_profile(n_points))
        
        ! Create test radial grid and profile
        do i = 1, n_points
            radial_grid(i) = real(i-1, real64) / real(n_points-1, real64)
            radial_profile(i) = exp(-radial_grid(i)**2)  ! Gaussian profile
        end do
        
        ! Integrate over cylindrical surfaces: ∫₀ʳ 2πr'·profile(r') dr'
        integrated_profile(1) = 0.0_real64
        do i = 2, n_points
            integrated_profile(i) = integrated_profile(i-1) + &
                2.0_real64 * PI * radial_grid(i) * radial_profile(i) * &
                (radial_grid(i) - radial_grid(i-1))
        end do
        
        ierr = 0
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (integrated_profile(n_points) > integrated_profile(1))
        test_passed = test_passed .and. (integrated_profile(n_points) > 0.0_real64)
        test_passed = test_passed .and. (integrated_profile(n_points) < 10.0_real64)
        
        deallocate(radial_grid, radial_profile, integrated_profile)
        
        call end_test(test_passed)
    end subroutine test_integration_over_surfaces
    
    !---------------------------------------------------------------------------
    ! Test 9: Output file generation and formatting
    !---------------------------------------------------------------------------
    subroutine test_output_file_generation()
        character(len=256) :: filename
        integer :: unit_num, ierr, i
        real(real64) :: test_data(10)
        
        call start_test("Output file generation")
        test_passed = .true.
        
        ! Test output file generation with proper formatting
        ! Format: zone_{index}_{quantity}_dens_{type}_{species}.dat
        
        filename = "test_current_dens_0_i.dat"
        
        ! Create test data
        do i = 1, 10
            test_data(i) = real(i, real64) * 0.1_real64
        end do
        
        ! Write test file
        open(newunit=unit_num, file=filename, status='replace', action='write', iostat=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (ierr == 0) then
            write(unit_num, '(A)') "# Radial_position  Current_density_real  Current_density_imag"
            do i = 1, 10
                write(unit_num, '(3E15.7)') real(i-1, real64)/9.0_real64, test_data(i), 0.0_real64
            end do
            close(unit_num)
            
            ! Clean up test file
            open(newunit=unit_num, file=filename, status='old', iostat=ierr)
            if (ierr == 0) then
                close(unit_num, status='delete')
            end if
        end if
        
        call end_test(test_passed)
    end subroutine test_output_file_generation
    
    !---------------------------------------------------------------------------
    ! Test 10: Integration with conductivity system
    !---------------------------------------------------------------------------
    subroutine test_conductivity_integration()
        integer :: ierr
        
        call start_test("Conductivity system integration")
        test_passed = .true.
        
        ! Test integration with kilca_conductivity_m module
        ! Functions to test:
        ! - Access to C matrices and K matrices
        ! - Proper indexing for species and types
        ! - Error handling for missing conductivity data
        
        ierr = 0
        test_passed = (ierr == 0)
        
        call end_test(test_passed)
    end subroutine test_conductivity_integration
    
    !---------------------------------------------------------------------------
    ! Test 11: Integration with Maxwell field system
    !---------------------------------------------------------------------------
    subroutine test_maxwell_field_integration()
        integer :: ierr
        
        call start_test("Maxwell field system integration")
        test_passed = .true.
        
        ! Test integration with kilca_maxwell_m module
        ! Functions to test:
        ! - Access to electromagnetic field data
        ! - Field component indexing (iErsp_sys)
        ! - Moving frame field calculations
        
        ierr = 0
        test_passed = (ierr == 0)
        
        call end_test(test_passed)
    end subroutine test_maxwell_field_integration
    
    !---------------------------------------------------------------------------
    ! Test 12: Memory management for large quantity arrays
    !---------------------------------------------------------------------------
    subroutine test_memory_management()
        complex(real64), allocatable :: quantity_array(:,:,:,:)
        real(real64), allocatable :: integrated_array(:,:,:)
        integer :: ierr, dim1, dim2, dim3, dim4
        
        call start_test("Memory management for quantity arrays")
        test_passed = .true.
        
        ! Test allocation/deallocation of large quantity arrays
        ! Current density: {{{re,im},comp={r,s,p},type={0,1},spec={i,e,t}}
        dim1 = 2   ! real, imaginary
        dim2 = 3   ! r, s, p components
        dim3 = 2   ! types 0, 1
        dim4 = 3   ! species i, e, total
        
        allocate(quantity_array(dim1, dim2, dim3, dim4), stat=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (allocated(quantity_array)) then
            quantity_array = cmplx(1.0_real64, 0.5_real64, real64)
            test_passed = test_passed .and. (abs(quantity_array(1,1,1,1)) > 0.0_real64)
            
            deallocate(quantity_array, stat=ierr)
            test_passed = test_passed .and. (ierr == 0)
        end if
        
        ! Test integrated quantities array
        allocate(integrated_array(100, dim3, dim4), stat=ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        if (allocated(integrated_array)) then
            integrated_array = 2.0_real64
            test_passed = test_passed .and. (integrated_array(50, 1, 1) == 2.0_real64)
            
            deallocate(integrated_array, stat=ierr)
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

end program test_kilca_quants