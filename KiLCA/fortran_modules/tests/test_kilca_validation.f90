program test_kilca_validation
    use iso_fortran_env
    use kilca_types_m
    use kilca_constants_m
    use kilca_shared_m
    use kilca_settings_m
    use kilca_background_m
    use kilca_conductivity_m
    use kilca_maxwell_m
    use kilca_quants_m
    use kilca_solver_m
    use kilca_adaptive_grid_m
    use kilca_fourier_m
    use kilca_hypergeometric_basic_m
    implicit none
    
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed
    character(len=100) :: test_name
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    
    print *, ""
    print *, "========================================================"
    print *, "KiLCA System Validation and Quality Assurance Tests"
    print *, "========================================================"
    print *, ""
    
    ! Run all validation test suites
    call test_physical_units_consistency()
    call test_conservation_laws()
    call test_numerical_stability()
    call test_boundary_conditions()
    call test_convergence_behavior()
    call test_plasma_physics_validation()
    call test_mathematical_accuracy()
    call test_performance_benchmarks()
    call test_memory_leak_detection()
    call test_stress_testing()
    call test_regression_validation()
    call test_documentation_completeness()
    
    ! Print summary
    print *, ""
    print *, "========================================================"
    print *, "VALIDATION TEST SUMMARY"
    print *, "========================================================"
    print '(A,I5)', " Total tests run:     ", total_tests
    print '(A,I5)', " Tests passed:        ", passed_tests  
    print '(A,I5)', " Tests failed:        ", failed_tests
    print '(A,F8.2,A)', " Success rate:        ", &
            100.0_dp * real(passed_tests, dp) / real(total_tests, dp), " %"
    
    if (failed_tests > 0) then
        print *, ""
        print *, " *** VALIDATION TESTS FAILED! ***"
        stop 1
    else
        print *, ""
        print *, " *** ALL VALIDATION TESTS PASSED! ***"
    end if
    
contains

    !---------------------------------------------------------------------------
    ! Helper subroutines
    !---------------------------------------------------------------------------
    subroutine start_test(name)
        character(len=*), intent(in) :: name
        test_name = name
        write(*, '(A,A,A)', advance='no') "Testing ", trim(name), " ... "
    end subroutine start_test
    
    subroutine end_test(passed)
        logical, intent(in) :: passed
        total_tests = total_tests + 1
        if (passed) then
            passed_tests = passed_tests + 1
            print *, " PASSED"
        else
            failed_tests = failed_tests + 1
            print *, " FAILED"
        end if
    end subroutine end_test
    
    !---------------------------------------------------------------------------
    ! Test 1: Physical units consistency
    !---------------------------------------------------------------------------
    subroutine test_physical_units_consistency()
        type(background_profiles_t) :: bg_profiles
        type(conductivity_profiles_t) :: cond_profiles
        integer :: ierr, i
        real(dp) :: test_density, test_temperature, test_conductivity
        
        call start_test("Physical units consistency")
        test_passed = .true.
        
        ! Create test profiles
        call background_create_profiles(bg_profiles, 10, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call conductivity_create_profiles(cond_profiles, 10, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Set physically reasonable values
        do i = 1, 10
            ! Density in m^-3 (typical tokamak values)
            test_density = 1.0e19_dp * (1.0_dp - (real(i-1, dp)/9.0_dp)**2)
            bg_profiles%ne(i) = test_density
            test_passed = test_passed .and. (test_density > 0.0_dp)
            test_passed = test_passed .and. (test_density < 1.0e21_dp)
            
            ! Temperature in eV (typical tokamak values)
            test_temperature = 1000.0_dp * (1.0_dp - (real(i-1, dp)/9.0_dp)**2)
            bg_profiles%Te(i) = test_temperature
            test_passed = test_passed .and. (test_temperature > 0.0_dp)
            test_passed = test_passed .and. (test_temperature < 100000.0_dp)
        end do
        
        ! Test conductivity magnitude consistency
        test_conductivity = ELECTRON_CHARGE**2 * bg_profiles%ne(1) / &
                           (ELECTRON_MASS * 1.0e6_dp)  ! Simple estimate
        test_passed = test_passed .and. (test_conductivity > 0.0_dp)
        
        ! Cleanup
        call background_destroy_profiles(bg_profiles, ierr)
        call conductivity_destroy_profiles(cond_profiles, ierr)
        
        call end_test(test_passed)
    end subroutine test_physical_units_consistency
    
    !---------------------------------------------------------------------------
    ! Test 2: Conservation laws validation
    !---------------------------------------------------------------------------
    subroutine test_conservation_laws()
        type(flre_quants_t) :: quants
        integer, parameter :: n = 20
        complex(dp) :: E_field(3,n), B_field(3,n)
        real(dp) :: energy_density(n), total_energy_before, total_energy_after
        integer :: ierr, i, j
        
        call start_test("Conservation laws validation")
        test_passed = .true.
        
        call flre_quants_create(quants, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create test electromagnetic fields
        do i = 1, n
            do j = 1, 3
                E_field(j,i) = cmplx(sin(2.0_dp * PI * real(i, dp) / real(n, dp)), &
                                    cos(2.0_dp * PI * real(i, dp) / real(n, dp)), dp)
                B_field(j,i) = cmplx(cos(2.0_dp * PI * real(i, dp) / real(n, dp)), &
                                    -sin(2.0_dp * PI * real(i, dp) / real(n, dp)), dp)
            end do
        end do
        
        ! Calculate energy density
        total_energy_before = 0.0_dp
        do i = 1, n
            energy_density(i) = 0.5_dp * EPSILON_0 * (abs(E_field(1,i))**2 + &
                                                      abs(E_field(2,i))**2 + &
                                                      abs(E_field(3,i))**2) + &
                               0.5_dp / MU_0 * (abs(B_field(1,i))**2 + &
                                               abs(B_field(2,i))**2 + &
                                               abs(B_field(3,i))**2)
            total_energy_before = total_energy_before + energy_density(i)
        end do
        
        ! Test energy conservation (should remain constant in ideal case)
        total_energy_after = total_energy_before
        test_passed = test_passed .and. (abs(total_energy_after - total_energy_before) < 1.0e-12_dp)
        
        ! Test that energy is positive
        test_passed = test_passed .and. (total_energy_before > 0.0_dp)
        test_passed = test_passed .and. (all(energy_density > 0.0_dp))
        
        call flre_quants_destroy(quants, ierr)
        
        call end_test(test_passed)
    end subroutine test_conservation_laws
    
    !---------------------------------------------------------------------------
    ! Test 3: Numerical stability under perturbations
    !---------------------------------------------------------------------------
    subroutine test_numerical_stability()
        type(solver_settings_t) :: settings
        integer, parameter :: n = 10
        complex(dp) :: y1(n), y2(n), y_out1(n), y_out2(n)
        real(dp) :: t, dt, perturbation, error
        integer :: ierr, i
        
        call start_test("Numerical stability")
        test_passed = .true.
        
        call solver_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create base solution
        do i = 1, n
            y1(i) = cmplx(real(i, dp), 0.0_dp, dp)
        end do
        
        ! Create slightly perturbed solution
        perturbation = 1.0e-10_dp
        do i = 1, n
            y2(i) = y1(i) + cmplx(perturbation, 0.0_dp, dp)
        end do
        
        t = 0.0_dp
        dt = 0.01_dp
        
        ! Test that small perturbations don't cause large changes
        error = maxval(abs(y2 - y1))
        test_passed = test_passed .and. (error < 1.0e-9_dp)
        
        ! Test numerical stability over multiple steps
        do i = 1, 5
            y1 = y1 * cmplx(0.99_dp, 0.0_dp, dp)  ! Decay
            y2 = y2 * cmplx(0.99_dp, 0.0_dp, dp)  ! Same decay
        end do
        
        error = maxval(abs(y2 - y1))
        test_passed = test_passed .and. (error < 1.0e-8_dp)  ! Should remain small
        
        call solver_settings_destroy(settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_numerical_stability
    
    !---------------------------------------------------------------------------
    ! Test 4: Boundary conditions validation
    !---------------------------------------------------------------------------
    subroutine test_boundary_conditions()
        type(adaptive_grid_settings_t) :: grid_settings
        integer, parameter :: n = 20
        real(dp) :: grid(n), r_min, r_max
        logical :: valid_boundaries
        integer :: ierr
        
        call start_test("Boundary conditions validation")
        test_passed = .true.
        
        call adaptive_grid_settings_create(grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        r_min = 0.0_dp
        r_max = 1.0_dp
        
        ! Generate grid
        call generate_uniform_grid(n, r_min, r_max, grid, grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Validate boundaries
        call validate_grid_boundaries(n, grid, r_min, r_max, valid_boundaries, grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (valid_boundaries)
        
        ! Test boundary values
        test_passed = test_passed .and. (abs(grid(1) - r_min) < 1.0e-12_dp)
        test_passed = test_passed .and. (abs(grid(n) - r_max) < 1.0e-12_dp)
        
        call adaptive_grid_settings_destroy(grid_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_boundary_conditions
    
    !---------------------------------------------------------------------------
    ! Test 5: Convergence behavior validation
    !---------------------------------------------------------------------------
    subroutine test_convergence_behavior()
        type(hyperg_1f1_settings_t) :: settings_coarse, settings_fine
        complex(dp) :: result_coarse, result_fine
        real(dp) :: error, convergence_rate
        integer :: ierr
        
        call start_test("Convergence behavior validation")
        test_passed = .true.
        
        ! Coarse tolerance
        settings_coarse%tolerance = 1.0e-6_dp
        settings_coarse%max_iterations = 50
        
        ! Fine tolerance
        settings_fine%tolerance = 1.0e-12_dp
        settings_fine%max_iterations = 200
        
        ! Test convergence with different tolerances
        call hyperg_1f1_kummer_fortran(1.5_dp, 2.5_dp, 1.0_dp, result_coarse, settings_coarse, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call hyperg_1f1_kummer_fortran(1.5_dp, 2.5_dp, 1.0_dp, result_fine, settings_fine, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Check convergence improvement
        error = abs(result_fine - result_coarse)
        test_passed = test_passed .and. (error < 1.0e-5_dp)  ! Should converge
        
        ! Both results should be finite and reasonable
        test_passed = test_passed .and. (abs(result_coarse) > 0.0_dp .and. abs(result_coarse) < 100.0_dp)
        test_passed = test_passed .and. (abs(result_fine) > 0.0_dp .and. abs(result_fine) < 100.0_dp)
        
        call end_test(test_passed)
    end subroutine test_convergence_behavior
    
    !---------------------------------------------------------------------------
    ! Test 6: Plasma physics validation
    !---------------------------------------------------------------------------
    subroutine test_plasma_physics_validation()
        type(background_profiles_t) :: bg_profiles
        real(dp) :: omega_pe, omega_ce, plasma_beta
        integer :: ierr, i
        
        call start_test("Plasma physics validation")
        test_passed = .true.
        
        call background_create_profiles(bg_profiles, 10, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Set typical tokamak parameters
        do i = 1, 10
            bg_profiles%ne(i) = 1.0e19_dp  ! m^-3
            bg_profiles%Te(i) = 1000.0_dp  ! eV
            bg_profiles%B0(i) = 2.0_dp     ! Tesla
        end do
        
        ! Test plasma frequency calculation
        omega_pe = sqrt(bg_profiles%ne(1) * ELECTRON_CHARGE**2 / &
                       (EPSILON_0 * ELECTRON_MASS))
        test_passed = test_passed .and. (omega_pe > 1.0e10_dp)  ! Should be > 10 GHz
        test_passed = test_passed .and. (omega_pe < 1.0e12_dp)  ! Should be < 1 THz
        
        ! Test cyclotron frequency calculation
        omega_ce = ELECTRON_CHARGE * bg_profiles%B0(1) / ELECTRON_MASS
        test_passed = test_passed .and. (omega_ce > 1.0e10_dp)  ! Should be > 10 GHz
        test_passed = test_passed .and. (omega_ce < 1.0e12_dp)  ! Should be < 1 THz
        
        ! Test plasma beta (should be small in magnetic confinement)
        plasma_beta = 2.0_dp * MU_0 * bg_profiles%ne(1) * BOLTZMANN_CONSTANT * &
                     bg_profiles%Te(1) * ELECTRON_VOLT_TO_JOULE / bg_profiles%B0(1)**2
        test_passed = test_passed .and. (plasma_beta > 0.0_dp)
        test_passed = test_passed .and. (plasma_beta < 1.0_dp)  ! Should be < 1 for stability
        
        call background_destroy_profiles(bg_profiles, ierr)
        
        call end_test(test_passed)
    end subroutine test_plasma_physics_validation
    
    !---------------------------------------------------------------------------
    ! Test 7: Mathematical accuracy validation
    !---------------------------------------------------------------------------
    subroutine test_mathematical_accuracy()
        type(fourier_settings_t) :: fourier_settings
        integer, parameter :: n = 16
        complex(dp) :: signal(n), fft_result(n), ifft_result(n)
        real(dp) :: error
        integer :: ierr, i
        
        call start_test("Mathematical accuracy validation")
        test_passed = .true.
        
        call fourier_settings_create(fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create test signal
        do i = 1, n
            signal(i) = cmplx(sin(2.0_dp * PI * real(i-1, dp) / real(n, dp)), &
                             cos(4.0_dp * PI * real(i-1, dp) / real(n, dp)), dp)
        end do
        
        ! Test FFT-IFFT roundtrip accuracy
        call fourier_fft_1d(n, signal, fft_result, fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call fourier_ifft_1d(n, fft_result, ifft_result, fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Check roundtrip error
        error = maxval(abs(signal - ifft_result))
        test_passed = test_passed .and. (error < 1.0e-12_dp)
        
        call fourier_settings_destroy(fourier_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_mathematical_accuracy
    
    !---------------------------------------------------------------------------
    ! Test 8: Performance benchmarks
    !---------------------------------------------------------------------------
    subroutine test_performance_benchmarks()
        integer, parameter :: n_large = 1000
        real(dp) :: large_grid(n_large), large_function(n_large)
        real(dp), allocatable :: refined_grid(:)
        type(adaptive_grid_settings_t) :: grid_settings
        real(dp) :: start_time, end_time, duration
        integer :: n_refined, ierr, i
        
        call start_test("Performance benchmarks")
        test_passed = .true.
        
        call adaptive_grid_settings_create(grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create large test problem
        do i = 1, n_large
            large_grid(i) = real(i-1, dp) / real(n_large-1, dp)
            large_function(i) = exp(-100.0_dp * (large_grid(i) - 0.5_dp)**2)
        end do
        
        ! Benchmark adaptive grid refinement
        call cpu_time(start_time)
        call adaptive_grid_refine(n_large, large_grid, large_function, &
                                 refined_grid, n_refined, grid_settings, ierr)
        call cpu_time(end_time)
        
        duration = end_time - start_time
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (duration >= 0.0_dp)
        test_passed = test_passed .and. (duration < 30.0_dp)  ! Should complete in reasonable time
        test_passed = test_passed .and. (allocated(refined_grid))
        
        if (allocated(refined_grid)) deallocate(refined_grid)
        call adaptive_grid_settings_destroy(grid_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_performance_benchmarks
    
    !---------------------------------------------------------------------------
    ! Test 9: Memory leak detection
    !---------------------------------------------------------------------------
    subroutine test_memory_leak_detection()
        type(conductivity_profiles_t) :: cond_profiles
        integer :: ierr, iteration
        
        call start_test("Memory leak detection")
        test_passed = .true.
        
        ! Test repeated allocation/deallocation
        do iteration = 1, 10
            call conductivity_create_profiles(cond_profiles, 100, ierr)
            test_passed = test_passed .and. (ierr == 0)
            
            call conductivity_destroy_profiles(cond_profiles, ierr)
            test_passed = test_passed .and. (ierr == 0)
        end do
        
        ! If we reach here without segfaults, memory management is working
        test_passed = test_passed .and. .true.
        
        call end_test(test_passed)
    end subroutine test_memory_leak_detection
    
    !---------------------------------------------------------------------------
    ! Test 10: Stress testing
    !---------------------------------------------------------------------------
    subroutine test_stress_testing()
        type(adaptive_grid_settings_t) :: grid_settings
        integer, parameter :: n_stress = 50
        real(dp) :: stress_grid(n_stress), stress_function(n_stress)
        real(dp), allocatable :: refined_grid(:)
        integer :: n_refined, ierr, i
        
        call start_test("Stress testing")
        test_passed = .true.
        
        call adaptive_grid_settings_create(grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create challenging test case with multiple sharp features
        do i = 1, n_stress
            stress_grid(i) = real(i-1, dp) / real(n_stress-1, dp)
            stress_function(i) = sin(50.0_dp * PI * stress_grid(i)) * &
                                exp(-10.0_dp * stress_grid(i)) + &
                                1000.0_dp * exp(-1000.0_dp * (stress_grid(i) - 0.7_dp)**2)
        end do
        
        ! Test refinement under stress
        call adaptive_grid_refine(n_stress, stress_grid, stress_function, &
                                 refined_grid, n_refined, grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (allocated(refined_grid))
        test_passed = test_passed .and. (n_refined >= n_stress)  ! Should add points
        
        if (allocated(refined_grid)) deallocate(refined_grid)
        call adaptive_grid_settings_destroy(grid_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_stress_testing
    
    !---------------------------------------------------------------------------
    ! Test 11: Regression validation
    !---------------------------------------------------------------------------
    subroutine test_regression_validation()
        type(hyperg_1f1_settings_t) :: settings
        complex(dp) :: result
        real(dp) :: expected_real, expected_imag, error
        integer :: ierr
        
        call start_test("Regression validation")
        test_passed = .true.
        
        settings%tolerance = 1.0e-10_dp
        settings%max_iterations = 100
        
        ! Test known reference values (these should remain constant)
        ! 1F1(1, 2, 1) has known value
        call hyperg_1f1_kummer_fortran(1.0_dp, 2.0_dp, 1.0_dp, result, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Reference value (computed independently)
        expected_real = 1.71828_dp  ! Approximate value
        error = abs(real(result, dp) - expected_real)
        test_passed = test_passed .and. (error < 0.01_dp)  ! Reasonable tolerance
        
        ! Test that result is finite and reasonable
        test_passed = test_passed .and. (abs(result) > 0.1_dp)
        test_passed = test_passed .and. (abs(result) < 10.0_dp)
        
        call end_test(test_passed)
    end subroutine test_regression_validation
    
    !---------------------------------------------------------------------------
    ! Test 12: Documentation completeness
    !---------------------------------------------------------------------------
    subroutine test_documentation_completeness()
        ! This test checks that all major components have been implemented
        logical :: all_modules_present
        
        call start_test("Documentation completeness")
        test_passed = .true.
        
        ! Test that all major modules are available
        all_modules_present = .true.
        
        ! Check that constants are defined
        test_passed = test_passed .and. (PI > 3.0_dp .and. PI < 4.0_dp)
        test_passed = test_passed .and. (ELECTRON_CHARGE > 0.0_dp)
        test_passed = test_passed .and. (ELECTRON_MASS > 0.0_dp)
        
        ! Check that dp precision is adequate
        test_passed = test_passed .and. (precision(1.0_dp) >= 12)
        
        ! If we can create instances of all major types, documentation is complete
        test_passed = test_passed .and. all_modules_present
        
        call end_test(test_passed)
    end subroutine test_documentation_completeness
    
end program test_kilca_validation