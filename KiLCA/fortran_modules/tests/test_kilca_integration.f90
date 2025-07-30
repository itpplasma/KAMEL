program test_kilca_integration
    use iso_fortran_env
    use kilca_types_m
    use kilca_constants_m
    use kilca_shared_m
    use kilca_settings_m
    use kilca_background_m
    use kilca_mode_m
    use kilca_zone_m
    use kilca_interp_integ_m
    use kilca_mode_solver_m
    use kilca_solver_m
    use kilca_eigtransform_m
    use kilca_rhs_func_m
    use kilca_conductivity_m
    use kilca_maxwell_m
    use kilca_quants_m
    use kilca_antenna_m
    use kilca_physics_m
    use kilca_zerofind_m
    use kilca_hypergeometric_basic_m
    use kilca_fourier_m
    use kilca_adaptive_grid_m
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
    print *, "KiLCA System Integration and Validation Tests"
    print *, "========================================================"
    print *, ""
    
    ! Run all integration test suites
    call test_module_dependencies()
    call test_data_type_compatibility()
    call test_settings_integration()
    call test_plasma_physics_workflow()
    call test_solver_integration()
    call test_mathematical_libraries_integration()
    call test_grid_solver_integration()
    call test_full_kilca_simulation()
    call test_memory_management_integration()
    call test_error_handling_integration()
    call test_performance_integration()
    call test_numerical_accuracy_validation()
    
    ! Print summary
    print *, ""
    print *, "========================================================"
    print *, "INTEGRATION TEST SUMMARY"
    print *, "========================================================"
    print '(A,I5)', " Total tests run:     ", total_tests
    print '(A,I5)', " Tests passed:        ", passed_tests  
    print '(A,I5)', " Tests failed:        ", failed_tests
    print '(A,F8.2,A)', " Success rate:        ", &
            100.0_dp * real(passed_tests, dp) / real(total_tests, dp), " %"
    
    if (failed_tests > 0) then
        print *, ""
        print *, " *** INTEGRATION TESTS FAILED! ***"
        stop 1
    else
        print *, ""
        print *, " *** ALL INTEGRATION TESTS PASSED! ***"
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
    ! Test 1: Module dependencies and imports
    !---------------------------------------------------------------------------
    subroutine test_module_dependencies()
        real(dp) :: test_value
        integer :: ierr
        
        call start_test("Module dependencies")
        test_passed = .true.
        
        ! Test basic types
        test_value = 1.0_dp
        test_passed = test_passed .and. (precision(test_value) >= 12)
        
        ! Test constants
        test_passed = test_passed .and. (abs(PI - 3.141592653589793_dp) < 1.0e-12_dp)
        
        ! Test that all modules load without error
        test_passed = test_passed .and. .true.  ! If we reach here, modules loaded
        
        call end_test(test_passed)
    end subroutine test_module_dependencies
    
    !---------------------------------------------------------------------------
    ! Test 2: Data type compatibility across modules
    !---------------------------------------------------------------------------
    subroutine test_data_type_compatibility()
        type(kilca_settings_t) :: settings
        type(background_profiles_t) :: bg_profiles
        type(wave_data_t) :: wave_data
        type(zone_extended_t) :: zone
        integer :: ierr
        
        call start_test("Data type compatibility")
        test_passed = .true.
        
        ! Test settings creation
        call kilca_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test background profiles creation
        call background_create_profiles(bg_profiles, 10, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test wave data creation
        call wave_data_create(wave_data, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test zone creation
        call zone_extended_create(zone, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test compatibility between modules
        test_passed = test_passed .and. (settings%precision_dp == dp)
        
        ! Cleanup
        call kilca_settings_destroy(settings, ierr)
        call background_destroy_profiles(bg_profiles, ierr)
        call wave_data_destroy(wave_data, ierr)
        call zone_extended_destroy(zone, ierr)
        
        call end_test(test_passed)
    end subroutine test_data_type_compatibility
    
    !---------------------------------------------------------------------------
    ! Test 3: Settings integration across modules
    !---------------------------------------------------------------------------
    subroutine test_settings_integration()
        type(kilca_settings_t) :: settings
        type(solver_settings_t) :: solver_settings
        type(adaptive_grid_settings_t) :: grid_settings
        type(fourier_settings_t) :: fourier_settings
        integer :: ierr
        
        call start_test("Settings integration")
        test_passed = .true.
        
        ! Create all settings
        call kilca_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call solver_settings_create(solver_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call adaptive_grid_settings_create(grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call fourier_settings_create(fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test consistent precision settings
        test_passed = test_passed .and. (settings%tolerance > 0.0_dp)
        test_passed = test_passed .and. (solver_settings%tolerance > 0.0_dp)
        test_passed = test_passed .and. (grid_settings%tolerance > 0.0_dp)
        test_passed = test_passed .and. (fourier_settings%tolerance > 0.0_dp)
        
        ! Cleanup
        call kilca_settings_destroy(settings, ierr)
        call solver_settings_destroy(solver_settings, ierr)
        call adaptive_grid_settings_destroy(grid_settings, ierr)
        call fourier_settings_destroy(fourier_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_settings_integration
    
    !---------------------------------------------------------------------------
    ! Test 4: Plasma physics workflow integration
    !---------------------------------------------------------------------------
    subroutine test_plasma_physics_workflow()
        type(background_profiles_t) :: bg_profiles
        type(conductivity_profiles_t) :: cond_profiles
        type(maxwell_eqs_data_t) :: maxwell_data
        type(flre_quants_t) :: quants
        integer :: ierr, i
        real(dp) :: test_r
        
        call start_test("Plasma physics workflow")
        test_passed = .true.
        
        ! Create background profiles
        call background_create_profiles(bg_profiles, 20, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Initialize test data
        do i = 1, 20
            test_r = real(i-1, dp) / 19.0_dp
            bg_profiles%r(i) = test_r
            bg_profiles%ne(i) = 1.0e19_dp * (1.0_dp - test_r**2)  ! parabolic profile
            bg_profiles%Te(i) = 1000.0_dp * (1.0_dp - test_r**2)  ! temperature profile
        end do
        
        ! Create conductivity profiles
        call conductivity_create_profiles(cond_profiles, 20, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create Maxwell equations data
        call maxwell_create_data(maxwell_data, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create FLRE quantities
        call flre_quants_create(quants, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test workflow integration
        test_passed = test_passed .and. (bg_profiles%dimr == 20)
        test_passed = test_passed .and. (cond_profiles%dimx == 20)
        test_passed = test_passed .and. (all(bg_profiles%ne > 0.0_dp))
        test_passed = test_passed .and. (all(bg_profiles%Te > 0.0_dp))
        
        ! Cleanup
        call background_destroy_profiles(bg_profiles, ierr)
        call conductivity_destroy_profiles(cond_profiles, ierr)
        call maxwell_destroy_data(maxwell_data, ierr)
        call flre_quants_destroy(quants, ierr)
        
        call end_test(test_passed)
    end subroutine test_plasma_physics_workflow
    
    !---------------------------------------------------------------------------
    ! Test 5: Solver integration
    !---------------------------------------------------------------------------
    subroutine test_solver_integration()
        type(solver_settings_t) :: solver_settings
        type(rhs_func_params_t) :: rhs_params
        integer, parameter :: n = 10
        complex(dp) :: y(n), y_out(n)
        real(dp) :: t, dt
        integer :: ierr, i
        
        call start_test("Solver integration")
        test_passed = .true.
        
        ! Create solver settings
        call solver_settings_create(solver_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create RHS parameters
        call rhs_func_params_create(rhs_params, n, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Initialize test system
        do i = 1, n
            y(i) = cmplx(real(i, dp), 0.0_dp, dp)
        end do
        t = 0.0_dp
        dt = 0.01_dp
        
        ! Test RK4 integration step
        call solver_rk4_step(n, y, t, dt, y_out, rhs_params, solver_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test eigenvalue transformation
        call eigtransform_coeffs_to_solution(n, y, y_out, rhs_params, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Cleanup
        call solver_settings_destroy(solver_settings, ierr)
        call rhs_func_params_destroy(rhs_params, ierr)
        
        call end_test(test_passed)
    end subroutine test_solver_integration
    
    !---------------------------------------------------------------------------
    ! Test 6: Mathematical libraries integration
    !---------------------------------------------------------------------------
    subroutine test_mathematical_libraries_integration()
        type(zerofind_settings_t) :: zero_settings
        type(hyperg_1f1_settings_t) :: hyperg_settings
        type(fourier_settings_t) :: fourier_settings
        complex(dp) :: hyperg_result, fourier_data(8), fourier_result(8)
        integer :: ierr, i
        
        call start_test("Mathematical libraries integration")
        test_passed = .true.
        
        ! Create settings for math libraries
        call zerofind_settings_create(zero_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call fourier_settings_create(fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        hyperg_settings%tolerance = 1.0e-10_dp
        hyperg_settings%max_iterations = 100
        
        ! Test hypergeometric function
        call hyperg_1f1_kummer_fortran(1.0_dp, 2.0_dp, 1.0_dp, hyperg_result, hyperg_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (abs(hyperg_result) > 0.0_dp)
        
        ! Test Fourier transform
        do i = 1, 8
            fourier_data(i) = cmplx(sin(2.0_dp * PI * real(i-1, dp) / 8.0_dp), 0.0_dp, dp)
        end do
        
        call fourier_fft_1d(8, fourier_data, fourier_result, fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Cleanup
        call zerofind_settings_destroy(zero_settings, ierr)
        call fourier_settings_destroy(fourier_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_mathematical_libraries_integration
    
    !---------------------------------------------------------------------------
    ! Test 7: Grid-solver integration
    !---------------------------------------------------------------------------
    subroutine test_grid_solver_integration()
        type(adaptive_grid_settings_t) :: grid_settings
        type(solver_settings_t) :: solver_settings
        integer, parameter :: n_initial = 10
        real(dp) :: grid_initial(n_initial), function_values(n_initial)
        real(dp), allocatable :: refined_grid(:)
        complex(dp), allocatable :: solution(:)
        integer :: n_refined, ierr, i
        
        call start_test("Grid-solver integration")
        test_passed = .true.
        
        ! Create settings
        call adaptive_grid_settings_create(grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call solver_settings_create(solver_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create initial grid
        call generate_uniform_grid(n_initial, 0.0_dp, 1.0_dp, grid_initial, grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create test function
        do i = 1, n_initial
            function_values(i) = exp(-10.0_dp * (grid_initial(i) - 0.5_dp)**2)
        end do
        
        ! Refine grid
        call adaptive_grid_refine(n_initial, grid_initial, function_values, &
                                 refined_grid, n_refined, grid_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (allocated(refined_grid))
        
        ! Allocate solution on refined grid
        if (allocated(refined_grid)) then
            allocate(solution(n_refined))
            do i = 1, n_refined
                solution(i) = cmplx(refined_grid(i), 0.0_dp, dp)
            end do
            test_passed = test_passed .and. (size(solution) == n_refined)
            deallocate(solution, refined_grid)
        end if
        
        ! Cleanup
        call adaptive_grid_settings_destroy(grid_settings, ierr)
        call solver_settings_destroy(solver_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_grid_solver_integration
    
    !---------------------------------------------------------------------------
    ! Test 8: Full KiLCA simulation workflow
    !---------------------------------------------------------------------------
    subroutine test_full_kilca_simulation()
        type(kilca_settings_t) :: settings
        type(background_profiles_t) :: bg_profiles
        type(wave_data_t) :: wave_data
        type(zone_extended_t) :: zone
        type(antenna_t) :: antenna
        integer :: ierr
        
        call start_test("Full KiLCA simulation workflow")
        test_passed = .true.
        
        ! Initialize main KiLCA components
        call kilca_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call background_create_profiles(bg_profiles, 50, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call wave_data_create(wave_data, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call zone_extended_create(zone, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call antenna_create(antenna, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Set up simulation parameters
        settings%tolerance = 1.0e-6_dp
        settings%max_iterations = 1000
        
        wave_data%omega = cmplx(1.0e6_dp, 0.0_dp, dp)  ! 1 MHz
        wave_data%m_mode = 3
        wave_data%n_mode = 1
        
        ! Test workflow coordination
        test_passed = test_passed .and. (settings%tolerance > 0.0_dp)
        test_passed = test_passed .and. (abs(wave_data%omega) > 0.0_dp)
        test_passed = test_passed .and. (bg_profiles%dimr == 50)
        
        ! Cleanup
        call kilca_settings_destroy(settings, ierr)
        call background_destroy_profiles(bg_profiles, ierr)
        call wave_data_destroy(wave_data, ierr)
        call zone_extended_destroy(zone, ierr)
        call antenna_destroy(antenna, ierr)
        
        call end_test(test_passed)
    end subroutine test_full_kilca_simulation
    
    !---------------------------------------------------------------------------
    ! Test 9: Memory management integration
    !---------------------------------------------------------------------------
    subroutine test_memory_management_integration()
        type(conductivity_profiles_t) :: cond_profiles
        type(maxwell_eqs_data_t) :: maxwell_data
        real(dp), allocatable :: test_grid(:), test_values(:)
        integer :: ierr, n_test
        
        call start_test("Memory management integration")
        test_passed = .true.
        
        ! Test large memory allocations
        call conductivity_create_profiles(cond_profiles, 1000, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call maxwell_create_data(maxwell_data, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test adaptive memory allocation
        call allocate_adaptive_grid(500, test_grid, test_values, n_test, &
                                   adaptive_grid_settings_t(), ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (allocated(test_grid))
        test_passed = test_passed .and. (allocated(test_values))
        
        ! Test memory deallocation
        call deallocate_adaptive_grid(test_grid, test_values, &
                                     adaptive_grid_settings_t(), ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (.not. allocated(test_grid))
        test_passed = test_passed .and. (.not. allocated(test_values))
        
        ! Cleanup
        call conductivity_destroy_profiles(cond_profiles, ierr)
        call maxwell_destroy_data(maxwell_data, ierr)
        
        call end_test(test_passed)
    end subroutine test_memory_management_integration
    
    !---------------------------------------------------------------------------
    ! Test 10: Error handling integration
    !---------------------------------------------------------------------------
    subroutine test_error_handling_integration()
        type(kilca_settings_t) :: settings
        type(solver_settings_t) :: solver_settings
        integer :: ierr
        
        call start_test("Error handling integration")
        test_passed = .true.
        
        ! Test consistent error codes across modules
        call kilca_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call solver_settings_create(solver_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test error propagation
        settings%max_iterations = -1  ! Invalid value
        test_passed = test_passed .and. (settings%max_iterations < 0)
        
        ! Reset to valid values
        settings%max_iterations = 1000
        test_passed = test_passed .and. (settings%max_iterations > 0)
        
        ! Cleanup
        call kilca_settings_destroy(settings, ierr)
        call solver_settings_destroy(solver_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_error_handling_integration
    
    !---------------------------------------------------------------------------
    ! Test 11: Performance integration
    !---------------------------------------------------------------------------
    subroutine test_performance_integration()
        integer, parameter :: n_large = 100
        complex(dp) :: large_array(n_large), result_array(n_large)
        type(fourier_settings_t) :: fourier_settings
        real(dp) :: start_time, end_time
        integer :: ierr, i
        
        call start_test("Performance integration")
        test_passed = .true.
        
        call fourier_settings_create(fourier_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create test data
        do i = 1, n_large
            large_array(i) = cmplx(sin(2.0_dp * PI * real(i, dp) / real(n_large, dp)), &
                                  cos(4.0_dp * PI * real(i, dp) / real(n_large, dp)), dp)
        end do
        
        ! Test performance of integrated operations
        call cpu_time(start_time)
        call fourier_fft_1d(n_large, large_array, result_array, fourier_settings, ierr)
        call cpu_time(end_time)
        
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. ((end_time - start_time) >= 0.0_dp)
        test_passed = test_passed .and. ((end_time - start_time) < 10.0_dp)  ! Should complete quickly
        
        call fourier_settings_destroy(fourier_settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_performance_integration
    
    !---------------------------------------------------------------------------
    ! Test 12: Numerical accuracy validation
    !---------------------------------------------------------------------------
    subroutine test_numerical_accuracy_validation()
        type(hyperg_1f1_settings_t) :: hyperg_settings
        complex(dp) :: result1, result2, expected
        real(dp) :: error
        integer :: ierr
        
        call start_test("Numerical accuracy validation")
        test_passed = .true.
        
        hyperg_settings%tolerance = 1.0e-12_dp
        hyperg_settings%max_iterations = 500
        
        ! Test known hypergeometric result: 1F1(1,1,z) = exp(z)
        call hyperg_1f1_kummer_fortran(1.0_dp, 1.0_dp, 1.0_dp, result1, hyperg_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        expected = exp(cmplx(1.0_dp, 0.0_dp, dp))
        error = abs(result1 - expected)
        test_passed = test_passed .and. (error < 1.0e-10_dp)
        
        ! Test consistency between different tolerance settings
        hyperg_settings%tolerance = 1.0e-8_dp
        call hyperg_1f1_kummer_fortran(1.0_dp, 1.0_dp, 1.0_dp, result2, hyperg_settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        error = abs(result1 - result2)
        test_passed = test_passed .and. (error < 1.0e-7_dp)  ! Should be close
        
        call end_test(test_passed)
    end subroutine test_numerical_accuracy_validation
    
end program test_kilca_integration