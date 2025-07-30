program test_kilca_adaptive_grid
    use iso_fortran_env
    use kilca_types_m
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
    print *, "KiLCA Adaptive Grid System - Unit Tests"
    print *, "========================================================"
    print *, ""
    
    ! Run all test suites
    call test_grid_settings_management()
    call test_uniform_grid_generation()
    call test_adaptive_refinement_basic()
    call test_resonance_layer_refinement()
    call test_gradient_based_refinement() 
    call test_curvature_based_refinement()
    call test_error_estimation()
    call test_grid_redistribution()
    call test_boundary_handling()
    call test_grid_validation()
    call test_memory_management()
    call test_performance_optimization()
    
    ! Print summary
    print *, ""
    print *, "========================================================"
    print *, "TEST SUMMARY"
    print *, "========================================================"
    print '(A,I5)', " Total tests run:     ", total_tests
    print '(A,I5)', " Tests passed:        ", passed_tests  
    print '(A,I5)', " Tests failed:        ", failed_tests
    print '(A,F8.2,A)', " Success rate:        ", &
            100.0_dp * real(passed_tests, dp) / real(total_tests, dp), " %"
    
    if (failed_tests > 0) then
        print *, ""
        print *, " *** SOME TESTS FAILED! ***"
        stop 1
    else
        print *, ""
        print *, " *** ALL TESTS PASSED! ***"
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
    ! Test 1: Grid settings management
    !---------------------------------------------------------------------------
    subroutine test_grid_settings_management()
        type(adaptive_grid_settings_t) :: settings
        integer :: ierr
        
        call start_test("Grid settings management")
        test_passed = .true.
        
        call adaptive_grid_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test default values
        test_passed = test_passed .and. (settings%method == GRID_METHOD_GRADIENT)
        test_passed = test_passed .and. (settings%max_points == 1000)
        test_passed = test_passed .and. (settings%min_spacing > 0.0_dp)
        test_passed = test_passed .and. (settings%tolerance > 0.0_dp)
        
        ! Test setting custom values
        settings%method = GRID_METHOD_CURVATURE
        settings%max_points = 2000
        settings%resonance_refinement = .true.
        settings%debug_level = 1
        
        test_passed = test_passed .and. (settings%method == GRID_METHOD_CURVATURE)
        test_passed = test_passed .and. (settings%max_points == 2000)
        test_passed = test_passed .and. (settings%resonance_refinement)
        
        call adaptive_grid_settings_destroy(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call end_test(test_passed)
    end subroutine test_grid_settings_management
    
    !---------------------------------------------------------------------------
    ! Test 2: Uniform grid generation
    !---------------------------------------------------------------------------
    subroutine test_uniform_grid_generation()
        integer, parameter :: n = 10
        real(dp) :: r_min, r_max, grid(n)
        type(adaptive_grid_settings_t) :: settings
        integer :: ierr, i
        
        call start_test("Uniform grid generation")
        test_passed = .true.
        
        call adaptive_grid_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        r_min = 0.0_dp
        r_max = 1.0_dp
        
        call generate_uniform_grid(n, r_min, r_max, grid, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test grid properties
        test_passed = test_passed .and. (grid(1) == r_min)
        test_passed = test_passed .and. (grid(n) == r_max)
        
        ! Test monotonicity
        do i = 2, n
            test_passed = test_passed .and. (grid(i) > grid(i-1))
        end do
        
        call end_test(test_passed)
    end subroutine test_uniform_grid_generation
    
    !---------------------------------------------------------------------------
    ! Test 3: Basic adaptive refinement
    !---------------------------------------------------------------------------
    subroutine test_adaptive_refinement_basic()
        integer, parameter :: n_initial = 10
        real(dp) :: r_min, r_max, grid_initial(n_initial)
        real(dp) :: function_values(n_initial)
        real(dp), allocatable :: refined_grid(:)
        type(adaptive_grid_settings_t) :: settings
        integer :: n_refined, ierr, i
        
        call start_test("Basic adaptive refinement")
        test_passed = .true.
        
        call adaptive_grid_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        r_min = 0.0_dp
        r_max = 1.0_dp
        
        ! Create initial uniform grid
        call generate_uniform_grid(n_initial, r_min, r_max, grid_initial, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create test function with sharp gradient at r=0.5
        do i = 1, n_initial
            function_values(i) = tanh(10.0_dp * (grid_initial(i) - 0.5_dp))
        end do
        
        ! Perform adaptive refinement
        call adaptive_grid_refine(n_initial, grid_initial, function_values, &
                                 refined_grid, n_refined, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (allocated(refined_grid))
        test_passed = test_passed .and. (n_refined > n_initial)
        
        if (allocated(refined_grid)) deallocate(refined_grid)
        
        call end_test(test_passed)
    end subroutine test_adaptive_refinement_basic
    
    !---------------------------------------------------------------------------
    ! Test 4: Resonance layer refinement
    !---------------------------------------------------------------------------
    subroutine test_resonance_layer_refinement()
        integer, parameter :: n = 50
        real(dp) :: r_min, r_max, grid(n), resonance_position
        real(dp), allocatable :: refined_grid(:)
        type(adaptive_grid_settings_t) :: settings
        integer :: n_refined, ierr
        
        call start_test("Resonance layer refinement")
        test_passed = .true.
        
        call adaptive_grid_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        r_min = 0.0_dp
        r_max = 1.0_dp
        resonance_position = 0.3_dp
        
        settings%resonance_refinement = .true.
        settings%resonance_width = 0.1_dp
        settings%eps_res = 1.0e-4_dp
        settings%eps_out = 1.0e-2_dp
        settings%dr_res = 0.01_dp  ! Larger step in resonance region
        settings%dr_out = 0.05_dp  ! Larger step outside resonance region
        
        call generate_uniform_grid(n, r_min, r_max, grid, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call refine_resonance_layer(n, grid, resonance_position, refined_grid, &
                                   n_refined, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (allocated(refined_grid))
        
        if (allocated(refined_grid)) deallocate(refined_grid)
        
        call end_test(test_passed)
    end subroutine test_resonance_layer_refinement
    
    !---------------------------------------------------------------------------
    ! Test 5: Gradient-based refinement
    !---------------------------------------------------------------------------
    subroutine test_gradient_based_refinement()
        integer, parameter :: n = 20
        real(dp) :: grid(n), function_values(n), gradients(n-1)
        type(adaptive_grid_settings_t) :: settings
        integer :: ierr, i
        logical, allocatable :: refine_flags(:)
        
        call start_test("Gradient-based refinement")
        test_passed = .true.
        
        call adaptive_grid_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create test grid and function
        do i = 1, n
            grid(i) = real(i-1, dp) / real(n-1, dp)
            function_values(i) = exp(-10.0_dp * (grid(i) - 0.5_dp)**2)
        end do
        
        call compute_gradients(n, grid, function_values, gradients, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call identify_refinement_regions(n-1, gradients, refine_flags, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (allocated(refine_flags))
        
        if (allocated(refine_flags)) deallocate(refine_flags)
        
        call end_test(test_passed)
    end subroutine test_gradient_based_refinement
    
    !---------------------------------------------------------------------------
    ! Test 6: Curvature-based refinement
    !---------------------------------------------------------------------------
    subroutine test_curvature_based_refinement()
        integer, parameter :: n = 25
        real(dp) :: grid(n), function_values(n), curvatures(n-2)
        type(adaptive_grid_settings_t) :: settings
        integer :: ierr, i
        
        call start_test("Curvature-based refinement")
        test_passed = .true.
        
        call adaptive_grid_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        settings%method = GRID_METHOD_CURVATURE
        
        ! Create test grid and function with curvature features
        do i = 1, n
            grid(i) = real(i-1, dp) / real(n-1, dp)
            function_values(i) = sin(8.0_dp * PI * grid(i)) * exp(-2.0_dp * grid(i))
        end do
        
        call compute_curvatures(n, grid, function_values, curvatures, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (size(curvatures) == n-2)
        
        call end_test(test_passed)
    end subroutine test_curvature_based_refinement
    
    !---------------------------------------------------------------------------
    ! Test 7: Error estimation
    !---------------------------------------------------------------------------
    subroutine test_error_estimation()
        integer, parameter :: n = 15
        real(dp) :: grid(n), function_values(n), errors(n-1)
        type(adaptive_grid_settings_t) :: settings
        integer :: ierr, i
        
        call start_test("Error estimation")
        test_passed = .true.
        
        call adaptive_grid_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create test data
        do i = 1, n
            grid(i) = real(i-1, dp) / real(n-1, dp)
            function_values(i) = grid(i)**3 - 2.0_dp * grid(i)**2 + grid(i)
        end do
        
        call estimate_interpolation_errors(n, grid, function_values, errors, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (all(errors >= 0.0_dp))
        
        call end_test(test_passed)
    end subroutine test_error_estimation
    
    !---------------------------------------------------------------------------
    ! Test 8: Grid redistribution
    !---------------------------------------------------------------------------
    subroutine test_grid_redistribution()
        integer, parameter :: n_old = 20, n_new = 30
        real(dp) :: old_grid(n_old), old_values(n_old)
        real(dp) :: new_grid(n_new), new_values(n_new)
        type(adaptive_grid_settings_t) :: settings
        integer :: ierr, i
        
        call start_test("Grid redistribution")
        test_passed = .true.
        
        call adaptive_grid_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create old grid and values
        do i = 1, n_old
            old_grid(i) = real(i-1, dp) / real(n_old-1, dp)
            old_values(i) = cos(2.0_dp * PI * old_grid(i))
        end do
        
        ! Create new grid
        do i = 1, n_new
            new_grid(i) = real(i-1, dp) / real(n_new-1, dp)
        end do
        
        call redistribute_function_on_grid(n_old, old_grid, old_values, &
                                          n_new, new_grid, new_values, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call end_test(test_passed)
    end subroutine test_grid_redistribution
    
    !---------------------------------------------------------------------------
    ! Test 9: Boundary handling
    !---------------------------------------------------------------------------
    subroutine test_boundary_handling()
        integer, parameter :: n = 12
        real(dp) :: grid(n), r_min, r_max
        type(adaptive_grid_settings_t) :: settings
        integer :: ierr
        logical :: valid_boundaries
        
        call start_test("Boundary handling")
        test_passed = .true.
        
        call adaptive_grid_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        r_min = 0.1_dp
        r_max = 0.9_dp
        
        call generate_uniform_grid(n, r_min, r_max, grid, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        call validate_grid_boundaries(n, grid, r_min, r_max, valid_boundaries, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (valid_boundaries)
        
        call end_test(test_passed)
    end subroutine test_boundary_handling
    
    !---------------------------------------------------------------------------
    ! Test 10: Grid validation
    !---------------------------------------------------------------------------
    subroutine test_grid_validation()
        integer, parameter :: n = 8
        real(dp) :: valid_grid(n), invalid_grid(n)
        type(adaptive_grid_settings_t) :: settings
        integer :: ierr, i
        logical :: is_valid
        
        call start_test("Grid validation")
        test_passed = .true.
        
        call adaptive_grid_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create valid monotonic grid
        do i = 1, n
            valid_grid(i) = real(i-1, dp) / real(n-1, dp)
        end do
        
        call validate_grid_monotonicity(n, valid_grid, is_valid, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (is_valid)
        
        ! Create invalid non-monotonic grid
        invalid_grid = valid_grid
        invalid_grid(n/2) = invalid_grid(n/2+1) + 0.1_dp  ! Break monotonicity
        
        call validate_grid_monotonicity(n, invalid_grid, is_valid, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (.not. is_valid)
        
        call end_test(test_passed)
    end subroutine test_grid_validation
    
    !---------------------------------------------------------------------------
    ! Test 11: Memory management
    !---------------------------------------------------------------------------
    subroutine test_memory_management()
        integer, parameter :: n_max = 100
        real(dp), allocatable :: grid(:), function_vals(:)
        type(adaptive_grid_settings_t) :: settings
        integer :: ierr, n_actual
        
        call start_test("Memory management")
        test_passed = .true.
        
        call adaptive_grid_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Test dynamic allocation
        call allocate_adaptive_grid(n_max, grid, function_vals, n_actual, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (allocated(grid))
        test_passed = test_passed .and. (allocated(function_vals))
        test_passed = test_passed .and. (n_actual <= n_max)
        
        ! Test deallocation
        call deallocate_adaptive_grid(grid, function_vals, settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. (.not. allocated(grid))
        test_passed = test_passed .and. (.not. allocated(function_vals))
        
        call end_test(test_passed)
    end subroutine test_memory_management
    
    !---------------------------------------------------------------------------
    ! Test 12: Performance optimization
    !---------------------------------------------------------------------------
    subroutine test_performance_optimization()
        integer, parameter :: n = 100
        real(dp) :: grid(n), function_values(n)
        type(adaptive_grid_settings_t) :: settings
        integer :: ierr, i
        real(dp) :: start_time, end_time
        
        call start_test("Performance optimization")
        test_passed = .true.
        
        call adaptive_grid_settings_create(settings, ierr)
        test_passed = test_passed .and. (ierr == 0)
        
        ! Create test data
        do i = 1, n
            grid(i) = real(i-1, dp) / real(n-1, dp)
            function_values(i) = sin(20.0_dp * PI * grid(i)) * exp(-5.0_dp * grid(i))
        end do
        
        ! Test performance
        call cpu_time(start_time)
        call optimize_grid_performance(n, grid, function_values, settings, ierr)
        call cpu_time(end_time)
        
        test_passed = test_passed .and. (ierr == 0)
        test_passed = test_passed .and. ((end_time - start_time) >= 0.0_dp)
        
        call end_test(test_passed)
    end subroutine test_performance_optimization
    
end program test_kilca_adaptive_grid