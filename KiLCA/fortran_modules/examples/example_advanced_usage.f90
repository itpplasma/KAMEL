!> @file example_advanced_usage.f90
!> @brief Advanced usage example for KiLCA settings module
!> @details This example demonstrates advanced usage patterns including
!>          error handling, recovery, parameter sweeps, and robust configuration.

program example_advanced_usage
    use iso_fortran_env, only: real64
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    print *, "=== KiLCA Settings Module - Advanced Usage Example ==="
    print *, ""
    
    ! Demonstrate various advanced patterns
    call demonstrate_error_handling()
    call demonstrate_parameter_sweep()
    call demonstrate_robust_configuration()
    call demonstrate_performance_patterns()
    
    print *, ""
    print *, "=== Advanced Usage Example Completed Successfully ==="
    
contains

    !> @brief Demonstrate comprehensive error handling and recovery
    subroutine demonstrate_error_handling()
        type(settings_t), pointer :: settings
        integer :: ierr
        character(len=1024) :: error_msg
        character(len=256) :: error_name
        logical :: recovered
        
        print *, "1. Demonstrating Error Handling and Recovery"
        print *, "   " // repeat("-", 50)
        
        ! Example 1: Handle creation failure
        print *, "   a) Testing creation with invalid path..."
        call settings_create(settings, "", ierr)  ! Empty path should fail
        
        if (ierr /= KILCA_SUCCESS) then
            call settings_get_error_name(ierr, error_name, ierr)
            call settings_format_error_message(ierr, "settings creation", "empty path", error_msg)
            
            print *, "      ✓ Caught expected error:"
            print *, "        Code:", ierr, "(" // trim(error_name) // ")"
            print *, "        Message:", trim(error_msg)
            
            ! Log the error
            call settings_log_error(ierr, "settings creation", "empty path", "error_demo.log", ierr)
            print *, "      ✓ Error logged to error_demo.log"
        end if
        
        ! Example 2: Create valid settings and test recovery
        print *, ""
        print *, "   b) Testing error recovery..."
        call settings_create(settings, "/tmp/kilca_advanced", ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "      ✗ Failed to create valid settings"
            return
        end if
        
        ! Create invalid antenna state
        settings%antenna_settings%ra = -10.0_dp     ! Invalid negative radius
        settings%antenna_settings%dma = 2           ! But no modes array allocated
        
        print *, "      - Created invalid antenna state (ra < 0, missing modes array)"
        
        ! Attempt recovery
        call settings_attempt_recovery(settings%antenna_settings, recovered, ierr)
        if (recovered) then
            print *, "      ✓ Successfully recovered from invalid state:"
            print *, "        Fixed ra =", settings%antenna_settings%ra
            if (allocated(settings%antenna_settings%modes)) then
                print *, "        Allocated modes array with size", size(settings%antenna_settings%modes)
            end if
        else
            print *, "      ✗ Could not recover from invalid state"
        end if
        
        ! Example 3: Detailed validation errors
        print *, ""
        print *, "   c) Testing detailed error reporting..."
        
        ! Set up invalid background settings
        settings%background_settings%rtor = -100.0_dp   ! Invalid
        settings%background_settings%rp = 1000.0_dp     ! Invalid (larger than rtor)
        settings%background_settings%B0 = -500.0_dp     ! Invalid
        settings%background_settings%N = 4              ! Invalid (must be odd)
        
        call settings_get_detailed_validation_errors(settings%background_settings, error_msg, ierr)
        if (ierr == KILCA_SUCCESS) then
            print *, "      ✓ Detailed validation errors:"
            print *, "        ", trim(error_msg)
        end if
        
        call settings_destroy(settings, ierr)
        print *, "   ✓ Error handling demonstration completed"
    end subroutine demonstrate_error_handling
    
    !> @brief Demonstrate parameter sweep functionality
    subroutine demonstrate_parameter_sweep()
        type(settings_t), pointer :: base_settings, sweep_settings
        integer :: ierr, i
        real(dp), parameter :: frequencies(5) = [25.0e6_dp, 30.0e6_dp, 35.0e6_dp, 40.0e6_dp, 45.0e6_dp]
        character(len=256) :: filename
        
        print *, ""
        print *, "2. Demonstrating Parameter Sweep"
        print *, "   " // repeat("-", 50)
        
        ! Create base configuration
        call settings_initialize_defaults(base_settings, "/tmp/kilca_sweep", ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "   ✗ Failed to create base settings"
            return
        end if
        
        ! Configure common parameters
        base_settings%antenna_settings%ra = 15.0_dp
        base_settings%antenna_settings%wa = 1.0_dp
        base_settings%antenna_settings%I0 = 800.0_dp
        base_settings%background_settings%rtor = 165.0_dp
        base_settings%background_settings%rp = 50.0_dp
        base_settings%background_settings%B0 = 25000.0_dp
        
        print *, "   Base configuration created"
        
        ! Perform frequency sweep
        print *, "   Performing frequency sweep:"
        do i = 1, size(frequencies)
            ! Deep copy base settings
            call settings_deep_copy(base_settings, sweep_settings, ierr)
            if (ierr /= KILCA_SUCCESS) then
                print *, "      ✗ Failed to copy settings for frequency", frequencies(i)/1.0e6_dp, "MHz"
                cycle
            end if
            
            ! Modify frequency
            sweep_settings%antenna_settings%flab = cmplx(frequencies(i), 0.0_dp)
            
            ! Update output filename
            write(filename, '("freq_sweep_",I0,"MHz.dat")') int(frequencies(i)/1.0e6_dp)
            sweep_settings%eigmode_settings%fname = trim(sweep_settings%path2project) // "/" // trim(filename)
            
            ! Set eigenmode search around this frequency
            sweep_settings%eigmode_settings%search_flag = 1
            sweep_settings%eigmode_settings%rfmin = frequencies(i) * 0.9_dp
            sweep_settings%eigmode_settings%rfmax = frequencies(i) * 1.1_dp
            
            print *, "      ✓ Configured for", frequencies(i)/1.0e6_dp, "MHz ->", trim(filename)
            
            ! In a real application, you would run the calculation here
            ! call run_kilca_calculation(sweep_settings)
            
            call settings_destroy(sweep_settings, ierr)
        end do
        
        call settings_destroy(base_settings, ierr)
        print *, "   ✓ Parameter sweep demonstration completed"
    end subroutine demonstrate_parameter_sweep
    
    !> @brief Demonstrate robust configuration with retries and validation
    subroutine demonstrate_robust_configuration()
        type(settings_t), pointer :: settings
        logical :: success
        integer :: ierr, attempt
        character(len=1024) :: error_msg
        logical :: is_valid
        
        print *, ""
        print *, "3. Demonstrating Robust Configuration"
        print *, "   " // repeat("-", 50)
        
        ! Robust setup with multiple attempts
        success = .false.
        
        do attempt = 1, 3
            print *, "   Attempt", attempt, "to create robust configuration..."
            
            ! Try to create and configure settings
            call settings_initialize_defaults(settings, "/tmp/kilca_robust", ierr)
            if (ierr /= KILCA_SUCCESS) then
                print *, "      ✗ Failed to initialize settings, retrying..."
                cycle
            end if
            
            ! Configure realistic ASDEX Upgrade parameters
            settings%antenna_settings%ra = 18.5_dp       ! Close to plasma boundary
            settings%antenna_settings%wa = 0.8_dp        ! Narrow current layer
            settings%antenna_settings%I0 = 800.0_dp      ! Moderate current
            settings%antenna_settings%flab = cmplx(30.0e6_dp, 0.0_dp)  ! 30 MHz ICRF
            
            ! ASDEX Upgrade dimensions
            settings%background_settings%rtor = 165.0_dp  ! Major radius
            settings%background_settings%rp = 50.0_dp     ! Minor radius
            settings%background_settings%B0 = 25000.0_dp  ! 2.5 T field
            
            ! Set up antenna modes for ICRF
            settings%antenna_settings%dma = 2
            allocate(settings%antenna_settings%modes(4))
            settings%antenna_settings%modes = [0, 15, 0, 16]  ! n=15,16
            
            ! Comprehensive validation
            call settings_validate_complete(settings, is_valid, error_msg, ierr)
            if (ierr /= KILCA_SUCCESS) then
                print *, "      ✗ Validation failed with error code:", ierr
                call settings_destroy(settings, ierr)
                cycle
            end if
            
            if (.not. is_valid) then
                print *, "      ✗ Settings invalid:", trim(error_msg)
                call settings_destroy(settings, ierr)
                cycle
            end if
            
            ! Check consistency
            call settings_validate_consistency(settings, is_valid, error_msg, ierr)
            if (.not. is_valid) then
                print *, "      ✗ Consistency check failed:", trim(error_msg)
                call settings_destroy(settings, ierr)
                cycle
            end if
            
            ! Success!
            success = .true.
            print *, "      ✓ Robust configuration created successfully"
            print *, "        - Antenna at r =", settings%antenna_settings%ra, "cm"
            print *, "        - ICRF frequency =", real(settings%antenna_settings%flab)/1.0e6_dp, "MHz"
            print *, "        - Major radius =", settings%background_settings%rtor, "cm"
            print *, "        - Toroidal field =", settings%background_settings%B0/10000.0_dp, "T"
            exit
        end do
        
        if (.not. success) then
            print *, "   ✗ Failed to create robust configuration after", attempt, "attempts"
        else
            print *, "   ✓ Robust configuration validated and ready for use"
            call settings_destroy(settings, ierr)
        end if
    end subroutine demonstrate_robust_configuration
    
    !> @brief Demonstrate performance-oriented patterns
    subroutine demonstrate_performance_patterns()
        type(settings_t), pointer :: template_settings, run_settings
        integer :: ierr, i
        real(dp) :: current_values(10)
        
        print *, ""
        print *, "4. Demonstrating Performance Patterns"
        print *, "   " // repeat("-", 50)
        
        ! Create template for multiple runs
        call settings_initialize_defaults(template_settings, "/tmp/kilca_perf", ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "   ✗ Failed to create template settings"
            return
        end if
        
        ! Configure template with common parameters
        template_settings%background_settings%rtor = 165.0_dp
        template_settings%background_settings%rp = 50.0_dp
        template_settings%background_settings%B0 = 25000.0_dp
        template_settings%antenna_settings%ra = 15.0_dp
        template_settings%antenna_settings%wa = 1.0_dp
        template_settings%antenna_settings%flab = cmplx(30.0e6_dp, 0.0_dp)
        
        print *, "   Template configuration created"
        
        ! Generate current values for parameter scan
        do i = 1, 10
            current_values(i) = 500.0_dp + i * 100.0_dp  ! 600-1500 statamps
        end do
        
        print *, "   Performing efficient parameter scan using template copying:"
        
        ! Efficient parameter scan using template
        do i = 1, 10
            ! Deep copy from template (much faster than re-initializing)
            call settings_deep_copy(template_settings, run_settings, ierr)
            if (ierr /= KILCA_SUCCESS) then
                print *, "      ✗ Failed to copy template for run", i
                cycle
            end if
            
            ! Modify run-specific parameters
            run_settings%antenna_settings%I0 = current_values(i)
            
            ! Update output paths
            block
                character(len=256) :: run_path
                write(run_path, '("run_",I0)') i
                run_settings%eigmode_settings%fname = &
                    trim(run_settings%path2project) // "/" // trim(run_path) // "_results.dat"
            end block
            
            print *, "      ✓ Run", i, ": I0 =", current_values(i), "statamps"
            
            ! In a real application, you would run the calculation here
            ! call run_kilca_calculation(run_settings)
            
            ! Clean up this run
            call settings_destroy(run_settings, ierr)
        end do
        
        ! Clean up template
        call settings_destroy(template_settings, ierr)
        
        print *, "   ✓ Performance pattern demonstration completed"
        print *, "     (Template copying is much faster than repeated initialization)"
    end subroutine demonstrate_performance_patterns

end program example_advanced_usage