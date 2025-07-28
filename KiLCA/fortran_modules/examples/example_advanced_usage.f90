!> @file example_advanced_usage.f90
!> @brief Advanced usage example for KiLCA settings module
!> @details Demonstrates advanced features including file I/O, detailed error
!>          handling, custom initialization, and integration with simulation workflows

program example_advanced_usage
    use iso_fortran_env, only: real64, int32, output_unit
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    type(settings_t), pointer :: sd_main, sd_backup
    integer :: ierr, file_unit
    logical :: is_valid, recovered
    character(len=2048) :: error_msg, detailed_errors
    character(len=256) :: error_name
    
    print *, "=== KiLCA Settings Advanced Usage Example ==="
    print *
    
    ! === Advanced Initialization Strategies ===
    print *, "1. Advanced initialization strategies..."
    
    ! Custom antenna initialization
    type(antenna_t) :: custom_antenna
    call antenna_initialize_custom(custom_antenna, &
                                  ra=45.0_dp, wa=1.8_dp, I0=1200.0_dp, ierr=ierr)
    if (ierr == KILCA_SUCCESS) then
        print *, "   Custom antenna initialized: ra=", custom_antenna%ra, &
                ", wa=", custom_antenna%wa, ", I0=", custom_antenna%I0
    end if
    
    ! === Comprehensive Error Handling ===
    print *, "2. Comprehensive error handling..."
    
    ! Create main settings
    call settings_create(sd_main, "/tmp/advanced_example/", ierr)
    call settings_initialize_defaults(sd_main, "/tmp/advanced_example/", ierr)
    
    ! Introduce multiple validation errors
    sd_main%background_settings%rtor = -100.0_dp    ! Invalid
    sd_main%background_settings%rp = 1000.0_dp      ! Invalid (larger than rtor)  
    sd_main%background_settings%B0 = -500.0_dp      ! Invalid
    sd_main%background_settings%N = 4               ! Invalid (must be odd)
    
    ! Get detailed validation errors
    call settings_get_detailed_validation_errors(sd_main%background_settings, &
                                                 detailed_errors, ierr)
    if (len_trim(detailed_errors) > 0) then
        print *, "   Detailed validation errors found:"
        print *, "   ", trim(detailed_errors)
    end if
    
    ! Format specific error messages
    call settings_format_error_message(KILCA_ERROR_INVALID_INPUT, &
                                      "background validation", &
                                      "magnetic field B0", error_msg)
    print *, "   Formatted error:", trim(error_msg)
    
    ! Get error code names
    call settings_get_error_name(KILCA_ERROR_INVALID_INPUT, error_name, ierr)
    print *, "   Error code name:", trim(error_name)
    
    ! === Settings Recovery ===
    print *, "3. Settings recovery mechanisms..."
    
    ! Attempt recovery from invalid background settings
    call settings_attempt_recovery(sd_main%background_settings, recovered, ierr)
    if (recovered) then
        print *, "   Background settings recovered successfully"
        print *, "   New rtor:", sd_main%background_settings%rtor
        print *, "   New B0:", sd_main%background_settings%B0
    else
        print *, "   Could not automatically recover background settings"
        ! Manual recovery
        sd_main%background_settings%rtor = 625.0_dp
        sd_main%background_settings%rp = 200.0_dp
        sd_main%background_settings%B0 = 25000.0_dp
        sd_main%background_settings%N = 5
        print *, "   Applied manual recovery"
    end if
    
    ! === Error Logging ===
    print *, "4. Error logging to file..."
    
    ! Log errors to file
    call settings_log_error(KILCA_ERROR_INVALID_INPUT, &
                           "advanced example", &
                           "background B0 field", &
                           "advanced_example_errors.log", ierr)
    
    logical :: log_exists
    call settings_check_error_log_exists("advanced_example_errors.log", log_exists, ierr)
    if (log_exists) then
        print *, "   Error log created successfully"
    end if
    
    ! === File I/O Operations ===
    print *, "5. File I/O operations..."
    
    ! Save settings to file
    open(newunit=file_unit, file="advanced_settings_output.txt", action="write")
    call antenna_print_settings_to_unit(sd_main%antenna_settings, file_unit, ierr)
    call back_sett_print_settings_to_unit(sd_main%background_settings, file_unit, ierr)
    call output_sett_print_settings_to_unit(sd_main%output_settings, file_unit, ierr)
    call eigmode_sett_print_settings_to_unit(sd_main%eigmode_settings, file_unit, ierr)
    close(file_unit)
    print *, "   Settings saved to advanced_settings_output.txt"
    
    ! === Advanced Antenna Configuration ===
    print *, "6. Advanced antenna configuration..."
    
    ! Configure multi-mode antenna
    sd_main%antenna_settings%ra = 48.5_dp
    sd_main%antenna_settings%wa = 2.2_dp
    sd_main%antenna_settings%I0 = 950.0_dp
    sd_main%antenna_settings%flab = (8.5e4_dp, 1.2e3_dp)  ! Complex frequency
    
    ! Set multiple modes: (1,1), (2,1), (3,1), (1,2)
    sd_main%antenna_settings%dma = 8
    if (allocated(sd_main%antenna_settings%modes)) &
        deallocate(sd_main%antenna_settings%modes)
    allocate(sd_main%antenna_settings%modes(8))
    sd_main%antenna_settings%modes = [1, 1, 2, 1, 3, 1, 1, 2]
    
    sd_main%antenna_settings%flag_eigmode = 1  ! Enable eigenmode search
    sd_main%antenna_settings%flag_debug = 1    ! Enable debug output
    
    print *, "   Multi-mode antenna configured with", sd_main%antenna_settings%dma/2, "modes"
    print *, "   Complex frequency: real=", real(sd_main%antenna_settings%flab), &
             ", imag=", aimag(sd_main%antenna_settings%flab)
    
    ! === Advanced Background Configuration ===
    print *, "7. Advanced background configuration..."
    
    ! Set particle species (deuterium plasma)
    if (allocated(sd_main%background_settings%mass)) &
        deallocate(sd_main%background_settings%mass)
    if (allocated(sd_main%background_settings%charge)) &
        deallocate(sd_main%background_settings%charge)
    
    allocate(sd_main%background_settings%mass(2))     ! ions, electrons
    allocate(sd_main%background_settings%charge(2))   ! ions, electrons
    
    sd_main%background_settings%mass = [2.0_dp, 1.0_dp/1836.0_dp]  ! Deuterium, electron
    sd_main%background_settings%charge = [1.0_dp, -1.0_dp]         ! Ion, electron charges
    
    ! Advanced background parameters
    sd_main%background_settings%V_gal_sys = 0.02_dp    ! Moving frame velocity
    sd_main%background_settings%V_scale = 0.95_dp      ! Velocity scale factor
    sd_main%background_settings%zele = 0.8_dp          ! Reduced electron collisions
    sd_main%background_settings%zion = 1.2_dp          ! Enhanced ion collisions
    sd_main%background_settings%huge_factor = 1.0e25_dp ! Reduced huge factor
    
    print *, "   Advanced background: Deuterium plasma with collision modifications"
    print *, "   Velocity scale:", sd_main%background_settings%V_scale
    print *, "   Collision factors: zele=", sd_main%background_settings%zele, &
             ", zion=", sd_main%background_settings%zion
    
    ! === Advanced Output Configuration ===
    print *, "8. Advanced output configuration..."
    
    ! Configure selective quantity output
    sd_main%output_settings%num_quants = 6
    if (allocated(sd_main%output_settings%flag_quants)) &
        deallocate(sd_main%output_settings%flag_quants)
    allocate(sd_main%output_settings%flag_quants(6))
    sd_main%output_settings%flag_quants = [1, 1, 0, 1, 0, 1]  ! Selective output
    
    sd_main%output_settings%flag_background = 2    ! Compute and store
    sd_main%output_settings%flag_emfield = 2       ! Compute and store
    sd_main%output_settings%flag_additional = 1    ! Compute additional quantities
    sd_main%output_settings%flag_dispersion = 2    ! Compute and store dispersion
    sd_main%output_settings%flag_debug = 1         ! Enable debug output
    
    print *, "   Advanced output: selective quantities with storage"
    print *, "   Quantity flags:", sd_main%output_settings%flag_quants
    
    ! === Advanced Eigenmode Configuration ===
    print *, "9. Advanced eigenmode configuration..."
    
    ! High-resolution eigenmode search
    sd_main%eigmode_settings%rdim = 300
    sd_main%eigmode_settings%idim = 300
    sd_main%eigmode_settings%rfmin = 5.0e4_dp       ! Focus on TAE range
    sd_main%eigmode_settings%rfmax = 1.5e5_dp
    sd_main%eigmode_settings%ifmin = -5.0e3_dp      ! Narrow imaginary range
    sd_main%eigmode_settings%ifmax = 5.0e3_dp
    
    ! Tight convergence criteria
    sd_main%eigmode_settings%eps_res = 1.0e-10_dp
    sd_main%eigmode_settings%eps_abs = 1.0e-12_dp
    sd_main%eigmode_settings%eps_rel = 1.0e-10_dp
    sd_main%eigmode_settings%delta = 1.0e-10_dp
    
    ! Advanced solver settings
    sd_main%eigmode_settings%test_roots = 1
    sd_main%eigmode_settings%kmin = 2
    sd_main%eigmode_settings%kmax = 15
    sd_main%eigmode_settings%n_zeros = 20
    sd_main%eigmode_settings%use_winding = 1
    
    ! Set starting frequency guesses
    sd_main%eigmode_settings%Nguess = 3
    if (allocated(sd_main%eigmode_settings%fstart)) &
        deallocate(sd_main%eigmode_settings%fstart)
    allocate(sd_main%eigmode_settings%fstart(3))
    sd_main%eigmode_settings%fstart = [(8.0e4_dp, 1.0e3_dp), &
                                       (1.2e5_dp, -5.0e2_dp), &
                                       (9.5e4_dp, 2.0e3_dp)]
    
    print *, "   High-resolution eigenmode search:", sd_main%eigmode_settings%rdim, &
             "x", sd_main%eigmode_settings%idim
    print *, "   Convergence tolerance:", sd_main%eigmode_settings%eps_res
    print *, "   Starting guesses:", sd_main%eigmode_settings%Nguess
    
    ! === Settings Backup and Comparison ===
    print *, "10. Settings backup and comparison..."
    
    ! Create backup copy
    call settings_create(sd_backup, sd_main%path2project, ierr)
    call settings_deep_copy(sd_main, sd_backup, ierr)
    print *, "    Backup copy created"
    
    ! Modify original
    sd_main%antenna_settings%I0 = 1000.0_dp
    
    ! Compare
    logical :: is_equal
    call settings_compare(sd_main, sd_backup, is_equal, ierr)
    if (.not. is_equal) then
        print *, "    Settings differ from backup (expected after modification)"
    end if
    
    ! === Final Validation and Summary ===
    print *, "11. Final validation and summary..."
    
    call settings_validate_complete(sd_main, ierr)
    if (ierr == KILCA_SUCCESS) then
        print *, "    All advanced settings validated successfully"
    else
        print *, "    WARNING: Advanced settings validation failed:", ierr
    end if
    
    ! Print summary with context
    call settings_validate_with_context(sd_main, "advanced example", error_msg, ierr)
    if (ierr == KILCA_SUCCESS) then
        print *, "    Context validation successful"
    end if
    
    ! === Performance Demonstration ===
    print *, "12. Performance demonstration..."
    
    ! Multiple copy operations
    integer :: i, start_time, end_time, count_rate
    type(settings_t), pointer :: temp_settings
    
    call system_clock(start_time, count_rate)
    do i = 1, 100
        call settings_create(temp_settings, "/tmp/perf_test/", ierr)
        call settings_deep_copy(sd_main, temp_settings, ierr)
        call settings_destroy(temp_settings, ierr)
    end do
    call system_clock(end_time)
    
    print *, "    100 copy operations completed in", &
             real(end_time - start_time) / real(count_rate), "seconds"
    
    ! === Cleanup ===
    print *, "13. Cleanup and finalization..."
    
    ! Clean up error log
    call settings_clear_error_log("advanced_example_errors.log", ierr)
    
    ! Destroy settings
    call settings_destroy(sd_main, ierr)
    call settings_destroy(sd_backup, ierr)
    
    if (ierr == KILCA_SUCCESS) then
        print *, "    All resources cleaned up successfully"
    end if
    
    print *
    print *, "=== Advanced Usage Example Completed Successfully ==="
    
end program example_advanced_usage