!> @file example_basic_usage.f90
!> @brief Basic usage example for KiLCA settings module
!> @details Demonstrates the fundamental workflow for creating, initializing,
!>          configuring, and using settings in a typical KiLCA simulation

program example_basic_usage
    use iso_fortran_env, only: real64, int32
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    type(settings_t), pointer :: sd, sd_copy
    integer :: ierr
    logical :: is_valid, is_valid_check, recovered, is_equal
    character(len=1024) :: error_msg, validation_error
    
    print *, "=== KiLCA Settings Basic Usage Example ==="
    print *
    
    ! Step 1: Create settings structure with project path
    print *, "1. Creating settings structure..."
    call settings_create(sd, "/tmp/kilca_example/", ierr)
    if (ierr /= KILCA_SUCCESS) then
        print *, "ERROR: Failed to create settings structure"
        stop 1
    end if
    print *, "   Settings structure created successfully"
    
    ! Step 2: Initialize with scientific defaults
    print *, "2. Initializing with default values..."
    call settings_initialize_defaults(sd, "/tmp/kilca_example/", ierr)
    if (ierr /= KILCA_SUCCESS) then
        print *, "ERROR: Failed to initialize defaults"
        call settings_destroy(sd, ierr)
        stop 1
    end if
    print *, "   Default initialization completed"
    
    ! Step 3: Configure antenna settings for a typical ASDEX Upgrade shot
    print *, "3. Configuring antenna settings..."
    sd%antenna_settings%ra = 50.0_dp          ! Antenna radius (cm)
    sd%antenna_settings%wa = 2.0_dp           ! Current layer width
    sd%antenna_settings%I0 = 800.0_dp         ! Current (statamp)
    sd%antenna_settings%flab = (1.2e5_dp, 0.0_dp)  ! Frequency 120 kHz
    
    ! Set antenna modes (m=1,n=1 and m=2,n=1)
    sd%antenna_settings%dma = 2  ! Number of mode pairs
    if (allocated(sd%antenna_settings%modes)) deallocate(sd%antenna_settings%modes)
    allocate(sd%antenna_settings%modes(4))  ! 2 * dma
    sd%antenna_settings%modes = [1, 1, 2, 1]  ! [m1, n1, m2, n2]
    print *, "   Antenna configured: ra=", sd%antenna_settings%ra, " cm"
    print *, "   Frequency:", real(sd%antenna_settings%flab), " Hz"
    print *, "   Modes: (1,1), (2,1)"
    
    ! Step 4: Configure background settings for ASDEX Upgrade
    print *, "4. Configuring background settings..."
    sd%background_settings%rtor = 625.0_dp    ! Major radius (cm)
    sd%background_settings%rp = 200.0_dp      ! Plasma radius (cm)
    sd%background_settings%B0 = 25000.0_dp    ! Magnetic field (G)
    sd%background_settings%m_i = 2.0_dp       ! Deuterium mass
    sd%background_settings%N = 5              ! Spline order
    
    ! Set background type
    if (allocated(sd%background_settings%flag_back)) &
        deallocate(sd%background_settings%flag_back)
    sd%background_settings%flag_back = "normal"
    print *, "   Background configured: rtor=", sd%background_settings%rtor, " cm"
    print *, "   Magnetic field:", sd%background_settings%B0, " G"
    print *, "   Ion mass:", sd%background_settings%m_i, " proton masses"
    
    ! Step 5: Configure output settings
    print *, "5. Configuring output settings..."
    sd%output_settings%flag_background = 1    ! Compute background
    sd%output_settings%flag_emfield = 1       ! Compute EM fields
    sd%output_settings%flag_additional = 0    ! Skip additional quantities
    sd%output_settings%flag_dispersion = 1    ! Compute dispersion
    print *, "   Output flags set: background=", sd%output_settings%flag_background, &
             ", emfield=", sd%output_settings%flag_emfield
    
    ! Step 6: Configure eigenmode settings  
    print *, "6. Configuring eigenmode settings..."
    sd%eigmode_settings%search_flag = 1       ! Enable eigenmode search
    sd%eigmode_settings%rdim = 150            ! Real frequency points
    sd%eigmode_settings%idim = 150            ! Imaginary frequency points
    sd%eigmode_settings%rfmin = 0.0_dp        ! Real frequency min (Hz)
    sd%eigmode_settings%rfmax = 2.0e5_dp      ! Real frequency max (Hz)
    sd%eigmode_settings%ifmin = -1.0e4_dp     ! Imaginary frequency min (Hz)
    sd%eigmode_settings%ifmax = 1.0e4_dp      ! Imaginary frequency max (Hz)
    sd%eigmode_settings%eps_res = 1.0e-8_dp   ! Residual tolerance
    
    ! Set eigenmode output file
    if (allocated(sd%eigmode_settings%fname)) &
        deallocate(sd%eigmode_settings%fname)
    sd%eigmode_settings%fname = "example_eigenmode_output.dat"
    print *, "   Eigenmode search configured with", sd%eigmode_settings%rdim, "x", &
             sd%eigmode_settings%idim, " frequency grid"
    
    ! Step 7: Validate all settings
    print *, "7. Validating settings..."
    call settings_validate_complete(sd, is_valid_check, validation_error, ierr)
    if (ierr /= KILCA_SUCCESS .or. .not. is_valid_check) then
        print *, "ERROR: Settings validation failed:", trim(validation_error)
        call settings_destroy(sd, ierr)
        stop 1
    end if
    print *, "   All settings validated successfully"
    
    ! Step 8: Print complete settings summary
    print *, "8. Settings summary:"
    print *, "   ================================"
    call settings_print_all(sd, ierr)
    
    ! Step 9: Demonstrate error handling
    print *, "9. Demonstrating error handling..."
    
    ! Test with invalid value
    sd%antenna_settings%ra = -10.0_dp  ! Invalid negative radius
    call antenna_validate(sd%antenna_settings, is_valid, error_msg, ierr)
    if (.not. is_valid) then
        print *, "   Detected invalid antenna radius:", trim(error_msg)
    end if
    
    ! Attempt recovery
    call settings_attempt_recovery(sd%antenna_settings, recovered, ierr)
    if (recovered) then
        print *, "   Successfully recovered from invalid settings"
        print *, "   New antenna radius:", sd%antenna_settings%ra
    end if
    
    ! Step 10: Demonstrate settings copying
    print *, "10. Demonstrating settings copying..."
    
    call settings_create(sd_copy, sd%path2project, ierr)
    call settings_deep_copy(sd, sd_copy, ierr)
    
    call settings_compare(sd, sd_copy, is_equal, ierr)
    if (is_equal) then
        print *, "    Copy created successfully - settings are identical"
    end if
    
    ! Clean up copy
    call settings_destroy(sd_copy, ierr)
    
    ! Step 11: Clean up
    print *, "11. Cleaning up..."
    call settings_destroy(sd, ierr)
    if (ierr /= KILCA_SUCCESS) then
        print *, "WARNING: Error during cleanup:", ierr
    else
        print *, "    Cleanup completed successfully"
    end if
    
    print *
    print *, "=== Basic Usage Example Completed Successfully ==="
    
end program example_basic_usage