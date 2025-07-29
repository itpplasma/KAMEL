!> @file example_basic_usage.f90
!> @brief Basic usage example for KiLCA settings module
!> @details This example demonstrates the basic usage patterns for the
!>          kilca_settings_m module, including initialization, customization,
!>          validation, and cleanup.

program example_basic_usage
    use iso_fortran_env, only: real64
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    ! Variables
    type(settings_t), pointer :: settings
    integer :: ierr
    logical :: is_valid
    character(len=1024) :: error_msg
    
    print *, "=== KiLCA Settings Module - Basic Usage Example ==="
    print *, ""
    
    ! Step 1: Initialize settings with defaults
    print *, "1. Initializing settings with defaults..."
    call settings_initialize_defaults(settings, "/tmp/kilca_example", ierr)
    if (ierr /= KILCA_SUCCESS) then
        print *, "ERROR: Failed to initialize settings, code:", ierr
        stop 1
    end if
    print *, "   ✓ Settings initialized successfully"
    
    ! Step 2: Customize antenna settings
    print *, ""
    print *, "2. Customizing antenna settings..."
    settings%antenna_settings%ra = 10.5_dp      ! Antenna radius (cm)
    settings%antenna_settings%wa = 1.2_dp       ! Current layer width (cm)
    settings%antenna_settings%I0 = 1000.0_dp    ! Antenna current (statamps)
    settings%antenna_settings%flab = cmplx(50.0e6_dp, 2000.0_dp)  ! Frequency (Hz)
    settings%antenna_settings%flag_debug = 1    ! Enable debugging
    
    ! Set up antenna modes (m,n pairs)
    settings%antenna_settings%dma = 2
    allocate(settings%antenna_settings%modes(4))  ! 2*dma elements
    settings%antenna_settings%modes = [1, 1, 2, 2]  ! (m=1,n=1), (m=2,n=2)
    
    print *, "   ✓ Antenna radius:", settings%antenna_settings%ra, "cm"
    print *, "   ✓ Antenna current:", settings%antenna_settings%I0, "statamps"
    print *, "   ✓ Laboratory frequency:", real(settings%antenna_settings%flab)/1.0e6_dp, "MHz"
    
    ! Step 3: Customize background plasma settings
    print *, ""
    print *, "3. Customizing background plasma settings..."
    settings%background_settings%rtor = 625.0_dp    ! Major radius (cm)
    settings%background_settings%rp = 200.0_dp      ! Minor radius (cm)
    settings%background_settings%B0 = 25000.0_dp    ! Toroidal field (G)
    settings%background_settings%calc_back = 1      ! Calculate background
    
    print *, "   ✓ Major radius:", settings%background_settings%rtor, "cm"
    print *, "   ✓ Minor radius:", settings%background_settings%rp, "cm"
    print *, "   ✓ Toroidal field:", settings%background_settings%B0/10000.0_dp, "T"
    
    ! Step 4: Configure output settings
    print *, ""
    print *, "4. Configuring output settings..."
    settings%output_settings%flag_background = 1   ! Output background data
    settings%output_settings%flag_emfield = 1      ! Output EM field data
    settings%output_settings%flag_additional = 1   ! Output additional quantities
    settings%output_settings%num_quants = 3
    
    allocate(settings%output_settings%flag_quants(3))
    settings%output_settings%flag_quants = [1, 0, 1]  ! Select specific quantities
    
    print *, "   ✓ Output flags configured"
    print *, "   ✓ Additional quantities:", settings%output_settings%num_quants
    
    ! Step 5: Set up eigenmode calculation (optional)
    print *, ""
    print *, "5. Setting up eigenmode calculation..."
    settings%eigmode_settings%search_flag = 1       ! Enable eigenmode search
    settings%eigmode_settings%rdim = 100           ! Real frequency mesh points
    settings%eigmode_settings%idim = 100           ! Imaginary frequency mesh points
    settings%eigmode_settings%rfmin = 45.0e6_dp    ! Min real freq (Hz)
    settings%eigmode_settings%rfmax = 55.0e6_dp    ! Max real freq (Hz)
    settings%eigmode_settings%n_zeros = 10         ! Number of zeros to find
    
    ! Set output filename
    settings%eigmode_settings%fname = "/tmp/kilca_example/eigenmode_results.dat"
    
    print *, "   ✓ Eigenmode search enabled"
    print *, "   ✓ Frequency range:", settings%eigmode_settings%rfmin/1.0e6_dp, "-", &
             settings%eigmode_settings%rfmax/1.0e6_dp, "MHz"
    
    ! Step 6: Validate the configuration
    print *, ""
    print *, "6. Validating configuration..."
    
    ! Validate individual components
    call antenna_validate(settings%antenna_settings, is_valid, error_msg, ierr)
    if (.not. is_valid) then
        print *, "   ✗ Antenna validation failed:", trim(error_msg)
        stop 1
    end if
    print *, "   ✓ Antenna settings valid"
    
    call back_sett_validate(settings%background_settings, is_valid, error_msg, ierr)
    if (.not. is_valid) then
        print *, "   ✗ Background validation failed:", trim(error_msg)
        stop 1
    end if
    print *, "   ✓ Background settings valid"
    
    ! Validate complete settings
    call settings_validate_complete(settings, is_valid, error_msg, ierr)
    if (.not. is_valid) then
        print *, "   ✗ Complete validation failed:", trim(error_msg)
        stop 1
    end if
    print *, "   ✓ Complete settings valid"
    
    ! Check consistency between components
    call settings_validate_consistency(settings, is_valid, error_msg, ierr)
    if (.not. is_valid) then
        print *, "   ✗ Consistency check failed:", trim(error_msg)
        stop 1
    end if
    print *, "   ✓ Settings are consistent"
    
    ! Step 7: Display the complete configuration
    print *, ""
    print *, "7. Complete configuration summary:"
    print *, "   " // repeat("=", 60)
    call settings_print_all(settings, ierr)
    print *, "   " // repeat("=", 60)
    
    ! Step 8: Demonstrate copying settings
    print *, ""
    print *, "8. Demonstrating settings copy..."
    block
        type(settings_t), pointer :: settings_copy
        logical :: are_equal
        
        call settings_deep_copy(settings, settings_copy, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "   ✗ Failed to copy settings"
        else
            print *, "   ✓ Settings copied successfully"
            
            ! Compare original and copy
            call settings_compare(settings, settings_copy, are_equal, ierr)
            if (are_equal) then
                print *, "   ✓ Copy is identical to original"
            else
                print *, "   ✗ Copy differs from original"
            end if
            
            ! Clean up copy
            call settings_destroy(settings_copy, ierr)
        end if
    end block
    
    ! Step 9: Clean up
    print *, ""
    print *, "9. Cleaning up..."
    call settings_destroy(settings, ierr)
    if (ierr /= KILCA_SUCCESS) then
        print *, "   ✗ Error during cleanup, code:", ierr
    else
        print *, "   ✓ Memory cleaned up successfully"
    end if
    
    print *, ""
    print *, "=== Basic Usage Example Completed Successfully ==="
    print *, ""
    print *, "Next steps:"
    print *, "  - See example_advanced_usage.f90 for more complex scenarios"
    print *, "  - Read SETTINGS_USAGE.md for detailed usage patterns"
    print *, "  - Check PROCEDURE_REFERENCE.md for complete API documentation"
    
end program example_basic_usage