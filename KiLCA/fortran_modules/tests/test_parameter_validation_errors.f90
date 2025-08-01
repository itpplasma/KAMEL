program test_parameter_validation_errors
    use kilca_types_m, only: dp, KILCA_SUCCESS, KILCA_ERROR_INVALID_PARAMETER, &
                             KILCA_ERROR_BOUNDS, KILCA_ERROR_FORMAT
    use kilca_settings_m, only: settings_t, settings_read_namelist, &
                                settings_validate_parameters, &
                                settings_validate_physics_constraints
    implicit none
    
    type(settings_t) :: settings
    integer :: ierr
    integer :: test_status = 0
    character(len=*), parameter :: test_path = "./"
    
    print *, "=========================================="
    print *, "Testing Parameter Validation Errors [RED PHASE - SHOULD FAIL]"
    print *, "=========================================="
    
    call test_parameter_validation_functionality()
    
    if (test_status == 0) then
        print *, ""
        print *, "Parameter validation tests PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "Parameter validation tests FAILED:", test_status, "failure(s) (expected in RED phase)"
        stop 1
    end if

contains

    !> Test parameter validation functionality
    subroutine test_parameter_validation_functionality()
        print *, "Testing parameter validation functionality..."
        
        ! Test 1: Negative radius values (physics constraint)
        print *, ""
        print *, "Test 1: Negative radius values..."
        call create_negative_radius_config()
        call test_physics_validation("negative_radius.conf", "NEGATIVE_RADIUS")
        
        ! Test 2: Invalid magnetic field values
        print *, ""
        print *, "Test 2: Invalid magnetic field values..."
        call create_invalid_bfield_config()
        call test_physics_validation("invalid_bfield.conf", "INVALID_BFIELD")
        
        ! Test 3: Inconsistent geometry (rtor < rp)
        print *, ""
        print *, "Test 3: Inconsistent geometry..."
        call create_inconsistent_geometry_config()
        call test_physics_validation("inconsistent_geometry.conf", "GEOMETRY_INCONSISTENT")
        
        ! Test 4: Out of range frequency values
        print *, ""
        print *, "Test 4: Out of range frequency values..."
        call create_invalid_frequency_config()
        call test_physics_validation("invalid_frequency.conf", "FREQUENCY_RANGE")
        
        ! Test 5: Invalid array dimensions
        print *, ""
        print *, "Test 5: Invalid array dimensions..."
        call create_invalid_array_config()
        call test_parameter_validation("invalid_array.conf", "ARRAY_BOUNDS")
        
        ! Test 6: Invalid mode numbers
        print *, ""
        print *, "Test 6: Invalid mode numbers..."
        call create_invalid_modes_config()
        call test_parameter_validation("invalid_modes.conf", "MODE_INVALID")
        
        ! Test 7: Temperature/density range validation
        print *, ""
        print *, "Test 7: Temperature/density ranges..."
        call create_invalid_profiles_config()
        call test_physics_validation("invalid_profiles.conf", "PROFILE_RANGE")
        
        ! Test 8: Numerical precision validation
        print *, ""
        print *, "Test 8: Numerical precision validation..."
        call create_precision_config()
        call test_parameter_validation("precision.conf", "PRECISION_INVALID")
        
        ! Test 9: File path validation
        print *, ""
        print *, "Test 9: File path validation..."
        call create_invalid_path_config()
        call test_parameter_validation("invalid_path.conf", "PATH_INVALID")
        
        ! Test 10: Complex validation across sections
        print *, ""
        print *, "Test 10: Cross-section validation..."
        call create_cross_section_config()
        call test_cross_validation("cross_section.conf", "CROSS_SECTION")
        
        ! Clean up test files
        call cleanup_test_files()
        
    end subroutine test_parameter_validation_functionality
    
    !> Test physics constraint validation
    subroutine test_physics_validation(filename, error_type)
        character(len=*), intent(in) :: filename, error_type
        integer :: validation_error
        
        ! Read the file successfully
        call settings_read_namelist(filename, settings, ierr)
        
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not read config file for", trim(error_type)
            test_status = test_status + 1
            return
        end if
        
        ! Now validate physics constraints - this should fail
        call settings_validate_physics_constraints(settings, validation_error)
        
        if (validation_error == KILCA_SUCCESS) then
            print *, "FAIL: Should have failed physics validation for", trim(error_type)
            test_status = test_status + 1
        else
            print *, "PASS: Correctly detected physics error for", trim(error_type), "error:", validation_error
        end if
        
    end subroutine test_physics_validation
    
    !> Test parameter range validation
    subroutine test_parameter_validation(filename, error_type)
        character(len=*), intent(in) :: filename, error_type
        integer :: validation_error
        
        ! Read the file successfully
        call settings_read_namelist(filename, settings, ierr)
        
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not read config file for", trim(error_type)
            test_status = test_status + 1
            return
        end if
        
        ! Now validate parameters - this should fail
        call settings_validate_parameters(settings, validation_error)
        
        if (validation_error == KILCA_SUCCESS) then
            print *, "FAIL: Should have failed parameter validation for", trim(error_type)
            test_status = test_status + 1
        else
            print *, "PASS: Correctly detected parameter error for", trim(error_type), "error:", validation_error
        end if
        
    end subroutine test_parameter_validation
    
    !> Test cross-section validation
    subroutine test_cross_validation(filename, error_type)
        character(len=*), intent(in) :: filename, error_type
        
        ! This functionality doesn't exist yet - should fail
        print *, "FAIL: Cross-section validation not implemented yet"
        test_status = test_status + 1
        
    end subroutine test_cross_validation
    
    ! ===== Test file creation routines =====
    
    !> Create config with negative radius values
    subroutine create_negative_radius_config()
        integer :: unit
        
        open(newunit=unit, file="negative_radius.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = -90.0"      ! Negative antenna radius
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  flab = (1.0e6, 0.0)"
        write(unit, '(a)') "  dma = 2"
        write(unit, '(a)') "  flag_debug_ant = 0"
        write(unit, '(a)') "  flag_eigmode = 1"
        write(unit, '(a)') "  modes = 1, 1, 2, -1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
        write(unit, '(a)') "  rp = -65.0"      ! Negative plasma radius
        write(unit, '(a)') "  B0 = 25000.0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "/"
        close(unit)
    end subroutine create_negative_radius_config
    
    !> Create config with invalid magnetic field
    subroutine create_invalid_bfield_config()
        integer :: unit
        
        open(newunit=unit, file="invalid_bfield.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  flab = (1.0e6, 0.0)"
        write(unit, '(a)') "  dma = 2"
        write(unit, '(a)') "  flag_debug_ant = 0"
        write(unit, '(a)') "  flag_eigmode = 1"
        write(unit, '(a)') "  modes = 1, 1, 2, -1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
        write(unit, '(a)') "  rp = 65.0"
        write(unit, '(a)') "  B0 = -25000.0"   ! Negative magnetic field
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "/"
        close(unit)
    end subroutine create_invalid_bfield_config
    
    !> Create config with inconsistent geometry
    subroutine create_inconsistent_geometry_config()
        integer :: unit
        
        open(newunit=unit, file="inconsistent_geometry.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  flab = (1.0e6, 0.0)"
        write(unit, '(a)') "  dma = 2"
        write(unit, '(a)') "  flag_debug_ant = 0"
        write(unit, '(a)') "  flag_eigmode = 1"
        write(unit, '(a)') "  modes = 1, 1, 2, -1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 50.0"     ! rtor < rp (physically impossible)
        write(unit, '(a)') "  rp = 100.0"
        write(unit, '(a)') "  B0 = 25000.0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "/"
        close(unit)
    end subroutine create_inconsistent_geometry_config
    
    !> Create config with invalid frequency values
    subroutine create_invalid_frequency_config()
        integer :: unit
        
        open(newunit=unit, file="invalid_frequency.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  flab = (-1.0e6, 0.0)"  ! Negative frequency
        write(unit, '(a)') "  dma = 2"
        write(unit, '(a)') "  flag_debug_ant = 0"
        write(unit, '(a)') "  flag_eigmode = 1"
        write(unit, '(a)') "  modes = 1, 1, 2, -1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
        write(unit, '(a)') "  rp = 65.0"
        write(unit, '(a)') "  B0 = 25000.0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "  rfmin = 1.0e10"    ! rfmin > rfmax
        write(unit, '(a)') "  rfmax = 1.0e8"
        write(unit, '(a)') "/"
        close(unit)
    end subroutine create_invalid_frequency_config
    
    !> Create config with invalid array dimensions
    subroutine create_invalid_array_config()
        integer :: unit
        
        open(newunit=unit, file="invalid_array.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  flab = (1.0e6, 0.0)"
        write(unit, '(a)') "  dma = 0"           ! Zero modes - invalid
        write(unit, '(a)') "  flag_debug_ant = 0"
        write(unit, '(a)') "  flag_eigmode = 1"
        write(unit, '(a)') "  modes = 1, 1, 2, -1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
        write(unit, '(a)') "  rp = 65.0"
        write(unit, '(a)') "  B0 = 25000.0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 1"
        write(unit, '(a)') "  num_quants = -5"   ! Negative array size
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "/"
        close(unit)
    end subroutine create_invalid_array_config
    
    !> Create config with invalid mode numbers
    subroutine create_invalid_modes_config()
        integer :: unit
        
        open(newunit=unit, file="invalid_modes.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  flab = (1.0e6, 0.0)"
        write(unit, '(a)') "  dma = 3"
        write(unit, '(a)') "  flag_debug_ant = 0"
        write(unit, '(a)') "  flag_eigmode = 1"
        write(unit, '(a)') "  modes = 0, 0, 2, -1, 3, -2"  ! m=0, n=0 is invalid mode
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
        write(unit, '(a)') "  rp = 65.0"
        write(unit, '(a)') "  B0 = 25000.0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "/"
        close(unit)
    end subroutine create_invalid_modes_config
    
    !> Create config with invalid profile ranges
    subroutine create_invalid_profiles_config()
        integer :: unit
        
        open(newunit=unit, file="invalid_profiles.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  flab = (1.0e6, 0.0)"
        write(unit, '(a)') "  dma = 2"
        write(unit, '(a)') "  flag_debug_ant = 0"
        write(unit, '(a)') "  flag_eigmode = 1"
        write(unit, '(a)') "  modes = 1, 1, 2, -1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
        write(unit, '(a)') "  rp = 65.0"
        write(unit, '(a)') "  B0 = 25000.0"
        write(unit, '(a)') "  mass = -2.0, 1.0, 4.0"        ! Negative mass
        write(unit, '(a)') "  charge = 1.0, -1.0, 5.0"      ! Invalid charge values
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "/"
        close(unit)
    end subroutine create_invalid_profiles_config
    
    !> Create config with precision issues
    subroutine create_precision_config()
        integer :: unit
        
        open(newunit=unit, file="precision.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  flab = (1.0e6, 0.0)"
        write(unit, '(a)') "  dma = 2"
        write(unit, '(a)') "  flag_debug_ant = 0"
        write(unit, '(a)') "  flag_eigmode = 1"
        write(unit, '(a)') "  modes = 1, 1, 2, -1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
        write(unit, '(a)') "  rp = 65.0"
        write(unit, '(a)') "  B0 = 25000.0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "  eps_res = 1.0e-20"     ! Precision too high
        write(unit, '(a)') "  eps_abs = 1.0"         ! Precision too low
        write(unit, '(a)') "/"
        close(unit)
    end subroutine create_precision_config
    
    !> Create config with invalid paths
    subroutine create_invalid_path_config()
        integer :: unit
        
        open(newunit=unit, file="invalid_path.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  flab = (1.0e6, 0.0)"
        write(unit, '(a)') "  dma = 2"
        write(unit, '(a)') "  flag_debug_ant = 0"
        write(unit, '(a)') "  flag_eigmode = 1"
        write(unit, '(a)') "  modes = 1, 1, 2, -1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
        write(unit, '(a)') "  rp = 65.0"
        write(unit, '(a)') "  B0 = 25000.0"
        write(unit, '(a)') "  path2profiles = '/invalid/path/that/does/not/exist/'"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "  fname = '/tmp/nonexistent/output.dat'"
        write(unit, '(a)') "/"
        close(unit)
    end subroutine create_invalid_path_config
    
    !> Create config with cross-section validation issues
    subroutine create_cross_section_config()
        integer :: unit
        
        open(newunit=unit, file="cross_section.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 30.0"        ! Antenna inside plasma
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  flab = (1.0e6, 0.0)"
        write(unit, '(a)') "  dma = 2"
        write(unit, '(a)') "  flag_debug_ant = 0"
        write(unit, '(a)') "  flag_eigmode = 1"
        write(unit, '(a)') "  modes = 1, 1, 2, -1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
        write(unit, '(a)') "  rp = 65.0"        ! rp > ra (antenna inside plasma)
        write(unit, '(a)') "  B0 = 25000.0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "/"
        close(unit)
    end subroutine create_cross_section_config
    
    !> Clean up test files
    subroutine cleanup_test_files()
        call system("rm -f negative_radius.conf")
        call system("rm -f invalid_bfield.conf")
        call system("rm -f inconsistent_geometry.conf")
        call system("rm -f invalid_frequency.conf")
        call system("rm -f invalid_array.conf")
        call system("rm -f invalid_modes.conf")
        call system("rm -f invalid_profiles.conf")
        call system("rm -f precision.conf")
        call system("rm -f invalid_path.conf")
        call system("rm -f cross_section.conf")
        print *, ""
        print *, "Test files cleaned up"
    end subroutine cleanup_test_files

end program test_parameter_validation_errors