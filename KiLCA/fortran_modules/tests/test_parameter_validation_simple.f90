program test_parameter_validation_simple
    use kilca_types_m, only: dp, KILCA_SUCCESS, KILCA_ERROR_INVALID_PARAMETER
    use kilca_settings_m, only: settings_t, settings_read_namelist, &
                                settings_validate_parameters, &
                                settings_validate_physics_constraints
    implicit none
    
    type(settings_t) :: settings
    integer :: ierr
    integer :: test_status = 0
    
    print *, "=========================================="
    print *, "Testing Parameter Validation Functions [GREEN PHASE]"
    print *, "=========================================="
    
    call test_validation_functions_exist()
    
    if (test_status == 0) then
        print *, ""
        print *, "Parameter validation functions test PASSED"
        stop 0
    else
        print *, ""
        print *, "Parameter validation functions test FAILED:", test_status, "failure(s)"
        stop 1
    end if

contains

    !> Test that validation functions exist and can be called
    subroutine test_validation_functions_exist()
        integer :: validation_error
        
        print *, "Testing validation function availability..."
        
        ! Create a minimal valid config
        call create_minimal_valid_config()
        
        ! Read the config
        call settings_read_namelist("minimal_valid.conf", settings, ierr)
        
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not read minimal config, error:", ierr
            test_status = test_status + 1
            call cleanup_test_files()
            return
        end if
        
        print *, "SUCCESS: Config read successfully"
        
        ! Test parameter validation function - should pass for valid config
        call settings_validate_parameters(settings, validation_error)
        
        if (validation_error == KILCA_SUCCESS) then
            print *, "PASS: Parameter validation function works"
        else
            print *, "FAIL: Parameter validation failed for valid config, error:", validation_error
            test_status = test_status + 1
        end if
        
        ! Test physics constraint validation function - should pass for valid config
        call settings_validate_physics_constraints(settings, validation_error)
        
        if (validation_error == KILCA_SUCCESS) then
            print *, "PASS: Physics constraint validation function works"
        else
            print *, "FAIL: Physics constraint validation failed for valid config, error:", validation_error
            test_status = test_status + 1
        end if
        
        ! Clean up
        call cleanup_test_files()
        
    end subroutine test_validation_functions_exist
    
    !> Create minimal valid configuration file
    subroutine create_minimal_valid_config()
        integer :: unit
        
        open(newunit=unit, file="minimal_valid.conf", status="replace")
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
        write(unit, '(a)') "/"
        close(unit)
    end subroutine create_minimal_valid_config
    
    !> Clean up test files
    subroutine cleanup_test_files()
        call system("rm -f minimal_valid.conf")
        print *, "Test files cleaned up"
    end subroutine cleanup_test_files

end program test_parameter_validation_simple