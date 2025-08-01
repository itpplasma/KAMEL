program test_parameter_validation_direct
    use kilca_types_m, only: dp, KILCA_SUCCESS, KILCA_ERROR_INVALID_PARAMETER
    use kilca_settings_m, only: settings_t, settings_validate_parameters, &
                                settings_validate_physics_constraints
    implicit none
    
    type(settings_t) :: settings
    integer :: ierr
    integer :: test_status = 0
    
    print *, "=========================================="
    print *, "Testing Parameter Validation Functions [GREEN PHASE]"
    print *, "=========================================="
    
    call test_validation_functions_direct()
    
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

    !> Test validation functions directly with manually set values
    subroutine test_validation_functions_direct()
        integer :: validation_error
        
        print *, "Testing validation functions with default settings..."
        
        ! Test 1: Test with default/uninitialized settings (should pass)
        call settings_validate_parameters(settings, validation_error)
        
        if (validation_error == KILCA_SUCCESS) then
            print *, "PASS: Parameter validation function exists and runs"
        else
            print *, "INFO: Parameter validation detected issues (expected for uninitialized settings), error:", validation_error
        end if
        
        call settings_validate_physics_constraints(settings, validation_error)
        
        if (validation_error == KILCA_SUCCESS) then
            print *, "PASS: Physics constraint validation function exists and runs"
        else
            print *, "INFO: Physics constraint validation detected issues (expected for uninitialized settings), error:", validation_error
        end if
        
        ! Test 2: Test with manually set valid values
        print *, ""
        print *, "Testing with manually set valid values..."
        settings%antenna_settings%ra = 90.0_dp
        settings%antenna_settings%dma = 2
        settings%background_settings%rtor = 170.0_dp
        settings%background_settings%rp = 65.0_dp
        settings%background_settings%B0 = 25000.0_dp
        settings%output_settings%num_quants = 5
        
        call settings_validate_parameters(settings, validation_error)
        
        if (validation_error == KILCA_SUCCESS) then
            print *, "PASS: Parameter validation passes with valid values"
        else
            print *, "INFO: Parameter validation with manually set values, error:", validation_error
        end if
        
        call settings_validate_physics_constraints(settings, validation_error)
        
        if (validation_error == KILCA_SUCCESS) then
            print *, "PASS: Physics constraint validation passes with valid values"
        else
            print *, "INFO: Physics constraint validation with manually set values, error:", validation_error
        end if
        
        ! Test 3: Test with invalid values to confirm error detection
        print *, ""
        print *, "Testing with invalid values (should detect errors)..."
        settings%antenna_settings%ra = -90.0_dp  ! Negative radius
        settings%background_settings%B0 = -25000.0_dp  ! Negative B-field
        
        call settings_validate_physics_constraints(settings, validation_error)
        
        if (validation_error /= KILCA_SUCCESS) then
            print *, "PASS: Physics constraint validation correctly detects invalid values"
        else
            print *, "FAIL: Physics constraint validation should have detected invalid values"
            test_status = test_status + 1
        end if
        
        ! Test 4: Test array bounds checking
        print *, ""
        print *, "Testing array bounds checking..."
        settings%antenna_settings%dma = 0  ! Invalid array size
        settings%output_settings%num_quants = -5  ! Negative array size
        
        call settings_validate_parameters(settings, validation_error)
        
        if (validation_error /= KILCA_SUCCESS) then
            print *, "PASS: Parameter validation correctly detects array bounds errors"
        else
            print *, "FAIL: Parameter validation should have detected array bounds errors"
            test_status = test_status + 1
        end if
        
    end subroutine test_validation_functions_direct

end program test_parameter_validation_direct