program test_output_settings_cpp_vars
    use kilca_types_m, only: dp, KILCA_SUCCESS
    use kilca_settings_m, only: output_sett_t, output_sett_initialize_defaults, &
                                output_sett_validate, output_settings_set_flag_quants, &
                                output_settings_get_flag_quants
    implicit none
    
    type(output_sett_t) :: os
    integer :: ierr
    integer :: test_status = 0
    logical :: is_valid
    character(len=1024) :: error_msg
    
    print *, "=========================================="
    print *, "Testing Output Settings Complete C++ Variables [RED PHASE - SHOULD FAIL]"
    print *, "=========================================="
    
    call test_all_cpp_variables_exist()
    call test_dynamic_array_operations()
    call test_array_size_validation()
    
    if (test_status == 0) then
        print *, ""
        print *, "All output settings tests PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "Output settings tests FAILED:", test_status, "failure(s) (expected in RED phase)"
        stop 1
    end if

contains

    !> Test that all C++ variables from output_sett class exist and work
    subroutine test_all_cpp_variables_exist()
        print *, "Testing all C++ variables exist..."
        
        ! Initialize with defaults
        call output_sett_initialize_defaults(os, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not initialize output settings"
            test_status = test_status + 1
            return
        end if
        
        ! Test all C++ variables are accessible and settable
        os%flag_background = 2
        os%flag_emfield = 2
        os%flag_additional = 2
        os%flag_dispersion = 0
        os%flag_debug = 0
        os%num_quants = 8
        
        ! This should work - dynamic array already exists in structure
        if (allocated(os%flag_quants)) deallocate(os%flag_quants)
        allocate(os%flag_quants(8))
        os%flag_quants = [1, 1, 1, 1, 1, 1, 1, 0]
        
        ! Validate the structure
        call output_sett_validate(os, is_valid, error_msg, ierr)
        if (.not. is_valid .or. ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Output settings validation failed:", trim(error_msg)
            test_status = test_status + 1
        else
            print *, "PASS: All output C++ variables exist and validation works"
        end if
    end subroutine test_all_cpp_variables_exist
    
    !> Test dynamic array operations that should fail initially
    subroutine test_dynamic_array_operations()
        print *, "Testing dynamic array operations..."
        
        call output_sett_initialize_defaults(os, ierr)
        os%num_quants = 5
        
        ! This should FAIL initially - advanced array management not implemented
        call output_settings_set_flag_quants(os, [1, 0, 1, 0, 1], ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not set flag_quants array, error:", ierr
            test_status = test_status + 1
        else
            print *, "PASS: Advanced array management works correctly"
        end if
        
        ! Test getting array back
        block
            integer, allocatable :: quants_out(:)
            call output_settings_get_flag_quants(os, quants_out, ierr)
            if (ierr /= KILCA_SUCCESS) then
                print *, "FAIL: Could not get flag_quants array, error:", ierr
                test_status = test_status + 1
            else if (.not. allocated(quants_out)) then
                print *, "FAIL: Returned array not allocated"
                test_status = test_status + 1
            else if (size(quants_out) /= 5) then
                print *, "FAIL: Returned array wrong size:", size(quants_out)
                test_status = test_status + 1
            else if (any(quants_out /= [1, 0, 1, 0, 1])) then
                print *, "FAIL: Returned array has wrong values"
                test_status = test_status + 1
            else
                print *, "PASS: Array get operation works correctly"
            end if
        end block
    end subroutine test_dynamic_array_operations
    
    !> Test array size validation against num_quants
    subroutine test_array_size_validation()
        print *, "Testing array size validation..."
        
        call output_sett_initialize_defaults(os, ierr)
        
        ! Test inconsistent num_quants vs array size
        os%num_quants = 5
        if (allocated(os%flag_quants)) deallocate(os%flag_quants)
        allocate(os%flag_quants(3))  ! Wrong size
        os%flag_quants = [1, 0, 1]
        
        call output_sett_validate(os, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Array size mismatch not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Array size mismatch correctly detected:", trim(error_msg)
        end if
        
        ! Test correct size
        os%num_quants = 3
        call output_sett_validate(os, is_valid, error_msg, ierr)
        if (.not. is_valid) then
            print *, "FAIL: Valid array size rejected:", trim(error_msg)
            test_status = test_status + 1
        else
            print *, "PASS: Valid array size accepted"
        end if
        
        ! Test enhanced validation for array values
        os%flag_quants = [1, 3, 0]  ! Invalid value 3
        call output_sett_validate(os, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Invalid array values not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Invalid array values correctly detected:", trim(error_msg)
        end if
    end subroutine test_array_size_validation

end program test_output_settings_cpp_vars