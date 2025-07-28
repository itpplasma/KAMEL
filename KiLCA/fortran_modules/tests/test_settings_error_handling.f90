program test_settings_error_handling
    use iso_fortran_env, only: real64, int32, int64
    use iso_c_binding
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    ! Test variables
    integer :: test_status = 0
    
    ! Test 1: Test error message generation
    call test_error_message_generation()
    
    ! Test 2: Test error code mapping
    call test_error_code_mapping()
    
    ! Test 3: Test error propagation
    call test_error_propagation()
    
    ! Test 4: Test error recovery
    call test_error_recovery()
    
    ! Test 5: Test validation error details
    call test_validation_error_details()
    
    ! Test 6: Test file operation errors
    call test_file_operation_errors()
    
    ! Test 7: Test memory allocation errors
    call test_memory_allocation_errors()
    
    ! Test 8: Test error logging functionality
    call test_error_logging()
    
    if (test_status == 0) then
        print *, "All settings error handling tests PASSED!"
    else
        print *, "Some tests FAILED. Status:", test_status
        stop 1
    end if
    
contains

    subroutine test_error_message_generation()
        character(len=1024) :: error_msg
        integer :: error_code
        character(len=256) :: context
        
        print *, "Testing error message generation..."
        
        ! Test basic error message generation
        error_code = KILCA_ERROR_INVALID_INPUT
        context = "antenna validation"
        call settings_format_error_message(error_code, context, "ra parameter", error_msg)
        
        if (len_trim(error_msg) == 0) then
            print *, "FAIL: Error message should not be empty"
            test_status = test_status + 1
        end if
        
        if (index(error_msg, "antenna validation") == 0) then
            print *, "FAIL: Error message should contain context"
            test_status = test_status + 1
        end if
        
        if (index(error_msg, "ra parameter") == 0) then
            print *, "FAIL: Error message should contain specific parameter"
            test_status = test_status + 1
        end if
        
        ! Test different error codes
        error_code = KILCA_ERROR_MEMORY
        call settings_format_error_message(error_code, "allocation", "antenna modes", error_msg)
        
        if (index(error_msg, "memory") == 0 .and. index(error_msg, "Memory") == 0) then
            print *, "FAIL: Memory error message should contain memory reference"
            test_status = test_status + 1
        end if
        
        print *, "test_error_message_generation completed"
    end subroutine test_error_message_generation
    
    subroutine test_error_code_mapping()
        character(len=256) :: error_name
        integer :: ierr
        
        print *, "Testing error code mapping..."
        
        ! Test mapping error codes to names
        call settings_get_error_name(KILCA_SUCCESS, error_name, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: get_error_name should succeed for valid code"
            test_status = test_status + 1
        end if
        
        if (index(error_name, "SUCCESS") == 0) then
            print *, "FAIL: Success code should map to SUCCESS name"
            test_status = test_status + 1
        end if
        
        call settings_get_error_name(KILCA_ERROR_INVALID_INPUT, error_name, ierr)
        if (index(error_name, "INVALID_INPUT") == 0) then
            print *, "FAIL: Invalid input code should map correctly"
            test_status = test_status + 1
        end if
        
        call settings_get_error_name(KILCA_ERROR_MEMORY, error_name, ierr)
        if (index(error_name, "MEMORY") == 0) then
            print *, "FAIL: Memory error code should map correctly"
            test_status = test_status + 1
        end if
        
        ! Test invalid error code
        call settings_get_error_name(-999, error_name, ierr)
        if (ierr == KILCA_SUCCESS) then
            print *, "FAIL: Invalid error code should return error"
            test_status = test_status + 1
        end if
        
        print *, "test_error_code_mapping completed"
    end subroutine test_error_code_mapping
    
    subroutine test_error_propagation()
        type(settings_t), pointer :: sd
        integer :: ierr
        character(len=1024) :: error_msg
        
        print *, "Testing error propagation..."
        
        ! Test that errors propagate correctly through the call stack
        call settings_create(sd, "", ierr)  ! Empty path should cause error
        if (ierr == KILCA_SUCCESS) then
            print *, "FAIL: Empty path should cause error"
            test_status = test_status + 1
        end if
        
        ! Test validation error propagation
        if (associated(sd)) then
            sd%antenna_settings%ra = -1.0_dp  ! Invalid value
            call settings_validate_with_context(sd, "test context", error_msg, ierr)
            if (ierr == KILCA_SUCCESS) then
                print *, "FAIL: Invalid antenna settings should propagate error"
                test_status = test_status + 1
            end if
            
            if (index(error_msg, "test context") == 0) then
                print *, "FAIL: Error context should be preserved"
                test_status = test_status + 1
            end if
            
            call settings_destroy(sd, ierr)
        end if
        
        print *, "test_error_propagation completed"
    end subroutine test_error_propagation
    
    subroutine test_error_recovery()
        type(antenna_t) :: ant
        integer :: ierr
        logical :: recovered
        
        print *, "Testing error recovery..."
        
        ! Set invalid antenna state
        ant%dma = 5
        ant%ra = -1.0_dp
        ! Don't allocate modes array
        
        ! Test recovery from invalid state
        call settings_attempt_recovery(ant, recovered, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Recovery attempt should not fail"
            test_status = test_status + 1
        end if
        
        if (recovered) then
            ! Check if recovery succeeded
            if (ant%ra < 0.0_dp) then
                print *, "FAIL: Recovery should fix invalid ra"
                test_status = test_status + 1
            end if
            
            if (ant%dma > 0 .and. .not. allocated(ant%modes)) then
                print *, "FAIL: Recovery should allocate modes if dma > 0"
                test_status = test_status + 1
            end if
        end if
        
        ! Clean up
        if (allocated(ant%modes)) deallocate(ant%modes)
        
        print *, "test_error_recovery completed"
    end subroutine test_error_recovery
    
    subroutine test_validation_error_details()
        type(back_sett_t) :: bs
        integer :: ierr
        character(len=2048) :: detailed_error
        
        print *, "Testing validation error details..."
        
        ! Set up multiple invalid values
        bs%rtor = -100.0_dp  ! Invalid
        bs%rp = 1000.0_dp    ! Invalid (larger than rtor)
        bs%B0 = -500.0_dp    ! Invalid
        bs%N = 4             ! Invalid (must be odd)
        
        ! Get detailed validation errors
        call settings_get_detailed_validation_errors(bs, detailed_error, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Detailed validation should succeed"
            test_status = test_status + 1
        end if
        
        ! Check that all errors are mentioned
        if (index(detailed_error, "rtor") == 0) then
            print *, "FAIL: Should mention rtor error"
            test_status = test_status + 1
        end if
        
        if (index(detailed_error, "rp") == 0) then
            print *, "FAIL: Should mention rp error"
            test_status = test_status + 1
        end if
        
        if (index(detailed_error, "B0") == 0) then
            print *, "FAIL: Should mention B0 error"
            test_status = test_status + 1
        end if
        
        if (index(detailed_error, "N") == 0) then
            print *, "FAIL: Should mention N error"
            test_status = test_status + 1
        end if
        
        print *, "test_validation_error_details completed"
    end subroutine test_validation_error_details
    
    subroutine test_file_operation_errors()
        type(antenna_t) :: ant
        integer :: ierr
        character(len=1024) :: error_msg
        
        print *, "Testing file operation errors..."
        
        ! Test reading from non-existent file
        call antenna_read_settings_with_error_context(ant, "/nonexistent/path/", error_msg, ierr)
        if (ierr /= KILCA_ERROR_FILE) then
            print *, "FAIL: Should return FILE error for non-existent path"
            test_status = test_status + 1
        end if
        
        if (index(error_msg, "file") == 0 .and. index(error_msg, "File") == 0) then
            print *, "FAIL: Error message should mention file issue"
            test_status = test_status + 1
        end if
        
        if (index(error_msg, "/nonexistent/path/") == 0) then
            print *, "FAIL: Error message should include attempted path"
            test_status = test_status + 1
        end if
        
        print *, "test_file_operation_errors completed"
    end subroutine test_file_operation_errors
    
    subroutine test_memory_allocation_errors()
        type(settings_t), pointer :: sd
        integer :: ierr
        character(len=1024) :: error_msg
        
        print *, "Testing memory allocation errors..."
        
        ! Test handling of allocation failures (simulated)
        call settings_test_allocation_failure(sd, error_msg, ierr)
        if (ierr /= KILCA_ERROR_MEMORY) then
            print *, "FAIL: Should return MEMORY error on allocation failure"
            test_status = test_status + 1
        end if
        
        if (index(error_msg, "memory") == 0 .and. index(error_msg, "Memory") == 0) then
            print *, "FAIL: Error message should mention memory"
            test_status = test_status + 1
        end if
        
        if (associated(sd)) then
            print *, "FAIL: Pointer should be null on allocation failure"
            test_status = test_status + 1
        end if
        
        print *, "test_memory_allocation_errors completed"
    end subroutine test_memory_allocation_errors
    
    subroutine test_error_logging()
        integer :: ierr
        character(len=256) :: log_file
        logical :: logged
        
        print *, "Testing error logging functionality..."
        
        log_file = "test_error.log"
        
        ! Test error logging
        call settings_log_error(KILCA_ERROR_INVALID_INPUT, "Test error logging", "ra parameter", log_file, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Error logging should succeed"
            test_status = test_status + 1
        end if
        
        ! Check if log was created
        call settings_check_error_log_exists(log_file, logged, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Log check should succeed"
            test_status = test_status + 1
        end if
        
        if (.not. logged) then
            print *, "FAIL: Error log should have been created"
            test_status = test_status + 1
        end if
        
        ! Clean up log file
        call settings_clear_error_log(log_file, ierr)
        
        print *, "test_error_logging completed"
    end subroutine test_error_logging

end program test_settings_error_handling