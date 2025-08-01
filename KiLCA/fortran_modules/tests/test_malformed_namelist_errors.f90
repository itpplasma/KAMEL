program test_malformed_namelist_errors
    use kilca_types_m, only: dp, KILCA_SUCCESS, KILCA_ERROR_FORMAT, &
                             KILCA_ERROR_NAMELIST_READ, KILCA_ERROR_INVALID_PARAMETER, &
                             KILCA_ERROR_NAMELIST_SYNTAX, KILCA_ERROR_TYPE_MISMATCH, &
                             KILCA_ERROR_MISSING_SECTION, KILCA_ERROR_DUPLICATE_SECTION, &
                             KILCA_ERROR_UNKNOWN_PARAMETER, KILCA_ERROR_EMPTY_FILE, &
                             KILCA_ERROR_COMPLEX_FORMAT, KILCA_ERROR_ARRAY_SIZE
    use kilca_settings_m, only: settings_t, settings_read_namelist, &
                                settings_read_namelist_with_details, &
                                settings_read_namelist_all_errors
    implicit none
    
    type(settings_t) :: settings
    integer :: ierr
    integer :: test_status = 0
    character(len=*), parameter :: test_path = "./"
    
    print *, "==========================================="
    print *, "Testing Malformed Namelist Error Handling [RED PHASE - SHOULD FAIL]"
    print *, "==========================================="
    
    call test_malformed_namelist_functionality()
    
    if (test_status == 0) then
        print *, ""
        print *, "Malformed namelist error tests PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "Malformed namelist error tests FAILED:", test_status, "failure(s) (expected in RED phase)"
        stop 1
    end if

contains

    !> Test malformed namelist error handling functionality
    subroutine test_malformed_namelist_functionality()
        print *, "Testing malformed namelist error handling..."
        
        ! Test 1: Missing closing slash
        print *, ""
        print *, "Test 1: Missing closing slash..."
        call create_missing_slash_namelist()
        call test_error_handling("missing_slash.conf", "NAMELIST_SYNTAX", KILCA_ERROR_NAMELIST_SYNTAX)
        
        ! Test 2: Wrong data type
        print *, ""
        print *, "Test 2: Wrong data type..."
        call create_wrong_type_namelist()
        call test_error_handling("wrong_type.conf", "TYPE_MISMATCH", KILCA_ERROR_TYPE_MISMATCH)
        
        ! Test 3: Missing required section
        print *, ""
        print *, "Test 3: Missing required section..."
        call create_missing_section_namelist()
        call test_error_handling("missing_section.conf", "MISSING_SECTION", KILCA_ERROR_MISSING_SECTION)
        
        ! Test 4: Duplicate section
        print *, ""
        print *, "Test 4: Duplicate section..."
        call create_duplicate_section_namelist()
        call test_error_handling("duplicate_section.conf", "DUPLICATE_SECTION", KILCA_ERROR_DUPLICATE_SECTION)
        
        ! Test 5: Invalid parameter names
        print *, ""
        print *, "Test 5: Invalid parameter names..."
        call create_invalid_parameter_namelist()
        call test_error_handling("invalid_param.conf", "UNKNOWN_PARAMETER", KILCA_ERROR_UNKNOWN_PARAMETER)
        
        ! Test 6: Empty file
        print *, ""
        print *, "Test 6: Empty file..."
        call create_empty_namelist()
        call test_error_handling("empty.conf", "EMPTY_FILE", KILCA_ERROR_EMPTY_FILE)
        
        ! Test 7: Corrupted complex number format
        print *, ""
        print *, "Test 7: Corrupted complex number format..."
        call create_bad_complex_namelist()
        call test_error_handling("bad_complex.conf", "COMPLEX_FORMAT", KILCA_ERROR_COMPLEX_FORMAT)
        
        ! Test 8: Array size mismatch
        print *, ""
        print *, "Test 8: Array size mismatch..."
        call create_array_mismatch_namelist()
        call test_error_handling("array_mismatch.conf", "ARRAY_SIZE", KILCA_ERROR_ARRAY_SIZE)
        
        ! Test 9: Syntax error with line number tracking
        print *, ""
        print *, "Test 9: Syntax error with line tracking..."
        call create_syntax_error_namelist()
        call test_detailed_error_reporting("syntax_error.conf")
        
        ! Test 10: Multiple errors in one file
        print *, ""
        print *, "Test 10: Multiple errors in one file..."
        call create_multiple_errors_namelist()
        call test_multiple_error_reporting("multiple_errors.conf")
        
        ! Clean up test files
        call cleanup_test_files()
        
    end subroutine test_malformed_namelist_functionality
    
    !> Test basic error handling and error code
    subroutine test_error_handling(filename, error_type, expected_code)
        character(len=*), intent(in) :: filename, error_type
        integer, intent(in) :: expected_code
        integer :: actual_code
        
        call settings_read_namelist(filename, settings, actual_code)
        
        if (actual_code == KILCA_SUCCESS) then
            print *, "FAIL: Should have failed for", trim(error_type), "but succeeded"
            test_status = test_status + 1
        else if (actual_code /= expected_code) then
            print *, "FAIL: Wrong error code for", trim(error_type)
            print *, "      Expected:", expected_code, "Got:", actual_code
            test_status = test_status + 1
        else
            print *, "PASS: Correct error code for", trim(error_type)
        end if
        
    end subroutine test_error_handling
    
    !> Test detailed error reporting with line numbers
    subroutine test_detailed_error_reporting(filename)
        character(len=*), intent(in) :: filename
        character(len=1024) :: error_msg
        integer :: error_line
        
        ! This functionality doesn't exist yet - should fail
        ! call settings_read_namelist_with_details(filename, settings, ierr, error_msg, error_line)
        
        print *, "FAIL: Detailed error reporting not implemented yet"
        test_status = test_status + 1
        
    end subroutine test_detailed_error_reporting
    
    !> Test multiple error reporting
    subroutine test_multiple_error_reporting(filename)
        character(len=*), intent(in) :: filename
        integer :: num_errors
        character(len=1024), dimension(10) :: error_list
        
        ! This functionality doesn't exist yet - should fail
        ! call settings_read_namelist_all_errors(filename, settings, ierr, num_errors, error_list)
        
        print *, "FAIL: Multiple error reporting not implemented yet"
        test_status = test_status + 1
        
    end subroutine test_multiple_error_reporting
    
    ! ===== Test file creation routines =====
    
    !> Create namelist with missing closing slash
    subroutine create_missing_slash_namelist()
        integer :: unit
        
        open(newunit=unit, file="missing_slash.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  ! Missing closing slash here"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
        write(unit, '(a)') "/"
        close(unit)
    end subroutine create_missing_slash_namelist
    
    !> Create namelist with wrong data types
    subroutine create_wrong_type_namelist()
        integer :: unit
        
        open(newunit=unit, file="wrong_type.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 'not_a_number'"  ! String instead of real
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  flab = (1.0e6, 0.0)"
        write(unit, '(a)') "  dma = 2.5"             ! Real instead of integer
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
    end subroutine create_wrong_type_namelist
    
    !> Create namelist with missing required section
    subroutine create_missing_section_namelist()
        integer :: unit
        
        open(newunit=unit, file="missing_section.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        ! Missing background section!
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "/"
        close(unit)
    end subroutine create_missing_section_namelist
    
    !> Create namelist with duplicate section
    subroutine create_duplicate_section_namelist()
        integer :: unit
        
        open(newunit=unit, file="duplicate_section.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&antenna"  ! Duplicate antenna section
        write(unit, '(a)') "  ra = 95.0"
        write(unit, '(a)') "  wa = 10.0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
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
    end subroutine create_duplicate_section_namelist
    
    !> Create namelist with invalid parameter names
    subroutine create_invalid_parameter_namelist()
        integer :: unit
        
        open(newunit=unit, file="invalid_param.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  unknown_param = 42"  ! Parameter doesn't exist
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
        write(unit, '(a)') "  rp = 65.0"
        write(unit, '(a)') "  B0 = 25000.0"
        write(unit, '(a)') "  bad_parameter = 'test'"  ! Unknown parameter
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
    end subroutine create_invalid_parameter_namelist
    
    !> Create empty namelist file
    subroutine create_empty_namelist()
        integer :: unit
        
        open(newunit=unit, file="empty.conf", status="replace")
        ! Write nothing
        close(unit)
    end subroutine create_empty_namelist
    
    !> Create namelist with bad complex number format
    subroutine create_bad_complex_namelist()
        integer :: unit
        
        open(newunit=unit, file="bad_complex.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  flab = (1.0e6 0.0)"  ! Missing comma
        write(unit, '(a)') "  dma = 2"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "  fstart = (1.0, 2.0, 3.0)"  ! Not a complex number
        write(unit, '(a)') "/"
        close(unit)
    end subroutine create_bad_complex_namelist
    
    !> Create namelist with array size mismatch
    subroutine create_array_mismatch_namelist()
        integer :: unit
        
        open(newunit=unit, file="array_mismatch.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  flab = (1.0e6, 0.0)"
        write(unit, '(a)') "  dma = 3"  ! Says 3 mode pairs
        write(unit, '(a)') "  flag_debug_ant = 0"
        write(unit, '(a)') "  flag_eigmode = 1"
        write(unit, '(a)') "  modes = 1, 1, 2, -1"  ! But only provides 2 pairs
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 1"
        write(unit, '(a)') "  num_quants = 5"  ! Says 5 quantities
        write(unit, '(a)') "  flag_quants = 1, 0, 1"  ! But only provides 3
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "/"
        close(unit)
    end subroutine create_array_mismatch_namelist
    
    !> Create namelist with syntax errors on specific lines
    subroutine create_syntax_error_namelist()
        integer :: unit
        
        open(newunit=unit, file="syntax_error.conf", status="replace")
        write(unit, '(a)') "! Configuration file with syntax errors"
        write(unit, '(a)') ""
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa == 5.0"  ! Double equals (line 5)
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  flab = (1.0e6, 0.0)"
        write(unit, '(a)') "  dma = 2"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor 170.0"  ! Missing equals (line 12)
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
    end subroutine create_syntax_error_namelist
    
    !> Create namelist with multiple errors
    subroutine create_multiple_errors_namelist()
        integer :: unit
        
        open(newunit=unit, file="multiple_errors.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = -90.0"  ! Negative radius (physics error)
        write(unit, '(a)') "  wa = 'wrong'"  ! Type error
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  unknown = 42"  ! Unknown parameter
        write(unit, '(a)') "  dma = 2"
        ! Missing closing slash
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 50.0"  ! rtor < rp (physics error)
        write(unit, '(a)') "  rp = 100.0"
        write(unit, '(a)') "  B0 = -1000.0"  ! Negative magnetic field
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        ! Missing output section
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  search_flag = 5"  ! Out of range value
        write(unit, '(a)') "/"
        close(unit)
    end subroutine create_multiple_errors_namelist
    
    !> Clean up test files
    subroutine cleanup_test_files()
        call system("rm -f missing_slash.conf")
        call system("rm -f wrong_type.conf")
        call system("rm -f missing_section.conf")
        call system("rm -f duplicate_section.conf")
        call system("rm -f invalid_param.conf")
        call system("rm -f empty.conf")
        call system("rm -f bad_complex.conf")
        call system("rm -f array_mismatch.conf")
        call system("rm -f syntax_error.conf")
        call system("rm -f multiple_errors.conf")
        print *, ""
        print *, "Test files cleaned up"
    end subroutine cleanup_test_files

end program test_malformed_namelist_errors