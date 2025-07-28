program test_settings_documentation
    use iso_fortran_env, only: real64, int32, int64
    use iso_c_binding
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    ! Test variables
    integer :: test_status = 0
    
    ! Test 1: Test that documentation files exist
    call test_documentation_files_exist()
    
    ! Test 2: Test that procedure reference is complete
    call test_procedure_reference_complete()
    
    ! Test 3: Test that usage examples are functional
    call test_usage_examples_functional()
    
    ! Test 4: Test that error code documentation exists
    call test_error_code_documentation()
    
    ! Test 5: Test that validation procedure docs exist
    call test_validation_documentation()
    
    ! Test 6: Test that initialization procedure docs exist
    call test_initialization_documentation()
    
    ! Test 7: Test that I/O procedure documentation exists
    call test_io_documentation()
    
    ! Test 8: Test that examples directory has proper structure
    call test_examples_directory_structure()
    
    if (test_status == 0) then
        print *, "All settings documentation tests PASSED!"
    else
        print *, "Some tests FAILED. Status:", test_status
        stop 1
    end if
    
contains

    subroutine test_documentation_files_exist()
        logical :: file_exists
        
        print *, "Testing documentation files exist..."
        
        ! Test PROCEDURE_REFERENCE.md exists
        inquire(file="PROCEDURE_REFERENCE.md", exist=file_exists)
        if (.not. file_exists) then
            print *, "FAIL: PROCEDURE_REFERENCE.md does not exist"
            test_status = test_status + 1
        end if
        
        ! Test SETTINGS_USAGE.md exists
        inquire(file="SETTINGS_USAGE.md", exist=file_exists)
        if (.not. file_exists) then
            print *, "FAIL: SETTINGS_USAGE.md does not exist"
            test_status = test_status + 1
        end if
        
        ! Test examples directory exists
        inquire(file="examples", exist=file_exists)
        if (.not. file_exists) then
            print *, "FAIL: examples directory does not exist"
            test_status = test_status + 1
        end if
        
        print *, "test_documentation_files_exist completed"
    end subroutine test_documentation_files_exist
    
    subroutine test_procedure_reference_complete()
        logical :: file_exists
        integer :: unit_num, ierr
        character(len=1024) :: line
        logical :: found_settings_create = .false.
        logical :: found_antenna_procedures = .false.
        logical :: found_error_handling = .false.
        
        print *, "Testing procedure reference completeness..."
        
        inquire(file="PROCEDURE_REFERENCE.md", exist=file_exists)
        if (.not. file_exists) then
            print *, "FAIL: Cannot test procedure reference - file missing"
            test_status = test_status + 1
            return
        end if
        
        ! Open and read the procedure reference file
        open(newunit=unit_num, file="PROCEDURE_REFERENCE.md", status='old', &
             action='read', iostat=ierr)
        if (ierr /= 0) then
            print *, "FAIL: Cannot open PROCEDURE_REFERENCE.md"
            test_status = test_status + 1
            return
        end if
        
        ! Search for key procedure documentation
        do
            read(unit_num, '(A)', iostat=ierr) line
            if (ierr /= 0) exit
            
            if (index(line, "settings_create") > 0) found_settings_create = .true.
            if (index(line, "antenna_") > 0) found_antenna_procedures = .true.
            if (index(line, "error") > 0 .or. index(line, "Error") > 0) found_error_handling = .true.
        end do
        close(unit_num)
        
        ! Check that key procedures are documented
        if (.not. found_settings_create) then
            print *, "FAIL: settings_create not documented in procedure reference"
            test_status = test_status + 1
        end if
        
        if (.not. found_antenna_procedures) then
            print *, "FAIL: antenna procedures not documented"
            test_status = test_status + 1
        end if
        
        if (.not. found_error_handling) then
            print *, "FAIL: error handling not documented"
            test_status = test_status + 1
        end if
        
        print *, "test_procedure_reference_complete completed"
    end subroutine test_procedure_reference_complete
    
    subroutine test_usage_examples_functional()
        logical :: file_exists
        integer :: unit_num, ierr
        character(len=1024) :: line
        logical :: found_basic_example = .false.
        logical :: found_advanced_example = .false.
        
        print *, "Testing usage examples functionality..."
        
        inquire(file="SETTINGS_USAGE.md", exist=file_exists)
        if (.not. file_exists) then
            print *, "FAIL: Cannot test usage examples - file missing"
            test_status = test_status + 1
            return
        end if
        
        ! Open and read the usage file
        open(newunit=unit_num, file="SETTINGS_USAGE.md", status='old', &
             action='read', iostat=ierr)
        if (ierr /= 0) then
            print *, "FAIL: Cannot open SETTINGS_USAGE.md"
            test_status = test_status + 1
            return
        end if
        
        ! Search for usage examples
        do
            read(unit_num, '(A)', iostat=ierr) line
            if (ierr /= 0) exit
            
            if (index(line, "Basic Usage") > 0) found_basic_example = .true.
            if (index(line, "Advanced") > 0) found_advanced_example = .true.
        end do
        close(unit_num)
        
        if (.not. found_basic_example) then
            print *, "FAIL: Basic usage example not found"
            test_status = test_status + 1
        end if
        
        if (.not. found_advanced_example) then
            print *, "FAIL: Advanced usage example not found"
            test_status = test_status + 1
        end if
        
        print *, "test_usage_examples_functional completed"
    end subroutine test_usage_examples_functional
    
    subroutine test_error_code_documentation()
        logical :: file_exists
        integer :: unit_num, ierr
        character(len=1024) :: line
        logical :: found_kilca_success = .false.
        logical :: found_kilca_error_memory = .false.
        logical :: found_error_handling_section = .false.
        
        print *, "Testing error code documentation..."
        
        inquire(file="PROCEDURE_REFERENCE.md", exist=file_exists)
        if (.not. file_exists) then
            print *, "FAIL: Cannot test error codes - file missing"
            test_status = test_status + 1
            return
        end if
        
        open(newunit=unit_num, file="PROCEDURE_REFERENCE.md", status='old', &
             action='read', iostat=ierr)
        if (ierr /= 0) then
            print *, "FAIL: Cannot open procedure reference"
            test_status = test_status + 1
            return
        end if
        
        do
            read(unit_num, '(A)', iostat=ierr) line
            if (ierr /= 0) exit
            
            if (index(line, "KILCA_SUCCESS") > 0) found_kilca_success = .true.
            if (index(line, "KILCA_ERROR_MEMORY") > 0) found_kilca_error_memory = .true.
            if (index(line, "Error Handling") > 0) found_error_handling_section = .true.
        end do
        close(unit_num)
        
        if (.not. found_kilca_success) then
            print *, "FAIL: KILCA_SUCCESS not documented"
            test_status = test_status + 1
        end if
        
        if (.not. found_kilca_error_memory) then
            print *, "FAIL: Error codes not properly documented"
            test_status = test_status + 1
        end if
        
        if (.not. found_error_handling_section) then
            print *, "FAIL: Error handling section not found"
            test_status = test_status + 1
        end if
        
        print *, "test_error_code_documentation completed"
    end subroutine test_error_code_documentation
    
    subroutine test_validation_documentation()
        logical :: file_exists
        integer :: unit_num, ierr
        character(len=1024) :: line
        logical :: found_validation_section = .false.
        logical :: found_settings_validate = .false.
        
        print *, "Testing validation documentation..."
        
        inquire(file="PROCEDURE_REFERENCE.md", exist=file_exists)
        if (.not. file_exists) then
            print *, "FAIL: Cannot test validation docs - file missing"
            test_status = test_status + 1
            return
        end if
        
        open(newunit=unit_num, file="PROCEDURE_REFERENCE.md", status='old', &
             action='read', iostat=ierr)
        if (ierr /= 0) then
            print *, "FAIL: Cannot open procedure reference"
            test_status = test_status + 1
            return
        end if
        
        do
            read(unit_num, '(A)', iostat=ierr) line
            if (ierr /= 0) exit
            
            if (index(line, "Validation") > 0) found_validation_section = .true.
            if (index(line, "settings_validate") > 0) found_settings_validate = .true.
        end do
        close(unit_num)
        
        if (.not. found_validation_section) then
            print *, "FAIL: Validation section not found"
            test_status = test_status + 1
        end if
        
        if (.not. found_settings_validate) then
            print *, "FAIL: settings_validate not documented"
            test_status = test_status + 1
        end if
        
        print *, "test_validation_documentation completed"
    end subroutine test_validation_documentation
    
    subroutine test_initialization_documentation()
        logical :: file_exists
        integer :: unit_num, ierr
        character(len=1024) :: line
        logical :: found_initialization_section = .false.
        logical :: found_initialize_defaults = .false.
        
        print *, "Testing initialization documentation..."
        
        inquire(file="PROCEDURE_REFERENCE.md", exist=file_exists)
        if (.not. file_exists) then
            print *, "FAIL: Cannot test initialization docs - file missing"
            test_status = test_status + 1
            return
        end if
        
        open(newunit=unit_num, file="PROCEDURE_REFERENCE.md", status='old', &
             action='read', iostat=ierr)
        if (ierr /= 0) then
            print *, "FAIL: Cannot open procedure reference"
            test_status = test_status + 1
            return
        end if
        
        do
            read(unit_num, '(A)', iostat=ierr) line
            if (ierr /= 0) exit
            
            if (index(line, "Initialization") > 0) found_initialization_section = .true.
            if (index(line, "initialize_defaults") > 0) found_initialize_defaults = .true.
        end do
        close(unit_num)
        
        if (.not. found_initialization_section) then
            print *, "FAIL: Initialization section not found"
            test_status = test_status + 1
        end if
        
        if (.not. found_initialize_defaults) then
            print *, "FAIL: initialize_defaults not documented"
            test_status = test_status + 1
        end if
        
        print *, "test_initialization_documentation completed"
    end subroutine test_initialization_documentation
    
    subroutine test_io_documentation()
        logical :: file_exists
        integer :: unit_num, ierr
        character(len=1024) :: line
        logical :: found_io_section = .false.
        logical :: found_print_settings = .false.
        
        print *, "Testing I/O documentation..."
        
        inquire(file="PROCEDURE_REFERENCE.md", exist=file_exists)
        if (.not. file_exists) then
            print *, "FAIL: Cannot test I/O docs - file missing"
            test_status = test_status + 1
            return
        end if
        
        open(newunit=unit_num, file="PROCEDURE_REFERENCE.md", status='old', &
             action='read', iostat=ierr)
        if (ierr /= 0) then
            print *, "FAIL: Cannot open procedure reference"
            test_status = test_status + 1
            return
        end if
        
        do
            read(unit_num, '(A)', iostat=ierr) line
            if (ierr /= 0) exit
            
            if (index(line, "Input/Output") > 0 .or. index(line, "I/O") > 0) found_io_section = .true.
            if (index(line, "print_settings") > 0) found_print_settings = .true.
        end do
        close(unit_num)
        
        if (.not. found_io_section) then
            print *, "FAIL: I/O section not found"
            test_status = test_status + 1
        end if
        
        if (.not. found_print_settings) then
            print *, "FAIL: print_settings not documented"
            test_status = test_status + 1
        end if
        
        print *, "test_io_documentation completed"
    end subroutine test_io_documentation
    
    subroutine test_examples_directory_structure()
        logical :: file_exists
        
        print *, "Testing examples directory structure..."
        
        ! Check for basic example
        inquire(file="examples/example_basic_usage.f90", exist=file_exists)
        if (.not. file_exists) then
            print *, "FAIL: example_basic_usage.f90 not found"
            test_status = test_status + 1
        end if
        
        ! Check for advanced example
        inquire(file="examples/example_advanced_usage.f90", exist=file_exists)
        if (.not. file_exists) then
            print *, "FAIL: example_advanced_usage.f90 not found"
            test_status = test_status + 1
        end if
        
        ! Check for Makefile
        inquire(file="examples/Makefile", exist=file_exists)
        if (.not. file_exists) then
            print *, "FAIL: examples/Makefile not found"
            test_status = test_status + 1
        end if
        
        print *, "test_examples_directory_structure completed"
    end subroutine test_examples_directory_structure

end program test_settings_documentation