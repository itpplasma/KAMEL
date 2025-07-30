!> Test suite for KiLCA Fortran main program
!! Tests the main program functionality including:
!! - Command line argument parsing
!! - Project path handling
!! - Core data initialization
!! - Mode-independent calculations
!! - Mode-dependent calculations (antenna and eigmode)
!! - Error handling and cleanup
program test_kilca_main
    use iso_fortran_env, only: real64, error_unit
    use kilca_types_m
    use kilca_settings_m
    use kilca_background_m
    use kilca_mode_m
    use kilca_core_m
    use kilca_main_utils_m
    implicit none
    
    logical :: all_tests_passed
    
    all_tests_passed = .true.
    
    write(*, '(A)') "Running KiLCA main program tests..."
    
    ! Test 1: Command line argument parsing
    if (.not. test_command_line_parsing()) then
        all_tests_passed = .false.
    end if
    
    ! Test 2: Project path validation
    if (.not. test_project_path_validation()) then
        all_tests_passed = .false.
    end if
    
    ! Test 3: Core data initialization
    if (.not. test_core_data_initialization()) then
        all_tests_passed = .false.
    end if
    
    ! Test 4: Mode-independent calculations
    if (.not. test_mode_independent_calculations()) then
        all_tests_passed = .false.
    end if
    
    ! Test 5: Mode-dependent antenna calculations
    if (.not. test_mode_dependent_antenna_calculations()) then
        all_tests_passed = .false.
    end if
    
    ! Test 6: Mode-dependent eigmode calculations
    if (.not. test_mode_dependent_eigmode_calculations()) then
        all_tests_passed = .false.
    end if
    
    ! Test 7: Error handling and cleanup
    if (.not. test_error_handling_and_cleanup()) then
        all_tests_passed = .false.
    end if
    
    ! Test 8: Memory management
    if (.not. test_memory_management()) then
        all_tests_passed = .false.
    end if
    
    if (all_tests_passed) then
        write(*, '(A)') "All KiLCA main program tests passed!"
    else
        write(error_unit, '(A)') "Some KiLCA main program tests failed!"
        stop 1
    end if

contains

    !> Test command line argument parsing functionality
    function test_command_line_parsing() result(passed)
        logical :: passed
        character(len=1024) :: test_path, parsed_path
        character(len=1024) :: argv_test(2)
        integer :: ierr
        
        passed = .true.
        write(*, '(A)') "  Testing command line argument parsing..."
        
        ! Test 1: No arguments (should use current directory)
        call kilca_main_parse_arguments(0, [character(len=1024) ::], parsed_path, ierr)
        if (ierr /= 0) then
            write(error_unit, '(A)') "    FAIL: No arguments should default to current directory"
            passed = .false.
        else
            write(*, '(A)') "    PASS: No arguments handled correctly"
        end if
        
        ! Test 2: Single argument (project path)
        test_path = "/test/project/path"
        call kilca_main_parse_arguments(1, [test_path], parsed_path, ierr)
        if (ierr /= 0 .or. trim(parsed_path) /= trim(test_path)) then
            write(error_unit, '(A)') "    FAIL: Single argument not parsed correctly"
            passed = .false.
        else
            write(*, '(A)') "    PASS: Single argument parsed correctly"
        end if
        
        ! Test 3: Multiple arguments (should use first)
        argv_test(1) = test_path
        argv_test(2) = "/another/path"
        call kilca_main_parse_arguments(2, argv_test, parsed_path, ierr)
        if (ierr /= 0 .or. trim(parsed_path) /= trim(test_path)) then
            write(error_unit, '(A)') "    FAIL: Multiple arguments not handled correctly"
            passed = .false.
        else
            write(*, '(A)') "    PASS: Multiple arguments handled correctly"
        end if
        
    end function test_command_line_parsing

    !> Test project path validation
    function test_project_path_validation() result(passed)
        logical :: passed
        character(len=1024) :: test_path
        integer :: ierr
        
        passed = .true.
        write(*, '(A)') "  Testing project path validation..."
        
        ! Test 1: Valid path format
        test_path = "/valid/project/path/"
        call kilca_main_validate_project_path(test_path, ierr)
        if (ierr /= 0) then
            write(error_unit, '(A)') "    FAIL: Valid path rejected"
            passed = .false.
        else
            write(*, '(A)') "    PASS: Valid path accepted"
        end if
        
        ! Test 2: Path without trailing slash (should be fixed)
        test_path = "/valid/project/path"
        call kilca_main_validate_project_path(test_path, ierr)
        if (ierr /= 0 .or. test_path(len_trim(test_path):len_trim(test_path)) /= '/') then
            write(error_unit, '(A)') "    FAIL: Path without trailing slash not fixed"
            passed = .false.
        else
            write(*, '(A)') "    PASS: Trailing slash added correctly"
        end if
        
        ! Test 3: Empty path (should fail)
        test_path = ""
        call kilca_main_validate_project_path(test_path, ierr)
        if (ierr == 0) then
            write(error_unit, '(A)') "    FAIL: Empty path should be rejected"
            passed = .false.
        else
            write(*, '(A)') "    PASS: Empty path rejected correctly"
        end if
        
    end function test_project_path_validation

    !> Test core data initialization
    function test_core_data_initialization() result(passed)
        logical :: passed
        type(core_data_t), pointer :: core => null()
        character(len=1024) :: test_path
        integer :: ierr
        
        passed = .true.
        write(*, '(A)') "  Testing core data initialization..."
        
        ! Test basic initialization
        test_path = "/test/project/"
        call kilca_main_initialize_core_data(core, test_path, ierr)
        
        if (ierr /= 0) then
            write(error_unit, '(A)') "    FAIL: Core data initialization failed"
            passed = .false.
        else
            write(*, '(A)') "    PASS: Core data initialized successfully"
            
            ! Test core data components are allocated
            if (.not. associated(core%sd)) then
                write(error_unit, '(A)') "    FAIL: Settings not allocated"
                passed = .false.
            else
                write(*, '(A)') "    PASS: Settings allocated"
            end if
            
            if (.not. associated(core%bp)) then
                write(error_unit, '(A)') "    FAIL: Background not allocated"
                passed = .false.
            else
                write(*, '(A)') "    PASS: Background allocated"
            end if
        end if
        
        ! Clean up
        call kilca_main_cleanup_core_data(core, ierr)
        
    end function test_core_data_initialization

    !> Test mode-independent calculations
    function test_mode_independent_calculations() result(passed)
        logical :: passed
        type(core_data_t), pointer :: core => null()
        character(len=1024) :: test_path
        integer :: ierr
        
        passed = .true.
        write(*, '(A)') "  Testing mode-independent calculations..."
        
        ! Initialize core data first
        test_path = "/test/project/"
        call kilca_main_initialize_core_data(core, test_path, ierr)
        
        if (ierr == 0) then
            ! Test mode-independent calculations
            call kilca_main_calc_mode_independent_data(core, ierr)
            
            if (ierr /= 0) then
                write(error_unit, '(A)') "    FAIL: Mode-independent calculations failed"
                passed = .false.
            else
                write(*, '(A)') "    PASS: Mode-independent calculations completed"
            end if
        else
            write(error_unit, '(A)') "    FAIL: Could not initialize core data for test"
            passed = .false.
        end if
        
        ! Clean up
        call kilca_main_cleanup_core_data(core, ierr)
        
    end function test_mode_independent_calculations

    !> Test mode-dependent antenna calculations
    function test_mode_dependent_antenna_calculations() result(passed)
        logical :: passed
        type(core_data_t), pointer :: core => null()
        character(len=1024) :: test_path
        integer :: ierr
        
        passed = .true.
        write(*, '(A)') "  Testing mode-dependent antenna calculations..."
        
        ! Initialize core data first
        test_path = "/test/project/"
        call kilca_main_initialize_core_data(core, test_path, ierr)
        
        if (ierr == 0) then
            ! Test mode-dependent antenna calculations
            call kilca_main_calc_mode_dependent_data_antenna(core, ierr)
            
            if (ierr /= 0) then
                write(error_unit, '(A)') "    FAIL: Mode-dependent antenna calculations failed"
                passed = .false.
            else
                write(*, '(A)') "    PASS: Mode-dependent antenna calculations completed"
            end if
        else
            write(error_unit, '(A)') "    FAIL: Could not initialize core data for test"
            passed = .false.
        end if
        
        ! Clean up
        call kilca_main_cleanup_core_data(core, ierr)
        
    end function test_mode_dependent_antenna_calculations

    !> Test mode-dependent eigmode calculations
    function test_mode_dependent_eigmode_calculations() result(passed)
        logical :: passed
        type(core_data_t), pointer :: core => null()
        character(len=1024) :: test_path
        integer :: ierr
        
        passed = .true.
        write(*, '(A)') "  Testing mode-dependent eigmode calculations..."
        
        ! Initialize core data first
        test_path = "/test/project/"
        call kilca_main_initialize_core_data(core, test_path, ierr)
        
        if (ierr == 0) then
            ! Test mode-dependent eigmode calculations
            call kilca_main_calc_mode_dependent_data_eigmode(core, ierr)
            
            if (ierr /= 0) then
                write(error_unit, '(A)') "    FAIL: Mode-dependent eigmode calculations failed"
                passed = .false.
            else
                write(*, '(A)') "    PASS: Mode-dependent eigmode calculations completed"
            end if
        else
            write(error_unit, '(A)') "    FAIL: Could not initialize core data for test"
            passed = .false.
        end if
        
        ! Clean up
        call kilca_main_cleanup_core_data(core, ierr)
        
    end function test_mode_dependent_eigmode_calculations

    !> Test error handling and cleanup
    function test_error_handling_and_cleanup() result(passed)
        logical :: passed
        type(core_data_t), pointer :: core => null()
        character(len=1024) :: invalid_path
        integer :: ierr
        
        passed = .true.
        write(*, '(A)') "  Testing error handling and cleanup..."
        
        ! Test 1: Invalid path handling
        invalid_path = "/nonexistent/invalid/path/"
        call kilca_main_initialize_core_data(core, invalid_path, ierr)
        
        if (ierr == 0) then
            write(error_unit, '(A)') "    FAIL: Invalid path should cause error"
            passed = .false.
        else
            write(*, '(A)') "    PASS: Invalid path handled correctly"
        end if
        
        ! Test 2: Cleanup after error
        call kilca_main_cleanup_core_data(core, ierr)
        if (ierr /= 0) then
            write(error_unit, '(A)') "    FAIL: Cleanup after error failed"
            passed = .false.
        else
            write(*, '(A)') "    PASS: Cleanup after error successful"
        end if
        
    end function test_error_handling_and_cleanup

    !> Test memory management
    function test_memory_management() result(passed)
        logical :: passed
        type(core_data_t), pointer :: core => null()
        character(len=1024) :: test_path
        integer :: ierr, i
        
        passed = .true.
        write(*, '(A)') "  Testing memory management..."
        
        ! Test multiple allocation/deallocation cycles
        test_path = "/test/project/"
        
        do i = 1, 3
            call kilca_main_initialize_core_data(core, test_path, ierr)
            if (ierr /= 0) then
                write(error_unit, '(A,I0)') "    FAIL: Initialization failed on cycle ", i
                passed = .false.
                exit
            end if
            
            call kilca_main_cleanup_core_data(core, ierr)
            if (ierr /= 0) then
                write(error_unit, '(A,I0)') "    FAIL: Cleanup failed on cycle ", i
                passed = .false.
                exit
            end if
        end do
        
        if (passed) then
            write(*, '(A)') "    PASS: Multiple allocation/deallocation cycles successful"
        end if
        
    end function test_memory_management

end program test_kilca_main