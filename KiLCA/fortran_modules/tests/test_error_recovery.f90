program test_error_recovery
    use kilca_types_m, only: dp, KILCA_SUCCESS, KILCA_ERROR_INVALID_PARAMETER
    use kilca_settings_m, only: settings_t, settings_read_namelist_with_recovery, &
                                settings_get_warning_count, settings_get_warnings, &
                                settings_validate_with_warnings
    implicit none
    
    type(settings_t) :: settings
    integer :: ierr
    integer :: test_status = 0
    character(len=*), parameter :: test_path = "./"
    
    print *, "=========================================="
    print *, "Testing Error Recovery and Partial Loading [REFACTOR PHASE]"
    print *, "=========================================="
    
    call test_error_recovery_functionality()
    
    if (test_status == 0) then
        print *, ""
        print *, "Error recovery tests PASSED"
        stop 0
    else
        print *, ""
        print *, "Error recovery tests FAILED:", test_status, "failure(s)"
        stop 1
    end if

contains

    !> Test error recovery functionality
    subroutine test_error_recovery_functionality()
        print *, "Testing error recovery functionality..."
        
        ! Test 1: Partial loading with minor errors
        print *, ""
        print *, "Test 1: Partial loading with recoverable errors..."
        call create_config_with_minor_errors()
        call test_partial_loading("minor_errors.conf")
        
        ! Test 2: Warning system for questionable parameters
        print *, ""
        print *, "Test 2: Warning system for questionable values..."
        call create_config_with_warnings()
        call test_warning_system("questionable.conf")
        
        ! Test 3: Recovery from missing optional sections
        print *, ""
        print *, "Test 3: Recovery from missing optional sections..."
        call create_incomplete_config()
        call test_missing_sections("incomplete.conf")
        
        ! Test 4: Recovery with default values
        print *, ""
        print *, "Test 4: Recovery with default values..."
        call create_minimal_config()
        call test_default_recovery("minimal.conf")
        
        ! Test 5: Multiple warnings aggregation
        print *, ""
        print *, "Test 5: Multiple warnings aggregation..."
        call create_config_multiple_warnings()
        call test_multiple_warnings("multiple_warnings.conf")
        
        ! Clean up test files
        call cleanup_test_files()
        
    end subroutine test_error_recovery_functionality
    
    !> Test partial loading with recoverable errors
    subroutine test_partial_loading(filename)
        character(len=*), intent(in) :: filename
        integer :: warning_count
        character(len=256), dimension(10) :: warnings
        integer :: i
        
        ! Read with recovery enabled
        call settings_read_namelist_with_recovery(filename, settings, ierr, warning_count, warnings)
        
        if (ierr == KILCA_SUCCESS) then
            print *, "PASS: Partial loading succeeded despite errors"
            
            ! Check if warnings were generated
            if (warning_count > 0) then
                print *, "PASS: Generated", warning_count, "warnings"
                do i = 1, warning_count
                    print *, "  Warning:", trim(warnings(i))
                end do
            else
                print *, "FAIL: Should have generated warnings"
                test_status = test_status + 1
            end if
        else
            print *, "FAIL: Should have recovered from minor errors"
            test_status = test_status + 1
        end if
        
    end subroutine test_partial_loading
    
    !> Test warning system
    subroutine test_warning_system(filename)
        character(len=*), intent(in) :: filename
        integer :: warning_count
        character(len=256), dimension(10) :: warnings
        
        ! Read and validate with warnings
        call settings_read_namelist_with_recovery(filename, settings, ierr, warning_count, warnings)
        
        if (ierr == KILCA_SUCCESS) then
            ! Validate with warning generation
            call settings_validate_with_warnings(settings, ierr, warning_count, warnings)
            
            if (warning_count > 0) then
                print *, "PASS: Warning system detected questionable values"
                print *, "      Generated", warning_count, "warnings"
            else
                print *, "FAIL: Should have detected questionable values"
                test_status = test_status + 1
            end if
        else
            print *, "FAIL: Reading should have succeeded with warnings"
            test_status = test_status + 1
        end if
        
    end subroutine test_warning_system
    
    !> Test recovery from missing sections
    subroutine test_missing_sections(filename)
        character(len=*), intent(in) :: filename
        integer :: warning_count
        character(len=256), dimension(10) :: warnings
        
        call settings_read_namelist_with_recovery(filename, settings, ierr, warning_count, warnings)
        
        if (ierr == KILCA_SUCCESS) then
            print *, "PASS: Recovered from missing sections with defaults"
            
            ! Check that defaults were applied
            if (settings%output_settings%flag_background == 1) then
                print *, "PASS: Default values applied correctly"
            else
                print *, "FAIL: Default values not applied"
                test_status = test_status + 1
            end if
        else
            print *, "FAIL: Should have recovered with defaults"
            test_status = test_status + 1
        end if
        
    end subroutine test_missing_sections
    
    !> Test default recovery
    subroutine test_default_recovery(filename)
        character(len=*), intent(in) :: filename
        integer :: warning_count
        character(len=256), dimension(10) :: warnings
        
        call settings_read_namelist_with_recovery(filename, settings, ierr, warning_count, warnings)
        
        if (ierr == KILCA_SUCCESS) then
            print *, "PASS: Minimal config loaded with extensive defaults"
            
            ! Check that sensible defaults were applied
            if (settings%eigmode_settings%eps_res > 0.0_dp) then
                print *, "PASS: Sensible defaults applied"
            else
                print *, "FAIL: Defaults not properly applied"
                test_status = test_status + 1
            end if
        else
            print *, "FAIL: Should load minimal config with defaults"
            test_status = test_status + 1
        end if
        
    end subroutine test_default_recovery
    
    !> Test multiple warnings aggregation
    subroutine test_multiple_warnings(filename)
        character(len=*), intent(in) :: filename
        integer :: warning_count
        character(len=256), dimension(10) :: warnings
        
        call settings_read_namelist_with_recovery(filename, settings, ierr, warning_count, warnings)
        
        if (ierr == KILCA_SUCCESS .and. warning_count >= 3) then
            print *, "PASS: Multiple warnings correctly aggregated"
            print *, "      Total warnings:", warning_count
        else
            print *, "FAIL: Should have aggregated multiple warnings"
            test_status = test_status + 1
        end if
        
    end subroutine test_multiple_warnings
    
    ! ===== Test file creation routines =====
    
    !> Create config with minor recoverable errors
    subroutine create_config_with_minor_errors()
        integer :: unit
        
        open(newunit=unit, file="minor_errors.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "  wa = 5.0"
        write(unit, '(a)') "  I0 = 1.0e12"
        write(unit, '(a)') "  flab = (1.0e6, 0.0)"
        write(unit, '(a)') "  dma = 2"
        write(unit, '(a)') "  flag_debug_ant = 0"
        write(unit, '(a)') "  flag_eigmode = 1"
        write(unit, '(a)') "  modes = 1, 1, 2, -1"
        write(unit, '(a)') "  unknown_param = 42"  ! Unknown parameter - should warn but continue
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
        write(unit, '(a)') "  rp = 65.0"
        write(unit, '(a)') "  B0 = 25000.0"
        write(unit, '(a)') "  typo_param = 1.0"   ! Typo - should warn but continue
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
    end subroutine create_config_with_minor_errors
    
    !> Create config with questionable values
    subroutine create_config_with_warnings()
        integer :: unit
        
        open(newunit=unit, file="questionable.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 0.1"          ! Very small radius - warning
        write(unit, '(a)') "  wa = 500.0"        ! Very large width - warning
        write(unit, '(a)') "  I0 = 1.0e20"       ! Extremely high intensity - warning
        write(unit, '(a)') "  flab = (1.0e3, 0.0)"  ! Very low frequency - warning
        write(unit, '(a)') "  dma = 2"
        write(unit, '(a)') "  flag_debug_ant = 0"
        write(unit, '(a)') "  flag_eigmode = 1"
        write(unit, '(a)') "  modes = 1, 1, 2, -1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
        write(unit, '(a)') "  rp = 169.9"        ! rp very close to rtor - warning
        write(unit, '(a)') "  B0 = 100.0"        ! Very low B field - warning
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&output"
        write(unit, '(a)') "  flag_background = 1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&eigenmode"
        write(unit, '(a)') "  search_flag = 1"
        write(unit, '(a)') "  eps_res = 1.0e-15"  ! Extremely tight tolerance - warning
        write(unit, '(a)') "/"
        close(unit)
    end subroutine create_config_with_warnings
    
    !> Create incomplete config (missing sections)
    subroutine create_incomplete_config()
        integer :: unit
        
        open(newunit=unit, file="incomplete.conf", status="replace")
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
        ! Missing output and eigenmode sections
        close(unit)
    end subroutine create_incomplete_config
    
    !> Create minimal config
    subroutine create_minimal_config()
        integer :: unit
        
        open(newunit=unit, file="minimal.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 90.0"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 170.0"
        write(unit, '(a)') "/"
        ! Minimal required parameters only
        close(unit)
    end subroutine create_minimal_config
    
    !> Create config with multiple warnings
    subroutine create_config_multiple_warnings()
        integer :: unit
        
        open(newunit=unit, file="multiple_warnings.conf", status="replace")
        write(unit, '(a)') "&antenna"
        write(unit, '(a)') "  ra = 0.001"        ! Warning 1: extremely small
        write(unit, '(a)') "  wa = 1000.0"       ! Warning 2: extremely large
        write(unit, '(a)') "  I0 = 1.0e25"       ! Warning 3: unrealistic intensity
        write(unit, '(a)') "  flab = (10.0, 0.0)" ! Warning 4: extremely low frequency
        write(unit, '(a)') "  dma = 50"          ! Warning 5: unusually many modes
        write(unit, '(a)') "  flag_debug_ant = 0"
        write(unit, '(a)') "  flag_eigmode = 1"
        write(unit, '(a)') "  modes = 1, 1, 2, -1"
        write(unit, '(a)') "/"
        write(unit, '(a)') ""
        write(unit, '(a)') "&background"
        write(unit, '(a)') "  rtor = 10.0"       ! Warning 6: very small tokamak
        write(unit, '(a)') "  rp = 9.99"         ! Warning 7: aspect ratio near 1
        write(unit, '(a)') "  B0 = 10.0"         ! Warning 8: extremely low field
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
    end subroutine create_config_multiple_warnings
    
    !> Clean up test files
    subroutine cleanup_test_files()
        call system("rm -f minor_errors.conf")
        call system("rm -f questionable.conf")
        call system("rm -f incomplete.conf")
        call system("rm -f minimal.conf")
        call system("rm -f multiple_warnings.conf")
        print *, ""
        print *, "Test files cleaned up"
    end subroutine cleanup_test_files

end program test_error_recovery