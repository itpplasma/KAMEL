program test_error_recovery_simple
    use kilca_types_m, only: dp, KILCA_SUCCESS
    use kilca_settings_m, only: settings_t, settings_validate_with_warnings
    implicit none
    
    type(settings_t) :: settings
    integer :: ierr
    integer :: test_status = 0
    
    print *, "=========================================="
    print *, "Testing Error Recovery Functions [REFACTOR PHASE]"
    print *, "=========================================="
    
    call test_warning_system_direct()
    
    if (test_status == 0) then
        print *, ""
        print *, "Error recovery functions test PASSED"
        stop 0
    else
        print *, ""
        print *, "Error recovery functions test FAILED:", test_status, "failure(s)"
        stop 1
    end if

contains

    !> Test warning system directly without file I/O
    subroutine test_warning_system_direct()
        integer :: warning_count
        character(len=256), dimension(10) :: warnings
        integer :: i
        
        print *, "Testing warning system with direct values..."
        
        ! Set questionable values
        settings%antenna_settings%ra = 0.5_dp      ! Very small
        settings%antenna_settings%wa = 200.0_dp    ! Very large
        settings%antenna_settings%I0 = 1.0e19_dp   ! Extremely high
        settings%antenna_settings%flab = cmplx(5000.0_dp, 0.0_dp, dp)  ! Very low
        settings%antenna_settings%dma = 25         ! Many modes
        
        settings%background_settings%rtor = 40.0_dp  ! Small tokamak
        settings%background_settings%rp = 38.0_dp    ! Low aspect ratio
        settings%background_settings%B0 = 500.0_dp   ! Low field
        
        settings%eigmode_settings%eps_res = 1.0e-15_dp  ! Very tight
        
        ! Validate with warnings
        call settings_validate_with_warnings(settings, ierr, warning_count, warnings)
        
        if (ierr == KILCA_SUCCESS) then
            print *, "PASS: Validation completed"
        else
            print *, "INFO: Validation returned error:", ierr
        end if
        
        if (warning_count > 0) then
            print *, "PASS: Generated", warning_count, "warnings:"
            do i = 1, warning_count
                print *, "  -", trim(warnings(i))
            end do
        else
            print *, "FAIL: Should have generated warnings for questionable values"
            test_status = test_status + 1
        end if
        
        ! Test with valid values (no warnings expected)
        print *, ""
        print *, "Testing with normal values (no warnings expected)..."
        settings%antenna_settings%ra = 90.0_dp
        settings%antenna_settings%wa = 5.0_dp
        settings%antenna_settings%I0 = 1.0e12_dp
        settings%antenna_settings%flab = cmplx(1.0e6_dp, 0.0_dp, dp)
        settings%antenna_settings%dma = 2
        
        settings%background_settings%rtor = 170.0_dp
        settings%background_settings%rp = 65.0_dp
        settings%background_settings%B0 = 25000.0_dp
        
        settings%eigmode_settings%eps_res = 1.0e-8_dp
        
        call settings_validate_with_warnings(settings, ierr, warning_count, warnings)
        
        if (warning_count == 0) then
            print *, "PASS: No warnings for normal values"
        else
            print *, "INFO: Generated", warning_count, "warnings for normal values"
            do i = 1, warning_count
                print *, "  -", trim(warnings(i))
            end do
        end if
        
    end subroutine test_warning_system_direct

end program test_error_recovery_simple