program test_antenna_validation_enhanced
    use kilca_types_m, only: dp, KILCA_SUCCESS
    use kilca_settings_m, only: antenna_t, antenna_validate, antenna_initialize_defaults, &
                                antenna_settings_set_modes
    implicit none
    
    type(antenna_t) :: ant
    logical :: is_valid
    character(len=1024) :: error_msg
    integer :: ierr
    integer :: test_status = 0
    
    print *, "=========================================="
    print *, "Testing Enhanced Antenna Settings Validation"
    print *, "=========================================="
    
    call test_basic_validation()
    call test_physics_based_validation()
    call test_modes_validation()
    call test_frequency_validation()
    
    if (test_status == 0) then
        print *, ""
        print *, "All enhanced validation tests PASSED"
        stop 0
    else
        print *, ""
        print *, "Enhanced validation tests FAILED:", test_status, "failure(s)"
        stop 1
    end if

contains

    !> Test basic parameter validation
    subroutine test_basic_validation()
        print *, "Testing basic parameter validation..."
        
        ! Valid configuration should pass
        call antenna_initialize_defaults(ant, ierr)
        ant%ra = 50.0_dp
        ant%wa = 1.0_dp
        ant%I0 = 1000.0_dp
        ant%flab = (30.0e6_dp, 0.0_dp)
        ant%dma = 2
        ant%flag_debug = 1
        ant%flag_eigmode = 0
        call antenna_settings_set_modes(ant, [1, 15, 2, 30], ierr)
        
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (.not. is_valid .or. ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Valid configuration rejected:", trim(error_msg)
            test_status = test_status + 1
        else
            print *, "PASS: Valid configuration accepted"
        end if
        
        ! Test invalid radius
        ant%ra = -10.0_dp
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Negative radius not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Negative radius correctly rejected"
        end if
    end subroutine test_basic_validation
    
    !> Test physics-based validation enhancements
    subroutine test_physics_based_validation()
        print *, "Testing physics-based validation..."
        
        ! Reset to valid state
        call antenna_initialize_defaults(ant, ierr)
        ant%ra = 50.0_dp
        ant%wa = 1.0_dp
        ant%I0 = 1000.0_dp
        ant%flab = (30.0e6_dp, 0.0_dp)
        ant%dma = 1
        call antenna_settings_set_modes(ant, [1, 15], ierr)
        
        ! Test unreasonable radius (too small)
        ant%ra = 0.5_dp  ! Less than 1 cm
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably small radius not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably small radius correctly rejected"
        end if
        
        ! Test unreasonable radius (too large)
        ant%ra = 1500.0_dp  ! More than 1000 cm
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably large radius not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably large radius correctly rejected"
        end if
        
        ! Test wa >= ra (invalid physics)
        ant%ra = 50.0_dp
        ant%wa = 60.0_dp  ! Larger than radius
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Width >= radius not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Width >= radius correctly rejected"
        end if
        
        ! Test unreasonably large current
        ant%wa = 1.0_dp  ! Fix width
        ant%I0 = 2.0e15_dp  ! Unreasonably large
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably large current not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably large current correctly rejected"
        end if
    end subroutine test_physics_based_validation
    
    !> Test frequency validation
    subroutine test_frequency_validation()
        print *, "Testing frequency validation..."
        
        ! Reset to valid state
        call antenna_initialize_defaults(ant, ierr)
        ant%ra = 50.0_dp
        ant%wa = 1.0_dp
        ant%I0 = 1000.0_dp
        ant%dma = 1
        call antenna_settings_set_modes(ant, [1, 15], ierr)
        
        ! Test unreasonably low frequency
        ant%flab = (100.0_dp, 0.0_dp)  ! 100 Hz, too low
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably low frequency not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably low frequency correctly rejected"
        end if
        
        ! Test unreasonably high frequency
        ant%flab = (5.0e9_dp, 0.0_dp)  ! 5 GHz, too high for ICRF
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably high frequency not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably high frequency correctly rejected"
        end if
        
        ! Test valid ICRF frequency
        ant%flab = (50.0e6_dp, 0.0_dp)  ! 50 MHz, typical ICRF
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (.not. is_valid) then
            print *, "FAIL: Valid ICRF frequency rejected:", trim(error_msg)
            test_status = test_status + 1
        else
            print *, "PASS: Valid ICRF frequency accepted"
        end if
    end subroutine test_frequency_validation
    
    !> Test modes array validation
    subroutine test_modes_validation()
        print *, "Testing modes array validation..."
        
        ! Reset to valid state
        call antenna_initialize_defaults(ant, ierr)
        ant%ra = 50.0_dp
        ant%wa = 1.0_dp
        ant%I0 = 1000.0_dp
        ant%flab = (30.0e6_dp, 0.0_dp)
        
        ! Test unreasonably large m number
        call antenna_settings_set_modes(ant, [25, 15], ierr)  ! m=25 > 20
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably large m number not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably large m number correctly rejected"
        end if
        
        ! Test unreasonably large n number
        call antenna_settings_set_modes(ant, [1, 150], ierr)  ! n=150 > 100
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably large n number not detected" 
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably large n number correctly rejected"
        end if
        
        ! Test valid mode numbers
        call antenna_settings_set_modes(ant, [1, 15, -2, 30, 3, -45], ierr)
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (.not. is_valid) then
            print *, "FAIL: Valid mode numbers rejected:", trim(error_msg)
            test_status = test_status + 1
        else
            print *, "PASS: Valid mode numbers accepted"
        end if
        
        ! Test enhanced debug flag validation
        ant%flag_debug = 2  ! Should be valid (0, 1, or 2)
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (.not. is_valid) then
            print *, "FAIL: Valid debug flag=2 rejected:", trim(error_msg)
            test_status = test_status + 1
        else
            print *, "PASS: Valid debug flag=2 accepted"
        end if
        
        ! Test invalid debug flag
        ant%flag_debug = 3  ! Invalid
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Invalid debug flag=3 not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Invalid debug flag=3 correctly rejected"
        end if
    end subroutine test_modes_validation

end program test_antenna_validation_enhanced