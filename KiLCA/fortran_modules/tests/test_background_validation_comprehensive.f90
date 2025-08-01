program test_background_validation_comprehensive
    use kilca_types_m, only: dp, KILCA_SUCCESS
    use kilca_settings_m, only: back_sett_t, back_sett_validate, back_sett_initialize_defaults
    implicit none
    
    type(back_sett_t) :: bs
    logical :: is_valid
    character(len=1024) :: error_msg
    integer :: ierr
    integer :: test_status = 0
    
    print *, "=========================================="
    print *, "Testing Comprehensive Background Settings Validation"
    print *, "=========================================="
    
    call test_physics_range_validation()
    call test_collision_coefficient_validation()
    call test_parameter_cross_validation()
    call test_path_validation()
    
    if (test_status == 0) then
        print *, ""
        print *, "All comprehensive validation tests PASSED"
        stop 0
    else
        print *, ""
        print *, "Comprehensive validation tests FAILED:", test_status, "failure(s)"
        stop 1
    end if

contains

    !> Test physics-based range validation
    subroutine test_physics_range_validation()
        print *, "Testing physics-based range validation..."
        
        ! Valid configuration should pass
        call back_sett_initialize_defaults(bs, ierr)
        bs%rtor = 170.0_dp
        bs%rp = 70.0_dp
        bs%B0 = 25000.0_dp
        bs%m_i = 2.0_dp
        bs%zele = 1.0_dp
        bs%zion = 1.0_dp
        bs%V_gal_sys = 1.0e9_dp
        bs%huge_factor = 1.0e20_dp
        
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (.not. is_valid .or. ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Valid configuration rejected:", trim(error_msg)
            test_status = test_status + 1
        else
            print *, "PASS: Valid configuration accepted"
        end if
        
        ! Test torus radius too small
        bs%rtor = 30.0_dp  ! Below 50 cm minimum
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably small torus radius not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably small torus radius correctly rejected"
        end if
        
        ! Test torus radius too large
        bs%rtor = 1200.0_dp  ! Above 1000 cm maximum
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably large torus radius not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably large torus radius correctly rejected"
        end if
        
        ! Reset and test plasma radius too large
        bs%rtor = 170.0_dp
        bs%rp = 250.0_dp  ! Above 200 cm maximum
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably large plasma radius not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably large plasma radius correctly rejected"
        end if
        
        ! Test magnetic field too small
        bs%rp = 70.0_dp
        bs%B0 = 500.0_dp  ! Below 1000 minimum (0.05 T)
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably small magnetic field not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably small magnetic field correctly rejected"
        end if
        
        ! Test magnetic field too large
        bs%B0 = 150000.0_dp  ! Above 100000 maximum (15 T)
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably large magnetic field not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably large magnetic field correctly rejected"
        end if
        
        ! Test ion mass too small
        bs%B0 = 25000.0_dp
        bs%m_i = 0.3_dp  ! Below 0.5 amu minimum
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably small ion mass not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably small ion mass correctly rejected"
        end if
        
        ! Test ion mass too large
        bs%m_i = 100.0_dp  ! Above 50 amu maximum
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably large ion mass not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably large ion mass correctly rejected"
        end if
    end subroutine test_physics_range_validation
    
    !> Test collision coefficient validation
    subroutine test_collision_coefficient_validation()
        print *, "Testing collision coefficient validation..."
        
        ! Reset to valid state
        call back_sett_initialize_defaults(bs, ierr)
        bs%rtor = 170.0_dp
        bs%rp = 70.0_dp
        bs%B0 = 25000.0_dp
        bs%m_i = 2.0_dp
        bs%V_gal_sys = 1.0e9_dp
        bs%huge_factor = 1.0e20_dp
        
        ! Test electron collision coefficient too large
        bs%zele = 150.0_dp  ! Above 100 maximum
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably large electron collision coefficient not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably large electron collision coefficient correctly rejected"
        end if
        
        ! Test ion collision coefficient too large
        bs%zele = 1.0_dp
        bs%zion = 150.0_dp  ! Above 100 maximum
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably large ion collision coefficient not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably large ion collision coefficient correctly rejected"
        end if
    end subroutine test_collision_coefficient_validation
    
    !> Test parameter cross-validation
    subroutine test_parameter_cross_validation()
        print *, "Testing parameter cross-validation..."
        
        ! Reset to valid state
        call back_sett_initialize_defaults(bs, ierr)
        bs%rtor = 170.0_dp
        bs%rp = 70.0_dp
        bs%B0 = 25000.0_dp
        bs%m_i = 2.0_dp
        bs%zele = 1.0_dp
        bs%zion = 1.0_dp
        bs%huge_factor = 1.0e20_dp
        
        ! Test galaxy system voltage too large
        bs%V_gal_sys = 2.0e13_dp  ! Above 1e12 maximum
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably large galaxy system voltage not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably large galaxy system voltage correctly rejected"
        end if
        
        ! Test huge factor zero
        bs%V_gal_sys = 1.0e9_dp
        bs%huge_factor = 0.0_dp
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Zero huge factor not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Zero huge factor correctly rejected"
        end if
        
        ! Test huge factor too large
        bs%huge_factor = 2.0e51_dp  ! Above 1e50 maximum
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Unreasonably large huge factor not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Unreasonably large huge factor correctly rejected"
        end if
    end subroutine test_parameter_cross_validation
    
    !> Test path validation
    subroutine test_path_validation()
        print *, "Testing path validation..."
        
        ! Reset to valid state
        call back_sett_initialize_defaults(bs, ierr)
        bs%rtor = 170.0_dp
        bs%rp = 70.0_dp
        bs%B0 = 25000.0_dp
        bs%m_i = 2.0_dp
        bs%zele = 1.0_dp
        bs%zion = 1.0_dp
        bs%V_gal_sys = 1.0e9_dp
        bs%huge_factor = 1.0e20_dp
        
        ! Test empty path string
        if (allocated(bs%path2profiles)) deallocate(bs%path2profiles)
        bs%path2profiles = ""  ! Empty string
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Empty path string not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Empty path string correctly rejected"
        end if
        
        ! Test valid path
        bs%path2profiles = "../profiles/"
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (.not. is_valid) then
            print *, "FAIL: Valid path string rejected:", trim(error_msg)
            test_status = test_status + 1
        else
            print *, "PASS: Valid path string accepted"
        end if
    end subroutine test_path_validation

end program test_background_validation_comprehensive