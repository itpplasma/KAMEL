program test_background_settings_cpp_vars
    use kilca_types_m, only: dp, KILCA_SUCCESS, m_proton, m_electron, e_charge
    use kilca_settings_m, only: back_sett_t, back_sett_validate, back_sett_initialize_defaults, &
                                background_settings_compute_derived
    implicit none
    
    type(back_sett_t) :: bs
    logical :: is_valid
    character(len=1024) :: error_msg
    integer :: ierr
    integer :: test_status = 0
    
    print *, "=========================================="
    print *, "Testing Background Settings Complete C++ Variables"
    print *, "=========================================="
    
    call test_all_cpp_variables_exist()
    call test_derived_values_computation()
    call test_physics_constraints()
    
    if (test_status == 0) then
        print *, ""
        print *, "All background tests PASSED"
        stop 0
    else
        print *, ""
        print *, "Background tests FAILED:", test_status, "failure(s)"
        stop 1
    end if

contains

    !> Test that all C++ variables from back_sett class exist and work
    subroutine test_all_cpp_variables_exist()
        print *, "Testing all C++ variables exist..."
        
        ! Initialize with defaults
        call back_sett_initialize_defaults(bs, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not initialize background settings"
            test_status = test_status + 1
            return
        end if
        
        ! Test all C++ variables are accessible and settable
        bs%rtor = 170.69_dp
        bs%rp = 70.0_dp
        bs%B0 = 23176.46_dp
        if (allocated(bs%path2profiles)) deallocate(bs%path2profiles)
        bs%path2profiles = "../profiles/"
        bs%calc_back = 1
        if (allocated(bs%flag_back)) deallocate(bs%flag_back)
        bs%flag_back = "f"
        bs%N = 9
        bs%V_gal_sys = 1.0e9_dp
        bs%V_scale = 1.0e0_dp
        bs%m_i = 2.0_dp
        bs%zele = 1.0e-0_dp
        bs%zion = 1.0e-0_dp
        bs%flag_debug = 0
        bs%huge_factor = 1.0e20_dp
        
        ! This should FAIL initially - derived computation not implemented
        call background_settings_compute_derived(bs, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not compute derived values, error:", ierr
            test_status = test_status + 1
        else
            print *, "PASS: All background C++ variables exist and derived computation works"
        end if
    end subroutine test_all_cpp_variables_exist
    
    !> Test derived values computation (mass and charge arrays)
    subroutine test_derived_values_computation()
        print *, "Testing derived values computation..."
        
        ! Set up test case
        call back_sett_initialize_defaults(bs, ierr)
        bs%m_i = 2.0_dp  ! Deuterium
        
        ! This should work after implementation
        call background_settings_compute_derived(bs, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Derived computation failed, error:", ierr
            test_status = test_status + 1
            return
        end if
        
        ! Verify mass array computed correctly
        if (.not. allocated(bs%mass) .or. size(bs%mass) /= 2) then
            print *, "FAIL: Mass array not properly allocated"
            test_status = test_status + 1
            return
        end if
        
        ! Check ion mass: should be m_i * proton_mass
        if (abs(bs%mass(1) - bs%m_i * m_proton) > 1.0e-15_dp) then
            print *, "FAIL: Ion mass incorrect:", bs%mass(1), "expected:", bs%m_i * m_proton
            test_status = test_status + 1
            return
        end if
        
        ! Check electron mass
        if (abs(bs%mass(2) - m_electron) > 1.0e-15_dp) then
            print *, "FAIL: Electron mass incorrect:", bs%mass(2), "expected:", m_electron
            test_status = test_status + 1
            return
        end if
        
        ! Verify charge array computed correctly
        if (.not. allocated(bs%charge) .or. size(bs%charge) /= 2) then
            print *, "FAIL: Charge array not properly allocated"
            test_status = test_status + 1
            return
        end if
        
        ! Check ion charge: should be +e
        if (abs(bs%charge(1) - e_charge) > 1.0e-15_dp) then
            print *, "FAIL: Ion charge incorrect:", bs%charge(1), "expected:", e_charge
            test_status = test_status + 1
            return
        end if
        
        ! Check electron charge: should be -e
        if (abs(bs%charge(2) - (-e_charge)) > 1.0e-15_dp) then
            print *, "FAIL: Electron charge incorrect:", bs%charge(2), "expected:", -e_charge
            test_status = test_status + 1
            return
        end if
        
        print *, "PASS: Derived values computed correctly"
    end subroutine test_derived_values_computation
    
    !> Test physics constraints and validation
    subroutine test_physics_constraints()
        print *, "Testing physics constraints..."
        
        ! Set up valid configuration
        call back_sett_initialize_defaults(bs, ierr)
        bs%rtor = 170.0_dp
        bs%rp = 70.0_dp
        bs%B0 = 25000.0_dp
        
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (.not. is_valid .or. ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Valid configuration rejected:", trim(error_msg)
            test_status = test_status + 1
        else
            print *, "PASS: Valid configuration accepted"
        end if
        
        ! Test rp >= rtor constraint
        bs%rp = 180.0_dp  ! Larger than rtor
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: rp >= rtor not detected"
            test_status = test_status + 1
        else
            print *, "PASS: rp >= rtor correctly rejected"
        end if
    end subroutine test_physics_constraints
    

end program test_background_settings_cpp_vars