program test_eigenmode_settings_cpp_vars
    use kilca_types_m, only: dp, KILCA_SUCCESS
    use kilca_settings_m, only: eigmode_sett_t, eigmode_sett_initialize_defaults, &
                                eigmode_sett_validate
    implicit none
    
    type(eigmode_sett_t) :: es
    integer :: ierr
    integer :: test_status = 0
    logical :: is_valid
    character(len=1024) :: error_msg
    
    print *, "=========================================="
    print *, "Testing Eigenmode Settings Complete C++ Variables [RED PHASE - SHOULD FAIL]"
    print *, "=========================================="
    
    call test_all_cpp_variables_exist()
    call test_dynamic_guess_array_operations()
    call test_eigenmode_advanced_validation()
    
    if (test_status == 0) then
        print *, ""
        print *, "All eigenmode settings tests PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "Eigenmode settings tests FAILED:", test_status, "failure(s) (expected in RED phase)"
        stop 1
    end if

contains

    !> Test that all C++ variables from eigmode_sett class exist and work
    subroutine test_all_cpp_variables_exist()
        print *, "Testing all eigenmode C++ variables exist..."
        
        ! Initialize with defaults
        call eigmode_sett_initialize_defaults(es, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not initialize eigenmode settings"
            test_status = test_status + 1
            return
        end if
        
        ! Test all C++ variables are accessible and settable
        es%search_flag = 1
        es%rdim = 50
        es%rfmin = 1.0e6_dp
        es%rfmax = 1.0e10_dp
        es%idim = 25
        es%ifmin = -1.0e7_dp
        es%ifmax = 1.0e7_dp
        es%stop_flag = 1
        es%eps_res = 1.0e-8_dp
        es%eps_abs = 1.0e-10_dp
        es%eps_rel = 1.0e-8_dp
        es%delta = 1.0e-8_dp
        es%test_roots = 1
        es%flag_debug = 1
        es%Nguess = 5
        es%kmin = 0
        es%kmax = 20
        es%n_zeros = 25
        es%use_winding = 1
        
        ! Dynamic array for starting frequencies
        if (allocated(es%fstart)) deallocate(es%fstart)
        allocate(es%fstart(5))
        es%fstart(1) = (1.0e6_dp, 0.0_dp)
        es%fstart(2) = (1.1e6_dp, 0.0_dp)
        es%fstart(3) = (1.2e6_dp, 0.0_dp)
        es%fstart(4) = (1.3e6_dp, 0.0_dp)
        es%fstart(5) = (1.4e6_dp, 0.0_dp)
        
        ! Validate the structure
        call eigmode_sett_validate(es, is_valid, error_msg, ierr)
        if (.not. is_valid .or. ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Eigenmode settings validation failed:", trim(error_msg)
            test_status = test_status + 1
        else
            print *, "PASS: All eigenmode C++ variables exist and validation works"
        end if
    end subroutine test_all_cpp_variables_exist
    
    !> Test dynamic array operations for guess frequencies that should fail initially
    subroutine test_dynamic_guess_array_operations()
        print *, "Testing dynamic guess array operations..."
        
        call eigmode_sett_initialize_defaults(es, ierr)
        es%Nguess = 3
        
        ! This should FAIL initially - advanced array management not implemented
        block
            complex(dp) :: test_freqs(3)
            test_freqs(1) = (1.5e6_dp, 0.1e6_dp)
            test_freqs(2) = (2.0e6_dp, 0.2e6_dp)
            test_freqs(3) = (2.5e6_dp, 0.3e6_dp)
            call eigenmode_settings_set_guess_frequencies(es, test_freqs, ierr)
        end block
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not set guess frequencies array, error:", ierr
            test_status = test_status + 1
        else
            print *, "PASS: Advanced guess frequency array management works correctly"
        end if
        
        ! Test getting array back
        block
            complex(dp), allocatable :: guess_freqs(:)
            call eigenmode_settings_get_guess_frequencies(es, guess_freqs, ierr)
            if (ierr /= KILCA_SUCCESS) then
                print *, "FAIL: Could not get guess frequencies array, error:", ierr
                test_status = test_status + 1
            else if (.not. allocated(guess_freqs)) then
                print *, "FAIL: Returned guess frequencies array not allocated"
                test_status = test_status + 1
            else if (size(guess_freqs) /= 3) then
                print *, "FAIL: Returned guess frequencies array wrong size:", size(guess_freqs)
                test_status = test_status + 1
            else
                print *, "PASS: Guess frequencies array get operation works correctly"
            end if
        end block
    end subroutine test_dynamic_guess_array_operations
    
    !> Test enhanced eigenmode validation that should fail initially
    subroutine test_eigenmode_advanced_validation()
        print *, "Testing enhanced eigenmode validation..."
        
        call eigmode_sett_initialize_defaults(es, ierr)
        
        ! Test inconsistent Nguess vs fstart array size
        es%Nguess = 5
        if (allocated(es%fstart)) deallocate(es%fstart)
        allocate(es%fstart(3))  ! Wrong size
        es%fstart(1) = (1.0e6_dp, 0.0_dp)
        es%fstart(2) = (2.0e6_dp, 0.0_dp)
        es%fstart(3) = (3.0e6_dp, 0.0_dp)
        
        call eigmode_sett_validate(es, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: Nguess vs fstart array size mismatch not detected"
            test_status = test_status + 1
        else
            print *, "PASS: Nguess vs fstart array size mismatch correctly detected:", trim(error_msg)
        end if
        
        ! Test correct size
        es%Nguess = 3
        call eigmode_sett_validate(es, is_valid, error_msg, ierr)
        if (.not. is_valid) then
            print *, "FAIL: Valid Nguess vs fstart array size rejected:", trim(error_msg)
            test_status = test_status + 1
        else
            print *, "PASS: Valid Nguess vs fstart array size accepted"
        end if
        
        ! Test enhanced logical consistency validation that should be missing
        es%kmin = 15
        es%kmax = 5  ! kmin > kmax should be invalid
        call eigmode_sett_validate(es, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: kmin > kmax logical inconsistency not detected"
            test_status = test_status + 1
        else
            print *, "PASS: kmin > kmax logical inconsistency correctly detected:", trim(error_msg)
        end if
        
        ! Test frequency range validation that should be missing
        es%kmin = 1
        es%kmax = 10
        es%rfmin = 1.0e9_dp
        es%rfmax = 1.0e6_dp  ! rfmin > rfmax should be invalid
        call eigmode_sett_validate(es, is_valid, error_msg, ierr)
        if (is_valid) then
            print *, "FAIL: rfmin > rfmax range inconsistency not detected"
            test_status = test_status + 1
        else
            print *, "PASS: rfmin > rfmax range inconsistency correctly detected:", trim(error_msg)
        end if
    end subroutine test_eigenmode_advanced_validation

end program test_eigenmode_settings_cpp_vars