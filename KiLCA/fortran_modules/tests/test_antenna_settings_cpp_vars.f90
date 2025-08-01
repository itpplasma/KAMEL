program test_antenna_settings_cpp_vars
    use kilca_types_m, only: dp, KILCA_SUCCESS
    use kilca_settings_m, only: antenna_t, antenna_settings_set_modes, antenna_settings_get_modes
    implicit none
    
    type(antenna_t) :: as
    integer :: ierr
    integer :: test_status = 0
    
    print *, "=========================================="
    print *, "Testing Antenna Settings Complete C++ Variables"
    print *, "=========================================="
    
    ! This test should FAIL initially - not all C++ variables are present
    call test_all_cpp_variables_exist()
    call test_array_operations()
    
    if (test_status == 0) then
        print *, ""
        print *, "All tests PASSED"
        stop 0
    else
        print *, ""
        print *, "Tests FAILED:", test_status, "failure(s)"
        stop 1
    end if

contains

    !> Test that all C++ variables from antenna class exist
    subroutine test_all_cpp_variables_exist()
        print *, "Testing all C++ variables exist..."
        
        ! Test all required C++ variables exist
        as%ra = 90.0_dp
        as%wa = 0.0_dp
        as%I0 = 1.0e13_dp
        as%flab = (1.0e0_dp, 0.0e0_dp)
        as%dma = 5
        as%flag_debug = 1
        as%flag_eigmode = 0
        
        ! This should FAIL - modes array not yet implemented properly for all C++ variables
        ! According to C++ analysis, we need:
        ! - modes array (already exists but needs proper management)
        ! - All variables match C++ exactly
        
        ! Try to allocate modes array
        if (allocated(as%modes)) deallocate(as%modes)
        allocate(as%modes(10))  ! This should work
        as%modes = [1, 15, 2, 30, 3, 45, 4, 60, 5, 75]
        
        call assert_antenna_settings_complete(as, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Antenna settings not complete, error:", ierr
            test_status = test_status + 1
        else
            print *, "PASS: All antenna C++ variables exist"
        end if
    end subroutine test_all_cpp_variables_exist
    
    !> Test array operations
    subroutine test_array_operations()
        integer, parameter :: test_modes(8) = [1, 15, 2, 30, 3, 45, 4, 60]
        integer, allocatable :: modes_out(:)
        
        print *, "Testing antenna settings array operations..."
        
        ! This should FAIL initially - array operations not implemented
        call antenna_settings_set_modes(as, test_modes, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not set modes array, error:", ierr
            test_status = test_status + 1
            return
        end if
        
        call antenna_settings_get_modes(as, modes_out, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not get modes array, error:", ierr
            test_status = test_status + 1
            return
        end if
        
        call assert_arrays_equal(modes_out, test_modes, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Modes arrays not equal, error:", ierr
            test_status = test_status + 1
        else
            print *, "PASS: Array operations work correctly"
        end if
    end subroutine test_array_operations
    
    !> Assert that antenna settings structure is complete with all C++ variables
    subroutine assert_antenna_settings_complete(ant, ierr)
        type(antenna_t), intent(in) :: ant
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        ! Verify all C++ variables are accessible and have reasonable values
        ! Based on C++ analysis: ra, wa, I0, flab, dma, flag_debug, flag_eigmode, modes
        
        ! These should all compile and work
        if (ant%ra < 0.0_dp) then
            ierr = -1
            return
        end if
        
        if (ant%dma < 0) then
            ierr = -2
            return
        end if
        
        ! Verify complex frequency works
        if (real(ant%flab) /= 1.0e0_dp .or. aimag(ant%flab) /= 0.0e0_dp) then
            ierr = -3
            return
        end if
        
        ! Verify modes array is allocated and has correct size
        if (.not. allocated(ant%modes)) then
            ierr = -4
            return
        end if
        
        if (size(ant%modes) /= 10) then
            ierr = -5
            return
        end if
        
    end subroutine assert_antenna_settings_complete
    
    !> Assert that two integer arrays are equal
    subroutine assert_arrays_equal(arr1, arr2, ierr)
        integer, dimension(:), intent(in) :: arr1, arr2
        integer, intent(out) :: ierr
        
        ierr = KILCA_SUCCESS
        
        if (size(arr1) /= size(arr2)) then
            ierr = -1
            return
        end if
        
        if (any(arr1 /= arr2)) then
            ierr = -2
            return
        end if
        
    end subroutine assert_arrays_equal
    

end program test_antenna_settings_cpp_vars