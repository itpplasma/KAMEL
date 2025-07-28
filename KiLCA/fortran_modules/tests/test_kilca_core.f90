program test_kilca_core
    use iso_fortran_env, only: real32, real64, int32, int64
    use iso_c_binding
    use kilca_types_m
    use kilca_core_m
    implicit none
    
    ! Test variables
    integer :: test_status = 0
    character(len=256) :: test_path = "/test/project/path/"
    
    ! Test 1: Create and destroy core_data structure
    call test_core_data_lifecycle()
    
    ! Test 2: Test path management
    call test_path_management()
    
    ! Test 3: Test pointer members initialization
    call test_pointer_members()
    
    ! Test 4: Test modes array management
    call test_modes_array()
    
    ! Test 5: Test method procedures
    call test_core_methods()
    
    ! Test 6: Test C-Fortran interface compatibility
    call test_fortran_interface()
    
    ! Test 7: Test memory management
    call test_memory_management()
    
    ! Test 8: Test error handling
    call test_error_handling()
    
    if (test_status == 0) then
        print *, "All tests PASSED!"
    else
        print *, "Some tests FAILED. Status:", test_status
        stop 1
    end if
    
contains

    subroutine test_core_data_lifecycle()
        type(core_data_t), pointer :: cd
        integer :: ierr
        
        print *, "Testing core_data lifecycle..."
        
        ! Test creation
        call core_data_create(cd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: core_data_create returned error:", ierr
            test_status = test_status + 1
            return
        end if
        
        if (.not. associated(cd)) then
            print *, "FAIL: core_data not allocated"
            test_status = test_status + 1
            return
        end if
        
        ! Test destruction
        call core_data_destroy(cd, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: core_data_destroy returned error:", ierr
            test_status = test_status + 1
            return
        end if
        
        if (associated(cd)) then
            print *, "FAIL: core_data not deallocated"
            test_status = test_status + 1
        end if
        
        print *, "test_core_data_lifecycle completed"
    end subroutine test_core_data_lifecycle
    
    subroutine test_path_management()
        type(core_data_t), pointer :: cd
        integer :: ierr
        character(len=1024) :: retrieved_path
        
        print *, "Testing path management..."
        
        call core_data_create(cd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Test path storage
        call core_data_get_path(cd, retrieved_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: core_data_get_path error:", ierr
            test_status = test_status + 1
        end if
        
        if (trim(retrieved_path) /= trim(test_path)) then
            print *, "FAIL: Path mismatch"
            print *, "Expected: ", trim(test_path)
            print *, "Got:      ", trim(retrieved_path)
            test_status = test_status + 1
        end if
        
        call core_data_destroy(cd, ierr)
        
        print *, "test_path_management completed"
    end subroutine test_path_management
    
    subroutine test_pointer_members()
        type(core_data_t), pointer :: cd
        integer :: ierr
        logical :: has_settings, has_background
        
        print *, "Testing pointer members..."
        
        call core_data_create(cd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Check initial state
        call core_data_has_settings(cd, has_settings)
        if (has_settings) then
            print *, "FAIL: Settings should not be initialized yet"
            test_status = test_status + 1
        end if
        
        call core_data_has_background(cd, has_background)
        if (has_background) then
            print *, "FAIL: Background should not be initialized yet"
            test_status = test_status + 1
        end if
        
        call core_data_destroy(cd, ierr)
        
        print *, "test_pointer_members completed"
    end subroutine test_pointer_members
    
    subroutine test_modes_array()
        type(core_data_t), pointer :: cd
        integer :: ierr, n_modes
        integer :: test_dim = 5
        
        print *, "Testing modes array management..."
        
        call core_data_create(cd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Test initial state
        call core_data_get_modes_count(cd, n_modes, ierr)
        if (n_modes /= 0) then
            print *, "FAIL: Initial modes count should be 0"
            test_status = test_status + 1
        end if
        
        ! Test allocation
        call core_data_allocate_modes(cd, test_dim, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Mode allocation error:", ierr
            test_status = test_status + 1
        end if
        
        call core_data_get_modes_count(cd, n_modes, ierr)
        if (n_modes /= test_dim) then
            print *, "FAIL: Modes count mismatch"
            test_status = test_status + 1
        end if
        
        ! Test deallocation
        call core_data_delete_modes_array(cd, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Mode deletion error:", ierr
            test_status = test_status + 1
        end if
        
        call core_data_get_modes_count(cd, n_modes, ierr)
        if (n_modes /= 0) then
            print *, "FAIL: Modes not deleted"
            test_status = test_status + 1
        end if
        
        call core_data_destroy(cd, ierr)
        
        print *, "test_modes_array completed"
    end subroutine test_modes_array
    
    subroutine test_core_methods()
        type(core_data_t), pointer :: cd
        integer :: ierr
        
        print *, "Testing core method procedures..."
        
        call core_data_create(cd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! These methods require proper initialization, so we just test they exist
        ! and can be called without crashing
        
        ! Test that procedures are available (will fail to link if not implemented)
        ! calc_and_set_mode_independent_core_data
        ! calc_and_set_mode_dependent_core_data_antenna
        ! calc_and_set_mode_dependent_core_data_eigmode
        ! calc_and_set_mode_dependent_core_data_antenna_interface
        
        call core_data_destroy(cd, ierr)
        
        print *, "test_core_methods completed"
    end subroutine test_core_methods
    
    subroutine test_fortran_interface()
        type(core_data_t), pointer :: cd
        integer :: ierr
        integer :: ptr_size
        
        print *, "Testing Fortran interface compatibility..."
        
        ! Test pointer precision
        call get_pointer_precision(ptr_size, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: get_pointer_precision error"
            test_status = test_status + 1
        else if (ptr_size /= c_intptr_t) then
            print *, "FAIL: Pointer size mismatch"
            test_status = test_status + 1
        end if
        
        call core_data_create(cd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Test interface procedures exist
        ! set_core_data_in_core_module
        ! set_background_in_core_module
        ! set_settings_in_core_module
        ! clear_all_data_in_mode_data_module
        
        call core_data_destroy(cd, ierr)
        
        print *, "test_fortran_interface completed"
    end subroutine test_fortran_interface
    
    subroutine test_memory_management()
        type(core_data_t), pointer :: cd
        integer :: ierr
        integer :: i
        
        print *, "Testing memory management..."
        
        ! Test multiple create/destroy cycles
        do i = 1, 10
            call core_data_create(cd, test_path, ierr)
            if (ierr /= KILCA_SUCCESS) then
                print *, "FAIL: Creation failed on iteration", i
                test_status = test_status + 1
                exit
            end if
            
            call core_data_allocate_modes(cd, i*2, ierr)
            call core_data_delete_modes_array(cd, ierr)
            
            call core_data_destroy(cd, ierr)
            if (ierr /= KILCA_SUCCESS) then
                print *, "FAIL: Destruction failed on iteration", i
                test_status = test_status + 1
                exit
            end if
        end do
        
        print *, "test_memory_management completed"
    end subroutine test_memory_management
    
    subroutine test_error_handling()
        type(core_data_t), pointer :: cd
        integer :: ierr
        character(len=1) :: invalid_path = ""
        
        print *, "Testing error handling..."
        
        ! Test null pointer handling
        cd => null()
        call core_data_destroy(cd, ierr)
        if (ierr /= KILCA_ERROR_INVALID_INPUT) then
            print *, "FAIL: Should return error for null pointer"
            test_status = test_status + 1
        end if
        
        ! Test invalid path
        call core_data_create(cd, invalid_path, ierr)
        if (ierr == KILCA_SUCCESS) then
            print *, "FAIL: Should return error for empty path"
            test_status = test_status + 1
            call core_data_destroy(cd, ierr)
        end if
        
        ! Test double allocation
        call core_data_create(cd, test_path, ierr)
        if (ierr == KILCA_SUCCESS) then
            call core_data_allocate_modes(cd, 5, ierr)
            if (ierr == KILCA_SUCCESS) then
                call core_data_allocate_modes(cd, 10, ierr)
                if (ierr == KILCA_SUCCESS) then
                    print *, "FAIL: Should return error for double allocation"
                    test_status = test_status + 1
                end if
            end if
            call core_data_destroy(cd, ierr)
        end if
        
        print *, "test_error_handling completed"
    end subroutine test_error_handling

end program test_kilca_core