program test_core_deep_copy
    use iso_fortran_env, only: real64, int32, int64
    use iso_c_binding
    use kilca_types_m
    use kilca_core_m
    use kilca_settings_m
    implicit none
    
    ! Test variables
    integer :: test_status = 0
    character(len=256) :: test_path = "/test/project/path/"
    
    ! Test 1: Deep copy empty core_data
    call test_deep_copy_empty()
    
    ! Test 2: Deep copy with settings
    call test_deep_copy_with_settings()
    
    ! Test 3: Deep copy with mode array
    call test_deep_copy_with_modes()
    
    ! Test 4: Deep copy full structure
    call test_deep_copy_full()
    
    ! Test 5: Test independence after copy
    call test_copy_independence()
    
    ! Test 6: Test copy assignment operator
    call test_copy_assignment()
    
    if (test_status == 0) then
        print *, "All deep copy tests PASSED!"
    else
        print *, "Some tests FAILED. Status:", test_status
        stop 1
    end if
    
contains

    subroutine test_deep_copy_empty()
        type(core_data_t), pointer :: cd_src, cd_dst
        integer :: ierr
        
        print *, "Testing deep copy of empty core_data..."
        
        ! Create source
        call core_data_create(cd_src, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not create source"
            test_status = test_status + 1
            return
        end if
        
        ! Deep copy
        call core_data_deep_copy(cd_src, cd_dst, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Deep copy failed"
            test_status = test_status + 1
            call core_data_destroy(cd_src, ierr)
            return
        end if
        
        ! Verify copy
        if (cd_dst%path2project /= cd_src%path2project) then
            print *, "FAIL: Path not copied correctly"
            test_status = test_status + 1
        end if
        
        if (cd_dst%dim /= cd_src%dim) then
            print *, "FAIL: Dimension not copied"
            test_status = test_status + 1
        end if
        
        ! Clean up
        call core_data_destroy(cd_src, ierr)
        call core_data_destroy(cd_dst, ierr)
        
        print *, "test_deep_copy_empty completed"
    end subroutine test_deep_copy_empty
    
    subroutine test_deep_copy_with_settings()
        type(core_data_t), pointer :: cd_src, cd_dst
        type(settings_t), pointer :: sd
        integer :: ierr
        
        print *, "Testing deep copy with settings..."
        
        ! Create source with settings
        call core_data_create(cd_src, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Create and attach settings
        call settings_create(sd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            call core_data_destroy(cd_src, ierr)
            return
        end if
        
        cd_src%sd => sd
        cd_src%owns_settings = .true.
        
        ! Deep copy
        call core_data_deep_copy(cd_src, cd_dst, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Deep copy with settings failed"
            test_status = test_status + 1
            call core_data_destroy(cd_src, ierr)
            return
        end if
        
        ! Verify settings were copied
        if (.not. associated(cd_dst%sd)) then
            print *, "FAIL: Settings not copied"
            test_status = test_status + 1
        else
            ! Verify settings are independent
            if (associated(cd_dst%sd, cd_src%sd)) then
                print *, "FAIL: Settings not deep copied (same pointer)"
                test_status = test_status + 1
            end if
        end if
        
        ! Clean up
        call core_data_destroy(cd_src, ierr)
        call core_data_destroy(cd_dst, ierr)
        
        print *, "test_deep_copy_with_settings completed"
    end subroutine test_deep_copy_with_settings
    
    subroutine test_deep_copy_with_modes()
        type(core_data_t), pointer :: cd_src, cd_dst
        integer :: ierr, i
        
        print *, "Testing deep copy with mode array..."
        
        ! Create source
        call core_data_create(cd_src, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Allocate mode array
        cd_src%dim = 5
        allocate(cd_src%mda(cd_src%dim))
        do i = 1, cd_src%dim
            cd_src%mda(i)%initialized = .true.
        end do
        
        ! Deep copy
        call core_data_deep_copy(cd_src, cd_dst, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Deep copy with modes failed"
            test_status = test_status + 1
            call core_data_destroy(cd_src, ierr)
            return
        end if
        
        ! Verify modes were copied
        if (.not. allocated(cd_dst%mda)) then
            print *, "FAIL: Mode array not allocated in copy"
            test_status = test_status + 1
        else
            if (size(cd_dst%mda) /= size(cd_src%mda)) then
                print *, "FAIL: Mode array size mismatch"
                test_status = test_status + 1
            end if
            
            ! Check all modes
            do i = 1, cd_dst%dim
                if (.not. cd_dst%mda(i)%initialized) then
                    print *, "FAIL: Mode", i, "not initialized in copy"
                    test_status = test_status + 1
                end if
            end do
        end if
        
        ! Clean up
        call core_data_destroy(cd_src, ierr)
        call core_data_destroy(cd_dst, ierr)
        
        print *, "test_deep_copy_with_modes completed"
    end subroutine test_deep_copy_with_modes
    
    subroutine test_deep_copy_full()
        type(core_data_t), pointer :: cd_src, cd_dst
        type(settings_t), pointer :: sd
        integer :: ierr
        
        print *, "Testing deep copy of full structure..."
        
        ! Create fully populated source
        call core_data_create(cd_src, "/test/vacuum/project/", ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Initialize mode independent data
        call calc_and_set_mode_independent_core_data(cd_src, ierr)
        
        ! Create settings with antenna data
        call settings_create(sd, cd_src%path2project, ierr)
        if (ierr == KILCA_SUCCESS) then
            sd%antenna_settings%dma = 3
            allocate(sd%antenna_settings%modes(6))
            sd%antenna_settings%modes = [1, 1, 2, 2, 3, 3]
            sd%antenna_settings%flab = cmplx(50.0e6_dp, 0.0_dp, dp)
            cd_src%sd => sd
            cd_src%owns_settings = .true.
            
            ! Initialize mode dependent data
            call calc_and_set_mode_dependent_core_data_antenna(cd_src, ierr)
        end if
        
        ! Deep copy
        call core_data_deep_copy(cd_src, cd_dst, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Full deep copy failed"
            test_status = test_status + 1
            call core_data_destroy(cd_src, ierr)
            return
        end if
        
        ! Verify everything was copied
        if (cd_dst%path2project /= cd_src%path2project) then
            print *, "FAIL: Path mismatch in full copy"
            test_status = test_status + 1
        end if
        
        if (cd_dst%dim /= cd_src%dim) then
            print *, "FAIL: Dimension mismatch in full copy"
            test_status = test_status + 1
        end if
        
        if (.not. associated(cd_dst%sd)) then
            print *, "FAIL: Settings not present in full copy"
            test_status = test_status + 1
        end if
        
        if (.not. associated(cd_dst%bp)) then
            print *, "FAIL: Background not present in full copy"
            test_status = test_status + 1
        end if
        
        ! Clean up
        call core_data_destroy(cd_src, ierr)
        call core_data_destroy(cd_dst, ierr)
        
        print *, "test_deep_copy_full completed"
    end subroutine test_deep_copy_full
    
    subroutine test_copy_independence()
        type(core_data_t), pointer :: cd_src, cd_dst
        character(len=256) :: new_path
        integer :: ierr
        
        print *, "Testing independence after copy..."
        
        ! Create source
        call core_data_create(cd_src, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Deep copy
        call core_data_deep_copy(cd_src, cd_dst, ierr)
        if (ierr /= KILCA_SUCCESS) then
            call core_data_destroy(cd_src, ierr)
            return
        end if
        
        ! Modify source
        new_path = "/modified/path/"
        cd_src%path2project = new_path
        cd_src%dim = 999
        
        ! Verify destination unchanged
        if (cd_dst%path2project == new_path) then
            print *, "FAIL: Destination affected by source modification"
            test_status = test_status + 1
        end if
        
        if (cd_dst%dim == 999) then
            print *, "FAIL: Destination dim affected by source"
            test_status = test_status + 1
        end if
        
        ! Clean up
        call core_data_destroy(cd_src, ierr)
        call core_data_destroy(cd_dst, ierr)
        
        print *, "test_copy_independence completed"
    end subroutine test_copy_independence
    
    subroutine test_copy_assignment()
        type(core_data_t), pointer :: cd_src, cd_dst
        integer :: ierr
        
        print *, "Testing copy assignment..."
        
        ! Create source and destination
        call core_data_create(cd_src, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        call core_data_create(cd_dst, "/old/path/", ierr)
        if (ierr /= KILCA_SUCCESS) then
            call core_data_destroy(cd_src, ierr)
            return
        end if
        
        ! Copy assignment
        call core_data_copy_assign(cd_dst, cd_src, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Copy assignment failed"
            test_status = test_status + 1
        else
            ! Verify assignment worked
            if (cd_dst%path2project /= cd_src%path2project) then
                print *, "FAIL: Copy assignment didn't update path"
                test_status = test_status + 1
            end if
        end if
        
        ! Clean up
        call core_data_destroy(cd_src, ierr)
        call core_data_destroy(cd_dst, ierr)
        
        print *, "test_copy_assignment completed"
    end subroutine test_copy_assignment

end program test_core_deep_copy