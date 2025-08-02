program test_disabled_files_cleanup
    implicit none
    
    integer :: test_status = 0
    logical :: file_exists
    character(len=256) :: disabled_files(9)
    integer :: i
    
    print *, "===========================================" 
    print *, "Testing Disabled Files Cleanup [MAINTENANCE - SHOULD PASS]"
    print *, "==========================================="
    
    ! List of disabled test files that should be removed
    disabled_files(1) = "test_kilca_bdf_solver.f90.disabled"
    disabled_files(2) = "test_kilca_integration.f90.disabled"
    disabled_files(3) = "test_kilca_integration_simple.f90.disabled"
    disabled_files(4) = "test_kilca_mode_solver.f90.disabled"
    disabled_files(5) = "test_kilca_mode_solver_simple.f90.disabled"
    disabled_files(6) = "test_kilca_solver_advanced.f90.disabled"
    disabled_files(7) = "test_kilca_validation.f90.disabled"
    disabled_files(8) = "test_kilca_zone.f90.disabled"
    disabled_files(9) = "test_memory_management.f90.disabled"
    
    print *, "Checking for presence of disabled test files..."
    print *, ""
    
    ! Check each disabled file
    do i = 1, 9
        inquire(file=trim(disabled_files(i)), exist=file_exists)
        if (file_exists) then
            print *, "FAIL: Disabled file still exists:", trim(disabled_files(i))
            test_status = test_status + 1
        else
            print *, "PASS: Disabled file properly removed:", trim(disabled_files(i))
        end if
    end do
    
    print *, ""
    if (test_status == 0) then
        print *, "Disabled files cleanup test PASSED - repository is clean"
        stop 0
    else
        print *, "Disabled files cleanup test FAILED:", test_status, "disabled file(s) still present"
        stop 1
    end if
    
end program test_disabled_files_cleanup