program test_unused_variable_warnings
    implicit none
    
    integer :: test_status = 0
    
    print *, "===========================================" 
    print *, "Testing Unused Variable Warnings [RED PHASE - SHOULD FAIL]"
    print *, "==========================================="
    
    call test_for_unused_variable_warnings()
    
    if (test_status == 0) then
        print *, ""
        print *, "Unused variable warnings test PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "Unused variable warnings test FAILED:", test_status, "warning(s) detected (expected in RED phase)"
        stop 1
    end if

contains

    !> Test for unused variable warnings
    subroutine test_for_unused_variable_warnings()
        integer :: unit, iostat
        character(len=1024) :: line
        character(len=*), parameter :: build_log = "build_warnings.log"
        logical :: file_exists
        
        print *, "Testing for unused variable warnings in build output..."
        print *, ""
        
        ! Build the project and capture warnings
        print *, "Building project to capture warnings..."
        call system("cd .. && make clean > /dev/null 2>&1")
        call system("cd .. && make all > build_warnings.log 2>&1")
        
        ! Check if build log exists
        inquire(file="../" // build_log, exist=file_exists)
        if (.not. file_exists) then
            print *, "FAIL: Build log file not found"
            test_status = test_status + 1
            return
        end if
        
        ! Parse build log for unused variable warnings
        open(newunit=unit, file="../" // build_log, status="old", iostat=iostat)
        if (iostat /= 0) then
            print *, "FAIL: Could not read build log"
            test_status = test_status + 1
            return
        end if
        
        ! Count unused variable warnings
        do
            read(unit, '(a)', iostat=iostat) line
            if (iostat /= 0) exit
            
            if (index(line, "Unused dummy argument") > 0 .or. &
                index(line, "unused-dummy-argument") > 0 .or. &
                index(line, "Unused variable") > 0) then
                
                print *, "FOUND UNUSED WARNING:", trim(line)
                test_status = test_status + 1
            end if
        end do
        
        close(unit)
        
        if (test_status > 0) then
            print *, ""
            print *, "Found", test_status, "unused variable/argument warnings"
        else
            print *, "No unused variable warnings found"
        end if
        
        ! Clean up build log
        call system("rm -f ../" // build_log)
        
    end subroutine test_for_unused_variable_warnings

end program test_unused_variable_warnings