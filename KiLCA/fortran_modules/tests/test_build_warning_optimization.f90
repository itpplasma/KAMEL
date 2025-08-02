program test_build_warning_optimization
    implicit none
    
    integer :: test_status = 0
    
    print *, "===========================================" 
    print *, "Testing Build Warning Optimization [MAINTENANCE - SHOULD PASS]"
    print *, "==========================================="
    
    call test_for_critical_build_warnings()
    
    if (test_status == 0) then
        print *, ""
        print *, "Build warning optimization test PASSED - critical warnings resolved"
        stop 0
    else
        print *, ""
        print *, "Build warning optimization test FAILED:", test_status, "critical warning(s) detected"
        stop 1
    end if

contains

    !> Test for critical build warnings that should be optimized
    subroutine test_for_critical_build_warnings()
        integer :: unit, iostat
        character(len=1024) :: line
        character(len=*), parameter :: build_log = "build_warnings_detailed.log"
        logical :: file_exists
        integer :: unused_warnings, line_truncation_warnings, format_warnings, intent_warnings
        
        print *, "Testing for critical build warnings..."
        print *, ""
        
        ! Initialize counters
        unused_warnings = 0
        line_truncation_warnings = 0
        format_warnings = 0
        intent_warnings = 0
        
        ! Build the project and capture warnings
        print *, "Building project to capture detailed warnings..."
        call system("cd .. && make clean > /dev/null 2>&1")
        call system("cd .. && make ninja > build_warnings_detailed.log 2>&1")
        
        ! Check if build log exists
        inquire(file="../" // build_log, exist=file_exists)
        if (.not. file_exists) then
            print *, "FAIL: Build log file not found"
            test_status = test_status + 1
            return
        end if
        
        ! Parse build log for different types of warnings
        open(newunit=unit, file="../" // build_log, status="old", iostat=iostat)
        if (iostat /= 0) then
            print *, "FAIL: Could not read build log"
            test_status = test_status + 1
            return
        end if
        
        ! Count different types of warnings
        do
            read(unit, '(a)', iostat=iostat) line
            if (iostat /= 0) exit
            
            ! Count unused variable/argument warnings
            if (index(line, "Unused dummy argument") > 0 .or. &
                index(line, "unused-dummy-argument") > 0 .or. &
                index(line, "Unused variable") > 0) then
                unused_warnings = unused_warnings + 1
            end if
            
            ! Count line truncation warnings (critical for readability)
            if (index(line, "Line truncated") > 0 .or. &
                index(line, "line-truncation") > 0) then
                line_truncation_warnings = line_truncation_warnings + 1
                print *, "CRITICAL LINE TRUNCATION:", trim(line)
            end if
            
            ! Count format warnings
            if (index(line, "format specifies type") > 0 .or. &
                index(line, "Wformat") > 0) then
                format_warnings = format_warnings + 1
                print *, "FORMAT WARNING:", trim(line)
            end if
            
            ! Count intent warnings (parameters not set)
            if (index(line, "INTENT(OUT)") > 0 .and. index(line, "was not set") > 0) then
                intent_warnings = intent_warnings + 1
                print *, "INTENT WARNING:", trim(line)
            end if
        end do
        
        close(unit)
        
        ! Report warning counts
        print *, ""
        print *, "=== BUILD WARNING SUMMARY ==="
        print *, "Unused variable/argument warnings:", unused_warnings
        print *, "Line truncation warnings:", line_truncation_warnings  
        print *, "Format warnings:", format_warnings
        print *, "Intent(OUT) not set warnings:", intent_warnings
        print *, ""
        
        ! Set failure conditions for critical warnings
        if (line_truncation_warnings > 0) then
            print *, "CRITICAL: Line truncation warnings must be fixed (readability)"
            test_status = test_status + line_truncation_warnings
        end if
        
        if (format_warnings > 0) then
            print *, "CRITICAL: Format warnings must be fixed (potential runtime errors)"
            test_status = test_status + format_warnings
        end if
        
        if (intent_warnings > 5) then
            print *, "CRITICAL: Too many intent(OUT) warnings (>5) - likely interface issues"
            test_status = test_status + 1
        end if
        
        if (unused_warnings > 200) then
            print *, "CRITICAL: Excessive unused warnings (>200) - code quality issue"
            test_status = test_status + 1
        end if
        
        ! Quality thresholds
        print *, "=== QUALITY ASSESSMENT ==="
        if (unused_warnings < 50) then
            print *, "GOOD: Unused warnings are manageable (<50)"
        else if (unused_warnings < 150) then
            print *, "ACCEPTABLE: Unused warnings are moderate (50-150)"
        else
            print *, "POOR: Unused warnings are excessive (>150)"
        end if
        
        print *, ""
        print *, "Total critical issues detected:", test_status
        
        ! Clean up build log
        call system("rm -f ../" // build_log)
        
    end subroutine test_for_critical_build_warnings

end program test_build_warning_optimization