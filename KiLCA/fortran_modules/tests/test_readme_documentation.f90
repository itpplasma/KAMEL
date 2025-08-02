program test_readme_documentation
    implicit none
    
    integer :: test_status = 0
    
    print *, "===========================================" 
    print *, "Testing README Documentation [RED PHASE - SHOULD FAIL]"
    print *, "==========================================="
    
    call test_readme_fortran_documentation()
    
    if (test_status == 0) then
        print *, ""
        print *, "README documentation test PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "README documentation test FAILED:", test_status, "requirement(s) not met (expected in RED phase)"
        stop 1
    end if

contains

    !> Test for README documentation of Fortran main program
    subroutine test_readme_fortran_documentation()
        integer :: unit, iostat
        character(len=1024) :: line
        character(len=*), parameter :: readme_file = "../README.md"
        logical :: file_exists
        logical :: has_fortran_section = .false.
        logical :: has_usage_section = .false.
        logical :: has_build_instructions = .false.
        logical :: has_kilca_main_mention = .false.
        logical :: has_command_examples = .false.
        logical :: has_requirements = .false.
        
        print *, "Testing README.md for Fortran program documentation..."
        print *, ""
        
        ! Check if README exists
        inquire(file=readme_file, exist=file_exists)
        if (.not. file_exists) then
            print *, "FAIL: README.md not found in fortran_modules directory"
            test_status = test_status + 1
            return
        end if
        
        print *, "PASS: README.md exists"
        
        ! Parse README for required content
        open(newunit=unit, file=readme_file, status="old", iostat=iostat)
        if (iostat /= 0) then
            print *, "FAIL: Could not read README.md"
            test_status = test_status + 1
            return
        end if
        
        ! Check for required sections and content
        do
            read(unit, '(a)', iostat=iostat) line
            if (iostat /= 0) exit
            
            ! Convert to lowercase for case-insensitive search
            call to_lower(line)
            
            ! Check for Fortran main program section
            if (index(line, "fortran main program") > 0 .or. &
                index(line, "kilca fortran") > 0 .or. &
                index(line, "fortran implementation") > 0) then
                has_fortran_section = .true.
            end if
            
            ! Check for usage section
            if (index(line, "## usage") > 0 .or. &
                index(line, "### usage") > 0 .or. &
                index(line, "how to use") > 0) then
                has_usage_section = .true.
            end if
            
            ! Check for build instructions
            if (index(line, "build") > 0 .and. &
                (index(line, "instruction") > 0 .or. index(line, "compile") > 0)) then
                has_build_instructions = .true.
            end if
            
            ! Check for kilca_main executable mention
            if (index(line, "kilca_main") > 0 .or. &
                index(line, "./kilca_main") > 0) then
                has_kilca_main_mention = .true.
            end if
            
            ! Check for command examples
            if (index(line, "./kilca_main") > 0 .and. &
                (index(line, "path") > 0 .or. index(line, "project") > 0)) then
                has_command_examples = .true.
            end if
            
            ! Check for requirements/dependencies
            if (index(line, "requirement") > 0 .or. &
                index(line, "dependenc") > 0 .or. &
                index(line, "prerequisite") > 0) then
                has_requirements = .true.
            end if
        end do
        
        close(unit)
        
        ! Report findings
        print *, ""
        print *, "=== README CONTENT ANALYSIS ==="
        
        if (.not. has_fortran_section) then
            print *, "FAIL: No Fortran main program section found"
            test_status = test_status + 1
        else
            print *, "PASS: Fortran main program section found"
        end if
        
        if (.not. has_usage_section) then
            print *, "FAIL: No usage section found"
            test_status = test_status + 1
        else
            print *, "PASS: Usage section found"
        end if
        
        if (.not. has_build_instructions) then
            print *, "FAIL: No build instructions found"
            test_status = test_status + 1
        else
            print *, "PASS: Build instructions found"
        end if
        
        if (.not. has_kilca_main_mention) then
            print *, "FAIL: No mention of kilca_main executable"
            test_status = test_status + 1
        else
            print *, "PASS: kilca_main executable mentioned"
        end if
        
        if (.not. has_command_examples) then
            print *, "FAIL: No command line examples found"
            test_status = test_status + 1
        else
            print *, "PASS: Command line examples found"
        end if
        
        if (.not. has_requirements) then
            print *, "FAIL: No requirements/dependencies section found"
            test_status = test_status + 1
        else
            print *, "PASS: Requirements/dependencies section found"
        end if
        
        print *, ""
        print *, "Total documentation requirements not met:", test_status
        
    end subroutine test_readme_fortran_documentation
    
    !> Convert string to lowercase
    subroutine to_lower(str)
        character(len=*), intent(inout) :: str
        integer :: i, ich
        
        do i = 1, len_trim(str)
            ich = iachar(str(i:i))
            if (ich >= 65 .and. ich <= 90) then
                str(i:i) = achar(ich + 32)
            end if
        end do
    end subroutine to_lower

end program test_readme_documentation