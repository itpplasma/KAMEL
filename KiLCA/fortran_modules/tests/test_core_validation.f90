program test_core_validation
    use iso_fortran_env, only: real64, int32, int64
    use iso_c_binding
    use kilca_types_m
    use kilca_core_m
    use kilca_settings_m
    implicit none
    
    ! Test variables
    integer :: test_status = 0
    character(len=256) :: test_path = "/test/project/path/"
    
    ! Test 1: Validate empty core_data
    call test_validate_empty()
    
    ! Test 2: Validate initialized core_data
    call test_validate_initialized()
    
    ! Test 3: Validate with invalid path
    call test_validate_invalid_path()
    
    ! Test 4: Validate with settings
    call test_validate_with_settings()
    
    ! Test 5: Validate with modes
    call test_validate_with_modes()
    
    ! Test 6: Validate consistency checks
    call test_validate_consistency()
    
    ! Test 7: Validate full structure
    call test_validate_full_structure()
    
    ! Test 8: Get validation report
    call test_validation_report()
    
    if (test_status == 0) then
        print *, "All validation tests PASSED!"
    else
        print *, "Some tests FAILED. Status:", test_status
        stop 1
    end if
    
contains

    subroutine test_validate_empty()
        type(core_data_t), pointer :: cd
        logical :: is_valid
        character(len=1024) :: error_msg
        integer :: ierr
        
        print *, "Testing validation of empty core_data..."
        
        ! Create empty core_data
        call core_data_create(cd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not create core_data"
            test_status = test_status + 1
            return
        end if
        
        ! Validate
        call core_data_validate(cd, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Validation procedure failed"
            test_status = test_status + 1
        end if
        
        ! Empty but properly created core_data should be valid
        if (.not. is_valid) then
            print *, "FAIL: Empty core_data marked as invalid"
            print *, "Error: ", trim(error_msg)
            test_status = test_status + 1
        end if
        
        ! Clean up
        call core_data_destroy(cd, ierr)
        
        print *, "test_validate_empty completed"
    end subroutine test_validate_empty
    
    subroutine test_validate_initialized()
        type(core_data_t), pointer :: cd
        logical :: is_valid
        character(len=1024) :: error_msg
        integer :: ierr
        
        print *, "Testing validation of initialized core_data..."
        
        ! Create and initialize
        call core_data_create(cd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        cd%initialized = .true.
        
        ! Validate
        call core_data_validate(cd, is_valid, error_msg, ierr)
        
        if (.not. is_valid) then
            print *, "FAIL: Initialized core_data marked as invalid"
            test_status = test_status + 1
        end if
        
        ! Clean up
        call core_data_destroy(cd, ierr)
        
        print *, "test_validate_initialized completed"
    end subroutine test_validate_initialized
    
    subroutine test_validate_invalid_path()
        type(core_data_t), pointer :: cd
        logical :: is_valid
        character(len=1024) :: error_msg
        integer :: ierr
        
        print *, "Testing validation with invalid path..."
        
        ! Allocate manually to bypass constructor checks
        allocate(cd)
        cd%path2project = ""  ! Invalid empty path
        cd%sd => null()
        cd%bp => null()
        cd%dim = 0
        cd%initialized = .false.
        
        ! Validate
        call core_data_validate(cd, is_valid, error_msg, ierr)
        
        if (is_valid) then
            print *, "FAIL: Invalid path not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected invalid path: ", trim(error_msg)
        end if
        
        ! Clean up
        deallocate(cd)
        
        print *, "test_validate_invalid_path completed"
    end subroutine test_validate_invalid_path
    
    subroutine test_validate_with_settings()
        type(core_data_t), pointer :: cd
        type(settings_t), pointer :: sd
        logical :: is_valid
        character(len=1024) :: error_msg
        integer :: ierr
        
        print *, "Testing validation with settings..."
        
        ! Create core_data with settings
        call core_data_create(cd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        call settings_create(sd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            call core_data_destroy(cd, ierr)
            return
        end if
        
        cd%sd => sd
        cd%owns_settings = .true.
        
        ! Validate
        call core_data_validate(cd, is_valid, error_msg, ierr)
        
        if (.not. is_valid) then
            print *, "FAIL: Valid settings marked as invalid"
            print *, "Error: ", trim(error_msg)
            test_status = test_status + 1
        end if
        
        ! Test with inconsistent settings
        cd%dim = 5
        sd%antenna_settings%dma = 3  ! Different from dim
        
        call core_data_validate(cd, is_valid, error_msg, ierr)
        
        if (is_valid) then
            print *, "FAIL: Inconsistent antenna settings not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected inconsistent settings: ", trim(error_msg)
        end if
        
        ! Clean up
        call core_data_destroy(cd, ierr)
        
        print *, "test_validate_with_settings completed"
    end subroutine test_validate_with_settings
    
    subroutine test_validate_with_modes()
        type(core_data_t), pointer :: cd
        logical :: is_valid
        character(len=1024) :: error_msg
        integer :: ierr
        
        print *, "Testing validation with modes..."
        
        ! Create core_data
        call core_data_create(cd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Allocate modes with inconsistent dimension
        cd%dim = 5
        allocate(cd%mda(3))  ! Wrong size!
        
        ! Validate
        call core_data_validate(cd, is_valid, error_msg, ierr)
        
        if (is_valid) then
            print *, "FAIL: Mode array size mismatch not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected mode array mismatch: ", trim(error_msg)
        end if
        
        ! Fix and revalidate
        deallocate(cd%mda)
        allocate(cd%mda(cd%dim))
        
        call core_data_validate(cd, is_valid, error_msg, ierr)
        
        if (.not. is_valid) then
            print *, "FAIL: Fixed mode array still marked invalid"
            test_status = test_status + 1
        end if
        
        ! Clean up
        call core_data_destroy(cd, ierr)
        
        print *, "test_validate_with_modes completed"
    end subroutine test_validate_with_modes
    
    subroutine test_validate_consistency()
        type(core_data_t), pointer :: cd
        logical :: is_consistent
        character(len=1024) :: report
        integer :: ierr
        
        print *, "Testing consistency validation..."
        
        ! Create core_data
        call core_data_create(cd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Check consistency
        call core_data_check_consistency(cd, is_consistent, report, ierr)
        
        if (.not. is_consistent) then
            print *, "WARN: Empty core_data marked as inconsistent"
            print *, "Report: ", trim(report)
        end if
        
        ! Create inconsistent state
        cd%owns_settings = .true.  ! But sd is null!
        
        call core_data_check_consistency(cd, is_consistent, report, ierr)
        
        if (is_consistent) then
            print *, "FAIL: Inconsistent ownership not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected inconsistency: ", trim(report)
        end if
        
        ! Clean up
        call core_data_destroy(cd, ierr)
        
        print *, "test_validate_consistency completed"
    end subroutine test_validate_consistency
    
    subroutine test_validate_full_structure()
        type(core_data_t), pointer :: cd
        type(settings_t), pointer :: sd
        logical :: is_valid
        character(len=1024) :: error_msg
        integer :: ierr
        
        print *, "Testing validation of full structure..."
        
        ! Create fully populated structure
        call core_data_create(cd, "/test/vacuum/project/", ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Initialize mode independent data
        call calc_and_set_mode_independent_core_data(cd, ierr)
        
        ! Create valid settings
        call settings_create(sd, cd%path2project, ierr)
        if (ierr == KILCA_SUCCESS) then
            ! Set valid antenna parameters
            sd%antenna_settings%ra = 10.0_dp
            sd%antenna_settings%wa = 1.0_dp
            sd%antenna_settings%I0 = 1000.0_dp
            sd%antenna_settings%dma = 2
            allocate(sd%antenna_settings%modes(4))
            sd%antenna_settings%modes = [1, 1, 2, 2]
            sd%antenna_settings%flab = cmplx(50.0e6_dp, 0.0_dp, dp)
            
            ! Set valid background parameters
            sd%background_settings%rtor = 625.0_dp
            sd%background_settings%rp = 99.0_dp
            sd%background_settings%B0 = 20000.0_dp
            
            cd%sd => sd
            cd%owns_settings = .true.
            
            ! Initialize modes
            call calc_and_set_mode_dependent_core_data_antenna(cd, ierr)
        end if
        
        ! Validate full structure
        call core_data_validate_full(cd, is_valid, error_msg, ierr)
        
        if (.not. is_valid) then
            print *, "FAIL: Valid full structure marked as invalid"
            print *, "Error: ", trim(error_msg)
            test_status = test_status + 1
        end if
        
        ! Clean up
        call core_data_destroy(cd, ierr)
        
        print *, "test_validate_full_structure completed"
    end subroutine test_validate_full_structure
    
    subroutine test_validation_report()
        type(core_data_t), pointer :: cd
        character(len=4096) :: report
        integer :: ierr
        
        print *, "Testing validation report generation..."
        
        ! Create core_data
        call core_data_create(cd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Get validation report
        call core_data_get_validation_report(cd, report, ierr)
        
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not generate validation report"
            test_status = test_status + 1
        else
            print *, "Validation report generated successfully"
            if (len_trim(report) == 0) then
                print *, "FAIL: Empty validation report"
                test_status = test_status + 1
            end if
        end if
        
        ! Clean up
        call core_data_destroy(cd, ierr)
        
        print *, "test_validation_report completed"
    end subroutine test_validation_report

end program test_core_validation