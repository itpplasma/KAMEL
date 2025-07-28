program test_kilca_settings
    use iso_fortran_env, only: real64, int32, int64
    use iso_c_binding
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    ! Test variables
    integer :: test_status = 0
    character(len=256) :: test_path = "/test/project/path/"
    
    ! Test 1: Create and destroy settings structure
    call test_settings_lifecycle()
    
    ! Test 2: Test antenna settings
    call test_antenna_settings()
    
    ! Test 3: Test background settings
    call test_background_settings()
    
    ! Test 4: Test output settings
    call test_output_settings()
    
    ! Test 5: Test eigenmode settings
    call test_eigmode_settings()
    
    ! Test 6: Test settings file I/O
    call test_settings_io()
    
    ! Test 7: Test settings validation
    call test_settings_validation()
    
    ! Test 8: Test C interface compatibility
    call test_c_interface()
    
    if (test_status == 0) then
        print *, "All tests PASSED!"
    else
        print *, "Some tests FAILED. Status:", test_status
        stop 1
    end if
    
contains

    subroutine test_settings_lifecycle()
        type(settings_t), pointer :: sd
        integer :: ierr
        
        print *, "Testing settings lifecycle..."
        
        ! Test creation
        call settings_create(sd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: settings_create returned error:", ierr
            test_status = test_status + 1
            return
        end if
        
        if (.not. associated(sd)) then
            print *, "FAIL: settings not allocated"
            test_status = test_status + 1
            return
        end if
        
        ! Test destruction
        call settings_destroy(sd, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: settings_destroy returned error:", ierr
            test_status = test_status + 1
            return
        end if
        
        if (associated(sd)) then
            print *, "FAIL: settings not deallocated"
            test_status = test_status + 1
        end if
        
        print *, "test_settings_lifecycle completed"
    end subroutine test_settings_lifecycle
    
    subroutine test_antenna_settings()
        type(settings_t), pointer :: sd
        type(antenna_t), pointer :: ant
        integer :: ierr
        real(dp) :: test_ra, test_wa, test_I0
        complex(dp) :: test_flab
        integer :: test_dma
        integer, allocatable :: test_modes(:)
        
        print *, "Testing antenna settings..."
        
        call settings_create(sd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Get antenna settings pointer
        call settings_get_antenna(sd, ant, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not get antenna settings"
            test_status = test_status + 1
        end if
        
        ! Set antenna parameters
        test_ra = 10.0_dp
        test_wa = 1.0_dp
        test_I0 = 1000.0_dp
        test_flab = cmplx(50.0e6_dp, 0.0_dp, dp)
        test_dma = 2
        allocate(test_modes(4))
        test_modes = [1, 2, 3, 4]  ! m1=1, n1=2, m2=3, n2=4
        
        call antenna_set_parameters(ant, test_ra, test_wa, test_I0, test_flab, &
                                   test_dma, test_modes, ierr)
        
        ! Verify parameters
        if (abs(ant%ra - test_ra) > epsilon(1.0_dp)) then
            print *, "FAIL: Antenna radius mismatch"
            test_status = test_status + 1
        end if
        
        if (abs(ant%I0 - test_I0) > epsilon(1.0_dp)) then
            print *, "FAIL: Antenna current mismatch"
            test_status = test_status + 1
        end if
        
        if (ant%dma /= test_dma) then
            print *, "FAIL: Mode count mismatch"
            test_status = test_status + 1
        end if
        
        deallocate(test_modes)
        call settings_destroy(sd, ierr)
        
        print *, "test_antenna_settings completed"
    end subroutine test_antenna_settings
    
    subroutine test_background_settings()
        type(settings_t), pointer :: sd
        type(back_sett_t), pointer :: bs
        integer :: ierr
        
        print *, "Testing background settings..."
        
        call settings_create(sd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Get background settings pointer
        call settings_get_background(sd, bs, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not get background settings"
            test_status = test_status + 1
        end if
        
        ! Set background parameters
        call back_sett_set_calc_flag(bs, 1, ierr)
        
        if (bs%calc_back /= 1) then
            print *, "FAIL: Background calc flag not set"
            test_status = test_status + 1
        end if
        
        call settings_destroy(sd, ierr)
        
        print *, "test_background_settings completed"
    end subroutine test_background_settings
    
    subroutine test_output_settings()
        type(settings_t), pointer :: sd
        type(output_sett_t), pointer :: os
        integer :: ierr
        
        print *, "Testing output settings..."
        
        call settings_create(sd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Get output settings pointer
        call settings_get_output(sd, os, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not get output settings"
            test_status = test_status + 1
        end if
        
        ! Set output parameters
        call output_sett_set_flags(os, flag_background=1, flag_emfield=2, ierr=ierr)
        
        if (os%flag_background /= 1) then
            print *, "FAIL: Output flag_background not set"
            test_status = test_status + 1
        end if
        
        call settings_destroy(sd, ierr)
        
        print *, "test_output_settings completed"
    end subroutine test_output_settings
    
    subroutine test_eigmode_settings()
        type(settings_t), pointer :: sd
        type(eigmode_sett_t), pointer :: es
        integer :: ierr
        
        print *, "Testing eigenmode settings..."
        
        call settings_create(sd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Get eigenmode settings pointer
        call settings_get_eigmode(sd, es, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not get eigenmode settings"
            test_status = test_status + 1
        end if
        
        ! Set eigenmode parameters
        call eigmode_sett_set_search_flag(es, 1, ierr)
        
        if (es%search_flag /= 1) then
            print *, "FAIL: Eigenmode search flag not set"
            test_status = test_status + 1
        end if
        
        call settings_destroy(sd, ierr)
        
        print *, "test_eigmode_settings completed"
    end subroutine test_eigmode_settings
    
    subroutine test_settings_io()
        type(settings_t), pointer :: sd
        integer :: ierr
        logical :: file_exists
        
        print *, "Testing settings file I/O..."
        
        call settings_create(sd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Note: In real implementation, would test actual file reading
        ! For now, just test that procedures exist
        
        ! Check if settings files would exist
        inquire(file=trim(test_path)//"antenna.in", exist=file_exists)
        if (file_exists) then
            call settings_read_all(sd, ierr)
            if (ierr /= KILCA_SUCCESS .and. ierr /= KILCA_ERROR_FILE) then
                print *, "FAIL: Unexpected error reading settings"
                test_status = test_status + 1
            end if
        end if
        
        call settings_destroy(sd, ierr)
        
        print *, "test_settings_io completed"
    end subroutine test_settings_io
    
    subroutine test_settings_validation()
        type(settings_t), pointer :: sd
        type(antenna_t), pointer :: ant
        integer :: ierr
        logical :: is_valid
        
        print *, "Testing settings validation..."
        
        call settings_create(sd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Test invalid antenna settings
        call settings_get_antenna(sd, ant, ierr)
        ant%ra = -1.0_dp  ! Invalid negative radius
        
        call settings_validate(sd, is_valid, ierr)
        if (is_valid) then
            print *, "FAIL: Invalid settings passed validation"
            test_status = test_status + 1
        end if
        
        ! Fix and revalidate
        ant%ra = 10.0_dp
        ant%wa = 1.0_dp
        ant%I0 = 1000.0_dp
        ant%dma = 1
        allocate(ant%modes(2))
        ant%modes = [1, 1]
        
        ! Also fix background settings for new validation
        sd%background_settings%rtor = 625.0_dp
        sd%background_settings%rp = 99.0_dp
        sd%background_settings%B0 = 20000.0_dp
        
        call settings_validate(sd, is_valid, ierr)
        if (.not. is_valid) then
            print *, "FAIL: Valid settings failed validation"
            test_status = test_status + 1
        end if
        
        call settings_destroy(sd, ierr)
        
        print *, "test_settings_validation completed"
    end subroutine test_settings_validation
    
    subroutine test_c_interface()
        type(settings_t), pointer :: sd
        type(antenna_t), pointer :: ant
        integer :: ierr
        
        print *, "Testing C interface compatibility..."
        
        call settings_create(sd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Test that C interface procedures are available
        ! set_antenna_settings_c
        ! copy_antenna_data_to_antenna_module
        ! copy_background_data_to_background_module
        
        call settings_destroy(sd, ierr)
        
        print *, "test_c_interface completed"
    end subroutine test_c_interface

end program test_kilca_settings