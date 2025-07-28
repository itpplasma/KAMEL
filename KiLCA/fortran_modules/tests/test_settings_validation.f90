program test_settings_validation
    use iso_fortran_env, only: real64, int32, int64
    use iso_c_binding
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    ! Test variables
    integer :: test_status = 0
    
    ! Test 1: Validate antenna settings
    call test_antenna_validation()
    
    ! Test 2: Validate background settings
    call test_background_validation()
    
    ! Test 3: Validate output settings
    call test_output_validation()
    
    ! Test 4: Validate eigenmode settings
    call test_eigenmode_validation()
    
    ! Test 5: Validate complete settings
    call test_settings_complete_validation()
    
    ! Test 6: Test validation error messages
    call test_validation_error_messages()
    
    ! Test 7: Test validation ranges and boundaries
    call test_validation_boundaries()
    
    ! Test 8: Test consistency validations
    call test_consistency_validations()
    
    if (test_status == 0) then
        print *, "All settings validation tests PASSED!"
    else
        print *, "Some tests FAILED. Status:", test_status
        stop 1
    end if
    
contains

    subroutine test_antenna_validation()
        type(antenna_t) :: ant
        logical :: is_valid
        character(len=1024) :: error_msg
        integer :: ierr
        
        print *, "Testing antenna validation..."
        
        ! Test 1: Valid antenna settings
        ant%ra = 10.5_dp
        ant%wa = 1.2_dp
        ant%I0 = 1000.0_dp
        ant%flab = cmplx(5.0e7_dp, 0.0_dp, dp)
        ant%dma = 2
        allocate(ant%modes(4))
        ant%modes = [1, 1, 2, 2]
        ant%flag_debug = 0
        ant%flag_eigmode = 0
        
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. .not. is_valid) then
            print *, "FAIL: Valid antenna marked as invalid"
            print *, "Error:", trim(error_msg)
            test_status = test_status + 1
        end if
        
        ! Test 2: Invalid ra (negative)
        ant%ra = -1.0_dp
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Negative ra not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected negative ra:", trim(error_msg)
        end if
        ant%ra = 10.5_dp  ! Reset
        
        ! Test 3: Invalid wa (zero)
        ant%wa = 0.0_dp
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Zero wa not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected zero wa:", trim(error_msg)
        end if
        ant%wa = 1.2_dp  ! Reset
        
        ! Test 4: Invalid I0 (negative)
        ant%I0 = -100.0_dp
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Negative I0 not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected negative I0:", trim(error_msg)
        end if
        ant%I0 = 1000.0_dp  ! Reset
        
        ! Test 5: Invalid dma (zero)
        ant%dma = 0
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Zero dma not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected zero dma:", trim(error_msg)
        end if
        ant%dma = 2  ! Reset
        
        ! Test 6: Modes array not allocated
        deallocate(ant%modes)
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Unallocated modes not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected unallocated modes:", trim(error_msg)
        end if
        
        ! Test 7: Modes array size mismatch
        allocate(ant%modes(2))  ! Should be 2*dma = 4
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Mode array size mismatch not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected mode array size mismatch:", trim(error_msg)
        end if
        
        ! Clean up
        if (allocated(ant%modes)) deallocate(ant%modes)
        
        print *, "test_antenna_validation completed"
    end subroutine test_antenna_validation
    
    subroutine test_background_validation()
        type(back_sett_t) :: bs
        logical :: is_valid
        character(len=1024) :: error_msg
        integer :: ierr
        
        print *, "Testing background validation..."
        
        ! Test 1: Valid background settings
        bs%rtor = 625.0_dp
        bs%rp = 99.0_dp
        bs%B0 = 20000.0_dp
        bs%calc_back = 1
        bs%flag_back = "normal"
        bs%N = 5
        bs%V_gal_sys = 0.0_dp
        bs%V_scale = 1.0_dp
        bs%m_i = 1.0_dp
        bs%zele = 1.0_dp
        bs%zion = 1.0_dp
        bs%flag_debug = 0
        allocate(bs%mass(2))
        allocate(bs%charge(2))
        bs%mass = [1.67e-24_dp, 9.11e-28_dp]
        bs%charge = [1.602e-19_dp, -1.602e-19_dp]
        bs%huge_factor = 1.0e30_dp
        
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. .not. is_valid) then
            print *, "FAIL: Valid background marked as invalid"
            print *, "Error:", trim(error_msg)
            test_status = test_status + 1
        end if
        
        ! Test 2: Invalid rtor (negative)
        bs%rtor = -100.0_dp
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Negative rtor not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected negative rtor:", trim(error_msg)
        end if
        bs%rtor = 625.0_dp  ! Reset
        
        ! Test 3: Invalid rp (zero)
        bs%rp = 0.0_dp
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Zero rp not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected zero rp:", trim(error_msg)
        end if
        bs%rp = 99.0_dp  ! Reset
        
        ! Test 4: Invalid B0 (negative)
        bs%B0 = -1000.0_dp
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Negative B0 not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected negative B0:", trim(error_msg)
        end if
        bs%B0 = 20000.0_dp  ! Reset
        
        ! Test 5: Invalid calc_back (zero)
        bs%calc_back = 0
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Zero calc_back not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected zero calc_back:", trim(error_msg)
        end if
        bs%calc_back = 1  ! Reset
        
        ! Test 6: Invalid N (even number)
        bs%N = 4  ! Must be odd
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Even N not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected even N:", trim(error_msg)
        end if
        bs%N = 5  ! Reset
        
        ! Test 7: Invalid m_i (negative)
        bs%m_i = -1.0_dp
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Negative m_i not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected negative m_i:", trim(error_msg)
        end if
        bs%m_i = 1.0_dp  ! Reset
        
        ! Clean up
        if (allocated(bs%mass)) deallocate(bs%mass)
        if (allocated(bs%charge)) deallocate(bs%charge)
        if (allocated(bs%flag_back)) deallocate(bs%flag_back)
        if (allocated(bs%path2profiles)) deallocate(bs%path2profiles)
        
        print *, "test_background_validation completed"
    end subroutine test_background_validation
    
    subroutine test_output_validation()
        type(output_sett_t) :: os
        logical :: is_valid
        character(len=1024) :: error_msg
        integer :: ierr
        
        print *, "Testing output validation..."
        
        ! Test 1: Valid output settings
        os%flag_background = 1
        os%flag_emfield = 1
        os%flag_additional = 0
        os%flag_dispersion = 0
        os%num_quants = 3
        allocate(os%flag_quants(3))
        os%flag_quants = [1, 0, 1]
        os%flag_debug = 0
        
        call output_sett_validate(os, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. .not. is_valid) then
            print *, "FAIL: Valid output marked as invalid"
            print *, "Error:", trim(error_msg)
            test_status = test_status + 1
        end if
        
        ! Test 2: Invalid flag values (out of range)
        os%flag_background = 3  ! Should be 0, 1, or 2
        call output_sett_validate(os, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Invalid flag_background not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected invalid flag_background:", trim(error_msg)
        end if
        os%flag_background = 1  ! Reset
        
        ! Test 3: Negative num_quants
        os%num_quants = -1
        call output_sett_validate(os, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Negative num_quants not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected negative num_quants:", trim(error_msg)
        end if
        os%num_quants = 3  ! Reset
        
        ! Test 4: flag_quants array size mismatch
        deallocate(os%flag_quants)
        allocate(os%flag_quants(2))  ! Should be num_quants = 3
        call output_sett_validate(os, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: flag_quants size mismatch not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected flag_quants size mismatch:", trim(error_msg)
        end if
        
        ! Clean up
        if (allocated(os%flag_quants)) deallocate(os%flag_quants)
        
        print *, "test_output_validation completed"
    end subroutine test_output_validation
    
    subroutine test_eigenmode_validation()
        type(eigmode_sett_t) :: es
        logical :: is_valid
        character(len=1024) :: error_msg
        integer :: ierr
        
        print *, "Testing eigenmode validation..."
        
        ! Test 1: Valid eigenmode settings
        es%fname = "output.dat"
        es%search_flag = 1
        es%rdim = 100
        es%rfmin = 0.0_dp
        es%rfmax = 1.0e9_dp
        es%idim = 100
        es%ifmin = -1.0e6_dp
        es%ifmax = 1.0e6_dp
        es%stop_flag = 0
        es%eps_res = 1.0e-6_dp
        es%eps_abs = 1.0e-8_dp
        es%eps_rel = 1.0e-6_dp
        es%delta = 1.0e-6_dp
        es%test_roots = 0
        es%flag_debug = 0
        es%Nguess = 2
        es%kmin = 1
        es%kmax = 10
        allocate(es%fstart(2))
        es%fstart = [cmplx(1.0e6_dp, 0.0_dp, dp), cmplx(2.0e6_dp, 0.0_dp, dp)]
        es%n_zeros = 10
        es%use_winding = 0
        
        call eigmode_sett_validate(es, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. .not. is_valid) then
            print *, "FAIL: Valid eigenmode marked as invalid"
            print *, "Error:", trim(error_msg)
            test_status = test_status + 1
        end if
        
        ! Test 2: Invalid dimensions (negative)
        es%rdim = -10
        call eigmode_sett_validate(es, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Negative rdim not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected negative rdim:", trim(error_msg)
        end if
        es%rdim = 100  ! Reset
        
        ! Test 3: Invalid frequency range (rfmin >= rfmax)
        es%rfmin = 1.0e9_dp
        es%rfmax = 0.0_dp
        call eigmode_sett_validate(es, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Invalid frequency range not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected invalid frequency range:", trim(error_msg)
        end if
        es%rfmin = 0.0_dp
        es%rfmax = 1.0e9_dp  ! Reset
        
        ! Test 4: Invalid tolerance (negative eps_res)
        es%eps_res = -1.0e-6_dp
        call eigmode_sett_validate(es, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Negative eps_res not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected negative eps_res:", trim(error_msg)
        end if
        es%eps_res = 1.0e-6_dp  ! Reset
        
        ! Test 5: fstart array size mismatch
        deallocate(es%fstart)
        allocate(es%fstart(1))  ! Should be Nguess = 2
        call eigmode_sett_validate(es, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: fstart size mismatch not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected fstart size mismatch:", trim(error_msg)
        end if
        
        ! Clean up
        if (allocated(es%fstart)) deallocate(es%fstart)
        if (allocated(es%fname)) deallocate(es%fname)
        
        print *, "test_eigenmode_validation completed"
    end subroutine test_eigenmode_validation
    
    subroutine test_settings_complete_validation()
        type(settings_t), pointer :: sd
        logical :: is_valid
        character(len=1024) :: error_msg
        integer :: ierr
        
        print *, "Testing complete settings validation..."
        
        call settings_create(sd, "/test/path/", ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not create settings"
            test_status = test_status + 1
            return
        end if
        
        ! Set up valid settings
        sd%antenna_settings%ra = 10.5_dp
        sd%antenna_settings%wa = 1.2_dp
        sd%antenna_settings%I0 = 1000.0_dp
        sd%antenna_settings%flab = cmplx(5.0e7_dp, 0.0_dp, dp)
        sd%antenna_settings%dma = 2
        allocate(sd%antenna_settings%modes(4))
        sd%antenna_settings%modes = [1, 1, 2, 2]
        
        sd%background_settings%rtor = 625.0_dp
        sd%background_settings%rp = 99.0_dp
        sd%background_settings%B0 = 20000.0_dp
        sd%background_settings%calc_back = 1
        sd%background_settings%N = 5
        sd%background_settings%m_i = 1.0_dp
        
        call settings_validate_complete(sd, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. .not. is_valid) then
            print *, "FAIL: Valid complete settings marked as invalid"
            print *, "Error:", trim(error_msg)
            test_status = test_status + 1
        end if
        
        ! Test with invalid antenna
        sd%antenna_settings%ra = -1.0_dp
        call settings_validate_complete(sd, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Invalid antenna in complete settings not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected invalid antenna in complete settings"
        end if
        
        call settings_destroy(sd, ierr)
        
        print *, "test_settings_complete_validation completed"
    end subroutine test_settings_complete_validation
    
    subroutine test_validation_error_messages()
        type(antenna_t) :: ant
        logical :: is_valid
        character(len=1024) :: error_msg
        integer :: ierr
        
        print *, "Testing validation error messages..."
        
        ! Test specific error message content
        ant%ra = -5.0_dp
        ant%wa = 1.0_dp
        ant%I0 = 100.0_dp
        ant%dma = 1
        
        call antenna_validate(ant, is_valid, error_msg, ierr)
        if (ierr == KILCA_SUCCESS .and. .not. is_valid) then
            if (index(error_msg, "ra") > 0 .and. index(error_msg, "positive") > 0) then
                print *, "Correctly generated specific error message for ra"
            else
                print *, "FAIL: Error message not specific enough:", trim(error_msg)
                test_status = test_status + 1
            end if
        else
            print *, "FAIL: Expected validation error not generated"
            test_status = test_status + 1
        end if
        
        print *, "test_validation_error_messages completed"
    end subroutine test_validation_error_messages
    
    subroutine test_validation_boundaries()
        type(back_sett_t) :: bs
        logical :: is_valid
        character(len=1024) :: error_msg
        integer :: ierr
        
        print *, "Testing validation boundaries..."
        
        ! Test boundary values for various parameters
        bs%rtor = 1.0e-10_dp  ! Very small positive
        bs%rp = 1.0e-10_dp
        bs%B0 = 1.0e-10_dp
        bs%calc_back = 1
        bs%N = 1  ! Minimum odd value
        bs%m_i = 1.0e-10_dp
        bs%zele = 1.0e-10_dp
        bs%zion = 1.0e-10_dp
        
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. .not. is_valid) then
            print *, "FAIL: Boundary values marked as invalid"
            print *, "Error:", trim(error_msg)
            test_status = test_status + 1
        else
            print *, "Correctly accepted boundary values"
        end if
        
        ! Test just below boundary
        bs%rtor = 0.0_dp  ! Zero (invalid)
        call back_sett_validate(bs, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Zero boundary value not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected zero boundary value"
        end if
        
        print *, "test_validation_boundaries completed"
    end subroutine test_validation_boundaries
    
    subroutine test_consistency_validations()
        type(settings_t), pointer :: sd
        logical :: is_valid
        character(len=1024) :: error_msg
        integer :: ierr
        
        print *, "Testing consistency validations..."
        
        call settings_create(sd, "/test/path/", ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Set up settings with inconsistent dimensions
        sd%antenna_settings%dma = 3
        allocate(sd%antenna_settings%modes(4))  ! Should be 6 for dma=3
        
        sd%output_settings%num_quants = 2
        allocate(sd%output_settings%flag_quants(3))  ! Should be 2
        
        call settings_validate_consistency(sd, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS .or. is_valid) then
            print *, "FAIL: Consistency errors not detected"
            test_status = test_status + 1
        else
            print *, "Correctly detected consistency errors:", trim(error_msg)
        end if
        
        call settings_destroy(sd, ierr)
        
        print *, "test_consistency_validations completed"
    end subroutine test_consistency_validations

end program test_settings_validation