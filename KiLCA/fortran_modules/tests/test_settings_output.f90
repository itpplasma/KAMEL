program test_settings_output
    use iso_fortran_env, only: real64, int32, int64, output_unit
    use iso_c_binding
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    ! Test variables
    integer :: test_status = 0
    
    ! Test 1: Test antenna print_settings
    call test_antenna_print_settings()
    
    ! Test 2: Test background print_settings
    call test_background_print_settings()
    
    ! Test 3: Test output print_settings
    call test_output_print_settings()
    
    ! Test 4: Test eigenmode print_settings
    call test_eigenmode_print_settings()
    
    ! Test 5: Test complete settings print
    call test_complete_settings_print()
    
    ! Test 6: Test print to file
    call test_print_to_file()
    
    ! Test 7: Test print with various formatting
    call test_print_formatting()
    
    if (test_status == 0) then
        print *, "All settings output tests PASSED!"
    else
        print *, "Some tests FAILED. Status:", test_status
        stop 1
    end if
    
contains

    subroutine test_antenna_print_settings()
        type(antenna_t) :: ant
        integer :: ierr
        
        print *, "Testing antenna print_settings..."
        
        ! Set up antenna with valid data
        ant%ra = 10.5_dp
        ant%wa = 1.2_dp
        ant%I0 = 1000.0_dp
        ant%flab = cmplx(5.0e7_dp, 2.0e3_dp, dp)
        ant%dma = 2
        allocate(ant%modes(4))
        ant%modes = [1, 1, 2, 2]
        ant%flag_debug = 1
        ant%flag_eigmode = 0
        
        ! Test print to stdout
        call antenna_print_settings(ant, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: antenna_print_settings failed"
            test_status = test_status + 1
        end if
        
        ! Clean up
        if (allocated(ant%modes)) deallocate(ant%modes)
        
        print *, "test_antenna_print_settings completed"
    end subroutine test_antenna_print_settings
    
    subroutine test_background_print_settings()
        type(back_sett_t) :: bs
        integer :: ierr
        
        print *, "Testing background print_settings..."
        
        ! Set up background with valid data
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
        bs%flag_debug = 1
        bs%path2profiles = "/test/path/profiles"
        allocate(bs%mass(2))
        allocate(bs%charge(2))
        bs%mass = [1.67e-24_dp, 9.11e-28_dp]
        bs%charge = [1.602e-19_dp, -1.602e-19_dp]
        bs%huge_factor = 1.0e30_dp
        
        ! Test print to stdout
        call back_sett_print_settings(bs, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: back_sett_print_settings failed"
            test_status = test_status + 1
        end if
        
        ! Clean up
        if (allocated(bs%mass)) deallocate(bs%mass)
        if (allocated(bs%charge)) deallocate(bs%charge)
        if (allocated(bs%flag_back)) deallocate(bs%flag_back)
        if (allocated(bs%path2profiles)) deallocate(bs%path2profiles)
        
        print *, "test_background_print_settings completed"
    end subroutine test_background_print_settings
    
    subroutine test_output_print_settings()
        type(output_sett_t) :: os
        integer :: ierr
        
        print *, "Testing output print_settings..."
        
        ! Set up output with valid data
        os%flag_background = 1
        os%flag_emfield = 1
        os%flag_additional = 0
        os%flag_dispersion = 0
        os%num_quants = 3
        allocate(os%flag_quants(3))
        os%flag_quants = [1, 0, 1]
        os%flag_debug = 1
        
        ! Test print to stdout
        call output_sett_print_settings(os, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: output_sett_print_settings failed"
            test_status = test_status + 1
        end if
        
        ! Clean up
        if (allocated(os%flag_quants)) deallocate(os%flag_quants)
        
        print *, "test_output_print_settings completed"
    end subroutine test_output_print_settings
    
    subroutine test_eigenmode_print_settings()
        type(eigmode_sett_t) :: es
        integer :: ierr
        
        print *, "Testing eigenmode print_settings..."
        
        ! Set up eigenmode with valid data
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
        es%flag_debug = 1
        es%Nguess = 2
        es%kmin = 1
        es%kmax = 10
        allocate(es%fstart(2))
        es%fstart = [cmplx(1.0e6_dp, 100.0_dp, dp), cmplx(2.0e6_dp, -50.0_dp, dp)]
        es%n_zeros = 10
        es%use_winding = 0
        
        ! Test print to stdout
        call eigmode_sett_print_settings(es, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: eigmode_sett_print_settings failed"
            test_status = test_status + 1
        end if
        
        ! Clean up
        if (allocated(es%fstart)) deallocate(es%fstart)
        if (allocated(es%fname)) deallocate(es%fname)
        
        print *, "test_eigenmode_print_settings completed"
    end subroutine test_eigenmode_print_settings
    
    subroutine test_complete_settings_print()
        type(settings_t), pointer :: sd
        integer :: ierr
        
        print *, "Testing complete settings print..."
        
        call settings_create(sd, "/test/path/", ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not create settings"
            test_status = test_status + 1
            return
        end if
        
        ! Set up minimal valid settings
        sd%antenna_settings%ra = 10.5_dp
        sd%antenna_settings%wa = 1.2_dp
        sd%antenna_settings%I0 = 1000.0_dp
        sd%antenna_settings%flab = cmplx(5.0e7_dp, 0.0_dp, dp)
        sd%antenna_settings%dma = 1
        allocate(sd%antenna_settings%modes(2))
        sd%antenna_settings%modes = [1, 1]
        
        sd%background_settings%rtor = 625.0_dp
        sd%background_settings%rp = 99.0_dp
        sd%background_settings%B0 = 20000.0_dp
        
        sd%output_settings%flag_background = 1
        sd%output_settings%num_quants = 1
        allocate(sd%output_settings%flag_quants(1))
        sd%output_settings%flag_quants = [1]
        
        sd%eigmode_settings%fname = "test.dat"
        sd%eigmode_settings%Nguess = 1
        allocate(sd%eigmode_settings%fstart(1))
        sd%eigmode_settings%fstart = [cmplx(1.0e6_dp, 0.0_dp, dp)]
        
        ! Test print all settings
        call settings_print_all(sd, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: settings_print_all failed"
            test_status = test_status + 1
        end if
        
        call settings_destroy(sd, ierr)
        
        print *, "test_complete_settings_print completed"
    end subroutine test_complete_settings_print
    
    subroutine test_print_to_file()
        type(antenna_t) :: ant
        integer :: ierr, file_unit
        character(len=256) :: filename
        logical :: file_exists
        
        print *, "Testing print to file..."
        
        filename = "test_antenna_output.txt"
        
        ! Set up antenna
        ant%ra = 10.5_dp
        ant%wa = 1.2_dp
        ant%I0 = 1000.0_dp
        ant%flab = cmplx(5.0e7_dp, 0.0_dp, dp)
        ant%dma = 1
        allocate(ant%modes(2))
        ant%modes = [1, 1]
        ant%flag_debug = 1
        ant%flag_eigmode = 0
        
        ! Open file for writing
        open(newunit=file_unit, file=filename, action='write', status='replace', iostat=ierr)
        if (ierr /= 0) then
            print *, "FAIL: Could not open test file"
            test_status = test_status + 1
            return
        end if
        
        ! Test print to file
        call antenna_print_settings_to_unit(ant, file_unit, ierr)
        close(file_unit)
        
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: antenna_print_settings_to_unit failed"
            test_status = test_status + 1
        end if
        
        ! Check if file was created and has content
        inquire(file=filename, exist=file_exists)
        if (.not. file_exists) then
            print *, "FAIL: Output file was not created"
            test_status = test_status + 1
        else
            ! Clean up test file
            open(newunit=file_unit, file=filename)
            close(file_unit, status='delete')
        end if
        
        ! Clean up
        if (allocated(ant%modes)) deallocate(ant%modes)
        
        print *, "test_print_to_file completed"
    end subroutine test_print_to_file
    
    subroutine test_print_formatting()
        type(back_sett_t) :: bs
        integer :: ierr
        
        print *, "Testing print formatting..."
        
        ! Set up background with edge case values to test formatting
        bs%rtor = 1.23456789e-10_dp  ! Very small
        bs%rp = 9.87654321e10_dp     ! Very large
        bs%B0 = 0.0_dp               ! Zero
        bs%calc_back = -1            ! Negative
        bs%flag_back = "very_long_flag_name_to_test_string_formatting"
        bs%N = 999                   ! Large integer
        bs%V_gal_sys = -1.23e-100_dp ! Very small negative
        bs%V_scale = 1.0_dp
        bs%m_i = 1.0_dp
        bs%zele = 1.0_dp
        bs%zion = 1.0_dp
        bs%flag_debug = 0
        bs%path2profiles = ""        ! Empty string
        bs%huge_factor = 1.0e30_dp
        
        ! Test that formatting handles edge cases
        call back_sett_print_settings(bs, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: back_sett_print_settings failed with edge case formatting"
            test_status = test_status + 1
        end if
        
        ! Clean up
        if (allocated(bs%flag_back)) deallocate(bs%flag_back)
        if (allocated(bs%path2profiles)) deallocate(bs%path2profiles)
        
        print *, "test_print_formatting completed"
    end subroutine test_print_formatting

end program test_settings_output