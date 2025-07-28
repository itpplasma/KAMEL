program test_settings_copy_compare
    use iso_fortran_env, only: real64, int32, int64
    use iso_c_binding
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    ! Test variables
    integer :: test_status = 0
    
    ! Test 1: Test antenna copy and comparison
    call test_antenna_copy_compare()
    
    ! Test 2: Test background copy and comparison
    call test_background_copy_compare()
    
    ! Test 3: Test output copy and comparison
    call test_output_copy_compare()
    
    ! Test 4: Test eigenmode copy and comparison
    call test_eigenmode_copy_compare()
    
    ! Test 5: Test complete settings deep copy (already implemented)
    call test_settings_deep_copy()
    
    ! Test 6: Test settings comparison
    call test_settings_comparison()
    
    ! Test 7: Test copy error handling
    call test_copy_error_handling()
    
    ! Test 8: Test comparison edge cases
    call test_comparison_edge_cases()
    
    if (test_status == 0) then
        print *, "All settings copy and comparison tests PASSED!"
    else
        print *, "Some tests FAILED. Status:", test_status
        stop 1
    end if
    
contains

    subroutine test_antenna_copy_compare()
        type(antenna_t) :: ant1, ant2
        logical :: is_equal
        integer :: ierr
        
        print *, "Testing antenna copy and comparison..."
        
        ! Set up original antenna
        ant1%ra = 10.5_dp
        ant1%wa = 1.2_dp
        ant1%I0 = 1000.0_dp
        ant1%flab = cmplx(5.0e7_dp, 2.0e3_dp, dp)
        ant1%dma = 2
        allocate(ant1%modes(4))
        ant1%modes = [1, 1, 2, 2]
        ant1%flag_debug = 1
        ant1%flag_eigmode = 0
        
        ! Test deep copy
        call antenna_deep_copy(ant1, ant2, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: antenna_deep_copy failed"
            test_status = test_status + 1
        end if
        
        ! Test comparison - should be equal after copy
        call antenna_compare(ant1, ant2, is_equal, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: antenna_compare failed"
            test_status = test_status + 1
        end if
        if (.not. is_equal) then
            print *, "FAIL: Copied antenna not equal to original"
            test_status = test_status + 1
        end if
        
        ! Test independence - modify copy, should not affect original
        ant2%ra = 20.0_dp
        call antenna_compare(ant1, ant2, is_equal, ierr)
        if (is_equal) then
            print *, "FAIL: Antenna copies are not independent"
            test_status = test_status + 1
        end if
        
        ! Test array independence
        if (allocated(ant2%modes)) then
            ant2%modes(1) = 999
            if (ant1%modes(1) == 999) then
                print *, "FAIL: Antenna mode arrays are not independent"
                test_status = test_status + 1
            end if
        end if
        
        ! Clean up
        if (allocated(ant1%modes)) deallocate(ant1%modes)
        if (allocated(ant2%modes)) deallocate(ant2%modes)
        
        print *, "test_antenna_copy_compare completed"
    end subroutine test_antenna_copy_compare
    
    subroutine test_background_copy_compare()
        type(back_sett_t) :: bs1, bs2
        logical :: is_equal
        integer :: ierr
        
        print *, "Testing background copy and comparison..."
        
        ! Set up original background
        bs1%rtor = 625.0_dp
        bs1%rp = 99.0_dp
        bs1%B0 = 20000.0_dp
        bs1%calc_back = 1
        bs1%flag_back = "normal"
        bs1%N = 5
        bs1%V_gal_sys = 0.0_dp
        bs1%V_scale = 1.0_dp
        bs1%m_i = 1.0_dp
        bs1%zele = 1.0_dp
        bs1%zion = 1.0_dp
        bs1%flag_debug = 1
        bs1%path2profiles = "/test/path/profiles"
        allocate(bs1%mass(2))
        allocate(bs1%charge(2))
        bs1%mass = [1.67e-24_dp, 9.11e-28_dp]
        bs1%charge = [1.602e-19_dp, -1.602e-19_dp]
        bs1%huge_factor = 1.0e30_dp
        
        ! Test deep copy
        call back_sett_deep_copy(bs1, bs2, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: back_sett_deep_copy failed"
            test_status = test_status + 1
        end if
        
        ! Test comparison - should be equal after copy
        call back_sett_compare(bs1, bs2, is_equal, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: back_sett_compare failed"
            test_status = test_status + 1
        end if
        if (.not. is_equal) then
            print *, "FAIL: Copied background not equal to original"
            test_status = test_status + 1
        end if
        
        ! Test independence - modify copy
        bs2%rtor = 1000.0_dp
        call back_sett_compare(bs1, bs2, is_equal, ierr)
        if (is_equal) then
            print *, "FAIL: Background copies are not independent"
            test_status = test_status + 1
        end if
        
        ! Test string independence
        if (allocated(bs2%flag_back)) then
            bs2%flag_back = "modified"
            if (bs1%flag_back == "modified") then
                print *, "FAIL: Background string fields are not independent"
                test_status = test_status + 1
            end if
        end if
        
        ! Clean up
        if (allocated(bs1%mass)) deallocate(bs1%mass)
        if (allocated(bs1%charge)) deallocate(bs1%charge)
        if (allocated(bs1%flag_back)) deallocate(bs1%flag_back)
        if (allocated(bs1%path2profiles)) deallocate(bs1%path2profiles)
        if (allocated(bs2%mass)) deallocate(bs2%mass)
        if (allocated(bs2%charge)) deallocate(bs2%charge)
        if (allocated(bs2%flag_back)) deallocate(bs2%flag_back)
        if (allocated(bs2%path2profiles)) deallocate(bs2%path2profiles)
        
        print *, "test_background_copy_compare completed"
    end subroutine test_background_copy_compare
    
    subroutine test_output_copy_compare()
        type(output_sett_t) :: os1, os2
        logical :: is_equal
        integer :: ierr
        
        print *, "Testing output copy and comparison..."
        
        ! Set up original output
        os1%flag_background = 1
        os1%flag_emfield = 1
        os1%flag_additional = 0
        os1%flag_dispersion = 0
        os1%num_quants = 3
        allocate(os1%flag_quants(3))
        os1%flag_quants = [1, 0, 1]
        os1%flag_debug = 1
        
        ! Test deep copy
        call output_sett_deep_copy(os1, os2, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: output_sett_deep_copy failed"
            test_status = test_status + 1
        end if
        
        ! Test comparison - should be equal after copy
        call output_sett_compare(os1, os2, is_equal, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: output_sett_compare failed"
            test_status = test_status + 1
        end if
        if (.not. is_equal) then
            print *, "FAIL: Copied output not equal to original"
            test_status = test_status + 1
        end if
        
        ! Test independence
        os2%flag_background = 2
        call output_sett_compare(os1, os2, is_equal, ierr)
        if (is_equal) then
            print *, "FAIL: Output copies are not independent"
            test_status = test_status + 1
        end if
        
        ! Test array independence
        if (allocated(os2%flag_quants)) then
            os2%flag_quants(1) = 999
            if (os1%flag_quants(1) == 999) then
                print *, "FAIL: Output flag arrays are not independent"
                test_status = test_status + 1
            end if
        end if
        
        ! Clean up
        if (allocated(os1%flag_quants)) deallocate(os1%flag_quants)
        if (allocated(os2%flag_quants)) deallocate(os2%flag_quants)
        
        print *, "test_output_copy_compare completed"
    end subroutine test_output_copy_compare
    
    subroutine test_eigenmode_copy_compare()
        type(eigmode_sett_t) :: es1, es2
        logical :: is_equal
        integer :: ierr
        
        print *, "Testing eigenmode copy and comparison..."
        
        ! Set up original eigenmode
        es1%fname = "output.dat"
        es1%search_flag = 1
        es1%rdim = 100
        es1%rfmin = 0.0_dp
        es1%rfmax = 1.0e9_dp
        es1%idim = 100
        es1%ifmin = -1.0e6_dp
        es1%ifmax = 1.0e6_dp
        es1%stop_flag = 0
        es1%eps_res = 1.0e-6_dp
        es1%eps_abs = 1.0e-8_dp
        es1%eps_rel = 1.0e-6_dp
        es1%delta = 1.0e-6_dp
        es1%test_roots = 0
        es1%flag_debug = 1
        es1%Nguess = 2
        es1%kmin = 1
        es1%kmax = 10
        allocate(es1%fstart(2))
        es1%fstart = [cmplx(1.0e6_dp, 100.0_dp, dp), cmplx(2.0e6_dp, -50.0_dp, dp)]
        es1%n_zeros = 10
        es1%use_winding = 0
        
        ! Test deep copy
        call eigmode_sett_deep_copy(es1, es2, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: eigmode_sett_deep_copy failed"
            test_status = test_status + 1
        end if
        
        ! Test comparison - should be equal after copy
        call eigmode_sett_compare(es1, es2, is_equal, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: eigmode_sett_compare failed"
            test_status = test_status + 1
        end if
        if (.not. is_equal) then
            print *, "FAIL: Copied eigenmode not equal to original"
            test_status = test_status + 1
        end if
        
        ! Test independence
        es2%rdim = 200
        call eigmode_sett_compare(es1, es2, is_equal, ierr)
        if (is_equal) then
            print *, "FAIL: Eigenmode copies are not independent"
            test_status = test_status + 1
        end if
        
        ! Test complex array independence
        if (allocated(es2%fstart)) then
            es2%fstart(1) = cmplx(999.0_dp, 888.0_dp, dp)
            if (real(es1%fstart(1)) == 999.0_dp) then
                print *, "FAIL: Eigenmode complex arrays are not independent"
                test_status = test_status + 1
            end if
        end if
        
        ! Clean up
        if (allocated(es1%fstart)) deallocate(es1%fstart)
        if (allocated(es2%fstart)) deallocate(es2%fstart)
        if (allocated(es1%fname)) deallocate(es1%fname)
        if (allocated(es2%fname)) deallocate(es2%fname)
        
        print *, "test_eigenmode_copy_compare completed"
    end subroutine test_eigenmode_copy_compare
    
    subroutine test_settings_deep_copy()
        type(settings_t), pointer :: sd1, sd2
        logical :: is_equal
        integer :: ierr
        
        print *, "Testing settings deep copy (already implemented)..."
        
        call settings_create(sd1, "/test/path/", ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not create first settings"
            test_status = test_status + 1
            return
        end if
        
        ! Set up original settings
        sd1%antenna_settings%ra = 10.5_dp
        sd1%antenna_settings%dma = 1
        allocate(sd1%antenna_settings%modes(2))
        sd1%antenna_settings%modes = [1, 1]
        
        sd1%background_settings%rtor = 625.0_dp
        allocate(sd1%background_settings%mass(2))
        sd1%background_settings%mass = [1.0_dp, 2.0_dp]
        
        ! Test deep copy (already implemented)
        call settings_deep_copy(sd1, sd2, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: settings_deep_copy failed"
            test_status = test_status + 1
        end if
        
        ! Test settings comparison
        call settings_compare(sd1, sd2, is_equal, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: settings_compare failed"
            test_status = test_status + 1
        end if
        if (.not. is_equal) then
            print *, "FAIL: Deep copied settings not equal to original"
            test_status = test_status + 1
        end if
        
        ! Test independence
        sd2%antenna_settings%ra = 999.0_dp
        call settings_compare(sd1, sd2, is_equal, ierr)
        if (is_equal) then
            print *, "FAIL: Settings copies are not independent"
            test_status = test_status + 1
        end if
        
        call settings_destroy(sd1, ierr)
        call settings_destroy(sd2, ierr)
        
        print *, "test_settings_deep_copy completed"
    end subroutine test_settings_deep_copy
    
    subroutine test_settings_comparison()
        type(settings_t), pointer :: sd1, sd2
        logical :: is_equal
        integer :: ierr
        
        print *, "Testing settings comparison functionality..."
        
        call settings_create(sd1, "/test/path1/", ierr)
        call settings_create(sd2, "/test/path2/", ierr)
        
        ! Test different paths
        call settings_compare(sd1, sd2, is_equal, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: settings_compare failed for different paths"
            test_status = test_status + 1
        end if
        if (is_equal) then
            print *, "FAIL: Settings with different paths should not be equal"
            test_status = test_status + 1
        end if
        
        ! Make paths same
        sd2%path2project = sd1%path2project
        
        ! Now they should be equal (default values)
        call settings_compare(sd1, sd2, is_equal, ierr)
        if (.not. is_equal) then
            print *, "FAIL: Default settings should be equal"
            test_status = test_status + 1
        end if
        
        call settings_destroy(sd1, ierr)
        call settings_destroy(sd2, ierr)
        
        print *, "test_settings_comparison completed"
    end subroutine test_settings_comparison
    
    subroutine test_copy_error_handling()
        type(antenna_t) :: ant1, ant2
        integer :: ierr
        
        print *, "Testing copy error handling..."
        
        ! Test copy with unallocated source array
        ant1%dma = 2
        ! Don't allocate modes array
        
        call antenna_deep_copy(ant1, ant2, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "Copy with unallocated array handled correctly"
        end if
        
        print *, "test_copy_error_handling completed"
    end subroutine test_copy_error_handling
    
    subroutine test_comparison_edge_cases()
        type(back_sett_t) :: bs1, bs2
        logical :: is_equal
        integer :: ierr
        
        print *, "Testing comparison edge cases..."
        
        ! Test comparison with different array sizes
        bs1%rtor = 100.0_dp
        bs2%rtor = 100.0_dp
        allocate(bs1%mass(2))
        allocate(bs2%mass(3))  ! Different size
        bs1%mass = [1.0_dp, 2.0_dp]
        bs2%mass = [1.0_dp, 2.0_dp, 3.0_dp]
        
        call back_sett_compare(bs1, bs2, is_equal, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: back_sett_compare failed for different array sizes"
            test_status = test_status + 1
        end if
        if (is_equal) then
            print *, "FAIL: Settings with different array sizes should not be equal"
            test_status = test_status + 1
        end if
        
        ! Test comparison with one array allocated, one not
        deallocate(bs2%mass)
        call back_sett_compare(bs1, bs2, is_equal, ierr)
        if (is_equal) then
            print *, "FAIL: Settings with different allocation states should not be equal"
            test_status = test_status + 1
        end if
        
        ! Clean up
        if (allocated(bs1%mass)) deallocate(bs1%mass)
        if (allocated(bs2%mass)) deallocate(bs2%mass)
        
        print *, "test_comparison_edge_cases completed"
    end subroutine test_comparison_edge_cases

end program test_settings_copy_compare