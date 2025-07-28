program test_kilca_settings_complete
    use iso_fortran_env, only: real64, int32, int64
    use iso_c_binding
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    ! Test variables
    integer :: test_status = 0
    character(len=256) :: test_path = "/test/project/path/"
    
    ! Test complete background settings structure
    call test_complete_background_settings()
    
    ! Test complete output settings structure  
    call test_complete_output_settings()
    
    ! Test complete eigenmode settings structure
    call test_complete_eigmode_settings()
    
    if (test_status == 0) then
        print *, "All complete settings tests PASSED!"
    else
        print *, "Some tests FAILED. Status:", test_status
        stop 1
    end if
    
contains

    subroutine test_complete_background_settings()
        type(settings_t), pointer :: sd
        type(back_sett_t), pointer :: bs
        integer :: ierr
        
        print *, "Testing complete background settings..."
        
        call settings_create(sd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Get background settings pointer
        call settings_get_background(sd, bs, ierr)
        
        ! Test all member variables exist and can be set
        bs%rtor = 625.0_dp  ! Machine radius
        bs%rp = 99.0_dp     ! Plasma radius
        bs%B0 = 20000.0_dp  ! Magnetic field
        
        bs%calc_back = 1
        bs%flag_back = "normal"
        bs%N = 5
        bs%V_gal_sys = 0.0_dp
        bs%V_scale = 1.0_dp
        bs%m_i = 2.0_dp
        bs%zele = 1.0_dp
        bs%zion = 1.0_dp
        bs%flag_debug = 0
        bs%huge_factor = 1.0e30_dp
        
        ! Verify settings
        if (abs(bs%rtor - 625.0_dp) > epsilon(1.0_dp)) then
            print *, "FAIL: rtor not set correctly"
            test_status = test_status + 1
        end if
        
        call settings_destroy(sd, ierr)
        
        print *, "test_complete_background_settings completed"
    end subroutine test_complete_background_settings
    
    subroutine test_complete_output_settings()
        type(settings_t), pointer :: sd
        type(output_sett_t), pointer :: os
        integer :: ierr
        
        print *, "Testing complete output settings..."
        
        call settings_create(sd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Get output settings pointer
        call settings_get_output(sd, os, ierr)
        
        ! Test all member variables
        os%flag_background = 1
        os%flag_emfield = 2
        os%flag_additional = 1
        os%flag_dispersion = 2
        os%flag_debug = 0
        
        ! Test flag_quants array
        os%num_quants = 5
        if (allocated(os%flag_quants)) deallocate(os%flag_quants)
        allocate(os%flag_quants(os%num_quants))
        os%flag_quants = [1, 1, 0, 1, 0]
        
        ! Verify
        if (os%flag_background /= 1) then
            print *, "FAIL: flag_background not set"
            test_status = test_status + 1
        end if
        
        if (.not. allocated(os%flag_quants)) then
            print *, "FAIL: flag_quants not allocated"
            test_status = test_status + 1
        end if
        
        call settings_destroy(sd, ierr)
        
        print *, "test_complete_output_settings completed"
    end subroutine test_complete_output_settings
    
    subroutine test_complete_eigmode_settings()
        type(settings_t), pointer :: sd
        type(eigmode_sett_t), pointer :: es
        integer :: ierr
        
        print *, "Testing complete eigenmode settings..."
        
        call settings_create(sd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Get eigenmode settings pointer
        call settings_get_eigmode(sd, es, ierr)
        
        ! Test all member variables
        es%fname = "eigenmodes.out"
        es%search_flag = 1
        
        ! Grid parameters
        es%rdim = 100
        es%rfmin = 0.0_dp
        es%rfmax = 1.0e9_dp
        es%idim = 50
        es%ifmin = -1.0e6_dp
        es%ifmax = 1.0e6_dp
        
        es%stop_flag = 0
        
        ! Accuracies
        es%eps_res = 1.0e-6_dp
        es%eps_abs = 1.0e-8_dp
        es%eps_rel = 1.0e-6_dp
        es%delta = 1.0e-6_dp
        
        es%test_roots = 1
        es%flag_debug = 0
        
        ! Starting values
        es%Nguess = 5
        es%kmin = 1
        es%kmax = 10
        
        ! Allocate and set fstart array
        if (allocated(es%fstart)) deallocate(es%fstart)
        allocate(es%fstart(es%Nguess))
        es%fstart = [(cmplx(1.0e6_dp, 1.0e3_dp, dp), ierr=1,es%Nguess)]
        
        es%n_zeros = 10
        es%use_winding = 1
        
        ! Verify
        if (es%rdim /= 100) then
            print *, "FAIL: rdim not set"
            test_status = test_status + 1
        end if
        
        if (.not. allocated(es%fstart)) then
            print *, "FAIL: fstart not allocated"
            test_status = test_status + 1
        end if
        
        call settings_destroy(sd, ierr)
        
        print *, "test_complete_eigmode_settings completed"
    end subroutine test_complete_eigmode_settings

end program test_kilca_settings_complete