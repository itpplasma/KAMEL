program test_settings_initialization
    use iso_fortran_env, only: real64, int32, int64
    use iso_c_binding
    use kilca_types_m
    use kilca_settings_m
    implicit none
    
    ! Test variables
    integer :: test_status = 0
    
    ! Test 1: Test antenna default initialization
    call test_antenna_default_initialization()
    
    ! Test 2: Test background default initialization
    call test_background_default_initialization()
    
    ! Test 3: Test output default initialization
    call test_output_default_initialization()
    
    ! Test 4: Test eigenmode default initialization
    call test_eigenmode_default_initialization()
    
    ! Test 5: Test complete settings default initialization
    call test_settings_default_initialization()
    
    ! Test 6: Test custom initialization with parameters
    call test_custom_initialization()
    
    ! Test 7: Test initialization reset/reinitialize
    call test_initialization_reset()
    
    ! Test 8: Test initialization validation
    call test_initialization_validation()
    
    if (test_status == 0) then
        print *, "All settings initialization tests PASSED!"
    else
        print *, "Some tests FAILED. Status:", test_status
        stop 1
    end if
    
contains

    subroutine test_antenna_default_initialization()
        type(antenna_t) :: ant
        integer :: ierr
        
        print *, "Testing antenna default initialization..."
        
        ! Test initialize with defaults
        call antenna_initialize_defaults(ant, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: antenna_initialize_defaults failed"
            test_status = test_status + 1
        end if
        
        ! Verify default values are set correctly
        if (ant%ra /= 0.0_dp) then
            print *, "FAIL: Default ra not set correctly, got:", ant%ra
            test_status = test_status + 1
        end if
        
        if (ant%wa /= 0.0_dp) then
            print *, "FAIL: Default wa not set correctly, got:", ant%wa
            test_status = test_status + 1
        end if
        
        if (ant%I0 /= 0.0_dp) then
            print *, "FAIL: Default I0 not set correctly, got:", ant%I0
            test_status = test_status + 1
        end if
        
        if (ant%flab /= cmplx_zero) then
            print *, "FAIL: Default flab not set correctly, got:", ant%flab
            test_status = test_status + 1
        end if
        
        if (ant%dma /= 0) then
            print *, "FAIL: Default dma not set correctly, got:", ant%dma
            test_status = test_status + 1
        end if
        
        if (ant%flag_debug /= 0) then
            print *, "FAIL: Default flag_debug not set correctly, got:", ant%flag_debug
            test_status = test_status + 1
        end if
        
        if (ant%flag_eigmode /= 0) then
            print *, "FAIL: Default flag_eigmode not set correctly, got:", ant%flag_eigmode
            test_status = test_status + 1
        end if
        
        ! Arrays should not be allocated by default
        if (allocated(ant%modes)) then
            print *, "FAIL: modes array should not be allocated by default"
            test_status = test_status + 1
        end if
        
        print *, "test_antenna_default_initialization completed"
    end subroutine test_antenna_default_initialization
    
    subroutine test_background_default_initialization()
        type(back_sett_t) :: bs
        integer :: ierr
        
        print *, "Testing background default initialization..."
        
        ! Test initialize with defaults
        call back_sett_initialize_defaults(bs, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: back_sett_initialize_defaults failed"
            test_status = test_status + 1
        end if
        
        ! Verify meaningful default values
        if (abs(bs%rtor - 625.0_dp) > epsilon(bs%rtor)) then
            print *, "FAIL: Default rtor not set correctly, got:", bs%rtor
            test_status = test_status + 1
        end if
        
        if (abs(bs%rp - 200.0_dp) > epsilon(bs%rp)) then
            print *, "FAIL: Default rp not set correctly, got:", bs%rp
            test_status = test_status + 1
        end if
        
        if (abs(bs%B0 - 20000.0_dp) > epsilon(bs%B0)) then
            print *, "FAIL: Default B0 not set correctly, got:", bs%B0
            test_status = test_status + 1
        end if
        
        if (bs%calc_back /= 1) then
            print *, "FAIL: Default calc_back not set correctly, got:", bs%calc_back
            test_status = test_status + 1
        end if
        
        if (bs%N /= 5) then
            print *, "FAIL: Default N not set correctly, got:", bs%N
            test_status = test_status + 1
        end if
        
        if (abs(bs%V_gal_sys - 0.0_dp) > epsilon(bs%V_gal_sys)) then
            print *, "FAIL: Default V_gal_sys not set correctly, got:", bs%V_gal_sys
            test_status = test_status + 1
        end if
        
        if (abs(bs%V_scale - 1.0_dp) > epsilon(bs%V_scale)) then
            print *, "FAIL: Default V_scale not set correctly, got:", bs%V_scale
            test_status = test_status + 1
        end if
        
        if (abs(bs%m_i - 1.0_dp) > epsilon(bs%m_i)) then
            print *, "FAIL: Default m_i not set correctly, got:", bs%m_i
            test_status = test_status + 1
        end if
        
        if (abs(bs%zele - 1.0_dp) > epsilon(bs%zele)) then
            print *, "FAIL: Default zele not set correctly, got:", bs%zele
            test_status = test_status + 1
        end if
        
        if (abs(bs%zion - 1.0_dp) > epsilon(bs%zion)) then
            print *, "FAIL: Default zion not set correctly, got:", bs%zion
            test_status = test_status + 1
        end if
        
        if (bs%flag_debug /= 0) then
            print *, "FAIL: Default flag_debug not set correctly, got:", bs%flag_debug
            test_status = test_status + 1
        end if
        
        if (abs(bs%huge_factor - 1.0e30_dp) > epsilon(bs%huge_factor)) then
            print *, "FAIL: Default huge_factor not set correctly, got:", bs%huge_factor
            test_status = test_status + 1
        end if
        
        ! Test default string allocation
        if (.not. allocated(bs%flag_back)) then
            print *, "FAIL: flag_back should be allocated with default"
            test_status = test_status + 1
        else
            if (bs%flag_back /= "normal") then
                print *, "FAIL: Default flag_back not set correctly, got:", bs%flag_back
                test_status = test_status + 1
            end if
        end if
        
        ! Clean up
        if (allocated(bs%flag_back)) deallocate(bs%flag_back)
        if (allocated(bs%path2profiles)) deallocate(bs%path2profiles)
        if (allocated(bs%mass)) deallocate(bs%mass)
        if (allocated(bs%charge)) deallocate(bs%charge)
        
        print *, "test_background_default_initialization completed"
    end subroutine test_background_default_initialization
    
    subroutine test_output_default_initialization()
        type(output_sett_t) :: os
        integer :: ierr
        
        print *, "Testing output default initialization..."
        
        ! Test initialize with defaults
        call output_sett_initialize_defaults(os, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: output_sett_initialize_defaults failed"
            test_status = test_status + 1
        end if
        
        ! Verify default values
        if (os%flag_background /= 1) then
            print *, "FAIL: Default flag_background not set correctly, got:", os%flag_background
            test_status = test_status + 1
        end if
        
        if (os%flag_emfield /= 1) then
            print *, "FAIL: Default flag_emfield not set correctly, got:", os%flag_emfield
            test_status = test_status + 1
        end if
        
        if (os%flag_additional /= 0) then
            print *, "FAIL: Default flag_additional not set correctly, got:", os%flag_additional
            test_status = test_status + 1
        end if
        
        if (os%flag_dispersion /= 0) then
            print *, "FAIL: Default flag_dispersion not set correctly, got:", os%flag_dispersion
            test_status = test_status + 1
        end if
        
        if (os%num_quants /= 0) then
            print *, "FAIL: Default num_quants not set correctly, got:", os%num_quants
            test_status = test_status + 1
        end if
        
        if (os%flag_debug /= 0) then
            print *, "FAIL: Default flag_debug not set correctly, got:", os%flag_debug
            test_status = test_status + 1
        end if
        
        ! Arrays should not be allocated by default
        if (allocated(os%flag_quants)) then
            print *, "FAIL: flag_quants array should not be allocated by default"
            test_status = test_status + 1
        end if
        
        print *, "test_output_default_initialization completed"
    end subroutine test_output_default_initialization
    
    subroutine test_eigenmode_default_initialization()
        type(eigmode_sett_t) :: es
        integer :: ierr
        
        print *, "Testing eigenmode default initialization..."
        
        ! Test initialize with defaults
        call eigmode_sett_initialize_defaults(es, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: eigmode_sett_initialize_defaults failed"
            test_status = test_status + 1
        end if
        
        ! Verify meaningful defaults for eigenmode calculations
        if (es%search_flag /= 0) then
            print *, "FAIL: Default search_flag not set correctly, got:", es%search_flag
            test_status = test_status + 1
        end if
        
        if (es%rdim /= 100) then
            print *, "FAIL: Default rdim not set correctly, got:", es%rdim
            test_status = test_status + 1
        end if
        
        if (abs(es%rfmin - 0.0_dp) > epsilon(es%rfmin)) then
            print *, "FAIL: Default rfmin not set correctly, got:", es%rfmin
            test_status = test_status + 1
        end if
        
        if (abs(es%rfmax - 1.0e9_dp) > epsilon(es%rfmax)) then
            print *, "FAIL: Default rfmax not set correctly, got:", es%rfmax
            test_status = test_status + 1
        end if
        
        if (es%idim /= 100) then
            print *, "FAIL: Default idim not set correctly, got:", es%idim
            test_status = test_status + 1
        end if
        
        if (abs(es%ifmin - (-1.0e6_dp)) > epsilon(es%ifmin)) then
            print *, "FAIL: Default ifmin not set correctly, got:", es%ifmin
            test_status = test_status + 1
        end if
        
        if (abs(es%ifmax - 1.0e6_dp) > epsilon(es%ifmax)) then
            print *, "FAIL: Default ifmax not set correctly, got:", es%ifmax
            test_status = test_status + 1
        end if
        
        if (es%stop_flag /= 0) then
            print *, "FAIL: Default stop_flag not set correctly, got:", es%stop_flag
            test_status = test_status + 1
        end if
        
        if (abs(es%eps_res - 1.0e-6_dp) > epsilon(es%eps_res)) then
            print *, "FAIL: Default eps_res not set correctly, got:", es%eps_res
            test_status = test_status + 1
        end if
        
        if (abs(es%eps_abs - 1.0e-8_dp) > epsilon(es%eps_abs)) then
            print *, "FAIL: Default eps_abs not set correctly, got:", es%eps_abs
            test_status = test_status + 1
        end if
        
        if (abs(es%eps_rel - 1.0e-6_dp) > epsilon(es%eps_rel)) then
            print *, "FAIL: Default eps_rel not set correctly, got:", es%eps_rel
            test_status = test_status + 1
        end if
        
        if (abs(es%delta - 1.0e-6_dp) > epsilon(es%delta)) then
            print *, "FAIL: Default delta not set correctly, got:", es%delta
            test_status = test_status + 1
        end if
        
        if (es%test_roots /= 0) then
            print *, "FAIL: Default test_roots not set correctly, got:", es%test_roots
            test_status = test_status + 1
        end if
        
        if (es%flag_debug /= 0) then
            print *, "FAIL: Default flag_debug not set correctly, got:", es%flag_debug
            test_status = test_status + 1
        end if
        
        if (es%Nguess /= 0) then
            print *, "FAIL: Default Nguess not set correctly, got:", es%Nguess
            test_status = test_status + 1
        end if
        
        if (es%kmin /= 1) then
            print *, "FAIL: Default kmin not set correctly, got:", es%kmin
            test_status = test_status + 1
        end if
        
        if (es%kmax /= 10) then
            print *, "FAIL: Default kmax not set correctly, got:", es%kmax
            test_status = test_status + 1
        end if
        
        if (es%n_zeros /= 10) then
            print *, "FAIL: Default n_zeros not set correctly, got:", es%n_zeros
            test_status = test_status + 1
        end if
        
        if (es%use_winding /= 0) then
            print *, "FAIL: Default use_winding not set correctly, got:", es%use_winding
            test_status = test_status + 1
        end if
        
        ! Test default string allocation
        if (.not. allocated(es%fname)) then
            print *, "FAIL: fname should be allocated with default"
            test_status = test_status + 1
        else
            if (es%fname /= "eigenmode_output.dat") then
                print *, "FAIL: Default fname not set correctly, got:", es%fname
                test_status = test_status + 1
            end if
        end if
        
        ! Arrays should not be allocated by default
        if (allocated(es%fstart)) then
            print *, "FAIL: fstart array should not be allocated by default"
            test_status = test_status + 1
        end if
        
        ! Clean up
        if (allocated(es%fname)) deallocate(es%fname)
        if (allocated(es%fstart)) deallocate(es%fstart)
        
        print *, "test_eigenmode_default_initialization completed"
    end subroutine test_eigenmode_default_initialization
    
    subroutine test_settings_default_initialization()
        type(settings_t), pointer :: sd
        integer :: ierr
        
        print *, "Testing complete settings default initialization..."
        
        ! Test initialize with defaults
        call settings_initialize_defaults(sd, "/test/path/", ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: settings_initialize_defaults failed"
            test_status = test_status + 1
            return
        end if
        
        ! Check that all sub-settings are properly initialized
        if (sd%path2project /= "/test/path/") then
            print *, "FAIL: Path not set correctly"
            test_status = test_status + 1
        end if
        
        ! Check antenna defaults
        if (sd%antenna_settings%ra /= 0.0_dp) then
            print *, "FAIL: Antenna defaults not set in complete settings"
            test_status = test_status + 1
        end if
        
        ! Check background defaults
        if (abs(sd%background_settings%rtor - 625.0_dp) > epsilon(sd%background_settings%rtor)) then
            print *, "FAIL: Background defaults not set in complete settings"
            test_status = test_status + 1
        end if
        
        ! Check output defaults
        if (sd%output_settings%flag_background /= 1) then
            print *, "FAIL: Output defaults not set in complete settings"
            test_status = test_status + 1
        end if
        
        ! Check eigenmode defaults
        if (sd%eigmode_settings%rdim /= 100) then
            print *, "FAIL: Eigenmode defaults not set in complete settings"
            test_status = test_status + 1
        end if
        
        call settings_destroy(sd, ierr)
        
        print *, "test_settings_default_initialization completed"
    end subroutine test_settings_default_initialization
    
    subroutine test_custom_initialization()
        type(antenna_t) :: ant
        integer :: ierr
        
        print *, "Testing custom initialization with parameters..."
        
        ! Test custom initialization with specific parameters
        call antenna_initialize_custom(ant, ra=10.5_dp, wa=1.2_dp, I0=1000.0_dp, ierr=ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: antenna_initialize_custom failed"
            test_status = test_status + 1
        end if
        
        ! Verify custom values were set
        if (abs(ant%ra - 10.5_dp) > epsilon(ant%ra)) then
            print *, "FAIL: Custom ra not set correctly, got:", ant%ra
            test_status = test_status + 1
        end if
        
        if (abs(ant%wa - 1.2_dp) > epsilon(ant%wa)) then
            print *, "FAIL: Custom wa not set correctly, got:", ant%wa
            test_status = test_status + 1
        end if
        
        if (abs(ant%I0 - 1000.0_dp) > epsilon(ant%I0)) then
            print *, "FAIL: Custom I0 not set correctly, got:", ant%I0
            test_status = test_status + 1
        end if
        
        ! Other values should still be defaults
        if (ant%flab /= cmplx_zero) then
            print *, "FAIL: Non-specified values should remain default"
            test_status = test_status + 1
        end if
        
        print *, "test_custom_initialization completed"
    end subroutine test_custom_initialization
    
    subroutine test_initialization_reset()
        type(back_sett_t) :: bs
        integer :: ierr
        
        print *, "Testing initialization reset..."
        
        ! Set some non-default values
        bs%rtor = 999.0_dp
        bs%calc_back = 99
        bs%flag_back = "modified"
        allocate(bs%mass(2))
        bs%mass = [1.0_dp, 2.0_dp]
        
        ! Reset to defaults
        call back_sett_initialize_defaults(bs, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: reset initialization failed"
            test_status = test_status + 1
        end if
        
        ! Verify values were reset
        if (abs(bs%rtor - 625.0_dp) > epsilon(bs%rtor)) then
            print *, "FAIL: Reset did not work for rtor, got:", bs%rtor
            test_status = test_status + 1
        end if
        
        if (bs%calc_back /= 1) then
            print *, "FAIL: Reset did not work for calc_back, got:", bs%calc_back
            test_status = test_status + 1
        end if
        
        ! Arrays should be deallocated and reallocated as needed
        if (allocated(bs%mass)) then
            print *, "FAIL: Arrays should be deallocated on reset"
            test_status = test_status + 1
        end if
        
        ! Clean up
        if (allocated(bs%flag_back)) deallocate(bs%flag_back)
        if (allocated(bs%path2profiles)) deallocate(bs%path2profiles)
        if (allocated(bs%mass)) deallocate(bs%mass)
        if (allocated(bs%charge)) deallocate(bs%charge)
        
        print *, "test_initialization_reset completed"
    end subroutine test_initialization_reset
    
    subroutine test_initialization_validation()
        type(eigmode_sett_t) :: es
        integer :: ierr
        logical :: is_valid
        character(len=1024) :: error_msg
        
        print *, "Testing initialization validation..."
        
        ! Initialize with defaults
        call eigmode_sett_initialize_defaults(es, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: initialization failed"
            test_status = test_status + 1
            return
        end if
        
        ! Validate that defaults result in valid settings
        call eigmode_sett_validate(es, is_valid, error_msg, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: validation failed after initialization"
            test_status = test_status + 1
        end if
        
        if (.not. is_valid) then
            print *, "FAIL: Default initialization should result in valid settings"
            print *, "Error:", trim(error_msg)
            test_status = test_status + 1
        end if
        
        ! Clean up
        if (allocated(es%fname)) deallocate(es%fname)
        if (allocated(es%fstart)) deallocate(es%fstart)
        
        print *, "test_initialization_validation completed"
    end subroutine test_initialization_validation

end program test_settings_initialization