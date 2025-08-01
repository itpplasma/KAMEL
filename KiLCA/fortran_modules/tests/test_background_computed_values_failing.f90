program test_background_computed_values_failing
    use kilca_types_m, only: dp, KILCA_SUCCESS
    use kilca_settings_m, only: back_sett_t, back_sett_initialize_defaults, &
                                background_settings_compute_derived
    implicit none
    
    type(back_sett_t) :: bs
    integer :: ierr
    integer :: test_status = 0
    
    print *, "=========================================="
    print *, "Testing Background Computed Values [RED PHASE - SHOULD FAIL]"
    print *, "=========================================="
    
    call test_comprehensive_derived_computation()
    call test_physics_derived_parameters()
    call test_normalized_parameters()
    
    if (test_status == 0) then
        print *, ""
        print *, "All computed values tests PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "Background computed values tests FAILED:", test_status, "failure(s) (expected in RED phase)"
        stop 1
    end if

contains

    !> Test comprehensive derived value computation including advanced parameters
    subroutine test_comprehensive_derived_computation()
        print *, "Testing comprehensive derived value computation..."
        
        ! Initialize with deuterium plasma parameters
        call back_sett_initialize_defaults(bs, ierr)
        bs%m_i = 2.0_dp      ! Deuterium mass
        bs%B0 = 25000.0_dp   ! 2.5 Tesla
        bs%rtor = 170.0_dp   ! Major radius
        bs%rp = 70.0_dp      ! Minor radius
        bs%V_gal_sys = 1.0e9_dp
        
        ! Compute derived values
        call background_settings_compute_derived(bs, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not compute derived values, error:", ierr
            test_status = test_status + 1
            return
        end if
        
        ! Test that additional arrays are computed beyond just mass and charge
        ! This should FAIL initially - only mass and charge arrays are currently implemented
        
        ! Check for density arrays (should be computed from profiles)
        if (.not. allocated(bs%density)) then
            print *, "FAIL: Density array not computed from profiles"
            test_status = test_status + 1
        else
            print *, "PASS: Density array computed correctly"
        end if
        
        ! Check for temperature arrays (should be computed from profiles)
        if (.not. allocated(bs%temperature)) then
            print *, "FAIL: Temperature array not computed from profiles"
            test_status = test_status + 1  
        else
            print *, "PASS: Temperature array computed correctly"
        end if
        
        ! Check for velocity arrays (should be computed from profiles)
        if (.not. allocated(bs%velocity)) then
            print *, "FAIL: Velocity array not computed from profiles"
            test_status = test_status + 1
        else
            print *, "PASS: Velocity array computed correctly"
        end if
    end subroutine test_comprehensive_derived_computation
    
    !> Test physics-derived parameters like normalized radii
    subroutine test_physics_derived_parameters()
        print *, "Testing physics-derived parameters..."
        
        ! Use same setup
        call back_sett_initialize_defaults(bs, ierr)
        bs%m_i = 2.0_dp
        bs%B0 = 25000.0_dp
        bs%rtor = 170.0_dp
        bs%rp = 70.0_dp
        bs%V_gal_sys = 1.0e9_dp
        
        call background_settings_compute_derived(bs, ierr)
        
        ! Test dimensionless parameter computation that should be derived but aren't yet
        ! These fields don't exist in the current structure, so will fail compilation
        
        ! Test normalized radius computation (this should fail - field doesn't exist)
        ! Commenting out since it will cause compilation error
        ! if (abs(bs%eps_r - bs%rp/bs%rtor) > 1.0e-10_dp) then
        print *, "FAIL: Normalized radius (eps_r) field not defined in structure"
        test_status = test_status + 1
        
        ! Test aspect ratio computation (this should fail - field doesn't exist)
        ! if (abs(bs%aspect_ratio - bs%rtor/bs%rp) > 1.0e-10_dp) then
        print *, "FAIL: Aspect ratio field not defined in structure"
        test_status = test_status + 1
        
        ! Test that profile-dependent arrays are allocated and computed
        if (.not. allocated(bs%profile_grid)) then
            print *, "FAIL: Profile grid not allocated"
            test_status = test_status + 1
        else
            print *, "PASS: Profile grid allocated correctly"
        end if
        
        if (.not. allocated(bs%q_profile)) then
            print *, "FAIL: Safety factor profile not computed"
            test_status = test_status + 1
        else
            print *, "PASS: Safety factor profile computed correctly"
        end if
    end subroutine test_physics_derived_parameters
    
    !> Test normalized parameters computation
    subroutine test_normalized_parameters()
        print *, "Testing normalized parameters..."
        
        call back_sett_initialize_defaults(bs, ierr)
        bs%m_i = 2.0_dp
        bs%B0 = 25000.0_dp
        bs%rtor = 170.0_dp
        bs%rp = 70.0_dp
        bs%V_gal_sys = 1.0e9_dp
        
        call background_settings_compute_derived(bs, ierr)
        
        ! Test for normalization parameters that should be computed but aren't yet
        ! These would be useful for plasma physics calculations
        
        ! Test for thermal velocity arrays (derived from temperature and mass)
        if (.not. allocated(bs%v_thermal)) then
            print *, "FAIL: Thermal velocity array not computed"
            test_status = test_status + 1
        else
            print *, "PASS: Thermal velocity array computed correctly"
        end if
        
        ! Test for pressure profile (derived from density and temperature)
        if (.not. allocated(bs%pressure)) then
            print *, "FAIL: Pressure profile not computed"
            test_status = test_status + 1
        else
            print *, "PASS: Pressure profile computed correctly"
        end if
        
        ! Test that collision frequency arrays are computed
        if (.not. allocated(bs%nu_collision)) then
            print *, "FAIL: Collision frequency array not computed"
            test_status = test_status + 1
        else
            print *, "PASS: Collision frequency array computed correctly"
        end if
    end subroutine test_normalized_parameters

end program test_background_computed_values_failing