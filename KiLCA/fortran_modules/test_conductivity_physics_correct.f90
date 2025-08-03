program test_conductivity_physics_correct
    ! Test program for CORRECT conductivity K-matrix physics
    ! This will initially FAIL with current broken implementation
    
    use iso_fortran_env, only: real64, error_unit
    use kilca_conductivity_m
    use kilca_plasma_physics_m
    implicit none
    
    integer, parameter :: N_TESTS = 8
    integer :: test_count = 0
    integer :: pass_count = 0
    
    write(*, '(A)') "Testing CORRECT Conductivity K-Matrix Physics"
    write(*, '(A)') "=============================================="
    
    ! Test 1: Z-function integration - must use plasma_z_function, not approximation
    call test_z_function_integration()
    
    ! Test 2: Complete tensor structure - all 9 elements non-zero in general
    call test_complete_tensor_structure()
    
    ! Test 3: Cyclotron resonance physics - proper frequency dependence
    call test_cyclotron_resonance_physics()
    
    ! Test 4: Species separation - electron vs ion contributions
    call test_species_separation()
    
    ! Test 5: FLRE contributions - proper Bessel function usage
    call test_flre_bessel_functions()
    
    ! Test 6: Stix tensor elements - verify known plasma physics relationships
    call test_stix_tensor_relationships()
    
    ! Test 7: Parallel vs perpendicular elements - K_zz physics
    call test_parallel_conductivity()
    
    ! Test 8: Known limits - cold plasma, high frequency limits
    call test_known_physics_limits()
    
    ! Summary
    write(*, '(A)') ""
    write(*, '(A, I0, A, I0, A)') "Results: ", pass_count, "/", test_count, " tests passed"
    
    if (pass_count == test_count) then
        write(*, '(A)') "SUCCESS: All conductivity physics tests passed!"
        stop 0
    else
        write(*, '(A)') "FAILURE: Conductivity physics implementation is broken!"
        stop 1
    end if
    
contains

    subroutine test_z_function_integration()
        ! Test that K-matrix calculation uses proper plasma Z-function
        type(cond_params_t) :: params
        real(real64) :: k_real_zz, k_imag_zz, k_real_xx, k_imag_xx
        integer :: ierr
        real(real64) :: zeta_expected
        complex(real64) :: z_expected
        
        call increment_test("Z-function integration in K-matrix")
        
        ! Setup test parameters with known values
        call setup_test_params(params)
        params%omega = 1.0e8_real64      ! Hz
        params%k_par = 100.0_real64      ! m^-1
        
        ! Calculate K_zz which should use Z-function: K_zz = (ω_p²/ω) * [1 + ζ*Z(ζ)]
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 2, 2, params, k_real_zz, k_imag_zz, ierr)
        
        ! Calculate K_xx which should use Z'(ζ): K_xx = (ω_p²/ω) * Γ_n * Z'(ζ)  
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 0, 0, params, k_real_xx, k_imag_xx, ierr)
        
        ! Test that both elements are non-zero (indicating Z-function is being used)
        if (abs(k_real_zz) < 1.0e-15_real64 .and. abs(k_imag_zz) < 1.0e-15_real64) then
            write(*, '(A)') "FAIL: K_zz must use plasma Z-function for parallel conductivity"
            return
        end if
        
        if (abs(k_real_xx) < 1.0e-15_real64 .and. abs(k_imag_xx) < 1.0e-15_real64) then
            write(*, '(A)') "FAIL: K_xx must use plasma Z-function derivative"
            return
        end if
        
        ! Test that K_zz and K_xx are different (they use different Z-function forms)
        if (abs(k_real_zz - k_real_xx) < 1.0e-10_real64) then
            write(*, '(A)') "FAIL: K_zz and K_xx must use different Z-function contributions"
            return
        end if
        
        call increment_pass("Z-function integration")
    end subroutine

    subroutine test_complete_tensor_structure()
        ! Test that all 9 tensor elements are properly calculated
        type(cond_params_t) :: params
        real(real64) :: k_elements(0:2, 0:2, 2)  ! [i,j,real/imag]
        integer :: ierr, i, j
        logical :: has_nonzero_offdiagonal
        
        call increment_test("Complete 3x3 tensor structure")
        
        call setup_test_params(params)
        
        ! Calculate all 9 tensor elements
        do i = 0, 2
            do j = 0, 2
                call calc_single_K_element(0.5_real64, 0, 0, 0, 0, i, j, params, &
                                          k_elements(i,j,1), k_elements(i,j,2), ierr)
            end do
        end do
        
        ! Check that off-diagonal elements are non-zero (current implementation sets many to 0)
        has_nonzero_offdiagonal = .false.
        if (abs(k_elements(0,2,1)) > 1.0e-10_real64) has_nonzero_offdiagonal = .true.  ! K_xz
        if (abs(k_elements(1,2,1)) > 1.0e-10_real64) has_nonzero_offdiagonal = .true.  ! K_yz
        if (abs(k_elements(2,0,1)) > 1.0e-10_real64) has_nonzero_offdiagonal = .true.  ! K_zx
        if (abs(k_elements(2,1,1)) > 1.0e-10_real64) has_nonzero_offdiagonal = .true.  ! K_zy
        
        ! This will FAIL because current implementation sets these to 0.0
        if (.not. has_nonzero_offdiagonal) then
            write(*, '(A)') "FAIL: Off-diagonal elements K_xz, K_yz, K_zx, K_zy incorrectly set to zero"
            return
        end if
        
        call increment_pass("Complete tensor structure")
    end subroutine

    subroutine test_cyclotron_resonance_physics()
        ! Test proper cyclotron frequency dependence
        type(cond_params_t) :: params
        real(real64) :: k_real_low, k_imag_low, k_real_high, k_imag_high
        integer :: ierr
        
        call increment_test("Cyclotron resonance physics")
        
        call setup_test_params(params)
        
        ! Test at low frequency (well below cyclotron)
        params%omega = 1.0e6_real64  ! Much less than typical cyclotron frequency
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 0, 1, params, &
                                  k_real_low, k_imag_low, ierr)
        
        ! Test at high frequency (well above cyclotron) 
        params%omega = 1.0e12_real64  ! Much greater than cyclotron frequency
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 0, 1, params, &
                                  k_real_high, k_imag_high, ierr)
        
        ! K_xy should have different behavior at different frequencies
        ! Current implementation uses hard-coded omega_test = 1.0e8, so this will FAIL
        if (abs(k_real_high - k_real_low) < 1.0e-10_real64) then
            write(*, '(A)') "FAIL: K_xy must depend on actual wave frequency, not hard-coded value"
            return
        end if
        
        call increment_pass("Cyclotron resonance physics")
    end subroutine

    subroutine test_species_separation()
        ! Test that electron and ion contributions are properly separated
        type(cond_params_t) :: params
        real(real64) :: k_electron, k_ion, k_dummy
        integer :: ierr
        
        call increment_test("Species separation (electron vs ion)")
        
        call setup_test_params(params)
        
        ! Get electron contribution (spec = 0)
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 0, 0, params, &
                                  k_electron, k_dummy, ierr)
        
        ! Get ion contribution (spec = 1) 
        call calc_single_K_element(0.5_real64, 1, 0, 0, 0, 0, 0, params, &
                                  k_ion, k_dummy, ierr)
        
        ! Should be different due to mass difference
        if (abs(k_electron - k_ion) < 1.0e-10_real64) then
            write(*, '(A)') "FAIL: Electron and ion contributions must be different"
            return  
        end if
        
        call increment_pass("Species separation")
    end subroutine

    subroutine test_flre_bessel_functions()
        ! Test that FLRE uses proper Bessel functions, not factorial approximation
        type(cond_params_t) :: params
        real(real64) :: k_00, k_11, k_22, k_dummy
        integer :: ierr
        
        call increment_test("FLRE Bessel function contributions")
        
        call setup_test_params(params)
        
        ! Test different FLRE orders (p,q combinations)
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 0, 0, params, k_00, k_dummy, ierr)  ! p=0,q=0
        call calc_single_K_element(0.5_real64, 0, 0, 1, 1, 0, 0, params, k_11, k_dummy, ierr)  ! p=1,q=1  
        call calc_single_K_element(0.5_real64, 0, 0, 2, 2, 0, 0, params, k_22, k_dummy, ierr)  ! p=2,q=2
        
        ! Current implementation uses factorial approximation, should use proper Bessel functions
        ! This creates wrong scaling between different orders
        if (abs(k_11/k_00 - 0.5_real64) < 0.1_real64) then
            write(*, '(A)') "FAIL: FLRE must use proper Bessel functions, not factorial approximation"
            return
        end if
        
        call increment_pass("FLRE Bessel functions")
    end subroutine

    subroutine test_stix_tensor_relationships()
        ! Test known relationships in Stix conductivity tensor
        type(cond_params_t) :: params
        real(real64) :: k_xx, k_yy, k_xy, k_yx, k_dummy
        integer :: ierr
        
        call increment_test("Stix tensor relationships")
        
        call setup_test_params(params)
        
        ! Get tensor elements
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 0, 0, params, k_xx, k_dummy, ierr)
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 1, 1, params, k_yy, k_dummy, ierr)
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 0, 1, params, k_xy, k_dummy, ierr)
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 1, 0, params, k_yx, k_dummy, ierr)
        
        ! In cold plasma limit: K_xx = K_yy and K_yx = -K_xy
        ! Current implementation doesn't enforce these relationships
        if (abs(k_xx - k_yy) > 1.0e-6_real64 * abs(k_xx)) then
            write(*, '(A)') "FAIL: K_xx should equal K_yy in isotropic plasma"
            return
        end if
        
        call increment_pass("Stix tensor relationships")
    end subroutine

    subroutine test_parallel_conductivity()
        ! Test K_zz (parallel) conductivity physics
        type(cond_params_t) :: params
        real(real64) :: k_zz, k_dummy
        integer :: ierr
        
        call increment_test("Parallel conductivity K_zz physics")
        
        call setup_test_params(params)
        
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 2, 2, params, k_zz, k_dummy, ierr)
        
        ! K_zz should include proper Z-function contribution: (1 + Z_func_re)
        ! Current implementation uses wrong Z-function approximation
        if (abs(k_zz) < 1.0e-10_real64) then
            write(*, '(A)') "FAIL: K_zz must include proper plasma Z-function contribution"
            return
        end if
        
        call increment_pass("Parallel conductivity")
    end subroutine

    subroutine test_known_physics_limits()
        ! Test known physics limits (cold plasma, high frequency)
        type(cond_params_t) :: params
        real(real64) :: k_cold, k_hot, k_dummy
        integer :: ierr
        
        call increment_test("Known physics limits")
        
        ! Cold plasma limit (very low temperature)
        call setup_test_params(params)
        params%temp_e = 1.0_real64  ! Very cold
        params%temp_i = 1.0_real64
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 0, 0, params, k_cold, k_dummy, ierr)
        
        ! Hot plasma (high temperature)
        params%temp_e = 1.0e6_real64  ! Very hot
        params%temp_i = 1.0e6_real64
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 0, 0, params, k_hot, k_dummy, ierr)
        
        ! Should show temperature dependence through Z-function
        if (abs(k_hot - k_cold) < 1.0e-10_real64) then
            write(*, '(A)') "FAIL: Conductivity must depend on temperature through Z-function"
            return
        end if
        
        call increment_pass("Known physics limits")
    end subroutine

    subroutine setup_test_params(params)
        type(cond_params_t), intent(out) :: params
        
        ! Realistic plasma parameters for testing
        params%density_e = 1.0e19_real64     ! m^-3
        params%density_i = 1.0e19_real64     ! m^-3
        params%temp_e = 1.0e4_real64         ! eV
        params%temp_i = 1.0e4_real64         ! eV
        params%mass_e = 9.109e-31_real64     ! kg
        params%mass_i = 1.673e-27_real64     ! kg (proton)
        params%charge_e = -1.602e-19_real64  ! C
        params%charge_i = 1.602e-19_real64   ! C
        params%B0 = 2.0_real64               ! T
        params%coll_freq_e = 1.0e6_real64    ! Hz
        params%coll_freq_i = 1.0e5_real64    ! Hz
        
        ! Wave parameters (should come from mode solver)
        params%omega = 1.0e8_real64          ! rad/s
        params%k_perp = 100.0_real64         ! m^-1
        params%k_par = 10.0_real64           ! m^-1
    end subroutine

    subroutine increment_test(test_name)
        character(len=*), intent(in) :: test_name
        test_count = test_count + 1
        write(*, '(A, I0, A, A, A)', advance='no') "Test ", test_count, ": ", test_name, " ... "
    end subroutine

    subroutine increment_pass(test_name)
        character(len=*), intent(in) :: test_name
        pass_count = pass_count + 1
        write(*, '(A)') "PASS"
    end subroutine

end program test_conductivity_physics_correct