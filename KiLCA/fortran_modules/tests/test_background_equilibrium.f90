program test_background_equilibrium
    use iso_fortran_env, only: real64
    use kilca_background_m
    use kilca_settings_m
    implicit none
    
    integer :: test_status = 0
    real(real64), parameter :: tolerance = 1.0e-8_real64
    
    print *, "===========================================" 
    print *, "Testing Background Equilibrium Calculations [RED PHASE - SHOULD FAIL]"
    print *, "==========================================="
    
    call test_background_structure()
    call test_equilibrium_physics()
    call test_profile_interpolation()
    call test_f0_parameters()
    
    if (test_status == 0) then
        print *, ""
        print *, "Background equilibrium tests PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "Background equilibrium tests FAILED:", test_status, "test(s) failed (expected in RED phase)"
        stop 1
    end if

contains

    !> Test background data structure and initialization
    subroutine test_background_structure()
        type(background_t) :: bg
        type(settings_t), pointer :: settings
        integer :: ierr, i
        integer :: expected_indices(5)
        
        print *, "Testing background equilibrium structure..."
        print *, ""
        
        ! Test 1: Background creation and settings link
        call settings_create(settings, "test_path", ierr)
        call check_error("Settings creation", ierr, 0)
        
        call background_create(bg, settings, ierr)
        call check_error("Background creation", ierr, 0)
        
        ! Test 2: Verify basic structure
        if (.not. associated(bg%sd)) then
            print *, "FAIL: Settings pointer not associated"
            test_status = test_status + 1
        else
            print *, "PASS: Settings linked to background"
        end if
        
        ! Test 3: Verify profile indices are set
        call background_set_profiles_indices(bg)
        
        ! Check key indices are assigned properly
        expected_indices = [bg%i_q, bg%i_n, bg%i_Er, bg%i_B, bg%i_dPhi0]
        if (any(expected_indices < 0)) then
            print *, "FAIL: Profile indices not set properly"
            test_status = test_status + 1
        else
            print *, "PASS: Profile indices assigned"
        end if
        
        ! Test 4: Check array allocations
        if (.not. allocated(bg%i_T) .or. .not. allocated(bg%i_Vth)) then
            print *, "FAIL: Species-dependent arrays not allocated"
            test_status = test_status + 1
        else
            print *, "PASS: Species arrays allocated"
        end if
        
        ! Test 5: Memory management
        call background_destroy(bg)
        
        call settings_destroy(settings, ierr)
        
    end subroutine test_background_structure
    
    !> Test equilibrium physics calculations
    subroutine test_equilibrium_physics()
        type(background_t) :: bg
        type(settings_t), pointer :: settings
        integer :: ierr
        integer :: dimx = 50
        real(real64) :: B0_test = 2.0_real64  ! Tesla
        real(real64) :: rtor_test = 1.0_real64  ! meters
        real(real64) :: r_test = 0.5_real64
        real(real64) :: q_val, n_val, Ti_val, Te_val, Vth_val, Vz_val, dPhi0_val
        integer :: i
        
        print *, ""
        print *, "Testing equilibrium physics calculations..."
        print *, ""
        
        ! Setup test configuration
        call settings_create(settings, "test_path", ierr)
        call check_error("Settings creation for physics test", ierr, 0)
        
        ! Set magnetic field parameters
        settings%background_settings%B0 = B0_test
        settings%background_settings%rtor = rtor_test
        
        call background_create(bg, settings, ierr)
        call check_error("Background creation for physics test", ierr, 0)
        
        ! Allocate profile arrays
        allocate(bg%x(dimx))
        allocate(bg%y(dimx * bg%dimy))
        bg%dimx = dimx
        
        ! Generate radial grid
        do i = 1, dimx
            bg%x(i) = real(i-1, real64) / real(dimx-1, real64)
        end do
        
        ! Test 1: Load input profiles
        call background_load_input_profiles(bg, ierr)
        call check_error("Load input profiles", ierr, 0)
        
        ! Test 2: Verify grid is generated
        if (bg%x(1) /= 0.0_real64 .or. abs(bg%x(dimx) - 1.0_real64) > tolerance) then
            print *, "FAIL: Radial grid not properly generated"
            test_status = test_status + 1
        else
            print *, "PASS: Radial grid generated"
        end if
        
        ! Test 3: Calculate equilibrium fields
        call background_calculate_equilibrium(bg, ierr)
        call check_error("Calculate equilibrium", ierr, 0)
        
        ! Test 4: Verify magnetic field calculations
        ! Check that B_z field is set to B0
        if (bg%i_Bz >= 0) then
            ! Get field at middle of domain
            i = dimx / 2
            if (abs(bg%y(bg%i_Bz * dimx + i) - B0_test) > tolerance) then
                print *, "FAIL: B_z field not set to B0:", bg%y(bg%i_Bz * dimx + i), "expected:", B0_test
                test_status = test_status + 1
            else
                print *, "PASS: B_z field correctly set to B0"
            end if
        end if
        
        ! Test 5: Check metric coefficients
        if (bg%i_hth >= 0) then
            i = dimx / 2
            ! h_theta should equal r in cylindrical coordinates
            if (abs(bg%y(bg%i_hth * dimx + i) - bg%x(i)) > tolerance) then
                print *, "FAIL: h_theta metric incorrect"
                test_status = test_status + 1
            else
                print *, "PASS: h_theta metric correctly calculated"
            end if
        end if
        
        ! Test 6: Spline profile calculations
        call background_spline_profiles(bg, ierr)
        call check_error("Background spline profiles", ierr, 0)
        
        ! Test 7: F0 parameter calculation
        call background_find_f0_parameters(bg, ierr)
        call check_error("F0 parameter calculation", ierr, 0)
        
        ! Test 8: Profile interpolation
        call background_interp_basic_profiles(bg, r_test, q_val, n_val, Ti_val, Te_val, Vth_val, Vz_val, dPhi0_val)
        
        ! Check that interpolated values are physical
        if (q_val <= 0.0_real64) then
            print *, "FAIL: Safety factor q <= 0:", q_val
            test_status = test_status + 1
        else
            print *, "PASS: Safety factor q is positive:", q_val
        end if
        
        if (n_val <= 0.0_real64) then
            print *, "FAIL: Density n <= 0:", n_val
            test_status = test_status + 1
        else
            print *, "PASS: Density is positive:", n_val
        end if
        
        if (Ti_val <= 0.0_real64 .or. Te_val <= 0.0_real64) then
            print *, "FAIL: Temperatures not positive: Ti=", Ti_val, "Te=", Te_val
            test_status = test_status + 1
        else
            print *, "PASS: Temperatures are positive"
        end if
        
        call background_destroy(bg)
        call settings_destroy(settings, ierr)
        
    end subroutine test_equilibrium_physics
    
    !> Test profile interpolation accuracy
    subroutine test_profile_interpolation()
        type(background_t) :: bg
        type(settings_t), pointer :: settings
        integer :: ierr, i
        integer :: dimx = 10
        real(real64) :: r1, r2, r3
        real(real64) :: q1, q2, q3, n1, n2, n3, Ti1, Ti2, Ti3, Te1, Te2, Te3
        real(real64) :: Vth1, Vth2, Vth3, Vz1, Vz2, Vz3, dPhi01, dPhi02, dPhi03
        logical :: interpolation_varies
        
        print *, ""
        print *, "Testing profile interpolation accuracy..."
        print *, ""
        
        call settings_create(settings, "test_path", ierr)
        call background_create(bg, settings, ierr)
        
        ! Set profile indices first
        call background_set_profiles_indices(bg)
        
        ! Setup profiles with known values
        allocate(bg%x(dimx))
        allocate(bg%y(dimx * bg%dimy))
        bg%dimx = dimx
        
        ! Create test radial grid
        do i = 1, dimx
            bg%x(i) = real(i-1, real64) / real(dimx-1, real64)
        end do
        
        call background_load_input_profiles(bg, ierr)
        call background_calculate_equilibrium(bg, ierr)
        call background_spline_profiles(bg, ierr)
        
        ! Test interpolation at different radial positions
        r1 = 0.0_real64  ! Edge
        r2 = 0.5_real64  ! Middle
        r3 = 1.0_real64  ! Center
        
        call background_interp_basic_profiles(bg, r1, q1, n1, Ti1, Te1, Vth1, Vz1, dPhi01)
        call background_interp_basic_profiles(bg, r2, q2, n2, Ti2, Te2, Vth2, Vz2, dPhi02)
        call background_interp_basic_profiles(bg, r3, q3, n3, Ti3, Te3, Vth3, Vz3, dPhi03)
        
        ! Test 1: Check if interpolation varies with radius
        interpolation_varies = .false.
        if (abs(q1 - q2) > tolerance .or. abs(q2 - q3) > tolerance) interpolation_varies = .true.
        if (abs(n1 - n2) > tolerance .or. abs(n2 - n3) > tolerance) interpolation_varies = .true.
        
        ! In current placeholder implementation, all values are constant
        if (.not. interpolation_varies) then
            print *, "FAIL: Profile interpolation returns constant values - placeholder implementation"
            test_status = test_status + 1
        else
            print *, "PASS: Profile interpolation varies with radius"
        end if
        
        ! Test 2: Check boundary behavior
        ! Values at boundaries should be reasonable
        if (q1 < 0.1_real64 .or. q1 > 10.0_real64) then
            print *, "FAIL: Safety factor at boundary unphysical:", q1
            test_status = test_status + 1
        else
            print *, "PASS: Safety factor at boundary is physical"
        end if
        
        ! Test 3: Check conservation properties
        ! For now, just check that values don't blow up
        if (n2 > 1.0e22_real64 .or. Ti2 > 1.0e5_real64) then
            print *, "FAIL: Interpolated values too large - numerical instability"
            test_status = test_status + 1
        else
            print *, "PASS: Interpolated values are bounded"
        end if
        
        call background_destroy(bg)
        call settings_destroy(settings, ierr)
        
    end subroutine test_profile_interpolation
    
    !> Test F0 distribution function parameters
    subroutine test_f0_parameters()
        type(background_t) :: bg
        type(settings_t), pointer :: settings
        integer :: ierr
        integer :: dimx = 20
        real(real64) :: dPhi0_before, dPhi0_after
        
        print *, ""
        print *, "Testing F0 distribution function parameters..."
        print *, ""
        
        call settings_create(settings, "test_path", ierr)
        call background_create(bg, settings, ierr)
        
        ! Set profile indices first
        call background_set_profiles_indices(bg)
        
        allocate(bg%x(dimx))
        allocate(bg%y(dimx * bg%dimy))
        bg%dimx = dimx
        
        call background_load_input_profiles(bg, ierr)
        call background_calculate_equilibrium(bg, ierr)
        
        ! Test 1: Check initial state of F0 parameters
        if (bg%i_dPhi0 >= 0) then
            dPhi0_before = bg%y(bg%i_dPhi0 * dimx + 1)
        else
            dPhi0_before = 0.0_real64
        end if
        
        ! Test 2: Calculate F0 parameters
        call background_find_f0_parameters(bg, ierr)
        call check_error("F0 parameter calculation", ierr, 0)
        
        ! Test 3: Verify parameters have been set
        if (bg%i_dPhi0 >= 0) then
            dPhi0_after = bg%y(bg%i_dPhi0 * dimx + 1)
            print *, "PASS: dPhi0 parameter calculated:", dPhi0_after
        else
            print *, "FAIL: dPhi0 index not set"
            test_status = test_status + 1
        end if
        
        ! Test 4: Check parameter self-consistency
        ! F0 parameters should satisfy quasi-neutrality and current conservation
        ! This is a placeholder test - full implementation would check physics constraints
        if (bg%flag_dPhi0_calc /= 0) then
            print *, "PASS: F0 calculation flag set"
        else
            print *, "FAIL: F0 calculation flag not set - may indicate placeholder implementation"
            test_status = test_status + 1
        end if
        
        ! Test 5: Check species-dependent parameters
        if (allocated(bg%i_n_p) .and. allocated(bg%i_Vp_p) .and. allocated(bg%i_Vt_p)) then
            print *, "PASS: Species-dependent F0 parameters allocated"
        else
            print *, "FAIL: Species-dependent F0 parameters not allocated"
            test_status = test_status + 1
        end if
        
        call background_destroy(bg)
        call settings_destroy(settings, ierr)
        
    end subroutine test_f0_parameters
    
    !> Check error codes
    subroutine check_error(test_name, actual, expected)
        character(len=*), intent(in) :: test_name
        integer, intent(in) :: actual, expected
        
        if (actual /= expected) then
            print *, "FAIL: ", test_name, " returned error ", actual, " (expected ", expected, ")"
            test_status = test_status + 1
        else
            print *, "PASS: ", test_name
        end if
    end subroutine check_error

end program test_background_equilibrium