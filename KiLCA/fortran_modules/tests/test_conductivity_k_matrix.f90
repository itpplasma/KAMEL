program test_conductivity_k_matrix
    use iso_fortran_env, only: real64
    use kilca_conductivity_m
    implicit none
    
    integer :: test_status = 0
    real(real64), parameter :: tolerance = 1.0e-8_real64
    
    print *, "===========================================" 
    print *, "Testing Conductivity K-Matrix Calculations [GREEN PHASE]"
    print *, "==========================================="
    
    call test_k_matrix_structure()
    call test_k_matrix_physics()
    call test_k_matrix_splines()
    call test_k_matrix_evaluation()
    
    if (test_status == 0) then
        print *, ""
        print *, "Conductivity K-matrix tests PASSED - GREEN phase successful"
        stop 0
    else
        print *, ""
        print *, "Conductivity K-matrix tests FAILED:", test_status, "test(s) failed"
        stop 1
    end if

contains

    !> Test K-matrix data structure and array indexing
    subroutine test_k_matrix_structure()
        type(cond_profiles_t) :: profiles
        type(cond_params_t) :: params
        integer :: ierr
        integer :: flreo = 2, dimt = 2, dimx = 10
        integer :: index1, index2, expected_dim
        
        print *, "Testing K-matrix structure and indexing..."
        print *, ""
        
        ! Test 1: Profile creation and dimension calculation
        call cond_profiles_create(profiles, flreo, dimt, dimx, .false., ierr)
        call check_error("Profile creation", ierr, 0)
        
        ! Test 2: Verify array dimensions
        ! K matrix: spec(2) * type(2) * p(3) * q(3) * i(3) * j(3) * parts(2) * nodes(10)
        expected_dim = 2 * 2 * 3 * 3 * 3 * 3 * 2 * 10  ! = 3240
        call check_integer("K array dimension", profiles%dim_K_array, expected_dim)
        
        ! Test 3: Index function consistency
        index1 = iKs_index(0, 0, 0, 0, 0, 0, 0, 0, flreo, dimt, dimx)
        index2 = iKa_index(0, 0, 0, 0, 0, 0, 0, 0, flreo, dimt, dimx)
        call check_integer("Index function consistency", index1, index2)
        
        ! Test 4: Index bounds checking
        index1 = iKs_index(1, 1, 2, 2, 2, 2, 1, 9, flreo, dimt, dimx)
        if (index1 < 1 .or. index1 > profiles%dim_K_array) then
            print *, "FAIL: Index out of bounds: ", index1
            test_status = test_status + 1
        else
            print *, "PASS: Index bounds check"
        end if
        
        call cond_profiles_destroy(profiles, ierr)
        
    end subroutine test_k_matrix_structure
    
    !> Test K-matrix physics calculations
    subroutine test_k_matrix_physics()
        type(cond_profiles_t) :: profiles
        type(cond_params_t) :: params
        integer :: ierr
        integer :: flreo = 1, dimt = 1, dimx = 5
        real(real64) :: k_real, k_imag
        real(real64) :: expected_plasma_freq, expected_cyclotron_freq
        real(real64) :: k_real_e, k_imag_e, k_real_i, k_imag_i
        real(real64) :: k_real_p0, k_real_p1
        real(real64) :: k_diag, k_offdiag
        
        print *, ""
        print *, "Testing K-matrix plasma physics calculations..."
        print *, ""
        
        ! Setup test parameters (typical tokamak plasma)
        params%B0 = 2.0_real64                    ! Tesla
        params%density_e = 1.0e19_real64          ! m^-3
        params%density_i = 1.0e19_real64          ! m^-3
        params%temp_e = 1000.0_real64 * 1.602e-19_real64  ! 1 keV in Joules
        params%temp_i = 1000.0_real64 * 1.602e-19_real64  ! 1 keV in Joules
        params%mass_e = 9.109e-31_real64          ! kg
        params%mass_i = 1.673e-27_real64          ! kg (proton)
        params%charge_e = -1.602e-19_real64       ! Coulomb
        params%charge_i = 1.602e-19_real64        ! Coulomb
        params%coll_freq_e = 1.0e4_real64         ! Hz
        params%coll_freq_i = 1.0e3_real64         ! Hz
        
        call cond_profiles_create(profiles, flreo, dimt, dimx, .false., ierr)
        call check_error("Profile creation for physics test", ierr, 0)
        
        ! Test 1: Electron cyclotron frequency calculation
        expected_cyclotron_freq = abs(params%charge_e) * params%B0 / params%mass_e
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 0, 0, params, k_real, k_imag, ierr)
        call check_error("K element calculation", ierr, 0)
        
        ! Test 2: Plasma frequency scaling
        ! K matrix should scale with density/temperature
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 0, 0, params, k_real, k_imag, ierr)
        if (abs(k_real) < 1.0e-20_real64) then
            print *, "FAIL: K matrix element too small - no plasma response: ", k_real
            test_status = test_status + 1
        else
            print *, "PASS: K matrix has finite plasma response: ", k_real
        end if
        
        ! Test 3: Species dependence (electrons vs ions)
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 0, 0, params, k_real_e, k_imag_e, ierr)
        call calc_single_K_element(0.5_real64, 1, 0, 0, 0, 0, 0, params, k_real_i, k_imag_i, ierr)
        
        ! Electrons should have different response than ions due to mass difference
        if (abs(k_real_e - k_real_i) < tolerance) then
            print *, "FAIL: No species dependence in K matrix: e=", k_real_e, " i=", k_real_i
            test_status = test_status + 1
        else
            print *, "PASS: Species dependence exists in K matrix"
        end if
        
        ! Test 4: Finite Larmor radius effects
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 0, 0, params, k_real_p0, k_imag, ierr)
        call calc_single_K_element(0.5_real64, 0, 0, 1, 0, 0, 0, params, k_real_p1, k_imag, ierr)
        
        ! Higher FLRE orders should show different behavior
        if (abs(k_real_p0 - k_real_p1) < tolerance) then
            print *, "FAIL: No FLRE order dependence: p=0:", k_real_p0, " p=1:", k_real_p1
            test_status = test_status + 1
        else
            print *, "PASS: FLRE order dependence exists"
        end if
        
        ! Test 5: Matrix structure (diagonal vs off-diagonal)
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 0, 0, params, k_diag, k_imag, ierr)
        call calc_single_K_element(0.5_real64, 0, 0, 0, 0, 0, 1, params, k_offdiag, k_imag, ierr)
        
        ! Diagonal elements should typically be larger than off-diagonal
        if (abs(k_diag) <= abs(k_offdiag)) then
            print *, "FAIL: Diagonal element not dominant: diag=", k_diag, " off=", k_offdiag
            test_status = test_status + 1
        else
            print *, "PASS: Diagonal element dominance"
        end if
        
        call cond_profiles_destroy(profiles, ierr)
        
    end subroutine test_k_matrix_physics
    
    !> Test K-matrix spline calculation and interpolation
    subroutine test_k_matrix_splines()
        type(cond_profiles_t) :: profiles
        type(cond_params_t) :: params
        integer :: ierr
        integer :: flreo = 1, dimt = 1, dimx = 5
        real(real64), allocatable :: k_values(:)
        integer :: i
        logical :: spline_varies
        real(real64), dimension(3) :: test_values
        
        print *, ""
        print *, "Testing K-matrix spline calculations..."
        print *, ""
        
        ! Setup minimal parameters
        params%B0 = 1.0_real64
        params%density_e = 1.0e19_real64
        params%density_i = 1.0e19_real64
        params%temp_e = 1.602e-16_real64  ! 1 keV
        params%temp_i = 1.602e-16_real64  ! 1 keV
        params%mass_e = 9.109e-31_real64
        params%mass_i = 1.673e-27_real64  ! proton
        params%charge_e = -1.602e-19_real64
        params%charge_i = 1.602e-19_real64
        params%coll_freq_e = 1.0e4_real64
        params%coll_freq_i = 1.0e3_real64
        
        call cond_profiles_create(profiles, flreo, dimt, dimx, .false., ierr)
        call check_error("Profile creation for spline test", ierr, 0)
        
        ! Test 1: Generate K matrix profiles
        call sample_cond_func(profiles, params, ierr)
        call check_error("K matrix sampling", ierr, 0)
        
        ! Test 2: Verify K matrix values are calculated
        if (all(abs(profiles%K) < 1.0e-20_real64)) then
            print *, "FAIL: All K matrix values are zero or too small"
            test_status = test_status + 1
        else
            print *, "PASS: K matrix values calculated"
        end if
        
        ! Test 3: Calculate splines
        call calc_splines_for_K(profiles, ierr)
        call check_error("K matrix spline calculation", ierr, 0)
        
        ! Test 4: Verify spline calculation actually worked
        if (profiles%sidK == 0) then
            print *, "FAIL: Spline ID not set after calculation"
            test_status = test_status + 1
        else
            print *, "PASS: Spline calculation completed"
        end if
        
        ! Test 5: Test spline evaluation at different radii
        allocate(k_values(profiles%dim_K_array))
        
        call eval_K_matrices(profiles, 0.0_real64, k_values, ierr)
        call check_error("K matrix evaluation at r=0", ierr, 0)
        
        call eval_K_matrices(profiles, 0.5_real64, k_values, ierr)
        call check_error("K matrix evaluation at r=0.5", ierr, 0)
        
        call eval_K_matrices(profiles, 1.0_real64, k_values, ierr)
        call check_error("K matrix evaluation at r=1", ierr, 0)
        
        ! Test 6: Verify interpolation provides smooth variation
        ! This test checks if the spline evaluation actually interpolates
        spline_varies = .false.
        call eval_K_matrices(profiles, 0.0_real64, k_values, ierr)
        test_values(1) = k_values(1)
        call eval_K_matrices(profiles, 0.5_real64, k_values, ierr)
        test_values(2) = k_values(1)
        call eval_K_matrices(profiles, 1.0_real64, k_values, ierr)
        test_values(3) = k_values(1)
        
        ! Check if interpolation gives smooth intermediate values
        if (test_values(2) /= test_values(1) .or. test_values(2) /= test_values(3)) then
            spline_varies = .true.
        end if
        
        ! In current placeholder implementation, this will fail (all return same values)
        if (.not. spline_varies) then
            print *, "FAIL: Spline interpolation not working - placeholder implementation"
            test_status = test_status + 1
        else
            print *, "PASS: Spline interpolation provides smooth variation"
        end if
        
        deallocate(k_values)
        call cond_profiles_destroy(profiles, ierr)
        
    end subroutine test_k_matrix_splines
    
    !> Test K-matrix evaluation and C-matrix derivation
    subroutine test_k_matrix_evaluation()
        type(cond_profiles_t) :: profiles
        type(cond_params_t) :: params
        integer :: ierr
        integer :: flreo = 1, dimt = 1, dimx = 5
        real(real64), allocatable :: k_values(:), c_values(:)
        real(real64) :: c_sum
        
        print *, ""
        print *, "Testing K-matrix evaluation and C-matrix derivation..."
        print *, ""
        
        ! Setup parameters
        params%B0 = 1.0_real64
        params%density_e = 1.0e19_real64
        params%density_i = 1.0e19_real64
        params%temp_e = 1.602e-16_real64
        params%temp_i = 1.602e-16_real64
        params%mass_e = 9.109e-31_real64
        params%mass_i = 1.673e-27_real64
        params%charge_e = -1.602e-19_real64
        params%charge_i = 1.602e-19_real64
        params%coll_freq_e = 1.0e4_real64
        params%coll_freq_i = 1.0e3_real64
        
        call cond_profiles_create(profiles, flreo, dimt, dimx, .false., ierr)
        call check_error("Profile creation for evaluation test", ierr, 0)
        
        ! Test 1: Full K-matrix calculation workflow
        call sample_cond_func(profiles, params, ierr)
        call check_error("K matrix sampling", ierr, 0)
        
        call calc_splines_for_K(profiles, ierr)
        call check_error("K matrix spline calculation", ierr, 0)
        
        ! Test 2: C-matrix derivation from K-matrices
        call calc_C_matrices(profiles, ierr)
        call check_error("C matrix calculation", ierr, 0)
        
        ! Test 3: Verify C-matrix has non-zero values
        if (all(abs(profiles%C) < 1.0e-20_real64)) then
            print *, "FAIL: All C matrix values are zero - derivation failed"
            test_status = test_status + 1
        else
            print *, "PASS: C matrix derived from K matrix"
        end if
        
        ! Test 4: C-matrix spline calculation
        call calc_splines_for_C(profiles, ierr)
        call check_error("C matrix spline calculation", ierr, 0)
        
        ! Test 5: Evaluate C-matrices
        allocate(c_values(profiles%dim_C_array))
        call eval_C_matrices(profiles, 0.5_real64, c_values, ierr)
        call check_error("C matrix evaluation", ierr, 0)
        
        ! Test 6: Conservation properties
        ! Sum of certain C-matrix elements should satisfy physics constraints
        ! This is a placeholder for more detailed physics validation
        c_sum = sum(abs(c_values))
        if (c_sum < 1.0e-20_real64) then
            print *, "FAIL: C matrix evaluation returned zero values"
            test_status = test_status + 1
        else
            print *, "PASS: C matrix evaluation has finite values: ", c_sum
        end if
        
        deallocate(c_values)
        call cond_profiles_destroy(profiles, ierr)
        
    end subroutine test_k_matrix_evaluation
    
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
    
    !> Check integer values
    subroutine check_integer(test_name, actual, expected)
        character(len=*), intent(in) :: test_name
        integer, intent(in) :: actual, expected
        
        if (actual /= expected) then
            print *, "FAIL: ", test_name, " = ", actual, " (expected ", expected, ")"
            test_status = test_status + 1
        else
            print *, "PASS: ", test_name, " = ", actual
        end if
    end subroutine check_integer

end program test_conductivity_k_matrix