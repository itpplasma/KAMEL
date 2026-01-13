program test_rhs_balance_integration
    !
    ! Integration test for rhs_balance module refactoring
    !
    ! This test verifies that the refactored rhs_balance_m module produces
    ! identical output to the original rhs_balance.f90 implementation.
    !
    ! Strategy:
    !   1. Set up realistic test state with npoic=100
    !   2. Call original implementation (rhs_balance_orig)
    !   3. Save outputs (nz, nsize, irow, icol, amat, rhsvec)
    !   4. Reset state
    !   5. Call refactored implementation (rhs_balance from rhs_balance_m)
    !   6. Compare outputs with tolerance 1e-12
    !
    use QLBalance_kinds, only: dp
    use grid_mod
    use plasma_parameters
    use baseparam_mod
    use wave_code_data
    use matrix_mod
    use rhs_balance_m, only: rhs_balance, initialize_rhs

    implicit none

    ! Saved outputs from original implementation
    integer :: nz_orig, nsize_orig
    integer, allocatable :: irow_orig(:), icol_orig(:)
    real(dp), allocatable :: amat_orig(:), rhsvec_orig(:)

    ! Test state vectors (use different names to avoid conflict with grid_mod)
    real(dp), allocatable :: y_test(:), dy_test(:)

    ! Test parameters
    integer, parameter :: test_npoic = 100
    real(dp), parameter :: tol = 1.0d-12

    ! Test result tracking
    integer :: num_passed, num_failed
    logical :: test_ok

    num_passed = 0
    num_failed = 0

    print *, "========================================"
    print *, "  RHS Balance Integration Test"
    print *, "========================================"
    print *, ""
    print *, "Testing with npoic = ", test_npoic
    print *, ""

    ! Set up test state
    call setup_test_state(test_npoic)

    ! Allocate state vectors
    allocate(y_test(neqset), dy_test(neqset))

    ! Initialize y vector from params
    call pack_y_from_params(y_test)

    print *, "1. Testing initialize_rhs + rhs_balance..."
    print *, ""

    !=========================================================================
    ! Step 1: Run original implementation
    !=========================================================================
    print *, "   Running original implementation..."

    call initialize_rhs_orig(y_test, dy_test)
    call rhs_balance_orig(0.0_dp, y_test, dy_test)

    ! Save original outputs
    nz_orig = nz
    nsize_orig = nsize
    allocate(irow_orig(nz), icol_orig(nz), amat_orig(nz), rhsvec_orig(nsize))
    irow_orig = irow(1:nz)
    icol_orig = icol(1:nz)
    amat_orig = amat(1:nz)
    rhsvec_orig = rhsvec(1:nsize)

    print *, "   Original: nz = ", nz_orig, ", nsize = ", nsize_orig

    ! Deallocate matrix arrays for fresh start
    deallocate(irow, icol, amat, rhsvec)

    !=========================================================================
    ! Step 2: Reset state and run refactored implementation
    !=========================================================================
    print *, "   Running refactored implementation..."

    ! Reset module state that may have been modified
    call reset_module_state()
    call pack_y_from_params(y_test)

    call initialize_rhs(y_test, dy_test)
    call rhs_balance(0.0_dp, y_test, dy_test)

    print *, "   Refactored: nz = ", nz, ", nsize = ", nsize

    !=========================================================================
    ! Step 3: Compare outputs
    !=========================================================================
    print *, ""
    print *, "   Comparing outputs..."

    ! Compare nz
    if (nz /= nz_orig) then
        print *, "   FAIL: nz mismatch"
        print *, "     Original:   ", nz_orig
        print *, "     Refactored: ", nz
        num_failed = num_failed + 1
    else
        print *, "   PASS: nz matches (", nz, ")"
        num_passed = num_passed + 1
    end if

    ! Compare nsize
    if (nsize /= nsize_orig) then
        print *, "   FAIL: nsize mismatch"
        print *, "     Original:   ", nsize_orig
        print *, "     Refactored: ", nsize
        num_failed = num_failed + 1
    else
        print *, "   PASS: nsize matches (", nsize, ")"
        num_passed = num_passed + 1
    end if

    ! Compare irow
    call compare_int_arrays(irow_orig, irow(1:nz), "irow", test_ok)
    if (test_ok) then
        num_passed = num_passed + 1
    else
        num_failed = num_failed + 1
    end if

    ! Compare icol
    call compare_int_arrays(icol_orig, icol(1:nz), "icol", test_ok)
    if (test_ok) then
        num_passed = num_passed + 1
    else
        num_failed = num_failed + 1
    end if

    ! Compare amat
    call compare_real_arrays(amat_orig, amat(1:nz), "amat", tol, test_ok)
    if (test_ok) then
        num_passed = num_passed + 1
    else
        num_failed = num_failed + 1
    end if

    ! Compare rhsvec
    call compare_real_arrays(rhsvec_orig, rhsvec(1:nsize), "rhsvec", tol, test_ok)
    if (test_ok) then
        num_passed = num_passed + 1
    else
        num_failed = num_failed + 1
    end if

    !=========================================================================
    ! Test Summary
    !=========================================================================
    print *, ""
    print *, "========================================"
    print *, "  Test Summary"
    print *, "========================================"
    print '(A,I3,A)', "  Passed: ", num_passed, " tests"
    print '(A,I3,A)', "  Failed: ", num_failed, " tests"
    print *, "========================================"

    ! Cleanup
    deallocate(y_test, dy_test)
    deallocate(irow_orig, icol_orig, amat_orig, rhsvec_orig)

    if (num_failed > 0) then
        print *, "INTEGRATION TEST FAILED"
        stop 1
    else
        print *, "ALL INTEGRATION TESTS PASSED"
    end if

contains

    subroutine setup_test_state(n_grid)
        !
        ! Set up realistic test state with given grid size
        !
        integer, intent(in) :: n_grid
        integer :: i
        real(dp) :: r_val, dr

        ! Grid parameters
        npoic = n_grid
        npoib = n_grid + 1
        nbaleqs = 4
        neqset = nbaleqs * npoic
        iboutype = 1  ! Boundary condition: npoi = npoic - 1
        npoi_der = 4
        mwind = 2

        ! Radial range (typical tokamak minor radius in cm)
        rmin = 50.0_dp   ! 50 cm
        rmax = 200.0_dp  ! 200 cm
        dr = (rmax - rmin) / real(n_grid, dp)

        ! Allocate grid arrays
        allocate(rb(npoib), rc(npoic), Sb(npoib), Sc(npoic))
        allocate(ipbeg(npoib), ipend(npoib))
        allocate(deriv_coef(2*mwind+1, npoib), reint_coef(2*mwind+1, npoib))

        ! Set up radial grids
        do i = 1, npoib
            rb(i) = rmin + (i-1) * dr
        end do
        do i = 1, npoic
            rc(i) = 0.5_dp * (rb(i) + rb(i+1))
        end do

        ! Surface areas (simplified 2*pi*r for cylindrical approximation)
        do i = 1, npoib
            Sb(i) = 2.0_dp * 3.14159265358979_dp * rb(i)
        end do
        do i = 1, npoic
            Sc(i) = 2.0_dp * 3.14159265358979_dp * rc(i)
        end do

        ! Interpolation stencil indices
        ! The stencil size must be exactly 2*mwind+1 to match deriv_coef
        ! At boundaries, we shift the stencil to stay within bounds
        do i = 1, npoib
            ! Default: centered stencil
            ipbeg(i) = i - mwind
            ipend(i) = i + mwind
            ! Shift stencil if it goes out of bounds
            if (ipbeg(i) < 1) then
                ipend(i) = ipend(i) + (1 - ipbeg(i))
                ipbeg(i) = 1
            end if
            if (ipend(i) > npoic) then
                ipbeg(i) = ipbeg(i) - (ipend(i) - npoic)
                ipend(i) = npoic
            end if
            ! Final safety check
            ipbeg(i) = max(1, ipbeg(i))
            ipend(i) = min(npoic, ipend(i))
        end do

        ! Interpolation coefficients
        ! For real use, these come from polynomial interpolation
        ! Here we use simplified 5-point stencil
        deriv_coef = 0.0_dp
        reint_coef = 0.0_dp
        do i = 1, npoib
            ! Use uniform weights for simplicity (centered difference at mwind+1 position)
            ! Index mapping: coefficient(k) maps to params(ipbeg + k - 1)
            ! For derivative: simple centered difference
            deriv_coef(mwind, i) = -1.0_dp / dr
            deriv_coef(mwind+2, i) = 1.0_dp / dr
            ! For interpolation: simple average of two central points
            reint_coef(mwind, i) = 0.5_dp
            reint_coef(mwind+2, i) = 0.5_dp
        end do

        ! Allocate flux arrays
        allocate(fluxes_dif(nbaleqs, npoib), fluxes_con(nbaleqs, npoib))
        allocate(fluxes_con_nl(nbaleqs, npoib))
        fluxes_dif = 0.0_dp
        fluxes_con = 0.0_dp
        fluxes_con_nl = 0.0_dp

        ! Allocate diffusion coefficient arrays
        allocate(dae11(npoib), dae12(npoib), dae22(npoib))
        allocate(dai11(npoib), dai12(npoib), dai22(npoib))
        allocate(dni22(npoib), visca(npoib))
        allocate(dqle11(npoib), dqle12(npoib), dqle21(npoib), dqle22(npoib))
        allocate(dqli11(npoib), dqli12(npoib), dqli21(npoib), dqli22(npoib))
        allocate(gpp_av(npoib))
        allocate(sqrt_g_times_B_theta_over_c(npoib))
        allocate(Ercov(npoib), Ercov_lin(npoib))
        allocate(polforce(npoib), polforce_ql(npoib), qlheat_e(npoib), qlheat_i(npoib))
        allocate(T_EM_phi_e(npoib), T_EM_phi_i(npoib))
        allocate(T_EM_phi_e_source(npoib), T_EM_phi_i_source(npoib))
        allocate(dery_equisource(neqset))

        ! Set realistic diffusion coefficients (CGS units)
        do i = 1, npoib
            r_val = rb(i)
            ! Anomalous transport coefficients (typical values)
            dae11(i) = 1.0e4_dp    ! cm^2/s
            dae12(i) = 0.5e4_dp
            dae22(i) = 2.0e4_dp
            dai11(i) = 0.8e4_dp
            dai12(i) = 0.4e4_dp
            dai22(i) = 1.5e4_dp
            dni22(i) = 0.5e4_dp
            visca(i) = 1.0e5_dp

            ! Quasilinear transport coefficients
            dqle11(i) = 2.0e3_dp
            dqle12(i) = 1.0e3_dp
            dqle21(i) = 1.5e3_dp
            dqle22(i) = 1.0e3_dp
            dqli11(i) = 1.5e3_dp
            dqli12(i) = 0.8e3_dp
            dqli21(i) = 1.2e3_dp
            dqli22(i) = 0.8e3_dp

            ! Geometric factors
            gpp_av(i) = 1.0_dp + 0.1_dp * (r_val - rmin) / (rmax - rmin)
            sqrt_g_times_B_theta_over_c(i) = 1.0e-3_dp * r_val
        end do

        Ercov = 0.0_dp
        Ercov_lin = 0.0_dp
        polforce = 0.0_dp
        polforce_ql = 0.0_dp
        qlheat_e = 0.0_dp
        qlheat_i = 0.0_dp
        T_EM_phi_e = 0.0_dp
        T_EM_phi_i = 0.0_dp
        T_EM_phi_e_source = 0.0_dp
        T_EM_phi_i_source = 0.0_dp
        dery_equisource = 0.0_dp

        ! Allocate and set up plasma parameter arrays
        allocate(params(nbaleqs, npoic), dot_params(nbaleqs, npoic))
        allocate(ddr_params(nbaleqs, npoib), ddr_params_nl(nbaleqs, npoib))
        allocate(params_lin(nbaleqs, npoic), params_b(nbaleqs, npoib))
        allocate(params_b_lin(nbaleqs, npoib))

        ! Initialize with realistic plasma profiles
        do i = 1, npoic
            r_val = rc(i)
            ! Density profile: peaked, ~1e13 cm^-3 at center
            params(1, i) = 1.0e13_dp * (1.0_dp - 0.5_dp * ((r_val - rmin)/(rmax - rmin))**2)
            ! Rotation frequency: ~1e4 rad/s
            params(2, i) = 1.0e4_dp * (1.0_dp - 0.3_dp * (r_val - rmin)/(rmax - rmin))
            ! Electron temperature: ~1 keV = 1.6e-9 erg at center
            params(3, i) = 1.6e-9_dp * (1.0_dp - 0.7_dp * ((r_val - rmin)/(rmax - rmin))**2)
            ! Ion temperature: similar to Te
            params(4, i) = 1.4e-9_dp * (1.0_dp - 0.7_dp * ((r_val - rmin)/(rmax - rmin))**2)
        end do

        dot_params = 0.0_dp
        ddr_params = 0.0_dp
        ddr_params_nl = 0.0_dp
        params_lin = 0.0_dp
        params_b = 0.0_dp
        params_b_lin = 0.0_dp

        ! Basic plasma parameters
        Z_i = 1.0_dp        ! Hydrogen
        am = 2.0_dp         ! Deuterium mass number
        btor = 2.0e4_dp     ! Toroidal field in Gauss (2 T)
        rtor = 165.0_dp     ! Major radius in cm

        ! Allocate and set wave_code_data arrays
        allocate(q(npoib), Vth(npoib))
        do i = 1, npoib
            r_val = rb(i)
            ! Safety factor profile
            q(i) = 1.0_dp + 2.0_dp * ((r_val - rmin)/(rmax - rmin))**2
            ! Toroidal velocity
            Vth(i) = 5.0e5_dp * (1.0_dp - 0.3_dp * (r_val - rmin)/(rmax - rmin))
        end do

        ! Allocate y and dery arrays in grid_mod
        allocate(y(neqset), dery(neqset))
        y = 0.0_dp
        dery = 0.0_dp

    end subroutine setup_test_state

    subroutine pack_y_from_params(y_vec)
        !
        ! Pack params array into y vector
        !
        real(dp), intent(out) :: y_vec(:)
        integer :: ipoi, ieq, i, n_poi

        if (iboutype == 1) then
            n_poi = npoic - 1
        else
            n_poi = npoic
        end if

        do ipoi = 1, n_poi
            do ieq = 1, nbaleqs
                i = nbaleqs * (ipoi - 1) + ieq
                y_vec(i) = params(ieq, ipoi)
            end do
        end do
    end subroutine pack_y_from_params

    subroutine reset_module_state()
        !
        ! Reset module arrays that may have been modified
        !
        integer :: ipoi, ieq, i, n_poi

        ! Reset intermediate arrays
        ddr_params = 0.0_dp
        ddr_params_nl = 0.0_dp
        params_lin = 0.0_dp
        params_b = 0.0_dp
        params_b_lin = 0.0_dp
        dot_params = 0.0_dp
        fluxes_dif = 0.0_dp
        fluxes_con = 0.0_dp
        fluxes_con_nl = 0.0_dp
        Ercov = 0.0_dp
        Ercov_lin = 0.0_dp
        polforce = 0.0_dp
        qlheat_e = 0.0_dp
        qlheat_i = 0.0_dp
        T_EM_phi_e = 0.0_dp
        T_EM_phi_i = 0.0_dp

        ! Reset matrix module state
        isw_rhs = 0
        nz = 0
        nsize = 0

    end subroutine reset_module_state

    subroutine compare_int_arrays(arr1, arr2, name, passed)
        !
        ! Compare two integer arrays
        !
        integer, intent(in) :: arr1(:), arr2(:)
        character(len=*), intent(in) :: name
        logical, intent(out) :: passed
        integer :: i, n_diff, first_diff

        passed = .true.
        n_diff = 0
        first_diff = 0

        if (size(arr1) /= size(arr2)) then
            print *, "   FAIL: ", name, " size mismatch"
            print *, "     Original size:   ", size(arr1)
            print *, "     Refactored size: ", size(arr2)
            passed = .false.
            return
        end if

        do i = 1, size(arr1)
            if (arr1(i) /= arr2(i)) then
                n_diff = n_diff + 1
                if (first_diff == 0) first_diff = i
            end if
        end do

        if (n_diff > 0) then
            print *, "   FAIL: ", name, " has ", n_diff, " mismatches out of ", size(arr1)
            print *, "     First mismatch at index ", first_diff
            print *, "       Original:   ", arr1(first_diff)
            print *, "       Refactored: ", arr2(first_diff)
            passed = .false.
        else
            print *, "   PASS: ", name, " matches (", size(arr1), " elements)"
        end if
    end subroutine compare_int_arrays

    subroutine compare_real_arrays(arr1, arr2, name, tolerance, passed)
        !
        ! Compare two real arrays with relative tolerance
        !
        real(dp), intent(in) :: arr1(:), arr2(:)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: tolerance
        logical, intent(out) :: passed
        integer :: i, n_diff, first_diff
        real(dp) :: rel_diff, max_rel_diff, ref_val

        passed = .true.
        n_diff = 0
        first_diff = 0
        max_rel_diff = 0.0_dp

        if (size(arr1) /= size(arr2)) then
            print *, "   FAIL: ", name, " size mismatch"
            print *, "     Original size:   ", size(arr1)
            print *, "     Refactored size: ", size(arr2)
            passed = .false.
            return
        end if

        do i = 1, size(arr1)
            ref_val = abs(arr1(i))
            if (ref_val > 1.0e-30_dp) then
                rel_diff = abs((arr1(i) - arr2(i)) / arr1(i))
            else
                rel_diff = abs(arr1(i) - arr2(i))
            end if

            if (rel_diff > max_rel_diff) max_rel_diff = rel_diff

            if (rel_diff > tolerance) then
                n_diff = n_diff + 1
                if (first_diff == 0) first_diff = i
            end if
        end do

        if (n_diff > 0) then
            print *, "   FAIL: ", name, " has ", n_diff, " mismatches (tol=", tolerance, ")"
            print *, "     Max relative diff: ", max_rel_diff
            print *, "     First mismatch at index ", first_diff
            print '(A,E22.15)', "       Original:   ", arr1(first_diff)
            print '(A,E22.15)', "       Refactored: ", arr2(first_diff)
            passed = .false.
        else
            print *, "   PASS: ", name, " matches (", size(arr1), " elements, max_rel_diff=", max_rel_diff, ")"
        end if
    end subroutine compare_real_arrays

end program test_rhs_balance_integration


!===========================================================================
! Mock for calc_equil_diffusion_coeffs
!===========================================================================
subroutine calc_equil_diffusion_coeffs()
    !
    ! Mock implementation - diffusion coefficients are already set up
    ! in setup_test_state, so this is a no-op
    !
end subroutine calc_equil_diffusion_coeffs
