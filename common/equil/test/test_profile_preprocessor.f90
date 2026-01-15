!> @file test_profile_preprocessor.f90
!> @brief Unit tests for the profile preprocessor module
!>
!> Tests the profile_preprocessor_m module functionality including:
!> - Cubic spline interpolation via the public interface
!> - Profile processing with synthetic data
!> - Coordinate type conversion (sqrt(psi_N) vs psi_N)

program test_profile_preprocessor
    use profile_preprocessor_m

    implicit none

    double precision, parameter :: tolerance = 1.0d-10

    logical :: all_passed

    all_passed = .true.

    print *, "========================================"
    print *, "Testing profile preprocessor"
    print *, "========================================"

    call test_cubic_spline_interpolation(all_passed)
    call test_profile_processing_sqrt_psin(all_passed)
    call test_profile_processing_psin(all_passed)
    call test_getter_functions(all_passed)

    call cleanup_test_files()

    print *, ""
    if (all_passed) then
        print *, "All profile preprocessor tests PASSED"
        stop 0
    else
        print *, "Some profile preprocessor tests FAILED"
        stop 1
    end if

contains

    !> Test cubic spline interpolation via the profile preprocessor
    !> Uses a simple quadratic function y = x^2 which should be reproduced
    !> exactly by a cubic spline (since quadratics are a subset of cubics)
    subroutine test_cubic_spline_interpolation(passed)
        implicit none
        logical, intent(inout) :: passed

        type(profile_preprocessor_t) :: pp
        double precision :: val, expected, error
        double precision :: test_r

        print *, ""
        print *, "Test: Cubic spline interpolation (via getter)"

        ! Create synthetic equilibrium and profile data
        call setup_test_case_quadratic()

        ! Initialize from the test namelist
        call pp%init('test_pp.nml')
        call pp%process()

        ! Test interpolation at various points
        ! Profile is y = x^2 where x = sqrt(psi_N) = r_eff (for this test)
        ! So n(r) = r^2

        ! Test at midpoint r = 0.5
        test_r = 0.5d0
        val = pp%get_n(test_r)
        expected = test_r**2
        error = abs(val - expected)

        print *, "  At r = ", test_r
        print *, "    Computed: ", val, "  Expected: ", expected, "  Error: ", error

        if (error > tolerance) then
            print *, "  FAILED: error exceeds tolerance"
            passed = .false.
        end if

        ! Test at r = 0.25
        test_r = 0.25d0
        val = pp%get_n(test_r)
        expected = test_r**2
        error = abs(val - expected)

        print *, "  At r = ", test_r
        print *, "    Computed: ", val, "  Expected: ", expected, "  Error: ", error

        if (error > tolerance) then
            print *, "  FAILED: error exceeds tolerance"
            passed = .false.
        end if

        ! Test at r = 0.75
        test_r = 0.75d0
        val = pp%get_n(test_r)
        expected = test_r**2
        error = abs(val - expected)

        print *, "  At r = ", test_r
        print *, "    Computed: ", val, "  Expected: ", expected, "  Error: ", error

        if (error > tolerance) then
            print *, "  FAILED: error exceeds tolerance"
            passed = .false.
        else
            print *, "  PASSED"
        end if

        call pp%cleanup()

    end subroutine test_cubic_spline_interpolation

    !> Test profile processing with sqrt(psi_N) input coordinate
    subroutine test_profile_processing_sqrt_psin(passed)
        implicit none
        logical, intent(inout) :: passed

        type(profile_preprocessor_t) :: pp
        double precision :: val, expected, error

        print *, ""
        print *, "Test: Profile processing with sqrt(psi_N) coordinate"

        ! Setup uses COORD_SQRT_PSIN by default
        call setup_test_case_linear()

        call pp%init('test_pp_linear.nml')
        call pp%process()

        ! For this test: input is sqrt(psi_N) from 0 to 1
        ! Te profile: Te(sqrt_psiN) = 1000 - 800*sqrt_psiN (linear)
        ! At r = 0.5, sqrt_psiN ~ 0.5, so Te ~ 1000 - 400 = 600

        val = pp%get_Te(0.5d0)
        expected = 600.d0
        error = abs(val - expected)

        print *, "  Te at r = 0.5: computed = ", val, "  expected = ", expected
        print *, "  Error = ", error

        if (error > 1.d0) then  ! Allow some interpolation error
            print *, "  FAILED: error too large"
            passed = .false.
        else
            print *, "  PASSED"
        end if

        call pp%cleanup()

    end subroutine test_profile_processing_sqrt_psin

    !> Test profile processing with psi_N input coordinate
    subroutine test_profile_processing_psin(passed)
        implicit none
        logical, intent(inout) :: passed

        type(profile_preprocessor_t) :: pp
        double precision :: val, expected, error

        print *, ""
        print *, "Test: Profile processing with psi_N coordinate"

        call setup_test_case_psin()

        call pp%init('test_pp_psin.nml')
        call pp%process()

        ! For this test: input is psi_N from 0 to 1
        ! Ti profile: Ti(psiN) = 500 (constant)
        val = pp%get_Ti(0.5d0)
        expected = 500.d0
        error = abs(val - expected)

        print *, "  Ti at r = 0.5: computed = ", val, "  expected = ", expected
        print *, "  Error = ", error

        if (error > 1.d-6) then
            print *, "  FAILED: error too large for constant profile"
            passed = .false.
        else
            print *, "  PASSED"
        end if

        call pp%cleanup()

    end subroutine test_profile_processing_psin

    !> Test getter functions with boundary conditions
    subroutine test_getter_functions(passed)
        implicit none
        logical, intent(inout) :: passed

        type(profile_preprocessor_t) :: pp
        double precision :: val_0, val_1
        logical :: test_ok

        print *, ""
        print *, "Test: Getter functions at boundaries"

        call setup_test_case_quadratic()

        call pp%init('test_pp.nml')
        call pp%process()

        test_ok = .true.

        ! Test at boundaries
        val_0 = pp%get_n(0.d0)
        val_1 = pp%get_n(1.d0)

        print *, "  n(r=0) = ", val_0, "  (expected ~ 0)"
        print *, "  n(r=1) = ", val_1, "  (expected ~ 1)"

        if (abs(val_0) > 1.d-6) then
            print *, "  FAILED: n(0) not close to 0"
            test_ok = .false.
        end if

        if (abs(val_1 - 1.d0) > 0.01d0) then
            print *, "  FAILED: n(1) not close to 1"
            test_ok = .false.
        end if

        if (test_ok) then
            print *, "  PASSED"
        else
            passed = .false.
        end if

        call pp%cleanup()

    end subroutine test_getter_functions

    !---------------------------------------------------------------------------
    ! Helper subroutines to create test data
    !---------------------------------------------------------------------------

    !> Setup test case with quadratic profile y = x^2
    subroutine setup_test_case_quadratic()
        implicit none
        integer :: iunit, i
        integer, parameter :: n = 21
        double precision :: sqrt_psin, val

        ! Create synthetic equilibrium file
        ! Format: r_eff, q, psi (we only need r_eff and psi)
        ! Use simple mapping: r_eff = sqrt_psiN, psi = psiN, psi_max = 1
        open(newunit=iunit, file='test_equil.dat', status='replace')
        write(iunit, '(A)') '# r_eff    q    psi'
        do i = 1, n
            sqrt_psin = dble(i-1) / dble(n-1)
            write(iunit, '(3ES20.12)') sqrt_psin, 1.5d0, sqrt_psin**2
        end do
        close(iunit)

        ! Create input profile: n(sqrt_psiN) = sqrt_psiN^2
        ! Since r_eff = sqrt_psiN in this test, we get n(r) = r^2
        open(newunit=iunit, file='test_n_of_psiN.dat', status='replace')
        write(iunit, '(A)') '# sqrt(psi_N)  n'
        do i = 1, n
            sqrt_psin = dble(i-1) / dble(n-1)
            val = sqrt_psin**2
            write(iunit, '(2ES20.12)') sqrt_psin, val
        end do
        close(iunit)

        ! Create namelist file
        open(newunit=iunit, file='test_pp.nml', status='replace')
        write(iunit, '(A)') '&profile_preprocessor'
        write(iunit, '(A)') "  equil_file = 'test_equil.dat'"
        write(iunit, '(A)') "  input_dir = '.'"
        write(iunit, '(A)') "  output_dir = '.'"
        write(iunit, '(A)') '  coord_type = 1'
        write(iunit, '(A)') "  n_input_file = 'test_n_of_psiN.dat'"
        write(iunit, '(A)') '/'
        close(iunit)

    end subroutine setup_test_case_quadratic

    !> Setup test case with linear Te profile
    subroutine setup_test_case_linear()
        implicit none
        integer :: iunit, i
        integer, parameter :: n = 21
        double precision :: sqrt_psin, val

        ! Create synthetic equilibrium file
        open(newunit=iunit, file='test_equil_linear.dat', status='replace')
        write(iunit, '(A)') '# r_eff    q    psi'
        do i = 1, n
            sqrt_psin = dble(i-1) / dble(n-1)
            write(iunit, '(3ES20.12)') sqrt_psin, 1.5d0, sqrt_psin**2
        end do
        close(iunit)

        ! Create Te profile: Te(sqrt_psiN) = 1000 - 800*sqrt_psiN
        open(newunit=iunit, file='test_Te_linear.dat', status='replace')
        write(iunit, '(A)') '# sqrt(psi_N)  Te [eV]'
        do i = 1, n
            sqrt_psin = dble(i-1) / dble(n-1)
            val = 1000.d0 - 800.d0 * sqrt_psin
            write(iunit, '(2ES20.12)') sqrt_psin, val
        end do
        close(iunit)

        ! Create namelist file
        open(newunit=iunit, file='test_pp_linear.nml', status='replace')
        write(iunit, '(A)') '&profile_preprocessor'
        write(iunit, '(A)') "  equil_file = 'test_equil_linear.dat'"
        write(iunit, '(A)') "  input_dir = '.'"
        write(iunit, '(A)') "  output_dir = '.'"
        write(iunit, '(A)') '  coord_type = 1'
        write(iunit, '(A)') "  Te_input_file = 'test_Te_linear.dat'"
        write(iunit, '(A)') '/'
        close(iunit)

    end subroutine setup_test_case_linear

    !> Setup test case with psi_N coordinate (not sqrt)
    subroutine setup_test_case_psin()
        implicit none
        integer :: iunit, i
        integer, parameter :: n = 21
        double precision :: sqrt_psin, psin

        ! Create synthetic equilibrium file
        open(newunit=iunit, file='test_equil_psin.dat', status='replace')
        write(iunit, '(A)') '# r_eff    q    psi'
        do i = 1, n
            sqrt_psin = dble(i-1) / dble(n-1)
            write(iunit, '(3ES20.12)') sqrt_psin, 1.5d0, sqrt_psin**2
        end do
        close(iunit)

        ! Create Ti profile as function of psi_N (not sqrt)
        ! Ti(psiN) = 500 (constant)
        open(newunit=iunit, file='test_Ti_psin.dat', status='replace')
        write(iunit, '(A)') '# psi_N  Ti [eV]'
        do i = 1, n
            psin = dble(i-1) / dble(n-1)
            write(iunit, '(2ES20.12)') psin, 500.d0
        end do
        close(iunit)

        ! Create namelist file with coord_type = 2 (COORD_PSIN)
        open(newunit=iunit, file='test_pp_psin.nml', status='replace')
        write(iunit, '(A)') '&profile_preprocessor'
        write(iunit, '(A)') "  equil_file = 'test_equil_psin.dat'"
        write(iunit, '(A)') "  input_dir = '.'"
        write(iunit, '(A)') "  output_dir = '.'"
        write(iunit, '(A)') '  coord_type = 2'
        write(iunit, '(A)') "  Ti_input_file = 'test_Ti_psin.dat'"
        write(iunit, '(A)') '/'
        close(iunit)

    end subroutine setup_test_case_psin

    !> Clean up test files
    subroutine cleanup_test_files()
        implicit none
        integer :: iunit
        logical :: exists

        ! Remove test files
        inquire(file='test_equil.dat', exist=exists)
        if (exists) then
            open(newunit=iunit, file='test_equil.dat', status='old')
            close(iunit, status='delete')
        end if

        inquire(file='test_equil_linear.dat', exist=exists)
        if (exists) then
            open(newunit=iunit, file='test_equil_linear.dat', status='old')
            close(iunit, status='delete')
        end if

        inquire(file='test_equil_psin.dat', exist=exists)
        if (exists) then
            open(newunit=iunit, file='test_equil_psin.dat', status='old')
            close(iunit, status='delete')
        end if

        inquire(file='test_n_of_psiN.dat', exist=exists)
        if (exists) then
            open(newunit=iunit, file='test_n_of_psiN.dat', status='old')
            close(iunit, status='delete')
        end if

        inquire(file='test_Te_linear.dat', exist=exists)
        if (exists) then
            open(newunit=iunit, file='test_Te_linear.dat', status='old')
            close(iunit, status='delete')
        end if

        inquire(file='test_Ti_psin.dat', exist=exists)
        if (exists) then
            open(newunit=iunit, file='test_Ti_psin.dat', status='old')
            close(iunit, status='delete')
        end if

        inquire(file='test_pp.nml', exist=exists)
        if (exists) then
            open(newunit=iunit, file='test_pp.nml', status='old')
            close(iunit, status='delete')
        end if

        inquire(file='test_pp_linear.nml', exist=exists)
        if (exists) then
            open(newunit=iunit, file='test_pp_linear.nml', status='old')
            close(iunit, status='delete')
        end if

        inquire(file='test_pp_psin.nml', exist=exists)
        if (exists) then
            open(newunit=iunit, file='test_pp_psin.nml', status='old')
            close(iunit, status='delete')
        end if

        ! Also remove output files that may have been created
        inquire(file='n.dat', exist=exists)
        if (exists) then
            open(newunit=iunit, file='n.dat', status='old')
            close(iunit, status='delete')
        end if

        inquire(file='Te.dat', exist=exists)
        if (exists) then
            open(newunit=iunit, file='Te.dat', status='old')
            close(iunit, status='delete')
        end if

        inquire(file='Ti.dat', exist=exists)
        if (exists) then
            open(newunit=iunit, file='Ti.dat', status='old')
            close(iunit, status='delete')
        end if

    end subroutine cleanup_test_files

end program test_profile_preprocessor
