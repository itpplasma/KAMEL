program test_neort_grid
    use iso_fortran_env, only: dp => real64
    use neort_interface, only: read_equil_file
    use control_mod, only: equil_path

    implicit none

    real(dp), parameter :: tolerance = 1.0e-12_dp
    integer :: num_failed

    num_failed = 0
    call test_equil_file_with_seven_columns()

    if (num_failed > 0) then
        print '(A,I0)', "FAILED assertions: ", num_failed
        stop 1
    end if
    print *, "ALL TESTS PASSED"

contains

    subroutine test_equil_file_with_seven_columns()
        ! Only r and psi_tor are needed; both the 7-column writer in
        ! equil_profiles.f90 and the 11-column preprocessor output must work.
        real(dp), dimension(:), allocatable :: r_eff, psi_tor
        integer :: iunit

        equil_path = "test_neort_grid_equil.dat"
        open (newunit=iunit, file=equil_path, status="replace", action="write")
        write (iunit, '(A)') " # nrad = 3"
        write (iunit, '(A)') " # r, q, psi_pol, psi_tor, dphi/dpsi, r_geom, volume"
        write (iunit, '(A)') " 1.0 0.9 -0.1 10.0 -9.0 1.1 5.0"
        write (iunit, '(A)') " 2.0 1.1 -0.4 40.0 -9.0 2.1 20.0"
        write (iunit, '(A)') " 3.0 1.5 -0.9 90.0 -9.0 3.1 45.0"
        close (iunit)

        call read_equil_file(r_eff=r_eff, psi_tor=psi_tor)

        call assert_true(size(r_eff) == 3, "three surfaces read")
        if (size(r_eff) == 3) then
            call assert_equal(r_eff(3), 3.0_dp, "last radius")
            call assert_equal(psi_tor(2), 40.0_dp, "middle toroidal flux")
        end if

        open (newunit=iunit, file=equil_path, status="old")
        close (iunit, status="delete")
    end subroutine test_equil_file_with_seven_columns

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label

        if (.not. condition) then
            print '(A,A)', "FAIL: ", label
            num_failed = num_failed + 1
        end if
    end subroutine assert_true

    subroutine assert_equal(actual, expected, label)
        real(dp), intent(in) :: actual, expected
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tolerance*max(abs(expected), 1.0_dp)) then
            print '(A,A)', "FAIL: ", label
            print '(A,ES22.14)', "  expected: ", expected
            print '(A,ES22.14)', "  actual:   ", actual
            num_failed = num_failed + 1
        end if
    end subroutine assert_equal

end program test_neort_grid
