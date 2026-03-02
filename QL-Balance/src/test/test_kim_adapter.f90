program test_kim_adapter
    !
    ! Unit tests for the KIM wave code adapter module.
    !
    ! These tests verify:
    !   1. Module compilation and linkage (the adapter + KIM_lib resolve)
    !   2. interp_complex_profile helper correctness
    !   3. wave_code selector variable default
    !
    ! Note: Full integration tests (calling kim_initialize) require
    ! KIM config files and profile data, and are covered by Task 10
    ! (end-to-end SingleStep validation).
    !
    use kim_wave_code_adapter_m, only: interp_complex_profile
    use control_mod, only: wave_code

    implicit none

    integer :: num_passed, num_failed

    num_passed = 0
    num_failed = 0

    print *, "========================================"
    print *, "  KIM Adapter Module Tests"
    print *, "========================================"
    print *, ""

    call test_wave_code_default()
    call test_interp_complex_constant()
    call test_interp_complex_linear()

    print *, ""
    print *, "========================================"
    print *, "  Test Summary"
    print *, "========================================"
    print '(A,I3,A)', "  Passed: ", num_passed, " tests"
    print '(A,I3,A)', "  Failed: ", num_failed, " tests"
    print *, "========================================"

    if (num_failed > 0) then
        print *, "TESTS FAILED"
        stop 1
    else
        print *, "ALL TESTS PASSED"
    end if

contains

    subroutine assert_equal_real(actual, expected, tol, label)
        real(8), intent(in) :: actual, expected, tol
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tol) then
            print '(A,A,A,ES15.8,A,ES15.8)', "  FAIL: ", label, &
                " got ", actual, " expected ", expected
            num_failed = num_failed + 1
        else
            print '(A,A)', "  PASS: ", label
            num_passed = num_passed + 1
        end if
    end subroutine

    subroutine assert_equal_complex(actual, expected, tol, label)
        complex(8), intent(in) :: actual, expected
        real(8), intent(in) :: tol
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tol) then
            print '(A,A,A,2ES15.8,A,2ES15.8)', "  FAIL: ", label, &
                " got (", real(actual), aimag(actual), &
                ") expected (", real(expected), aimag(expected), ")"
            num_failed = num_failed + 1
        else
            print '(A,A)', "  PASS: ", label
            num_passed = num_passed + 1
        end if
    end subroutine

    ! ------------------------------------------------------------------
    ! Test: wave_code defaults to 'KiLCA'
    ! ------------------------------------------------------------------
    subroutine test_wave_code_default()
        print *, "--- test_wave_code_default ---"

        if (trim(wave_code) == 'KiLCA') then
            print '(A)', "  PASS: wave_code default is 'KiLCA'"
            num_passed = num_passed + 1
        else
            print '(A,A,A)', "  FAIL: wave_code default is '", &
                trim(wave_code), "', expected 'KiLCA'"
            num_failed = num_failed + 1
        end if
    end subroutine

    ! ------------------------------------------------------------------
    ! Test: interp_complex_profile with constant function
    ! A constant complex function should interpolate exactly.
    ! ------------------------------------------------------------------
    subroutine test_interp_complex_constant()
        integer, parameter :: n_old = 5, n_new = 3
        real(8) :: r_old(n_old), r_new(n_new)
        complex(8) :: z_old(n_old), z_new(n_new)
        complex(8), parameter :: val = (2.5d0, -1.3d0)
        real(8), parameter :: tol = 1.0d-10
        integer :: i

        print *, "--- test_interp_complex_constant ---"

        ! Old grid: [10, 20, 30, 40, 50]
        do i = 1, n_old
            r_old(i) = 10.0d0 * dble(i)
            z_old(i) = val
        end do

        ! New grid: [15, 25, 35] (within old grid bounds)
        r_new(1) = 15.0d0
        r_new(2) = 25.0d0
        r_new(3) = 35.0d0

        call interp_complex_profile(n_old, r_old, z_old, n_new, r_new, z_new)

        do i = 1, n_new
            call assert_equal_complex(z_new(i), val, tol, &
                "constant interp at new grid point")
        end do
    end subroutine

    ! ------------------------------------------------------------------
    ! Test: interp_complex_profile with linear function
    ! z(r) = (1+2i)*r should interpolate exactly with polynomial interp.
    ! ------------------------------------------------------------------
    subroutine test_interp_complex_linear()
        integer, parameter :: n_old = 10, n_new = 5
        real(8) :: r_old(n_old), r_new(n_new)
        complex(8) :: z_old(n_old), z_new(n_new)
        complex(8) :: expected
        real(8), parameter :: tol = 1.0d-8
        integer :: i

        print *, "--- test_interp_complex_linear ---"

        ! Old grid: uniform [10, 20, 30, ..., 100]
        do i = 1, n_old
            r_old(i) = 10.0d0 * dble(i)
            z_old(i) = cmplx(1.0d0, 2.0d0, kind=8) * r_old(i)
        end do

        ! New grid: non-coincident points
        r_new(1) = 15.0d0
        r_new(2) = 25.0d0
        r_new(3) = 55.0d0
        r_new(4) = 75.0d0
        r_new(5) = 95.0d0

        call interp_complex_profile(n_old, r_old, z_old, n_new, r_new, z_new)

        do i = 1, n_new
            expected = cmplx(1.0d0, 2.0d0, kind=8) * r_new(i)
            call assert_equal_complex(z_new(i), expected, tol, &
                "linear interp at new grid point")
        end do
    end subroutine

end program test_kim_adapter
