module test_periodization_helpers_m
    !> Helper module hosting the toy aperiodic-function callback used by the
    !> periodization regression test. Defining it as a module procedure (rather
    !> than a contained/internal one) makes it a safe actual argument for the
    !> procedure(aperfuns_i) dummy across all Fortran compilers.
    use KIM_kinds_m, only: dp
    implicit none
    private
    public :: toy_aperfuns

contains

    !> Aperiodic sample functions: funs(1)=x, funs(2)=x**2. Matches aperfuns_i.
    subroutine toy_aperfuns(nfuns, x, funs)
        integer, intent(in) :: nfuns
        real(dp), intent(in) :: x
        real(dp), intent(out) :: funs(nfuns)
        funs(1) = x
        funs(2) = x**2
    end subroutine toy_aperfuns

end module test_periodization_helpers_m

program test_periodization
    !> Regression test for the vendored periodization utility
    !> (KIM/src/background_equilibrium/periodization.f90). Pins make_periodic
    !> against the known reference-program output (fort.3) for
    !> aperfuns=[x, x**2], x_mid=0.5, dx_asis=0.25, dx_tr=0.5 (period L=1.5),
    !> plus periodicity and as-is-region fidelity checks.
    use KIM_kinds_m, only: dp
    use periodization_m, only: make_periodic
    use test_periodization_helpers_m, only: toy_aperfuns
    implicit none

    logical :: all_passed
    all_passed = .true.

    call test_reference_rows(all_passed)
    call test_periodicity(all_passed)
    call test_as_is_region(all_passed)

    if (all_passed) then
        print *, 'All tests PASSED'
        stop 0
    else
        print *, 'Some periodization tests FAILED'
        stop 1
    end if

contains

    !> Known reference rows from the reference test_driver output (fort.3).
    subroutine test_reference_rows(passed)
        logical, intent(inout) :: passed
        call check_row(-0.99_dp, 0.51_dp, 0.2601_dp, passed)
        call check_row(-0.98_dp, 0.52_dp, 0.2704_dp, passed)
        call check_row(-0.95_dp, 0.55_dp, 0.3025_dp, passed)
        call check_row(0.30_dp, 0.30_dp, 0.09_dp, passed)
        call check_row(0.70_dp, 0.70_dp, 0.49_dp, passed)
        call check_row(1.99_dp, 0.49_dp, 0.2401_dp, passed)
        call check_row(2.00_dp, 0.50_dp, 0.25_dp, passed)
    end subroutine test_reference_rows

    !> make_periodic must be L-periodic in x with L = 2*(dx_asis+dx_tr) = 1.5.
    subroutine test_periodicity(passed)
        logical, intent(inout) :: passed
        real(dp), parameter :: L = 1.5_dp
        real(dp) :: funs0(2), funsL(2)
        real(dp) :: xs(4)
        integer :: i
        xs = [-0.87_dp, 0.13_dp, 0.62_dp, 1.41_dp]
        do i = 1, size(xs)
            call make_periodic(toy_aperfuns, 2, 0.5_dp, 0.25_dp, 0.5_dp, xs(i), funs0)
            call make_periodic(toy_aperfuns, 2, 0.5_dp, 0.25_dp, 0.5_dp, xs(i) + L, funsL)
            call report('periodicity: f(x) == f(x+L)', &
                        maxval(abs(funs0 - funsL)) < 1.0e-12_dp, passed)
        end do
    end subroutine test_periodicity

    !> In the as-is region |x-x_mid| <= dx_asis the periodized functions must
    !> equal the aperiodic ones exactly ([x, x**2]).
    subroutine test_as_is_region(passed)
        logical, intent(inout) :: passed
        real(dp) :: funs(2)
        real(dp) :: xs(3)
        integer :: i
        real(dp) :: x
        xs = [0.30_dp, 0.50_dp, 0.70_dp]
        do i = 1, size(xs)
            x = xs(i)
            call make_periodic(toy_aperfuns, 2, 0.5_dp, 0.25_dp, 0.5_dp, x, funs)
            call report('as-is fidelity: funs == [x, x**2]', &
                        abs(funs(1) - x) < 1.0e-12_dp .and. &
                        abs(funs(2) - x**2) < 1.0e-12_dp, passed)
        end do
    end subroutine test_as_is_region

    !> Evaluate make_periodic at x and compare against the reference values.
    subroutine check_row(x, ref1, ref2, passed)
        real(dp), intent(in) :: x, ref1, ref2
        logical, intent(inout) :: passed
        real(dp) :: funs(2)
        character(len=64) :: label
        call make_periodic(toy_aperfuns, 2, 0.5_dp, 0.25_dp, 0.5_dp, x, funs)
        write (label, '(a,f6.2)') 'reference row x=', x
        call report(trim(label), &
                    abs(funs(1) - ref1) < 1.0e-12_dp .and. &
                    abs(funs(2) - ref2) < 1.0e-12_dp, passed)
    end subroutine check_row

    subroutine report(name, ok, passed)
        character(*), intent(in) :: name
        logical, intent(in) :: ok
        logical, intent(inout) :: passed
        if (ok) then
            print *, 'PASS: ', name
        else
            print *, 'FAIL: ', name
            passed = .false.
        end if
    end subroutine report

end program test_periodization
