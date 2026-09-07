!> Unit test for the Fortran Neville port (kilca_neville_m).
!>
!> A degree-deg Neville interpolant reproduces a degree-deg polynomial exactly,
!> value and derivatives, at any point. This checks both the ready wrapper
!> (eval_neville_polynom_ready_) and the grid/search wrapper
!> (eval_neville_polynom_) against the analytic polynomial and its derivatives.
program test_neville
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use kilca_neville_m, only: eval_neville_polynom, eval_neville_polynom_ready
    implicit none

    real(c_double), parameter :: tol = 1.0d-9
    real(c_double) :: xa(0:4), ya(0:4), R(0:2)
    real(c_double) :: xg(0:8), yg(0:8), Rg(0:0)
    integer(c_int) :: i, ind, failures
    real(c_double) :: xq

    failures = 0

    xa = [0.0d0, 0.5d0, 1.0d0, 1.7d0, 2.5d0]
    do i = 0, 4
        ya(i) = p(xa(i))
    end do
    xq = 0.8d0
    call eval_neville_polynom_ready(xa, ya, 4_c_int, xq, 0_c_int, 2_c_int, R)
    call check("ready value", R(0), p(xq))
    call check("ready d1", R(1), dp(xq))
    call check("ready d2", R(2), ddp(xq))

    do i = 0, 8
        xg(i) = 0.3d0*i
        yg(i) = p(xg(i))
    end do
    xq = 1.31d0
    ind = 4
    call eval_neville_polynom(9_c_int, xg, yg, 4_c_int, xq, 0_c_int, 0_c_int, ind, Rg)
    call check("grid value", Rg(0), p(xq))

    if (failures == 0) then
        write (*, '(a)') "PASS: kilca_neville_m reproduces degree-4 polynomial"
    else
        write (*, '(a,i0)') "FAILED: ", failures
        stop 1
    end if

contains

    real(c_double) function p(x)
        real(c_double), intent(in) :: x
        p = 1.0d0 - 2.0d0*x + 0.5d0*x*x - 0.3d0*x*x*x + 0.1d0*x*x*x*x
    end function p

    real(c_double) function dp(x)
        real(c_double), intent(in) :: x
        dp = -2.0d0 + x - 0.9d0*x*x + 0.4d0*x*x*x
    end function dp

    real(c_double) function ddp(x)
        real(c_double), intent(in) :: x
        ddp = 1.0d0 - 1.8d0*x + 1.2d0*x*x
    end function ddp

    subroutine check(label, got, want)
        character(*), intent(in) :: label
        real(c_double), intent(in) :: got, want
        if (abs(got - want) > tol) then
            write (*, '(a,a,2(1x,es23.16))') "FAIL: ", label, got, want
            failures = failures + 1
        end if
    end subroutine check

end program test_neville
