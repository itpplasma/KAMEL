!> Unit test for the Fortran spline port (kilca_spline_m).
!>
!> Any odd-degree spline with natural boundary conditions reproduces linear data
!> exactly, value and first derivative, on the whole interval. This exercises
!> spline_alloc/spline_calc (dgesv boundary + dgbsv band solve)/spline_eval and
!> spline_eval_d, plus the opaque-handle round trip, independently of the heavy
!> end-to-end golden record.
program test_spline
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_intptr_t, c_ptr, &
        c_loc, c_null_ptr
    use kilca_spline_m, only: spline_alloc, spline_calc, spline_eval, &
        spline_eval_d, spline_free
    implicit none

    integer(c_int), parameter :: N = 3, stype = 1, dimx = 11, dimz = 5
    real(c_double), parameter :: a = 2.0d0, b = 3.0d0, tol = 1.0d-10
    real(c_double), target :: x(dimx), y(dimx), Carr((N + 1)*dimx)
    real(c_double), target :: z(dimz), R(dimz*2), Rd(dimz*2), Rnode(dimx)
    integer(c_intptr_t) :: sid
    integer(c_int) :: i, ierr, failures
    real(c_double) :: val, der, want

    failures = 0

    do i = 1, dimx
        x(i) = dble(i - 1)/dble(dimx - 1)
        y(i) = a + b*x(i)
    end do
    z = [0.05d0, 0.27d0, 0.5d0, 0.73d0, 0.94d0]

    call spline_alloc(N, stype, dimx, c_loc(x), c_loc(Carr), sid)
    if (sid == 0_c_intptr_t) then
        write (*, '(a)') "FAIL: spline_alloc returned null handle"
        stop 1
    end if

    call spline_calc(sid, c_loc(y), 0_c_int, 0_c_int, c_null_ptr, ierr)
    if (ierr /= 0) then
        write (*, '(a,i0)') "FAIL: spline_calc ierr=", ierr
        failures = failures + 1
    end if

    ! spline_eval: R laid out [value, deriv] per z point (D2 = 2)
    call spline_eval(sid, dimz, c_loc(z), 0_c_int, 1_c_int, 0_c_int, 0_c_int, c_loc(R))
    do i = 1, dimz
        val = R(2*(i - 1) + 1)
        der = R(2*(i - 1) + 2)
        want = a + b*z(i)
        if (abs(val - want) > tol) then
            write (*, '(a,i0,a,es16.8,a,es16.8)') &
                "FAIL: eval value at z(", i, ") got ", val, " want ", want
            failures = failures + 1
        end if
        if (abs(der - b) > tol) then
            write (*, '(a,i0,a,es16.8)') &
                "FAIL: eval deriv at z(", i, ") got ", der, " want 3.0"
            failures = failures + 1
        end if
    end do

    ! spline_eval_d: with a single function and single derivative slot per pair
    ! the layout matches spline_eval for D1=1; verify it agrees with the analytic
    ! values too.
    call spline_eval_d(sid, dimz, c_loc(z), 0_c_int, 1_c_int, 0_c_int, 0_c_int, c_loc(Rd))
    do i = 1, dimz
        val = Rd(2*(i - 1) + 1)
        der = Rd(2*(i - 1) + 2)
        want = a + b*z(i)
        if (abs(val - want) > tol) then
            write (*, '(a,i0)') "FAIL: eval_d value at z=", i
            failures = failures + 1
        end if
        if (abs(der - b) > tol) then
            write (*, '(a,i0)') "FAIL: eval_d deriv at z=", i
            failures = failures + 1
        end if
    end do

    ! interpolation condition at the nodes must hold exactly for the value
    call spline_eval(sid, dimx, c_loc(x), 0_c_int, 0_c_int, 0_c_int, 0_c_int, c_loc(Rnode))
    do i = 1, dimx
        if (abs(Rnode(i) - y(i)) > tol) then
            write (*, '(a,i0)') "FAIL: node interpolation at i=", i
            failures = failures + 1
        end if
    end do

    call spline_free(sid)

    call check_two_functions(failures)

    if (failures == 0) then
        write (*, '(a)') "PASS: kilca_spline_m linear reproduction and node interpolation"
    else
        write (*, '(a,i0,a)') "FAILED with ", failures, " errors"
        stop 1
    end if

contains

    !> Two functions at once (dimy=2) with nonlinear data, checking the node
    !> interpolation condition for both. This exercises the multi-function
    !> strides (j*len in C, j*dimx in y, j*(N+1) in boundary blocks) that the
    !> single-function case leaves untested.
    subroutine check_two_functions(fails)
        integer(c_int), intent(inout) :: fails
        integer(c_int), parameter :: nf = 2
        real(c_double), target :: xx(dimx), yy(dimx*nf), CC((N + 1)*dimx*nf)
        real(c_double), target :: Rn(dimx*nf)
        integer(c_intptr_t) :: h
        integer(c_int) :: ii, ie, jf

        do ii = 1, dimx
            xx(ii) = dble(ii - 1)/dble(dimx - 1)
            yy(ii) = xx(ii)*xx(ii)                 ! f0 = x^2
            yy(dimx + ii) = cos(3.0d0*xx(ii))      ! f1 = cos(3x)
        end do

        call spline_alloc(N, stype, dimx, c_loc(xx), c_loc(CC), h)
        call spline_calc(h, c_loc(yy), 0_c_int, 1_c_int, c_null_ptr, ie)
        if (ie /= 0) then
            write (*, '(a,i0)') "FAIL: two-func spline_calc ierr=", ie
            fails = fails + 1
        end if

        ! spline_eval_d output ordering: R[(j-Imin) + D1*(n-Dmin) + k*D2],
        ! D1 = Imax-Imin+1 = 2, Dmax=Dmin=0 -> R index = j + 2*k for node k.
        call spline_eval_d(h, dimx, c_loc(xx), 0_c_int, 0_c_int, 0_c_int, 1_c_int, c_loc(Rn))
        do ii = 1, dimx
            do jf = 0, nf - 1
                if (abs(Rn(jf + nf*(ii - 1) + 1) - yy(jf*dimx + ii)) > tol) then
                    write (*, '(a,i0,a,i0)') "FAIL: two-func node interp f=", jf, " i=", ii
                    fails = fails + 1
                end if
            end do
        end do

        call spline_free(h)
    end subroutine check_two_functions

end program test_spline
