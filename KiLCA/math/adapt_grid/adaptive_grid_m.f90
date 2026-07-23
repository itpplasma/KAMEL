!> Scalar adaptive 1D grid refinement, formerly
!> KiLCA/math/adapt_grid/adaptive_grid.cpp. Only the vector-valued entry
!> calc_adaptive_1D_grid_4vector_ is live: its sole caller is the Fortran
!> sysmat_profiles_create (kilca_sysmat_profiles_m), which reaches it through
!> a bind(C, name="calc_adaptive_1D_grid_4vector_") interface block, so the
!> C symbol name is preserved and that caller needs no change.
!>
!> The companion scalar entry calc_adaptive_1D_grid_ and its only helper
!> check_and_remove_grid_condensations_ had zero callers anywhere in KiLCA or
!> QL-Balance (verified by grep across all C++, header and Fortran sources);
!> they are dropped, matching this port's dead-code precedent. set_interval,
!> add_new_interval_to_the_array, update_max_err_interval and eval_error were
!> shared by both entries and survive as private helpers of the live one.
module adaptive_grid_m
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr, c_funptr, &
                                           c_f_procpointer
    implicit none
    private

    public :: calc_adaptive_1D_grid_4vector

    type :: interval5_t
        real(c_double) :: xs, ys
        real(c_double) :: xl, yl
        real(c_double) :: xm, ym
        real(c_double) :: xr, yr
        real(c_double) :: xe, ye
    end type interval5_t

    abstract interface
        subroutine sample_cb(r, fval, ctx) bind(C)
            import :: c_double, c_ptr
            real(c_double), intent(in) :: r
            real(c_double), intent(out) :: fval
            type(c_ptr), value :: ctx
        end subroutine sample_cb
    end interface

contains

    !> The callback f writes one scalar accuracy proxy per evaluated abscissa
    !> (it stores the full vector through p); x and y carry the initial grid on
    !> input, and on output dimx holds the count of evaluated grid points while
    !> eps holds the error reached. The grid itself is collected through p by
    !> the callback, so this routine never sorts or emits it (unlike the
    !> dropped scalar entry).
    subroutine calc_adaptive_1D_grid_4vector(f, p, max_dimx, eps, dimx, x, y) &
        bind(C, name="calc_adaptive_1D_grid_4vector_")
        type(c_funptr), value :: f
        type(c_ptr), value :: p
        integer(c_int), intent(in) :: max_dimx
        real(c_double), intent(inout) :: eps
        integer(c_int), intent(inout) :: dimx
        real(c_double), intent(in) :: x(0:*)
        real(c_double), intent(in) :: y(0:*)

        procedure(sample_cb), pointer :: cb
        type(interval5_t), allocatable :: Iarr(:)
        real(c_double), allocatable :: err(:)
        integer(c_int) :: max_dimI, dimI, i, max_ind
        real(c_double) :: max_err

        call c_f_procpointer(f, cb)

        max_dimI = (max_dimx - 1)/4

        allocate (Iarr(0:max_dimI - 1))
        allocate (err(0:max_dimI - 1))

        do i = 0, dimx - 2
            call set_interval(Iarr(i), x(i), x(i + 1), y(i), y(i + 1), cb, p)
            call eval_error(Iarr(i), err(i))
        end do

        dimI = dimx - 1
        max_err = 0.0_c_double

        do while (dimI < max_dimI)
            max_err = err(0)
            max_ind = 0

            do i = 1, dimI - 1
                if (err(i) > max_err) then
                    max_err = err(i)
                    max_ind = i
                end if
            end do

            if (max_err < eps) exit

            call add_new_interval_to_the_array(Iarr, max_ind, dimI, cb, p)
            call eval_error(Iarr(dimI), err(dimI))
            dimI = dimI + 1

            call update_max_err_interval(Iarr, max_ind, cb, p)
            call eval_error(Iarr(max_ind), err(max_ind))
        end do

        eps = max_err
        dimx = 4*dimI + 1
    end subroutine calc_adaptive_1D_grid_4vector

    subroutine set_interval(iv, x1, x2, y1, y2, cb, p)
        type(interval5_t), intent(inout) :: iv
        real(c_double), intent(in) :: x1, x2, y1, y2
        procedure(sample_cb) :: cb
        type(c_ptr), value :: p

        iv%xs = x1
        iv%ys = y1

        iv%xe = x2
        iv%ye = y2

        iv%xm = 0.5_c_double*(x1 + x2)
        call cb(iv%xm, iv%ym, p)

        iv%xl = (3.0_c_double*x1 + x2)/4.0_c_double
        call cb(iv%xl, iv%yl, p)

        iv%xr = (x1 + 3.0_c_double*x2)/4.0_c_double
        call cb(iv%xr, iv%yr, p)
    end subroutine set_interval

    subroutine add_new_interval_to_the_array(Iarr, max_ind, Icount, cb, p)
        type(interval5_t), intent(inout) :: Iarr(0:*)
        integer(c_int), intent(in) :: max_ind, Icount
        procedure(sample_cb) :: cb
        type(c_ptr), value :: p

        Iarr(Icount)%xs = Iarr(max_ind)%xm
        Iarr(Icount)%ys = Iarr(max_ind)%ym

        Iarr(Icount)%xe = Iarr(max_ind)%xe
        Iarr(Icount)%ye = Iarr(max_ind)%ye

        Iarr(Icount)%xm = Iarr(max_ind)%xr
        Iarr(Icount)%ym = Iarr(max_ind)%yr

        Iarr(Icount)%xl = (3.0_c_double*Iarr(Icount)%xs + Iarr(Icount)%xe)/4.0_c_double
        call cb(Iarr(Icount)%xl, Iarr(Icount)%yl, p)

        Iarr(Icount)%xr = (Iarr(Icount)%xs + 3.0_c_double*Iarr(Icount)%xe)/4.0_c_double
        call cb(Iarr(Icount)%xr, Iarr(Icount)%yr, p)
    end subroutine add_new_interval_to_the_array

    subroutine update_max_err_interval(Iarr, max_ind, cb, p)
        type(interval5_t), intent(inout) :: Iarr(0:*)
        integer(c_int), intent(in) :: max_ind
        procedure(sample_cb) :: cb
        type(c_ptr), value :: p

        Iarr(max_ind)%xe = Iarr(max_ind)%xm
        Iarr(max_ind)%ye = Iarr(max_ind)%ym

        Iarr(max_ind)%xm = Iarr(max_ind)%xl
        Iarr(max_ind)%ym = Iarr(max_ind)%yl

        Iarr(max_ind)%xl = (3.0_c_double*Iarr(max_ind)%xs + Iarr(max_ind)%xe) &
                           /4.0_c_double
        call cb(Iarr(max_ind)%xl, Iarr(max_ind)%yl, p)

        Iarr(max_ind)%xr = (Iarr(max_ind)%xs + 3.0_c_double*Iarr(max_ind)%xe) &
                           /4.0_c_double
        call cb(Iarr(max_ind)%xr, Iarr(max_ind)%yr, p)
    end subroutine update_max_err_interval

    subroutine eval_error(iv, err)
        type(interval5_t), intent(in) :: iv
        real(c_double), intent(out) :: err

        err = abs((iv%xe - iv%xs)*abs(iv%ys - 4.0_c_double*iv%yl + 6.0_c_double*iv%ym &
                                      - 4.0_c_double*iv%yr + iv%ye))
    end subroutine eval_error

end module adaptive_grid_m
