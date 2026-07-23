!> Fakes for test_sysmat_profiles, isolating the M-array sort/rearrangement
!> logic in kilca_sysmat_profiles_m (translated from the C++
!> sort_index_doubles + the i/j flat-index loop) from the real adaptive-grid
!> refinement and spline machinery.
!>
!> calc_adaptive_1D_grid_4vector_ calls back through the real callback with
!> two extra points deliberately out of sorted order (mimicking what the
!> production adaptive-grid algorithm would do), so the create() pipeline's
!> argsort + M rearrangement has real work to verify, not just the identity
!> permutation of the 3 already-sorted boundary points. calc_diff_sys_matrix_
!> encodes r into its output; the fake spline_calc_ (which receives the fully
!> rearranged M array as its `y` argument) checks every grid point landed at
!> the expected flat offset and exposes the failure count via
!> get_sysmat_test_failures_.
module sysmat_test_state
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    integer(c_int) :: failures = 0
end module

subroutine calc_diff_sys_matrix_c_fake(r, flagback, Rarr, fb_len) &
    bind(C, name="calc_diff_sys_matrix_")
    use, intrinsic :: iso_c_binding, only: c_double, c_char, c_int
    real(c_double), intent(in) :: r
    character(kind=c_char), intent(in) :: flagback(*)
    real(c_double), intent(out) :: Rarr(*)
    integer(c_int), value :: fb_len
    integer(c_int) :: j

    ! dimM is 8 in the test (Nwaves=2 -> 2*Nwaves*Nwaves=8). Encode r into
    ! every element so each grid point's full row is identifiable.
    do j = 1, 8
        Rarr(j) = 1000.0d0*r + j
    end do
end subroutine

subroutine calc_adaptive_1D_grid_4vector_(f, p, max_dimx, eps, dimx, x, y) &
    bind(C, name="calc_adaptive_1D_grid_4vector_")
    use, intrinsic :: iso_c_binding, only: c_funptr, c_ptr, c_int, c_double, &
        c_f_procpointer
    type(c_funptr), value :: f
    type(c_ptr), value :: p
    integer(c_int), intent(in) :: max_dimx
    real(c_double), intent(inout) :: eps
    integer(c_int), intent(inout) :: dimx
    real(c_double), intent(in) :: x(*)
    real(c_double), intent(in) :: y(*)

    abstract interface
        subroutine sample_cb(r, fval, ctx) bind(C)
            import :: c_double, c_ptr
            real(c_double), intent(in) :: r
            real(c_double), intent(out) :: fval
            type(c_ptr), value :: ctx
        end subroutine sample_cb
    end interface
    procedure(sample_cb), pointer :: cb
    real(c_double) :: fval

    call c_f_procpointer(f, cb)
    call cb(2.5d0, fval, p)
    call cb(1.5d0, fval, p)
    dimx = dimx + 2
end subroutine

subroutine spline_alloc_(N, styp, dimx, x, Carr, sid) bind(C, name="spline_alloc_")
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_intptr_t
    integer(c_int), value :: N, styp, dimx
    real(c_double), intent(in) :: x(*)
    real(c_double), intent(inout) :: Carr(*)
    integer(c_intptr_t), intent(out) :: sid
    sid = 42_c_intptr_t
end subroutine

!> Receives the fully sorted/rearranged M array (y) and verifies, for the
!> known input r values {1.0, 2.0, 3.0, 2.5, 1.5} sampled in that order, that
!> after sorting (expected x = {1.0, 1.5, 2.0, 2.5, 3.0}) every M row equals
!> 1000*x(i) + (j+1) at flat offset i + j*dimx (1-based i, 0-based j).
subroutine spline_calc_(sid, y, Imin, Imax, W, ierr) bind(C, name="spline_calc_")
    use, intrinsic :: iso_c_binding, only: c_intptr_t, c_double, c_int, c_ptr
    use sysmat_test_state, only: failures
    integer(c_intptr_t), value :: sid
    real(c_double), intent(in) :: y(*)
    integer(c_int), value :: Imin, Imax
    type(c_ptr), value :: W
    integer(c_int), intent(out) :: ierr
    real(c_double), parameter :: x_sorted(5) = [1.0d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0]
    integer(c_int), parameter :: dimx = 5
    integer(c_int) :: i, j, dimM

    dimM = Imax - Imin + 1
    do i = 1, dimx
        do j = 0, dimM - 1
            if (abs(y(i + j*dimx) - (1000.0d0*x_sorted(i) + (j + 1))) > 1.0d-9) &
                failures = failures + 1
        end do
    end do
    ierr = 0
end subroutine

function save_cmplx_matrix(Nrows, Ncols, Npoints, xgrid, arr, full_name) &
    result(ierr) bind(C, name="save_cmplx_matrix")
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char
    integer(c_int), value :: Nrows, Ncols, Npoints
    real(c_double), intent(in) :: xgrid(*)
    real(c_double), intent(in) :: arr(*)
    character(kind=c_char), intent(in) :: full_name(*)
    integer(c_int) :: ierr
    ierr = 0
end function

integer(c_int) function get_sysmat_test_failures() bind(C, name="get_sysmat_test_failures_")
    use, intrinsic :: iso_c_binding, only: c_int
    use sysmat_test_state, only: failures
    get_sysmat_test_failures = failures
end function
