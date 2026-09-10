!> Minimal domain parameters for this isolated native-interface test.
module constants
    use iso_c_binding, only: c_double, c_double_complex, c_intptr_t
    implicit none
    integer, parameter :: dp = c_double, dpc = c_double_complex, pp = c_intptr_t
end module constants

module flre_sett
    implicit none
    integer :: nwaves = 2, flre_order = 1
end module flre_sett

module sysmat_test_state
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    integer(c_int) :: failures = 0
    integer(c_int) :: spline_frees = 0, spline_calculations = 0
contains
    integer(c_int) function get_sysmat_test_failures()
        if (spline_frees /= 1) then
            print *, "Expected one spline release, got", spline_frees
            failures = failures + 1
        end if
        if (spline_calculations /= 1) then
            print *, "Expected one spline calculation, got", spline_calculations
            failures = failures + 1
        end if
        get_sysmat_test_failures = failures
    end function
end module sysmat_test_state

subroutine calc_diff_sys_matrix(r, flagback, Rarr)
    use constants, only: dp, dpc
    use flre_sett, only: nwaves
    real(dp), intent(in) :: r
    character(len=*), intent(in) :: flagback
    complex(dpc), intent(out) :: Rarr(nwaves, nwaves)
    integer :: i, j, k

    ! Encode each real and imaginary element into the original packed layout.
    do j = 1, nwaves
        do i = 1, nwaves
            k = 2 * ((j - 1) * nwaves + i) - 1
            Rarr(i, j) = cmplx(1000.0_dp * r + k, 1000.0_dp * r + k + 1, dpc)
        end do
    end do
end subroutine

module adaptive_grid_m
    use iso_c_binding, only: c_ptr, c_int, c_double
    implicit none
    abstract interface
        subroutine sample_cb(r, fval, ctx)
            import :: c_double, c_ptr
            real(c_double), intent(in) :: r
            real(c_double), intent(out) :: fval
            type(c_ptr), value :: ctx
        end subroutine sample_cb
    end interface
contains
    subroutine calc_adaptive_1D_grid_4vector(f, p, max_dimx, eps, dimx, x, y)
        procedure(sample_cb) :: f
        type(c_ptr), value :: p
        integer(c_int), intent(in) :: max_dimx
        real(c_double), intent(inout) :: eps
        integer(c_int), intent(inout) :: dimx
        real(c_double), intent(in) :: x(0:*), y(0:*)
        real(c_double) :: fval

        call f(2.5d0, fval, p)
        call f(1.5d0, fval, p)
        dimx = dimx + 2
    end subroutine
end module adaptive_grid_m

module kilca_spline_m
    use iso_c_binding, only: c_ptr, c_int, c_double, c_intptr_t, c_f_pointer
    use sysmat_test_state, only: failures, spline_frees, spline_calculations
    implicit none
contains

    subroutine spline_alloc(N, styp, dimx, x, Carr, sid)
        use, intrinsic :: iso_c_binding, only: c_int, c_double, c_intptr_t
        integer(c_int), value :: N, styp, dimx
        type(c_ptr), value :: x, Carr
        integer(c_intptr_t), intent(out) :: sid
        sid = 42_c_intptr_t
    end subroutine

!> Receives the fully sorted/rearranged M array (y) and verifies, for the
!> known input r values {1.0, 2.0, 3.0, 2.5, 1.5} sampled in that order, that
!> after sorting (expected x = {1.0, 1.5, 2.0, 2.5, 3.0}) every M row equals
!> 1000*x(i) + (j+1) at flat offset i + j*dimx (1-based i, 0-based j).
    subroutine spline_calc(sid, y, Imin, Imax, W, ierr)
        use, intrinsic :: iso_c_binding, only: c_intptr_t, c_double, c_int, c_ptr
        use sysmat_test_state, only: failures
        integer(c_intptr_t), value :: sid
        type(c_ptr), value :: y
        real(c_double), pointer :: yp(:)
        integer(c_int), value :: Imin, Imax
        type(c_ptr), value :: W
        integer(c_int), intent(out) :: ierr
        real(c_double), parameter :: x_sorted(5) = [1.0d0, 1.5d0, 2.0d0, 2.5d0, &
                                                    3.0d0]
        integer(c_int), parameter :: dimx = 5
        integer(c_int) :: i, j, dimM

        spline_calculations = spline_calculations + 1
        dimM = Imax - Imin + 1
        call c_f_pointer(y, yp, [dimx * dimM])
        do i = 1, dimx
            do j = 0, dimM - 1
                if (abs(yp(i + j * dimx) - (1000.0d0 * x_sorted(i) + (j + 1))) > 1.0d-9) &
                    failures = failures + 1
            end do
        end do
        ierr = 0
    end subroutine

    subroutine spline_free(sid)
        use iso_c_binding, only: c_intptr_t
        use sysmat_test_state, only: failures, spline_frees, spline_calculations
        integer(c_intptr_t), value :: sid
        if (sid /= 42_c_intptr_t) failures = failures + 1
        spline_frees = spline_frees + 1
    end subroutine

end module kilca_spline_m

module kilca_inout_m
    implicit none
contains
    function save_cmplx_matrix_(Nrows, Ncols, Npoints, xgrid, arr, full_name) &
        result(ierr)
        use, intrinsic :: iso_c_binding, only: c_int, c_double
        integer(c_int), value :: Nrows, Ncols, Npoints
        real(c_double), intent(in) :: xgrid(0:Npoints - 1)
        real(c_double), intent(in) :: arr(0:2 * Nrows * Ncols * Npoints - 1)
        character(len=*), intent(in) :: full_name
        integer(c_int) :: ierr
        ierr = 0
    end function

end module kilca_inout_m
