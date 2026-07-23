!> Fake calc_dispersion and save_cmplx_matrix_to_one_file for
!> test_disp_profiles, so the per-grid-point array offset arithmetic in
!> kilca_disp_profiles_m (k+dimk*i / p+dimp*i in the original C++, translated
!> to Fortran array slicing) can be checked with known values, independent of
!> real zone/background module setup. calc_dispersion encodes the grid index
!> into its output; save_cmplx_matrix_to_one_file (which receives the full k
!> array disp_profiles_save_ forwards) verifies every grid point landed at the
!> expected flat offset, and exposes the failure count via get_disp_test_failures_.
module disp_profiles_test_state
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    integer(c_int) :: failures = 0
end module

subroutine calc_dispersion(r, flagback, flagprint, kval, polvec)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    real(dp), intent(in) :: r
    character(*), intent(in) :: flagback
    integer, intent(in) :: flagprint
    real(dp), intent(out) :: kval(*)
    real(dp), intent(out) :: polvec(*)
    integer :: idx

    ! r holds the 0-based grid index (the test sets x(i) = i-1).
    idx = nint(r)
    kval(1) = 1000.0_dp + idx
    kval(2) = 2000.0_dp + idx
    polvec(1) = 3000.0_dp + idx
end subroutine

function save_cmplx_matrix_to_one_file(Nrows, Ncols, Npoints, xgrid, arr, full_name) &
    result(ierr) bind(C, name="save_cmplx_matrix_to_one_file")
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char
    use disp_profiles_test_state, only: failures
    integer(c_int), value :: Nrows, Ncols, Npoints
    real(c_double), intent(in) :: xgrid(*)
    real(c_double), intent(in) :: arr(*)
    character(kind=c_char), intent(in) :: full_name(*)
    integer(c_int) :: ierr
    integer(c_int) :: i, dimk

    dimk = 2*Nrows
    do i = 0, Npoints - 1
        if (abs(arr(dimk*i + 1) - (1000.0d0 + i)) > 1.0d-12) failures = failures + 1
        if (abs(arr(dimk*i + 2) - (2000.0d0 + i)) > 1.0d-12) failures = failures + 1
    end do
    ierr = 0
end function

integer(c_int) function get_disp_test_failures() bind(C, name="get_disp_test_failures_")
    use, intrinsic :: iso_c_binding, only: c_int
    use disp_profiles_test_state, only: failures
    get_disp_test_failures = failures
end function
