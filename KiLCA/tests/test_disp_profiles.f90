!> Unit test for kilca_disp_profiles_m, focused on the per-grid-point array
!> offset arithmetic (k+dimk*i / p+dimp*i in the original C++, translated to
!> Fortran array slicing). The fake calc_dispersion (test_disp_profiles_fakes)
!> encodes the grid index into kval; the fake save_cmplx_matrix_to_one_file
!> verifies every grid point landed at the expected flat offset in the
!> assembled k array and reports the failure count via get_disp_test_failures.
program test_disp_profiles
    use kilca_disp_profiles_m, only: disp_profiles_calculate
    use kilca_disp_profiles_m, only: disp_profiles_create
    use kilca_disp_profiles_m, only: disp_profiles_destroy
    use kilca_disp_profiles_m, only: disp_profiles_save
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_double, c_loc
    use disp_profiles_test_state, only: get_disp_test_failures
    implicit none

    integer(c_int), parameter :: Nw = 2, dimx = 4
    real(c_double), target :: x(dimx)
    integer(c_intptr_t) :: handle
    integer :: i, failures

    do i = 1, dimx
        x(i) = real(i - 1, c_double)
    end do

    handle = disp_profiles_create(Nw, dimx, c_loc(x), 'f')
    if (handle == 0_c_intptr_t) then
        write (*, '(a)') "FAIL: disp_profiles_create returned null handle"
        stop 1
    end if

    call disp_profiles_calculate(handle)
    call disp_profiles_save(handle, 'out.d')
    call disp_profiles_destroy(handle)

    failures = get_disp_test_failures()
    if (failures == 0) then
        write (*, '(a)') "PASS: kilca_disp_profiles_m offset arithmetic matches expected layout"
    else
        write (*, '(a,i0)') "FAILED: ", failures
        stop 1
    end if
end program test_disp_profiles
