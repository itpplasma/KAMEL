!> Unit test for kilca_disp_profiles_m, focused on the per-grid-point array
!> offset arithmetic (k+dimk*i / p+dimp*i in the original C++, translated to
!> Fortran array slicing). The fake calc_dispersion (test_disp_profiles_fakes)
!> encodes the grid index into kval; the fake save_cmplx_matrix_to_one_file
!> verifies every grid point landed at the expected flat offset in the
!> assembled k array and reports the failure count via get_disp_test_failures_.
program test_disp_profiles
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_double, c_char, &
        c_null_char
    implicit none

    interface
        function disp_profiles_create(Nw, dimx_p, x_p, flag_back_p) result(handle) &
            bind(C, name="disp_profiles_create_")
            import :: c_int, c_intptr_t, c_double, c_char
            integer(c_int), value :: Nw, dimx_p
            real(c_double), intent(in) :: x_p(*)
            character(kind=c_char), intent(in) :: flag_back_p(*)
            integer(c_intptr_t) :: handle
        end function
        subroutine disp_profiles_destroy(handle) bind(C, name="disp_profiles_destroy_")
            import :: c_intptr_t
            integer(c_intptr_t), value :: handle
        end subroutine
        subroutine disp_profiles_calculate(handle) bind(C, name="disp_profiles_calculate_")
            import :: c_intptr_t
            integer(c_intptr_t), value :: handle
        end subroutine
        subroutine disp_profiles_save(handle, filename) bind(C, name="disp_profiles_save_")
            import :: c_intptr_t, c_char
            integer(c_intptr_t), value :: handle
            character(kind=c_char), intent(in) :: filename(*)
        end subroutine
        integer(c_int) function get_disp_test_failures() bind(C, name="get_disp_test_failures_")
            import :: c_int
        end function
    end interface

    integer(c_int), parameter :: Nw = 2, dimx = 4
    real(c_double) :: x(dimx)
    character(kind=c_char) :: flag_back(2)
    character(kind=c_char) :: fname(16)
    integer(c_intptr_t) :: handle
    integer :: i, failures

    do i = 1, dimx
        x(i) = real(i - 1, c_double)
    end do
    flag_back = ['f', c_null_char]
    fname = c_null_char
    fname(1:5) = ['o', 'u', 't', '.', 'd']

    handle = disp_profiles_create(Nw, dimx, x, flag_back)
    if (handle == 0_c_intptr_t) then
        write (*, '(a)') "FAIL: disp_profiles_create returned null handle"
        stop 1
    end if

    call disp_profiles_calculate(handle)
    call disp_profiles_save(handle, fname)
    call disp_profiles_destroy(handle)

    failures = get_disp_test_failures()
    if (failures == 0) then
        write (*, '(a)') "PASS: kilca_disp_profiles_m offset arithmetic matches expected layout"
    else
        write (*, '(a,i0)') "FAILED: ", failures
        stop 1
    end if
end program test_disp_profiles
