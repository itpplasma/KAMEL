!> Dispersion-relation profiles for a flre_zone, formerly the C++ disp_profiles
!> class. Per-instance handle (same pattern as kilca_spline_m/
!> kilca_maxwell_eqs_data_m), since each flre_zone owns its own instance.
!>
!> calculate_dispersion_profiles drives the existing Fortran calc_dispersion
!> subroutine (called natively here, exactly as it is from other Fortran code,
!> so the gfortran hidden-string-length convention for the flagback argument is
!> handled by the compiler automatically -- not via the old C ABI's manual
!> trailing-length argument).
!>
!> sort_dispersion_profiles is NOT translated: SORT_DISPERSION_PROFILES is
!> hard-coded to 0 in CMakeLists.txt, so the C++ #if SORT_DISPERSION_PROFILES
!> == 1 branch calling it was dead code in the active build.
module kilca_disp_profiles_m
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_char, c_double, &
        c_ptr, c_null_ptr, c_loc, c_f_pointer
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: disp_profiles_create, disp_profiles_destroy
    public :: disp_profiles_calculate, disp_profiles_save

    type :: disp_profiles_t
        integer(c_int) :: Nwaves
        character(len=1) :: flag_back
        integer(c_int) :: dimx
        type(c_ptr) :: x_ptr = c_null_ptr
        integer(c_int) :: dimk
        real(c_double), allocatable :: k(:)
        integer(c_int) :: dimp
        real(c_double), allocatable :: p(:)
    end type disp_profiles_t

    interface
        subroutine calc_dispersion(r, flagback, flagprint, kval, polvec)
            import :: dp
            real(dp), intent(in) :: r
            character(*), intent(in) :: flagback
            integer, intent(in) :: flagprint
            real(dp), intent(out) :: kval(*)
            real(dp), intent(out) :: polvec(*)
        end subroutine calc_dispersion

        function save_cmplx_matrix_to_one_file(Nrows, Ncols, Npoints, xgrid, arr, full_name) &
            result(ierr) bind(C, name="save_cmplx_matrix_to_one_file")
            import :: c_int, c_double, c_char
            integer(c_int), value :: Nrows, Ncols, Npoints
            real(c_double), intent(in) :: xgrid(*)
            real(c_double), intent(in) :: arr(*)
            character(kind=c_char), intent(in) :: full_name(*)
            integer(c_int) :: ierr
        end function save_cmplx_matrix_to_one_file
    end interface

contains

    function disp_profiles_create(Nw, dimx_p, x_p, flag_back_p) result(handle) &
        bind(C, name="disp_profiles_create_")
        integer(c_int), value :: Nw, dimx_p
        type(c_ptr), value :: x_p
        type(c_ptr), value :: flag_back_p
        integer(c_intptr_t) :: handle
        type(disp_profiles_t), pointer :: d
        character(kind=c_char), pointer :: fb

        allocate (d)
        d%Nwaves = Nw
        d%dimx = dimx_p
        d%x_ptr = x_p

        call c_f_pointer(flag_back_p, fb)
        d%flag_back = fb

        d%dimk = 2*Nw
        allocate (d%k(d%dimk*dimx_p))

        d%dimp = 2*Nw*Nw
        allocate (d%p(d%dimp*dimx_p))

        handle = transfer(c_loc(d), handle)
    end function disp_profiles_create

    subroutine disp_profiles_destroy(handle) bind(C, name="disp_profiles_destroy_")
        integer(c_intptr_t), value :: handle
        type(disp_profiles_t), pointer :: d

        if (handle == 0_c_intptr_t) return
        call handle_to_d(handle, d)
        if (allocated(d%k)) deallocate (d%k)
        if (allocated(d%p)) deallocate (d%p)
        deallocate (d)
    end subroutine disp_profiles_destroy

    subroutine disp_profiles_calculate(handle) bind(C, name="disp_profiles_calculate_")
        integer(c_intptr_t), value :: handle
        type(disp_profiles_t), pointer :: d
        real(c_double), pointer :: x(:)
        integer(c_int) :: i

        call handle_to_d(handle, d)
        call c_f_pointer(d%x_ptr, x, [d%dimx])

        do i = 0, d%dimx - 1
            call calc_dispersion(x(i + 1), d%flag_back, 0, &
                                 d%k(d%dimk*i + 1:), d%p(d%dimp*i + 1:))
        end do
    end subroutine disp_profiles_calculate

    subroutine disp_profiles_save(handle, filename) bind(C, name="disp_profiles_save_")
        integer(c_intptr_t), value :: handle
        type(c_ptr), value :: filename
        type(disp_profiles_t), pointer :: d
        real(c_double), pointer :: x(:)
        character(kind=c_char), pointer :: fname(:)
        integer(c_int) :: ierr

        call handle_to_d(handle, d)
        call c_f_pointer(d%x_ptr, x, [d%dimx])
        call c_f_pointer(filename, fname, [1024])

        ierr = save_cmplx_matrix_to_one_file(d%Nwaves, 1, d%dimx, x, d%k, fname)
    end subroutine disp_profiles_save

    subroutine handle_to_d(handle, d)
        integer(c_intptr_t), value :: handle
        type(disp_profiles_t), pointer, intent(out) :: d
        type(c_ptr) :: cp
        cp = transfer(handle, cp)
        call c_f_pointer(cp, d)
    end subroutine handle_to_d

end module kilca_disp_profiles_m
