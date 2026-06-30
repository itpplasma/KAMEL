!> ODE system matrix (u' = A*u) evaluated on an adaptive radial grid, formerly
!> the C++ sysmat_profiles class. Per-instance handle (same pattern as
!> kilca_spline_m/kilca_maxwell_eqs_data_m/kilca_disp_profiles_m), addressed
!> the same way the existing mode_data Fortran module already stores it
!> (integer(pp) sp_ptr) -- that bridge is agnostic to whether the handle
!> points at a C++ object or a Fortran one, so no change was needed there.
!>
!> calc_diff_sys_matrix_ and the adaptive-grid callback registration
!> (calc_adaptive_1D_grid_4vector_, a C++ function taking a C function
!> pointer) are called via their existing C ABI entry points, exactly as the
!> former C++ methods did, so the values and grid-refinement behavior are
!> unchanged. Sorting uses fortnum's argsort directly (1-based, native
!> Fortran) instead of the former sort_index_doubles C++ wrapper.
module kilca_sysmat_profiles_m
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_double, c_char, &
        c_ptr, c_funptr, c_funloc, c_loc, c_f_pointer, c_null_ptr
    use fortnum_multiroot, only: argsort
    implicit none
    private

    public :: sysmat_profiles_create, sysmat_profiles_destroy
    public :: get_sysmat_sidm, get_sysmat_dimm, get_sysmat_dimx, get_sysmat_x_ptr
    public :: get_sysmat_flag_back

    type :: sysmat_profiles_t
        character(len=8) :: flag_back
        character(len=1024) :: path2linear
        integer(c_int) :: N, ind, dimx, Nwaves, dimM
        real(c_double), allocatable :: x(:)
        real(c_double), allocatable :: M(:)
        real(c_double), allocatable :: C(:)
        real(c_double), allocatable :: R(:)
        integer(c_intptr_t) :: sidM = 0
        real(c_double), allocatable :: xt(:)
        real(c_double), allocatable :: yt(:)
    end type sysmat_profiles_t

    interface
        subroutine calc_diff_sys_matrix_c(r, flagback, Dmat, fb_len) &
            bind(C, name="calc_diff_sys_matrix_")
            import :: c_double, c_char, c_int
            real(c_double), intent(in) :: r
            character(kind=c_char), intent(in) :: flagback(*)
            real(c_double), intent(out) :: Dmat(*)
            integer(c_int), value :: fb_len
        end subroutine calc_diff_sys_matrix_c

        subroutine calc_adaptive_1d_grid_4vector(f, p, max_dimx, eps, dimx, x, y) &
            bind(C, name="calc_adaptive_1D_grid_4vector_")
            import :: c_funptr, c_ptr, c_int, c_double
            type(c_funptr), value :: f
            type(c_ptr), value :: p
            integer(c_int), intent(in) :: max_dimx
            real(c_double), intent(inout) :: eps
            integer(c_int), intent(inout) :: dimx
            real(c_double), intent(in) :: x(*)
            real(c_double), intent(in) :: y(*)
        end subroutine calc_adaptive_1d_grid_4vector

        subroutine spline_alloc_c(N, styp, dimx, x, Carr, sid) bind(C, name="spline_alloc_")
            import :: c_int, c_double, c_intptr_t
            integer(c_int), value :: N, styp, dimx
            real(c_double), intent(in) :: x(*)
            real(c_double), intent(inout) :: Carr(*)
            integer(c_intptr_t), intent(out) :: sid
        end subroutine spline_alloc_c

        subroutine spline_calc_c(sid, y, Imin, Imax, W, ierr) bind(C, name="spline_calc_")
            import :: c_intptr_t, c_double, c_int, c_ptr
            integer(c_intptr_t), value :: sid
            real(c_double), intent(in) :: y(*)
            integer(c_int), value :: Imin, Imax
            type(c_ptr), value :: W
            integer(c_int), intent(out) :: ierr
        end subroutine spline_calc_c

        function save_cmplx_matrix(Nrows, Ncols, Npoints, xgrid, arr, full_name) &
            result(ierr) bind(C, name="save_cmplx_matrix")
            import :: c_int, c_double, c_char
            integer(c_int), value :: Nrows, Ncols, Npoints
            real(c_double), intent(in) :: xgrid(*)
            real(c_double), intent(in) :: arr(*)
            character(kind=c_char), intent(in) :: full_name(*)
            integer(c_int) :: ierr
        end function save_cmplx_matrix
    end interface

contains

    !> Mirrors sysmat_profiles::calc_and_spline_sysmatrix_profiles, with the
    !> former zone->cp->{flag_back,path2linear,NC} and zone->{Nwaves,max_dim_c,
    !> eps_out,flag_debug,r1,r2,wd->r_res} reads moved to the caller (flre_zone
    !> is still a C++ class; only it can read those fields), passed in directly.
    function sysmat_profiles_create(Nwaves, flag_back_p, path2linear_p, NC, &
        max_dim, eps_out, flag_debug, r1, r2, rm) result(handle) &
        bind(C, name="sysmat_profiles_create_")
        integer(c_int), value :: Nwaves, NC, max_dim, flag_debug
        type(c_ptr), value :: flag_back_p, path2linear_p
        real(c_double), value :: eps_out, r1, r2, rm
        integer(c_intptr_t) :: handle

        type(sysmat_profiles_t), pointer :: sp
        real(c_double) :: xa(3), ya(3), eps
        integer(c_int) :: i, dimxa, ierr
        integer(c_int), allocatable :: perm(:)

        allocate (sp)
        sp%Nwaves = Nwaves

        call read_c_string(flag_back_p, sp%flag_back)
        call read_c_string(path2linear_p, sp%path2linear)

        sp%N = NC

        sp%dimM = 2*Nwaves*Nwaves
        allocate (sp%R((sp%N + 1)*sp%dimM))

        sp%ind = 0

        allocate (sp%xt(max_dim))
        allocate (sp%yt(sp%dimM*max_dim))

        xa(1) = r1
        xa(3) = r2
        if (rm /= 0.0d0 .and. (rm > r1 .and. rm < r2)) then
            xa(2) = 1.01d0*rm
        else
            xa(2) = 0.5d0*(xa(1) + xa(3))
        end if

        dimxa = 3
        do i = 1, dimxa
            call sample_sysmat_func(xa(i), ya(i), c_loc(sp))
        end do

        eps = eps_out
        call calc_adaptive_1d_grid_4vector(c_funloc(sample_sysmat_func), c_loc(sp), &
                                           max_dim, eps, dimxa, xa, ya)

        sp%dimx = dimxa
        allocate (sp%x(sp%dimx))
        allocate (sp%M(sp%dimx*sp%dimM))
        allocate (sp%C((sp%N + 1)*sp%dimx*sp%dimM))

        block
            integer(c_int) :: j
            allocate (perm(sp%dimx))
            call argsort(sp%xt(1:sp%dimx), perm)

            do i = 1, sp%dimx
                sp%x(i) = sp%xt(perm(i))
                do j = 0, sp%dimM - 1
                    sp%M(i + j*sp%dimx) = sp%yt(j + (perm(i) - 1)*sp%dimM + 1)
                end do
            end do
        end block
        deallocate (perm)
        deallocate (sp%xt)
        deallocate (sp%yt)

        call spline_alloc_c(sp%N, 1, sp%dimx, sp%x, sp%C, sp%sidM)

        call spline_calc_c(sp%sidM, sp%M, 0, sp%dimM - 1, c_null_ptr, ierr)

        if (flag_debug > 1) call sysmat_profiles_save_m(sp, 10)

        handle = transfer(c_loc(sp), handle)
    end function sysmat_profiles_create

    subroutine sysmat_profiles_destroy(handle) bind(C, name="sysmat_profiles_destroy_")
        integer(c_intptr_t), value :: handle
        type(sysmat_profiles_t), pointer :: sp

        if (handle == 0_c_intptr_t) return
        call handle_to_sp(handle, sp)
        if (allocated(sp%x)) deallocate (sp%x)
        if (allocated(sp%M)) deallocate (sp%M)
        if (allocated(sp%C)) deallocate (sp%C)
        if (allocated(sp%R)) deallocate (sp%R)
        deallocate (sp)
    end subroutine sysmat_profiles_destroy

    integer(c_intptr_t) function get_sysmat_sidm(handle) bind(C, name="get_sysmat_sidm_")
        integer(c_intptr_t), value :: handle
        type(sysmat_profiles_t), pointer :: sp
        call handle_to_sp(handle, sp)
        get_sysmat_sidm = sp%sidM
    end function get_sysmat_sidm

    integer(c_int) function get_sysmat_dimm(handle) bind(C, name="get_sysmat_dimm_")
        integer(c_intptr_t), value :: handle
        type(sysmat_profiles_t), pointer :: sp
        call handle_to_sp(handle, sp)
        get_sysmat_dimm = sp%dimM
    end function get_sysmat_dimm

    integer(c_int) function get_sysmat_dimx(handle) bind(C, name="get_sysmat_dimx_")
        integer(c_intptr_t), value :: handle
        type(sysmat_profiles_t), pointer :: sp
        call handle_to_sp(handle, sp)
        get_sysmat_dimx = sp%dimx
    end function get_sysmat_dimx

    function get_sysmat_x_ptr(handle) result(ptr) bind(C, name="get_sysmat_x_ptr_")
        integer(c_intptr_t), value :: handle
        type(c_ptr) :: ptr
        type(sysmat_profiles_t), pointer :: sp
        call handle_to_sp(handle, sp)
        ptr = c_loc(sp%x(1))
    end function get_sysmat_x_ptr

    function get_sysmat_flag_back(handle) result(ch) bind(C, name="get_sysmat_flag_back_")
        integer(c_intptr_t), value :: handle
        character(kind=c_char) :: ch
        type(sysmat_profiles_t), pointer :: sp
        call handle_to_sp(handle, sp)
        ch = sp%flag_back(1:1)
    end function get_sysmat_flag_back

    subroutine handle_to_sp(handle, sp)
        integer(c_intptr_t), value :: handle
        type(sysmat_profiles_t), pointer, intent(out) :: sp
        type(c_ptr) :: cp
        cp = transfer(handle, cp)
        call c_f_pointer(cp, sp)
    end subroutine handle_to_sp

    !> bind(C) callback matching void(double *r, double *f, void *p): mirrors
    !> sample_sysmat_func, storing the grid point and evaluated diff-sys-matrix
    !> values, with the convergence target being log(1+sum(yt^2)).
    subroutine sample_sysmat_func(r, f, p) bind(C)
        real(c_double), intent(in) :: r
        real(c_double), intent(out) :: f
        type(c_ptr), value :: p
        type(sysmat_profiles_t), pointer :: sp
        integer(c_int) :: j

        call c_f_pointer(p, sp)

        sp%xt(sp%ind + 1) = r

        call calc_diff_sys_matrix_c(r, sp%flag_back, &
                                    sp%yt(sp%dimM*sp%ind + 1:sp%dimM*sp%ind + sp%dimM), 1_c_int)

        f = 0.0d0
        do j = 0, sp%dimM - 1
            f = f + log(1.0d0 + sp%yt(j + sp%dimM*sp%ind + 1)**2)
        end do

        sp%ind = sp%ind + 1
    end subroutine sample_sysmat_func

    subroutine sysmat_profiles_save_m(sp, dimf)
        type(sysmat_profiles_t), pointer, intent(in) :: sp
        integer(c_int), intent(in) :: dimf
        real(c_double), allocatable :: grid(:), vals(:)
        integer(c_int) :: i, k, idx, dimt, ierr_unused
        real(c_double) :: r
        character(kind=c_char) :: filename(1024)
        character(len=1024) :: fname_f

        dimt = dimf*(sp%dimx - 1)
        allocate (grid(dimt))
        allocate (vals(sp%dimM*dimt))

        do i = 0, sp%dimx - 2
            do k = 0, dimf - 1
                idx = k + dimf*i
                r = sp%x(i + 1) + k*(sp%x(i + 2) - sp%x(i + 1))/dimf
                grid(idx + 1) = r
                call calc_diff_sys_matrix_c(r, sp%flag_back, &
                                            vals(sp%dimM*idx + 1:sp%dimM*idx + sp%dimM), 1_c_int)
            end do
        end do

        fname_f = trim(sp%path2linear)//'debug-data/amat'
        call to_c_string(fname_f, filename)
        ierr_unused = save_cmplx_matrix(sp%Nwaves, sp%Nwaves, dimt, grid, vals, filename)

        deallocate (grid)
        deallocate (vals)
    end subroutine sysmat_profiles_save_m

    subroutine to_c_string(f, c)
        use, intrinsic :: iso_c_binding, only: c_null_char
        character(*), intent(in) :: f
        character(kind=c_char), intent(out) :: c(*)
        integer :: i, n
        n = len_trim(f)
        do i = 1, n
            c(i) = f(i:i)
        end do
        c(n + 1) = c_null_char
    end subroutine to_c_string

    subroutine read_c_string(cp, out)
        use, intrinsic :: iso_c_binding, only: c_null_char
        type(c_ptr), value :: cp
        character(*), intent(out) :: out
        character(kind=c_char), pointer :: chars(:)
        integer :: i, n

        n = len(out)
        call c_f_pointer(cp, chars, [n])
        out = ''
        do i = 1, n
            if (chars(i) == c_null_char) exit
            out(i:i) = chars(i)
        end do
    end subroutine read_c_string

end module kilca_sysmat_profiles_m
