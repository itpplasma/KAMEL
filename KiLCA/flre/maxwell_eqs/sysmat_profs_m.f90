!> ODE system matrix (u' = A*u) evaluated on an adaptive radial grid, formerly
!> the C++ sysmat_profiles class. Per-instance handle (same pattern as
!> kilca_spline_m/kilca_maxwell_eqs_data_m/kilca_disp_profiles_m), stored in
!> the legacy mode_data module as integer(pp) sp_ptr.
!>
!> System-matrix evaluation uses a native complex-array interface; adaptive-grid
!> sampling uses native procedure callbacks. Sorting uses fortnum's 1-based argsort.
module kilca_sysmat_profiles_m
    use kilca_legacy_interfaces_m, &
        only: calc_diff_sys_matrix_c => calc_diff_sys_matrix
    use adaptive_grid_m, only: calc_adaptive_1d_grid_4vector
    use kilca_inout_m, only: save_cmplx_matrix => save_cmplx_matrix_
    use kilca_spline_m, only: spline_alloc_c => spline_alloc
    use kilca_spline_m, only: spline_calc_c => spline_calc
    use kilca_spline_m, only: spline_free_c => spline_free
    use, intrinsic :: iso_c_binding, only: &
        c_int, c_intptr_t, c_double, c_char, c_ptr, c_loc, c_f_pointer, c_null_ptr
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
    public :: sample_sysmat_func
contains

    !> Construct profiles from the owning zone's settings and radial bounds.
    function sysmat_profiles_create(Nwaves, flag_back_p, path2linear_p, NC, &
                                    max_dim, eps_out, flag_debug, r1, r2, rm) result(handle)
        integer(c_int), value :: Nwaves, NC, max_dim, flag_debug
        character(len=*), intent(in) :: flag_back_p, path2linear_p
        real(c_double), value :: eps_out, r1, r2, rm
        integer(c_intptr_t) :: handle

        type(sysmat_profiles_t), pointer :: sp
        real(c_double) :: xa(3), ya(3), eps
        integer(c_int) :: i, dimxa, ierr
        integer(c_int), allocatable :: perm(:)

        allocate (sp)
        sp%Nwaves = Nwaves

        sp%flag_back = flag_back_p
        sp%path2linear = path2linear_p

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
        call calc_adaptive_1d_grid_4vector(sample_sysmat_func, c_loc(sp), &
                                           max_dim, eps, dimxa, xa, ya)

        sp%dimx = dimxa
        allocate (sp%x(sp%dimx))
        allocate (sp%M(sp%dimx*sp%dimM))
        allocate (sp%C((sp%N + 1)*sp%dimx*sp%dimM))

        block
            integer(c_int) :: j
            allocate (perm(sp%dimx))
            call argsort(sp%xt(1:sp%dimx), perm)

            ! fortnum argsort is an unstable heapsort; the sysmat grid shares
            ! zone-boundary nodes (duplicated x), so break exact-key ties by
            ! ascending original index, reproducing the oracle's
            ! sort_index_doubles (KiLCA/core/shared.cpp) so the spline is
            ! well-defined and matches the C++ oracle bit-for-bit.
            block
                integer :: ta, tb, kk, tmp
                ta = 1
                do while (ta <= sp%dimx)
                    tb = ta + 1
                    do while (tb <= sp%dimx)
                        if (sp%xt(perm(tb)) /= sp%xt(perm(ta))) exit
                        tb = tb + 1
                    end do
                    do kk = ta + 1, tb - 1
                        tmp = perm(kk)
                        i = kk - 1
                        do while (i >= ta)
                            if (perm(i) <= tmp) exit
                            perm(i + 1) = perm(i)
                            i = i - 1
                        end do
                        perm(i + 1) = tmp
                    end do
                    ta = tb
                end do
            end block

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

        call spline_alloc_c(sp%N, 1, sp%dimx, c_loc(sp%x), c_loc(sp%C), sp%sidM)

        call spline_calc_c(sp%sidM, c_loc(sp%M), 0, sp%dimM - 1, c_null_ptr, ierr)

        if (flag_debug > 1) call sysmat_profiles_save_m(sp, 10)

        handle = transfer(c_loc(sp), handle)
    end function sysmat_profiles_create

    subroutine sysmat_profiles_destroy(handle)
        integer(c_intptr_t), value :: handle
        type(sysmat_profiles_t), pointer :: sp

        if (handle == 0_c_intptr_t) return
        call handle_to_sp(handle, sp)
        if (sp%sidM /= 0_c_intptr_t) call spline_free_c(sp%sidM)
        if (allocated(sp%x)) deallocate (sp%x)
        if (allocated(sp%M)) deallocate (sp%M)
        if (allocated(sp%C)) deallocate (sp%C)
        if (allocated(sp%R)) deallocate (sp%R)
        deallocate (sp)
    end subroutine sysmat_profiles_destroy

    integer(c_intptr_t) function get_sysmat_sidm(handle)
        integer(c_intptr_t), value :: handle
        type(sysmat_profiles_t), pointer :: sp
        call handle_to_sp(handle, sp)
        get_sysmat_sidm = sp%sidM
    end function get_sysmat_sidm

    integer(c_int) function get_sysmat_dimm(handle)
        integer(c_intptr_t), value :: handle
        type(sysmat_profiles_t), pointer :: sp
        call handle_to_sp(handle, sp)
        get_sysmat_dimm = sp%dimM
    end function get_sysmat_dimm

    integer(c_int) function get_sysmat_dimx(handle)
        integer(c_intptr_t), value :: handle
        type(sysmat_profiles_t), pointer :: sp
        call handle_to_sp(handle, sp)
        get_sysmat_dimx = sp%dimx
    end function get_sysmat_dimx

    function get_sysmat_x_ptr(handle) result(ptr)
        integer(c_intptr_t), value :: handle
        type(c_ptr) :: ptr
        type(sysmat_profiles_t), pointer :: sp
        call handle_to_sp(handle, sp)
        ptr = c_loc(sp%x(1))
    end function get_sysmat_x_ptr

    function get_sysmat_flag_back(handle) result(ch)
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

    !> Native adaptive-grid callback storing the grid point and system matrix,
    !> with the convergence target log(1+sum(yt^2)).
    subroutine sample_sysmat_func(r, f, p)
        real(c_double), intent(in) :: r
        real(c_double), intent(out) :: f
        type(c_ptr), value :: p
        type(sysmat_profiles_t), pointer :: sp
        integer(c_int) :: j

        call c_f_pointer(p, sp)

        sp%xt(sp%ind + 1) = r

        block
            complex(c_double) :: matrix(sp%Nwaves, sp%Nwaves)
            call calc_diff_sys_matrix_c(r, sp%flag_back, matrix)
            sp%yt(sp%dimM * sp%ind + 1:sp%dimM * (sp%ind + 1)) = &
                transfer(matrix, sp%yt(1), sp%dimM)
        end block

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
        character(len=1024) :: filename
        character(len=1024) :: fname_f

        dimt = dimf*(sp%dimx - 1)
        allocate (grid(dimt))
        allocate (vals(sp%dimM*dimt))

        do i = 0, sp%dimx - 2
            do k = 0, dimf - 1
                idx = k + dimf*i
                r = sp%x(i + 1) + k*(sp%x(i + 2) - sp%x(i + 1))/dimf
                grid(idx + 1) = r
                block
                    complex(c_double) :: matrix(sp%Nwaves, sp%Nwaves)
                    call calc_diff_sys_matrix_c(r, sp%flag_back, matrix)
                    vals(sp%dimM * idx + 1:sp%dimM * (idx + 1)) = &
                        transfer(matrix, vals(1), sp%dimM)
                end block
            end do
        end do

        fname_f = trim(sp%path2linear)//'debug-data/amat'
        filename = fname_f
        ierr_unused = save_cmplx_matrix(sp%Nwaves, sp%Nwaves, dimt, grid, vals, filename)

        deallocate (grid)
        deallocate (vals)
    end subroutine sysmat_profiles_save_m

end module kilca_sysmat_profiles_m
