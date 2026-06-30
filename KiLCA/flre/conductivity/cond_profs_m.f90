!> Conductivity (K and C matrix) profiles for a flre_zone, formerly the C++
!> cond_profiles class. Per-instance handle (same pattern as
!> kilca_spline_m/kilca_maxwell_eqs_data_m/kilca_sysmat_profiles_m), since each
!> flre_zone owns its own instance.
!>
!> Only the live (polynomial-grid) path is translated: the "exact" per-point
!> fallback (calc_and_spline_conductivity_for_point_/alloc_conductivity_
!> profiles_/the second branch of ctensor/kmatrices in conductivity.f90) is
!> provably unreachable -- background's flag_back is read once at startup
!> (background_m.f90) and never reassigned, so the `backflag == flag_back`
!> guard in ctensor/kmatrices is always true. Likewise the non-_polynom
!> siblings (sample_cond_func, calc_splines_for_K, set_arrays_for_K,
!> smooth_arrays_for_K) and the *_fine debug dumpers are dead (only reachable
!> via the inactive cond_profs.cpp or commented-out call sites).
!>
!> sd/bp (settings*/background*) are write-only in the C++ constructor (read
!> nowhere else in cond_profiles), so they are not translated at all; wd
!> (wave_data*) is only read for path2linear/omov/r_res at construction time,
!> so the caller (flre_zone.cpp, still C++) passes those directly instead of
!> an opaque wave_data handle.
module kilca_cond_profiles_m
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_double, c_char, &
        c_ptr, c_funptr, c_funloc, c_loc, c_f_pointer, c_null_ptr, c_null_char
    implicit none
    private

    public :: cond_profiles_create, cond_profiles_destroy
    public :: get_cond_nk, get_cond_dimk, get_cond_nc, get_cond_dimc
    public :: get_cond_dimx, get_cond_x_ptr, get_cond_k_ptr, get_cond_iks
    public :: get_cond_flag_back, get_cond_path2linear
    public :: eval_all_k_matrices, eval_all_c_matrices
    public :: eval_c_matrices_f, eval_k_matrices_f
    public :: calc_and_spline_conductivity_for_point, delete_conductivity_profiles_f

    real(c_double), parameter :: pi = 3.141592653589793238462643383279502884197d0

    type :: cond_profiles_t
        character(len=1) :: flag_back
        character(len=1024) :: path2linear
        integer(c_int) :: flreo, gal_corr, dimt
        integer(c_int) :: NK, NC
        complex(c_double) :: omov
        integer(c_int) :: dimx
        real(c_double), allocatable :: x(:)
        integer(c_int) :: dimK
        real(c_double), allocatable :: K(:)
        real(c_double), allocatable :: CK(:)
        real(c_double), allocatable :: RK(:)
        integer(c_intptr_t) :: sidK = 0
        integer(c_int) :: dimC
        real(c_double), allocatable :: C(:)
        real(c_double), allocatable :: CC(:)
        real(c_double), allocatable :: RC(:)
        integer(c_intptr_t) :: sidC = 0
        real(c_double), allocatable :: bico(:)
        real(c_double), allocatable :: xt(:)
        real(c_double), allocatable :: yt(:)
    end type cond_profiles_t

    interface
        function adaptive_grid_polynom_res(f, p, a, b, dimy, deg, xdim, eps, &
            r_res, D, eps_res, dim_err, ind_err, x1, y1) result(stat) &
            bind(C, name="adaptive_grid_polynom_res")
            import :: c_funptr, c_ptr, c_int, c_double
            type(c_funptr), value :: f
            type(c_ptr), value :: p
            real(c_double), value :: a, b, r_res, D, eps_res
            real(c_double), intent(inout) :: eps
            integer(c_int), value :: dimy, deg, dim_err
            integer(c_int), intent(inout) :: xdim
            integer(c_int), intent(in) :: ind_err(*)
            real(c_double), intent(inout) :: x1(*), y1(*)
            integer(c_int) :: stat
        end function adaptive_grid_polynom_res

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

        subroutine spline_eval_d_c(sid, dimz, z, Dmin, Dmax, Imin, Imax, R) &
            bind(C, name="spline_eval_d_")
            import :: c_intptr_t, c_double, c_int
            integer(c_intptr_t), value :: sid
            integer(c_int), value :: dimz, Dmin, Dmax, Imin, Imax
            real(c_double), intent(in) :: z(*)
            real(c_double), intent(out) :: R(*)
        end subroutine spline_eval_d_c

        subroutine spline_free_c(sid) bind(C, name="spline_free_")
            import :: c_intptr_t
            integer(c_intptr_t), value :: sid
        end subroutine spline_free_c

        subroutine binomial_coefficients_c(N, BC) bind(C, name="binomial_coefficients_")
            import :: c_int, c_double
            integer(c_int), value :: N
            real(c_double), intent(out) :: BC(*)
        end subroutine binomial_coefficients_c

        subroutine eval_and_set_background_parameters_spec_independent_c(r, flagback, fb_len) &
            bind(C, name="eval_and_set_background_parameters_spec_independent_")
            import :: c_double, c_char, c_int
            real(c_double), intent(in) :: r
            character(kind=c_char), intent(in) :: flagback(*)
            integer(c_int), value :: fb_len
        end subroutine eval_and_set_background_parameters_spec_independent_c

        subroutine eval_and_set_wave_parameters_c(r, flagback, fb_len) &
            bind(C, name="eval_and_set_wave_parameters_")
            import :: c_double, c_char, c_int
            real(c_double), intent(in) :: r
            character(kind=c_char), intent(in) :: flagback(*)
            integer(c_int), value :: fb_len
        end subroutine eval_and_set_wave_parameters_c

        subroutine eval_a_matrix_c() bind(C, name="eval_a_matrix_")
        end subroutine eval_a_matrix_c

        subroutine calc_dem_djmi_arrays_c(r) bind(C, name="calc_dem_djmi_arrays_")
            import :: c_double
            real(c_double), intent(in) :: r
        end subroutine calc_dem_djmi_arrays_c

        subroutine eval_and_set_background_parameters_spec_dependent_c(r, spec, flagback, fb_len) &
            bind(C, name="eval_and_set_background_parameters_spec_dependent_")
            import :: c_double, c_int, c_char
            real(c_double), intent(in) :: r
            integer(c_int), intent(in) :: spec
            character(kind=c_char), intent(in) :: flagback(*)
            integer(c_int), value :: fb_len
        end subroutine eval_and_set_background_parameters_spec_dependent_c

        subroutine eval_and_set_f0_parameters_nu_and_derivs_c(r, spec, flagback, fb_len) &
            bind(C, name="eval_and_set_f0_parameters_nu_and_derivs_")
            import :: c_double, c_int, c_char
            real(c_double), intent(in) :: r
            integer(c_int), intent(in) :: spec
            character(kind=c_char), intent(in) :: flagback(*)
            integer(c_int), value :: fb_len
        end subroutine eval_and_set_f0_parameters_nu_and_derivs_c

        subroutine eval_electric_drift_velocities_c() bind(C, name="eval_electric_drift_velocities_")
        end subroutine eval_electric_drift_velocities_c

        subroutine eval_fgi_arrays_c() bind(C, name="eval_fgi_arrays_")
        end subroutine eval_fgi_arrays_c

        subroutine calc_w2_array_c(spec) bind(C, name="calc_w2_array_")
            import :: c_int
            integer(c_int), intent(in) :: spec
        end subroutine calc_w2_array_c

        subroutine calc_d_array_c() bind(C, name="calc_d_array_")
        end subroutine calc_d_array_c

        subroutine calc_k_matrices_c(K) bind(C, name="calc_k_matrices_")
            import :: c_double
            real(c_double), intent(out) :: K(*)
        end subroutine calc_k_matrices_c

        subroutine calc_k1_matrices_c(K) bind(C, name="calc_k1_matrices_")
            import :: c_double
            real(c_double), intent(out) :: K(*)
        end subroutine calc_k1_matrices_c

        subroutine calc_and_add_galilelian_correction_c(r, spec, flagback, ct, fb_len) &
            bind(C, name="calc_and_add_galilelian_correction_")
            import :: c_double, c_int, c_char
            real(c_double), intent(in) :: r
            integer(c_int), intent(in) :: spec
            character(kind=c_char), intent(in) :: flagback(*)
            real(c_double), intent(inout) :: ct(*)
            integer(c_int), value :: fb_len
        end subroutine calc_and_add_galilelian_correction_c

        real(c_double) function get_background_huge_factor_c() bind(C, name="get_background_huge_factor_")
            import :: c_double
        end function get_background_huge_factor_c

        real(c_double) function get_background_charge_c(i) bind(C, name="get_background_charge_")
            import :: c_int, c_double
            integer(c_int), value :: i
        end function get_background_charge_c

        function get_background_flag_back_c() result(ch) bind(C, name="get_background_flag_back_")
            import :: c_char
            character(kind=c_char) :: ch
        end function get_background_flag_back_c

        integer(c_int) function get_background_n_c() bind(C, name="get_background_N_")
            import :: c_int
        end function get_background_n_c

        subroutine get_flre_order_c(flre_order) bind(C, name="get_flre_order_")
            import :: c_int
            integer(c_int), intent(out) :: flre_order
        end subroutine get_flre_order_c

        subroutine get_gal_corr_c(gal_corr) bind(C, name="get_gal_corr_")
            import :: c_int
            integer(c_int), intent(out) :: gal_corr
        end subroutine get_gal_corr_c

        integer(c_int) function get_background_obj_dimx_c(bp) bind(C, name="get_background_obj_dimx_")
            import :: c_int, c_ptr
            type(c_ptr), value :: bp
        end function get_background_obj_dimx_c

        real(c_double) function get_background_obj_x0_c(bp) bind(C, name="get_background_obj_x0_")
            import :: c_double, c_ptr
            type(c_ptr), value :: bp
        end function get_background_obj_x0_c

        real(c_double) function get_background_obj_xlast_c(bp) bind(C, name="get_background_obj_xlast_")
            import :: c_double, c_ptr
            type(c_ptr), value :: bp
        end function get_background_obj_xlast_c

        real(c_double) function get_wave_data_obj_omov_re_c(wd) bind(C, name="get_wave_data_obj_omov_re_")
            import :: c_double, c_ptr
            type(c_ptr), value :: wd
        end function get_wave_data_obj_omov_re_c

        real(c_double) function get_wave_data_obj_omov_im_c(wd) bind(C, name="get_wave_data_obj_omov_im_")
            import :: c_double, c_ptr
            type(c_ptr), value :: wd
        end function get_wave_data_obj_omov_im_c

        function adaptive_grid_polynom_err(f, p, a, b, dimy, deg, xdim, eps, &
            dim_err, ind_err, x1, y1) result(stat) &
            bind(C, name="adaptive_grid_polynom_err")
            import :: c_funptr, c_ptr, c_int, c_double
            type(c_funptr), value :: f
            type(c_ptr), value :: p
            real(c_double), value :: a, b
            integer(c_int), value :: dimy, deg, dim_err
            integer(c_int), intent(inout) :: xdim
            real(c_double), intent(inout) :: eps
            integer(c_int), intent(in) :: ind_err(*)
            real(c_double), intent(inout) :: x1(*), y1(*)
            integer(c_int) :: stat
        end function adaptive_grid_polynom_err
    end interface

contains

    integer(c_int) function iKs(cp, spec, ttype, p, q, i, j, part, node) result(idx)
        type(cond_profiles_t), intent(in) :: cp
        integer(c_int), intent(in) :: spec, ttype, p, q, i, j, part, node
        idx = node + cp%dimx*(part + 2*(j + 3*(i + 3*(q + (cp%flreo + 1)*(p + &
              (cp%flreo + 1)*(ttype + cp%dimt*spec))))))
    end function iKs

    integer(c_int) function iKa(cp, spec, ttype, p, q, i, j, part, node) result(idx)
        type(cond_profiles_t), intent(in) :: cp
        integer(c_int), intent(in) :: spec, ttype, p, q, i, j, part, node
        idx = part + 2*(j + 3*(i + 3*(q + (cp%flreo + 1)*(p + (cp%flreo + 1)*( &
              ttype + cp%dimt*(spec + 2*node))))))
    end function iKa

    integer(c_int) function iCs(cp, spec, ttype, s, i, j, part, node) result(idx)
        type(cond_profiles_t), intent(in) :: cp
        integer(c_int), intent(in) :: spec, ttype, s, i, j, part, node
        idx = node + cp%dimx*(part + 2*(j + 3*(i + 3*(s + (2*cp%flreo + 1)*( &
              ttype + cp%dimt*spec)))))
    end function iCs

    integer(c_int) function iCa(s, i, j, part) result(idx)
        integer(c_int), intent(in) :: s, i, j, part
        idx = part + 2*(j + 3*(i + 3*s))
    end function iCa

    !> path2linear_p: caller-computed (via eval_path_to_linear_data, still
    !> C++) null-terminated string. a/b: caller-computed adaptive-grid
    !> bounds (max(r1-1,bp->x[0])/min(r2+1,bp->x[dimx-1])) since bp (the C++
    !> background physics instance) is not Fortran-resident. r_res/omov_re/
    !> omov_im: from the zone's still-C++ wave_data instance.
    function cond_profiles_create(path2linear_p, flreo, gal_corr, N, max_dim_c, &
        r1, r2, D, eps_out, eps_res, a, b, r_res, omov_re, omov_im, flag_debug, flag) &
        result(handle) bind(C, name="cond_profiles_create_")
        type(c_ptr), value :: path2linear_p
        integer(c_int), value :: flreo, gal_corr, N, max_dim_c, flag_debug, flag
        real(c_double), value :: r1, r2, D, eps_out, eps_res, a, b, r_res, omov_re, omov_im
        integer(c_intptr_t) :: handle

        type(cond_profiles_t), pointer :: cp
        integer(c_int) :: l, dim_err
        integer(c_int), allocatable :: ind_err(:)
        integer(c_int) :: spec, ttype, p, q, i, j, part
        real(c_double) :: epso, epsi
        integer(c_int) :: stat_unused

        allocate (cp)
        cp%flag_back = get_background_flag_back_c()
        call read_c_string(path2linear_p, cp%path2linear)

        cp%flreo = flreo
        cp%gal_corr = gal_corr
        cp%omov = cmplx(omov_re, omov_im, c_double)
        cp%NK = N + flreo + 1

        allocate (cp%bico(0:(flreo + 1)*(flreo + 1) - 1))
        call binomial_coefficients_c(flreo, cp%bico)

        cp%dimt = 2
        cp%dimK = 2*cp%dimt*(flreo + 1)*(flreo + 1)*3*3*2

        cp%dimx = max_dim_c
        allocate (cp%xt(0:cp%dimx - 1))
        allocate (cp%yt(0:cp%dimK*cp%dimx - 1))

        if (a > r1 .or. b < r2) then
            write (*, '(a)') 'warning: a or b is inside the zone.'
            stop 1
        end if

        dim_err = 2*3*3*2
        allocate (ind_err(0:dim_err - 1))
        l = 0
        do spec = 0, 1
            do ttype = 0, 0
                do p = 0, 0
                    do q = 0, 0
                        do i = 0, 2
                            do j = 0, 2
                                do part = 0, 1
                                    ind_err(l) = iKa(cp, spec, ttype, p, q, i, j, part, 0)
                                    l = l + 1
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do

        epso = eps_out
        epsi = eps_res
        stat_unused = adaptive_grid_polynom_res(c_funloc(sample_cond_func_polynom), c_loc(cp), &
                                                a, b, cp%dimK, cp%NK, cp%dimx, epso, &
                                                r_res, D, epsi, dim_err, ind_err, cp%xt, cp%yt)
        deallocate (ind_err)

        allocate (cp%x(0:cp%dimx - 1))
        allocate (cp%K(0:cp%dimx*cp%dimK - 1))
        allocate (cp%CK(0:cp%dimx*cp%dimK*(cp%NK + 1) - 1))
        allocate (cp%RK(0:(cp%NK + 1)*cp%dimK - 1))

        if (flag == 0) then
            call calc_splines_for_K_polynom(cp)
        else
            call set_arrays_for_K_polynom(cp)
        end if

        deallocate (cp%xt)
        deallocate (cp%yt)

        if (flag_debug > 1) call save_K_matrices(cp, -1, 0)

        if (flag == 0) then
            cp%NC = N
            cp%dimC = 2*cp%dimt*(2*flreo + 1)*3*3*2
            allocate (cp%C(0:cp%dimx*cp%dimC - 1))
            allocate (cp%CC(0:(cp%NC + 1)*cp%dimx*cp%dimC - 1))
            allocate (cp%RC(0:(cp%NC + 1)*cp%dimC - 1))

            call calc_splines_for_C(cp)

            if (flag_debug > 1) call save_C_matrices(cp, -1, 0)
        end if

        handle = transfer(c_loc(cp), handle)
    end function cond_profiles_create

    subroutine cond_profiles_destroy(handle) bind(C, name="cond_profiles_destroy_")
        integer(c_intptr_t), value :: handle
        type(cond_profiles_t), pointer :: cp

        if (handle == 0_c_intptr_t) return
        call handle_to_cp(handle, cp)
        if (cp%sidK /= 0) call spline_free_c(cp%sidK)
        if (cp%sidC /= 0) call spline_free_c(cp%sidC)
        deallocate (cp)
    end subroutine cond_profiles_destroy

    !> bind(C) replacement for calc_and_spline_conductivity_for_point_,
    !> reachable only via the provably-dead branch in ctensor/kmatrices
    !> (conductivity.f90), but translated in full since the symbol is
    !> statically referenced by that still-live Fortran code. sd_ptr/bp_ptr/
    !> wd_ptr/cp_ptr are passed BY REFERENCE (non-VALUE), matching the
    !> original `settings **sd_ptr`/`cond_profiles **cp_ptr` C ABI that
    !> conductivity.f90's implicit-interface call expects (same convention as
    !> eval_c_matrices_f/eval_k_matrices_f above). bp/wd are still C++
    !> instances (per-zone, not the Fortran-resident globals), so their
    !> needed fields are read through small additive C++ accessors
    !> (get_background_obj_*_/get_wave_data_obj_omov_*_).
    subroutine calc_and_spline_conductivity_for_point(sd_ptr, bp_ptr, wd_ptr, &
        flag_back_p, r, cp_ptr) bind(C, name="calc_and_spline_conductivity_for_point_")
        integer(c_intptr_t), intent(in) :: sd_ptr, bp_ptr, wd_ptr
        character(kind=c_char), intent(in) :: flag_back_p(*)
        real(c_double), intent(in) :: r
        integer(c_intptr_t), intent(out) :: cp_ptr

        type(cond_profiles_t), pointer :: cp
        type(c_ptr) :: bp_cptr, wd_cptr
        integer(c_int) :: l, dim_err
        integer(c_int), allocatable :: ind_err(:)
        integer(c_int) :: spec, ttype, p, q, i, j, part
        real(c_double) :: delta, a, b, eps
        real(c_double) :: bp_x0, bp_xlast, omov_re, omov_im
        integer(c_int) :: stat_unused

        allocate (cp)

        bp_cptr = transfer(bp_ptr, bp_cptr)
        wd_cptr = transfer(wd_ptr, wd_cptr)

        cp%flag_back = flag_back_p(1)
        cp%path2linear = ''

        call get_flre_order_c(cp%flreo)
        call get_gal_corr_c(cp%gal_corr)

        allocate (cp%bico(0:(cp%flreo + 1)*(cp%flreo + 1) - 1))
        call binomial_coefficients_c(cp%flreo, cp%bico)

        cp%NK = get_background_n_c() - (cp%flreo + 1)

        cp%dimt = 2
        cp%dimK = 2*cp%dimt*(cp%flreo + 1)*(cp%flreo + 1)*3*3*2

        cp%dimx = 3*(cp%NK + 1)

        allocate (cp%xt(0:cp%dimx - 1))
        allocate (cp%yt(0:cp%dimx*cp%dimK - 1))

        bp_x0 = get_background_obj_x0_c(bp_cptr)
        bp_xlast = get_background_obj_xlast_c(bp_cptr)
        omov_re = get_wave_data_obj_omov_re_c(wd_cptr)
        omov_im = get_wave_data_obj_omov_im_c(wd_cptr)
        cp%omov = cmplx(omov_re, omov_im, c_double)

        delta = 0.01d0*(bp_xlast - bp_x0)
        a = max(bp_x0, r - delta)
        b = min(bp_xlast, r + delta)

        dim_err = 2*3*3*2
        allocate (ind_err(0:dim_err - 1))
        l = 0
        do spec = 0, 1
            do ttype = 0, 0
                do p = 0, 0
                    do q = 0, 0
                        do i = 0, 2
                            do j = 0, 2
                                do part = 0, 1
                                    ind_err(l) = iKa(cp, spec, ttype, p, q, i, j, part, 0)
                                    l = l + 1
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do

        eps = 0.0d0
        stat_unused = adaptive_grid_polynom_err(c_funloc(sample_cond_func_polynom), c_loc(cp), &
                                                a, b, cp%dimK, cp%NK, cp%dimx, eps, &
                                                dim_err, ind_err, cp%xt, cp%yt)
        deallocate (ind_err)

        allocate (cp%x(0:cp%dimx - 1))
        allocate (cp%K(0:cp%dimx*cp%dimK - 1))
        allocate (cp%CK(0:cp%dimx*cp%dimK*(cp%NK + 1) - 1))
        allocate (cp%RK(0:(cp%NK + 1)*cp%dimK - 1))

        call calc_splines_for_K_polynom(cp)

        deallocate (cp%xt)
        deallocate (cp%yt)

        cp%NC = cp%NK - (cp%flreo + 1)
        cp%dimC = 2*cp%dimt*(2*cp%flreo + 1)*3*3*2
        allocate (cp%C(0:cp%dimx*cp%dimC - 1))
        allocate (cp%CC(0:(cp%NC + 1)*cp%dimx*cp%dimC - 1))
        allocate (cp%RC(0:(cp%NC + 1)*cp%dimC - 1))

        call calc_splines_for_C(cp)

        cp_ptr = transfer(c_loc(cp), cp_ptr)
    end subroutine calc_and_spline_conductivity_for_point

    !> bind(C) replacement for delete_conductivity_profiles_f_ (same
    !> reference-passing convention as above).
    subroutine delete_conductivity_profiles_f(cp_ptr) bind(C, name="delete_conductivity_profiles_f_")
        integer(c_intptr_t), intent(in) :: cp_ptr
        type(cond_profiles_t), pointer :: cp

        call handle_to_cp(cp_ptr, cp)
        if (cp%sidK /= 0) call spline_free_c(cp%sidK)
        if (cp%sidC /= 0) call spline_free_c(cp%sidC)
        deallocate (cp)
    end subroutine delete_conductivity_profiles_f

    !> bind(C) callback matching void(double *r, double *f, void *p): mirrors
    !> sample_cond_func_polynom exactly (target buffer is f, not *cp->yt).
    subroutine sample_cond_func_polynom(r, f, p) bind(C)
        real(c_double), intent(in) :: r
        real(c_double), intent(out) :: f(*)
        type(c_ptr), value :: p
        type(cond_profiles_t), pointer :: cp
        integer(c_int) :: spec, ttype

        call c_f_pointer(p, cp)

        call eval_and_set_background_parameters_spec_independent_c(r, cp%flag_back, 1_c_int)
        call eval_and_set_wave_parameters_c(r, cp%flag_back, 1_c_int)
        call eval_a_matrix_c()
        call calc_dem_djmi_arrays_c(r)

        do spec = 0, 1
            call eval_and_set_background_parameters_spec_dependent_c(r, spec, cp%flag_back, 1_c_int)
            call eval_and_set_f0_parameters_nu_and_derivs_c(r, spec, cp%flag_back, 1_c_int)
            call eval_electric_drift_velocities_c()
            call eval_fgi_arrays_c()
            call calc_w2_array_c(spec)
            call calc_d_array_c()

            ttype = 0
            call calc_k_matrices_c(f(iKa(cp, spec, ttype, 0, 0, 0, 0, 0, 0) + 1))

            ttype = 1
            call calc_k1_matrices_c(f(iKa(cp, spec, ttype, 0, 0, 0, 0, 0, 0) + 1))
        end do
    end subroutine sample_cond_func_polynom

    subroutine set_arrays_for_K_polynom(cp)
        type(cond_profiles_t), intent(inout) :: cp
        integer(c_int) :: node, spec, ttype, p, q, i, j, part

        do node = 0, cp%dimx - 1
            cp%x(node) = cp%xt(node)
            do spec = 0, 1
                do ttype = 0, cp%dimt - 1
                    do p = 0, cp%flreo
                        do q = 0, cp%flreo
                            do i = 0, 2
                                do j = 0, 2
                                    do part = 0, 1
                                        cp%K(iKs(cp, spec, ttype, p, q, i, j, part, node)) = &
                                            cp%yt(iKa(cp, spec, ttype, p, q, i, j, part, node))
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine set_arrays_for_K_polynom

    subroutine calc_splines_for_K_polynom(cp)
        type(cond_profiles_t), intent(inout) :: cp
        integer(c_int) :: ierr

        call set_arrays_for_K_polynom(cp)

        call spline_alloc_c(cp%NK, 1, cp%dimx, cp%x, cp%CK, cp%sidK)
        call spline_calc_c(cp%sidK, cp%K, 0, cp%dimK - 1, c_null_ptr, ierr)
    end subroutine calc_splines_for_K_polynom

    subroutine calc_splines_for_C(cp)
        type(cond_profiles_t), intent(inout) :: cp
        integer(c_int) :: spec, ttype, node, s, i, j, part, ierr

        do spec = 0, 1
            do ttype = 0, cp%dimt - 1
                do node = 0, cp%dimx - 1
                    call calc_C_matrices(cp, spec, ttype, cp%x(node), cp%RC)
                    do s = 0, 2*cp%flreo
                        do i = 0, 2
                            do j = 0, 2
                                do part = 0, 1
                                    cp%C(iCs(cp, spec, ttype, s, i, j, part, node)) = &
                                        cp%RC(iCa(s, i, j, part))
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do

        call spline_alloc_c(cp%NC, 1, cp%dimx, cp%x, cp%CC, cp%sidC)
        call spline_calc_c(cp%sidC, cp%C, 0, cp%dimC - 1, c_null_ptr, ierr)
    end subroutine calc_splines_for_C

    !> Mirrors calc_C_matrices exactly: evaluates K-matrices and derivatives,
    !> combines them via the binomial-coefficient sum into the C matrix, scales
    !> by the conductivity prefactor, and applies the Galilean correction.
    subroutine calc_C_matrices(cp, spec, ttype, r, Cout)
        type(cond_profiles_t), intent(inout) :: cp
        integer(c_int), intent(in) :: spec, ttype
        real(c_double), intent(in) :: r
        real(c_double), intent(out) :: Cout(*)

        integer(c_int) :: flreo, dimc, dimk, p, nmin, nmax, n, m, i, j, ind, k
        real(c_double) :: coeff, scale_fac
        complex(c_double) :: Cm(0:9*(2*cp%flreo + 1) - 1)
        complex(c_double) :: cft

        flreo = cp%flreo
        dimc = 9*(2*flreo + 1)

        call eval_K_matrices(cp, spec, ttype, 0, flreo, r, cp%RK)

        Cm = (0.0d0, 0.0d0)

        dimk = 2*9*(flreo + 1)*(flreo + 1)

        do p = 0, 2*flreo
            nmin = max(0, p - flreo)
            nmax = min(p, flreo)
            do n = nmin, nmax
                do m = 0, flreo - (p - n)
                    coeff = (-1.0d0)**(m + p - n)*cp%bico(m + p - n + (p - n)*(flreo + 1))
                    do i = 0, 2
                        do j = 0, 2
                            ind = 2*(j + 3*(i + 3*(n + (flreo + 1)*(m + p - n))))
                            Cm(j + 3*(i + 3*p)) = Cm(j + 3*(i + 3*p)) + &
                                coeff*cmplx(cp%RK(ind + m*dimk), cp%RK(ind + 1 + m*dimk), c_double)
                        end do
                    end do
                end do
            end do
        end do

        if (cp%flag_back(1:1) /= 'f') then
            scale_fac = get_background_huge_factor_c()
        else
            scale_fac = 1.0d0
        end if

        cft = 2.0d0*pi*(0.0d0, 1.0d0)*get_background_charge_c(spec)*get_background_charge_c(spec)/ &
              cp%omov/(r*scale_fac)

        do k = 0, dimc - 1
            Cout(2*k + 1) = real(cft*Cm(k), c_double)
            Cout(2*k + 2) = aimag(cft*Cm(k))
        end do

        if (cp%gal_corr == 1 .and. cp%flag_back(1:1) == 'f' .and. ttype == 0 .and. flreo == 1) then
            call calc_and_add_galilelian_correction_c(r, spec, cp%flag_back, Cout, 1_c_int)
        end if
    end subroutine calc_C_matrices

    !> Matches the original eval_K_matrices: evaluates the spline (with
    !> derivatives Dmin..Dmax) at r and undoes the huge_factor scaling applied
    !> when flag_back != 'f'.
    subroutine eval_K_matrices(cp, spec, ttype, Dmin, Dmax, r, Kout)
        type(cond_profiles_t), intent(in) :: cp
        integer(c_int), intent(in) :: spec, ttype, Dmin, Dmax
        real(c_double), intent(in) :: r
        real(c_double), intent(out) :: Kout(*)
        integer(c_int) :: dimk, ind_ka, n, ind, jj
        real(c_double) :: scale_fac, rloc(1)

        dimk = 2*9*(cp%flreo + 1)*(cp%flreo + 1)
        ind_ka = iKa(cp, spec, ttype, 0, 0, 0, 0, 0, 0)

        rloc(1) = r
        call spline_eval_d_c(cp%sidK, 1, rloc, Dmin, Dmax, ind_ka, ind_ka + dimk - 1, Kout)

        if (cp%flag_back(1:1) /= 'f') then
            do n = Dmin, Dmax
                scale_fac = get_background_huge_factor_c()**n
                ind = dimk*(n - Dmin)
                do jj = 0, dimk - 1
                    Kout(jj + ind + 1) = Kout(jj + ind + 1)/scale_fac
                end do
            end do
        end if
    end subroutine eval_K_matrices

    !> Matches eval_C_matrices (per spec/type spline evaluation; dimc here
    !> already counts (re,im) pairs, unlike calc_C_matrices' dimc).
    subroutine eval_C_matrices(cp, spec, ttype, Dmin, Dmax, r, Cout)
        type(cond_profiles_t), intent(in) :: cp
        integer(c_int), intent(in) :: spec, ttype, Dmin, Dmax
        real(c_double), intent(in) :: r
        real(c_double), intent(out) :: Cout(*)
        integer(c_int) :: dimc, ind_ca, n, ind, jj
        real(c_double) :: scale_fac, rloc(1)

        dimc = 18*(2*cp%flreo + 1)
        ind_ca = dimc*(ttype + cp%dimt*spec)

        rloc(1) = r
        call spline_eval_d_c(cp%sidC, 1, rloc, Dmin, Dmax, ind_ca, ind_ca + dimc - 1, Cout)

        if (cp%flag_back(1:1) /= 'f') then
            do n = Dmin, Dmax
                scale_fac = get_background_huge_factor_c()**n
                ind = dimc*(n - Dmin)
                do jj = 0, dimc - 1
                    Cout(jj + ind + 1) = Cout(jj + ind + 1)/scale_fac
                end do
            end do
        end if
    end subroutine eval_C_matrices

    !> bind(C) replacement for eval_c_matrices_f_, called from the still-live
    !> Fortran ctensor (conductivity.f90), itself called every RHS evaluation
    !> from diff_sys.f90. Handle and index args are passed BY REFERENCE,
    !> matching the original `cond_profiles **cp_ptr`/`int *spec` C ABI that
    !> conductivity.f90's implicit-interface call expects.
    subroutine eval_c_matrices_f(cp_ptr, spec, ttype, Dmin, Dmax, r, ct) &
        bind(C, name="eval_c_matrices_f_")
        integer(c_intptr_t), intent(in) :: cp_ptr
        integer(c_int), intent(in) :: spec, ttype, Dmin, Dmax
        real(c_double), intent(in) :: r
        real(c_double), intent(out) :: ct(*)
        type(cond_profiles_t), pointer :: cp

        call handle_to_cp(cp_ptr, cp)
        call eval_C_matrices(cp, spec, ttype, Dmin, Dmax, r, ct)
    end subroutine eval_c_matrices_f

    !> bind(C) replacement for eval_k_matrices_f_ (same calling convention).
    subroutine eval_k_matrices_f(cp_ptr, spec, ttype, Dmin, Dmax, r, km) &
        bind(C, name="eval_k_matrices_f_")
        integer(c_intptr_t), intent(in) :: cp_ptr
        integer(c_int), intent(in) :: spec, ttype, Dmin, Dmax
        real(c_double), intent(in) :: r
        real(c_double), intent(out) :: km(*)
        type(cond_profiles_t), pointer :: cp

        call handle_to_cp(cp_ptr, cp)
        call eval_K_matrices(cp, spec, ttype, Dmin, Dmax, r, km)
    end subroutine eval_k_matrices_f

    !> Matches eval_all_K_matrices: full dimK spline evaluation (no per-spec
    !> index offset), same huge_factor unscaling. Handle by VALUE (called
    !> from still-C++ flre_quants.cpp, normal C call convention).
    subroutine eval_all_k_matrices(handle, Dmin, Dmax, r, Kout) bind(C, name="eval_all_k_matrices_")
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: Dmin, Dmax
        real(c_double), value :: r
        real(c_double), intent(out) :: Kout(*)
        type(cond_profiles_t), pointer :: cp
        integer(c_int) :: dimk, n, ind, jj
        real(c_double) :: scale_fac, rloc(1)

        call handle_to_cp(handle, cp)
        dimk = cp%dimK
        rloc(1) = r
        call spline_eval_d_c(cp%sidK, 1, rloc, Dmin, Dmax, 0, dimk - 1, Kout)

        if (cp%flag_back(1:1) /= 'f') then
            do n = Dmin, Dmax
                scale_fac = get_background_huge_factor_c()**n
                ind = dimk*(n - Dmin)
                do jj = 0, dimk - 1
                    Kout(jj + ind + 1) = Kout(jj + ind + 1)/scale_fac
                end do
            end do
        end if
    end subroutine eval_all_k_matrices

    !> Matches eval_all_C_matrices.
    subroutine eval_all_c_matrices(handle, Dmin, Dmax, r, Cout) bind(C, name="eval_all_c_matrices_")
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: Dmin, Dmax
        real(c_double), value :: r
        real(c_double), intent(out) :: Cout(*)
        type(cond_profiles_t), pointer :: cp
        integer(c_int) :: dimc, n, ind, jj
        real(c_double) :: scale_fac, rloc(1)

        call handle_to_cp(handle, cp)
        dimc = cp%dimC
        rloc(1) = r
        call spline_eval_d_c(cp%sidC, 1, rloc, Dmin, Dmax, 0, dimc - 1, Cout)

        if (cp%flag_back(1:1) /= 'f') then
            do n = Dmin, Dmax
                scale_fac = get_background_huge_factor_c()**n
                ind = dimc*(n - Dmin)
                do jj = 0, dimc - 1
                    Cout(jj + ind + 1) = Cout(jj + ind + 1)/scale_fac
                end do
            end do
        end if
    end subroutine eval_all_c_matrices

    !> Debug dump (flag_debug > 1 only); numeric formatting is an
    !> ES24.16/uppercase-E approximation of the original "%.16le", not
    !> byte-identical -- acceptable since these files are not part of any
    !> correctness-critical comparison.
    subroutine save_K_matrices(cp, spec, ttype)
        type(cond_profiles_t), intent(in) :: cp
        integer(c_int), intent(in) :: spec, ttype
        character(len=2) :: reim(0:1)
        character(len=1024) :: fname
        integer(c_int) :: p, q, part, k, i, j, u
        real(c_double) :: k_m

        reim(0) = 're'; reim(1) = 'im'

        do p = 0, cp%flreo
            do q = 0, cp%flreo
                do part = 0, 1
                    write (fname, '(a,a,i0,i0,a,a,a)') trim(cp%path2linear), &
                        'debug-data/kt_', p, q, '_', trim(reim(part)), '.dat'
                    open (newunit=u, file=trim(fname), status='replace', action='write')
                    do k = 0, cp%dimx - 1
                        do i = 0, 2
                            do j = 0, 2
                                if (spec == -1) then
                                    k_m = cp%K(iKs(cp, 0, ttype, p, q, i, j, part, k)) + &
                                          cp%K(iKs(cp, 1, ttype, p, q, i, j, part, k))
                                else
                                    k_m = cp%K(iKs(cp, spec, ttype, p, q, i, j, part, k))
                                end if
                                write (u, '(es24.16,1x)', advance='no') k_m
                            end do
                        end do
                        write (u, '(es24.16)') cp%x(k)
                    end do
                    close (u)
                end do
            end do
        end do
    end subroutine save_K_matrices

    subroutine save_C_matrices(cp, spec, ttype)
        type(cond_profiles_t), intent(in) :: cp
        integer(c_int), intent(in) :: spec, ttype
        character(len=2) :: reim(0:1)
        character(len=1024) :: fname
        integer(c_int) :: p, part, k, i, j, u
        real(c_double) :: c_m

        reim(0) = 're'; reim(1) = 'im'

        do p = 0, 2*cp%flreo
            do part = 0, 1
                write (fname, '(a,a,i0,a,a,a)') trim(cp%path2linear), &
                    'debug-data/ct_', p, '_', trim(reim(part)), '.dat'
                open (newunit=u, file=trim(fname), status='replace', action='write')
                do k = 0, cp%dimx - 1
                    do i = 0, 2
                        do j = 0, 2
                            if (spec == -1) then
                                c_m = cp%C(iCs(cp, 0, ttype, p, i, j, part, k)) + &
                                      cp%C(iCs(cp, 1, ttype, p, i, j, part, k))
                            else
                                c_m = cp%C(iCs(cp, spec, ttype, p, i, j, part, k))
                            end if
                            write (u, '(es24.16,1x)', advance='no') c_m
                        end do
                    end do
                    write (u, '(es24.16)') cp%x(k)
                end do
                close (u)
            end do
        end do
    end subroutine save_C_matrices

    integer(c_int) function get_cond_nk(handle) bind(C, name="get_cond_nk_")
        integer(c_intptr_t), value :: handle
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        get_cond_nk = cp%NK
    end function get_cond_nk

    integer(c_int) function get_cond_dimk(handle) bind(C, name="get_cond_dimk_")
        integer(c_intptr_t), value :: handle
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        get_cond_dimk = cp%dimK
    end function get_cond_dimk

    integer(c_int) function get_cond_nc(handle) bind(C, name="get_cond_nc_")
        integer(c_intptr_t), value :: handle
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        get_cond_nc = cp%NC
    end function get_cond_nc

    integer(c_int) function get_cond_dimc(handle) bind(C, name="get_cond_dimc_")
        integer(c_intptr_t), value :: handle
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        get_cond_dimc = cp%dimC
    end function get_cond_dimc

    integer(c_int) function get_cond_dimx(handle) bind(C, name="get_cond_dimx_")
        integer(c_intptr_t), value :: handle
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        get_cond_dimx = cp%dimx
    end function get_cond_dimx

    function get_cond_x_ptr(handle) result(ptr) bind(C, name="get_cond_x_ptr_")
        integer(c_intptr_t), value :: handle
        type(c_ptr) :: ptr
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        ptr = c_loc(cp%x(0))
    end function get_cond_x_ptr

    function get_cond_k_ptr(handle) result(ptr) bind(C, name="get_cond_k_ptr_")
        integer(c_intptr_t), value :: handle
        type(c_ptr) :: ptr
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        ptr = c_loc(cp%K(0))
    end function get_cond_k_ptr

    integer(c_int) function get_cond_iks(handle, spec, ttype, p, q, i, j, part, node) &
        bind(C, name="get_cond_iks_")
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: spec, ttype, p, q, i, j, part, node
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        get_cond_iks = iKs(cp, spec, ttype, p, q, i, j, part, node)
    end function get_cond_iks

    function get_cond_flag_back(handle) result(ch) bind(C, name="get_cond_flag_back_")
        integer(c_intptr_t), value :: handle
        character(kind=c_char) :: ch
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        ch = cp%flag_back(1:1)
    end function get_cond_flag_back

    subroutine get_cond_path2linear(handle, buf) bind(C, name="get_cond_path2linear_")
        integer(c_intptr_t), value :: handle
        character(kind=c_char), intent(out) :: buf(*)
        type(cond_profiles_t), pointer :: cp
        integer :: i, n

        call handle_to_cp(handle, cp)
        n = len_trim(cp%path2linear)
        do i = 1, n
            buf(i) = cp%path2linear(i:i)
        end do
        buf(n + 1) = c_null_char
    end subroutine get_cond_path2linear

    subroutine handle_to_cp(handle, cp)
        integer(c_intptr_t), value :: handle
        type(cond_profiles_t), pointer, intent(out) :: cp
        type(c_ptr) :: ccp
        ccp = transfer(handle, ccp)
        call c_f_pointer(ccp, cp)
    end subroutine handle_to_cp

    subroutine read_c_string(cp_, out)
        type(c_ptr), value :: cp_
        character(len=*), intent(out) :: out
        character(kind=c_char), pointer :: chars(:)
        integer :: i

        call c_f_pointer(cp_, chars, [len(out)])
        out = ''
        do i = 1, len(out)
            if (chars(i) == c_null_char) exit
            out(i:i) = chars(i)
        end do
    end subroutine read_c_string

end module kilca_cond_profiles_m
