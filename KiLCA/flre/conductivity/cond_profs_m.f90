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
!> Construction receives the output path, frequency, and resonance position
!> directly from the owning FLRE zone; no settings or background object is owned.
module kilca_cond_profiles_m
    use kilca_legacy_interfaces_m, &
        only: binomial_coefficients_c => binomial_coefficients
    use kilca_legacy_interfaces_m, only: &
        eval_and_set_background_parameters_spec_independent_c => &
            eval_and_set_background_parameters_spec_independent
    use kilca_legacy_interfaces_m, only: &
        eval_and_set_wave_parameters_c => &
            eval_and_set_wave_parameters
    use kilca_legacy_interfaces_m, only: eval_a_matrix_c => eval_a_matrix
    use kilca_legacy_interfaces_m, &
        only: calc_dem_djmi_arrays_c => calc_dem_djmi_arrays
    use kilca_legacy_interfaces_m, only: &
        eval_and_set_background_parameters_spec_dependent_c => &
            eval_and_set_background_parameters_spec_dependent
    use kilca_legacy_interfaces_m, only: &
        eval_and_set_f0_parameters_nu_and_derivs_c => &
            eval_and_set_f0_parameters_nu_and_derivs
    use kilca_legacy_interfaces_m, only: &
        eval_electric_drift_velocities_c => &
            eval_electric_drift_velocities
    use kilca_legacy_interfaces_m, only: eval_fgi_arrays_c => eval_fgi_arrays
    use kilca_legacy_interfaces_m, only: calc_w2_array_c => calc_w2_array
    use kilca_legacy_interfaces_m, only: calc_d_array_c => calc_d_array
    use kilca_legacy_interfaces_m, only: calc_k_matrices_c => calc_k_matrices
    use kilca_legacy_interfaces_m, only: calc_k1_matrices_c => calc_k1_matrices
    use kilca_legacy_interfaces_m, only: &
        calc_and_add_galilelian_correction_c => &
            calc_and_add_galilelian_correction
    use kilca_legacy_interfaces_m, only: get_flre_order_c => get_flre_order
    use kilca_legacy_interfaces_m, only: get_gal_corr_c => get_gal_corr
    use adaptive_grid_pol_m, only: adaptive_grid_polynom_err
    use adaptive_grid_pol_m, only: adaptive_grid_polynom_res
    use kilca_background_data_m, only: get_background_obj_x0_c => get_background_x0
    use kilca_background_data_m, &
        only: get_background_obj_xlast_c => get_background_xlast
    use kilca_background_settings_m, &
        only: get_background_charge_c => get_background_charge
    use kilca_background_settings_m, &
        only: get_background_flag_back_c => get_background_flag_back
    use kilca_background_settings_m, only: &
        get_background_huge_factor_c => &
            get_background_huge_factor
    use kilca_background_settings_m, only: get_background_n_c => get_background_n
    use kilca_spline_m, only: spline_alloc_c => spline_alloc
    use kilca_spline_m, only: spline_calc_c => spline_calc
    use kilca_spline_m, only: spline_eval_d_c => spline_eval_d
    use kilca_spline_m, only: spline_free_c => spline_free
    use kilca_wave_data_m, &
        only: get_wave_data_obj_omov_im_c => get_wave_data_obj_omov_im
    use kilca_wave_data_m, &
        only: get_wave_data_obj_omov_re_c => get_wave_data_obj_omov_re
    use, intrinsic :: iso_c_binding, only: &
        c_int, c_intptr_t, c_double, c_char, c_ptr, c_loc, c_f_pointer, c_null_ptr
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
    public :: sample_cond_func_polynom
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

    !> The owning zone supplies a native output path, adaptive-grid bounds a/b,
    !> and the wave's resonance position and complex frequency components.
    function cond_profiles_create(path2linear_p, flreo, gal_corr, N, max_dim_c, &
         r1, r2, D, eps_out, eps_res, a, b, r_res, omov_re, omov_im, flag_debug, &
        flag) &
         result(handle)
        character(len=*), intent(in) :: path2linear_p
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
        cp%path2linear = path2linear_p

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
                                    ind_err(l) = iKa(cp, spec, ttype, p, q, i, j, &
        part, 0)
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
        stat_unused = adaptive_grid_polynom_res(sample_cond_func_polynom, c_loc(cp), &
                                                a, b, cp%dimK, cp%NK, cp%dimx, epso, &
                                                r_res, D, epsi, dim_err, ind_err, &
        cp%xt, cp%yt)
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

    subroutine cond_profiles_destroy(handle)
        integer(c_intptr_t), value :: handle
        type(cond_profiles_t), pointer :: cp

        if (handle == 0_c_intptr_t) return
        call handle_to_cp(handle, cp)
        if (cp%sidK /= 0) call spline_free_c(cp%sidK)
        if (cp%sidC /= 0) call spline_free_c(cp%sidC)
        deallocate (cp)
    end subroutine cond_profiles_destroy

    !> Per-point fallback retained for the statically referenced branch in
    !> conductivity.f90. Handles are passed by reference through native interfaces;
    !> background bounds and wave parameters come from their owning modules.
    subroutine calc_and_spline_conductivity_for_point(sd_ptr, bp_ptr, wd_ptr, &
         flag_back_p, r, cp_ptr)
        integer(c_intptr_t), intent(in) :: sd_ptr, bp_ptr, wd_ptr
        character(len=*), intent(in) :: flag_back_p
        real(c_double), intent(in) :: r
        integer(c_intptr_t), intent(out) :: cp_ptr

        type(cond_profiles_t), pointer :: cp
        type(c_ptr) :: wd_cptr
        integer(c_int) :: l, dim_err
        integer(c_int), allocatable :: ind_err(:)
        integer(c_int) :: spec, ttype, p, q, i, j, part
        real(c_double) :: delta, a, b, eps
        real(c_double) :: bp_x0, bp_xlast, omov_re, omov_im
        integer(c_int) :: stat_unused

        allocate (cp)

        wd_cptr = transfer(wd_ptr, wd_cptr)

        cp%flag_back = flag_back_p(1:1)
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

        bp_x0 = get_background_obj_x0_c()
        bp_xlast = get_background_obj_xlast_c()
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
                                    ind_err(l) = iKa(cp, spec, ttype, p, q, i, j, &
        part, 0)
                                    l = l + 1
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do

        eps = 0.0d0
        stat_unused = adaptive_grid_polynom_err(sample_cond_func_polynom, c_loc(cp), &
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

    !> Release the conductivity instance referenced by the caller's handle.
    subroutine delete_conductivity_profiles_f(cp_ptr)
        integer(c_intptr_t), intent(in) :: cp_ptr
        type(cond_profiles_t), pointer :: cp

        call handle_to_cp(cp_ptr, cp)
        if (cp%sidK /= 0) call spline_free_c(cp%sidK)
        if (cp%sidC /= 0) call spline_free_c(cp%sidC)
        deallocate (cp)
    end subroutine delete_conductivity_profiles_f

    !> Native adaptive-grid callback; f receives the sampled convergence vector.
    subroutine sample_cond_func_polynom(r, f, p)
        real(c_double), intent(in) :: r
        real(c_double), intent(out) :: f(*)
        type(c_ptr), value :: p
        type(cond_profiles_t), pointer :: cp
        integer(c_int) :: spec, ttype

        call c_f_pointer(p, cp)

        call eval_and_set_background_parameters_spec_independent_c(r, cp%flag_back)
        call eval_and_set_wave_parameters_c(r, cp%flag_back)
        call eval_a_matrix_c()
        call calc_dem_djmi_arrays_c(r)

        do spec = 0, 1
            call eval_and_set_background_parameters_spec_dependent_c(r, spec, &
        cp%flag_back)
            call eval_and_set_f0_parameters_nu_and_derivs_c(r, spec, cp%flag_back)
            call eval_electric_drift_velocities_c()
            call eval_fgi_arrays_c()
            call calc_w2_array_c(spec)
            call calc_d_array_c()

            ttype = 0
            block
                complex(c_double) :: kmat(3, 3, 0:cp%flreo, 0:cp%flreo)
                integer :: offset, count
                call calc_k_matrices_c(kmat)
                offset = iKa(cp, spec, ttype, 0, 0, 0, 0, 0, 0)
                count = 2*size(kmat)
                f(offset + 1:offset + count) = transfer(kmat, f(1), count)
            end block

            ttype = 1
            block
                complex(c_double) :: kmat(3, 3, 0:cp%flreo, 0:cp%flreo)
                integer :: offset, count
                call calc_k1_matrices_c(kmat)
                offset = iKa(cp, spec, ttype, 0, 0, 0, 0, 0, 0)
                count = 2*size(kmat)
                f(offset + 1:offset + count) = transfer(kmat, f(1), count)
            end block
        end do
    end subroutine sample_cond_func_polynom

    subroutine set_arrays_for_K_polynom(cp)
        type(cond_profiles_t), target, intent(inout) :: cp
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
        type(cond_profiles_t), target, intent(inout) :: cp
        integer(c_int) :: ierr

        call set_arrays_for_K_polynom(cp)

        call spline_alloc_c(cp%NK, 1, cp%dimx, c_loc(cp%x), c_loc(cp%CK), cp%sidK)
        call spline_calc_c(cp%sidK, c_loc(cp%K), 0, cp%dimK - 1, c_null_ptr, ierr)
    end subroutine calc_splines_for_K_polynom

    subroutine calc_splines_for_C(cp)
        type(cond_profiles_t), target, intent(inout) :: cp
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

        call spline_alloc_c(cp%NC, 1, cp%dimx, c_loc(cp%x), c_loc(cp%CC), cp%sidC)
        call spline_calc_c(cp%sidC, c_loc(cp%C), 0, cp%dimC - 1, c_null_ptr, ierr)
    end subroutine calc_splines_for_C

    !> Mirrors calc_C_matrices exactly: evaluates K-matrices and derivatives,
    !> combines them via the binomial-coefficient sum into the C matrix, scales
    !> by the conductivity prefactor, and applies the Galilean correction.
    subroutine calc_C_matrices(cp, spec, ttype, r, Cout)
        type(cond_profiles_t), target, intent(inout) :: cp
        integer(c_int), intent(in) :: spec, ttype
        real(c_double), intent(in) :: r
        real(c_double), target, intent(out) :: Cout(*)

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
            block
                complex(c_double) :: ct(3, 3, 0:2*flreo)
                ct = reshape(transfer(Cout(1:2*size(ct)), [(0.0d0, 0.0d0)]), &
        shape(ct))
                call calc_and_add_galilelian_correction_c(r, spec, cp%flag_back, ct)
                Cout(1:2*size(ct)) = transfer(ct, Cout(1), 2*size(ct))
            end block
        end if
    end subroutine calc_C_matrices

    !> Matches the original eval_K_matrices: evaluates the spline (with
    !> derivatives Dmin..Dmax) at r and undoes the huge_factor scaling applied
    !> when flag_back != 'f'.
    subroutine eval_K_matrices(cp, spec, ttype, Dmin, Dmax, r, Kout)
        type(cond_profiles_t), intent(in) :: cp
        integer(c_int), intent(in) :: spec, ttype, Dmin, Dmax
        real(c_double), intent(in) :: r
        real(c_double), target, intent(out) :: Kout(*)
        integer(c_int) :: dimk, ind_ka, n, ind, jj
        real(c_double) :: scale_fac
        real(c_double), target :: rloc(1)

        dimk = 2*9*(cp%flreo + 1)*(cp%flreo + 1)
        ind_ka = iKa(cp, spec, ttype, 0, 0, 0, 0, 0, 0)

        rloc(1) = r
        call spline_eval_d_c(cp%sidK, 1, c_loc(rloc), Dmin, Dmax, ind_ka, &
        ind_ka + dimk - 1, c_loc(Kout))

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
        real(c_double), target, intent(out) :: Cout(*)
        integer(c_int) :: dimc, ind_ca, n, ind, jj
        real(c_double) :: scale_fac
        real(c_double), target :: rloc(1)

        dimc = 18*(2*cp%flreo + 1)
        ind_ca = dimc*(ttype + cp%dimt*spec)

        rloc(1) = r
        call spline_eval_d_c(cp%sidC, 1, c_loc(rloc), Dmin, Dmax, ind_ca, &
        ind_ca + dimc - 1, c_loc(Cout))

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

    !> Native complex-matrix interface used by ctensor during RHS evaluation.
    subroutine eval_c_matrices_f(cp_ptr, spec, ttype, Dmin, Dmax, r, ct)
        integer(c_intptr_t), intent(in) :: cp_ptr
        integer(c_int), intent(in) :: spec, ttype, Dmin, Dmax
        real(c_double), intent(in) :: r
        complex(c_double), target, intent(out) :: ct(*)
        real(c_double), pointer :: flat(:)
        type(cond_profiles_t), pointer :: cp

        call handle_to_cp(cp_ptr, cp)
        call c_f_pointer(c_loc(ct), flat, [2*9*(2*cp%flreo + 1)*(Dmax - Dmin + 1)])
        call eval_C_matrices(cp, spec, ttype, Dmin, Dmax, r, flat)
    end subroutine eval_c_matrices_f

    !> Native complex-matrix interface used by kmatrices.
    subroutine eval_k_matrices_f(cp_ptr, spec, ttype, Dmin, Dmax, r, km)
        integer(c_intptr_t), intent(in) :: cp_ptr
        integer(c_int), intent(in) :: spec, ttype, Dmin, Dmax
        real(c_double), intent(in) :: r
        complex(c_double), target, intent(out) :: km(*)
        real(c_double), pointer :: flat(:)
        type(cond_profiles_t), pointer :: cp

        call handle_to_cp(cp_ptr, cp)
        call c_f_pointer(c_loc(km), flat, [2*9*(cp%flreo + 1)**2*(Dmax - Dmin + 1)])
        call eval_K_matrices(cp, spec, ttype, Dmin, Dmax, r, flat)
    end subroutine eval_k_matrices_f

    !> Matches eval_all_K_matrices: full dimK spline evaluation (no per-spec
    !> index offset), with the same huge_factor unscaling and a handle passed by value.
    subroutine eval_all_k_matrices(handle, Dmin, Dmax, r, Kout)
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: Dmin, Dmax
        real(c_double), value :: r
        real(c_double), target, intent(out) :: Kout(*)
        type(cond_profiles_t), pointer :: cp
        integer(c_int) :: dimk, n, ind, jj
        real(c_double) :: scale_fac
        real(c_double), target :: rloc(1)

        call handle_to_cp(handle, cp)
        dimk = cp%dimK
        rloc(1) = r
        call spline_eval_d_c(cp%sidK, 1, c_loc(rloc), Dmin, Dmax, 0, dimk - 1, &
        c_loc(Kout))

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
    subroutine eval_all_c_matrices(handle, Dmin, Dmax, r, Cout)
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: Dmin, Dmax
        real(c_double), value :: r
        real(c_double), target, intent(out) :: Cout(*)
        type(cond_profiles_t), pointer :: cp
        integer(c_int) :: dimc, n, ind, jj
        real(c_double) :: scale_fac
        real(c_double), target :: rloc(1)

        call handle_to_cp(handle, cp)
        dimc = cp%dimC
        rloc(1) = r
        call spline_eval_d_c(cp%sidC, 1, c_loc(rloc), Dmin, Dmax, 0, dimc - 1, &
        c_loc(Cout))

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
                    open (newunit=u, file=trim(fname), status='replace', &
        action='write')
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

    integer(c_int) function get_cond_nk(handle)
        integer(c_intptr_t), value :: handle
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        get_cond_nk = cp%NK
    end function get_cond_nk

    integer(c_int) function get_cond_dimk(handle)
        integer(c_intptr_t), value :: handle
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        get_cond_dimk = cp%dimK
    end function get_cond_dimk

    integer(c_int) function get_cond_nc(handle)
        integer(c_intptr_t), value :: handle
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        get_cond_nc = cp%NC
    end function get_cond_nc

    integer(c_int) function get_cond_dimc(handle)
        integer(c_intptr_t), value :: handle
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        get_cond_dimc = cp%dimC
    end function get_cond_dimc

    integer(c_int) function get_cond_dimx(handle)
        integer(c_intptr_t), value :: handle
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        get_cond_dimx = cp%dimx
    end function get_cond_dimx

    function get_cond_x_ptr(handle) result(ptr)
        integer(c_intptr_t), value :: handle
        type(c_ptr) :: ptr
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        ptr = c_loc(cp%x(0))
    end function get_cond_x_ptr

    function get_cond_k_ptr(handle) result(ptr)
        integer(c_intptr_t), value :: handle
        type(c_ptr) :: ptr
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        ptr = c_loc(cp%K(0))
    end function get_cond_k_ptr

    integer(c_int) function get_cond_iks(handle, spec, ttype, p, q, i, j, part, node)
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: spec, ttype, p, q, i, j, part, node
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        get_cond_iks = iKs(cp, spec, ttype, p, q, i, j, part, node)
    end function get_cond_iks

    function get_cond_flag_back(handle) result(ch)
        integer(c_intptr_t), value :: handle
        character(kind=c_char) :: ch
        type(cond_profiles_t), pointer :: cp
        call handle_to_cp(handle, cp)
        ch = cp%flag_back(1:1)
    end function get_cond_flag_back

    subroutine get_cond_path2linear(handle, buf)
        integer(c_intptr_t), value :: handle
        character(len=*), intent(out) :: buf
        type(cond_profiles_t), pointer :: cp
        integer :: i, n

        call handle_to_cp(handle, cp)
        buf = cp%path2linear
    end subroutine get_cond_path2linear

    subroutine handle_to_cp(handle, cp)
        integer(c_intptr_t), value :: handle
        type(cond_profiles_t), pointer, intent(out) :: cp
        type(c_ptr) :: ccp
        ccp = transfer(handle, ccp)
        call c_f_pointer(ccp, cp)
    end subroutine handle_to_cp

end module kilca_cond_profiles_m
