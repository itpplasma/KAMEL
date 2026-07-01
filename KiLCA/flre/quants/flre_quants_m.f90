!> Auxiliary FLRE quantities (current/power/flux/density/torque), formerly the
!> C++ flre_quants class. Per-instance handle (same pattern as
!> kilca_cond_profiles_m/kilca_sysmat_profiles_m), since each flre_zone owns
!> its own instance.
!>
!> The C++ class held a FIXED set of 8 quantities (num_tot=8, "add more
!> quants if needed" never exercised) dispatched through calc_quant[]/
!> save_quant[] function-pointer arrays populated once in the constructor.
!> Since the set is fixed and known at translation time, the dispatch is
!> rendered as direct conditional calls (matching exactly which function
!> fires for which quantity index and loop) instead of literal Fortran
!> procedure-pointer arrays -- behaviorally identical, far less error-prone
!> to hand-translate. For the same reason, qloc[i]/qint[i] (a C++ array of
!> per-quantity heap buffers, sized dimq[i]*dimx and allocated only for
!> requested quantities) becomes 8 separately-named, always-allocated flat
!> arrays: never reading an unallocated buffer was already guaranteed by the
!> oracle's own get_output_flag_quants_/flag[] gating, so always-allocating
!> changes memory footprint in the (untested) all-quantities-off corner case
!> but not observable output.
!>
!> All `if (DEBUG_FLAG) fprintf(...)` blocks from the oracle are dropped:
!> DEBUG_FLAG is a compile-time macro hard-coded 0 in code_settings.h, so
!> these blocks provably never execute (matching the SORT_DISPERSION_PROFILES/
!> USE_SPLINES_IN_RHS_EVALUATION precedent). save_total_absorbed_power
!> (declared in calc_flre_quants.h) and the file-local `binary_search` are
!> dropped: both have zero definitions/callers anywhere in the live tree.
!>
!> calc_quant/save_quant member functions are called ONLY from within this
!> class's own dispatch methods (never directly by flre_zone.cpp or other
!> files), so they are plain private module subroutines, not bind(C) -- only
!> the handful of true external entry points (create/destroy/
!> calculate_local_profiles/calculate_integrated_profiles/save_profiles/
!> calculate_JaE/transform_quants_to_lab_cyl_frame/interp_diss_power_density/
!> interp_current_density) need bind(C) names.
module kilca_flre_quants_m
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_double, c_char, &
        c_ptr, c_loc, c_f_pointer, c_null_ptr, c_null_char
    use kilca_cond_profiles_m, only: eval_all_k_matrices, eval_all_c_matrices, &
        get_cond_nc
    use kilca_maxwell_eqs_data_m, only: get_me_iersp_sys, get_me_ibrsp_sys, &
        get_me_der_order
    use kilca_neville_m, only: eval_neville_polynom
    implicit none
    private

    public :: flre_quants_create, flre_quants_destroy
    public :: flre_quants_calculate_local_profiles
    public :: flre_quants_calculate_integrated_profiles
    public :: flre_quants_save_profiles
    public :: flre_quants_calculate_jae
    public :: flre_quants_transform_quants_to_lab_cyl_frame
    public :: flre_quants_interp_diss_power_density
    public :: flre_quants_interp_current_density

    real(c_double), parameter :: pi = 3.141592653589793238462643383279502884197d0
    complex(c_double), parameter :: ii = (0.0d0, 1.0d0)
    real(c_double), parameter :: cspeed = 2.9979245800d10
    real(c_double), parameter :: echarge = 4.8032d-10
    integer(c_int), parameter :: boundary_antenna = 4

    !quants indices, 0-based to match get_output_flag_quants_(i):
    integer(c_int), parameter :: current_dens_q = 0
    integer(c_int), parameter :: abs_power_dens_q = 1
    integer(c_int), parameter :: diss_power_dens_q = 2
    integer(c_int), parameter :: kin_flux_q = 3
    integer(c_int), parameter :: poy_flux_q = 4
    integer(c_int), parameter :: tot_flux_q = 5
    integer(c_int), parameter :: number_dens_q = 6
    integer(c_int), parameter :: lor_torque_dens_q = 7

    type :: flre_quants_t
        integer(c_intptr_t) :: zone_cp = 0
        integer(c_intptr_t) :: zone_me = 0
        type(c_ptr) :: zone_bp = c_null_ptr
        character(len=1024) :: path2linear
        integer(c_int) :: zone_index
        integer(c_int) :: bc1, bc2
        real(c_double) :: vol_fac
        real(c_double), allocatable :: bico(:)
        integer(c_int) :: flreo
        integer(c_int) :: dimx
        real(c_double), pointer :: x(:) => null()
        integer(c_int) :: ncomps
        real(c_double), pointer :: eb_mov(:) => null()
        complex(c_double) :: wd_omov

        integer(c_int) :: numq
        integer(c_int) :: n_spline !spline degree, = NC of cond_profiles

        !per-node state:
        real(c_double) :: r
        integer(c_int) :: node
        logical :: flag_computed(0:7)
        logical :: flagC, flagK, flagS
        integer(c_int) :: dmax
        real(c_double), allocatable :: cmat(:)
        real(c_double), allocatable :: kmat(:)

        !spline for current density (jr) vs r:
        integer(c_intptr_t) :: sidY = 0
        integer(c_int) :: ny
        real(c_double), allocatable :: yarr(:), sarr(:)

        real(c_double) :: jae
        real(c_double), allocatable :: jae_arr(:), jaei_arr(:)
        real(c_double), allocatable :: cdlab(:)

        !per-quantity local/integrated storage (always allocated; see module doc):
        real(c_double), allocatable :: current_dens(:)
        real(c_double), allocatable :: abs_pow_dens(:), abs_pow_int(:)
        real(c_double), allocatable :: diss_pow_dens(:), diss_pow_int(:)
        real(c_double), allocatable :: kin_flux(:)
        real(c_double), allocatable :: poy_flux(:)
        real(c_double), allocatable :: tot_flux(:)
        real(c_double), allocatable :: number_dens(:)
        real(c_double), allocatable :: lor_torque_dens(:), lor_torque_int(:)
    end type flre_quants_t

    interface
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

        integer(c_int) function save_real_array_c(dim, xgrid, arr, full_name) &
            bind(C, name="save_real_array")
            import :: c_int, c_double, c_char
            integer(c_int), value :: dim
            real(c_double), intent(in) :: xgrid(*), arr(*)
            character(kind=c_char), intent(in) :: full_name(*)
        end function save_real_array_c

        subroutine current_density_c(jsurf) bind(C, name="current_density_")
            import :: c_double
            real(c_double), intent(out) :: jsurf(*)
        end subroutine current_density_c

        subroutine cyl2rsp_c(ra, jr, jt, js, jp) bind(C, name="cyl2rsp_")
            import :: c_double
            real(c_double), intent(in) :: ra, jr(*), jt(*)
            real(c_double), intent(out) :: js(*), jp(*)
        end subroutine cyl2rsp_c

        subroutine calc_current_density_r_s_p_c(r, jrsp) bind(C, name="calc_current_density_r_s_p_")
            import :: c_double
            real(c_double), intent(in) :: r
            real(c_double), intent(out) :: jrsp(*)
        end subroutine calc_current_density_r_s_p_c

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

        subroutine get_wave_parameters_c(kvals) bind(C, name="get_wave_parameters_")
            import :: c_double
            real(c_double), intent(out) :: kvals(*)
        end subroutine get_wave_parameters_c

        subroutine get_magnetic_field_parameters_c(hvals) bind(C, name="get_magnetic_field_parameters_")
            import :: c_double
            real(c_double), intent(out) :: hvals(*)
        end subroutine get_magnetic_field_parameters_c

        subroutine eval_hthz_c(rval, omin, omax, bp, hout) bind(C, name="eval_hthz_")
            import :: c_double, c_int, c_ptr
            real(c_double), intent(in) :: rval
            integer(c_int), intent(in) :: omin, omax
            type(c_ptr), intent(in) :: bp
            real(c_double), intent(out) :: hout(*)
        end subroutine eval_hthz_c

        function get_antenna_wa_c() result(wa) bind(C, name="get_antenna_wa_")
            import :: c_double
            real(c_double) :: wa
        end function get_antenna_wa_c

        function get_antenna_ra_c() result(ra) bind(C, name="get_antenna_ra_")
            import :: c_double
            real(c_double) :: ra
        end function get_antenna_ra_c

        function get_background_rtor_c() result(rtor) bind(C, name="get_background_rtor_")
            import :: c_double
            real(c_double) :: rtor
        end function get_background_rtor_c

        function get_background_v_gal_sys_c() result(v) bind(C, name="get_background_V_gal_sys_")
            import :: c_double
            real(c_double) :: v
        end function get_background_v_gal_sys_c

        function get_background_charge_c(i) result(ch) bind(C, name="get_background_charge_")
            import :: c_int, c_double
            integer(c_int), value :: i
            real(c_double) :: ch
        end function get_background_charge_c

        function get_background_flag_back_c() result(ch) bind(C, name="get_background_flag_back_")
            import :: c_char
            character(kind=c_char) :: ch
        end function get_background_flag_back_c

        function get_output_flag_quants_c(i) result(flag) bind(C, name="get_output_flag_quants_")
            import :: c_int
            integer(c_int), value :: i
            integer(c_int) :: flag
        end function get_output_flag_quants_c

        function get_output_flag_emfield_c() result(flag) bind(C, name="get_output_flag_emfield_")
            import :: c_int
            integer(c_int) :: flag
        end function get_output_flag_emfield_c

        function get_output_flag_additional_c() result(flag) bind(C, name="get_output_flag_additional_")
            import :: c_int
            integer(c_int) :: flag
        end function get_output_flag_additional_c

        function get_output_num_quants_c() result(n) bind(C, name="get_output_num_quants_")
            import :: c_int
            integer(c_int) :: n
        end function get_output_num_quants_c

        subroutine binomial_coefficients_c(N, BC) bind(C, name="binomial_coefficients_")
            import :: c_int, c_double
            integer(c_int), value :: N
            real(c_double), intent(out) :: BC(*)
        end subroutine binomial_coefficients_c
    end interface

contains

    !flat-index helpers, row-major (rightmost fastest), matching the C
    !typedef-cast layouts in calc_flre_quants.cpp exactly (verified against
    !the underlying spline-channel formulas in cond_profs_m's iKs/iCs and
    !shared.cpp's binomial_coefficients):

    integer(c_int) function idx_f(ncomps, node, comp, part) result(idx)
        integer(c_int), intent(in) :: ncomps, node, comp, part
        idx = part + 2*(comp + ncomps*node)
    end function idx_f

    integer(c_int) function idx_cd(dimx, spec, type_, i, part, node) result(idx)
        integer(c_int), intent(in) :: dimx, spec, type_, i, part, node
        idx = node + dimx*(part + 2*(i + 3*(type_ + 2*spec)))
    end function idx_cd

    integer(c_int) function idx_2(dimx, spec, type_, node) result(idx)
        integer(c_int), intent(in) :: dimx, spec, type_, node
        idx = node + dimx*type_ + dimx*2*spec
    end function idx_2

    integer(c_int) function idx_3(dimx, spec, type_, node) result(idx)
        integer(c_int), intent(in) :: dimx, spec, type_, node
        idx = node + dimx*type_ + dimx*3*spec
    end function idx_3

    integer(c_int) function idx_nd(dimx, spec, part, node) result(idx)
        integer(c_int), intent(in) :: dimx, spec, part, node
        idx = node + dimx*(part + 2*spec)
    end function idx_nd

    integer(c_int) function idx_kmat(dimk, flreo, deriv, spec, type_, p, q, i, j, part) result(idx)
        integer(c_int), intent(in) :: dimk, flreo, deriv, spec, type_, p, q, i, j, part
        idx = deriv*dimk + part + 2*(j + 3*(i + 3*(q + (flreo + 1)*(p + (flreo + 1)*( &
              type_ + 2*spec)))))
    end function idx_kmat

    integer(c_int) function idx_cmat(flreo, spec, type_, s, i, j, part) result(idx)
        integer(c_int), intent(in) :: flreo, spec, type_, s, i, j, part
        idx = part + 2*(j + 3*(i + 3*(s + (2*flreo + 1)*(type_ + 2*spec))))
    end function idx_cmat

    integer(c_int) function idx_bico(flreo, k, n) result(idx)
        integer(c_int), intent(in) :: flreo, k, n
        idx = n + k*(flreo + 1)
    end function idx_bico

    subroutine save_arr(dim, xgrid, arr, fname)
        integer(c_int), intent(in) :: dim
        real(c_double), intent(in) :: xgrid(*), arr(*)
        character(len=*), intent(in) :: fname
        character(kind=c_char) :: cbuf(len(fname) + 1)
        integer :: i, n, ierr_unused

        n = len_trim(fname)
        do i = 1, n
            cbuf(i) = fname(i:i)
        end do
        cbuf(n + 1) = c_null_char
        ierr_unused = save_real_array_c(dim, xgrid, arr, cbuf)
    end subroutine save_arr

    function itoa(n) result(s)
        integer(c_int), intent(in) :: n
        character(len=12) :: s
        write (s, '(i0)') n
    end function itoa

    function flre_quants_create(zone_cp, zone_me, zone_bp, path2linear_p, &
        flreo, dimx, x_ptr, ncomps, eb_mov_ptr, wd_omov_re, wd_omov_im, &
        bc1, bc2, zone_index) result(handle) bind(C, name="flre_quants_create_")
        integer(c_intptr_t), value :: zone_cp, zone_me
        type(c_ptr), value :: zone_bp, path2linear_p, x_ptr, eb_mov_ptr
        integer(c_int), value :: flreo, dimx, ncomps, bc1, bc2, zone_index
        real(c_double), value :: wd_omov_re, wd_omov_im
        integer(c_intptr_t) :: handle

        type(flre_quants_t), pointer :: qp
        integer(c_int) :: i, num_quants, ierr

        allocate (qp)
        qp%zone_cp = zone_cp
        qp%zone_me = zone_me
        qp%zone_bp = zone_bp
        call read_c_string(path2linear_p, qp%path2linear)
        qp%bc1 = bc1
        qp%bc2 = bc2
        qp%zone_index = zone_index
        qp%flreo = flreo
        qp%dimx = dimx
        call c_f_pointer(x_ptr, qp%x, [dimx])
        qp%ncomps = ncomps
        call c_f_pointer(eb_mov_ptr, qp%eb_mov, [2*ncomps*dimx])
        qp%wd_omov = cmplx(wd_omov_re, wd_omov_im, c_double)

        qp%vol_fac = (2.0d0*pi)*(2.0d0*pi)*get_background_rtor_c()

        allocate (qp%bico(0:(flreo + 1)*(flreo + 1) - 1))
        call binomial_coefficients_c(flreo, qp%bico)

        num_quants = get_output_num_quants_c()

        qp%numq = 0
        do i = 0, num_quants - 1
            if (get_output_flag_quants_c(i) /= 0) qp%numq = qp%numq + 1
        end do

        if (qp%numq == 0) then
            handle = transfer(c_loc(qp), handle)
            return
        end if

        ! 1-based (indices 1..N), NOT 0-based: every access site throughout
        ! this module uses a uniform idx_2/idx_3/idx_cd/idx_nd(...) + 1
        ! shift (the idx_* helpers themselves return a 0-based flat index
        ! in [0, N-1]; the +1 converts that into this array's own 1-based
        ! Fortran indexing). Allocating these 0-based instead put the
        ! valid range one past the array's declared upper bound.
        allocate (qp%current_dens(1:36*dimx))
        allocate (qp%abs_pow_dens(1:6*dimx), qp%abs_pow_int(1:6*dimx))
        allocate (qp%diss_pow_dens(1:9*dimx), qp%diss_pow_int(1:9*dimx))
        allocate (qp%kin_flux(1:9*dimx))
        allocate (qp%poy_flux(1:dimx))
        allocate (qp%tot_flux(1:dimx))
        allocate (qp%number_dens(1:6*dimx))
        allocate (qp%lor_torque_dens(1:9*dimx), qp%lor_torque_int(1:9*dimx))

        qp%n_spline = get_cond_nc(zone_cp)

        qp%flagC = .false.
        qp%flagK = .false.
        allocate (qp%cmat(0:(qp%n_spline + 1)*2*2*(2*flreo + 1)*3*3*2 - 1))
        ! 1-based: idx_kmat(...) + 1 is used at every kmat access site
        ! (unlike cmat/bico, which are accessed via idx_cmat/idx_bico with
        ! NO +1, correctly matching their own 0-based allocation).
        allocate (qp%kmat(1:(qp%n_spline + 1)*2*2*2*(flreo + 1)*(flreo + 1)*3*3*2))

        qp%ny = 2*1*1*3
        ! 1-based: calc_splines_for_current_density's ind+k+1 reaches up to
        ! 6*dimx (spec=2, part=1, k=dimx-1), one past a 0-based dimx*ny-1
        ! upper bound.
        allocate (qp%yarr(1:dimx*qp%ny))
        allocate (qp%sarr(0:(qp%n_spline + 1)*dimx*qp%ny - 1))
        call spline_alloc_c(qp%n_spline, 1, dimx, qp%x, qp%sarr, qp%sidY)
        qp%flagS = .false.

        ! 1-based (indices 1..dimx), NOT 0-based: every access site
        ! (calculate_jae, calculate_jae_delta/_distributed,
        ! flre_quants_calculate_jae) uses a uniform k+1 shift matching the
        ! oracle's 0-based qp->jaE[k]/qp->jaEi[qp->dimx-1] - the array
        ! itself must hold dimx elements at indices 1..dimx, not 0..dimx-1.
        allocate (qp%jae_arr(1:dimx))
        allocate (qp%jaei_arr(1:dimx))

        allocate (qp%cdlab(0:2*3*2*3*dimx - 1))

        handle = transfer(c_loc(qp), handle)
    end function flre_quants_create

    subroutine flre_quants_destroy(handle) bind(C, name="flre_quants_destroy_")
        integer(c_intptr_t), value :: handle
        type(flre_quants_t), pointer :: qp

        if (handle == 0_c_intptr_t) return
        call handle_to_qp(handle, qp)
        if (qp%sidY /= 0) call spline_free_c(qp%sidY)
        deallocate (qp)
    end subroutine flre_quants_destroy

    subroutine set_null_node_state(qp, k)
        type(flre_quants_t), intent(inout) :: qp
        integer(c_int), intent(in) :: k

        qp%node = k
        qp%r = qp%x(k + 1)
        qp%flag_computed = .false.
        qp%flagC = .false.
        qp%flagK = .false.
    end subroutine set_null_node_state

    subroutine flre_quants_calculate_local_profiles(handle) &
        bind(C, name="flre_quants_calculate_local_profiles_")
        integer(c_intptr_t), value :: handle
        type(flre_quants_t), pointer :: qp
        integer(c_int) :: k

        call handle_to_qp(handle, qp)
        if (qp%numq == 0) return

        qp%dmax = qp%flreo

        do k = 0, qp%dimx - 1
            call set_null_node_state(qp, k)

            call eval_all_k_matrices(qp%zone_cp, 0, qp%dmax, qp%r, qp%kmat)
            qp%flagK = .true.

            call eval_all_c_matrices(qp%zone_cp, 0, 0, qp%r, qp%cmat)
            qp%flagC = .true.

            if (get_output_flag_quants_c(current_dens_q) /= 0) call calc_current_density(qp)
            if (get_output_flag_quants_c(abs_power_dens_q) /= 0) call calc_absorbed_power_density(qp)
            if (get_output_flag_quants_c(diss_power_dens_q) /= 0) call calc_dissipated_power_density(qp)
            if (get_output_flag_quants_c(kin_flux_q) /= 0) call calc_kinetic_flux(qp)
            if (get_output_flag_quants_c(poy_flux_q) /= 0) call calc_poynting_flux(qp)
            if (get_output_flag_quants_c(tot_flux_q) /= 0) call calc_total_flux(qp)
        end do

        call calc_splines_for_current_density(qp)

        call calculate_field_profiles_poy_test(qp)

        do k = 0, qp%dimx - 1
            call set_null_node_state(qp, k)

            if (get_output_flag_quants_c(number_dens_q) /= 0) call calc_number_density(qp)
            if (get_output_flag_quants_c(lor_torque_dens_q) /= 0) call calc_lorentz_torque_density(qp)
        end do
    end subroutine flre_quants_calculate_local_profiles

    subroutine flre_quants_calculate_integrated_profiles(handle) &
        bind(C, name="flre_quants_calculate_integrated_profiles_")
        integer(c_intptr_t), value :: handle
        type(flre_quants_t), pointer :: qp

        call handle_to_qp(handle, qp)

        if (get_output_flag_quants_c(abs_power_dens_q) /= 0) call calc_absorbed_power_in_cylinder(qp)
        if (get_output_flag_quants_c(diss_power_dens_q) /= 0) call calc_dissipated_power_in_cylinder(qp)
        if (get_output_flag_quants_c(lor_torque_dens_q) /= 0) call calc_lorentz_torque_on_cylinder(qp)
    end subroutine flre_quants_calculate_integrated_profiles

    subroutine flre_quants_save_profiles(handle) bind(C, name="flre_quants_save_profiles_")
        integer(c_intptr_t), value :: handle
        type(flre_quants_t), pointer :: qp

        call handle_to_qp(handle, qp)

        if (get_output_flag_quants_c(current_dens_q) == 2) call save_current_density_basic(qp)
        if (get_output_flag_quants_c(abs_power_dens_q) == 2) call save_absorbed_power(qp)
        if (get_output_flag_quants_c(diss_power_dens_q) == 2) call save_dissipated_power(qp)
        if (get_output_flag_quants_c(kin_flux_q) == 2) call save_kinetic_flux(qp)
        if (get_output_flag_quants_c(poy_flux_q) == 2) call save_poynting_flux(qp)
        if (get_output_flag_quants_c(tot_flux_q) == 2) call save_total_flux(qp)
        if (get_output_flag_quants_c(number_dens_q) == 2) call save_number_density(qp)
        if (get_output_flag_quants_c(lor_torque_dens_q) == 2) call save_lorentz_torque(qp)
    end subroutine flre_quants_save_profiles

    !==================================================================
    ! current density
    !==================================================================

    subroutine calc_current_density(qp)
        type(flre_quants_t), intent(inout) :: qp
        integer(c_int) :: i, j, spec, type_, order
        complex(c_double) :: cd, cm, ef

        if (.not. qp%flagC) return

        do spec = 0, 1
            do type_ = 0, 1
                do i = 0, 2
                    cd = (0.0d0, 0.0d0)
                    do j = 0, 2
                        do order = 0, get_me_der_order(qp%zone_me, i, j)
                            cm = cmplx(qp%cmat(idx_cmat(qp%flreo, spec, type_, order, i, j, 0)), &
                                       qp%cmat(idx_cmat(qp%flreo, spec, type_, order, i, j, 1)), c_double)
                            ef = cmplx(qp%eb_mov(idx_f(qp%ncomps, qp%node, get_me_iersp_sys(qp%zone_me, j) + order, 0) + 1), &
                                       qp%eb_mov(idx_f(qp%ncomps, qp%node, &
                                                       get_me_iersp_sys(qp%zone_me, j) + order, 1) + 1), c_double)
                            cd = cd + cm*ef
                        end do
                    end do
                    qp%current_dens(idx_cd(qp%dimx, spec, type_, i, 0, qp%node) + 1) = real(cd, c_double)
                    qp%current_dens(idx_cd(qp%dimx, spec, type_, i, 1, qp%node) + 1) = aimag(cd)
                end do
            end do
        end do

        do type_ = 0, 1
            do i = 0, 2
                qp%current_dens(idx_cd(qp%dimx, 2, type_, i, 0, qp%node) + 1) = &
                    qp%current_dens(idx_cd(qp%dimx, 0, type_, i, 0, qp%node) + 1) + &
                    qp%current_dens(idx_cd(qp%dimx, 1, type_, i, 0, qp%node) + 1)
                qp%current_dens(idx_cd(qp%dimx, 2, type_, i, 1, qp%node) + 1) = &
                    qp%current_dens(idx_cd(qp%dimx, 0, type_, i, 1, qp%node) + 1) + &
                    qp%current_dens(idx_cd(qp%dimx, 1, type_, i, 1, qp%node) + 1)
            end do
        end do

        qp%flag_computed(current_dens_q) = .true.
    end subroutine calc_current_density

    subroutine save_current_density_basic(qp)
        type(flre_quants_t), intent(in) :: qp
        character(len=1) :: sort(0:2), comp(0:2)
        integer(c_int) :: i, spec, type_

        sort = ['i', 'e', 't']
        comp = ['r', 's', 'p']

        do spec = 0, 2
            do type_ = 0, 1
                do i = 0, 2
                    call save_arr(qp%dimx, qp%x, qp%current_dens(idx_cd(qp%dimx, spec, type_, i, 0, 0) + 1:), &
                                  trim(qp%path2linear)//'zone_'//trim(itoa(qp%zone_index))// &
                                  '_current_dens_'//comp(i)//'_'//trim(itoa(type_))//'_'//sort(spec)//'.dat')
                end do
            end do
        end do
    end subroutine save_current_density_basic

    !==================================================================
    ! absorbed power density
    !==================================================================

    subroutine calc_absorbed_power_density(qp)
        type(flre_quants_t), intent(inout) :: qp
        integer(c_int) :: spec, type_, i
        complex(c_double) :: cd, ef
        real(c_double) :: apd

        if (.not. qp%flag_computed(current_dens_q)) return

        do spec = 0, 2
            do type_ = 0, 1
                apd = 0.0d0
                do i = 0, 2
                    cd = cmplx(qp%current_dens(idx_cd(qp%dimx, spec, type_, i, 0, qp%node) + 1), &
                               qp%current_dens(idx_cd(qp%dimx, spec, type_, i, 1, qp%node) + 1), c_double)
                    ef = cmplx(qp%eb_mov(idx_f(qp%ncomps, qp%node, get_me_iersp_sys(qp%zone_me, i), 0) + 1), &
                               qp%eb_mov(idx_f(qp%ncomps, qp%node, get_me_iersp_sys(qp%zone_me, i), 1) + 1), c_double)
                    apd = apd + 0.5d0*real(cd*conjg(ef), c_double)
                end do
                qp%abs_pow_dens(idx_2(qp%dimx, spec, type_, qp%node) + 1) = apd
            end do
        end do

        qp%flag_computed(abs_power_dens_q) = .true.
    end subroutine calc_absorbed_power_density

    subroutine calc_absorbed_power_in_cylinder(qp)
        type(flre_quants_t), intent(inout) :: qp
        integer(c_int) :: spec, type_

        if (get_output_flag_quants_c(abs_power_dens_q) == 0) return

        do spec = 0, 2
            do type_ = 0, 1
                call integrate_over_cylinder(qp%dimx, qp%x, qp%abs_pow_dens(idx_2(qp%dimx, spec, type_, 0) + 1:), &
                                             qp%vol_fac, qp%abs_pow_int(idx_2(qp%dimx, spec, type_, 0) + 1:))
            end do
        end do
    end subroutine calc_absorbed_power_in_cylinder

    subroutine save_absorbed_power(qp)
        type(flre_quants_t), intent(in) :: qp
        character(len=1) :: sort(0:2)
        integer(c_int) :: spec, type_

        sort = ['i', 'e', 't']

        do spec = 0, 2
            do type_ = 0, 1
                call save_arr(qp%dimx, qp%x, qp%abs_pow_dens(idx_2(qp%dimx, spec, type_, 0) + 1:), &
                              trim(qp%path2linear)//'zone_'//trim(itoa(qp%zone_index))// &
                              '_abs_pow_dens_'//trim(itoa(type_))//'_'//sort(spec)//'.dat')
                call save_arr(qp%dimx, qp%x, qp%abs_pow_int(idx_2(qp%dimx, spec, type_, 0) + 1:), &
                              trim(qp%path2linear)//'zone_'//trim(itoa(qp%zone_index))// &
                              '_abs_pow_int_'//trim(itoa(type_))//'_'//sort(spec)//'.dat')
            end do
        end do
    end subroutine save_absorbed_power

    !==================================================================
    ! dissipated power density
    !==================================================================

    subroutine calc_dissipated_power_density(qp)
        type(flre_quants_t), intent(inout) :: qp
        integer(c_int) :: i, j, n1, n2, spec, type_, dimk
        complex(c_double) :: dpd, km, ef1, ef2, fac

        if (.not. qp%flagK) return

        dimk = 2*2*(qp%flreo + 1)*(qp%flreo + 1)*3*3*2

        do spec = 0, 1
            fac = pi*get_background_charge_c(spec)*get_background_charge_c(spec)/qp%wd_omov/qp%r

            do type_ = 0, 1
                dpd = (0.0d0, 0.0d0)
                do n1 = 0, qp%flreo
                    do n2 = 0, qp%flreo
                        do i = 0, 2
                            do j = 0, 2
                                km = cmplx(qp%kmat(idx_kmat(dimk, qp%flreo, 0, spec, type_, n1, n2, i, j, 0) + 1), &
                                           qp%kmat(idx_kmat(dimk, qp%flreo, 0, spec, type_, n1, n2, i, j, 1) + 1), c_double)
                                ef1 = cmplx(qp%eb_mov(idx_f(qp%ncomps, qp%node, get_me_iersp_sys(qp%zone_me, i) + n1, 0) + 1), &
                                            qp%eb_mov(idx_f(qp%ncomps, qp%node, &
                                                            get_me_iersp_sys(qp%zone_me, i) + n1, 1) + 1), c_double)
                                ef2 = cmplx(qp%eb_mov(idx_f(qp%ncomps, qp%node, get_me_iersp_sys(qp%zone_me, j) + n2, 0) + 1), &
                                            qp%eb_mov(idx_f(qp%ncomps, qp%node, &
                                                            get_me_iersp_sys(qp%zone_me, j) + n2, 1) + 1), c_double)
                                dpd = dpd + km*conjg(ef1)*ef2
                            end do
                        end do
                    end do
                end do
                qp%diss_pow_dens(idx_3(qp%dimx, spec, type_, qp%node) + 1) = real(fac*ii*dpd, c_double)
            end do
        end do

        do type_ = 0, 1
            qp%diss_pow_dens(idx_3(qp%dimx, 2, type_, qp%node) + 1) = &
                qp%diss_pow_dens(idx_3(qp%dimx, 0, type_, qp%node) + 1) + &
                qp%diss_pow_dens(idx_3(qp%dimx, 1, type_, qp%node) + 1)
        end do

        do spec = 0, 2
            qp%diss_pow_dens(idx_3(qp%dimx, spec, 2, qp%node) + 1) = &
                qp%diss_pow_dens(idx_3(qp%dimx, spec, 1, qp%node) + 1) - &
                qp%diss_pow_dens(idx_3(qp%dimx, spec, 0, qp%node) + 1)
        end do

        qp%flag_computed(diss_power_dens_q) = .true.
    end subroutine calc_dissipated_power_density

    subroutine calc_dissipated_power_in_cylinder(qp)
        type(flre_quants_t), intent(inout) :: qp
        integer(c_int) :: spec, type_

        if (get_output_flag_quants_c(diss_power_dens_q) == 0) return

        do spec = 0, 2
            do type_ = 0, 2
                call integrate_over_cylinder(qp%dimx, qp%x, qp%diss_pow_dens(idx_3(qp%dimx, spec, type_, 0) + 1:), &
                                             qp%vol_fac, qp%diss_pow_int(idx_3(qp%dimx, spec, type_, 0) + 1:))
            end do
        end do
    end subroutine calc_dissipated_power_in_cylinder

    subroutine save_dissipated_power(qp)
        type(flre_quants_t), intent(in) :: qp
        character(len=1) :: sort(0:2)
        integer(c_int) :: spec, type_

        sort = ['i', 'e', 't']

        do spec = 0, 2
            do type_ = 0, 2
                call save_arr(qp%dimx, qp%x, qp%diss_pow_dens(idx_3(qp%dimx, spec, type_, 0) + 1:), &
                              trim(qp%path2linear)//'zone_'//trim(itoa(qp%zone_index))// &
                              '_dis_pow_dens_'//trim(itoa(type_))//'_'//sort(spec)//'.dat')
                call save_arr(qp%dimx, qp%x, qp%diss_pow_int(idx_3(qp%dimx, spec, type_, 0) + 1:), &
                              trim(qp%path2linear)//'zone_'//trim(itoa(qp%zone_index))// &
                              '_dis_pow_int_'//trim(itoa(type_))//'_'//sort(spec)//'.dat')
            end do
        end do
    end subroutine save_dissipated_power

    !==================================================================
    ! kinetic flux
    !==================================================================

    subroutine calc_kinetic_flux(qp)
        type(flre_quants_t), intent(inout) :: qp
        integer(c_int) :: i, j, n1, n2, spec, type_, p, s, dimk
        complex(c_double) :: kf, km, ef1, ef2, fac
        real(c_double) :: coeff

        if (.not. qp%flagK) return

        dimk = 2*2*(qp%flreo + 1)*(qp%flreo + 1)*3*3*2

        do spec = 0, 1
            fac = pi*get_background_charge_c(spec)*get_background_charge_c(spec)/qp%wd_omov/qp%r

            do type_ = 0, 1
                kf = (0.0d0, 0.0d0)
                do n1 = 0, qp%flreo
                    do n2 = 0, qp%flreo
                        do p = 0, n1 - 1
                            do s = 0, n1 - p - 1
                                coeff = (-1.0d0)**(n1 + p)*qp%bico(idx_bico(qp%flreo, s, n1 - p - 1))
                                do i = 0, 2
                                    do j = 0, 2
                                        km = cmplx(qp%kmat(idx_kmat(dimk, qp%flreo, n1 - p - s - 1, &
                                                                    spec, type_, n1, n2, i, j, 0) + 1), &
                                                   qp%kmat(idx_kmat(dimk, qp%flreo, n1 - p - s - 1, &
                                                                    spec, type_, n1, n2, i, j, 1) + 1), c_double)
                                        ef1 = cmplx(qp%eb_mov(idx_f(qp%ncomps, qp%node, &
                                                                    get_me_iersp_sys(qp%zone_me, i) + p, 0) + 1), &
                                                    qp%eb_mov(idx_f(qp%ncomps, qp%node, &
                                                                    get_me_iersp_sys(qp%zone_me, i) + p, 1) + 1), c_double)
                                        ef2 = cmplx(qp%eb_mov(idx_f(qp%ncomps, qp%node, &
                                                                    get_me_iersp_sys(qp%zone_me, j) + n2 + s, 0) + 1), &
                                                    qp%eb_mov(idx_f(qp%ncomps, qp%node, &
                                                                    get_me_iersp_sys(qp%zone_me, j) + n2 + s, 1) + 1), c_double)
                                        kf = kf + coeff*km*conjg(ef1)*ef2
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
                qp%kin_flux(idx_3(qp%dimx, spec, type_, qp%node) + 1) = real(fac*ii*kf, c_double)*qp%vol_fac*qp%r
            end do
        end do

        do type_ = 0, 1
            qp%kin_flux(idx_3(qp%dimx, 2, type_, qp%node) + 1) = &
                qp%kin_flux(idx_3(qp%dimx, 0, type_, qp%node) + 1) + &
                qp%kin_flux(idx_3(qp%dimx, 1, type_, qp%node) + 1)
        end do

        do spec = 0, 2
            qp%kin_flux(idx_3(qp%dimx, spec, 2, qp%node) + 1) = &
                qp%kin_flux(idx_3(qp%dimx, spec, 1, qp%node) + 1) - &
                qp%kin_flux(idx_3(qp%dimx, spec, 0, qp%node) + 1)
        end do

        qp%flag_computed(kin_flux_q) = .true.
    end subroutine calc_kinetic_flux

    subroutine save_kinetic_flux(qp)
        type(flre_quants_t), intent(in) :: qp
        character(len=1) :: sort(0:2)
        integer(c_int) :: spec, type_

        sort = ['i', 'e', 't']

        do spec = 0, 2
            do type_ = 0, 2
                call save_arr(qp%dimx, qp%x, qp%kin_flux(idx_3(qp%dimx, spec, type_, 0) + 1:), &
                              trim(qp%path2linear)//'zone_'//trim(itoa(qp%zone_index))// &
                              '_kin_flux_'//trim(itoa(type_))//'_'//sort(spec)//'.dat')
            end do
        end do
    end subroutine save_kinetic_flux

    !==================================================================
    ! poynting flux, total flux
    !==================================================================

    subroutine calc_poynting_flux(qp)
        type(flre_quants_t), intent(inout) :: qp
        complex(c_double) :: es, ep, bs, bp_

        es = cmplx(qp%eb_mov(idx_f(qp%ncomps, qp%node, get_me_iersp_sys(qp%zone_me, 1), 0) + 1), &
                   qp%eb_mov(idx_f(qp%ncomps, qp%node, get_me_iersp_sys(qp%zone_me, 1), 1) + 1), c_double)
        ep = cmplx(qp%eb_mov(idx_f(qp%ncomps, qp%node, get_me_iersp_sys(qp%zone_me, 2), 0) + 1), &
                   qp%eb_mov(idx_f(qp%ncomps, qp%node, get_me_iersp_sys(qp%zone_me, 2), 1) + 1), c_double)
        bs = cmplx(qp%eb_mov(idx_f(qp%ncomps, qp%node, get_me_ibrsp_sys(qp%zone_me, 1), 0) + 1), &
                   qp%eb_mov(idx_f(qp%ncomps, qp%node, get_me_ibrsp_sys(qp%zone_me, 1), 1) + 1), c_double)
        bp_ = cmplx(qp%eb_mov(idx_f(qp%ncomps, qp%node, get_me_ibrsp_sys(qp%zone_me, 2), 0) + 1), &
                    qp%eb_mov(idx_f(qp%ncomps, qp%node, get_me_ibrsp_sys(qp%zone_me, 2), 1) + 1), c_double)

        qp%poy_flux(qp%node + 1) = qp%vol_fac*qp%r*cspeed/(4.0d0*pi)*0.5d0* &
                                   real(conjg(es)*bp_ - conjg(ep)*bs, c_double)

        qp%flag_computed(poy_flux_q) = .true.
    end subroutine calc_poynting_flux

    subroutine save_poynting_flux(qp)
        type(flre_quants_t), intent(in) :: qp

        call save_arr(qp%dimx, qp%x, qp%poy_flux, &
                      trim(qp%path2linear)//'zone_'//trim(itoa(qp%zone_index))//'_poy_flux.dat')
    end subroutine save_poynting_flux

    subroutine calc_total_flux(qp)
        type(flre_quants_t), intent(inout) :: qp

        if (.not. (qp%flag_computed(kin_flux_q) .and. qp%flag_computed(poy_flux_q))) return

        qp%tot_flux(qp%node + 1) = qp%poy_flux(qp%node + 1) + qp%kin_flux(idx_3(qp%dimx, 2, 0, qp%node) + 1)

        qp%flag_computed(tot_flux_q) = .true.
    end subroutine calc_total_flux

    subroutine save_total_flux(qp)
        type(flre_quants_t), intent(in) :: qp

        call save_arr(qp%dimx, qp%x, qp%tot_flux, &
                      trim(qp%path2linear)//'zone_'//trim(itoa(qp%zone_index))//'_tot_flux.dat')
    end subroutine save_total_flux

    subroutine calculate_field_profiles_poy_test(qp)
        type(flre_quants_t), intent(inout) :: qp
        integer(c_intptr_t) :: sidpf
        real(c_double), allocatable :: spf(:), err(:)
        real(c_double) :: divpfd(1), avrg_err, max_err
        integer(c_int) :: k, ierr

        if (.not. (qp%flag_computed(abs_power_dens_q) .and. qp%flag_computed(poy_flux_q))) return

        allocate (spf(0:(qp%n_spline + 1)*qp%dimx - 1))
        call spline_alloc_c(qp%n_spline, 1, qp%dimx, qp%x, spf, sidpf)
        call spline_calc_c(sidpf, qp%poy_flux, 0, 0, c_null_ptr, ierr)

        allocate (err(0:qp%dimx - 1))

        avrg_err = 0.0d0
        max_err = 0.0d0

        do k = 0, qp%dimx - 2
            qp%node = k
            qp%r = qp%x(k + 1)

            call spline_eval_d_c(sidpf, 1, qp%x(k + 1:k + 1), 1, 1, 0, 0, divpfd)

            divpfd(1) = divpfd(1)/(-qp%r*qp%vol_fac)

            err(k) = abs(divpfd(1) - (qp%abs_pow_dens(idx_2(qp%dimx, 2, 0, k) + 1) + qp%jae_arr(k + 1)))/ &
                     max(1.0d0, abs(divpfd(1)))

            if (err(k) > max_err) max_err = err(k)

            avrg_err = avrg_err + err(k)
        end do

        avrg_err = avrg_err/(qp%dimx - 1)

        if (get_output_flag_additional_c() > 1) then
            call save_arr(qp%dimx - 1, qp%x, err, &
                          trim(qp%path2linear)//'zone_'//trim(itoa(qp%zone_index))//'_poy_test_err.dat')
        end if

        call spline_free_c(sidpf)
        deallocate (spf)
        deallocate (err)
    end subroutine calculate_field_profiles_poy_test

    !==================================================================
    ! current density splines, number density, lorentz torque
    !==================================================================

    subroutine calc_splines_for_current_density(qp)
        type(flre_quants_t), intent(inout) :: qp
        integer(c_int) :: spec, part, ind, k, ierr
        integer(c_int), parameter :: type_ = 0, comp = 0

        if (get_output_flag_quants_c(current_dens_q) == 0) return

        do spec = 0, 2
            do part = 0, 1
                ind = qp%dimx*(part + 2*spec)
                do k = 0, qp%dimx - 1
                    qp%yarr(ind + k + 1) = qp%current_dens(idx_cd(qp%dimx, spec, type_, comp, part, k) + 1)
                end do
            end do
        end do

        call spline_calc_c(qp%sidY, qp%yarr, 0, qp%ny - 1, c_null_ptr, ierr)

        qp%flagS = .true.
    end subroutine calc_splines_for_current_density

    subroutine calc_number_density(qp)
        type(flre_quants_t), intent(inout) :: qp
        character(kind=c_char) :: flag_back_buf(1)
        real(c_double) :: kvals(3), djr(2)
        complex(c_double) :: j(0:2), dj, nd
        integer(c_int) :: spec, i, imin, imax

        if (.not. qp%flagS) return

        flag_back_buf(1) = get_background_flag_back_c()

        call eval_and_set_background_parameters_spec_independent_c(qp%r, flag_back_buf, 1_c_int)
        call eval_and_set_wave_parameters_c(qp%r, flag_back_buf, 1_c_int)
        call get_wave_parameters_c(kvals)

        do spec = 0, 1
            do i = 0, 2
                j(i) = cmplx(qp%current_dens(idx_cd(qp%dimx, spec, 0, i, 0, qp%node) + 1), &
                             qp%current_dens(idx_cd(qp%dimx, spec, 0, i, 1, qp%node) + 1), c_double)
            end do

            imin = 2*spec
            imax = imin + 1

            call spline_eval_d_c(qp%sidY, 1, qp%x(qp%node + 1:qp%node + 1), 1, 1, imin, imax, djr)

            dj = cmplx(djr(1), djr(2), c_double)

            nd = -ii/qp%wd_omov/get_background_charge_c(spec)* &
                 (j(0)/qp%r + dj + ii*(kvals(2)*j(1) + kvals(3)*j(2)))

            qp%number_dens(idx_nd(qp%dimx, spec, 0, qp%node) + 1) = real(nd, c_double)
            qp%number_dens(idx_nd(qp%dimx, spec, 1, qp%node) + 1) = aimag(nd)
        end do

        qp%number_dens(idx_nd(qp%dimx, 2, 0, qp%node) + 1) = &
            qp%number_dens(idx_nd(qp%dimx, 0, 0, qp%node) + 1)*(get_background_charge_c(0)/echarge) + &
            qp%number_dens(idx_nd(qp%dimx, 1, 0, qp%node) + 1)*(get_background_charge_c(1)/echarge)
        qp%number_dens(idx_nd(qp%dimx, 2, 1, qp%node) + 1) = &
            qp%number_dens(idx_nd(qp%dimx, 0, 1, qp%node) + 1)*(get_background_charge_c(0)/echarge) + &
            qp%number_dens(idx_nd(qp%dimx, 1, 1, qp%node) + 1)*(get_background_charge_c(1)/echarge)

        qp%flag_computed(number_dens_q) = .true.
    end subroutine calc_number_density

    subroutine save_number_density(qp)
        type(flre_quants_t), intent(in) :: qp
        character(len=1) :: sort(0:2)
        integer(c_int) :: spec

        sort = ['i', 'e', 't']

        do spec = 0, 2
            call save_arr(qp%dimx, qp%x, qp%number_dens(idx_nd(qp%dimx, spec, 0, 0) + 1:), &
                          trim(qp%path2linear)//'zone_'//trim(itoa(qp%zone_index))//'_density_re_'//sort(spec)//'.dat')
        end do
    end subroutine save_number_density

    subroutine calc_lorentz_torque_density(qp)
        type(flre_quants_t), intent(inout) :: qp
        real(c_double) :: h(3)
        character(kind=c_char) :: flag_back_buf(1)
        complex(c_double) :: j_rsp(0:2), e_rsp(0:2), b_rsp(0:2)
        complex(c_double) :: j_cyl(0:2), e_cyl(0:2), bc_cyl(0:2), jxbc(0:2)
        complex(c_double) :: nd
        integer(c_int) :: spec, i

        if (.not. qp%flag_computed(number_dens_q)) return

        flag_back_buf(1) = get_background_flag_back_c()
        call eval_and_set_background_parameters_spec_independent_c(qp%r, flag_back_buf, 1_c_int)
        call get_magnetic_field_parameters_c(h)

        do spec = 0, 1
            do i = 0, 2
                j_rsp(i) = cmplx(qp%current_dens(idx_cd(qp%dimx, spec, 0, i, 0, qp%node) + 1), &
                                 qp%current_dens(idx_cd(qp%dimx, spec, 0, i, 1, qp%node) + 1), c_double)
                e_rsp(i) = cmplx(qp%eb_mov(idx_f(qp%ncomps, qp%node, get_me_iersp_sys(qp%zone_me, i), 0) + 1), &
                                 qp%eb_mov(idx_f(qp%ncomps, qp%node, get_me_iersp_sys(qp%zone_me, i), 1) + 1), c_double)
                b_rsp(i) = cmplx(qp%eb_mov(idx_f(qp%ncomps, qp%node, get_me_ibrsp_sys(qp%zone_me, i), 0) + 1), &
                                 qp%eb_mov(idx_f(qp%ncomps, qp%node, get_me_ibrsp_sys(qp%zone_me, i), 1) + 1), c_double)
            end do

            j_cyl(0) = j_rsp(0)
            j_cyl(1) = h(3)*j_rsp(1) + h(2)*j_rsp(2)
            j_cyl(2) = h(3)*j_rsp(2) - h(2)*j_rsp(1)

            e_cyl(0) = e_rsp(0)
            e_cyl(1) = h(3)*e_rsp(1) + h(2)*e_rsp(2)
            e_cyl(2) = h(3)*e_rsp(2) - h(2)*e_rsp(1)

            bc_cyl(0) = conjg(b_rsp(0))
            bc_cyl(1) = conjg(h(3)*b_rsp(1) + h(2)*b_rsp(2))
            bc_cyl(2) = conjg(h(3)*b_rsp(2) - h(2)*b_rsp(1))

            call vec_product_3d(j_cyl, bc_cyl, jxbc)

            nd = cmplx(qp%number_dens(idx_nd(qp%dimx, spec, 0, qp%node) + 1), &
                       qp%number_dens(idx_nd(qp%dimx, spec, 1, qp%node) + 1), c_double)

            do i = 0, 2
                qp%lor_torque_dens(idx_3(qp%dimx, spec, i, qp%node) + 1) = &
                    0.5d0*real(get_background_charge_c(spec)*nd*conjg(e_cyl(i)) + &
                              echarge/cspeed*jxbc(i), c_double)
            end do

            qp%lor_torque_dens(idx_3(qp%dimx, spec, 1, qp%node) + 1) = &
                qp%lor_torque_dens(idx_3(qp%dimx, spec, 1, qp%node) + 1)*qp%r
            qp%lor_torque_dens(idx_3(qp%dimx, spec, 2, qp%node) + 1) = &
                qp%lor_torque_dens(idx_3(qp%dimx, spec, 2, qp%node) + 1)*get_background_rtor_c()
        end do

        do i = 0, 2
            qp%lor_torque_dens(idx_3(qp%dimx, 2, i, qp%node) + 1) = &
                qp%lor_torque_dens(idx_3(qp%dimx, 0, i, qp%node) + 1) + &
                qp%lor_torque_dens(idx_3(qp%dimx, 1, i, qp%node) + 1)
        end do

        qp%flag_computed(lor_torque_dens_q) = .true.
    end subroutine calc_lorentz_torque_density

    subroutine calc_lorentz_torque_on_cylinder(qp)
        type(flre_quants_t), intent(inout) :: qp
        integer(c_int) :: spec, i
        real(c_double) :: factor

        if (get_output_flag_quants_c(lor_torque_dens_q) == 0) return

        do spec = 0, 2
            do i = 0, 2
                if (i == 0) then
                    factor = 0.0d0
                else
                    factor = qp%vol_fac
                end if
                call integrate_over_cylinder(qp%dimx, qp%x, qp%lor_torque_dens(idx_3(qp%dimx, spec, i, 0) + 1:), &
                                             factor, qp%lor_torque_int(idx_3(qp%dimx, spec, i, 0) + 1:))
            end do
        end do
    end subroutine calc_lorentz_torque_on_cylinder

    subroutine save_lorentz_torque(qp)
        type(flre_quants_t), intent(in) :: qp
        character(len=1) :: sort(0:2), comp(0:2)
        integer(c_int) :: spec, i

        sort = ['i', 'e', 't']
        comp = ['r', 't', 'z']

        do spec = 0, 2
            do i = 0, 2
                call save_arr(qp%dimx, qp%x, qp%lor_torque_dens(idx_3(qp%dimx, spec, i, 0) + 1:), &
                              trim(qp%path2linear)//'zone_'//trim(itoa(qp%zone_index))// &
                              '_torque_dens_'//comp(i)//'_'//sort(spec)//'.dat')
                call save_arr(qp%dimx, qp%x, qp%lor_torque_int(idx_3(qp%dimx, spec, i, 0) + 1:), &
                              trim(qp%path2linear)//'zone_'//trim(itoa(qp%zone_index))// &
                              '_torque_int_'//comp(i)//'_'//sort(spec)//'.dat')
            end do
        end do
    end subroutine save_lorentz_torque

    !==================================================================
    ! helpers
    !==================================================================

    subroutine vec_product_3d(a, b, res)
        complex(c_double), intent(in) :: a(0:2), b(0:2)
        complex(c_double), intent(out) :: res(0:2)
        res(0) = a(1)*b(2) - a(2)*b(1)
        res(1) = -(a(0)*b(2) - a(2)*b(0))
        res(2) = a(0)*b(1) - a(1)*b(0)
    end subroutine vec_product_3d

    subroutine integrate_over_cylinder(dim, x, q, vol_fac, qi)
        integer(c_int), intent(in) :: dim
        real(c_double), intent(in) :: x(*), q(*), vol_fac
        real(c_double), intent(out) :: qi(*)
        integer(c_int) :: i

        qi(1) = 0.0d0
        do i = 2, dim
            qi(i) = qi(i - 1) + vol_fac*0.5d0*(x(i) - x(i - 1))*(x(i - 1)*q(i - 1) + x(i)*q(i))
        end do
    end subroutine integrate_over_cylinder

    !==================================================================
    ! lab-frame / cylindrical-system current density transform
    !==================================================================

    subroutine flre_quants_transform_quants_to_lab_cyl_frame(handle) &
        bind(C, name="flre_quants_transform_quants_to_lab_cyl_frame_")
        integer(c_intptr_t), value :: handle
        type(flre_quants_t), pointer :: qp
        real(c_double), allocatable :: cd(:)

        call handle_to_qp(handle, qp)

        if (get_output_flag_quants_c(current_dens_q) > 0) then
            call eval_current_dens_in_lab_frame(qp, qp%cdlab)
        end if

        if (get_output_flag_quants_c(current_dens_q) > 1) then
            allocate (cd(0:2*3*2*3*qp%dimx - 1))

            call eval_current_dens_in_lab_frame(qp, cd)
            call save_current_density_ext(qp, cd, 'lab', 'rsp')

            call eval_current_dens_in_cyl_sys(qp, cd)
            call save_current_density_ext(qp, cd, 'lab', 'rtz')

            deallocate (cd)
        end if
    end subroutine flre_quants_transform_quants_to_lab_cyl_frame

    subroutine eval_current_dens_in_lab_frame(qp, cd)
        type(flre_quants_t), intent(in) :: qp
        real(c_double), intent(inout) :: cd(0:)
        complex(c_double) :: jj
        real(c_double) :: htz(2), vel(0:2)
        integer(c_int) :: k, type_, i, spec

        do k = 0, qp%dimx - 1
            call eval_hthz_c(qp%x(k + 1), 0, 0, qp%zone_bp, htz)

            vel(0) = 0.0d0
            vel(1) = -htz(1)*get_background_v_gal_sys_c()
            vel(2) = htz(2)*get_background_v_gal_sys_c()

            do type_ = 0, 1
                do i = 0, 2
                    do spec = 0, 1
                        jj = cmplx(qp%current_dens(idx_cd(qp%dimx, spec, type_, i, 0, k) + 1), &
                                   qp%current_dens(idx_cd(qp%dimx, spec, type_, i, 1, k) + 1), c_double) + &
                             get_background_charge_c(spec)*vel(i)* &
                             cmplx(qp%number_dens(idx_nd(qp%dimx, spec, 0, k) + 1), &
                                   qp%number_dens(idx_nd(qp%dimx, spec, 1, k) + 1), c_double)

                        cd(idx_cd(qp%dimx, spec, type_, i, 0, k)) = real(jj, c_double)
                        cd(idx_cd(qp%dimx, spec, type_, i, 1, k)) = aimag(jj)
                    end do
                    cd(idx_cd(qp%dimx, 2, type_, i, 0, k)) = cd(idx_cd(qp%dimx, 0, type_, i, 0, k)) + &
                                                             cd(idx_cd(qp%dimx, 1, type_, i, 0, k))
                    cd(idx_cd(qp%dimx, 2, type_, i, 1, k)) = cd(idx_cd(qp%dimx, 0, type_, i, 1, k)) + &
                                                             cd(idx_cd(qp%dimx, 1, type_, i, 1, k))
                end do
            end do
        end do
    end subroutine eval_current_dens_in_lab_frame

    subroutine eval_current_dens_in_cyl_sys(qp, cd)
        type(flre_quants_t), intent(in) :: qp
        real(c_double), intent(inout) :: cd(0:)
        complex(c_double) :: jrsp(0:2), jcyl(0:2)
        real(c_double) :: htz(2)
        integer(c_int) :: k, type_, spec, i

        do k = 0, qp%dimx - 1
            call eval_hthz_c(qp%x(k + 1), 0, 0, qp%zone_bp, htz)

            do type_ = 0, 1
                do spec = 0, 2
                    do i = 0, 2
                        jrsp(i) = cmplx(cd(idx_cd(qp%dimx, spec, type_, i, 0, k)), &
                                        cd(idx_cd(qp%dimx, spec, type_, i, 1, k)), c_double)
                    end do

                    jcyl(0) = jrsp(0)
                    jcyl(1) = htz(2)*jrsp(1) + htz(1)*jrsp(2)
                    jcyl(2) = -htz(1)*jrsp(1) + htz(2)*jrsp(2)

                    do i = 0, 2
                        cd(idx_cd(qp%dimx, spec, type_, i, 0, k)) = real(jcyl(i), c_double)
                        cd(idx_cd(qp%dimx, spec, type_, i, 1, k)) = aimag(jcyl(i))
                    end do
                end do
            end do
        end do
    end subroutine eval_current_dens_in_cyl_sys

    subroutine save_current_density_ext(qp, cd, frame, comp_set)
        type(flre_quants_t), intent(in) :: qp
        real(c_double), intent(in) :: cd(0:)
        character(len=*), intent(in) :: frame
        character(len=3), intent(in) :: comp_set
        character(len=1) :: sort(0:2)
        integer(c_int) :: spec, type_, i

        sort = ['i', 'e', 't']

        do spec = 0, 2
            do type_ = 0, 1
                do i = 0, 2
                    call save_arr(qp%dimx, qp%x, cd(idx_cd(qp%dimx, spec, type_, i, 0, 0):), &
                                  trim(qp%path2linear)//'zone_'//trim(itoa(qp%zone_index))// &
                                  '_current_dens_'//comp_set(i + 1:i + 1)//'_'//trim(itoa(type_))// &
                                  '_'//sort(spec)//'_'//trim(frame)//'.dat')
                end do
            end do
        end do
    end subroutine save_current_density_ext

    !==================================================================
    ! interpolation entry points (called externally from flre_zone.cpp)
    !==================================================================

    subroutine flre_quants_interp_diss_power_density(handle, x, type_, spec, dpd) &
        bind(C, name="flre_quants_interp_diss_power_density_")
        integer(c_intptr_t), value :: handle
        real(c_double), value :: x
        integer(c_int), value :: type_, spec
        real(c_double), intent(out) :: dpd(1)
        type(flre_quants_t), pointer :: qp
        integer(c_int) :: ind, deg

        call handle_to_qp(handle, qp)

        deg = 5
        ind = qp%dimx/2

        call eval_neville_polynom(qp%dimx, qp%x, qp%diss_pow_dens(idx_3(qp%dimx, spec, type_, 0) + 1:), &
                                  deg, x, 0, 0, ind, dpd)
    end subroutine flre_quants_interp_diss_power_density

    subroutine flre_quants_interp_current_density(handle, x, type_, spec, comp, jout) &
        bind(C, name="flre_quants_interp_current_density_")
        integer(c_intptr_t), value :: handle
        real(c_double), value :: x
        integer(c_int), value :: type_, spec, comp
        real(c_double), intent(out) :: jout(2)
        type(flre_quants_t), pointer :: qp
        integer(c_int) :: ind, deg

        call handle_to_qp(handle, qp)

        deg = 5
        ind = qp%dimx/2

        call eval_neville_polynom(qp%dimx, qp%x, qp%cdlab(idx_cd(qp%dimx, spec, type_, comp, 0, 0) + 1:), &
                                  deg, x, 0, 0, ind, jout(1))
        ind = qp%dimx/2
        call eval_neville_polynom(qp%dimx, qp%x, qp%cdlab(idx_cd(qp%dimx, spec, type_, comp, 1, 0) + 1:), &
                                  deg, x, 0, 0, ind, jout(2))
    end subroutine flre_quants_interp_current_density

    !==================================================================
    ! absorbed energy from antenna current (work of E field on Ja)
    !==================================================================

    subroutine flre_quants_calculate_jae(handle) bind(C, name="flre_quants_calculate_jae_")
        integer(c_intptr_t), value :: handle
        type(flre_quants_t), pointer :: qp
        integer(c_int) :: k

        call handle_to_qp(handle, qp)

        do k = 0, qp%dimx - 1
            qp%jae_arr(k + 1) = 0.0d0
            qp%jaei_arr(k + 1) = 0.0d0
        end do

        if (get_antenna_wa_c() == 0.0d0) then
            call calculate_jae_delta(qp)
        else
            call calculate_jae_distributed(qp)
        end if
    end subroutine flre_quants_calculate_jae

    subroutine calculate_jae_delta(qp)
        type(flre_quants_t), intent(inout) :: qp
        integer(c_int) :: ia
        real(c_double) :: jsurf(4), jsurft(4), antenna_ra
        complex(c_double) :: ja(0:1), ef(0:1)
        integer(c_int) :: iersp_sys(0:2)
        real(c_double) :: es_re, es_im, ep_re, ep_im

        if (qp%bc1 == boundary_antenna) then
            ia = 0
        else if (qp%bc2 == boundary_antenna) then
            ia = qp%dimx - 1
        else
            return
        end if

        call current_density_c(jsurf)
        antenna_ra = get_antenna_ra_c()
        call cyl2rsp_c(antenna_ra, jsurf, jsurf(3:), jsurft, jsurft(3:))

        ja(0) = cmplx(jsurft(1), jsurft(2), c_double)
        ja(1) = cmplx(jsurft(3), jsurft(4), c_double)

        iersp_sys(0) = get_me_iersp_sys(qp%zone_me, 0)
        iersp_sys(1) = get_me_iersp_sys(qp%zone_me, 1)
        iersp_sys(2) = get_me_iersp_sys(qp%zone_me, 2)

        es_re = qp%eb_mov(idx_f(qp%ncomps, ia, iersp_sys(1), 0) + 1)
        es_im = qp%eb_mov(idx_f(qp%ncomps, ia, iersp_sys(1), 1) + 1)
        ep_re = qp%eb_mov(idx_f(qp%ncomps, ia, iersp_sys(2), 0) + 1)
        ep_im = qp%eb_mov(idx_f(qp%ncomps, ia, iersp_sys(2), 1) + 1)

        ef(0) = cmplx(es_re, es_im, c_double)
        ef(1) = cmplx(ep_re, ep_im, c_double)

        qp%jae = 0.5d0*qp%vol_fac*qp%x(ia + 1)*real(ja(0)*conjg(ef(0)) + ja(1)*conjg(ef(1)), c_double)

        if (get_output_flag_emfield_c() > 1) then
            call save_arr(1, qp%x(ia + 1:ia + 1), [qp%jae], &
                          trim(qp%path2linear)//'zone_'//trim(itoa(qp%zone_index))//'_JaE.dat')
        end if
    end subroutine calculate_jae_delta

    subroutine calculate_jae_distributed(qp)
        type(flre_quants_t), intent(inout) :: qp
        real(c_double) :: ja_rsp(4)
        complex(c_double) :: ja(0:1), ef(0:1)
        integer(c_int) :: iersp_sys(0:2)
        integer(c_int) :: k
        real(c_double) :: es_re, es_im, ep_re, ep_im

        iersp_sys(0) = get_me_iersp_sys(qp%zone_me, 0)
        iersp_sys(1) = get_me_iersp_sys(qp%zone_me, 1)
        iersp_sys(2) = get_me_iersp_sys(qp%zone_me, 2)

        do k = 0, qp%dimx - 1
            call calc_current_density_r_s_p_c(qp%x(k + 1), ja_rsp)

            ja(0) = cmplx(ja_rsp(1), ja_rsp(2), c_double)
            ja(1) = cmplx(ja_rsp(3), ja_rsp(4), c_double)

            es_re = qp%eb_mov(idx_f(qp%ncomps, k, iersp_sys(1), 0) + 1)
            es_im = qp%eb_mov(idx_f(qp%ncomps, k, iersp_sys(1), 1) + 1)
            ep_re = qp%eb_mov(idx_f(qp%ncomps, k, iersp_sys(2), 0) + 1)
            ep_im = qp%eb_mov(idx_f(qp%ncomps, k, iersp_sys(2), 1) + 1)

            ef(0) = cmplx(es_re, es_im, c_double)
            ef(1) = cmplx(ep_re, ep_im, c_double)

            qp%jae_arr(k + 1) = 0.5d0*real(ja(0)*conjg(ef(0)) + ja(1)*conjg(ef(1)), c_double)
        end do

        call integrate_over_cylinder(qp%dimx, qp%x, qp%jae_arr, qp%vol_fac, qp%jaei_arr)

        qp%jae = qp%jaei_arr(qp%dimx)

        if (get_output_flag_emfield_c() > 1) then
            call save_arr(1, qp%x(qp%dimx:qp%dimx), [qp%jae], &
                          trim(qp%path2linear)//'zone_'//trim(itoa(qp%zone_index))//'_JaE.dat')
            call save_arr(qp%dimx, qp%x, qp%jae_arr, &
                          trim(qp%path2linear)//'zone_'//trim(itoa(qp%zone_index))//'_jaE_dens.dat')
            call save_arr(qp%dimx, qp%x, qp%jaei_arr, &
                          trim(qp%path2linear)//'zone_'//trim(itoa(qp%zone_index))//'_jaE_int.dat')
        end if
    end subroutine calculate_jae_distributed

    !==================================================================
    ! handle plumbing
    !==================================================================

    subroutine handle_to_qp(handle, qp)
        integer(c_intptr_t), value :: handle
        type(flre_quants_t), pointer, intent(out) :: qp
        type(c_ptr) :: ccp
        ccp = transfer(handle, ccp)
        call c_f_pointer(ccp, qp)
    end subroutine handle_to_qp

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

end module kilca_flre_quants_m
