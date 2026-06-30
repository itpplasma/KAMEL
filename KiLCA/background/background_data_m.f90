!> Background equilibrium profiles and physics, formerly the C++ background
!> class (background.{h,cpp} + calc_back.cpp + eval_back.cpp). Translated as
!> a Fortran SINGLETON module, not a per-instance handle: the only live
!> `new background(...)` call site anywhere in the tree is core.cpp's
!> calc_and_set_mode_independent_core_data, called exactly once per run (the
!> other `new background` hit, post_processing/torque_facs.cpp, is dead -
!> not in any CMake target, like transf_quants.cpp/quants_profs.h which it
!> belongs to).
!>
!> bp_ptr (core_m.f90) is read by ~10 existing live Fortran files
!> (conduct_parameters.f90's eval_and_set_background_parameters_spec_
!> independent/_spec_dependent/eval_and_set_f0_parameters_nu_and_derivs, and
!> transitively gcorr.f90/kmatrices*.f90/diff_sys.f90) as an opaque handle
!> passed to the (formerly C++) eval_background_spec_independent_/
!> eval_f0_parameters_nu_and_derivs_/vs_0_f_/eval_hthz_ extern "C" functions,
!> which dereferenced it as `background*`. Since background is now a
!> singleton with no heap address to hand out, bp_ptr is kept only as a
!> harmless nonzero sentinel (set once by background_create_) so none of
!> those ~10 existing Fortran files need editing; the bind(C) functions
!> below still accept a handle argument to match those callers' existing
!> argument lists exactly, but ignore its value and operate on this
!> module's own state.
module kilca_background_data_m
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use constants, only: dp, pi, boltz, c, mp, me, e
    use kilca_spline_m, only: spline_alloc, spline_calc, spline_eval, spline_free
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_ode_rk8pd, only: rk8pd_state_t, rk8pd_evolve_init, rk8pd_evolve_apply
    use background, only: flag_back, calc_back, N, path2profiles, &
        rtor, B0, V_gal_sys, V_scale, mass, charge, zion, zele, flag_debug
    implicit none
    private

    public :: background_create, background_set_profiles_from_files
    public :: background_set_profiles_from_interface
    public :: background_interp_in_lab_frame
    public :: eval_background_spec_independent, eval_f0_parameters_nu_and_derivs
    public :: vs_0_f
    public :: eval_hthz_c, eval_hthz
    public :: eval_Bt_dBt_Bz_dBz, eval_p_dp, eval_mass_density, eval_Bt_Bz
    public :: eval_Bt, eval_Bz, eval_Vt, eval_Vz, eval_B0_ht_hz_n0_Vz
    public :: get_background_x0, get_background_xlast, get_background_dimx

    integer(c_int) :: Nprofiles
    character(len=8), allocatable :: profile_names(:)
    character(len=1024) :: path2background

    integer(c_int) :: dimx
    real(dp), allocatable, target :: x(:)

    integer(c_int) :: dimy
    real(dp), allocatable, target :: y(:)

    integer(c_int) :: ind_search
    real(dp), allocatable, target :: Carr(:)
    real(dp), allocatable, target :: Rarr(:)

    integer(c_intptr_t) :: sid = 0

    integer(c_int) :: i_q, i_n
    integer(c_int) :: i_T(0:1)
    integer(c_int) :: i_Vth(0:2), i_Vz(0:2)
    integer(c_int) :: i_Er
    integer(c_int) :: i_hth, i_hz, i_Bth, i_Bz, i_B
    integer(c_int) :: i_J0th(0:2), i_J0z(0:2)
    integer(c_int) :: i_dPhi0
    integer(c_int) :: i_n_p(0:1), i_Vp_p(0:1), i_Vt_p(0:1), i_nu(0:1)

    integer(c_int) :: flag_dPhi0_calc

    interface
        integer(c_int) function get_output_flag_background() bind(C, name="get_output_flag_background_")
            import :: c_int
        end function get_output_flag_background

        function eval_path_to_linear_data_c(buf) result(n) bind(C)
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: buf(*)
            integer(c_int) :: n
        end function eval_path_to_linear_data_c

        function count_lines_in_file_c(filename) result(n) bind(C, name="count_lines_in_file_")
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: filename(*)
            integer(c_int) :: n
        end function count_lines_in_file_c

        function load_data_file_c(file_name, dim_, ncols, rgrid, qgrid) result(ierr) &
            bind(C, name="load_data_file")
            import :: c_int, c_char, c_double
            character(kind=c_char), intent(in) :: file_name(*)
            integer(c_int), value :: dim_, ncols
            real(c_double), intent(in) :: rgrid(*)
            real(c_double), intent(out) :: qgrid(*)
            integer(c_int) :: ierr
        end function load_data_file_c

        integer(c_int) function save_real_array_c(dim_, xgrid, arr, full_name) &
            bind(C, name="save_real_array")
            import :: c_int, c_double, c_char
            integer(c_int), value :: dim_
            real(c_double), intent(in) :: xgrid(*), arr(*)
            character(kind=c_char), intent(in) :: full_name(*)
        end function save_real_array_c

        function get_background_dimension_from_balance_c(dimx_) result(stat) &
            bind(C, name="get_background_dimension_from_balance_")
            import :: c_int
            integer(c_int), intent(out) :: dimx_
            integer(c_int) :: stat
        end function get_background_dimension_from_balance_c
    end interface

    interface
        subroutine get_background_profiles_from_balance_x(dimx_, x_, q_, n_, Ti_, Te_, Vth_, Vz_, Er_) &
            bind(C, name="get_background_profiles_from_balance_")
            import :: c_int, c_double
            integer(c_int), intent(in) :: dimx_
            real(c_double), intent(out) :: x_(*), q_(*), n_(*), Ti_(*), Te_(*), Vth_(*), Vz_(*), Er_(*)
        end subroutine get_background_profiles_from_balance_x
    end interface

contains

    integer(c_int) function count_lines(fname) result(n)
        character(len=*), intent(in) :: fname
        character(kind=c_char) :: cbuf(len(fname) + 1)
        call to_cstr(fname, cbuf)
        n = count_lines_in_file_c(cbuf)
    end function count_lines

    subroutine load_data(fname, dim_, ncols, rgrid, qgrid)
        character(len=*), intent(in) :: fname
        integer(c_int), intent(in) :: dim_, ncols
        real(c_double), intent(in) :: rgrid(*)
        real(c_double), intent(out) :: qgrid(*)
        character(kind=c_char) :: cbuf(len(fname) + 1)
        integer(c_int) :: ierr_unused
        call to_cstr(fname, cbuf)
        ierr_unused = load_data_file_c(cbuf, dim_, ncols, rgrid, qgrid)
    end subroutine load_data

    subroutine save_arr(dim_, xgrid, arr, fname)
        integer(c_int), intent(in) :: dim_
        real(c_double), intent(in) :: xgrid(*), arr(*)
        character(len=*), intent(in) :: fname
        character(kind=c_char) :: cbuf(len(fname) + 1)
        integer :: ierr_unused
        call to_cstr(fname, cbuf)
        ierr_unused = save_real_array_c(dim_, xgrid, arr, cbuf)
    end subroutine save_arr

    subroutine to_cstr(str, cbuf)
        character(len=*), intent(in) :: str
        character(kind=c_char), intent(out) :: cbuf(*)
        integer :: i, n
        n = len_trim(str)
        do i = 1, n
            cbuf(i) = str(i:i)
        end do
        cbuf(n + 1) = c_null_char
    end subroutine to_cstr

    !> Mirrors background::background(const settings *s): sets up the index
    !> table and the output directory. path2project comes from the caller
    !> (core.cpp, still C++, has sd->path2project directly). Returns a fixed
    !> nonzero sentinel handle for compatibility with core.h's intptr_t bp
    !> field and the existing bp_ptr plumbing - the value itself is never
    !> read back by anything in this module.
    function background_create(path2project_p) result(handle) bind(C, name="background_create_")
        type(c_ptr), value :: path2project_p
        integer(c_intptr_t) :: handle
        character(len=1024) :: path2project
        integer(c_int) :: k, ierr_unused
        character(len=1024) :: sys_command
        external :: set_back_aliases_in_conductivity_parameters

        call read_c_string(path2project_p, path2project)

        call set_profiles_indices()

        call set_back_aliases_in_conductivity_parameters()

        path2background = trim(path2project)//'background-data/'

        if (get_output_flag_background() > 1) then
            sys_command = 'mkdir -p '//trim(path2background)
            call execute_command_line(trim(sys_command), exitstat=ierr_unused)
        end if

        handle = 1_c_intptr_t
    end function background_create

    !> Mirrors background::set_profiles_indices exactly: the k++ assignment
    !> order is load-bearing (the oracle's own comment: "DO NOT CHANGE, ONLY
    !> ADD NEW IF NEEDED") since every spline channel index downstream
    !> depends on it.
    subroutine set_profiles_indices()
        integer(c_int) :: k, acalc

        acalc = abs(calc_back)

        if (acalc == 1 .or. acalc == 3) then
            Nprofiles = 7
        else if (acalc == 2) then
            Nprofiles = 11
        else
            write (error_unit, '(a)') 'warning: set_profiles_indices: unknown flag!'
            stop 1
        end if

        allocate (profile_names(0:Nprofiles - 1))

        if (acalc == 1 .or. acalc == 3) then
            profile_names(0) = 'q'; profile_names(1) = 'n'; profile_names(2) = 'Ti'
            profile_names(3) = 'Te'; profile_names(4) = 'Vth'; profile_names(5) = 'Vz'
            profile_names(6) = 'Er'
        else if (acalc == 2) then
            profile_names(0) = 'q'; profile_names(1) = 'n'; profile_names(2) = 'Ti'
            profile_names(3) = 'Te'; profile_names(4) = 'Vth'; profile_names(5) = 'Vz'
            profile_names(6) = 'Er'; profile_names(7) = 'Bth'; profile_names(8) = 'Bz'
            profile_names(9) = 'Jth'; profile_names(10) = 'Jz'
        end if

        k = 0
        i_q = k; k = k + 1
        i_n = k; k = k + 1
        i_T(0) = k; k = k + 1
        i_T(1) = k; k = k + 1
        i_Vth(2) = k; k = k + 1
        i_Vz(2) = k; k = k + 1
        i_Er = k; k = k + 1
        i_Bth = k; k = k + 1
        i_Bz = k; k = k + 1
        i_B = k; k = k + 1
        i_hth = k; k = k + 1
        i_hz = k; k = k + 1
        i_dPhi0 = k; k = k + 1
        i_n_p(0) = k; k = k + 1
        i_Vp_p(0) = k; k = k + 1
        i_Vt_p(0) = k; k = k + 1
        i_nu(0) = k; k = k + 1
        i_n_p(1) = k; k = k + 1
        i_Vp_p(1) = k; k = k + 1
        i_Vt_p(1) = k; k = k + 1
        i_nu(1) = k; k = k + 1
        i_J0th(2) = k; k = k + 1
        i_J0z(2) = k; k = k + 1
        i_J0th(1) = k; k = k + 1
        i_J0z(1) = k; k = k + 1
        i_J0th(0) = k; k = k + 1
        i_J0z(0) = k; k = k + 1
        i_Vth(1) = k; k = k + 1
        i_Vz(1) = k; k = k + 1
        i_Vth(0) = k; k = k + 1
        i_Vz(0) = k; k = k + 1

        dimy = k
    end subroutine set_profiles_indices

    subroutine background_set_profiles_from_files() bind(C, name="background_set_profiles_from_files_")
        character(len=1024) :: file_name, path2profiles_buf
        integer(c_int) :: Nspl, ierr, acalc

        call get_path2profiles(path2profiles_buf)
        file_name = trim(path2profiles_buf)//'n.dat'

        dimx = count_lines(file_name)

        Nspl = N
        ind_search = int(0.5d0*dimx)

        allocate (x(0:dimx - 1))
        allocate (y(0:dimx*dimy - 1))
        allocate (Carr(0:(Nspl + 1)*dimx*dimy - 1))
        allocate (Rarr(0:(Nspl + 1)*dimy - 1))

        call spline_alloc(Nspl, 1, dimx, c_loc(x), c_loc(Carr), sid)

        call load_input_background_profiles()
        call spline_input_profiles()

        acalc = abs(calc_back)
        if (acalc == 1) then
            call calculate_equilibrium()
        else if (acalc == 2) then
            call check_and_spline_equilibrium()
        else
            write (error_unit, '(a)') &
                'set_background_profiles_from_files: this feature is not implemented yet!'
            stop 1
        end if

        flag_dPhi0_calc = 1

        call find_f0_parameters()

        if (get_output_flag_background() > 1) then
            call save_background()
            call eval_and_save_f0_moments()
        end if
    end subroutine background_set_profiles_from_files

    subroutine background_set_profiles_from_interface() &
        bind(C, name="background_set_profiles_from_interface_")
        integer(c_int) :: Nspl, ierr, acalc, stat_unused
        real(dp), allocatable :: q_(:), n_(:), Ti_(:), Te_(:), Vth_(:), Vz_(:), Er_(:)
        integer(c_int) :: i

        stat_unused = get_background_dimension_from_balance_c(dimx)

        ind_search = int(0.5d0*dimx)
        Nspl = N

        allocate (x(0:dimx - 1))
        allocate (y(0:dimx*dimy - 1))
        allocate (Carr(0:(Nspl + 1)*dimx*dimy - 1))
        allocate (Rarr(0:(Nspl + 1)*dimy - 1))

        call spline_alloc(Nspl, 1, dimx, c_loc(x), c_loc(Carr), sid)

        allocate (q_(0:dimx - 1), n_(0:dimx - 1), Ti_(0:dimx - 1), Te_(0:dimx - 1))
        allocate (Vth_(0:dimx - 1), Vz_(0:dimx - 1), Er_(0:dimx - 1))

        call get_background_profiles_from_balance_x(dimx, x, q_, n_, Ti_, Te_, Vth_, Vz_, Er_)

        do i = 0, dimx - 1
            y(i_q*dimx + i) = q_(i)
            y(i_n*dimx + i) = n_(i)
            y(i_T(0)*dimx + i) = Ti_(i)
            y(i_T(1)*dimx + i) = Te_(i)
            y(i_Vth(2)*dimx + i) = Vth_(i)
            Vz_(i) = Vz_(i) - V_gal_sys
            y(i_Vz(2)*dimx + i) = Vz_(i)
            y(i_Er*dimx + i) = Er_(i)
        end do

        deallocate (q_, n_, Ti_, Te_, Vth_, Vz_, Er_)

        call spline_input_profiles()

        flag_dPhi0_calc = 1

        acalc = abs(calc_back)
        if (acalc == 1) then
            call calculate_equilibrium()
        else if (acalc == 2) then
            write (error_unit, '(a)') &
                'set_background_profiles_from_interface: this feature is not implemented yet!'
            stop 1
        else if (acalc == 3) then
            flag_dPhi0_calc = 1
            call calculate_equilibrium()
        else
            write (error_unit, '(a)') &
                'set_background_profiles_from_files: this feature is not implemented yet!'
            stop 1
        end if

        call find_f0_parameters()

        if (get_output_flag_background() > 1) then
            call save_background()
            call eval_and_save_f0_moments()
        end if
    end subroutine background_set_profiles_from_interface

    subroutine get_path2profiles(buf)
        character(len=1024), intent(out) :: buf
        interface
            subroutine get_background_path2profiles_c(out) bind(C, name="get_background_path2profiles_")
                import :: c_char
                character(kind=c_char), intent(out) :: out(*)
            end subroutine get_background_path2profiles_c
        end interface
        character(kind=c_char) :: cbuf(1024)
        integer :: i

        call get_background_path2profiles_c(cbuf)
        buf = ''
        do i = 1, 1024
            if (cbuf(i) == c_null_char) exit
            buf(i:i) = cbuf(i)
        end do
    end subroutine get_path2profiles

    !> Mirrors background::load_input_background_profiles: all profiles are
    !> assumed to be in the lab frame; Vz is shifted to the moving frame here.
    subroutine load_input_background_profiles()
        character(len=1024) :: file_name, sys_command, path2profiles_buf
        integer(c_int) :: ind(0:10)
        integer(c_int) :: i, j
        real(dp) :: tmp1, tmp2

        call get_path2profiles(path2profiles_buf)

        ind = [i_q, i_n, i_T(0), i_T(1), i_Vth(2), i_Vz(2), i_Er, i_Bth, i_Bz, i_J0th(2), i_J0z(2)]

        tmp1 = 0.0d0
        do i = 0, Nprofiles - 1
            file_name = trim(path2profiles_buf)//trim(profile_names(i))//'.dat'
            call load_data(file_name, dimx, 1, x, y(ind(i)*dimx:))

            if (i == 5) then
                do j = 0, dimx - 1
                    y(ind(i)*dimx + j) = y(ind(i)*dimx + j)*V_scale
                    y(ind(i)*dimx + j) = y(ind(i)*dimx + j) - V_gal_sys
                end do
            end if

            tmp2 = 0.0d0
            do j = 0, dimx - 1
                tmp2 = tmp2 + x(j)
            end do
            if (tmp2 /= tmp1 .and. i /= 0) then
                write (output_unit, '(a,es25.16,a,es25.16,a,i0,a)') &
                    'warning: load_input_profiles: r grids look different: tmp1=', tmp1, &
                    ' tmp2=', tmp2, ' i=', i, '!'
            end if
            tmp1 = tmp2

            if (get_output_flag_background() > 1) then
                sys_command = 'cp '//trim(file_name)//' '//trim(path2background)// &
                              trim(profile_names(i))//'_i.dat'
                call execute_command_line(trim(sys_command))
            end if
        end do

        if (x(0) == 0.0d0) then
            write (error_unit, '(a)') &
                'warning: load_input_profiles: first point of profiles grid must be not zero!'
        end if
    end subroutine load_input_background_profiles

    !> Splines the first 7 profiles: q, n, Ti, Te, Vth, Vz (moving frame), Er
    !> (lab frame).
    subroutine spline_input_profiles()
        integer(c_int) :: ierr
        call spline_calc(sid, c_loc(y(i_q*dimx:)), i_q, i_Er, c_null_ptr, ierr)
    end subroutine spline_input_profiles

    !> ODE right-hand side for the equilibrium u-function (calc_back.cpp's
    !> rhs_back): no ctx needed, background is a singleton.
    subroutine rhs_back(t, yv, dydt, ctx)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: yv(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        real(dp), target :: rloc(1)
        real(dp) :: qq, n_, dn, Ti, dTi, Te, dTe, dpress, g

        rloc(1) = t
        call spline_eval(sid, 1, c_loc(rloc), 0, 1, i_q, i_T(1), c_loc(Rarr))

        qq = Rarr(0)
        n_ = Rarr(2)
        dn = Rarr(3)
        Ti = Rarr(4)
        dTi = Rarr(5)
        Te = Rarr(6)
        dTe = Rarr(7)

        dpress = boltz*((Ti + Te)*dn + (dTi + dTe)*n_)

        g = 1.0d0 + t*t/rtor/rtor/qq/qq

        dydt(1) = -2.0d0*t*yv(1)/(qq*qq*g*rtor*rtor) - 8.0d0*pi*dpress
    end subroutine rhs_back

    !> Mirrors background::calculate_equilibrium exactly, including the
    !> single continuously-carried rk8pd evolve across the whole grid (per
    !> the oracle's own comment, splitting this into independent per-segment
    !> integrations or changing the tolerances drifts from the recorded
    !> golden - do not "clean up" these magic numbers).
    subroutine calculate_equilibrium()
        integer(c_int), parameter :: Neq = 1
        real(dp) :: rc, g, rcur
        real(dp), allocatable :: u(:)
        real(dp) :: uval(1)
        real(dp), target :: rloc(1)
        type(rk8pd_state_t) :: ode_state
        type(fortnum_status_t) :: status
        integer(c_int) :: nfev, i, ierr
        real(dp), pointer :: qarr(:), bth(:), bz(:), b0arr(:), hth(:), hz(:), jth(:), jz(:)
        real(dp) :: dbz, dbth

        rc = x(0)

        rloc(1) = rc
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_q, i_q, c_loc(Rarr))

        g = 1.0d0 + rc*rc/rtor/rtor/Rarr(0)/Rarr(0)

        allocate (u(0:dimx - 1))
        u(0) = B0*B0*g
        uval(1) = u(0)

        call rk8pd_evolve_init(ode_state, Neq, x(1) - rc, status)
        if (status%code /= FORTNUM_OK) then
            write (error_unit, '(a)') 'calculate_equilibrium: failed to create ODE evolver'
            stop 1
        end if

        rcur = rc
        nfev = 0

        do i = 1, dimx - 1
            call rk8pd_evolve_apply(rhs_back, ode_state, rcur, x(i), uval, &
                                    1.0d-16, 1.0d-16, 100000, nfev, status)

            if (status%code /= FORTNUM_OK) then
                write (error_unit, '(a,es12.4,a,i0)') &
                    'calculate_equilibrium: ODE solver failed at r = ', rcur, ' (status=', status%code
                stop 1
            end if

            u(i) = uval(1)
        end do

        if (flag_debug > 1) then
            call save_arr(dimx, x, u, trim(path2background)//'u.dat')
        end if

        qarr(0:dimx - 1) => y(i_q*dimx:i_q*dimx + dimx - 1)
        bth(0:dimx - 1) => y(i_Bth*dimx:i_Bth*dimx + dimx - 1)
        bz(0:dimx - 1) => y(i_Bz*dimx:i_Bz*dimx + dimx - 1)
        b0arr(0:dimx - 1) => y(i_B*dimx:i_B*dimx + dimx - 1)

        do i = 0, dimx - 1
            bz(i) = sign(1.0d0, B0)*sqrt(u(i)/(1.0d0 + x(i)*x(i)/rtor/rtor/qarr(i)/qarr(i)))
            bth(i) = bz(i)*x(i)/qarr(i)/rtor
            b0arr(i) = sqrt(bth(i)*bth(i) + bz(i)*bz(i))
        end do

        hth(0:dimx - 1) => y(i_hth*dimx:i_hth*dimx + dimx - 1)
        hz(0:dimx - 1) => y(i_hz*dimx:i_hz*dimx + dimx - 1)

        do i = 0, dimx - 1
            hth(i) = bth(i)/b0arr(i)
            hz(i) = bz(i)/b0arr(i)
        end do

        call spline_calc(sid, c_loc(y(i_Bth*dimx:)), i_Bth, i_hz, c_null_ptr, ierr)

        jth(0:dimx - 1) => y(i_J0th(2)*dimx:i_J0th(2)*dimx + dimx - 1)
        jz(0:dimx - 1) => y(i_J0z(2)*dimx:i_J0z(2)*dimx + dimx - 1)

        do i = 0, dimx - 1
            rloc(1) = x(i)
            call spline_eval(sid, 1, c_loc(rloc), 1, 1, i_Bth, i_Bz, c_loc(Rarr))
            dbth = Rarr(0)
            dbz = Rarr(1)

            jth(i) = -c/4.0d0/pi*dbz
            jz(i) = c/4.0d0/pi*(bth(i)/x(i) + dbth)
        end do

        call spline_calc(sid, c_loc(y(i_J0th(2)*dimx:)), i_J0th(2), i_J0z(2), c_null_ptr, ierr)

        deallocate (u)
    end subroutine calculate_equilibrium

    !> Mirrors background::check_and_spline_equilibrium exactly.
    subroutine check_and_spline_equilibrium()
        integer(c_int) :: i, ierr
        real(dp), pointer :: bth(:), bz(:), b0arr(:), hth(:), hz(:), jth(:)
        real(dp) :: dev, Bt, dBt, Bz_, dBz_, p, dpp
        real(dp), target :: rloc(1)

        bth(0:dimx - 1) => y(i_Bth*dimx:i_Bth*dimx + dimx - 1)
        bz(0:dimx - 1) => y(i_Bz*dimx:i_Bz*dimx + dimx - 1)
        b0arr(0:dimx - 1) => y(i_B*dimx:i_B*dimx + dimx - 1)
        hth(0:dimx - 1) => y(i_hth*dimx:i_hth*dimx + dimx - 1)
        hz(0:dimx - 1) => y(i_hz*dimx:i_hz*dimx + dimx - 1)

        do i = 0, dimx - 1
            b0arr(i) = sqrt(bth(i)*bth(i) + bz(i)*bz(i))
            hth(i) = bth(i)/b0arr(i)
            hz(i) = bz(i)/b0arr(i)
        end do

        call spline_calc(sid, c_loc(y(i_Bth*dimx:)), i_Bth, i_hz, c_null_ptr, ierr)
        call spline_calc(sid, c_loc(y(i_J0th(2)*dimx:)), i_J0th(2), i_J0z(2), c_null_ptr, ierr)

        do i = 0, dimx - 1
            rloc(1) = x(i)
            call spline_eval(sid, 1, c_loc(rloc), 0, 1, i_Bth, i_Bth, c_loc(Rarr))
            Bt = Rarr(0); dBt = Rarr(1)
            call spline_eval(sid, 1, c_loc(rloc), 0, 1, i_Bz, i_Bz, c_loc(Rarr))
            Bz_ = Rarr(0); dBz_ = Rarr(1)

            call eval_p_dp_local(x(i), Rarr)
            p = Rarr(0); dpp = Rarr(1)

            dev = (4.0d0*pi*dpp + Bt*dBt + Bz_*dBz_)/(Bt*Bt/x(i)) + 1.0d0

            if (abs(dev) > 1.0d-6) then
                write (error_unit, '(a,es12.4,a,es12.4)') &
                    'warning: check_and_spline_equilibrium: r = ', x(i), ' dev = ', dev
            end if
        end do
    end subroutine check_and_spline_equilibrium

    !> Mirrors background::find_f0_parameters exactly.
    subroutine find_f0_parameters()
        integer(c_int) :: i, spec, ierr
        real(dp), pointer :: n_(:), Ti_(:), Te_(:), Vth_(:), Vz_(:), Er_(:)
        real(dp), pointer :: b0_(:), bth_(:), hth_(:), hz_(:), jth_(:), jz_(:), dPhi0(:)
        real(dp), pointer :: nui(:), nue(:), n_i_p(:), n_e_p(:)
        real(dp), pointer :: Vp_i_p(:), Vp_e_p(:), Vt_i_p(:), Vt_e_p(:)
        real(dp) :: Vs_tot, Vp_tot, Js_tot, Jp_tot
        real(dp) :: dpress, dn, dTi, dTe, Vsi_, Vse_
        real(dp) :: nuee, nuei, nuie, nuii, Lee, Lei, Lie, Lii
        real(dp), parameter :: vf = 1.0d0
        real(dp), target :: rloc(1)
        external :: eval_and_set_background_parameters_spec_independent
        external :: eval_and_set_background_parameters_spec_dependent
        external :: eval_and_set_f0_parameters_nu_and_derivs
        external :: dens_par

        n_(0:dimx - 1) => y(i_n*dimx:i_n*dimx + dimx - 1)
        Ti_(0:dimx - 1) => y(i_T(0)*dimx:i_T(0)*dimx + dimx - 1)
        Te_(0:dimx - 1) => y(i_T(1)*dimx:i_T(1)*dimx + dimx - 1)
        Vth_(0:dimx - 1) => y(i_Vth(2)*dimx:i_Vth(2)*dimx + dimx - 1)
        Vz_(0:dimx - 1) => y(i_Vz(2)*dimx:i_Vz(2)*dimx + dimx - 1)
        Er_(0:dimx - 1) => y(i_Er*dimx:i_Er*dimx + dimx - 1)

        b0_(0:dimx - 1) => y(i_B*dimx:i_B*dimx + dimx - 1)
        bth_(0:dimx - 1) => y(i_Bth*dimx:i_Bth*dimx + dimx - 1)
        hth_(0:dimx - 1) => y(i_hth*dimx:i_hth*dimx + dimx - 1)
        hz_(0:dimx - 1) => y(i_hz*dimx:i_hz*dimx + dimx - 1)

        jth_(0:dimx - 1) => y(i_J0th(2)*dimx:i_J0th(2)*dimx + dimx - 1)
        jz_(0:dimx - 1) => y(i_J0z(2)*dimx:i_J0z(2)*dimx + dimx - 1)

        dPhi0(0:dimx - 1) => y(i_dPhi0*dimx:i_dPhi0*dimx + dimx - 1)

        nui(0:dimx - 1) => y(i_nu(0)*dimx:i_nu(0)*dimx + dimx - 1)
        nue(0:dimx - 1) => y(i_nu(1)*dimx:i_nu(1)*dimx + dimx - 1)

        n_i_p(0:dimx - 1) => y(i_n_p(0)*dimx:i_n_p(0)*dimx + dimx - 1)
        n_e_p(0:dimx - 1) => y(i_n_p(1)*dimx:i_n_p(1)*dimx + dimx - 1)

        Vp_i_p(0:dimx - 1) => y(i_Vp_p(0)*dimx:i_Vp_p(0)*dimx + dimx - 1)
        Vp_e_p(0:dimx - 1) => y(i_Vp_p(1)*dimx:i_Vp_p(1)*dimx + dimx - 1)

        Vt_i_p(0:dimx - 1) => y(i_Vt_p(0)*dimx:i_Vt_p(0)*dimx + dimx - 1)
        Vt_e_p(0:dimx - 1) => y(i_Vt_p(1)*dimx:i_Vt_p(1)*dimx + dimx - 1)

        do i = 0, dimx - 1
            Lee = 23.5d0 - log(sqrt(n_(i))/Te_(i)**1.25d0) - (log(Te_(i)) - 2.0d0)/4.0d0
            Lei = 30.0d0 - log(sqrt(n_(i))/Ti_(i)**1.5d0*(charge(0)/e)**2/(mass(0)/mp))
            Lie = Lei
            Lii = 23.0d0 - log((charge(0)/e)**3/Ti_(i)**1.5d0*sqrt(2.0d0*n_(i)))

            nuee = (5.8d-6)*n_(i)*Lee/sqrt(Te_(i))/(vf*Te_(i))
            nuei = (7.7d-6)*n_(i)*Lei*(charge(0)/e)**2/(vf*Te_(i))**1.5d0
            nuie = (3.2d-9)*n_(i)*Lie*(charge(0)/e)**2/(mass(0)/mp)/sqrt(Te_(i))/(vf*Ti_(i))
            nuii = (1.4d-7)*n_(i)*Lii*(charge(0)/e)**4/sqrt(mass(0)/mp)/sqrt(Ti_(i))/(vf*Ti_(i))

            nui(i) = zion*(nuie + nuii + 10.0d0)
            nue(i) = zele*(nuee + nuei + 10.0d0)

            n_i_p(i) = n_(i)
            n_e_p(i) = n_(i)

            Vt_i_p(i) = sqrt(boltz*Ti_(i)/mass(0))
            Vt_e_p(i) = sqrt(boltz*Te_(i)/mass(1))

            Vp_tot = hth_(i)*Vth_(i) + hz_(i)*Vz_(i)
            Jp_tot = hth_(i)*jth_(i) + hz_(i)*jz_(i)

            Vp_i_p(i) = Vp_tot + Jp_tot*mass(1)/mass(0)/e/n_(i)/(1.0d0 + mass(1)/mass(0))
            Vp_e_p(i) = Vp_tot - Jp_tot/e/n_(i)/(1.0d0 + mass(1)/mass(0))

            if (flag_dPhi0_calc > 0) then
                Vs_tot = hz_(i)*Vth_(i) - hth_(i)*Vz_(i)
                Js_tot = hz_(i)*jth_(i) - hth_(i)*jz_(i)

                Vsi_ = Vs_tot + Js_tot*mass(1)/mass(0)/e/n_(i)/(1.0d0 + mass(1)/mass(0))
                Vse_ = Vs_tot - Js_tot/e/n_(i)/(1.0d0 + mass(1)/mass(0))

                rloc(1) = x(i)
                call spline_eval(sid, 1, c_loc(rloc), 1, 1, i_n, i_T(1), c_loc(Rarr))
                dn = Rarr(0); dTi = Rarr(1); dTe = Rarr(2)

                dpress = boltz*(Ti_(i)*dn + dTi*n_(i))
                dPhi0(i) = Vsi_/c*b0_(i) - dpress/charge(0)/n_(i)
            else
                dPhi0(i) = -Er_(i) + bth_(i)*V_gal_sys/c
            end if
        end do

        call spline_calc(sid, c_loc(y(i_dPhi0*dimx:)), i_dPhi0, i_nu(1), c_null_ptr, ierr)

        do i = 0, dimx - 1
            call eval_and_set_background_parameters_spec_independent(x(i), flag_back)

            do spec = 0, 1
                call eval_and_set_background_parameters_spec_dependent(x(i), spec, flag_back)
                call eval_and_set_f0_parameters_nu_and_derivs(x(i), spec, flag_back)
                call dens_par(n_(i), y(i_n_p(spec)*dimx + i))
            end do
        end do

        call spline_calc(sid, c_loc(y(i_n_p(0)*dimx:)), i_n_p(0), i_n_p(0), c_null_ptr, ierr)
        call spline_calc(sid, c_loc(y(i_n_p(1)*dimx:)), i_n_p(1), i_n_p(1), c_null_ptr, ierr)
    end subroutine find_f0_parameters

    !> Mirrors background::eval_and_save_f0_moments (debug dump, gated
    !> behind the RUNTIME flag get_output_flag_background_()>1, not a
    !> compile-time DEBUG_FLAG - keep computing it).
    subroutine eval_and_save_f0_moments()
        integer(c_int), parameter :: num_moms = 16
        integer(c_int) :: u_(0:num_moms - 1)
        character(len=8) :: moms_name(0:num_moms - 1)
        character(len=1024) :: full_name
        integer(c_int) :: i, k, spec
        real(dp) :: n_m(0:1), Vth_m(0:2), Vz_m(0:2), Vs_m(0:2), Vs_0_m(0:2), Vp_m(0:2), T_m(0:1)
        real(dp) :: j0th_m(0:2), j0z_m(0:2), j0s_m(0:2), j0p_m(0:2)
        external :: eval_and_set_background_parameters_spec_independent
        external :: eval_and_set_background_parameters_spec_dependent
        external :: eval_and_set_f0_parameters_nu_and_derivs
        external :: dens_mom, vth_mom, vz_mom, vs_mom, vp_mom, eterm_mom

        moms_name = ['ni_m    ', 'ne_m    ', 'Vth_m   ', 'Vz_m    ', &
                     'Vsi_m   ', 'Vse_m   ', 'Vpi_m   ', 'Vpe_m   ', &
                     'Vs0i_m  ', 'Vs0e_m  ', 'Ti_m    ', 'Te_m    ', &
                     'j0s_m   ', 'j0p_m   ', 'j0th_m  ', 'j0z_m   ']

        do k = 0, num_moms - 1
            full_name = trim(path2background)//trim(moms_name(k))//'.dat'
            open (newunit=u_(k), file=trim(full_name), status='replace', action='write')
        end do

        do i = 0, dimx - 1
            call eval_and_set_background_parameters_spec_independent(x(i), flag_back)

            do spec = 0, 1
                call eval_and_set_background_parameters_spec_dependent(x(i), spec, flag_back)
                call eval_and_set_f0_parameters_nu_and_derivs(x(i), spec, flag_back)

                call dens_mom(n_m(spec))
                call vth_mom(Vth_m(spec))
                Vth_m(spec) = Vth_m(spec)*x(i)
                call vz_mom(Vz_m(spec))
                call vs_mom(Vs_m(spec))
                call vs_0(x(i), spec, Vs_0_m(spec))
                call vp_mom(Vp_m(spec))
                call eterm_mom(T_m(spec))
                T_m(spec) = 2.0d0*T_m(spec)/3.0d0/n_m(spec)
            end do

            Vth_m(2) = (mass(0)*Vth_m(0) + mass(1)*Vth_m(1))/(mass(0) + mass(1))
            Vz_m(2) = (mass(0)*Vz_m(0) + mass(1)*Vz_m(1))/(mass(0) + mass(1))

            j0p_m(2) = e*(Vp_m(0) - Vp_m(1))
            j0s_m(2) = e*(Vs_m(0) - Vs_m(1))
            j0th_m(2) = e*(Vth_m(0) - Vth_m(1))
            j0z_m(2) = e*(Vz_m(0) - Vz_m(1))

            write (u_(0), '(es24.15,1x,es24.15)') x(i), n_m(0)
            write (u_(1), '(es24.15,1x,es24.15)') x(i), n_m(1)
            write (u_(2), '(es24.15,1x,es24.15)') x(i), Vth_m(2)/n_m(0)
            write (u_(3), '(es24.15,1x,es24.15)') x(i), Vz_m(2)/n_m(0)
            write (u_(4), '(es24.15,1x,es24.15)') x(i), Vs_m(0)/n_m(0)
            write (u_(5), '(es24.15,1x,es24.15)') x(i), Vs_m(1)/n_m(1)
            write (u_(6), '(es24.15,1x,es24.15)') x(i), Vp_m(0)/n_m(0)
            write (u_(7), '(es24.15,1x,es24.15)') x(i), Vp_m(1)/n_m(1)
            write (u_(8), '(es24.15,1x,es24.15)') x(i), Vs_0_m(0)
            write (u_(9), '(es24.15,1x,es24.15)') x(i), Vs_0_m(1)
            write (u_(10), '(es24.15,1x,es24.15)') x(i), T_m(0)
            write (u_(11), '(es24.15,1x,es24.15)') x(i), T_m(1)
            write (u_(12), '(es24.15,1x,es24.15)') x(i), j0s_m(2)
            write (u_(13), '(es24.15,1x,es24.15)') x(i), j0p_m(2)
            write (u_(14), '(es24.15,1x,es24.15)') x(i), j0th_m(2)
            write (u_(15), '(es24.15,1x,es24.15)') x(i), j0z_m(2)
        end do

        do k = 0, num_moms - 1
            close (u_(k))
        end do
    end subroutine eval_and_save_f0_moments

    subroutine save_background()
        call save_arr(dimx, x, y(i_Bth*dimx:), trim(path2background)//'b0th.dat')
        call save_arr(dimx, x, y(i_Bz*dimx:), trim(path2background)//'b0z.dat')
        call save_arr(dimx, x, y(i_B*dimx:), trim(path2background)//'b0.dat')
        call save_arr(dimx, x, y(i_hth*dimx:), trim(path2background)//'hth.dat')
        call save_arr(dimx, x, y(i_hz*dimx:), trim(path2background)//'hz.dat')
        call save_arr(dimx, x, y(i_J0th(2)*dimx:), trim(path2background)//'j0th.dat')
        call save_arr(dimx, x, y(i_J0z(2)*dimx:), trim(path2background)//'j0z.dat')
        call save_arr(dimx, x, y(i_nu(0)*dimx:), trim(path2background)//'nui.dat')
        call save_arr(dimx, x, y(i_nu(1)*dimx:), trim(path2background)//'nue.dat')
        call save_arr(dimx, x, y(i_dPhi0*dimx:), trim(path2background)//'dPhi0.dat')
        call save_arr(dimx, x, y(i_n_p(0)*dimx:), trim(path2background)//'ni_p.dat')
        call save_arr(dimx, x, y(i_n_p(1)*dimx:), trim(path2background)//'ne_p.dat')
        call save_arr(dimx, x, y(i_Vp_p(0)*dimx:), trim(path2background)//'Vpi_p.dat')
        call save_arr(dimx, x, y(i_Vp_p(1)*dimx:), trim(path2background)//'Vpe_p.dat')
        call save_arr(dimx, x, y(i_Vt_p(0)*dimx:), trim(path2background)//'Vti_p.dat')
        call save_arr(dimx, x, y(i_Vt_p(1)*dimx:), trim(path2background)//'Vte_p.dat')
    end subroutine save_background

    subroutine interp_basic_background_profiles(r, qq, n_, Ti, Te, Vth, Vz, dPhi0v)
        real(dp), intent(in) :: r
        real(dp), intent(out) :: qq, n_, Ti, Te, Vth, Vz, dPhi0v
        real(dp), target :: rloc(1)
        rloc(1) = r
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_q, i_Vz(2), c_loc(Rarr))
        qq = Rarr(0); n_ = Rarr(1); Ti = Rarr(2); Te = Rarr(3); Vth = Rarr(4); Vz = Rarr(5)
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_dPhi0, i_dPhi0, c_loc(Rarr))
        dPhi0v = Rarr(0)
    end subroutine interp_basic_background_profiles

    subroutine transform_basic_background_profiles_to_lab_frame(r, Vz, dPhi0v)
        real(dp), intent(in) :: r
        real(dp), intent(inout) :: Vz, dPhi0v
        real(dp) :: Bth
        real(dp), target :: rloc(1)
        rloc(1) = r
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_Bth, i_Bth, c_loc(Rarr))
        Bth = Rarr(0)
        Vz = Vz + V_gal_sys
        dPhi0v = dPhi0v - Bth*V_gal_sys/c
    end subroutine transform_basic_background_profiles_to_lab_frame

    subroutine background_interp_in_lab_frame(dim_, rv, qq, n_, Ti, Te, Vth, Vz, dPhi0v) &
        bind(C, name="interp_basic_background_profiles_in_lab_frame_")
        integer(c_int), value :: dim_
        real(c_double), intent(in) :: rv(0:dim_ - 1)
        real(c_double), intent(out) :: qq(0:dim_ - 1), n_(0:dim_ - 1), Ti(0:dim_ - 1), Te(0:dim_ - 1)
        real(c_double), intent(inout) :: Vth(0:dim_ - 1), Vz(0:dim_ - 1), dPhi0v(0:dim_ - 1)
        integer(c_int) :: i
        do i = 0, dim_ - 1
            call interp_basic_background_profiles(rv(i), qq(i), n_(i), Ti(i), Te(i), Vth(i), Vz(i), dPhi0v(i))
            call transform_basic_background_profiles_to_lab_frame(rv(i), Vz(i), dPhi0v(i))
        end do
    end subroutine background_interp_in_lab_frame

    real(c_double) function get_background_x0() bind(C, name="get_background_x0_")
        get_background_x0 = x(0)
    end function get_background_x0

    real(c_double) function get_background_xlast() bind(C, name="get_background_xlast_")
        get_background_xlast = x(dimx - 1)
    end function get_background_xlast

    integer(c_int) function get_background_dimx() bind(C, name="get_background_dimx_")
        get_background_dimx = dimx
    end function get_background_dimx

    !==================================================================
    ! eval_back.cpp live functions.
    !
    ! Two distinct calling conventions, matching the oracle exactly:
    !  - eval_background_spec_independent_/eval_f0_parameters_nu_and_derivs_
    !    (trailing underscore) are called ONLY from the pre-existing,
    !    unchanged Fortran conduct_parameters.f90 via an implicit interface
    !    (F77-style), so every argument is BY REFERENCE - matching the
    !    oracle's `double *r`/`int *spec`/`background **bpro` pointer args.
    !    The handle argument is accepted for ABI compatibility (so
    !    conduct_parameters.f90 needs zero edits) but never read, since
    !    background is a singleton.
    !  - eval_hthz/eval_Bt_dBt_Bz_dBz/eval_p_dp/eval_mass_density/eval_Bt_Bz/
    !    eval_Bt/eval_Bz/eval_Vt/eval_Vz/eval_B0_ht_hz_n0_Vz are called from
    !    still-C++ files (imhd/compressible_flow.cpp, imhd/incompressible.cpp,
    !    imhd/imhd_zone.cpp, flre/flre_zone.cpp, mode/transforms.cpp) as
    !    ordinary C++ function calls, so r/bp are BY VALUE (matching the
    !    oracle's plain `double r, const background *bp` signatures) - bp
    !    is accepted (background.h now forward-declares `class background;`
    !    with no body, so these C++ files keep compiling unchanged, passing
    !    around an opaque, never-dereferenced pointer value) but never read.
    !  eval_hthz_ (trailing underscore) is the BY-REFERENCE sibling of
    !  eval_hthz, called only from this port's own flre_quants_m.f90.
    !==================================================================

    subroutine eval_background_spec_independent(rval, handle, Rout) &
        bind(C, name="eval_background_spec_independent_")
        real(c_double), intent(in) :: rval
        integer(c_intptr_t), intent(in) :: handle
        real(c_double), intent(out), target :: Rout(*)
        real(dp), target :: rloc(1)
        rloc(1) = rval
        call spline_eval(sid, 1, c_loc(rloc), 0, 2, i_B, i_dPhi0, c_loc(Rout))
    end subroutine eval_background_spec_independent

    subroutine eval_f0_parameters_nu_and_derivs(rval, spec, handle, Rout) &
        bind(C, name="eval_f0_parameters_nu_and_derivs_")
        real(c_double), intent(in) :: rval
        integer(c_int), intent(in) :: spec
        integer(c_intptr_t), intent(in) :: handle
        real(c_double), intent(out), target :: Rout(*)
        real(dp), target :: rloc(1)
        rloc(1) = rval
        call spline_eval(sid, 1, c_loc(rloc), 0, 2, i_n_p(spec), i_nu(spec), c_loc(Rout))
    end subroutine eval_f0_parameters_nu_and_derivs

    !> Internal helper, called from eval_and_save_f0_moments and (via the
    !> bind(C) sibling vs_0_f below) gcorr.f90's still-live
    !> eval_and_set_params_for_additional_current_ - a genuine miss in an
    !> earlier "0 callers" grep pass (which only checked C++ files; gcorr.f90
    !> is Fortran), caught only by the link step. Re-confirms: always check
    !> Fortran (.f90) callers too, not just C++ ones, before calling
    !> something dead.
    subroutine vs_0(r, spec, res)
        real(dp), intent(in) :: r
        integer(c_int), intent(in) :: spec
        real(dp), intent(out) :: res
        real(dp), target :: rloc(1)
        real(dp) :: n_, dn, T, dT, dpress, Bv, dPhi0v

        rloc(1) = r
        call spline_eval(sid, 1, c_loc(rloc), 0, 1, i_n, i_n, c_loc(Rarr))
        n_ = Rarr(0); dn = Rarr(1)
        call spline_eval(sid, 1, c_loc(rloc), 0, 1, i_T(spec), i_T(spec), c_loc(Rarr))
        T = Rarr(0); dT = Rarr(1)

        dpress = boltz*(dn*T + n_*dT)

        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_B, i_B, c_loc(Rarr))
        Bv = Rarr(0)
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_dPhi0, i_dPhi0, c_loc(Rarr))
        dPhi0v = Rarr(0)

        res = c/Bv*(dPhi0v + dpress/charge(spec)/n_)
    end subroutine vs_0

    !> bind(C) replacement for vs_0_f_, called via bp_ptr from gcorr.f90 -
    !> same by-reference convention as eval_background_spec_independent_.
    subroutine vs_0_f(rval, spec, handle, res) bind(C, name="vs_0_f_")
        real(c_double), intent(in) :: rval
        integer(c_int), intent(in) :: spec
        integer(c_intptr_t), intent(in) :: handle
        real(c_double), intent(out) :: res
        real(dp) :: res_local
        call vs_0(rval, spec, res_local)
        res = res_local
    end subroutine vs_0_f

    subroutine eval_hthz_c(rval, omin, omax, bp, hout) bind(C, name="eval_hthz_")
        real(c_double), intent(in) :: rval
        integer(c_int), intent(in) :: omin, omax
        integer(c_intptr_t), intent(in) :: bp
        real(c_double), intent(out), target :: hout(*)
        real(dp), target :: rloc(1)
        rloc(1) = rval
        call spline_eval(sid, 1, c_loc(rloc), omin, omax, i_hth, i_hz, c_loc(hout))
    end subroutine eval_hthz_c

    subroutine eval_hthz(rval, omin, omax, bp, hout) bind(C, name="eval_hthz")
        real(c_double), value :: rval
        integer(c_int), value :: omin, omax
        type(c_ptr), value :: bp
        real(c_double), intent(out), target :: hout(*)
        real(dp), target :: rloc(1)
        rloc(1) = rval
        call spline_eval(sid, 1, c_loc(rloc), omin, omax, i_hth, i_hz, c_loc(hout))
    end subroutine eval_hthz

    subroutine eval_Bt_dBt_Bz_dBz(rval, bp, R) bind(C, name="eval_Bt_dBt_Bz_dBz")
        real(c_double), value :: rval
        type(c_ptr), value :: bp
        real(c_double), intent(out), target :: R(*)
        real(dp), target :: rloc(1)
        rloc(1) = rval
        call spline_eval(sid, 1, c_loc(rloc), 0, 1, i_Bth, i_Bz, c_loc(R))
    end subroutine eval_Bt_dBt_Bz_dBz

    subroutine eval_p_dp(rval, bp, R) bind(C, name="eval_p_dp")
        real(c_double), value :: rval
        type(c_ptr), value :: bp
        real(c_double), intent(out), target :: R(*)
        real(dp), target :: rloc(1)
        rloc(1) = rval
        call eval_p_dp_local(rloc(1), R)
    end subroutine eval_p_dp

    !> Shared by eval_p_dp (external C++ callers) and
    !> check_and_spline_equilibrium (internal use).
    subroutine eval_p_dp_local(r, Rout)
        real(dp), intent(in) :: r
        real(dp), intent(out) :: Rout(0:1)
        real(dp), target :: S(0:3)
        real(dp) :: n_, dn, Ti, dTi, Te, dTe
        real(dp), target :: rloc(1)
        rloc(1) = r
        call spline_eval(sid, 1, c_loc(rloc), 0, 1, i_n, i_n, c_loc(S))
        n_ = S(0); dn = S(1)
        call spline_eval(sid, 1, c_loc(rloc), 0, 1, i_T(0), i_T(1), c_loc(S))
        Ti = S(0); dTi = S(1); Te = S(2); dTe = S(3)
        Rout(0) = boltz*(n_*(Ti + Te))
        Rout(1) = boltz*(dn*(Ti + Te) + n_*(dTi + dTe))
    end subroutine eval_p_dp_local

    subroutine eval_mass_density(rval, bp, R) bind(C, name="eval_mass_density")
        real(c_double), value :: rval
        type(c_ptr), value :: bp
        real(c_double), intent(out), target :: R(*)
        real(dp), target :: rloc(1)
        rloc(1) = rval
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_n, i_n, c_loc(R))
        R(1) = R(1)*mass(0)
    end subroutine eval_mass_density

    subroutine eval_Bt_Bz(rval, bp, R) bind(C, name="eval_Bt_Bz")
        real(c_double), value :: rval
        type(c_ptr), value :: bp
        real(c_double), intent(out), target :: R(*)
        real(dp), target :: rloc(1)
        rloc(1) = rval
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_Bth, i_Bz, c_loc(R))
    end subroutine eval_Bt_Bz

    subroutine eval_Bt(rval, bp, R) bind(C, name="eval_Bt")
        real(c_double), value :: rval
        type(c_ptr), value :: bp
        real(c_double), intent(out), target :: R(*)
        real(dp), target :: rloc(1)
        rloc(1) = rval
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_Bth, i_Bth, c_loc(R))
    end subroutine eval_Bt

    subroutine eval_Bz(rval, bp, R) bind(C, name="eval_Bz")
        real(c_double), value :: rval
        type(c_ptr), value :: bp
        real(c_double), intent(out), target :: R(*)
        real(dp), target :: rloc(1)
        rloc(1) = rval
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_Bz, i_Bz, c_loc(R))
    end subroutine eval_Bz

    subroutine eval_Vt(rval, bp, R) bind(C, name="eval_Vt")
        real(c_double), value :: rval
        type(c_ptr), value :: bp
        real(c_double), intent(out), target :: R(*)
        real(dp), target :: rloc(1)
        rloc(1) = rval
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_Vth(2), i_Vth(2), c_loc(R))
    end subroutine eval_Vt

    subroutine eval_Vz(rval, bp, R) bind(C, name="eval_Vz")
        real(c_double), value :: rval
        type(c_ptr), value :: bp
        real(c_double), intent(out), target :: R(*)
        real(dp), target :: rloc(1)
        rloc(1) = rval
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_Vz(2), i_Vz(2), c_loc(R))
    end subroutine eval_Vz

    subroutine eval_B0_ht_hz_n0_Vz(rval, bp, R) bind(C, name="eval_B0_ht_hz_n0_Vz")
        real(c_double), value :: rval
        type(c_ptr), value :: bp
        real(c_double), intent(out), target :: R(0:4)
        real(dp), target :: rloc(1)
        rloc(1) = rval
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_B, i_hz, c_loc(R(0:2)))
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_n, i_n, c_loc(R(3:3)))
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_Vz(2), i_Vz(2), c_loc(R(4:4)))
    end subroutine eval_B0_ht_hz_n0_Vz

    !> Live (confirmed via calc_mode.cpp:30,41,70 - missed in an earlier
    !> "0 callers" grep pass because the call sites have no space before
    !> the opening paren, q(x,...), which the original pattern match
    !> assumed; always re-verify a "dead" claim by reading consumer code,
    !> not just trusting one grep pattern).
    !> Matches wave_code_interface.cpp's get_background_magnetic_fields_
    !> from_wave_code_ inline spline_eval_ (i_Bth..i_B, 3 consecutive
    !> channels) and get_collision_frequences_from_wave_code_'s two
    !> single-channel calls (i_nu(0), i_nu(1) - NOT consecutive in the
    !> index table, so two separate evaluations).
    subroutine get_background_magnetic_fields(rval, Bt, Bz, B0) &
        bind(C, name="get_background_magnetic_fields_")
        real(c_double), value :: rval
        real(c_double), intent(out), target :: Bt(1), Bz(1), B0(1)
        real(dp), target :: rloc(1), Rl(0:2)
        rloc(1) = rval
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_Bth, i_B, c_loc(Rl))
        Bt(1) = Rl(0); Bz(1) = Rl(1); B0(1) = Rl(2)
    end subroutine get_background_magnetic_fields

    subroutine get_background_collision_freqs(rval, nui, nue) &
        bind(C, name="get_background_collision_freqs_")
        real(c_double), value :: rval
        real(c_double), intent(out), target :: nui(1), nue(1)
        real(dp), target :: rloc(1)
        rloc(1) = rval
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_nu(0), i_nu(0), c_loc(nui))
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_nu(1), i_nu(1), c_loc(nue))
    end subroutine get_background_collision_freqs

    real(c_double) function eval_q(rval, bp) bind(C, name="q")
        real(c_double), value :: rval
        type(c_ptr), value :: bp
        real(dp), target :: rloc(1), ans(1)
        rloc(1) = rval
        call spline_eval(sid, 1, c_loc(rloc), 0, 0, i_q, i_q, c_loc(ans))
        eval_q = ans(1)
    end function eval_q

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

end module kilca_background_data_m
