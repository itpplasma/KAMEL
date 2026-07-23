!> Eigenmode root-search strategies, formerly eigmode/calc_eigmode.cpp and
!> eigmode/find_eigmodes.cpp. Three runtime-selected paths, dispatched by
!> get_eigmode_search_flag_ from core_data_m's
!> core_data_calc_and_set_mode_dependent_eigmode_:
!>
!>   find_det_zeros       - multi-start hybrid Newton solve via fortnum's
!>                          fortnum_multiroot_hybrid, with an optional winding
!>                          number sanity check (calc_circle_integral/inv_det).
!>   loop_over_frequences - brute-force 2D (re f, im f) determinant table dump.
!>   find_eigmodes        - reliable zeros search driving the vendored ZerSol
!>                          C++ library through its extern "C" wrapper
!>                          (zersolc.h); the templated zerosolver.hpp stays C++.
!>
!> The shared determinant evaluation (the oracle's eval_det and templated
!> determinant<T>, near-identical bodies) lives in compute_det; the two
!> callback ABIs (fortnum_vector_fn and ZerSol's complex_function) are thin
!> bind(C) shims over eval_determinant_core.
module kilca_eigmode_solve_m
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex, &
        c_intptr_t, c_ptr, c_funptr, c_char, c_null_char, c_loc, c_funloc, &
        c_f_pointer, c_null_funptr
    use, intrinsic :: iso_fortran_env, only: output_unit
    use constants, only: dp, pi, im
    use eigmode_sett_data, only: es_fname => fname, es_rdim => rdim, &
        es_idim => idim, es_rfmin => rfmin, es_rfmax => rfmax, &
        es_ifmin => ifmin, es_ifmax => ifmax, es_eps_res => eps_res, &
        es_eps_abs => eps_abs, es_eps_rel => eps_rel, es_delta => delta, &
        es_test_roots => test_roots, es_nguess => Nguess, es_kmin => kmin, &
        es_kmax => kmax, es_fstart => fstart, es_n_zeros => n_zeros, &
        es_use_winding => use_winding
    use kilca_core_data_m, only: core_data_get_sd_, core_data_get_bp_, &
        core_data_set_mda_element_, core_data_get_mda_element_
    use kilca_settings_m, only: settings_get_path2project_
    use kilca_mode_data_m, only: mode_data_create_, mode_data_destroy_, &
        mode_data_calc_all_mode_data_, mode_data_get_wd_
    use kilca_wave_data_m, only: wave_data_get_det_re, wave_data_get_det_im
    implicit none
    private

    public :: find_det_zeros, loop_over_frequences, find_eigmodes
    public :: calc_det

    !> Mirrors the C struct det_params { int ind, m, n; double delta;
    !> intptr_t cd; } field-for-field. Field order and interoperable kinds
    !> must match the C layout exactly. delta is carried for layout fidelity
    !> but read by neither callback (the oracle never read it either).
    type, bind(C) :: det_params_t
        integer(c_int) :: ind
        integer(c_int) :: m
        integer(c_int) :: n
        real(c_double) :: delta
        integer(c_intptr_t) :: cd
    end type det_params_t

    interface
        !> Legacy free Fortran subroutine in mode/mode_m.f90 (single trailing
        !> underscore), not a module procedure - declared bind(C) here.
        subroutine clear_all_data_in_mode_data_module() &
                bind(C, name="clear_all_data_in_mode_data_module_")
        end subroutine clear_all_data_in_mode_data_module

        !> fortnum multidimensional hybrid Newton solve, F(x)=0 in R^n. fdf is
        !> the residual-only callback (fortnum_vector_fn); the wrapper builds
        !> the Jacobian by central differences.
        function fortnum_multiroot_hybrid(fdf, n, x0, xtol, ftol, max_iter, x, ctx) &
                result(stat) bind(C, name="fortnum_multiroot_hybrid")
            import :: c_int, c_double, c_funptr, c_ptr
            type(c_funptr), value :: fdf
            integer(c_int), value :: n
            real(c_double), intent(in) :: x0(*)
            real(c_double), value :: xtol, ftol
            integer(c_int), value :: max_iter
            real(c_double), intent(out) :: x(*)
            type(c_ptr), value :: ctx
            integer(c_int) :: stat
        end function fortnum_multiroot_hybrid

        !> ZerSol extern "C" wrapper (math/zersol-0.0.0/src/zersolc.h),
        !> already specialized to _Real = double.
        function create_solver(f, df, data, description, xmin, xmax, ymin, ymax) &
                result(solver) bind(C, name="create_solver")
            import :: c_funptr, c_ptr, c_char, c_double
            type(c_funptr), value :: f, df
            type(c_ptr), value :: data
            character(kind=c_char), intent(in) :: description(*)
            real(c_double), value :: xmin, xmax, ymin, ymax
            type(c_ptr) :: solver
        end function create_solver

        function find_zeros(solver, n, z, v, n_zeros) result(stat) &
                bind(C, name="find_zeros")
            import :: c_ptr, c_int, c_double_complex
            type(c_ptr), value :: solver
            integer(c_int), value :: n
            complex(c_double_complex), intent(out) :: z(*), v(*)
            integer(c_int), intent(out) :: n_zeros
            integer(c_int) :: stat
        end function find_zeros

        subroutine print_status(solver, file_name) bind(C, name="print_status")
            import :: c_ptr, c_char
            type(c_ptr), value :: solver
            character(kind=c_char), intent(in) :: file_name(*)
        end subroutine print_status

        subroutine print_dump(solver, file_name) bind(C, name="print_dump")
            import :: c_ptr, c_char
            type(c_ptr), value :: solver
            character(kind=c_char), intent(in) :: file_name(*)
        end subroutine print_dump

        subroutine free_solver(solver) bind(C, name="free_solver")
            import :: c_ptr
            type(c_ptr), value :: solver
        end subroutine free_solver

        subroutine set_nd_method(solver, method) bind(C, name="set_nd_method")
            import :: c_ptr, c_int
            type(c_ptr), value :: solver
            integer(c_int), value :: method
        end subroutine set_nd_method

        subroutine set_nd_step(solver, h) bind(C, name="set_nd_step")
            import :: c_ptr, c_double
            type(c_ptr), value :: solver
            real(c_double), value :: h
        end subroutine set_nd_step

        subroutine set_min_rec_lev(solver, min_rec_lev) bind(C, name="set_min_rec_lev")
            import :: c_ptr, c_int
            type(c_ptr), value :: solver
            integer(c_int), value :: min_rec_lev
        end subroutine set_min_rec_lev

        subroutine set_n_target(solver, n_target) bind(C, name="set_n_target")
            import :: c_ptr, c_int
            type(c_ptr), value :: solver
            integer(c_int), value :: n_target
        end subroutine set_n_target

        subroutine set_start_array(solver, n_start, start) &
                bind(C, name="set_start_array")
            import :: c_ptr, c_int, c_double_complex
            type(c_ptr), value :: solver
            integer(c_int), value :: n_start
            complex(c_double_complex), intent(in) :: start(*)
        end subroutine set_start_array

        subroutine set_n_split_x(solver, n_split_x) bind(C, name="set_n_split_x")
            import :: c_ptr, c_int
            type(c_ptr), value :: solver
            integer(c_int), value :: n_split_x
        end subroutine set_n_split_x

        subroutine set_n_split_y(solver, n_split_y) bind(C, name="set_n_split_y")
            import :: c_ptr, c_int
            type(c_ptr), value :: solver
            integer(c_int), value :: n_split_y
        end subroutine set_n_split_y

        subroutine set_max_iter_num(solver, max_iter_num) &
                bind(C, name="set_max_iter_num")
            import :: c_ptr, c_int
            type(c_ptr), value :: solver
            integer(c_int), value :: max_iter_num
        end subroutine set_max_iter_num

        subroutine set_eps_for_arg(solver, abs_eps_z, rel_eps_z) &
                bind(C, name="set_eps_for_arg")
            import :: c_ptr, c_double
            type(c_ptr), value :: solver
            real(c_double), value :: abs_eps_z, rel_eps_z
        end subroutine set_eps_for_arg

        subroutine set_eps_for_func(solver, abs_eps_f, rel_eps_f) &
                bind(C, name="set_eps_for_func")
            import :: c_ptr, c_double
            type(c_ptr), value :: solver
            real(c_double), value :: abs_eps_f, rel_eps_f
        end subroutine set_eps_for_func

        subroutine set_use_winding(solver, use_winding) bind(C, name="set_use_winding")
            import :: c_ptr, c_int
            type(c_ptr), value :: solver
            integer(c_int), value :: use_winding
        end subroutine set_use_winding

        subroutine set_debug_level(solver, debug_level) bind(C, name="set_debug_level")
            import :: c_ptr, c_int
            type(c_ptr), value :: solver
            integer(c_int), value :: debug_level
        end subroutine set_debug_level

        subroutine set_print_level(solver, print_level) bind(C, name="set_print_level")
            import :: c_ptr, c_int
            type(c_ptr), value :: solver
            integer(c_int), value :: print_level
        end subroutine set_print_level
    end interface

contains

    !> --- multi-start hybrid Newton search (oracle find_det_zeros) ---

    function find_det_zeros(ind, m, n, cd) bind(C, name="find_det_zeros") result(stat)
        integer(c_int), value :: ind, m, n
        integer(c_intptr_t), value :: cd
        integer(c_int) :: stat

        type(det_params_t), target :: p
        character(len=1024) :: full_name
        integer :: u, k
        real(c_double) :: x_init(2), x_root(2), det_out(2)
        real(dp) :: radius
        complex(dp) :: center, quad1, quad2

        stat = 0

        p%ind = ind
        p%m = m
        p%n = n
        p%delta = es_delta
        p%cd = cd

        call build_output_name(cd, full_name)

        open (newunit=u, file=trim(full_name), status='replace', action='write')
        write (u, '(a)') '%Re(f)'//tab()//tab()//tab()//'Im(f)'//tab()//tab()//tab()// &
            'Re(det)'//tab()//tab()//tab()//'Im(det)'
        close (u)

        do k = es_kmin, es_kmax
            open (newunit=u, file=trim(full_name), position='append', action='write')

            x_init(1) = real(es_fstart(k), c_double)
            x_init(2) = aimag(es_fstart(k))

            stat = fortnum_multiroot_hybrid(c_funloc(eval_det_residual), 2_c_int, &
                x_init, real(es_eps_abs, c_double), real(es_eps_res, c_double), &
                1000000_c_int, x_root, c_loc(p))

            call eval_determinant_core(p, real(x_root(1), dp), real(x_root(2), dp), &
                det_out(1), det_out(2))

            write (u, '(es27.20e2, 2x, es27.20e2, a, es27.20e2, 2x, es27.20e2)') &
                x_root(1), x_root(2), tab(), det_out(1), det_out(2)
            write (u, '(a, i0)') '%status = ', stat

            if (es_test_roots == 1) then
                radius = 1.0e1_dp

                center = cmplx(x_root(1), x_root(2), dp)
                quad1 = calc_circle_integral(p, center, radius)

                center = center + 5.0_dp*radius*(1.0_dp + im)
                quad2 = calc_circle_integral(p, center, radius)

                write (output_unit, &
                    '(a, es27.20e2, 1x, es27.20e2, a, es27.20e2, 1x, es27.20e2)') &
                    'quad1 = ', real(quad1, dp), aimag(quad1), tab()//'quad2 = ', &
                    real(quad2, dp), aimag(quad2)

                write (u, &
                    '(a, es27.20e2, 1x, es27.20e2, a, es27.20e2, 1x, '// &
                    'es27.20e2, a, es27.20e2)') &
                    '%quad1 = ', real(quad1, dp), aimag(quad1), tab()//'quad2 = ', &
                    real(quad2, dp), aimag(quad2), tab()//'err = ', abs(quad1/quad2)
            end if

            close (u)
        end do
    end function find_det_zeros

    !> Closed circle quadrature of inv_det around center (oracle
    !> calc_circle_integral). The oracle's function-pointer argument had a
    !> single possible value (inv_det), so it is called directly here.
    function calc_circle_integral(p, center, radius) result(res)
        type(det_params_t), intent(in) :: p
        complex(dp), intent(in) :: center
        real(dp), intent(in) :: radius
        complex(dp) :: res, expon
        real(dp) :: dphi
        integer, parameter :: n = 100
        integer :: k

        dphi = 2.0_dp*pi/n
        res = (0.0_dp, 0.0_dp)
        do k = 0, n - 1
            expon = exp(k*dphi*im)
            res = res + inv_det(p, center + radius*expon)*expon
        end do
        res = dphi*im*radius*res
    end function calc_circle_integral

    !> Reciprocal of the determinant at z (oracle inv_det, E/det with E = 1).
    function inv_det(p, z) result(res)
        type(det_params_t), intent(in) :: p
        complex(dp), intent(in) :: z
        complex(dp) :: res, det
        real(dp) :: det_re, det_im

        call eval_determinant_core(p, real(z, dp), aimag(z), det_re, det_im)
        det = cmplx(det_re, det_im, dp)
        res = (1.0_dp, 0.0_dp)/det
    end function inv_det

    !> Imaginary part of the determinant on the imaginary frequency axis
    !> (oracle calc_det). Kept as an exported entry point, matching the
    !> oracle header, though no caller exists in the codebase.
    subroutine calc_det(p, freq, absdet)
        type(det_params_t), intent(in) :: p
        real(dp), intent(in) :: freq
        real(dp), intent(out) :: absdet
        real(dp) :: det_re, det_im

        call eval_determinant_core(p, 0.0_dp, freq, det_re, det_im)
        absdet = det_im
    end subroutine calc_det

    !> --- brute-force 2D frequency scan (oracle loop_over_frequences) ---

    function loop_over_frequences(ind, m, n, cd) bind(C, name="loop_over_frequences") &
            result(stat)
        integer(c_int), value :: ind, m, n
        integer(c_intptr_t), value :: cd
        integer(c_int) :: stat

        character(len=1024) :: full_name
        integer :: u, i, k
        real(dp) :: fre, fim, det_re, det_im

        call build_output_name(cd, full_name)

        open (newunit=u, file=trim(full_name), status='replace', action='write')
        write (u, '(a)') '%iter'//tab()//'Re(f)'//tab()//tab()//tab()//'Im(f)'// &
            tab()//tab()//tab()//'Re(det)'//tab()//tab()//tab()//'Im(det)'

        do i = 0, es_rdim - 1
            fre = es_rfmin + real(i, dp)*(es_rfmax - es_rfmin) &
                /real(max(es_rdim - 1, 1), dp)

            do k = 0, es_idim - 1
                fim = es_ifmin + real(k, dp)*(es_ifmax - es_ifmin) &
                    /real(max(es_idim - 1, 1), dp)

                call compute_det(ind, m, n, cd, fre, fim, det_re, det_im)

                write (u, '(i6, a, es27.20e2, 2x, es27.20e2, a, '// &
                    'es27.20e2, 2x, es27.20e2)') &
                    i*es_idim + k, tab(), fre, fim, tab(), det_re, det_im
            end do
        end do

        close (u)
        stat = 0
    end function loop_over_frequences

    !> --- reliable zeros search via the ZerSol C wrapper (oracle find_eigmodes) ---

    function find_eigmodes(ind, m, n, cd) bind(C, name="find_eigmodes") result(stat)
        integer(c_int), value :: ind, m, n
        integer(c_intptr_t), value :: cd
        integer(c_int) :: stat

        integer(c_int), parameter :: max_n_zeros = 128
        type(det_params_t), target :: p
        type(c_ptr) :: solver
        complex(c_double_complex) :: z(0:max_n_zeros - 1), v(0:max_n_zeros - 1)
        complex(c_double_complex), allocatable :: es_fstart_loc(:)
        integer(c_int) :: n_zeros
        character(len=1024) :: full_name
        integer :: u, i, k

        p%ind = ind
        p%m = m
        p%n = n
        p%delta = es_delta
        p%cd = cd

        solver = create_solver(c_funloc(determinant_cb), c_null_funptr, c_loc(p), &
            'determinant'//c_null_char, real(es_rfmin, c_double), &
            real(es_rfmax, c_double), real(es_ifmin, c_double), &
            real(es_ifmax, c_double))

        call set_nd_method(solver, 1_c_int)
        if (es_delta > 0.0_dp) call set_nd_step(solver, real(es_delta, c_double))

        call set_min_rec_lev(solver, 8_c_int)

        allocate (es_fstart_loc(0:es_nguess - 1))
        do k = 0, es_nguess - 1
            es_fstart_loc(k) = cmplx(real(es_fstart(k), c_double), &
                aimag(es_fstart(k)), c_double_complex)
        end do

        call set_n_target(solver, int(es_n_zeros, c_int))
        call set_start_array(solver, int(es_nguess, c_int), es_fstart_loc)
        call set_n_split_x(solver, 4_c_int)
        call set_n_split_y(solver, 4_c_int)
        call set_max_iter_num(solver, 24_c_int)
        call set_eps_for_arg(solver, real(es_eps_abs, c_double), &
            real(es_eps_rel, c_double))
        call set_eps_for_func(solver, 0.0_c_double, real(es_eps_res, c_double))
        call set_use_winding(solver, int(es_use_winding, c_int))
        call set_debug_level(solver, 1_c_int)
        call set_print_level(solver, 1_c_int)

        n_zeros = 0
        stat = find_zeros(solver, max_n_zeros, z, v, n_zeros)

        if (stat /= 0) call print_dump(solver, c_null_char)

        call print_status(solver, c_null_char)
        call print_dump(solver, 'search.dump'//c_null_char)

        call free_solver(solver)

        deallocate (es_fstart_loc)

        call build_output_name(cd, full_name)

        open (newunit=u, file=trim(full_name), status='replace', action='write')
        write (u, '(a)') '%#'//tab()//'Re(f)'//tab()//tab()//tab()//tab()//'Im(f)'// &
            tab()//tab()//tab()//tab()//'Re(det)'//tab()//tab()//tab()//tab()//'Im(det)'

        do i = 0, n_zeros - 1
            write (u, '(i5, a, es27.20e2, 2x, es27.20e2, a, '// &
                'es27.20e2, 2x, es27.20e2)') &
                i, tab(), real(z(i), dp), aimag(z(i)), tab(), &
                real(v(i), dp), aimag(v(i))
        end do

        close (u)
        stat = 0
    end function find_eigmodes

    !> --- shared determinant evaluation ---

    !> fortnum_vector_fn residual shim: x = (Re f, Im f), f = (Re det, Im det).
    subroutine eval_det_residual(n, x, f, params) bind(C)
        integer(c_int), value :: n
        real(c_double), intent(in) :: x(*)
        real(c_double), intent(out) :: f(*)
        type(c_ptr), value :: params
        type(det_params_t), pointer :: p

        call c_f_pointer(params, p)
        call eval_determinant_core(p, real(x(1), dp), real(x(2), dp), f(1), f(2))
    end subroutine eval_det_residual

    !> ZerSol complex_function shim: z (complex frequency) -> determinant value.
    function determinant_cb(z, params) result(res) bind(C)
        complex(c_double_complex), value :: z
        type(c_ptr), value :: params
        complex(c_double_complex) :: res
        type(det_params_t), pointer :: p
        real(dp) :: det_re, det_im

        call c_f_pointer(params, p)
        call eval_determinant_core(p, real(z, dp), aimag(z), det_re, det_im)
        res = cmplx(det_re, det_im, c_double_complex)
    end function determinant_cb

    !> Determinant at (fre, fim) plus the oracle's det.dat debug append, the
    !> body shared by eval_det and the templated determinant<T>.
    subroutine eval_determinant_core(p, fre, fim, det_re, det_im)
        type(det_params_t), intent(in) :: p
        real(dp), intent(in) :: fre, fim
        real(dp), intent(out) :: det_re, det_im
        integer :: u

        call compute_det(p%ind, p%m, p%n, p%cd, fre, fim, det_re, det_im)

        open (newunit=u, file='det.dat', position='append', action='write')
        write (u, '(es27.20e2, 1x, es27.20e2, a, es27.20e2, 1x, es27.20e2)') &
            fre, fim, tab(), det_re, det_im
        close (u)
    end subroutine eval_determinant_core

    !> Build mode_data at olab = 2 pi (fre + i fim), evaluate the dispersion
    !> determinant, then tear the mode_data back down - the create / calc /
    !> read / destroy sequence common to all three search paths.
    subroutine compute_det(ind, m, n, cd, fre, fim, det_re, det_im)
        integer(c_int), intent(in) :: ind, m, n
        integer(c_intptr_t), intent(in) :: cd
        real(dp), intent(in) :: fre, fim
        real(dp), intent(out) :: det_re, det_im
        real(c_double) :: olab_re, olab_im
        integer(c_intptr_t) :: sd, bp, mdh, wdh
        type(c_ptr) :: wd_cptr
        character(kind=c_char) :: path_c(1024)

        olab_re = 2.0_dp*pi*fre
        olab_im = 2.0_dp*pi*fim

        sd = core_data_get_sd_(cd)
        bp = core_data_get_bp_(cd)

        call settings_get_path2project_(sd, path_c)

        mdh = mode_data_create_(m, n, olab_re, olab_im, sd, bp, path_c)
        call core_data_set_mda_element_(cd, ind, mdh)

        call mode_data_calc_all_mode_data_(core_data_get_mda_element_(cd, ind), 0_c_int)

        wdh = mode_data_get_wd_(core_data_get_mda_element_(cd, ind))
        wd_cptr = transfer(wdh, wd_cptr)
        det_re = wave_data_get_det_re(wd_cptr)
        det_im = wave_data_get_det_im(wd_cptr)

        call mode_data_destroy_(core_data_get_mda_element_(cd, ind))
        call core_data_set_mda_element_(cd, ind, 0_c_intptr_t)
        call clear_all_data_in_mode_data_module()
    end subroutine compute_det

    !> --- helpers ---

    !> full output path: settings path2project followed by the eigmode file
    !> name (oracle sprintf("%s%s", path2project, es_fname)).
    subroutine build_output_name(cd, full_name)
        integer(c_intptr_t), intent(in) :: cd
        character(len=*), intent(out) :: full_name
        character(kind=c_char) :: path_c(1024)
        character(len=1024) :: path
        integer :: i

        call settings_get_path2project_(core_data_get_sd_(cd), path_c)

        path = ''
        do i = 1, 1024
            if (path_c(i) == c_null_char) exit
            path(i:i) = path_c(i)
        end do

        full_name = trim(path)//trim(es_fname)
    end subroutine build_output_name

    pure function tab() result(c)
        character(len=1) :: c
        c = achar(9)
    end function tab

end module kilca_eigmode_solve_m
