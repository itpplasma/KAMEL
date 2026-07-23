!> Per-mode orchestrator, formerly the C++ mode_data class (mode.{h,cpp})
!> and calc_mode.cpp. Owns one perturbation mode's zone chain, drives basis
!> calculation, builds and solves the stitching-equations linear system at
!> zone boundaries, and combines/interpolates the final lab-frame fields.
!>
!> zones: stored as a plain integer(c_intptr_t) handle array (NOT a
!> class(zone_t) array -- Fortran cannot mix dynamic types in one array, and
!> NOT the private zone_box_t pool type), dispatched via handle_to_zone right
!> before each use, mirroring the oracle's `zone *code = zones[iz]` pattern.
!>
!> EB/EB_int collapse: the oracle kept two differently-laid-out flat arrays,
!> EB (node-major, from iFFM) for storage/output and EB_int (node-minor, from
!> iFint) purely so Neville interpolation could see a node-contiguous y-grid
!> per (component, re/im). A native complex(dp) EB(6,dim) array makes that
!> reshuffle unnecessary: for fixed comp, real(EB(comp,a:b))/aimag(EB(comp,
!> a:b)) is exactly the node-contiguous slice eval_neville_polynom needs
!> (Fortran copies non-contiguous array-valued actual arguments transparently
!> for an INTENT(IN) dummy), so setup_wave_fields_for_interpolation becomes a
!> no-op and is not translated as a separate procedure; its call-sequence
!> point is kept as a comment in mode_data_calc_all_mode_data_.
!>
!> copy_mode_paths_to_mode_data_module_/copy_mode_paths_from_mode_data_struct_:
!> the oracle round-tripped path2linear/path2dispersion/path2poincare into
!> the legacy `mode_data` Fortran module (mode_m.f90) via a C++ method that
!> dereferenced `(mode_data **)this`. With mode_data itself now Fortran there
!> is no such C++ object, so this module writes those three module variables
!> directly (`use mode_data, only: ... => path2linear, ...`) instead of
!> calling the now-impossible round trip; the legacy module ends up holding
!> the exact same strings at the exact same point in construction.
!>
!> Directory scanning (allocate_and_setup_zones' zone_*.in discovery) is
!> translated via direct POSIX opendir/readdir/closedir/fnmatch C interop
!> (glibc x86_64 struct dirent layout, matching this build's platform
!> exactly, same as the oracle's own dirent.h/fnmatch.h includes), preserving
!> readdir's filesystem-order-dependent first-match semantics byte for byte
!> rather than substituting a sorted listing.
module kilca_mode_data_m
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_double, c_char, &
        c_ptr, c_funptr, c_funloc, c_loc, c_f_pointer, c_null_char, c_null_ptr, &
        c_associated, c_long, c_short, c_signed_char
    use constants, only: dp, twopi
    use kilca_zone_m, only: zone_t, handle_to_zone, zone_destroy_c, &
        PLASMA_MODEL_VACUUM, PLASMA_MODEL_MEDIUM, PLASMA_MODEL_IMHD, &
        PLASMA_MODEL_RMHD, PLASMA_MODEL_FLRE, &
        BOUNDARY_CENTER, BOUNDARY_INFINITY, BOUNDARY_IDEALWALL, &
        BOUNDARY_INTERFACE, BOUNDARY_ANTENNA, &
        Nmed, med_str, skip_line, read_real_before_hash, read_token_before_hash
    use kilca_hmedium_zone_m, only: hmedium_zone_create
    use kilca_imhd_zone_m, only: imhd_zone_create
    use kilca_flre_zone_m, only: flre_zone_create_
    use kilca_wave_data_m, only: wave_data_t, wave_data_create, wave_data_destroy
    use kilca_background_data_m, only: get_background_x0, get_background_xlast
    use kilca_neville_m, only: eval_neville_polynom
    use fortnum_capi, only: fortnum_root_brent
    use fortnum_status, only: FORTNUM_OK
    use kilca_progs_common_m, only: fmt_g, fmt_e
    implicit none
    private

    public :: mode_data_t
    public :: mode_data_create_, mode_data_destroy_
    public :: mode_data_calc_all_mode_data_
    public :: mode_data_get_wd_, mode_data_get_zone_handle_
    public :: mode_data_eval_EB_fields_, mode_data_eval_diss_power_density_
    public :: mode_data_eval_current_density_

    type :: mode_data_t
        integer(c_intptr_t) :: bp = 1_c_intptr_t
        character(len=1024) :: path2project = ''
        type(wave_data_t), pointer :: wd => null()
        character(len=1024) :: path2linear = ''
        character(len=1024) :: path2dispersion = ''
        character(len=1024) :: path2poincare = ''
        integer :: Nzones = 0
        integer(c_intptr_t), allocatable :: zone_handles(:)
        integer :: dim = 0
        real(dp), allocatable :: r(:)
        complex(dp), allocatable :: EB(:, :)
        integer, allocatable :: index(:)
        integer :: Nc = 0
        complex(dp), allocatable :: A(:, :)
        complex(dp), allocatable :: B(:)
        complex(dp), allocatable :: S(:)
    end type mode_data_t

    !> glibc x86_64 struct dirent (default 64-bit ino_t/off_t on LP64 Linux):
    !> 8+8+2+1 bytes of header then a NUL-terminated name, no padding before
    !> d_name since char has alignment 1. offsetof(d_name) == 19, confirmed
    !> against this platform's <dirent.h> via offsetof().
    type, bind(c) :: dirent_t
        integer(c_long) :: d_ino
        integer(c_long) :: d_off
        integer(c_short) :: d_reclen
        integer(c_signed_char) :: d_type
        character(kind=c_char) :: d_name(256)
    end type dirent_t

    interface
        integer(c_int) function get_output_flag_background() &
                bind(C, name="get_output_flag_background_")
            import :: c_int
        end function get_output_flag_background

        integer(c_int) function get_output_flag_emfield() &
                bind(C, name="get_output_flag_emfield_")
            import :: c_int
        end function get_output_flag_emfield

        integer(c_int) function get_output_flag_additional() &
                bind(C, name="get_output_flag_additional_")
            import :: c_int
        end function get_output_flag_additional

        integer(c_int) function get_output_flag_dispersion() &
                bind(C, name="get_output_flag_dispersion_")
            import :: c_int
        end function get_output_flag_dispersion

        real(c_double) function get_background_rtor() bind(C, name="get_background_rtor_")
            import :: c_double
        end function get_background_rtor

        real(c_double) function get_background_V_gal_sys() &
                bind(C, name="get_background_V_gal_sys_")
            import :: c_double
        end function get_background_V_gal_sys

        real(c_double) function eval_q_for_resonance(rval, bp) bind(C, name="q")
            import :: c_double, c_ptr
            real(c_double), value :: rval
            type(c_ptr), value :: bp
        end function eval_q_for_resonance

        function opendir_c(path) bind(C, name="opendir") result(dirp)
            import :: c_char, c_ptr
            character(kind=c_char), intent(in) :: path(*)
            type(c_ptr) :: dirp
        end function opendir_c

        function readdir_c(dirp) bind(C, name="readdir") result(ep)
            import :: c_ptr
            type(c_ptr), value :: dirp
            type(c_ptr) :: ep
        end function readdir_c

        integer(c_int) function closedir_c(dirp) bind(C, name="closedir")
            import :: c_ptr, c_int
            type(c_ptr), value :: dirp
        end function closedir_c

        integer(c_int) function fnmatch_c(pattern, str, flags) bind(C, name="fnmatch")
            import :: c_char, c_int
            character(kind=c_char), intent(in) :: pattern(*), str(*)
            integer(c_int), value :: flags
        end function fnmatch_c
    end interface

    !> Pre-existing bare (non-module) legacy Fortran subroutines: zone
    !> boundary-equation builders (hmedium.f90/imhd.f90/flre.f90) and the
    !> stitching-system assembly/solve routines (stitching.f90). Called via
    !> implicit interface, exactly as their own dummy argument lists declare.
    external :: center_equations_hommed, infinity_equations_hommed, &
        ideal_wall_equations_hommed, stitching_equations_hommed_hommed
    external :: center_equations_imhd, infinity_equations_imhd, &
        ideal_wall_equations_imhd, stitching_equations_imhd_hommed, &
        stitching_equations_imhd_imhd
    external :: center_equations_flre, infinity_equations_flre, &
        ideal_wall_equations_flre, stitching_equations_flre_hommed, &
        stitching_equations_flre_flre
    external :: update_system_matrix_and_rhs_vector, calc_system_determinant, &
        find_superposition_coeffs

    !> Pre-existing legacy mode_data module setters (mode_m.f90), bare
    !> external subroutines (not module procedures), called natively.
    external :: set_settings_in_mode_data_module, set_back_profiles_in_mode_data_module, &
        set_wave_data_in_mode_data_module, set_wave_parameters_in_mode_data_module, &
        set_resonance_location_in_mode_data_module

contains

    !> ---- construction / destruction ----

    function mode_data_create_(m, n, olab_re, olab_im, sd_ptr, bp_ptr, path2project) &
            result(handle) bind(C, name="mode_data_create_")
        integer(c_int), value :: m, n
        real(c_double), value :: olab_re, olab_im
        integer(c_intptr_t), value :: sd_ptr, bp_ptr
        character(kind=c_char), intent(in) :: path2project(*)
        integer(c_intptr_t) :: handle

        type(mode_data_t), pointer :: md
        complex(dp) :: olab, omov
        integer(c_intptr_t) :: wd_handle
        type(c_ptr) :: wd_cptr
        real(dp) :: r_res_local

        allocate (md)
        md%bp = bp_ptr
        md%path2project = c_string_to_fortran_local(path2project)

        call set_settings_in_mode_data_module(sd_ptr)
        call set_back_profiles_in_mode_data_module(bp_ptr)

        olab = cmplx(olab_re, olab_im, dp)
        omov = olab - real(n, dp)*get_background_V_gal_sys()/get_background_rtor()

        wd_handle = wave_data_create(m, n, olab_re, olab_im, real(omov, dp), aimag(omov))
        wd_cptr = transfer(wd_handle, wd_cptr)
        call c_f_pointer(wd_cptr, md%wd)

        call set_wave_data_in_mode_data_module(wd_handle)

        call set_wave_parameters_in_mode_data_module(int(m), int(n), olab_re, olab_im, &
            real(omov, dp), aimag(omov))

        if (get_output_flag_background() > 0) then
            call find_resonance_location(md)
            r_res_local = md%wd%r_res
            call set_resonance_location_in_mode_data_module(r_res_local)
        end if

        call allocate_and_setup_zones(md, sd_ptr, bp_ptr, wd_handle)

        call eval_path_to_linear_data(md%path2project, int(m), int(n), olab, md%path2linear)
        call eval_path_to_dispersion_data(md%path2project, int(m), int(n), olab, md%path2dispersion)
        call eval_path_to_poincare_data(md%path2project, int(m), int(n), olab, md%path2poincare)

        call set_and_make_mode_data_directories(md)

        call copy_mode_paths_to_mode_data_module_native(md)

        handle = mode_data_register(md)
    end function mode_data_create_

    subroutine mode_data_destroy_(handle) bind(C, name="mode_data_destroy_")
        integer(c_intptr_t), value :: handle
        type(mode_data_t), pointer :: md
        integer :: iz
        integer(c_intptr_t) :: wdh

        if (handle == 0_c_intptr_t) return
        call handle_to_mode_data(handle, md)

        wdh = transfer(c_loc(md%wd), wdh)
        call wave_data_destroy(wdh)

        if (allocated(md%r)) deallocate (md%r)
        if (allocated(md%EB)) deallocate (md%EB)
        if (allocated(md%index)) deallocate (md%index)
        if (allocated(md%A)) deallocate (md%A)
        if (allocated(md%B)) deallocate (md%B)
        if (allocated(md%S)) deallocate (md%S)

        if (allocated(md%zone_handles)) then
            do iz = 1, md%Nzones
                if (md%zone_handles(iz) /= 0_c_intptr_t) call zone_destroy_c(md%zone_handles(iz))
            end do
            deallocate (md%zone_handles)
        end if

        deallocate (md)
    end subroutine mode_data_destroy_

    function mode_data_register(md) result(handle)
        type(mode_data_t), pointer, intent(in) :: md
        integer(c_intptr_t) :: handle
        handle = transfer(c_loc(md), handle)
    end function mode_data_register

    subroutine handle_to_mode_data(handle, md)
        integer(c_intptr_t), value :: handle
        type(mode_data_t), pointer, intent(out) :: md
        type(c_ptr) :: cp
        cp = transfer(handle, cp)
        call c_f_pointer(cp, md)
    end subroutine handle_to_mode_data

    !> ---- main driver ----

    subroutine mode_data_calc_all_mode_data_(handle, flag) &
            bind(C, name="mode_data_calc_all_mode_data_")
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: flag
        type(mode_data_t), pointer :: md

        call handle_to_mode_data(handle, md)

        call calc_basis_fields_in_zones(md, int(flag))

        if (flag /= 0) return

        call calc_stitching_equations(md)

        call calc_stitching_equations_determinant(md)

        if (get_output_flag_emfield() > 1) call save_mode_det_data(md)

        call solve_stitching_equations(md)

        call calc_superposition_of_basis_fields_in_zones(md)

        call space_out_fields_in_zones(md)

        call combine_final_wave_fields(md)

        !> setup_wave_fields_for_interpolation: no-op here. The oracle's
        !> EB_int (node-minor layout) only existed to give Neville
        !> interpolation a node-contiguous y-grid; EB(comp,:) array sections
        !> already provide that directly (see module header), so there is
        !> nothing left to set up.

        if (get_output_flag_emfield() > 1) call save_final_wave_fields(md)

        !> calc_and_save_divEB: DEBUG_FLAG is hard-coded 0 in code_settings.h
        !> (KiLCA/code_settings.h), so this call is dead code in the oracle
        !> and is not invoked here either (matches the kilca_flre_zone_m
        !> precedent for the same macro).

        if (get_output_flag_additional() > 0) then
            call calc_quants_in_zones(md)

            if (get_output_flag_additional() > 1) call save_quants_in_zones(md)

            !> combine_final_quants / save_final_quants are no-ops in the
            !> oracle (empty method bodies), so there is nothing to call.
        end if
    end subroutine mode_data_calc_all_mode_data_

    !> ---- accessors for still-C++ callers ----

    function mode_data_get_wd_(handle) result(res) bind(C, name="mode_data_get_wd_")
        integer(c_intptr_t), value :: handle
        integer(c_intptr_t) :: res
        type(mode_data_t), pointer :: md
        call handle_to_mode_data(handle, md)
        res = transfer(c_loc(md%wd), res)
    end function mode_data_get_wd_

    function mode_data_get_zone_handle_(handle, zone_ind) result(res) &
            bind(C, name="mode_data_get_zone_handle_")
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: zone_ind
        integer(c_intptr_t) :: res
        type(mode_data_t), pointer :: md
        call handle_to_mode_data(handle, md)
        res = md%zone_handles(int(zone_ind) + 1)
    end function mode_data_get_zone_handle_

    subroutine mode_data_eval_EB_fields_(handle, x, EB_out) &
            bind(C, name="mode_data_eval_EB_fields_")
        integer(c_intptr_t), value :: handle
        real(c_double), value :: x
        real(c_double), intent(out) :: EB_out(*)
        type(mode_data_t), pointer :: md
        complex(dp) :: EBc(6)
        integer :: comp

        call handle_to_mode_data(handle, md)
        call eval_EB_fields(md, real(x, dp), EBc)
        do comp = 1, 6
            EB_out(2*comp - 1) = real(EBc(comp), dp)
            EB_out(2*comp) = aimag(EBc(comp))
        end do
    end subroutine mode_data_eval_EB_fields_

    subroutine mode_data_eval_diss_power_density_(handle, x, ttype, spec, dpd) &
            bind(C, name="mode_data_eval_diss_power_density_")
        integer(c_intptr_t), value :: handle
        real(c_double), value :: x
        integer(c_int), value :: ttype, spec
        real(c_double), intent(out) :: dpd(*)
        type(mode_data_t), pointer :: md
        class(zone_t), pointer :: z
        integer :: izz

        call handle_to_mode_data(handle, md)
        izz = determine_zone_index_for_point(md, real(x, dp))
        call handle_to_zone(md%zone_handles(izz), z)
        call z%eval_diss_power_density(real(x, dp), int(ttype), int(spec), dpd)
    end subroutine mode_data_eval_diss_power_density_

    subroutine mode_data_eval_current_density_(handle, x, ttype, spec, comp, J) &
            bind(C, name="mode_data_eval_current_density_")
        integer(c_intptr_t), value :: handle
        real(c_double), value :: x
        integer(c_int), value :: ttype, spec, comp
        real(c_double), intent(out) :: J(*)
        type(mode_data_t), pointer :: md
        class(zone_t), pointer :: z
        integer :: izz

        call handle_to_mode_data(handle, md)
        izz = determine_zone_index_for_point(md, real(x, dp))
        call handle_to_zone(md%zone_handles(izz), z)
        call z%eval_current_density(real(x, dp), int(ttype), int(spec), int(comp), J)
    end subroutine mode_data_eval_current_density_

    !> ---- legacy mode_data (mode_m.f90) module bridging ----

    !> Replaces the oracle's copy_mode_paths_to_mode_data_module_ /
    !> copy_mode_paths_from_mode_data_struct_ round trip (the latter
    !> dereferenced a C++ `mode_data **`, which no longer exists): the
    !> legacy mode_data module's path variables are public, so this writes
    !> them directly. Same end state, same point in construction.
    subroutine copy_mode_paths_to_mode_data_module_native(md)
        use mode_data, only: md_path2linear => path2linear, &
            md_path2dispersion => path2dispersion, md_path2poincare => path2poincare
        type(mode_data_t), intent(in) :: md
        md_path2linear = md%path2linear
        md_path2dispersion = md%path2dispersion
        md_path2poincare = md%path2poincare
    end subroutine copy_mode_paths_to_mode_data_module_native

    !> ---- mode_data::save_mode_det_data (mode.cpp) ----

    subroutine save_mode_det_data(md)
        type(mode_data_t), intent(in) :: md
        character(len=1024) :: fname
        integer :: unit, ios

        write (fname, '(a,a,i0,a,i0,a,a,a,a,a)') trim(md%path2linear), &
            'mode_', md%wd%m, '_', md%wd%n, '_[', fmt_e(real(md%wd%olab, dp), 16), ',', &
            fmt_e(aimag(md%wd%olab), 16), '].dat'

        open (newunit=unit, file=trim(fname), status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write (*, '(a,a)') 'Failed to open file ', trim(fname)
        end if

        ! Mirrors the oracle's fprintf format exactly:
        ! "%.16le  %.16le\t%.16le  %.16le\n" (two spaces / tab / two spaces).
        write (unit, '(a,a,a,a,a,a,a)') &
            fmt_e(real(md%wd%olab, dp), 16), '  ', fmt_e(aimag(md%wd%olab), 16), char(9), &
            fmt_e(real(md%wd%det, dp), 16), '  ', fmt_e(aimag(md%wd%det), 16)

        close (unit)
    end subroutine save_mode_det_data

    !> ---- mode_data::find_resonance_location (calc_mode.cpp) ----

    subroutine find_resonance_location(md)
        type(mode_data_t), intent(inout) :: md
        real(dp) :: q_res, r1v, r2v, root
        real(dp), target :: q_res_ctx
        integer(c_int) :: status

        q_res = -real(md%wd%m, dp)/real(md%wd%n, dp)
        r1v = get_background_x0()
        r2v = get_background_xlast()

        if ((eval_q_for_resonance(r1v, c_null_ptr) - q_res)* &
            (eval_q_for_resonance(r2v, c_null_ptr) - q_res) > 0.0_dp) then
            md%wd%r_res = 0.0_dp
            return
        end if

        q_res_ctx = q_res
        status = fortnum_root_brent(c_funloc(qminusq0_cb), r1v, r2v, 0.0_dp, 0.0_dp, &
            100, root, c_loc(q_res_ctx))

        if (status /= FORTNUM_OK) then
            md%wd%r_res = 0.0_dp
            return
        end if

        md%wd%r_res = root
    end subroutine find_resonance_location

    !> q(x) - q_res, the resonance-location root-finding callback. The
    !> background* parameter of the oracle's q() is unused by eval_q
    !> (background_data_m.f90), so c_null_ptr is passed for it, matching the
    !> sentinel-bp precedent already established for zone_t.
    function qminusq0_cb(x, ctx) result(y) bind(C)
        real(c_double), value :: x
        type(c_ptr), value :: ctx
        real(c_double) :: y
        real(dp), pointer :: q_res_p
        call c_f_pointer(ctx, q_res_p)
        y = eval_q_for_resonance(x, c_null_ptr) - q_res_p
    end function qminusq0_cb

    !> ---- path-string helpers (mode.cpp free functions) ----

    subroutine eval_path_to_linear_data(path2project, m, n, olab, path2linear)
        character(len=*), intent(in) :: path2project
        integer, intent(in) :: m, n
        complex(dp), intent(in) :: olab
        character(len=*), intent(out) :: path2linear
        complex(dp) :: flab
        flab = olab/twopi
        write (path2linear, '(a,a,i0,a,i0,a,a,a,a,a)') trim(path2project), &
            'linear-data/m_', m, '_n_', n, '_flab_[', fmt_g(real(flab, dp), 15), ',', &
            fmt_g(aimag(flab), 15), ']/'
    end subroutine eval_path_to_linear_data

    subroutine eval_path_to_dispersion_data(path2project, m, n, olab, path2dispersion)
        character(len=*), intent(in) :: path2project
        integer, intent(in) :: m, n
        complex(dp), intent(in) :: olab
        character(len=*), intent(out) :: path2dispersion
        complex(dp) :: flab
        flab = olab/twopi
        write (path2dispersion, '(a,a,i0,a,i0,a,a,a,a,a)') trim(path2project), &
            'dispersion-data/m_', m, '_n_', n, '_flab_[', fmt_g(real(flab, dp), 15), ',', &
            fmt_g(aimag(flab), 15), ']/'
    end subroutine eval_path_to_dispersion_data

    subroutine eval_path_to_poincare_data(path2project, m, n, olab, path2poincare)
        character(len=*), intent(in) :: path2project
        integer, intent(in) :: m, n
        complex(dp), intent(in) :: olab
        character(len=*), intent(out) :: path2poincare
        complex(dp) :: flab
        flab = olab/twopi
        write (path2poincare, '(a,a,i0,a,i0,a,a,a,a,a)') trim(path2project), &
            'poincare-data/m_', m, '_n_', n, '_flab_[', fmt_g(real(flab, dp), 15), ',', &
            fmt_g(aimag(flab), 15), ']/'
    end subroutine eval_path_to_poincare_data

    !> ---- mode_data::set_and_make_mode_data_directories (mode.cpp) ----

    subroutine set_and_make_mode_data_directories(md)
        type(mode_data_t), intent(in) :: md
        character(len=1024) :: sys_command, fname
        integer :: unit, ierr_unused

        if (get_output_flag_emfield() > 1) then
            sys_command = 'mkdir -p '//trim(md%path2linear)
            call execute_command_line(trim(sys_command), exitstat=ierr_unused)

            sys_command = 'mkdir -p '//trim(md%path2linear)//'debug-data/'
            call execute_command_line(trim(sys_command), exitstat=ierr_unused)

            fname = trim(md%path2linear)//'mode_data.dat'
            open (newunit=unit, file=trim(fname), status='replace', action='write')
            write (unit, '(a)') '%m n'
            write (unit, '(a)') '%Re(flab) Im(flab)'
            write (unit, '(a)') '%Re(fmov) Im(fmov)'
            write (unit, '(a)') '%r_res 0.0'
            write (unit, '(i0,1x,i0)') md%wd%m, md%wd%n
            write (unit, '(a,a,a)') fmt_g(real(md%wd%olab, dp)/twopi, 16), ' ', &
                fmt_g(aimag(md%wd%olab)/twopi, 16)
            write (unit, '(a,a,a)') fmt_g(real(md%wd%omov, dp)/twopi, 16), ' ', &
                fmt_g(aimag(md%wd%omov)/twopi, 16)
            write (unit, '(a,a,a)') fmt_g(md%wd%r_res, 16), ' ', fmt_g(0.0_dp, 16)
            close (unit)
        end if

        if (get_output_flag_dispersion() > 1) then
            sys_command = 'mkdir -p '//trim(md%path2dispersion)
            call execute_command_line(trim(sys_command), exitstat=ierr_unused)
        end if
    end subroutine set_and_make_mode_data_directories

    !> ---- mode_data::allocate_and_setup_zones (calc_mode.cpp) ----

    subroutine allocate_and_setup_zones(md, sd_ptr, bp_ptr, wd_handle)
        type(mode_data_t), intent(inout) :: md
        integer(c_intptr_t), intent(in) :: sd_ptr, bp_ptr, wd_handle
        integer :: k, ztype
        character(len=1024) :: filename
        character(kind=c_char), allocatable :: cpath(:)
        class(zone_t), pointer :: z

        md%Nzones = determine_number_of_zones(md%path2project)

        if (md%Nzones < 2) then
            write (*, '(a,i0)') 'allocate_and_setup_zones: Nzones must be >= 2: ', md%Nzones
            stop 1
        end if

        if (allocated(md%zone_handles)) deallocate (md%zone_handles)
        allocate (md%zone_handles(md%Nzones))

        cpath = to_cstr(md%path2project)

        do k = 0, md%Nzones - 1
            filename = get_zone_file_name(md%path2project, k)
            ztype = determine_zone_type(filename)

            if (ztype > 1 .and. get_output_flag_background() == 0) then
                write (*, '(a)') 'error: allocate_and_setup_zones: background is not set!'
                stop 1
            end if

            select case (ztype)
            case (PLASMA_MODEL_VACUUM, PLASMA_MODEL_MEDIUM)
                md%zone_handles(k + 1) = hmedium_zone_create(sd_ptr, bp_ptr, wd_handle, &
                    cpath, int(k, c_int))
            case (PLASMA_MODEL_IMHD)
                md%zone_handles(k + 1) = imhd_zone_create(sd_ptr, bp_ptr, wd_handle, &
                    cpath, int(k, c_int))
            case (PLASMA_MODEL_RMHD)
                write (*, '(a,i0,a)') 'The plasma model for zone ', k, ' is not implemented.'
                stop 1
            case (PLASMA_MODEL_FLRE)
                md%zone_handles(k + 1) = flre_zone_create_(sd_ptr, bp_ptr, wd_handle, &
                    cpath, int(k, c_int))
            case default
                write (*, '(a,i0,a)') 'The plasma model for zone ', k, ' is unknown.'
                stop 1
            end select

            call handle_to_zone(md%zone_handles(k + 1), z)
            call z%read_settings(trim(filename))
        end do

        call check_zones_parameters(md)
    end subroutine allocate_and_setup_zones

    subroutine check_zones_parameters(md)
        type(mode_data_t), intent(in) :: md
        integer :: k, iz
        class(zone_t), pointer :: z1, z2, zends

        do k = 1, md%Nzones - 1
            call handle_to_zone(md%zone_handles(k), z1)
            call handle_to_zone(md%zone_handles(k + 1), z2)

            if (z1%get_r2() /= z2%get_r1() .or. z1%bc2 /= z2%bc1) then
                write (*, '(a,i0)') 'check_zones: boundaries are different: k = ', k - 1
                stop 1
            end if

            if (.not. (z1%bc2 == BOUNDARY_INTERFACE .or. z1%bc2 == BOUNDARY_ANTENNA)) then
                write (*, '(a,i0,a,i0)') &
                    'check_zones: improper type of the right boundary for zone ', k - 1, &
                    ': ', z1%bc2
                stop 1
            end if
        end do

        iz = 1
        call handle_to_zone(md%zone_handles(iz), zends)
        if (.not. (zends%bc1 == BOUNDARY_CENTER .or. zends%bc1 == BOUNDARY_IDEALWALL)) then
            write (*, '(a,i0)') 'check_zones: improper type of the first boundary: ', zends%bc1
            stop 1
        end if

        iz = md%Nzones
        call handle_to_zone(md%zone_handles(iz), zends)
        if (.not. (zends%bc2 == BOUNDARY_INFINITY .or. zends%bc2 == BOUNDARY_IDEALWALL)) then
            write (*, '(a,i0)') 'check_zones: improper type of the last boundary: ', zends%bc2
            stop 1
        end if
    end subroutine check_zones_parameters

    !> ---- zone settings-file discovery (calc_mode.cpp statics) ----

    function determine_number_of_zones(path2project) result(nz)
        character(len=*), intent(in) :: path2project
        integer :: nz
        character(kind=c_char), allocatable :: cpath(:), cpattern(:)
        type(c_ptr) :: dirp, ep_cptr
        type(dirent_t), pointer :: ep
        integer(c_int) :: ret

        cpath = to_cstr(path2project)
        cpattern = to_cstr('zone_*.in')
        nz = 0

        dirp = opendir_c(cpath)
        if (.not. c_associated(dirp)) then
            write (*, '(a,a)') &
                'determine_number_of_zones: faled to open the project directory ', &
                trim(path2project)
            stop 1
        end if

        do
            ep_cptr = readdir_c(dirp)
            if (.not. c_associated(ep_cptr)) exit
            call c_f_pointer(ep_cptr, ep)
            if (fnmatch_c(cpattern, ep%d_name, 0_c_int) == 0) nz = nz + 1
        end do

        ret = closedir_c(dirp)
    end function determine_number_of_zones

    function get_zone_file_name(path2project, zone_index) result(fname)
        character(len=*), intent(in) :: path2project
        integer, intent(in) :: zone_index
        character(len=1024) :: fname
        character(kind=c_char), allocatable :: cpath(:), cpattern(:)
        character(len=32) :: pattern_str
        character(len=256) :: dname
        type(c_ptr) :: dirp, ep_cptr
        type(dirent_t), pointer :: ep
        integer :: found
        integer(c_int) :: ret

        cpath = to_cstr(path2project)
        write (pattern_str, '(a,i0,a)') '*zone_', zone_index + 1, '*.in'
        cpattern = to_cstr(trim(pattern_str))

        fname = ''
        found = 0

        dirp = opendir_c(cpath)
        if (.not. c_associated(dirp)) then
            write (*, '(a,a)') 'get_zone_file_name: faled to open the project directory ', &
                trim(path2project)
            stop 1
        end if

        do
            ep_cptr = readdir_c(dirp)
            if (.not. c_associated(ep_cptr)) exit
            call c_f_pointer(ep_cptr, ep)
            if (fnmatch_c(cpattern, ep%d_name, 0_c_int) == 0) then
                call dirent_name_to_fortran(ep%d_name, dname)
                fname = trim(path2project)//trim(dname)
                found = 1
                exit
            end if
        end do

        ret = closedir_c(dirp)

        if (found == 0) then
            write (*, '(a,i0,a)') &
                'get_zone_file_name: failed to find the file name for zone ', zone_index, '!'
            stop 1
        end if
    end function get_zone_file_name

    function determine_zone_type(file) result(medium)
        character(len=*), intent(in) :: file
        integer :: medium
        integer :: unit, ios, k
        real(dp) :: dummy_r1
        character(len=64) :: bstr1, mstr

        open (newunit=unit, file=trim(file), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write (*, '(a,a)') 'error: determine_zone_type: failed to open file ', trim(file)
            stop 0
        end if

        call skip_line(unit)
        call read_real_before_hash(unit, dummy_r1)
        call read_token_before_hash(unit, bstr1)
        call read_token_before_hash(unit, mstr)

        close (unit)

        medium = -1
        do k = 0, Nmed - 1
            if (trim(mstr) == trim(med_str(k))) then
                medium = k
                exit
            end if
        end do

        if (medium == -1) then
            write (*, '(a,a)') 'error: determine_zone_type: medium type is unknown: ', trim(mstr)
            stop 1
        end if
    end function determine_zone_type

    subroutine dirent_name_to_fortran(cname, fname)
        character(kind=c_char), intent(in) :: cname(:)
        character(len=*), intent(out) :: fname
        integer :: i
        fname = ''
        do i = 1, min(size(cname), len(fname))
            if (cname(i) == c_null_char) exit
            fname(i:i) = cname(i)
        end do
    end subroutine dirent_name_to_fortran

    function to_cstr(s) result(c)
        character(len=*), intent(in) :: s
        character(kind=c_char), allocatable :: c(:)
        integer :: i, n
        n = len_trim(s)
        allocate (c(n + 1))
        do i = 1, n
            c(i) = s(i:i)
        end do
        c(n + 1) = c_null_char
    end function to_cstr

    function c_string_to_fortran_local(cstr) result(fstr)
        character(kind=c_char), intent(in) :: cstr(*)
        character(len=1024) :: fstr
        integer :: i
        fstr = ''
        do i = 1, 1024
            if (cstr(i) == c_null_char) exit
            fstr(i:i) = cstr(i)
        end do
    end function c_string_to_fortran_local

    !> ---- mode_data::calc_basis_fields_in_zones (calc_mode.cpp) ----

    subroutine calc_basis_fields_in_zones(md, flag)
        type(mode_data_t), intent(inout) :: md
        integer, intent(in) :: flag
        integer :: iz
        class(zone_t), pointer :: z
        do iz = 1, md%Nzones
            call handle_to_zone(md%zone_handles(iz), z)
            call z%calc_basis_fields(flag)
        end do
    end subroutine calc_basis_fields_in_zones

    !> ---- mode_data::calc_stitching_equations (calc_mode.cpp) ----

    subroutine calc_stitching_equations(md)
        type(mode_data_t), intent(inout) :: md
        integer :: iz, neq, ieq, nvar, ivar, Nw, len_b, Nw1, Nw2, len1, len2, flg_ant
        class(zone_t), pointer :: z, z1, z2
        complex(dp), allocatable :: M(:, :), J(:)

        md%Nc = 0
        do iz = 1, md%Nzones
            call handle_to_zone(md%zone_handles(iz), z)
            md%Nc = md%Nc + z%get_dim_of_basis()
        end do

        if (allocated(md%A)) deallocate (md%A)
        allocate (md%A(md%Nc, md%Nc))
        if (allocated(md%B)) deallocate (md%B)
        allocate (md%B(md%Nc))

        allocate (M(md%Nc, md%Nc))
        allocate (J(md%Nc))

        neq = 1; ieq = 0; nvar = 1; ivar = 0

        ! zeroes A, B (the legacy routine's ieq == 0 branch is the zeroing call)
        call update_system_matrix_and_rhs_vector(md%Nc, md%A, md%B, neq, ieq, nvar, ivar, M, J)

        ieq = 1
        ivar = 1

        ! first boundary (usually, but not necessarily, a center)
        iz = 1
        call handle_to_zone(md%zone_handles(iz), z)
        Nw = z%get_dim_of_basis()
        len_b = z%get_dim_of_basis_vector()

        if (z%medium == PLASMA_MODEL_VACUUM .or. z%medium == PLASMA_MODEL_MEDIUM) then
            select case (z%bc1)
            case (BOUNDARY_CENTER)
                call center_equations_hommed(Nw, len_b, md%zone_handles(iz), z%basis(:, :, 1), &
                    neq, nvar, M, J)
            case (BOUNDARY_IDEALWALL)
                call ideal_wall_equations_hommed(Nw, len_b, md%zone_handles(iz), &
                    z%basis(:, :, 1), neq, nvar, M, J)
            case default
                write (*, '(a)') 'error: calc_stitching_equations: first boundary is unknown.'
                stop 1
            end select
        else if (z%medium == PLASMA_MODEL_IMHD) then
            select case (z%bc1)
            case (BOUNDARY_CENTER)
                call center_equations_imhd(Nw, len_b, md%zone_handles(iz), z%basis(:, :, 1), &
                    neq, nvar, M, J)
            case (BOUNDARY_IDEALWALL)
                call ideal_wall_equations_imhd(Nw, len_b, md%zone_handles(iz), &
                    z%basis(:, :, 1), neq, nvar, M, J)
            case default
                write (*, '(a)') 'error: calc_stitching_equations: first boundary is unknown.'
                stop 1
            end select
        else if (z%medium == PLASMA_MODEL_RMHD) then
            write (*, '(a)') 'error: calc_stitching_equations: not implemented.'
            stop 1
        else if (z%medium == PLASMA_MODEL_FLRE) then
            select case (z%bc1)
            case (BOUNDARY_CENTER)
                call center_equations_flre(Nw, len_b, md%zone_handles(iz), z%basis(:, :, 1), &
                    neq, nvar, M, J)
            case (BOUNDARY_IDEALWALL)
                call ideal_wall_equations_flre(Nw, len_b, md%zone_handles(iz), &
                    z%basis(:, :, 1), neq, nvar, M, J)
            case default
                write (*, '(a)') 'error: calc_stitching_equations: first boundary is unknown.'
                stop 1
            end select
        else
            write (*, '(a)') 'error: calc_stitching_equations: unknown plasma model.'
            stop 1
        end if

        call update_system_matrix_and_rhs_vector(md%Nc, md%A, md%B, neq, ieq, nvar, ivar, M, J)

        ieq = ieq + neq
        ! ivar unchanged: same zone & variables (oracle's `ivar += 0`)

        ! internal boundaries
        do iz = 1, md%Nzones - 1
            call handle_to_zone(md%zone_handles(iz), z1)
            call handle_to_zone(md%zone_handles(iz + 1), z2)

            Nw1 = z1%get_dim_of_basis()
            Nw2 = z2%get_dim_of_basis()
            len1 = z1%get_dim_of_basis_vector()
            len2 = z2%get_dim_of_basis_vector()

            if (z1%bc2 == BOUNDARY_ANTENNA) then
                flg_ant = 1
            else
                flg_ant = 0
            end if

            if ((z1%medium == PLASMA_MODEL_VACUUM .or. z1%medium == PLASMA_MODEL_MEDIUM) .and. &
                (z2%medium == PLASMA_MODEL_VACUUM .or. z2%medium == PLASMA_MODEL_MEDIUM)) then
                call stitching_equations_hommed_hommed(Nw1, len1, md%zone_handles(iz), &
                    z1%basis(:, :, z1%dim), Nw2, len2, md%zone_handles(iz + 1), &
                    z2%basis(:, :, 1), flg_ant, neq, nvar, M, J)
            else if (z1%medium == PLASMA_MODEL_IMHD .and. &
                    (z2%medium == PLASMA_MODEL_VACUUM .or. z2%medium == PLASMA_MODEL_MEDIUM)) then
                call stitching_equations_imhd_hommed(Nw1, len1, md%zone_handles(iz), &
                    z1%basis(:, :, z1%dim), Nw2, len2, md%zone_handles(iz + 1), &
                    z2%basis(:, :, 1), flg_ant, neq, nvar, M, J)
            else if (z1%medium == PLASMA_MODEL_IMHD .and. z2%medium == PLASMA_MODEL_IMHD) then
                call stitching_equations_imhd_imhd(Nw1, len1, md%zone_handles(iz), &
                    z1%basis(:, :, z1%dim), Nw2, len2, md%zone_handles(iz + 1), &
                    z2%basis(:, :, 1), flg_ant, neq, nvar, M, J)
            else if (z1%medium == PLASMA_MODEL_RMHD .and. &
                    (z2%medium == PLASMA_MODEL_VACUUM .or. z2%medium == PLASMA_MODEL_MEDIUM)) then
                write (*, '(a)') 'error: calc_stitching_equations: not implemented.'
                stop 1
            else if (z1%medium == PLASMA_MODEL_RMHD .and. z2%medium == PLASMA_MODEL_RMHD) then
                write (*, '(a)') 'error: calc_stitching_equations: not implemented.'
                stop 1
            else if (z1%medium == PLASMA_MODEL_FLRE .and. &
                    (z2%medium == PLASMA_MODEL_VACUUM .or. z2%medium == PLASMA_MODEL_MEDIUM)) then
                call stitching_equations_flre_hommed(Nw1, len1, md%zone_handles(iz), &
                    z1%basis(:, :, z1%dim), Nw2, len2, md%zone_handles(iz + 1), &
                    z2%basis(:, :, 1), flg_ant, neq, nvar, M, J)
            else if (z1%medium == PLASMA_MODEL_FLRE .and. z2%medium == PLASMA_MODEL_FLRE) then
                call stitching_equations_flre_flre(Nw1, len1, md%zone_handles(iz), &
                    z1%basis(:, :, z1%dim), Nw2, len2, md%zone_handles(iz + 1), &
                    z2%basis(:, :, 1), flg_ant, neq, nvar, M, J)
            else
                write (*, '(a)') 'error: calc_stitching_equations: unknown type of boundary.'
                stop 1
            end if

            call update_system_matrix_and_rhs_vector(md%Nc, md%A, md%B, neq, ieq, nvar, ivar, M, J)

            ieq = ieq + neq
            ivar = ivar + Nw1
        end do

        ! last boundary (usually, but not necessarily, infinity)
        iz = md%Nzones
        call handle_to_zone(md%zone_handles(iz), z)
        Nw = z%get_dim_of_basis()
        len_b = z%get_dim_of_basis_vector()

        if (z%medium == PLASMA_MODEL_VACUUM .or. z%medium == PLASMA_MODEL_MEDIUM) then
            select case (z%bc2)
            case (BOUNDARY_INFINITY)
                call infinity_equations_hommed(Nw, len_b, md%zone_handles(iz), &
                    z%basis(:, :, z%dim), neq, nvar, M, J)
            case (BOUNDARY_IDEALWALL)
                call ideal_wall_equations_hommed(Nw, len_b, md%zone_handles(iz), &
                    z%basis(:, :, z%dim), neq, nvar, M, J)
            case default
                write (*, '(a)') 'error: calc_stitching_equations: last boundary is unknown.'
                stop 1
            end select
        else if (z%medium == PLASMA_MODEL_IMHD) then
            select case (z%bc2)
            case (BOUNDARY_INFINITY)
                call infinity_equations_imhd(Nw, len_b, md%zone_handles(iz), &
                    z%basis(:, :, z%dim), neq, nvar, M, J)
            case (BOUNDARY_IDEALWALL)
                call ideal_wall_equations_imhd(Nw, len_b, md%zone_handles(iz), &
                    z%basis(:, :, z%dim), neq, nvar, M, J)
            case default
                write (*, '(a)') 'error: calc_stitching_equations: last boundary is unknown.'
                stop 1
            end select
        else if (z%medium == PLASMA_MODEL_RMHD) then
            write (*, '(a)') 'error: calc_stitching_equations: not implemented.'
            stop 1
        else if (z%medium == PLASMA_MODEL_FLRE) then
            select case (z%bc2)
            case (BOUNDARY_INFINITY)
                call infinity_equations_flre(Nw, len_b, md%zone_handles(iz), &
                    z%basis(:, :, z%dim), neq, nvar, M, J)
            case (BOUNDARY_IDEALWALL)
                call ideal_wall_equations_flre(Nw, len_b, md%zone_handles(iz), &
                    z%basis(:, :, z%dim), neq, nvar, M, J)
            case default
                write (*, '(a)') 'error: calc_stitching_equations: last boundary is unknown.'
                stop 1
            end select
        else
            write (*, '(a)') 'error: calc_stitching_equations: unknown plasma model.'
            stop 1
        end if

        call update_system_matrix_and_rhs_vector(md%Nc, md%A, md%B, neq, ieq, nvar, ivar, M, J)

        ieq = ieq + neq
        ivar = ivar + Nw

        if (ieq - 1 /= md%Nc .or. ivar - 1 /= md%Nc) then
            write (*, '(a,i0,a,i0)') &
                'error: calc_stitching_equations: wrong final indices: ieq = ', ieq, &
                ', ivar =', ivar
            stop 1
        end if

        deallocate (M, J)
    end subroutine calc_stitching_equations

    subroutine calc_stitching_equations_determinant(md)
        type(mode_data_t), intent(inout) :: md
        complex(dp) :: det
        call calc_system_determinant(md%Nc, md%A, det)
        md%wd%det = det
    end subroutine calc_stitching_equations_determinant

    subroutine solve_stitching_equations(md)
        type(mode_data_t), intent(inout) :: md
        if (allocated(md%S)) deallocate (md%S)
        allocate (md%S(md%Nc))
        call find_superposition_coeffs(md%Nc, md%A, md%B, md%S)
    end subroutine solve_stitching_equations

    subroutine calc_superposition_of_basis_fields_in_zones(md)
        type(mode_data_t), intent(inout) :: md
        integer :: iz, ind, nw, i
        class(zone_t), pointer :: z
        real(dp), allocatable :: S_p(:)

        ind = 0
        do iz = 1, md%Nzones
            call handle_to_zone(md%zone_handles(iz), z)
            nw = z%get_dim_of_basis()
            allocate (S_p(2*nw))
            do i = 1, nw
                S_p(2*i - 1) = real(md%S(ind + i), dp)
                S_p(2*i) = aimag(md%S(ind + i))
            end do
            call z%calc_superposition_of_basis_fields(S_p)
            deallocate (S_p)
            ind = ind + nw
        end do
    end subroutine calc_superposition_of_basis_fields_in_zones

    subroutine space_out_fields_in_zones(md)
        type(mode_data_t), intent(inout) :: md
        integer :: iz
        class(zone_t), pointer :: z
        do iz = 1, md%Nzones
            call handle_to_zone(md%zone_handles(iz), z)
            call z%calc_final_fields()
            !> zone_save_final_fields_ here is gated behind DEBUG_FLAG, hard
            !> -coded 0 (KiLCA/code_settings.h); dead in the oracle, omitted.
        end do
    end subroutine space_out_fields_in_zones

    subroutine combine_final_wave_fields(md)
        type(mode_data_t), intent(inout) :: md
        integer :: iz, ind, zdim, node, comp
        class(zone_t), pointer :: z
        real(dp), allocatable :: ebflat(:)

        md%dim = 0
        do iz = 1, md%Nzones
            call handle_to_zone(md%zone_handles(iz), z)
            md%dim = md%dim + z%get_radial_grid_dimension()
        end do

        if (allocated(md%r)) deallocate (md%r)
        allocate (md%r(md%dim))
        if (allocated(md%EB)) deallocate (md%EB)
        allocate (md%EB(6, md%dim))
        if (allocated(md%index)) deallocate (md%index)
        allocate (md%index(md%Nzones))

        ind = 0
        do iz = 1, md%Nzones
            call handle_to_zone(md%zone_handles(iz), z)
            md%index(iz) = ind
            zdim = z%get_radial_grid_dimension()

            call z%copy_radial_grid(md%r(ind + 1:ind + zdim))

            allocate (ebflat(12*zdim))
            call z%copy_E_and_B_fields(ebflat)
            do node = 0, zdim - 1
                do comp = 0, 5
                    md%EB(comp + 1, ind + node + 1) = cmplx(ebflat(2*(comp + 6*node) + 1), &
                        ebflat(2*(comp + 6*node) + 2), dp)
                end do
            end do
            deallocate (ebflat)

            ind = ind + zdim
        end do
    end subroutine combine_final_wave_fields

    subroutine save_final_wave_fields(md)
        type(mode_data_t), intent(in) :: md
        character(len=1024) :: fname
        integer :: unit, i, comp

        fname = trim(md%path2linear)//'EB.dat'
        open (newunit=unit, file=trim(fname), status='replace', action='write')

        do i = 1, md%dim
            write (unit, '(es24.16e3)', advance='no') md%r(i)
            do comp = 1, 6
                write (unit, '(a,es24.16e3,a,es24.16e3)', advance='no') &
                    char(9), real(md%EB(comp, i), dp), char(9), aimag(md%EB(comp, i))
            end do
            write (unit, *)
        end do

        close (unit)
    end subroutine save_final_wave_fields

    subroutine calc_quants_in_zones(md)
        type(mode_data_t), intent(inout) :: md
        integer :: iz
        class(zone_t), pointer :: z
        do iz = 1, md%Nzones
            call handle_to_zone(md%zone_handles(iz), z)
            call z%calc_all_quants()
        end do
    end subroutine calc_quants_in_zones

    subroutine save_quants_in_zones(md)
        type(mode_data_t), intent(inout) :: md
        integer :: iz
        class(zone_t), pointer :: z
        do iz = 1, md%Nzones
            call handle_to_zone(md%zone_handles(iz), z)
            call z%save_all_quants()
        end do
    end subroutine save_quants_in_zones

    !> ---- point lookups (mode_data::determine_zone_index_for_point etc) ----

    function determine_zone_index_for_point(md, x) result(izz)
        type(mode_data_t), intent(in) :: md
        real(dp), intent(in) :: x
        integer :: izz
        integer :: iz
        class(zone_t), pointer :: z
        logical :: hit

        izz = -1
        do iz = 1, md%Nzones
            call handle_to_zone(md%zone_handles(iz), z)
            hit = x >= z%get_r1() .and. x <= z%get_r2()
            if (hit) then
                izz = iz
                exit
            end if
        end do

        if (izz == -1) then
            write (*, '(a,es12.4)') &
                'warning: determine_zone_index_for_point: x is outside the grid: x = ', x
            izz = 1
        end if
    end function determine_zone_index_for_point

    subroutine eval_EB_fields(md, x, EB_out)
        type(mode_data_t), intent(in) :: md
        real(dp), intent(in) :: x
        complex(dp), intent(out) :: EB_out(6)
        integer :: izz, D, comp, start
        integer(c_int) :: neville_ind
        class(zone_t), pointer :: z
        real(dp) :: re_v(0:0), im_v(0:0)

        izz = determine_zone_index_for_point(md, x)
        call handle_to_zone(md%zone_handles(izz), z)
        D = z%get_radial_grid_dimension()
        start = md%index(izz)

        neville_ind = 0
        do comp = 1, 6
            call eval_neville_polynom(D, md%r(start + 1:start + D), &
                real(md%EB(comp, start + 1:start + D), dp), 5, x, 0, 0, neville_ind, re_v)
            call eval_neville_polynom(D, md%r(start + 1:start + D), &
                aimag(md%EB(comp, start + 1:start + D)), 5, x, 0, 0, neville_ind, im_v)
            EB_out(comp) = cmplx(re_v(0), im_v(0), dp)
        end do
    end subroutine eval_EB_fields

    !> ---- divEB / calc_and_save_divEB: translated, never invoked (dead in
    !> the oracle, gated behind DEBUG_FLAG hard-coded 0) ----

    subroutine divEB(md, x, div)
        type(mode_data_t), intent(in) :: md
        real(dp), intent(in) :: x
        complex(dp), intent(out) :: div(2)
        integer :: izz, D, comp, start
        integer(c_int) :: neville_ind
        class(zone_t), pointer :: z
        real(dp) :: re_v(0:1), im_v(0:1)
        complex(dp) :: EBc(6), dEBc(6)
        real(dp) :: kt, kz

        izz = determine_zone_index_for_point(md, x)
        call handle_to_zone(md%zone_handles(izz), z)
        D = z%get_radial_grid_dimension()
        start = md%index(izz)

        neville_ind = 0
        do comp = 1, 6
            call eval_neville_polynom(D, md%r(start + 1:start + D), &
                real(md%EB(comp, start + 1:start + D), dp), 5, x, 0, 1, neville_ind, re_v)
            call eval_neville_polynom(D, md%r(start + 1:start + D), &
                aimag(md%EB(comp, start + 1:start + D)), 5, x, 0, 1, neville_ind, im_v)
            EBc(comp) = cmplx(re_v(0), im_v(0), dp)
            dEBc(comp) = cmplx(re_v(1), im_v(1), dp)
        end do

        kt = real(md%wd%m, dp)/x
        kz = real(md%wd%n, dp)/get_background_rtor()

        div(1) = EBc(1)/x + dEBc(1) + cmplx(0.0_dp, 1.0_dp, dp)*kt*EBc(2) + &
            cmplx(0.0_dp, 1.0_dp, dp)*kz*EBc(3)
        div(2) = EBc(4)/x + dEBc(4) + cmplx(0.0_dp, 1.0_dp, dp)*kt*EBc(5) + &
            cmplx(0.0_dp, 1.0_dp, dp)*kz*EBc(6)
    end subroutine divEB

    subroutine calc_and_save_divEB(md)
        type(mode_data_t), intent(in) :: md
        character(len=1024) :: fname
        integer :: unit, i
        complex(dp) :: div(2)

        fname = trim(md%path2linear)//'divEB.dat'
        open (newunit=unit, file=trim(fname), status='replace', action='write')

        do i = 1, md%dim
            call divEB(md, md%r(i), div)
            write (unit, '(es24.16e3,a,es24.16e3,a,es24.16e3,a,es24.16e3,a,es24.16e3)') &
                md%r(i), char(9), real(div(1), dp), char(9), aimag(div(1)), char(9), &
                real(div(2), dp), char(9), aimag(div(2))
        end do

        close (unit)
    end subroutine calc_and_save_divEB

end module kilca_mode_data_m

!> Plain (non-module) Fortran subprogram, deliberately kept OUTSIDE
!> kilca_mode_data_m: its only caller, mode_m.f90's own legacy
!> copy_mode_paths_to_mode_data_module, is itself Fortran and reaches it
!> via an implicit external interface (plain `name_` mangling), not a
!> bind(C) name or a module-procedure `__module_MOD_name` mangling - so
!> this must be a free subprogram, not something `contains`ed in a module
!> (a module-contained version compiles fine but produces the wrong
!> symbol name and leaves the real `name_` symbol undefined at link time).
!> That caller (copy_mode_paths_to_mode_data_module) itself has zero
!> callers anywhere in this codebase now (kilca_mode_data_m's own
!> copy_mode_paths_to_mode_data_module_native replaced it for the only
!> call site that mattered), so this body never executes - but gfortran
!> still compiles mode_m.f90's dead subroutine, and the linker still needs
!> this external symbol resolved (it used to be provided by mode.cpp, now
!> deleted). Implemented for real rather than stubbed, since the cost is
!> the same either way and a real implementation can't silently drift
!> from intent.
subroutine copy_mode_paths_from_mode_data_struct(md, path2linear, path2dispersion, &
    path2poincare)
    use, intrinsic :: iso_c_binding, only: c_intptr_t, c_ptr, c_f_pointer
    use kilca_mode_data_m, only: mode_data_t
    implicit none
    integer(c_intptr_t), intent(in) :: md
    character(len=*), intent(out) :: path2linear, path2dispersion, path2poincare
    type(mode_data_t), pointer :: mdp
    type(c_ptr) :: cp
    cp = transfer(md, cp)
    call c_f_pointer(cp, mdp)
    path2linear = mdp%path2linear
    path2dispersion = mdp%path2dispersion
    path2poincare = mdp%path2poincare
end subroutine copy_mode_paths_from_mode_data_struct
