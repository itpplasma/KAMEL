!> The top-level KiLCA orchestrator, formerly the C++ core_data class
!> (core.{h,cpp}). Per-instance handle (same c_loc/transfer pattern as
!> cond_profiles etc - core_data is genuinely instantiated multiple times,
!> once per wave_code_interface_64bit.f90/QL-Balance call).
!>
!> eigmode/calc_eigmode.cpp and eigmode/find_eigmodes.cpp (the zersol-based
!> complex-root search, still C++, a separate and substantially larger
!> translation unit) call back into this module's getters with a handle
!> instead of the oracle's `core_data *this` - only their call SIGNATURES
!> were updated for this step, not their own internal root-finding logic.
module kilca_core_data_m
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_double, c_char, &
        c_ptr, c_loc, c_f_pointer, c_null_char
    use constants, only: dp, pi
    use kilca_settings_m, only: settings_create_, settings_read_settings_, &
        settings_get_path2project_
    use kilca_mode_data_m, only: mode_data_create_, mode_data_destroy_, &
        mode_data_calc_all_mode_data_
    implicit none
    private

    public :: core_data_create_, core_data_destroy_
    public :: core_data_calc_and_set_mode_independent_
    public :: core_data_calc_and_set_mode_dependent_antenna_
    public :: core_data_calc_and_set_mode_dependent_eigmode_
    public :: core_data_calc_and_set_mode_dependent_antenna_interface_
    public :: core_data_calc_and_set_mode_dependent_antenna_interface_mn_
    public :: core_data_get_dim_, core_data_get_mda_element_, core_data_set_mda_element_
    public :: core_data_get_bp_, core_data_get_sd_
    public :: core_data_get_path2project_

    type :: core_data_t
        character(len=1024) :: path2project = ''
        integer(c_intptr_t) :: sd = 0
        integer(c_intptr_t) :: bp = 1_c_intptr_t
        integer :: dim = 0
        integer(c_intptr_t), allocatable :: mda(:)
    end type core_data_t

    !> Static settings caches (one per project-type substring match),
    !> mirroring core_data::calc_and_set_mode_independent_core_data's own
    !> static_settings_vacuum/static_settings_flre - settings reading has
    !> always been effectively global (see kilca_settings_m), so this
    !> caching avoids re-parsing the same .in files across core_data
    !> instances for the same project type.
    integer(c_intptr_t) :: static_settings_vacuum = 0
    integer(c_intptr_t) :: static_settings_flre = 0

    interface
        subroutine get_pointer_precision(ppp) bind(C, name="get_pointer_precision_")
            import :: c_int
            integer(c_int), intent(out) :: ppp
        end subroutine get_pointer_precision

        function background_create(path2project_p) bind(C, name="background_create_") result(handle)
            import :: c_ptr, c_intptr_t
            type(c_ptr), value :: path2project_p
            integer(c_intptr_t) :: handle
        end function background_create

        integer(c_int) function get_background_calc_back() bind(C, name="get_background_calc_back_")
            import :: c_int
        end function get_background_calc_back

        subroutine background_set_profiles_from_files() &
            bind(C, name="background_set_profiles_from_files_")
        end subroutine background_set_profiles_from_files

        subroutine background_set_profiles_from_interface() &
            bind(C, name="background_set_profiles_from_interface_")
        end subroutine background_set_profiles_from_interface

        integer(c_int) function get_antenna_dma() bind(C, name="get_antenna_dma_")
            import :: c_int
        end function get_antenna_dma

        subroutine get_antenna_flab(re, im) bind(C, name="get_antenna_flab_")
            import :: c_double
            real(c_double), intent(out) :: re, im
        end subroutine get_antenna_flab

        subroutine get_antenna_mode(ind, m, n) bind(C, name="get_antenna_mode_")
            import :: c_int
            integer(c_int), value :: ind
            integer(c_int), intent(out) :: m, n
        end subroutine get_antenna_mode

        integer(c_int) function get_eigmode_search_flag() bind(C, name="get_eigmode_search_flag_")
            import :: c_int
        end function get_eigmode_search_flag

        subroutine clear_all_data_in_mode_data_module() &
            bind(C, name="clear_all_data_in_mode_data_module_")
        end subroutine clear_all_data_in_mode_data_module

        !> Pre-existing legacy Fortran (core_m.f90), stores the settings
        !> handle into the `core` module's sd_ptr - by-reference, matching
        !> its `integer(pp), intent(in) :: sd` dummy.
        subroutine set_settings_in_core_module(sd) bind(C, name="set_settings_in_core_module_")
            import :: c_intptr_t
            integer(c_intptr_t), intent(in) :: sd
        end subroutine set_settings_in_core_module

        !> Still C++ (eigmode/find_eigmodes.cpp/calc_eigmode.cpp), zersol-
        !> based complex-root search - takes the new opaque core_data
        !> handle in place of the oracle's `core_data *cd`.
        function loop_over_frequences(ind, m, n, cd) bind(C, name="loop_over_frequences") result(stat)
            import :: c_int, c_intptr_t
            integer(c_int), value :: ind, m, n
            integer(c_intptr_t), value :: cd
            integer(c_int) :: stat
        end function loop_over_frequences

        function find_det_zeros(ind, m, n, cd) bind(C, name="find_det_zeros") result(stat)
            import :: c_int, c_intptr_t
            integer(c_int), value :: ind, m, n
            integer(c_intptr_t), value :: cd
            integer(c_int) :: stat
        end function find_det_zeros

        function find_eigmodes(ind, m, n, cd) bind(C, name="find_eigmodes") result(stat)
            import :: c_int, c_intptr_t
            integer(c_int), value :: ind, m, n
            integer(c_intptr_t), value :: cd
            integer(c_int) :: stat
        end function find_eigmodes
    end interface

contains

    function core_data_create_(path) result(handle) bind(C, name="core_data_create_")
        character(kind=c_char), intent(in) :: path(*)
        integer(c_intptr_t) :: handle

        type(core_data_t), pointer :: cd
        integer(c_int) :: pp_size

        !> Mirrors core_data::core_data's sizeof(uintptr_t) != sizeof(pp)
        !> sanity check: pp (constants_m_64bit.f90) and handle are both
        !> meant to be 8-byte (64-bit) quantities on this platform.
        call get_pointer_precision(pp_size)
        if (pp_size /= storage_size(handle)/8) then
            write (*, '(a,i0,a,i0)') &
                'warning: core_data: pointer-precision mismatch: pp=', pp_size, &
                ' sizeof(handle)=', storage_size(handle)/8
            stop 1
        end if

        allocate (cd)
        cd%path2project = c_string_to_fortran(path)

        handle = transfer(c_loc(cd), handle)
    end function core_data_create_

    subroutine core_data_destroy_(handle) bind(C, name="core_data_destroy_")
        integer(c_intptr_t), value :: handle
        type(core_data_t), pointer :: cd
        integer :: ind

        if (handle == 0_c_intptr_t) return
        call handle_to_core_data(handle, cd)

        ! sd intentionally not destroyed: it is one of the static, reused
        ! settings caches above, matching the oracle's own comment ("Do
        ! NOT delete sd intentionally, as the static object ... will be
        ! reused"). bp is a singleton sentinel, nothing to free.

        if (allocated(cd%mda)) then
            do ind = 1, cd%dim
                if (cd%mda(ind) /= 0_c_intptr_t) call mode_data_destroy_(cd%mda(ind))
            end do
            deallocate (cd%mda)
        end if

        deallocate (cd)
    end subroutine core_data_destroy_

    subroutine core_data_calc_and_set_mode_independent_(handle) &
        bind(C, name="core_data_calc_and_set_mode_independent_")
        integer(c_intptr_t), value :: handle
        type(core_data_t), pointer :: cd
        character(kind=c_char), allocatable, target :: cpath(:)
        character(len=1024) :: sd_path2project

        call handle_to_core_data(handle, cd)

        if (index(cd%path2project, 'vacuum') > 0) then
            if (static_settings_vacuum == 0_c_intptr_t) then
                cpath = to_cstr(cd%path2project)
                static_settings_vacuum = settings_create_(cpath)
                call settings_read_settings_(static_settings_vacuum)
            end if
            cd%sd = static_settings_vacuum
        else if (index(cd%path2project, 'flre') > 0) then
            if (static_settings_flre == 0_c_intptr_t) then
                cpath = to_cstr(cd%path2project)
                static_settings_flre = settings_create_(cpath)
                call settings_read_settings_(static_settings_flre)
            end if
            cd%sd = static_settings_flre
        else
            write (*, '(a,a)') &
                'Error: calc_and_set_mode_independent_core_data: unknown project type in path: ', &
                trim(cd%path2project)
            stop 1
        end if

        call set_settings_in_core_module(cd%sd)

        call settings_get_path2project_(cd%sd, sd_path2project)
        cpath = to_cstr(sd_path2project)
        cd%bp = background_create(c_loc(cpath))

        if (get_background_calc_back() > 0) then
            call background_set_profiles_from_files()
        else if (get_background_calc_back() < 0) then
            call background_set_profiles_from_interface()
        else
            write (*, '(a)') &
                'warning: calc_and_set_mode_independent_core_data: unknown flag in background.in!'
            stop 1
        end if
    end subroutine core_data_calc_and_set_mode_independent_

    subroutine core_data_calc_and_set_mode_dependent_antenna_(handle) &
        bind(C, name="core_data_calc_and_set_mode_dependent_antenna_")
        integer(c_intptr_t), value :: handle
        type(core_data_t), pointer :: cd
        real(dp) :: flab_re, flab_im
        complex(dp) :: olab
        character(len=1024) :: sd_path2project
        integer :: ind, m, n

        call handle_to_core_data(handle, cd)

        cd%dim = get_antenna_dma()
        if (allocated(cd%mda)) deallocate (cd%mda)
        allocate (cd%mda(cd%dim))

        call get_antenna_flab(flab_re, flab_im)
        olab = (2.0_dp*pi)*cmplx(flab_re, flab_im, dp)

        do ind = 1, cd%dim
            call get_antenna_mode(int(ind - 1, c_int), m, n)

            call settings_get_path2project_(cd%sd, sd_path2project)
            cd%mda(ind) = mode_data_create_(m, n, real(olab, dp), aimag(olab), cd%sd, cd%bp, &
                trim(sd_path2project)//c_null_char)

            call mode_data_calc_all_mode_data_(cd%mda(ind), 0_c_int)

            call mode_data_destroy_(cd%mda(ind))
            cd%mda(ind) = 0

            call clear_all_data_in_mode_data_module()
        end do
    end subroutine core_data_calc_and_set_mode_dependent_antenna_

    subroutine core_data_calc_and_set_mode_dependent_eigmode_(handle) &
        bind(C, name="core_data_calc_and_set_mode_dependent_eigmode_")
        integer(c_intptr_t), value :: handle
        type(core_data_t), pointer :: cd
        integer :: ind, m, n, stat_unused

        call handle_to_core_data(handle, cd)

        cd%dim = get_antenna_dma()
        if (allocated(cd%mda)) deallocate (cd%mda)
        allocate (cd%mda(cd%dim))
        cd%mda = 0

        do ind = 1, cd%dim
            call get_antenna_mode(int(ind - 1, c_int), m, n)

            select case (get_eigmode_search_flag())
            case (1)
                stat_unused = loop_over_frequences(int(ind - 1, c_int), m, n, handle)
            case (0)
                stat_unused = find_det_zeros(int(ind - 1, c_int), m, n, handle)
            case (-1)
                stat_unused = find_eigmodes(int(ind - 1, c_int), m, n, handle)
            case default
                write (*, '(a,i0,a)') 'Error: unknown search_flag in eigmode options file: ', &
                    get_eigmode_search_flag(), '.'
                stop 1
            end select
        end do
    end subroutine core_data_calc_and_set_mode_dependent_eigmode_

    subroutine core_data_calc_and_set_mode_dependent_antenna_interface_(handle) &
        bind(C, name="core_data_calc_and_set_mode_dependent_antenna_interface_")
        integer(c_intptr_t), value :: handle
        type(core_data_t), pointer :: cd
        real(dp) :: flab_re, flab_im
        complex(dp) :: olab
        character(len=1024) :: sd_path2project
        integer :: ind, m, n

        call handle_to_core_data(handle, cd)

        cd%dim = get_antenna_dma()
        if (allocated(cd%mda)) deallocate (cd%mda)
        allocate (cd%mda(cd%dim))

        call get_antenna_flab(flab_re, flab_im)
        olab = (2.0_dp*pi)*cmplx(flab_re, flab_im, dp)

        do ind = 1, cd%dim
            call get_antenna_mode(int(ind - 1, c_int), m, n)

            call settings_get_path2project_(cd%sd, sd_path2project)
            cd%mda(ind) = mode_data_create_(m, n, real(olab, dp), aimag(olab), cd%sd, cd%bp, &
                trim(sd_path2project)//c_null_char)

            call mode_data_calc_all_mode_data_(cd%mda(ind), 0_c_int)

            call clear_all_data_in_mode_data_module()
        end do
    end subroutine core_data_calc_and_set_mode_dependent_antenna_interface_

    subroutine core_data_calc_and_set_mode_dependent_antenna_interface_mn_(handle, m, n, flag) &
        bind(C, name="core_data_calc_and_set_mode_dependent_antenna_interface_mn_")
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: m, n, flag
        type(core_data_t), pointer :: cd
        real(dp) :: flab_re, flab_im
        complex(dp) :: olab
        character(len=1024) :: sd_path2project
        integer :: ind

        call handle_to_core_data(handle, cd)

        cd%dim = 1
        if (allocated(cd%mda)) deallocate (cd%mda)
        allocate (cd%mda(cd%dim))

        call get_antenna_flab(flab_re, flab_im)
        olab = (2.0_dp*pi)*cmplx(flab_re, flab_im, dp)

        do ind = 1, cd%dim
            call settings_get_path2project_(cd%sd, sd_path2project)
            cd%mda(ind) = mode_data_create_(int(m), int(n), real(olab, dp), aimag(olab), cd%sd, &
                cd%bp, trim(sd_path2project)//c_null_char)

            call mode_data_calc_all_mode_data_(cd%mda(ind), flag)

            call clear_all_data_in_mode_data_module()
        end do
    end subroutine core_data_calc_and_set_mode_dependent_antenna_interface_mn_

    !> ---- accessors for still-C++ callers (wave_code_interface.cpp,
    !> eigmode/calc_eigmode.cpp, eigmode/find_eigmodes.cpp) ----

    integer(c_int) function core_data_get_dim_(handle) bind(C, name="core_data_get_dim_") result(res)
        integer(c_intptr_t), value :: handle
        type(core_data_t), pointer :: cd
        call handle_to_core_data(handle, cd)
        res = cd%dim
    end function core_data_get_dim_

    function core_data_get_mda_element_(handle, ind) bind(C, name="core_data_get_mda_element_") &
        result(res)
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: ind
        integer(c_intptr_t) :: res
        type(core_data_t), pointer :: cd
        call handle_to_core_data(handle, cd)
        res = cd%mda(int(ind) + 1)
    end function core_data_get_mda_element_

    !> Lets still-C++ callers (eigmode/calc_eigmode.cpp's eval_det/
    !> loop_over_frequences) write a freshly-created mode_data handle into
    !> cd->mda[ind] the same way the oracle did via direct field
    !> assignment.
    subroutine core_data_set_mda_element_(handle, ind, val) &
        bind(C, name="core_data_set_mda_element_")
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: ind
        integer(c_intptr_t), value :: val
        type(core_data_t), pointer :: cd
        call handle_to_core_data(handle, cd)
        cd%mda(int(ind) + 1) = val
    end subroutine core_data_set_mda_element_

    function core_data_get_bp_(handle) bind(C, name="core_data_get_bp_") result(res)
        integer(c_intptr_t), value :: handle
        integer(c_intptr_t) :: res
        type(core_data_t), pointer :: cd
        call handle_to_core_data(handle, cd)
        res = cd%bp
    end function core_data_get_bp_

    function core_data_get_sd_(handle) bind(C, name="core_data_get_sd_") result(res)
        integer(c_intptr_t), value :: handle
        integer(c_intptr_t) :: res
        type(core_data_t), pointer :: cd
        call handle_to_core_data(handle, cd)
        res = cd%sd
    end function core_data_get_sd_

    subroutine core_data_get_path2project_(handle, buf) bind(C, name="core_data_get_path2project_")
        integer(c_intptr_t), value :: handle
        character(kind=c_char), intent(out) :: buf(*)
        type(core_data_t), pointer :: cd
        integer :: i, n
        call handle_to_core_data(handle, cd)
        n = len_trim(cd%path2project)
        do i = 1, n
            buf(i) = cd%path2project(i:i)
        end do
        buf(n + 1) = c_null_char
    end subroutine core_data_get_path2project_

    subroutine handle_to_core_data(handle, cd)
        integer(c_intptr_t), value :: handle
        type(core_data_t), pointer, intent(out) :: cd
        type(c_ptr) :: cp
        cp = transfer(handle, cp)
        call c_f_pointer(cp, cd)
    end subroutine handle_to_core_data

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

    function c_string_to_fortran(cstr) result(fstr)
        character(kind=c_char), intent(in) :: cstr(*)
        character(len=1024) :: fstr
        integer :: i
        fstr = ''
        i = 0
        do
            if (cstr(i + 1) == c_null_char .or. i >= 1024) exit
            fstr(i + 1:i + 1) = cstr(i + 1)
            i = i + 1
        end do
    end function c_string_to_fortran

end module kilca_core_data_m
