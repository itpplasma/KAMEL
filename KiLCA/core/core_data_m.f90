!> The top-level KiLCA orchestrator, formerly the C++ core_data class
!> (core.{h,cpp}). Per-instance handle (same c_loc/transfer pattern as
!> cond_profiles etc - core_data is genuinely instantiated multiple times,
!> once per wave_code_interface_64bit.f90/QL-Balance call).
!>
!> Eigenmode orchestration is implemented in a submodule so the eigenmode
!> solver can import this module's native accessors without a module cycle.
module kilca_core_data_m
    use kilca_legacy_interfaces_m, only: get_pointer_precision
    use kilca_legacy_interfaces_m, only: clear_all_data_in_mode_data_module
    use kilca_legacy_interfaces_m, only: set_settings_in_core_module
    use kilca_antenna_settings_m, only: get_antenna_dma
    use kilca_antenna_settings_m, only: get_antenna_flab
    use kilca_antenna_settings_m, only: get_antenna_mode
    use kilca_background_data_m, only: background_create
    use kilca_background_data_m, only: background_set_profiles_from_files
    use kilca_background_data_m, only: background_set_profiles_from_interface
    use kilca_background_settings_m, only: get_background_calc_back
    use kilca_eigmode_settings_m, only: get_eigmode_search_flag
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_ptr, c_loc, &
                                                                              c_f_pointer
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
        module subroutine core_data_calc_and_set_mode_dependent_eigmode_(handle)
            integer(c_intptr_t), value :: handle
        end subroutine
    end interface

contains

    function core_data_create_(path) result(handle)
        character(len=*), intent(in) :: path
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
        cd%path2project = path

        handle = transfer(c_loc(cd), handle)
    end function core_data_create_

    subroutine core_data_destroy_(handle)
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

    subroutine core_data_calc_and_set_mode_independent_(handle)
        integer(c_intptr_t), value :: handle
        type(core_data_t), pointer :: cd
        character(len=1024) :: sd_path2project

        call handle_to_core_data(handle, cd)

        if (index(cd%path2project, 'vacuum') > 0) then
            if (static_settings_vacuum == 0_c_intptr_t) then

                static_settings_vacuum = settings_create_(cd%path2project)
                call settings_read_settings_(static_settings_vacuum)
            end if
            cd%sd = static_settings_vacuum
        else if (index(cd%path2project, 'flre') > 0) then
            if (static_settings_flre == 0_c_intptr_t) then

                static_settings_flre = settings_create_(cd%path2project)
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

        cd%bp = background_create(sd_path2project)

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

    subroutine core_data_calc_and_set_mode_dependent_antenna_(handle)
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
        olab = (2.0_dp * pi) * cmplx(flab_re, flab_im, dp)

        do ind = 1, cd%dim
            call get_antenna_mode(int(ind - 1, c_int), m, n)

            call settings_get_path2project_(cd%sd, sd_path2project)
            cd%mda(ind) = mode_data_create_(m, n, real(olab, dp), aimag(olab), &
                                            cd%sd, cd%bp, &
                                            trim(sd_path2project))

            call mode_data_calc_all_mode_data_(cd%mda(ind), 0_c_int)

            call mode_data_destroy_(cd%mda(ind))
            cd%mda(ind) = 0

            call clear_all_data_in_mode_data_module()
        end do
    end subroutine core_data_calc_and_set_mode_dependent_antenna_

    subroutine core_data_calc_and_set_mode_dependent_antenna_interface_(handle)
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
        olab = (2.0_dp * pi) * cmplx(flab_re, flab_im, dp)

        do ind = 1, cd%dim
            call get_antenna_mode(int(ind - 1, c_int), m, n)

            call settings_get_path2project_(cd%sd, sd_path2project)
            cd%mda(ind) = mode_data_create_(m, n, real(olab, dp), aimag(olab), &
                                            cd%sd, cd%bp, &
                                            trim(sd_path2project))

            call mode_data_calc_all_mode_data_(cd%mda(ind), 0_c_int)

            call clear_all_data_in_mode_data_module()
        end do
    end subroutine core_data_calc_and_set_mode_dependent_antenna_interface_

    subroutine core_data_calc_and_set_mode_dependent_antenna_interface_mn_(handle, &
                                                                           m, n, flag)
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
        olab = (2.0_dp * pi) * cmplx(flab_re, flab_im, dp)

        do ind = 1, cd%dim
            call settings_get_path2project_(cd%sd, sd_path2project)
            cd%mda(ind) = mode_data_create_(int(m), int(n), real(olab, dp), &
                                            aimag(olab), cd%sd, &
                                            cd%bp, trim(sd_path2project))

            call mode_data_calc_all_mode_data_(cd%mda(ind), flag)

            call clear_all_data_in_mode_data_module()
        end do
    end subroutine core_data_calc_and_set_mode_dependent_antenna_interface_mn_

    !> ---- accessors for the wave-code interface and eigenmode solver ----

    integer(c_int) function core_data_get_dim_(handle) result(res)
        integer(c_intptr_t), value :: handle
        type(core_data_t), pointer :: cd
        call handle_to_core_data(handle, cd)
        res = cd%dim
    end function core_data_get_dim_

    function core_data_get_mda_element_(handle, ind) &
        result(res)
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: ind
        integer(c_intptr_t) :: res
        type(core_data_t), pointer :: cd
        call handle_to_core_data(handle, cd)
        res = cd%mda(int(ind) + 1)
    end function core_data_get_mda_element_

    !> Lets callers (eigmode/eigmode_solve_m.f90's determinant evaluation)
    !> write a freshly-created mode_data handle into cd->mda[ind] the same
    !> way the oracle did via direct field assignment.
    subroutine core_data_set_mda_element_(handle, ind, val)
        integer(c_intptr_t), value :: handle
        integer(c_int), value :: ind
        integer(c_intptr_t), value :: val
        type(core_data_t), pointer :: cd
        call handle_to_core_data(handle, cd)
        cd%mda(int(ind) + 1) = val
    end subroutine core_data_set_mda_element_

    function core_data_get_bp_(handle) result(res)
        integer(c_intptr_t), value :: handle
        integer(c_intptr_t) :: res
        type(core_data_t), pointer :: cd
        call handle_to_core_data(handle, cd)
        res = cd%bp
    end function core_data_get_bp_

    function core_data_get_sd_(handle) result(res)
        integer(c_intptr_t), value :: handle
        integer(c_intptr_t) :: res
        type(core_data_t), pointer :: cd
        call handle_to_core_data(handle, cd)
        res = cd%sd
    end function core_data_get_sd_

    subroutine core_data_get_path2project_(handle, buf)
        integer(c_intptr_t), value :: handle
        character(len=*), intent(out) :: buf
        type(core_data_t), pointer :: cd
        integer :: i, n
        call handle_to_core_data(handle, cd)
        buf = cd%path2project
    end subroutine core_data_get_path2project_

    subroutine handle_to_core_data(handle, cd)
        integer(c_intptr_t), value :: handle
        type(core_data_t), pointer, intent(out) :: cd
        type(c_ptr) :: cp
        cp = transfer(handle, cp)
        call c_f_pointer(cp, cd)
    end subroutine handle_to_core_data

end module kilca_core_data_m
