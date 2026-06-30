!> Aggregates the four already-Fortran per-domain settings readers, formerly
!> the C++ settings class (settings.{h,cpp}). Per-instance handle (same
!> c_loc/transfer pattern as cond_profiles etc): core_data caches up to two
!> static settings instances (one per project-type path), each remembering
!> its own path2project, even though read_settings' actual side effect
!> (populating antenna_data/background/output_data/eigmode_sett_data module
!> state) is global, exactly as it already was for the C++ oracle - those
!> four reader modules were already global singletons before this port
!> began.
module kilca_settings_m
    use, intrinsic :: iso_c_binding, only: c_intptr_t, c_char, c_ptr, c_loc, &
        c_f_pointer, c_null_char
    implicit none
    private

    public :: settings_create_, settings_destroy_, settings_read_settings_
    public :: settings_get_path2project_

    type :: settings_t
        character(len=1024) :: path2project = ''
    end type settings_t

    interface
        subroutine read_antenna_settings(path) bind(C, name="read_antenna_settings_")
            import :: c_char
            character(kind=c_char), intent(in) :: path(*)
        end subroutine read_antenna_settings

        subroutine read_background_settings(path) bind(C, name="read_background_settings_")
            import :: c_char
            character(kind=c_char), intent(in) :: path(*)
        end subroutine read_background_settings

        subroutine read_output_settings(path) bind(C, name="read_output_settings_")
            import :: c_char
            character(kind=c_char), intent(in) :: path(*)
        end subroutine read_output_settings

        subroutine read_eigmode_settings(path) bind(C, name="read_eigmode_settings_")
            import :: c_char
            character(kind=c_char), intent(in) :: path(*)
        end subroutine read_eigmode_settings
    end interface

contains

    function settings_create_(path) result(handle) bind(C, name="settings_create_")
        character(kind=c_char), intent(in) :: path(*)
        integer(c_intptr_t) :: handle
        type(settings_t), pointer :: sd

        allocate (sd)
        sd%path2project = c_string_to_fortran(path)
        handle = transfer(c_loc(sd), handle)
    end function settings_create_

    subroutine settings_destroy_(handle) bind(C, name="settings_destroy_")
        integer(c_intptr_t), value :: handle
        type(settings_t), pointer :: sd

        if (handle == 0_c_intptr_t) return
        call handle_to_settings(handle, sd)
        deallocate (sd)
    end subroutine settings_destroy_

    !> Mirrors settings::read_settings exactly: a status line, the four
    !> domain readers in the same order, a completion line.
    subroutine settings_read_settings_(handle) bind(C, name="settings_read_settings_")
        integer(c_intptr_t), value :: handle
        type(settings_t), pointer :: sd
        character(kind=c_char), allocatable :: cpath(:)

        call handle_to_settings(handle, sd)

        write (*, '(a,a)') '>> KiLCA: Reading settings from ', trim(sd%path2project)

        cpath = to_cstr(sd%path2project)

        call read_antenna_settings(cpath)
        call read_background_settings(cpath)
        call read_output_settings(cpath)
        call read_eigmode_settings(cpath)

        write (*, '(a)') '>> KiLCA: Settings read successfully.'
    end subroutine settings_read_settings_

    subroutine settings_get_path2project_(handle, buf) bind(C, name="settings_get_path2project_")
        integer(c_intptr_t), value :: handle
        character(kind=c_char), intent(out) :: buf(*)
        type(settings_t), pointer :: sd
        integer :: i, n

        call handle_to_settings(handle, sd)
        n = len_trim(sd%path2project)
        do i = 1, n
            buf(i) = sd%path2project(i:i)
        end do
        buf(n + 1) = c_null_char
    end subroutine settings_get_path2project_

    subroutine handle_to_settings(handle, sd)
        integer(c_intptr_t), value :: handle
        type(settings_t), pointer, intent(out) :: sd
        type(c_ptr) :: cp
        cp = transfer(handle, cp)
        call c_f_pointer(cp, sd)
    end subroutine handle_to_settings

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

end module kilca_settings_m
