!> Coordinates the four domain settings readers through native module interfaces.
!> Each handle owns its project path; reading populates shared settings state in
!> antenna_data, background, output_data, and eigmode_sett_data. The core caches
!> separate handles for the vacuum and FLRE project paths.
module kilca_settings_m
    use, intrinsic :: iso_c_binding, only: c_intptr_t, c_ptr, c_loc, c_f_pointer
    use kilca_antenna_settings_m, only: read_antenna_settings
    use kilca_background_settings_m, only: read_background_settings
    use kilca_output_settings_m, only: read_output_settings
    use kilca_eigmode_settings_m, only: read_eigmode_settings
    implicit none
    private

    public :: settings_create_, settings_destroy_, settings_read_settings_
    public :: settings_get_path2project_

    type :: settings_t
        character(len=1024) :: path2project = ''
    end type settings_t

contains

    function settings_create_(path) result(handle)
        character(len=*), intent(in) :: path
        integer(c_intptr_t) :: handle
        type(settings_t), pointer :: sd

        allocate (sd)
        sd%path2project = path
        handle = transfer(c_loc(sd), handle)
    end function settings_create_

    subroutine settings_destroy_(handle)
        integer(c_intptr_t), value :: handle
        type(settings_t), pointer :: sd

        if (handle == 0_c_intptr_t) return
        call handle_to_settings(handle, sd)
        deallocate (sd)
    end subroutine settings_destroy_

    !> Read antenna, background, output, and eigenmode settings in that order.
    subroutine settings_read_settings_(handle)
        integer(c_intptr_t), value :: handle
        type(settings_t), pointer :: sd

        call handle_to_settings(handle, sd)

        write (*, '(a,a)') '>> KiLCA: Reading settings from ', trim(sd%path2project)

        call read_antenna_settings(trim(sd%path2project))
        call read_background_settings(trim(sd%path2project))
        call read_output_settings(trim(sd%path2project))
        call read_eigmode_settings(trim(sd%path2project))

        write (*, '(a)') '>> KiLCA: Settings read successfully.'
    end subroutine settings_read_settings_

    subroutine settings_get_path2project_(handle, buf)
        integer(c_intptr_t), value :: handle
        character(len=*), intent(out) :: buf
        type(settings_t), pointer :: sd

        call handle_to_settings(handle, sd)
        buf = sd%path2project
    end subroutine settings_get_path2project_

    subroutine handle_to_settings(handle, sd)
        integer(c_intptr_t), value :: handle
        type(settings_t), pointer, intent(out) :: sd
        type(c_ptr) :: cp
        cp = transfer(handle, cp)
        call c_f_pointer(cp, sd)
    end subroutine handle_to_settings

end module kilca_settings_m
