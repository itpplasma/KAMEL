!> Output settings, formerly the C++ output_sett class. The Fortran output_data
!> module parses output.in and serves the flags to the Fortran callers through module
!> procedures. Layout matches output_sett::read_settings exactly: a header line,
!> the four run flags, two skipped lines, flag_debug, two skipped lines, then
!> num_quants(=8) per-quantity flags.
module output_data

use constants, only: dp

integer :: flag_background
integer :: flag_emfield
integer :: flag_additional
integer :: flag_dispersion
integer :: flag_debug
integer :: num_quants
integer, dimension(:), allocatable :: flag_quants

end module

!------------------------------------------------------------------------------

module kilca_output_settings_m
    implicit none
    private

    public :: read_output_settings
    public :: get_output_flag_background
    public :: get_output_flag_emfield
    public :: get_output_flag_additional
    public :: get_output_flag_dispersion
    public :: get_output_num_quants
    public :: get_output_flag_quants

contains

    subroutine read_output_settings(path)

use, intrinsic :: iso_fortran_env, only: error_unit
use output_data

        character(len=*), intent(in) :: path

        character(len=1024) :: fpath, fname, buf, before
        integer :: u, ios, pos, k

        fpath = path

        fname = trim(fpath)//'/output.in'
        open (newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write (error_unit, '(a,a)') "error: read_output_settings: cannot open ", trim(fname)
            error stop
        end if

        call skip(u)
        call get_int(u, flag_background)
        call get_int(u, flag_emfield)
        call get_int(u, flag_additional)
        call get_int(u, flag_dispersion)
        call skip(u)
        call skip(u)
        call get_int(u, flag_debug)
        call skip(u)
        call skip(u)

        num_quants = 8
        if (allocated(flag_quants)) deallocate (flag_quants)
        allocate (flag_quants(num_quants))
        do k = 1, num_quants
            call get_int(u, flag_quants(k))
        end do

        close (u)

    contains

        subroutine skip(unit)
            integer, intent(in) :: unit
            integer :: jos
            read (unit, '(a)', iostat=jos) buf
            if (jos /= 0) then
                write (error_unit, '(a)') "error: read_output_settings: read error"
                error stop
            end if
        end subroutine skip

        subroutine get_int(unit, val)
            integer, intent(in) :: unit
            integer, intent(out) :: val
            integer :: jos
            read (unit, '(a)', iostat=jos) buf
            if (jos /= 0) then
                write (error_unit, '(a)') "error: read_output_settings: read error"
                error stop
            end if
            pos = index(buf, '#')
            if (pos > 0) then
                before = buf(1:pos - 1)
            else
                before = buf
            end if
            read (before, *) val
        end subroutine get_int

    end subroutine read_output_settings

!------------------------------------------------------------------------------

    integer(c_int) function get_output_flag_background()
    use, intrinsic :: iso_c_binding, only: c_int
    use output_data, only: flag_background
        get_output_flag_background = flag_background
    end function

    integer(c_int) function get_output_flag_emfield()
    use, intrinsic :: iso_c_binding, only: c_int
    use output_data, only: flag_emfield
        get_output_flag_emfield = flag_emfield
    end function

    integer(c_int) function get_output_flag_additional()
    use, intrinsic :: iso_c_binding, only: c_int
    use output_data, only: flag_additional
        get_output_flag_additional = flag_additional
    end function

    integer(c_int) function get_output_flag_dispersion()
    use, intrinsic :: iso_c_binding, only: c_int
    use output_data, only: flag_dispersion
        get_output_flag_dispersion = flag_dispersion
    end function

    integer(c_int) function get_output_num_quants()
    use, intrinsic :: iso_c_binding, only: c_int
    use output_data, only: num_quants
        get_output_num_quants = num_quants
    end function

    integer(c_int) function get_output_flag_quants(i)
    use, intrinsic :: iso_c_binding, only: c_int
    use output_data, only: flag_quants
        integer(c_int), value :: i
        get_output_flag_quants = flag_quants(i + 1)
    end function

end module kilca_output_settings_m
