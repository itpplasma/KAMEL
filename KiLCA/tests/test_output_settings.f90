!> Unit test for the Fortran output settings reader (output_data module).
!>
!> Writes a known output.in matching output_sett::read_settings' layout (header,
!> four run flags, two skips, flag_debug, two skips, then 8 flag_quants), runs
!> read_output_settings_, and checks every getter. Locks the skip positions and
!> the flag_quants array indexing.
program test_output_settings
    use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
    implicit none

    interface
        subroutine read_output_settings(path) bind(C, name="read_output_settings_")
            import :: c_char
            character(kind=c_char), dimension(*), intent(in) :: path
        end subroutine
        integer(c_int) function get_bg() bind(C, name="get_output_flag_background_")
            import :: c_int
        end function
        integer(c_int) function get_em() bind(C, name="get_output_flag_emfield_")
            import :: c_int
        end function
        integer(c_int) function get_add() bind(C, name="get_output_flag_additional_")
            import :: c_int
        end function
        integer(c_int) function get_disp() bind(C, name="get_output_flag_dispersion_")
            import :: c_int
        end function
        integer(c_int) function get_nq() bind(C, name="get_output_num_quants_")
            import :: c_int
        end function
        integer(c_int) function get_fq(i) bind(C, name="get_output_flag_quants_")
            import :: c_int
            integer(c_int), value :: i
        end function
    end interface

    character(kind=c_char), dimension(2) :: cpath
    integer :: u, failures, k

    failures = 0

    open (newunit=u, file='output.in', status='replace', action='write')
    write (u, '(a)') '#Run settings:'
    write (u, '(a)') '1   #flag_background'
    write (u, '(a)') '2   #flag_emfield'
    write (u, '(a)') '0   #flag_additional'
    write (u, '(a)') '1   #flag_dispersion'
    write (u, '(a)') '#skip'
    write (u, '(a)') '#Debug:'
    write (u, '(a)') '0   #flag_debug'
    write (u, '(a)') '#skip'
    write (u, '(a)') '#Quantity flags:'
    do k = 0, 7
        write (u, '(i0,a,i0)') k, '   #flag_quants ', k
    end do
    close (u)

    cpath(1) = '.'
    cpath(2) = c_null_char
    call read_output_settings(cpath)

    call check("flag_background", get_bg(), 1)
    call check("flag_emfield", get_em(), 2)
    call check("flag_additional", get_add(), 0)
    call check("flag_dispersion", get_disp(), 1)
    call check("num_quants", get_nq(), 8)
    do k = 0, 7
        call check("flag_quants", get_fq(k), k)
    end do

    if (failures == 0) then
        write (*, '(a)') "PASS: output settings reader matches expected values"
    else
        write (*, '(a,i0)') "FAILED: ", failures
        stop 1
    end if

contains

    subroutine check(label, got, want)
        character(*), intent(in) :: label
        integer, intent(in) :: got, want
        if (got /= want) then
            write (*, '(a,a,2(1x,i0))') "FAIL ", label, got, want
            failures = failures + 1
        end if
    end subroutine

end program test_output_settings
