!> Unit test for the Fortran antenna settings reader (antenna_data module).
!>
!> Writes a known antenna.in / modes.in in the "value #comment" format the C++
!> antenna::read_settings used, runs read_antenna_settings_, and checks every
!> exposed field through the bind(C) getters. This locks the parser behavior the
!> former C++ class provided (doubles, the (re,im) complex, ints, and the
!> "(m,n)" mode list).
program test_antenna_settings
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_null_char
    implicit none

    interface
        subroutine read_antenna_settings(path) bind(C, name="read_antenna_settings_")
            import :: c_char
            character(kind=c_char), dimension(*), intent(in) :: path
        end subroutine
        integer(c_int) function get_antenna_dma() bind(C, name="get_antenna_dma_")
            import :: c_int
        end function
        subroutine get_antenna_flab(re, im) bind(C, name="get_antenna_flab_")
            import :: c_double
            real(c_double), intent(out) :: re, im
        end subroutine
        subroutine get_antenna_mode(ind, m, n) bind(C, name="get_antenna_mode_")
            import :: c_int
            integer(c_int), value :: ind
            integer(c_int), intent(out) :: m, n
        end subroutine
        real(c_double) function get_antenna_ra() bind(C, name="get_antenna_ra_")
            import :: c_double
        end function
        real(c_double) function get_antenna_wa() bind(C, name="get_antenna_wa_")
            import :: c_double
        end function
        integer(c_int) function get_antenna_flag_eigmode() bind(C, name="get_antenna_flag_eigmode_")
            import :: c_int
        end function
    end interface

    real(c_double), parameter :: tol = 1.0d-12
    character(kind=c_char), dimension(2) :: cpath
    integer :: u, failures, m, n
    real(c_double) :: re, im

    failures = 0

    open (newunit=u, file='antenna.in', status='replace', action='write')
    write (u, '(a)') '#Antenna settings:'
    write (u, '(a)') '70.0         #small radius'
    write (u, '(a)') '0.5          #current density layer width'
    write (u, '(a)') '1.0e13       #current in the coils'
    write (u, '(a)') '(1.0e0, 2.0e0)   #complex frequency'
    write (u, '(a)') '3            #number of antenna modes'
    write (u, '(a)') '0            #flag for debugging'
    write (u, '(a)') '1            #flag to solve an eigenmode problem'
    write (u, '(a)') '#footer'
    close (u)

    open (newunit=u, file='modes.in', status='replace', action='write')
    write (u, '(a)') '(3,2)'
    write (u, '(a)') '(4,2)'
    write (u, '(a)') '(5,2)'
    close (u)

    cpath(1) = '.'
    cpath(2) = c_null_char
    call read_antenna_settings(cpath)

    call check_d("ra", get_antenna_ra(), 70.0d0)
    call check_d("wa", get_antenna_wa(), 0.5d0)
    call check_i("dma", get_antenna_dma(), 3)
    call check_i("flag_eigmode", get_antenna_flag_eigmode(), 1)

    call get_antenna_flab(re, im)
    call check_d("flab_re", re, 1.0d0)
    call check_d("flab_im", im, 2.0d0)

    call get_antenna_mode(0_c_int, m, n)
    call check_i("mode0_m", m, 3); call check_i("mode0_n", n, 2)
    call get_antenna_mode(2_c_int, m, n)
    call check_i("mode2_m", m, 5); call check_i("mode2_n", n, 2)

    if (failures == 0) then
        write (*, '(a)') "PASS: antenna settings reader matches expected values"
    else
        write (*, '(a,i0)') "FAILED: ", failures
        stop 1
    end if

contains

    subroutine check_d(label, got, want)
        character(*), intent(in) :: label
        real(c_double), intent(in) :: got, want
        if (abs(got - want) > tol) then
            write (*, '(a,a,2(1x,es23.16))') "FAIL ", label, got, want
            failures = failures + 1
        end if
    end subroutine

    subroutine check_i(label, got, want)
        character(*), intent(in) :: label
        integer, intent(in) :: got, want
        if (got /= want) then
            write (*, '(a,a,2(1x,i0))') "FAIL ", label, got, want
            failures = failures + 1
        end if
    end subroutine

end program test_antenna_settings
