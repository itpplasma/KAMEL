!> Eigenmode search settings, formerly the C++ eigmode_sett class. The Fortran
!> eigmode_sett_data module parses eigmode.in and serves the fields to the
!> remaining C++ via bind(C) getters. Layout matches eigmode_sett::read_settings
!> exactly.
module eigmode_sett_data

use constants, only: dp, dpc

character(len=1024) :: fname
integer :: search_flag
integer :: rdim
real(dp) :: rfmin, rfmax
integer :: idim
real(dp) :: ifmin, ifmax
integer :: stop_flag
real(dp) :: eps_res, eps_abs, eps_rel
real(dp) :: delta
integer :: test_roots
integer :: flag_debug
integer :: Nguess
integer :: kmin, kmax
complex(dpc), dimension(:), allocatable :: fstart
integer :: n_zeros
integer :: use_winding

end module

!------------------------------------------------------------------------------

subroutine read_eigmode_settings(path) bind(C, name="read_eigmode_settings_")

use, intrinsic :: iso_c_binding, only: c_char, c_null_char
use, intrinsic :: iso_fortran_env, only: error_unit
use eigmode_sett_data

character(kind=c_char), dimension(*), intent(in) :: path

character(len=1024) :: fpath, fullname, line, before
integer :: i, u, ios, k

fpath = ''
i = 1
do
    if (path(i) == c_null_char .or. i > 1024) exit
    fpath(i:i) = path(i)
    i = i + 1
end do

fullname = trim(fpath)//'/eigmode.in'
open (newunit=u, file=trim(fullname), status='old', action='read', iostat=ios)
if (ios /= 0) then
    write (error_unit, '(a,a)') "error: read_eigmode_settings: cannot open ", trim(fullname)
    error stop
end if

call skip(u)
call value_before_hash(u, before); fname = trim(adjustl(before))
call skip(u)

call skip(u)
call value_before_hash(u, before); read (before, *) search_flag
call skip(u)

call skip(u)
call value_before_hash(u, before); read (before, *) rdim
call value_before_hash(u, before); read (before, *) rfmin
call value_before_hash(u, before); read (before, *) rfmax
call value_before_hash(u, before); read (before, *) idim
call value_before_hash(u, before); read (before, *) ifmin
call value_before_hash(u, before); read (before, *) ifmax
call skip(u)

call skip(u)
call value_before_hash(u, before); read (before, *) stop_flag
call value_before_hash(u, before); read (before, *) eps_res
call value_before_hash(u, before); read (before, *) eps_abs
call value_before_hash(u, before); read (before, *) eps_rel
call skip(u)

call skip(u)
call value_before_hash(u, before); read (before, *) delta
call skip(u)

call skip(u)
call value_before_hash(u, before); read (before, *) test_roots
call value_before_hash(u, before); read (before, *) flag_debug
call skip(u)

call skip(u)
call value_before_hash(u, before); read (before, *) n_zeros
call value_before_hash(u, before); read (before, *) use_winding
call skip(u)

call skip(u)
call value_before_hash(u, before); read (before, *) Nguess
call value_before_hash(u, before); read (before, *) kmin
call value_before_hash(u, before); read (before, *) kmax
call skip(u)

if (allocated(fstart)) deallocate (fstart)
allocate (fstart(0:Nguess - 1))

call skip(u)
do k = 0, Nguess - 1
    call value_before_hash(u, before)
    read (before, *) fstart(k)
end do

close (u)

contains

    subroutine skip(unit)
        integer, intent(in) :: unit
        integer :: jos
        read (unit, '(a)', iostat=jos) line
        if (jos /= 0) then
            write (error_unit, '(a)') "error: read_eigmode_settings: read error"
            error stop
        end if
    end subroutine skip

    subroutine value_before_hash(unit, out)
        integer, intent(in) :: unit
        character(len=*), intent(out) :: out
        character(len=1024) :: buf
        integer :: pos, jos
        read (unit, '(a)', iostat=jos) buf
        if (jos /= 0) then
            write (error_unit, '(a)') "error: read_eigmode_settings: read error"
            error stop
        end if
        pos = index(buf, '#')
        if (pos > 0) then
            out = buf(1:pos - 1)
        else
            out = buf
        end if
    end subroutine value_before_hash

end subroutine read_eigmode_settings

!------------------------------------------------------------------------------

integer(c_int) function get_eigmode_search_flag() bind(C, name="get_eigmode_search_flag_")
    use, intrinsic :: iso_c_binding, only: c_int
    use eigmode_sett_data, only: search_flag
    get_eigmode_search_flag = search_flag
end function

real(c_double) function get_eigmode_delta() bind(C, name="get_eigmode_delta_")
    use, intrinsic :: iso_c_binding, only: c_double
    use eigmode_sett_data, only: delta
    get_eigmode_delta = delta
end function

real(c_double) function get_eigmode_rfmin() bind(C, name="get_eigmode_rfmin_")
    use, intrinsic :: iso_c_binding, only: c_double
    use eigmode_sett_data, only: rfmin
    get_eigmode_rfmin = rfmin
end function

real(c_double) function get_eigmode_rfmax() bind(C, name="get_eigmode_rfmax_")
    use, intrinsic :: iso_c_binding, only: c_double
    use eigmode_sett_data, only: rfmax
    get_eigmode_rfmax = rfmax
end function

real(c_double) function get_eigmode_ifmin() bind(C, name="get_eigmode_ifmin_")
    use, intrinsic :: iso_c_binding, only: c_double
    use eigmode_sett_data, only: ifmin
    get_eigmode_ifmin = ifmin
end function

real(c_double) function get_eigmode_ifmax() bind(C, name="get_eigmode_ifmax_")
    use, intrinsic :: iso_c_binding, only: c_double
    use eigmode_sett_data, only: ifmax
    get_eigmode_ifmax = ifmax
end function

integer(c_int) function get_eigmode_rdim() bind(C, name="get_eigmode_rdim_")
    use, intrinsic :: iso_c_binding, only: c_int
    use eigmode_sett_data, only: rdim
    get_eigmode_rdim = rdim
end function

integer(c_int) function get_eigmode_idim() bind(C, name="get_eigmode_idim_")
    use, intrinsic :: iso_c_binding, only: c_int
    use eigmode_sett_data, only: idim
    get_eigmode_idim = idim
end function

integer(c_int) function get_eigmode_n_zeros() bind(C, name="get_eigmode_n_zeros_")
    use, intrinsic :: iso_c_binding, only: c_int
    use eigmode_sett_data, only: n_zeros
    get_eigmode_n_zeros = n_zeros
end function

integer(c_int) function get_eigmode_use_winding() bind(C, name="get_eigmode_use_winding_")
    use, intrinsic :: iso_c_binding, only: c_int
    use eigmode_sett_data, only: use_winding
    get_eigmode_use_winding = use_winding
end function

integer(c_int) function get_eigmode_Nguess() bind(C, name="get_eigmode_Nguess_")
    use, intrinsic :: iso_c_binding, only: c_int
    use eigmode_sett_data, only: Nguess
    get_eigmode_Nguess = Nguess
end function

integer(c_int) function get_eigmode_kmin() bind(C, name="get_eigmode_kmin_")
    use, intrinsic :: iso_c_binding, only: c_int
    use eigmode_sett_data, only: kmin
    get_eigmode_kmin = kmin
end function

integer(c_int) function get_eigmode_kmax() bind(C, name="get_eigmode_kmax_")
    use, intrinsic :: iso_c_binding, only: c_int
    use eigmode_sett_data, only: kmax
    get_eigmode_kmax = kmax
end function

integer(c_int) function get_eigmode_test_roots() bind(C, name="get_eigmode_test_roots_")
    use, intrinsic :: iso_c_binding, only: c_int
    use eigmode_sett_data, only: test_roots
    get_eigmode_test_roots = test_roots
end function

real(c_double) function get_eigmode_eps_abs() bind(C, name="get_eigmode_eps_abs_")
    use, intrinsic :: iso_c_binding, only: c_double
    use eigmode_sett_data, only: eps_abs
    get_eigmode_eps_abs = eps_abs
end function

real(c_double) function get_eigmode_eps_rel() bind(C, name="get_eigmode_eps_rel_")
    use, intrinsic :: iso_c_binding, only: c_double
    use eigmode_sett_data, only: eps_rel
    get_eigmode_eps_rel = eps_rel
end function

real(c_double) function get_eigmode_eps_res() bind(C, name="get_eigmode_eps_res_")
    use, intrinsic :: iso_c_binding, only: c_double
    use eigmode_sett_data, only: eps_res
    get_eigmode_eps_res = eps_res
end function

subroutine get_eigmode_fname(out) bind(C, name="get_eigmode_fname_")
    use, intrinsic :: iso_c_binding, only: c_char, c_null_char
    use eigmode_sett_data, only: fname
    character(kind=c_char), dimension(*), intent(out) :: out
    integer :: i, n
    n = len_trim(fname)
    do i = 1, n
        out(i) = fname(i:i)
    end do
    out(n + 1) = c_null_char
end subroutine

real(c_double) function get_eigmode_fstart_re(k) bind(C, name="get_eigmode_fstart_re_")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use eigmode_sett_data, only: fstart
    integer(c_int), value :: k
    get_eigmode_fstart_re = real(fstart(k))
end function

real(c_double) function get_eigmode_fstart_im(k) bind(C, name="get_eigmode_fstart_im_")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use eigmode_sett_data, only: fstart
    integer(c_int), value :: k
    get_eigmode_fstart_im = aimag(fstart(k))
end function
