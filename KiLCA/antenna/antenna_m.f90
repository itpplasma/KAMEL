module antenna_data

!Contains antenna data and current density coefficients Jth(m,n), Jz(m,n) for the current mode

use constants, only: dp, dpc

real (dp) :: ra;       !small radius (cm) of antenna location
real (dp) :: wa;       !current density layer width
real (dp) :: I0;       !current in antenna coils (statamp)
complex (dpc) :: flab; !frequency (Hz) in the laboratory frame
integer :: dma;        !dimension of modes array
integer, dimension(:), allocatable :: modes; !array of modes (m,n)
integer :: flag_debug; !debug flag
integer :: flag_eigmode; !flag for eigmode search

complex(dpc), dimension(2) :: Ja_cyl

end module

!------------------------------------------------------------------------------

!> Read antenna.in and modes.in into the antenna_data module, replacing the
!> former C++ antenna::read_settings. Each scalar line is "value #comment"; the
!> value before '#' is read list-directed, so the complex flab parses from its
!> (re, im) form natively. modes.in holds dma lines "(m,n)".
subroutine read_antenna_settings(path) bind(C, name="read_antenna_settings_")

use, intrinsic :: iso_c_binding, only: c_char, c_null_char
use, intrinsic :: iso_fortran_env, only: error_unit
use antenna_data

character(kind=c_char), dimension(*), intent(in) :: path

character(len=1024) :: fpath, fname, line, before
integer :: i, u, ios, hp, k

fpath = ''
i = 1
do
    if (path(i) == c_null_char .or. i > 1024) exit
    fpath(i:i) = path(i)
    i = i + 1
end do

fname = trim(fpath)//'/antenna.in'
open (newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
if (ios /= 0) then
    write (error_unit, '(a,a)') "error: read_antenna_settings: cannot open ", trim(fname)
    error stop
end if

read (u, '(a)', iostat=ios) line              ! header line
call value_before_hash(u, before); read (before, *) ra
call value_before_hash(u, before); read (before, *) wa
call value_before_hash(u, before); read (before, *) I0
call value_before_hash(u, before); read (before, *) flab
call value_before_hash(u, before); read (before, *) dma
call value_before_hash(u, before); read (before, *) flag_debug
call value_before_hash(u, before); read (before, *) flag_eigmode
close (u)

fname = trim(fpath)//'/modes.in'
open (newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
if (ios /= 0) then
    write (error_unit, '(a,a)') "error: read_antenna_settings: cannot open ", trim(fname)
    error stop
end if

if (allocated(modes)) deallocate (modes)
allocate (modes(2*dma))
do k = 0, dma - 1
    read (u, '(a)', iostat=ios) line
    if (ios /= 0) then
        write (error_unit, '(a)') "error: read_antenna_settings: modes.in read error"
        error stop
    end if
    do i = 1, len(line)
        if (line(i:i) == '(' .or. line(i:i) == ')' .or. line(i:i) == ',') line(i:i) = ' '
    end do
    read (line, *) modes(2*k + 1), modes(2*k + 2)
end do
close (u)

if (flag_debug > 0) call print_antenna_module_data()

contains

    subroutine value_before_hash(unit, out)
        integer, intent(in) :: unit
        character(len=*), intent(out) :: out
        character(len=1024) :: buf
        integer :: pos, jos
        read (unit, '(a)', iostat=jos) buf
        if (jos /= 0) then
            write (error_unit, '(a)') "error: read_antenna_settings: antenna.in read error"
            error stop
        end if
        pos = index(buf, '#')
        if (pos > 0) then
            out = buf(1:pos - 1)
        else
            out = buf
        end if
    end subroutine value_before_hash

end subroutine read_antenna_settings

!------------------------------------------------------------------------------

integer(c_int) function get_antenna_dma() bind(C, name="get_antenna_dma_")
    use, intrinsic :: iso_c_binding, only: c_int
    use antenna_data, only: dma
    get_antenna_dma = dma
end function

subroutine get_antenna_flab(re, im) bind(C, name="get_antenna_flab_")
    use, intrinsic :: iso_c_binding, only: c_double
    use antenna_data, only: flab
    real(c_double), intent(out) :: re, im
    re = real(flab, c_double)
    im = aimag(flab)
end subroutine

subroutine get_antenna_mode(ind, m, n) bind(C, name="get_antenna_mode_")
    use, intrinsic :: iso_c_binding, only: c_int
    use antenna_data, only: modes
    integer(c_int), value :: ind
    integer(c_int), intent(out) :: m, n
    m = modes(2*ind + 1)
    n = modes(2*ind + 2)
end subroutine

real(c_double) function get_antenna_ra() bind(C, name="get_antenna_ra_")
    use, intrinsic :: iso_c_binding, only: c_double
    use antenna_data, only: ra
    get_antenna_ra = ra
end function

real(c_double) function get_antenna_wa() bind(C, name="get_antenna_wa_")
    use, intrinsic :: iso_c_binding, only: c_double
    use antenna_data, only: wa
    get_antenna_wa = wa
end function

integer(c_int) function get_antenna_flag_eigmode() bind(C, name="get_antenna_flag_eigmode_")
    use, intrinsic :: iso_c_binding, only: c_int
    use antenna_data, only: flag_eigmode
    get_antenna_flag_eigmode = flag_eigmode
end function

!------------------------------------------------------------------------------

subroutine set_current_density_in_antenna_module ()

use antenna_data, only: Ja_cyl;

implicit none;

call current_density (Ja_cyl); !Fourier amplutudes of the (Jth, Jz)

end subroutine

!------------------------------------------------------------------------------

subroutine calc_current_density_r_s_p (r, Ja_rsp)

use constants, only: dp, dpc, isqrt2pi
use antenna_data, only: ra, wa, Ja_cyl;

real(dp), intent(in) :: r
complex(dpc), dimension(2), intent(out) :: Ja_rsp

real(dp) :: delta

!(r,s,p) antenna current density Fourier amplutudes: optimize!
call cyl2rsp (r, Ja_cyl(1), Ja_cyl(2), Ja_rsp(1), Ja_rsp(2));

!delta-function approximation:
delta = (isqrt2pi/wa)*exp(-((r-ra)**2)/(2.0*wa**2));

Ja_rsp(1) = Ja_rsp(1)*delta;
Ja_rsp(2) = Ja_rsp(2)*delta;

end subroutine

!------------------------------------------------------------------------------

subroutine print_antenna_module_data ()

use antenna_data;

implicit none;

print *

!Antenna settings:
print *, "antenna radius: ", ra;
print *, "antenna layer width: ", wa;
print *, "antenna coils current: ", I0;
print *, "antenna lab frequency: ", flab;
print *, "dimension of modes array: ", dma;
print *, "flag for debugging mode: ", flag_debug;

end subroutine

!--------------------------------------------------------------------
