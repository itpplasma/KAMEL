!<The module for some of the background data

module background

use constants, only: dp, dpc;

implicit none;

!Machine settings:
real (dp) :: rtor;  !big torus radius (cm) of the machine
real (dp) :: rp;    !plasma radius (cm)
real (dp) :: B0;    !toroidal magnetic field (G) at the center

character (1) :: flag_back; !flag for background

real (dp) :: V_gal_sys; !velocity (cm/c) of a moving frame
real (dp) :: V_scale;   !scale of the Vz velocity profile: Vz = V_scale*Vz - V_gal_sys

real (dp) :: zele;  !collision coefficient for electrons
real (dp) :: zion;  !collision coefficient for ions

integer :: flag_debug;  !flag for debugging mode (additional checks are performed in the

!Particles settings:
real (dp), dimension (0:1) :: mass;
real (dp), dimension (0:1) :: charge;

!Other misc parameters:
real (dp) :: huge_factor;

integer :: calc_back;   !flag selecting the background calculation method
integer :: N;           !spline degree for background profiles
character (len=1024) :: path2profiles; !directory holding the profile .dat files

end module

!------------------------------------------------------------------------------

!> Read background.in into the background module, replacing the former C++
!> back_sett::read_settings. Layout: header, rtor/rp/B0, skip, header,
!> path2profiles, calc_back, flag_back, N, V_gal_sys, V_scale, m_i, zele, zion,
!> skip, header, flag_debug, skip; then mass/charge/huge_factor are derived
!> exactly as the C++ constructor did.
subroutine read_background_settings(path) bind(C, name="read_background_settings_")

use, intrinsic :: iso_c_binding, only: c_char, c_null_char
use, intrinsic :: iso_fortran_env, only: error_unit
use constants, only: mp, me, e
use background

character(kind=c_char), dimension(*), intent(in) :: path

character(len=1024) :: fpath, fname, line, before
real(dp) :: m_i
integer :: i, u, ios

fpath = ''
i = 1
do
    if (path(i) == c_null_char .or. i > 1024) exit
    fpath(i:i) = path(i)
    i = i + 1
end do

fname = trim(fpath)//'/background.in'
open (newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
if (ios /= 0) then
    write (error_unit, '(a,a)') "error: read_background_settings: cannot open ", trim(fname)
    error stop
end if

call skip(u)
call value_before_hash(u, before); read (before, *) rtor
call value_before_hash(u, before); read (before, *) rp
call value_before_hash(u, before); read (before, *) B0
call skip(u)

call skip(u)
call value_before_hash(u, before); path2profiles = trim(adjustl(before))
call value_before_hash(u, before); read (before, *) calc_back
call value_before_hash(u, before); flag_back = trim(adjustl(before))
call value_before_hash(u, before); read (before, *) N
call value_before_hash(u, before); read (before, *) V_gal_sys
call value_before_hash(u, before); read (before, *) V_scale
call value_before_hash(u, before); read (before, *) m_i
call value_before_hash(u, before); read (before, *) zele
call value_before_hash(u, before); read (before, *) zion
call skip(u)

call skip(u)
call value_before_hash(u, before); read (before, *) flag_debug
call skip(u)

close (u)

mass(0) = m_i*mp
mass(1) = me
charge(0) = e
charge(1) = -e
huge_factor = 1.0d20

if (flag_debug > 0) call print_background_module_data()

contains

    subroutine skip(unit)
        integer, intent(in) :: unit
        integer :: jos
        read (unit, '(a)', iostat=jos) line
        if (jos /= 0) then
            write (error_unit, '(a)') "error: read_background_settings: read error"
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
            write (error_unit, '(a)') "error: read_background_settings: read error"
            error stop
        end if
        pos = index(buf, '#')
        if (pos > 0) then
            out = buf(1:pos - 1)
        else
            out = buf
        end if
    end subroutine value_before_hash

end subroutine read_background_settings

!------------------------------------------------------------------------------

real(c_double) function get_background_rtor() bind(C, name="get_background_rtor_")
    use, intrinsic :: iso_c_binding, only: c_double
    use background, only: rtor
    get_background_rtor = rtor
end function

real(c_double) function get_background_rp() bind(C, name="get_background_rp_")
    use, intrinsic :: iso_c_binding, only: c_double
    use background, only: rp
    get_background_rp = rp
end function

real(c_double) function get_background_B0() bind(C, name="get_background_B0_")
    use, intrinsic :: iso_c_binding, only: c_double
    use background, only: B0
    get_background_B0 = B0
end function

real(c_double) function get_background_V_gal_sys() bind(C, name="get_background_V_gal_sys_")
    use, intrinsic :: iso_c_binding, only: c_double
    use background, only: V_gal_sys
    get_background_V_gal_sys = V_gal_sys
end function

real(c_double) function get_background_V_scale() bind(C, name="get_background_V_scale_")
    use, intrinsic :: iso_c_binding, only: c_double
    use background, only: V_scale
    get_background_V_scale = V_scale
end function

real(c_double) function get_background_zele() bind(C, name="get_background_zele_")
    use, intrinsic :: iso_c_binding, only: c_double
    use background, only: zele
    get_background_zele = zele
end function

real(c_double) function get_background_zion() bind(C, name="get_background_zion_")
    use, intrinsic :: iso_c_binding, only: c_double
    use background, only: zion
    get_background_zion = zion
end function

integer(c_int) function get_background_flag_debug() bind(C, name="get_background_flag_debug_")
    use, intrinsic :: iso_c_binding, only: c_int
    use background, only: flag_debug
    get_background_flag_debug = flag_debug
end function

real(c_double) function get_background_huge_factor() bind(C, name="get_background_huge_factor_")
    use, intrinsic :: iso_c_binding, only: c_double
    use background, only: huge_factor
    get_background_huge_factor = huge_factor
end function

integer(c_int) function get_background_calc_back() bind(C, name="get_background_calc_back_")
    use, intrinsic :: iso_c_binding, only: c_int
    use background, only: calc_back
    get_background_calc_back = calc_back
end function

integer(c_int) function get_background_N() bind(C, name="get_background_N_")
    use, intrinsic :: iso_c_binding, only: c_int
    use background, only: N
    get_background_N = N
end function

real(c_double) function get_background_mass(i) bind(C, name="get_background_mass_")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use background, only: mass
    integer(c_int), value :: i
    get_background_mass = mass(i)
end function

real(c_double) function get_background_charge(i) bind(C, name="get_background_charge_")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use background, only: charge
    integer(c_int), value :: i
    get_background_charge = charge(i)
end function

function get_background_flag_back() result(ch) bind(C, name="get_background_flag_back_")
    use, intrinsic :: iso_c_binding, only: c_char
    use background, only: flag_back
    character(kind=c_char) :: ch
    ch = flag_back(1:1)
end function

subroutine get_background_path2profiles(out) bind(C, name="get_background_path2profiles_")
    use, intrinsic :: iso_c_binding, only: c_char, c_null_char
    use background, only: path2profiles
    character(kind=c_char), dimension(*), intent(out) :: out
    integer :: i, n
    n = len_trim(path2profiles)
    do i = 1, n
        out(i) = path2profiles(i:i)
    end do
    out(n + 1) = c_null_char
end subroutine

!--------------------------------------------------------------------

subroutine print_background_module_data ()

use background;

implicit none;

print *;

!Machine vessel settings:
print *, "torus big radius: ", rtor;
print *, "plasma radius: ", rp;
print *, "toroidal magnetic field at the center: ", B0;

!Backround field and plasma settings:
print *, "flag for background: ", flag_back;
print *, "velocity of the moving frame: ", V_gal_sys;
print *, "scale factor for the Vz velocity profile: ", V_scale;
print *, "collision coefficient for electrons: ", zele;
print *, "collision coefficient for ions: ", zion;
print *, "flag for debugging mode: ", flag_debug;

!Particles settings:
print *, "particle masses: ", mass;
print *, "particle charges: ", charge;

print *, "huge_factor: ", huge_factor;

end subroutine

!--------------------------------------------------------------------
