!> Unit test for the Fortran background settings reader (background module).
!>
!> Writes a known background.in matching back_sett::read_settings' layout and
!> checks every getter, including the derived mass/charge/huge_factor values and
!> the single-character flag_back and path2profiles string fields that risked
!> subtle truncation bugs during translation.
program test_back_sett
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_c_binding, only: c_double
    use kilca_background_settings_m, only: &
        read_background_settings, get_rtor => &
        get_background_rtor, get_rp => &
        get_background_rp, get_B0 => &
        get_background_B0, get_Vgal => &
        get_background_V_gal_sys, get_Vscale => &
        get_background_V_scale, get_zele => &
        get_background_zele, get_zion => &
        get_background_zion, get_flag_debug => &
        get_background_flag_debug, get_huge_factor => &
        get_background_huge_factor, get_calc_back => &
        get_background_calc_back, get_N => &
        get_background_N, get_mass => &
        get_background_mass, get_charge => &
        get_background_charge, get_flag_back => &
        get_background_flag_back, get_path2profiles => &
        get_background_path2profiles
    implicit none

    real(c_double), parameter :: tol = 1.0d-10
    real(c_double), parameter :: mp = 1.67262158d-24, me = mp/1.8361526675d3, e = 4.8032d-10
    integer :: u, failures
    character(len=1024) :: path2profiles

    failures = 0

    open (newunit=u, file='background.in', status='replace', action='write')
    write (u, '(a)') '#Machine settings:'
    write (u, '(a)') '170.05    #rtor'
    write (u, '(a)') '67.0      #rp'
    write (u, '(a)') '-17563.3704  #B0'
    write (u, '(a)') '#skip'
    write (u, '(a)') '#Background settings:'
    write (u, '(a)') achar(9)//'./profiles/'//achar(9)//' #path'
    write (u, '(a)') '1         #calc_back'
    write (u, '(a)') 'f         #flag_back'
    write (u, '(a)') '9         #N'
    write (u, '(a)') '1.e9      #V_gal_sys'
    write (u, '(a)') '1.0e0     #V_scale'
    write (u, '(a)') '2.0       #m_i'
    write (u, '(a)') '1.0e-0    #zele'
    write (u, '(a)') '1.0e-0    #zion'
    write (u, '(a)') '#skip'
    write (u, '(a)') '#Checkings setting:'
    write (u, '(a)') '0         #flag_debug'
    write (u, '(a)') '#skip'
    close (u)

    call read_background_settings('.')

    call check_d("rtor", get_rtor(), 170.05d0)
    call check_d("rp", get_rp(), 67.0d0)
    call check_d("B0", get_B0(), -17563.3704d0)
    call check_i("calc_back", get_calc_back(), 1)
    call check_i("N", get_N(), 9)
    call check_d("V_gal_sys", get_Vgal(), 1.0d9)
    call check_d("V_scale", get_Vscale(), 1.0d0)
    call check_d("zele", get_zele(), 1.0d0)
    call check_d("zion", get_zion(), 1.0d0)
    call check_i("flag_debug", get_flag_debug(), 0)

    call check_d("mass0", get_mass(0_c_int), 2.0d0*mp)
    call check_d("mass1", get_mass(1_c_int), me)
    call check_d("charge0", get_charge(0_c_int), e)
    call check_d("charge1", get_charge(1_c_int), -e)
    call check_d("huge_factor", get_huge_factor(), 1.0d20)

    if (get_flag_back() /= 'f') then
        write (*, '(a)') "FAIL flag_back"
        failures = failures + 1
    end if

    call get_path2profiles(path2profiles)
    if (trim(path2profiles) /= './profiles/') then
        write (*, '(a,a)') "FAIL path2profiles: ", trim(path2profiles)
        failures = failures + 1
    end if

    if (failures == 0) then
        write (*, '(a)') "PASS: background settings reader matches expected values"
    else
        write (*, '(a,i0)') "FAILED: ", failures
        stop 1
    end if

contains

    subroutine check_d(label, got, want)
        character(*), intent(in) :: label
        real(c_double), intent(in) :: got, want
        if (abs(got - want) > tol*max(1.0d0, abs(want))) then
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

end program test_back_sett
