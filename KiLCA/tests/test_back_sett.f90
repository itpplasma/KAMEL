!> Unit test for the Fortran background settings reader (background module).
!>
!> Writes a known background.in matching back_sett::read_settings' layout and
!> checks every getter, including the derived mass/charge/huge_factor values and
!> the single-character flag_back and path2profiles string fields that risked
!> subtle truncation bugs during translation.
program test_back_sett
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_null_char
    implicit none

    interface
        subroutine read_background_settings(path) bind(C, name="read_background_settings_")
            import :: c_char
            character(kind=c_char), dimension(*), intent(in) :: path
        end subroutine
        real(c_double) function get_rtor() bind(C, name="get_background_rtor_")
            import :: c_double
        end function
        real(c_double) function get_rp() bind(C, name="get_background_rp_")
            import :: c_double
        end function
        real(c_double) function get_B0() bind(C, name="get_background_B0_")
            import :: c_double
        end function
        real(c_double) function get_Vgal() bind(C, name="get_background_V_gal_sys_")
            import :: c_double
        end function
        real(c_double) function get_Vscale() bind(C, name="get_background_V_scale_")
            import :: c_double
        end function
        real(c_double) function get_zele() bind(C, name="get_background_zele_")
            import :: c_double
        end function
        real(c_double) function get_zion() bind(C, name="get_background_zion_")
            import :: c_double
        end function
        integer(c_int) function get_flag_debug() bind(C, name="get_background_flag_debug_")
            import :: c_int
        end function
        real(c_double) function get_huge_factor() bind(C, name="get_background_huge_factor_")
            import :: c_double
        end function
        integer(c_int) function get_calc_back() bind(C, name="get_background_calc_back_")
            import :: c_int
        end function
        integer(c_int) function get_N() bind(C, name="get_background_N_")
            import :: c_int
        end function
        real(c_double) function get_mass(i) bind(C, name="get_background_mass_")
            import :: c_int, c_double
            integer(c_int), value :: i
        end function
        real(c_double) function get_charge(i) bind(C, name="get_background_charge_")
            import :: c_int, c_double
            integer(c_int), value :: i
        end function
        function get_flag_back() result(ch) bind(C, name="get_background_flag_back_")
            import :: c_char
            character(kind=c_char) :: ch
        end function
        subroutine get_path2profiles(out) bind(C, name="get_background_path2profiles_")
            import :: c_char
            character(kind=c_char), dimension(*), intent(out) :: out
        end subroutine
    end interface

    real(c_double), parameter :: tol = 1.0d-10
    real(c_double), parameter :: mp = 1.67262158d-24, me = mp/1.8361526675d3, e = 4.8032d-10
    character(kind=c_char), dimension(2) :: cpath
    character(kind=c_char), dimension(1024) :: pbuf
    integer :: u, failures, i
    character(len=1024) :: path2profiles

    failures = 0

    open (newunit=u, file='background.in', status='replace', action='write')
    write (u, '(a)') '#Machine settings:'
    write (u, '(a)') '170.05    #rtor'
    write (u, '(a)') '67.0      #rp'
    write (u, '(a)') '-17563.3704  #B0'
    write (u, '(a)') '#skip'
    write (u, '(a)') '#Background settings:'
    write (u, '(a)') './profiles/   #path'
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

    cpath(1) = '.'
    cpath(2) = c_null_char
    call read_background_settings(cpath)

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

    call get_path2profiles(pbuf)
    path2profiles = ''
    do i = 1, 1024
        if (pbuf(i) == c_null_char) exit
        path2profiles(i:i) = pbuf(i)
    end do
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
