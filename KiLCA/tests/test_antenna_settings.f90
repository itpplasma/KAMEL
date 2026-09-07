!> Unit test for the Fortran antenna settings reader (antenna_data module).
!>
!> Writes a known antenna.in / modes.in in the "value #comment" format the C++
!> antenna::read_settings used, runs read_antenna_settings, and checks every
!> exposed field through the module getters. This locks the parser behavior the
!> former C++ class provided (doubles, the (re,im) complex, ints, and the
!> "(m,n)" mode list).
program test_antenna_settings
    use, intrinsic :: iso_c_binding, only: c_int
    use, intrinsic :: iso_c_binding, only: c_double, c_intptr_t
    use kilca_antenna_settings_m, only: &
        read_antenna_settings, get_antenna_dma, get_antenna_flab, get_antenna_mode, &
        get_antenna_ra, get_antenna_wa, get_antenna_flag_eigmode
    use kilca_settings_m, only: settings_create_, settings_destroy_, &
                                settings_get_path2project_
    implicit none

    real(c_double), parameter :: tol = 1.0d-12
    integer :: u, failures, m, n
    real(c_double) :: re, im

    failures = 0
    call check_native_settings_path()

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

    call read_antenna_settings('.')

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

    subroutine check_native_settings_path()
        character(len=*), parameter :: project_path = 'project with spaces/case'
        character(len=1024) :: result_path
        character(len=12) :: short_path
        integer(c_intptr_t) :: handle

        handle = settings_create_(project_path)
        call settings_get_path2project_(handle, result_path)
        if (result_path /= project_path) then
            write (*, '(a,a)') 'FAIL native settings path roundtrip: ', &
                trim(result_path)
            failures = failures + 1
        end if
        call settings_get_path2project_(handle, short_path)
        if (short_path /= project_path(:len(short_path))) then
            write (*, '(a,a)') 'FAIL native settings path truncation: ', short_path
            failures = failures + 1
        end if
        call settings_destroy_(handle)
    end subroutine check_native_settings_path

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
