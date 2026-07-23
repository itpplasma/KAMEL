!> Unit test for the Fortran eigmode settings reader (eigmode_sett_data module).
!>
!> Writes a known eigmode.in matching eigmode_sett::read_settings' layout and
!> checks every getter, including the fname string and the fstart complex array.
program test_eigmode_sett
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_null_char
    implicit none

    interface
        subroutine read_eigmode_settings(path) bind(C, name="read_eigmode_settings_")
            import :: c_char
            character(kind=c_char), dimension(*), intent(in) :: path
        end subroutine
        integer(c_int) function get_search_flag() bind(C, name="get_eigmode_search_flag_")
            import :: c_int
        end function
        real(c_double) function get_delta() bind(C, name="get_eigmode_delta_")
            import :: c_double
        end function
        real(c_double) function get_rfmin() bind(C, name="get_eigmode_rfmin_")
            import :: c_double
        end function
        real(c_double) function get_rfmax() bind(C, name="get_eigmode_rfmax_")
            import :: c_double
        end function
        real(c_double) function get_ifmin() bind(C, name="get_eigmode_ifmin_")
            import :: c_double
        end function
        real(c_double) function get_ifmax() bind(C, name="get_eigmode_ifmax_")
            import :: c_double
        end function
        integer(c_int) function get_rdim() bind(C, name="get_eigmode_rdim_")
            import :: c_int
        end function
        integer(c_int) function get_idim() bind(C, name="get_eigmode_idim_")
            import :: c_int
        end function
        integer(c_int) function get_n_zeros() bind(C, name="get_eigmode_n_zeros_")
            import :: c_int
        end function
        integer(c_int) function get_use_winding() bind(C, name="get_eigmode_use_winding_")
            import :: c_int
        end function
        integer(c_int) function get_Nguess() bind(C, name="get_eigmode_Nguess_")
            import :: c_int
        end function
        integer(c_int) function get_kmin() bind(C, name="get_eigmode_kmin_")
            import :: c_int
        end function
        integer(c_int) function get_kmax() bind(C, name="get_eigmode_kmax_")
            import :: c_int
        end function
        integer(c_int) function get_test_roots() bind(C, name="get_eigmode_test_roots_")
            import :: c_int
        end function
        real(c_double) function get_eps_abs() bind(C, name="get_eigmode_eps_abs_")
            import :: c_double
        end function
        real(c_double) function get_eps_rel() bind(C, name="get_eigmode_eps_rel_")
            import :: c_double
        end function
        real(c_double) function get_eps_res() bind(C, name="get_eigmode_eps_res_")
            import :: c_double
        end function
        subroutine get_fname(out) bind(C, name="get_eigmode_fname_")
            import :: c_char
            character(kind=c_char), dimension(*), intent(out) :: out
        end subroutine
        real(c_double) function get_fstart_re(k) bind(C, name="get_eigmode_fstart_re_")
            import :: c_int, c_double
            integer(c_int), value :: k
        end function
        real(c_double) function get_fstart_im(k) bind(C, name="get_eigmode_fstart_im_")
            import :: c_int, c_double
            integer(c_int), value :: k
        end function
    end interface

    real(c_double), parameter :: tol = 1.0d-10
    character(kind=c_char), dimension(2) :: cpath
    character(kind=c_char), dimension(1024) :: fbuf
    character(len=1024) :: fname
    integer :: u, failures, i

    failures = 0

    open (newunit=u, file='eigmode.in', status='replace', action='write')
    write (u, '(a)') '#Output:'
    write (u, '(a)') 'eigmode_search.dat  #fname'
    write (u, '(a)') '#skip'
    write (u, '(a)') '#frequency scan or root search:'
    write (u, '(a)') '-1        #search_flag'
    write (u, '(a)') '#skip'
    write (u, '(a)') '#frequency grid:'
    write (u, '(a)') '10        #rdim'
    write (u, '(a)') '1.0e3     #rfmin'
    write (u, '(a)') '2.0e3     #rfmax'
    write (u, '(a)') '5         #idim'
    write (u, '(a)') '-1.0e2    #ifmin'
    write (u, '(a)') '1.0e2     #ifmax'
    write (u, '(a)') '#skip'
    write (u, '(a)') '#Stopping criteria:'
    write (u, '(a)') '1         #stop_flag'
    write (u, '(a)') '1.0e-8    #eps_res'
    write (u, '(a)') '1.0e-9    #eps_abs'
    write (u, '(a)') '1.0e-10   #eps_rel'
    write (u, '(a)') '#skip'
    write (u, '(a)') '#For derivative:'
    write (u, '(a)') '1.0e-3    #delta'
    write (u, '(a)') '#skip'
    write (u, '(a)') '#For testing:'
    write (u, '(a)') '1         #test_roots'
    write (u, '(a)') '0         #flag_debug'
    write (u, '(a)') '#skip'
    write (u, '(a)') '#ZerSol parameters:'
    write (u, '(a)') '5         #n_zeros'
    write (u, '(a)') '1         #use_winding'
    write (u, '(a)') '#skip'
    write (u, '(a)') '#Starting points:'
    write (u, '(a)') '3         #Nguess'
    write (u, '(a)') '0         #kmin'
    write (u, '(a)') '2         #kmax'
    write (u, '(a)') '#skip'
    write (u, '(a)') '#fstart array:'
    write (u, '(a)') '(1.0e3, 0.5e2)'
    write (u, '(a)') '(1.2e3, -0.3e2)'
    write (u, '(a)') '(1.5e3, 0.0)'
    close (u)

    cpath(1) = '.'
    cpath(2) = c_null_char
    call read_eigmode_settings(cpath)

    call check_i("search_flag", get_search_flag(), -1)
    call check_d("delta", get_delta(), 1.0d-3)
    call check_d("rfmin", get_rfmin(), 1.0d3)
    call check_d("rfmax", get_rfmax(), 2.0d3)
    call check_d("ifmin", get_ifmin(), -1.0d2)
    call check_d("ifmax", get_ifmax(), 1.0d2)
    call check_i("rdim", get_rdim(), 10)
    call check_i("idim", get_idim(), 5)
    call check_i("n_zeros", get_n_zeros(), 5)
    call check_i("use_winding", get_use_winding(), 1)
    call check_i("Nguess", get_Nguess(), 3)
    call check_i("kmin", get_kmin(), 0)
    call check_i("kmax", get_kmax(), 2)
    call check_i("test_roots", get_test_roots(), 1)
    call check_d("eps_abs", get_eps_abs(), 1.0d-9)
    call check_d("eps_rel", get_eps_rel(), 1.0d-10)
    call check_d("eps_res", get_eps_res(), 1.0d-8)

    call get_fname(fbuf)
    fname = ''
    do i = 1, 1024
        if (fbuf(i) == c_null_char) exit
        fname(i:i) = fbuf(i)
    end do
    if (trim(fname) /= 'eigmode_search.dat') then
        write (*, '(a,a)') "FAIL fname: ", trim(fname)
        failures = failures + 1
    end if

    call check_d("fstart_re0", get_fstart_re(0_c_int), 1.0d3)
    call check_d("fstart_im0", get_fstart_im(0_c_int), 0.5d2)
    call check_d("fstart_re2", get_fstart_re(2_c_int), 1.5d3)
    call check_d("fstart_im2", get_fstart_im(2_c_int), 0.0d0)

    if (failures == 0) then
        write (*, '(a)') "PASS: eigmode settings reader matches expected values"
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

end program test_eigmode_sett
