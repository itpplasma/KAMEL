! taken from: https://stackoverflow.com/a/17598866/16527499

subroutine read_namelist( &
    rtor, rp, B0, path2profiles_C, calc_back, flag_back, N, V_gal_sys, V_scale, m_i, &
    zele, zion, & ! background
    fname_C, search_flag, rdim, rfmin, rfmax, idim, ifmin, ifmax, stop_flag, eps_res, &
    eps_abs, eps_rel, delta, test_roots, Nguess, kmin, kmax, fstart_C, & ! eigmode
    ra, wa, I0, flab, dma, modes_C, flag_eigmode, & ! antenna
    flag_debug) bind(C, name="read_namelist")

    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_null_char, &
                                           c_bool, c_double_complex
    implicit none

    ! --- background ---
    ! Machine settings
    real(kind=c_double), intent(out) :: rtor, rp, B0

    ! Background settings
    character(kind=c_char), intent(out) :: path2profiles_C(*)
    integer(kind=c_int), intent(out) :: calc_back
    character(kind=c_char), intent(out) :: flag_back
    integer(kind=c_int), intent(out) :: N
    real(kind=c_double), intent(out) :: V_gal_sys, V_scale
    real(kind=c_double), intent(out) :: m_i, zele, zion

    ! --- eigmode ---
    ! Output
    character(kind=c_char), intent(out) :: fname_C(*)

    ! Search settings
    logical(kind=c_bool), intent(out) :: search_flag

    ! Frequency grid settings
    integer(kind=c_int), intent(out) :: rdim
    real(kind=c_double), intent(out) :: rfmin, rfmax
    integer(kind=c_int), intent(out) :: idim
    real(kind=c_double), intent(out) :: ifmin, ifmax

    ! Stopping criteria
    logical(kind=c_bool), intent(out) :: stop_flag
    real(kind=c_double), intent(out) :: eps_res, eps_abs, eps_rel

    ! Omega derivative settings
    real(kind=c_double), intent(out) :: delta

    ! Test roots
    logical(kind=c_bool), intent(out) :: test_roots

    ! Starting points
    integer(kind=c_int), intent(out) :: Nguess, kmin, kmax
    complex(kind=c_double_complex), intent(out) :: fstart_C(*)

    ! --- antenna ---
    real(kind=c_double), intent(out) :: ra, wa, I0
    complex(kind=c_double_complex), intent(out) :: flab
    integer(kind=c_int), intent(out) :: dma
    integer(kind=c_int), intent(out) :: modes_C(*)
    logical(kind=c_bool), intent(out) :: flag_eigmode

    ! --- debuggroup ---
    logical(kind=c_bool), intent(out) :: flag_debug


    ! --- INTERNAL VARIABLES ---
    logical :: print_content_of_this_file
    integer, parameter :: maxlen = 256 ! for strings
    integer, parameter :: maxsize = 100 ! for arrays
    character(len=maxlen) :: path2profiles
    character(len=maxlen) :: fname
    complex(kind=c_double_complex), dimension(maxsize) :: fstart
    integer(kind=c_int), dimension(maxsize) :: modes
    integer :: unit, i, strlen

    namelist /background/ rtor, rp, B0, path2profiles, calc_back, flag_back, &
                          N, V_gal_sys, V_scale, m_i, zele, zion
    namelist /eigmode/ fname, search_flag, rdim, rfmin, rfmax, idim, ifmin, ifmax, &
                       stop_flag, eps_res, eps_abs, eps_rel, delta, test_roots, &
                       Nguess, kmin, kmax, fstart
    namelist /antenna/ ra, wa, I0, flab, dma, modes, flag_eigmode
    namelist /debuggroup/ flag_debug, print_content_of_this_file

    unit = 10
    open (unit, file="kilca_config.nml", status="old")
    read (unit, nml=background)
    read (unit, nml=eigmode)
    read (unit, nml=antenna)
    read (unit, nml=debuggroup)
    close (unit)

    if (print_content_of_this_file) then
      write(*,*) "Content of the namelist file:"
      write(*,*) "--------------------------------------------------------------------"
      write(*, nml=background)
      write(*, nml=eigmode)
      write(*, nml=antenna)
      write(*, nml=debuggroup)
      write(*,*) "--------------------------------------------------------------------"
    end if

    ! Copy the string buffers to the output variables
    strlen = min(len_trim(path2profiles), maxlen - 1)
    do i = 1, strlen
      path2profiles_C(i) = transfer(path2profiles(i:i), c_null_char)
    end do
    path2profiles_C(strlen + 1) = c_null_char

    strlen = min(len_trim(fname), maxlen - 1)
    do i = 1, strlen
      fname_C(i) = transfer(fname(i:i), c_null_char)
    end do
    fname_C(strlen + 1) = c_null_char

    ! Copy the arrays to the output variables
    do i = 1, maxsize
      fstart_C(i) = fstart(i)
    end do

    do i = 1, maxsize
      modes_C(i) = modes(i)
    end do

end subroutine read_namelist

subroutine read_namelist_unit_test(a, b, c, d_C, &
                                   e_C, f_C) bind(C, name="read_namelist_unit_test")

    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_null_char, &
                                           c_double_complex
    implicit none

    integer(kind=c_int), intent(out) :: a
    real(kind=c_double), intent(out) :: b
    character(kind=c_char), intent(out) :: c
    character(kind=c_char), intent(out) :: d_C(*)
    integer(kind=c_int), intent(out) :: e_C(*)
    complex(kind=c_double_complex), intent(out) :: f_C(*)

    ! internal variables
    character(len=16) :: d
    integer(kind=c_int), dimension(3) :: e
    complex(kind=c_double_complex), dimension(3) :: f
    integer :: unit, i, strlen

    namelist /testnml/ a, b, c, d, e, f

    unit = 100
    open (unit, file="simplified_namelist.nml", status="old")
    read (unit, nml=testnml)
    close (unit)

    ! Copy the string buffer to the output variable
    strlen = min(len_trim(d), 15)
    do i = 1, strlen
      d_C(i) = transfer(d(i:i), c_null_char)
    end do
    d_C(strlen + 1) = c_null_char

    ! Copy the integer array to the output variable
    do i = 1, 3
      e_C(i) = e(i)
    end do

    ! Copy the complex array to the output variable
    do i = 1, 3
      f_C(i) = f(i)
    end do

end subroutine read_namelist_unit_test
