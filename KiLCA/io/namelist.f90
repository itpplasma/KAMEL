! taken from: https://stackoverflow.com/a/17598866/16527499

subroutine read_namelist(rtor, rp, B0, path2profiles_C, calc_back, flag_back, N, &
                         V_gal_sys, V_scale, m_i, zele, zion, &
                         flag_debug) bind(C, name="read_namelist")

    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_bool, c_size_t, &
                                           c_null_char
    implicit none

    ! background
    ! Machine settings
    real(kind=c_double), intent(out) :: rtor, rp, B0

    ! Background settings
    character(kind=c_char), intent(out) :: path2profiles_C(*)
    integer(kind=c_int), intent(out) :: calc_back
    character(kind=c_char), intent(out) :: flag_back
    integer(kind=c_int), intent(out) :: N
    real(kind=c_double), intent(out) :: V_gal_sys, V_scale
    real(kind=c_double), intent(out) :: m_i, zele, zion

    ! Debug settings
    logical(kind=c_bool), intent(out) :: flag_debug

    ! internal variables
    character(len=256) :: path2profiles
    integer :: i, strlen

    namelist /background/ rtor, rp, B0, path2profiles, calc_back, flag_back, &
                          N, V_gal_sys, V_scale, m_i, zele, zion, flag_debug

    open (unit=100, file="kilca_config.nml", status="old")
    read (unit=100, nml=background)
    close (unit=100)

    ! Copy the string buffer to the output variable
    strlen = min(len_trim(path2profiles), 255)
    do i = 1, strlen
      path2profiles_C(i) = transfer(path2profiles(i:i), c_null_char)
    end do
    path2profiles_C(strlen + 1) = c_null_char
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
