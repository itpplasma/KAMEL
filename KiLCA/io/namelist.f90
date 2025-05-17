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

subroutine read_namelist_unit_test(a, b, c) bind(C, name="read_namelist_unit_test")
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char
    implicit none

    integer(kind=c_int), intent(inout) :: a
    real(kind=c_double), intent(inout) :: b
    character(kind=c_char), intent(inout) :: c

    namelist /testnml/ a, b, c

    open (unit=100, file="test_namelist.nml", status="old")
    read (unit=100, nml=testnml)
    close (unit=100)
end subroutine read_namelist_unit_test
