! taken from: https://stackoverflow.com/a/17598866/16527499

! TODO:
subroutine read_namelist(n, m, l) bind(C, name="read_namelist")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none

    real(kind=c_double), intent(inout) :: n
    real(kind=c_double), intent(inout) :: m
    integer(kind=c_int), intent(inout) :: l

    namelist /kilcaConfig/ n, m, l

    open (unit=100, file="kilcaConfig.nml", status="old")
    read (unit=100, nml=kilcaConfig)
    close (unit=100)
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
