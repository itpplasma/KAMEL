program test_integration

    use iso_fortran_env, only: dp => real64
    use integration, only: simpson_nonequi

    implicit none

    real(dp) :: approx
    real(dp), dimension(:), allocatable :: x, f
    integer :: n, i
    real(dp), parameter :: pi = 3.14159265358979323846_dp

    ! Example: integrate sin(x) from 0 to π on a non-equidistant grid
    n = 55
    allocate (x(n), f(n))

    do i = 1, n
        x(i) = 0.5_dp * pi * (1.0_dp - cos((i - 1.0_dp) * pi / (n - 1.0_dp)))
        f(i) = sin(x(i))
    end do

    call simpson_nonequi(approx, x, f)

    if (abs(approx - 2.0_dp) > 1.0e-6_dp) then
        print *, "Approximate integral: ", approx
        print *, "Expected integral: 2.0"
        error stop "Integration test failed: Approximate integral does not match expected value."
    end if

end program test_integration
