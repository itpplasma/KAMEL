program test_integration
    use iso_fortran_env, only: dp => real64
    use integration, only: simpson_nonequi
    implicit none

    integer :: number_of_fails
    real(dp), parameter :: tol_exact = 1.0e-12_dp
    real(dp), parameter :: tol_gen = 5.0e-7_dp
    real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp

    number_of_fails = 0

    call test_quadratic(number_of_fails, tol_exact)
    call test_cubic_convergence(number_of_fails)
    call test_sin_clustered(number_of_fails, 3.0e-6_dp)
    call test_even_points_leftover(number_of_fails, 4.0e-6_dp)
    call test_exp(number_of_fails, 2.0e-7_dp)
    call test_inv_quadratic(number_of_fails, 2.0e-7_dp)

    if (number_of_fails == 0) then
        print *, "All advanced integration tests PASSED."
    else
        print *, "FAILED tests = ", number_of_fails
        error stop "Advanced integration tests failed."
    end if

contains

    subroutine test_quadratic(fails, tol)
        integer, intent(inout) :: fails
        real(dp), intent(in) :: tol
        real(dp), allocatable :: x(:), f(:)
        integer :: n, i
        real(dp) :: approx, exact

        ! Integrate p(x)=2 - 3x + 5x^2 over [0,1]; exact = ∫ (2 -3x +5x^2) dx = 2 - 3/2 + 5/3
        exact = 2.0_dp - 1.5_dp + 5.0_dp / 3.0_dp

        n = 9 ! odd number of points (8 intervals)
        allocate (x(n), f(n))
        do i = 1, n
            x(i) = (real(i - 1, dp) / (n - 1))**1.3_dp ! non-uniform monotone
            f(i) = 2.0_dp - 3.0_dp * x(i) + 5.0_dp * x(i)**2
        end do
        call simpson_nonequi(approx, x, f)
        if (abs(approx - exact) > tol) then
            print *, "test_quadratic failed: err=", abs(approx - exact)
            fails = fails + 1
        end if
        deallocate (x, f)
    end subroutine

    subroutine test_cubic_convergence(fails)
        integer, intent(inout) :: fails
        integer :: n1, n2, i
        real(dp), allocatable :: x1(:), f1(:), x2(:), f2(:)
        real(dp) :: a, b, approx1, approx2, exact, err1, err2, ratio

        a = 0.0_dp; b = 1.0_dp
        exact = 0.25_dp ! ∫_0^1 x^3 dx

        n1 = 15 ! odd
        n2 = 29 ! finer (step ~ half)
        allocate (x1(n1), f1(n1), x2(n2), f2(n2))

        do i = 1, n1
            x1(i) = a + (b - a) * ((real(i - 1, dp) / (n1 - 1))**1.4_dp)
            f1(i) = x1(i)**3
        end do
        do i = 1, n2
            x2(i) = a + (b - a) * ((real(i - 1, dp) / (n2 - 1))**1.4_dp)
            f2(i) = x2(i)**3
        end do

        call simpson_nonequi(approx1, x1, f1)
        call simpson_nonequi(approx2, x2, f2)
        err1 = abs(approx1 - exact)
        err2 = abs(approx2 - exact)
        ! Expected ratio ~ (h1/h2)^4 ~ ( (n2-1)/(n1-1) )^4 ≈ (28/14)^4 = 16
        if (err2 == 0.0_dp) then
            ratio = 1.0e6_dp
        else
            ratio = err1 / err2
        end if
        if (ratio < 8.0_dp) then
            print *, "test_cubic_convergence failed: ratio = ", ratio, " err1=", err1, &
                " err2=", err2
            fails = fails + 1
        end if

        deallocate (x1, f1, x2, f2)
    end subroutine

    subroutine test_sin_clustered(fails, tol)
        integer, intent(inout) :: fails
        real(dp), intent(in) :: tol
        integer :: n, i
        real(dp), allocatable :: x(:), f(:)
        real(dp) :: approx, exact
        exact = 2.0_dp ! ∫_0^π sin(x) dx

        n = 41
        allocate (x(n), f(n))
        do i = 1, n
            ! Chebyshev-like clustering near endpoints
            x(i) = 0.5_dp * pi * (1.0_dp - cos((real(i - 1, dp) / (n - 1)) * pi))
            f(i) = sin(x(i))
        end do
        call simpson_nonequi(approx, x, f)
        if (abs(approx - exact) > tol) then
            print *, "test_sin_clustered failed: err=", abs(approx - exact)
            fails = fails + 1
        end if
        deallocate (x, f)
    end subroutine

    subroutine test_even_points_leftover(fails, tol)
        integer, intent(inout) :: fails
        real(dp), intent(in) :: tol
        integer :: n, i
        real(dp), allocatable :: x(:), f(:)
        real(dp) :: approx, exact
        exact = 2.0_dp ! same sin integral

        n = 40 ! EVEN number of points -> leftover last interval triggers trapezoid
        allocate (x(n), f(n))
        do i = 1, n
            x(i) = pi * (real(i - 1, dp) / (n - 1))**1.2_dp
            f(i) = sin(x(i))
        end do
        call simpson_nonequi(approx, x, f)
        if (abs(approx - exact) > tol) then
            print *, "test_even_points_leftover failed: err=", abs(approx - exact)
            fails = fails + 1
        end if
        deallocate (x, f)
    end subroutine

    subroutine test_exp(fails, tol)
        integer, intent(inout) :: fails
        real(dp), intent(in) :: tol
        integer :: n, i
        real(dp), allocatable :: x(:), f(:)
        real(dp) :: approx, exact
        exact = exp(1.0_dp) - 1.0_dp ! ∫_0^1 e^x dx
        n = 37
        allocate (x(n), f(n))
        do i = 1, n
            x(i) = (real(i - 1, dp) / (n - 1))**1.5_dp
            f(i) = exp(x(i))
        end do
        call simpson_nonequi(approx, x, f)
        if (abs(approx - exact) > tol) then
            print *, "test_exp failed: err=", abs(approx - exact)
            fails = fails + 1
        end if
        deallocate (x, f)
    end subroutine

    subroutine test_inv_quadratic(fails, tol)
        integer, intent(inout) :: fails
        real(dp), intent(in) :: tol
        integer :: n, i
        real(dp), allocatable :: x(:), f(:)
        real(dp) :: approx, exact
        exact = atan(1.0_dp) ! π/4, ∫_0^1 1/(1+x^2) dx
        n = 45
        allocate (x(n), f(n))
        do i = 1, n
            x(i) = (real(i - 1, dp) / (n - 1))**1.1_dp
            f(i) = 1.0_dp / (1.0_dp + x(i)**2)
        end do
        call simpson_nonequi(approx, x, f)
        if (abs(approx - exact) > tol) then
            print *, "test_inv_quadratic failed: err=", abs(approx - exact)
            fails = fails + 1
        end if
        deallocate (x, f)
    end subroutine

end program test_integration
