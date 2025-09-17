module integration

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

    !> Composite Simpson’s rule on a non-equidistant grid.
    !! This subroutine computes the numerical integral of a function using the
    !! composite Simpson's rule on a non-equidistant grid.
    !!
    !! @param[out] integral Computed approximation to ∫ x(1)->x(n) f(x) dx.
    !! @param[in] x Array of abscissae, must be sorted in ascending order.
    !! @param[in] f Array of function values at the abscissae x.
    !! @pre The number of entries in x and f must be odd (=> even number of intervals).
    !! @note The method uses Lagrange basis polynomials for integration.
    !! @see https://doi.org/10.1007/978-3-319-27265-8, p. 38f.
    subroutine simpson_nonequi(integral, x, f)
        real(dp), intent(out) :: integral
        real(dp), dimension(:), intent(in) :: x, f

        integer :: i, npts, nmax
        real(dp) :: xi, xi1, xi2
        real(dp) :: d0, d1, d2
        real(dp) :: w0, w1, w2
        real(dp) :: hlast

        npts = size(x)
        integral = 0.0_dp

        if (size(f) /= npts) then
            print *, "simpson_nonequi: size(x) /= size(f)"
            error stop
        end if
        if (npts < 3) then
            print *, "simpson_nonequi: need at least 3 points"
            error stop
        end if

        ! check monotonic increasing grid
        do i = 1, npts - 1
            if (x(i + 1) <= x(i)) then
                print *, "simpson_nonequi: x must be strictly increasing"
                error stop
            end if
        end do

        ! number of points to use with composite Simpson (make it odd)
        if (mod(npts - 1, 2) == 0) then
            nmax = npts
        else
            nmax = npts - 1 ! leave last interval for trapezoid below
        end if

        do i = 1, nmax - 2, 2
            xi = x(i)
            xi1 = x(i + 1)
            xi2 = x(i + 2)

            ! Denominators for Lagrange basis polynomials
            d0 = (xi - xi1) * (xi - xi2)
            d1 = (xi1 - xi) * (xi1 - xi2)
            d2 = (xi2 - xi) * (xi2 - xi1)

            ! Compute weights by integrating Lagrange basis on [xi,xi2]
            ! L0(x) = (x-xi1)(x-xi2) / d0
            ! Integral from xi to xi2: expand and integrate term by term
            w0 = ((xi2**3 - xi**3) / 3.0_dp - (xi1 + xi2) * (xi2**2 - xi**2) / 2.0_dp + &
                  xi1 * xi2 * (xi2 - xi)) / d0

            ! L1(x) = (x-xi)(x-xi2) / d1
            w1 = ((xi2**3 - xi**3) / 3.0_dp - (xi + xi2) * (xi2**2 - xi**2) / 2.0_dp + &
                  xi * xi2 * (xi2 - xi)) / d1

            ! L2(x) = (x-xi)(x-xi1) / d2
            w2 = ((xi2**3 - xi**3) / 3.0_dp - (xi + xi1) * (xi2**2 - xi**2) / 2.0_dp + &
                  xi * xi1 * (xi2 - xi)) / d2

            integral = integral + w0 * f(i) + w1 * f(i + 1) + w2 * f(i + 2)
        end do

        ! if there was an extra last interval (npts even), integrate it with trapezoid
        if (nmax /= npts) then
            hlast = x(npts) - x(npts - 1)
            integral = integral + 0.5_dp * (f(npts - 1) + f(npts)) * hlast
        end if

    end subroutine simpson_nonequi

end module integration
