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

        integer :: i, n
        real(dp) :: xi, xi1, xi2
        real(dp) :: d0, d1, d2
        real(dp) :: w0, w1, w2

        n = size(x)

        if (mod(n - 1, 2) /= 0) then
            write (*, "(a, a, i0, a)") "Warning: The supplied arrays for Simpson integration ", &
                "should have an odd number of elements but have ", n, &
                ". Skipping the last element."
            n = n - 1
        end if

        integral = 0.0d0

        ! Loop over each pair of intervals: nodes (i,i+1,i+2)
        do i = 1, n - 2, 2
            xi = x(i)
            xi1 = x(i + 1)
            xi2 = x(i + 2)

            ! Denominators for Lagrange basis polynomials
            d0 = (xi - xi1) * (xi - xi2)
            d1 = (xi1 - xi) * (xi1 - xi2)
            d2 = (xi2 - xi) * (xi2 - xi1)

            ! Compute weights by integrating each Lagrange basis:
            ! w_k = ∫_{xi}^{xi2} L_k(x) dx, where
            !       L0(x) = (x-xi1)(x-xi2)/d0, L1(x) = (x-xi)(x-xi2)/d1, L2(x) = (x-xi)(x-xi1)/d2
            ! ∫ (x^2 + px + q) dx = x^3/3 + px^2/2 + qx + c
            w0 = (xi2**3 / 3.d0 - (xi1 + xi2) * xi2**2 / 2.d0 + xi1 * xi2 * xi2 &
                  - (xi**3 / 3.d0 - (xi1 + xi2) * xi**2 / 2.d0 + xi1 * xi2 * xi)) / d0

            w1 = (xi2**3 / 3.d0 - (xi + xi2) * xi2**2 / 2.d0 + xi * xi2 * xi2 &
                  - (xi**3 / 3.d0 - (xi + xi2) * xi**2 / 2.d0 + xi * xi2 * xi)) / d1

            w2 = (xi2**3 / 3.d0 - (xi + xi1) * xi2**2 / 2.d0 + xi * xi1 * xi2 &
                  - (xi**3 / 3.d0 - (xi + xi1) * xi**2 / 2.d0 + xi * xi1 * xi)) / d2

            integral = integral + w0 * f(i) + w1 * f(i + 1) + w2 * f(i + 2)
        end do

    end subroutine simpson_nonequi

end module integration
