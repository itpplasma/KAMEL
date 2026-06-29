module gauss_quadrature_m
    !> Shared Gauss-Legendre quadrature helpers for the integration paths.
    !> Previously duplicated byte-for-byte: compute_nodes_weights in
    !> integrals_gauss_m / integrals_rkf45_m, and calc_xbar in
    !> integrands_gauss_m / integrands_rkf45_m.

    use KIM_kinds_m, only: dp

    implicit none
    private

    public :: compute_nodes_weights, calc_xbar

contains

    subroutine compute_nodes_weights(n, x, w)

        implicit none

        integer, intent(in) :: n
        real(dp), intent(out) :: x(n), w(n)

        integer :: i, j, m
        real(dp) :: z, z1, p1, p2, p3, pp
        real(dp), parameter :: EPS = 1.0d-14

        m = (n + 1) / 2  ! number of roots to compute (use symmetry)

        do i = 1, m
            ! Initial guess (Chebyshev-Gauss nodes)
            z = cos(3.14159265358979d0 * (i - 0.25d0) / (n + 0.5d0))

            do
                p1 = 1.0d0
                p2 = 0.0d0
                ! Evaluate Legendre polynomial using recurrence
                do j = 1, n
                    p3 = p2
                    p2 = p1
                    p1 = ((2.0d0 * j - 1.0d0) * z * p2 - (j - 1.0d0) * p3) / j
                end do

                ! Derivative of P_n
                pp = n * (z * p1 - p2) / (z*z - 1.0d0)

                z1 = z
                z = z1 - p1 / pp  ! Newton-Raphson step

                if (abs(z - z1) < EPS) exit
            end do

            x(i) = -z
            x(n + 1 - i) = z
            w(i) = 2.0d0 / ((1.0d0 - z*z) * pp * pp)
            w(n + 1 - i) = w(i)
        end do

    end subroutine compute_nodes_weights

    function calc_xbar(x, xp) result(xbar)

        implicit none

        real(dp), intent(in) :: x, xp
        real(dp) :: xbar

        xbar = 0.5d0 * (xp + x)

    end function calc_xbar

end module gauss_quadrature_m
