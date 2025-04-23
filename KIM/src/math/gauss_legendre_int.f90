module gauss_quad

    use KIM_kinds, only: dp

    implicit none

    type :: gauss_config_t
        integer :: n  ! Number of nodes
        real(dp), allocatable, dimension(:) :: x, w  ! Integration limits
    end type gauss_config_t

    contains

    subroutine init_gauss_int(gauss_conf)

        use KIM_kinds, only: dp

        implicit none

        type(gauss_config_t), intent(inout) :: gauss_conf

        allocate(gauss_conf%x(gauss_conf%n), gauss_conf%w(gauss_conf%n))
        call compute_nodes_weights(gauss_conf%n, gauss_conf%x, gauss_conf%w)

    end subroutine

    subroutine gauss_integrate_B0(int_B, a, b, result, gauss_conf)

        use KIM_kinds, only: dp
        use reduced_integrands, only: int_B0_rho_phi_t

        implicit none

        real(dp), intent(in) :: a, b
        type(gauss_config_t), intent(in) :: gauss_conf
        real(dp), intent(out) :: result
        type(int_B0_rho_phi_t), intent(in) :: int_B

        integer :: i
        real(dp) :: xm, xr

        result = 0.0d0
        xm = 0.5d0 * (b + a)
        xr = 0.5d0 * (b - a)

        do i = 1, gauss_conf%n
            result = result + gauss_conf%w(i) * int_B%f(xr * gauss_conf%x(i) + xm)
        end do
        result = result * xr

    end subroutine gauss_integrate_B0

    subroutine gauss_integrate_B1(int_B1, result, gauss_conf)

        use KIM_kinds, only: dp
        use reduced_integrands, only: int_B1_rho_phi_t
        use grid, only: xl_grid
        use constants, only: pi

        implicit none

        class(int_B1_rho_phi_t), intent(in) :: int_B1

        type(gauss_config_t), intent(in) :: gauss_conf
        real(dp), intent(out) :: result
        real(dp) :: x_mapped, xp_mapped, theta_mapped
        integer :: i,j,k


        do i=1,gauss_conf%n ! theta
            theta_mapped = 0.5d0 * (pi * gauss_conf%x(i) + pi)

            do j=1,gauss_conf%n ! xp 
                !xp_mapped = 0.5d0 * ((xl_grid%xb(int_B1%lp+1) - xl_grid%xb(int_B1%lp-1)) * gauss_conf%x(j) + &
                !    xl_grid%xb(int_B1%lp+1) + xl_grid%xb(int_B1%lp-1))
                xp_mapped = 0.5d0 * ((int_B1%xlpp1 - int_B1%xlpm1) * gauss_conf%x(j) + &
                    int_B1%xlpp1 + int_B1%xlpm1)


                do k=1,gauss_conf%n !x
                    !x_mapped = 0.5d0 * ((xl_grid%xb(int_B1%l+1) - xl_grid%xb(int_B1%l-1)) * gauss_conf%x(k) + &
                        !xl_grid%xb(int_B1%l+1) + xl_grid%xb(int_B1%l-1))
                    x_mapped = 0.5d0 * ((int_B1%xlp1 - int_B1%xlm1) * gauss_conf%x(j) + &
                        int_B1%xlp1 + int_B1%xlm1)

                    result = result + gauss_conf%w(i) * gauss_conf%w(j) * gauss_conf%w(k) &
                        * int_B1%f(x_mapped, xp_mapped, theta_mapped) &
                        * pi * (int_B1%xlp1 - int_B1%xlm1) & ! normalization due to integral range shift
                        * (int_B1%xlpp1 - int_B1%xlp-1) / 8.0d0
                        !* pi * (xl_grid%xb(int_B1%l+1) - xl_grid%xb(int_B1%l-1)) & ! normalization due to integral range shift
                        !* (xl_grid%xb(int_B1%lp+1) - xl_grid%xb(int_B1%lp-1)) / 8.0d0
                end do
            end do
        end do

    end subroutine

    subroutine compute_nodes_weights(n, x, w)

        use KIM_kinds, only: dp

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

end module gauss_quad
