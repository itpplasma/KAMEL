module electrostatic_integrals_rkf45_mod
    ! module for integrals of the electrostatic problem
    ! uses gauss legendre quadrature for integration

    use KIM_kinds, only: dp
    use constants, only: pi

    implicit none

    type :: rkf45_config_t
        integer :: Nx, Nxp, Ntheta ! Number of nodes in each dimension
        real(dp), allocatable, dimension(:) :: x_x, w_x
        real(dp), allocatable, dimension(:) :: x_xp, w_xp
    end type rkf45_config_t

    real(dp) :: f0 = 0.0d0
    real(dp) :: theta_0 = 0.01d0 ! add small epsilon to avoid singularity
    real(dp) :: theta_max = pi - 0.01d0
    real(dp) :: h0 = 0.1d0

    contains

    subroutine init_rkf45_int(rkf45_conf)

        implicit none

        type(rkf45_config_t), intent(inout) :: rkf45_conf

        allocate(rkf45_conf%x_x(rkf45_conf%Nx), rkf45_conf%w_x(rkf45_conf%Nx))
        call compute_nodes_weights(rkf45_conf%Nx, rkf45_conf%x_x, rkf45_conf%w_x)

        allocate(rkf45_conf%x_xp(rkf45_conf%Nxp), rkf45_conf%w_xp(rkf45_conf%Nxp))
        call compute_nodes_weights(rkf45_conf%Nxp, rkf45_conf%x_xp, rkf45_conf%w_xp)

    end subroutine

    subroutine rkf45_integrate_F0(result, rkf45_conf, context)

        use electrostatic_integrands_rkf45_mod, only: rkf45_integrand_context_t, rkf45_integrand_F0

        implicit none

        type(rkf45_integrand_context_t), intent(inout) :: context
        type(rkf45_config_t), intent(in) :: rkf45_conf
        real(dp), intent(out) :: result

        integer :: i
        real(dp) :: xm, xr

        result = 0.0d0
        xm = 0.5d0 * (context%xlp1 + context%xlm1)
        xr = 0.5d0 * (context%xlp1 - context%xlm1)
        !$omp parallel do collapse(1) private(i) firstprivate(context) reduction(+:result)
        do i = 1, rkf45_conf%Nx
            context%x = xr * rkf45_conf%x_x(i) + xm
            result = result + rkf45_conf%w_x(i) * rkf45_integrand_F0(context)
        end do
        !$omp end parallel do
        result = result * xr

    end subroutine

    subroutine rkf45_integrate_F1(result, rkf45_conf, context)

        use electrostatic_integrands_rkf45_mod, only: rkf45_integrand_context_t, rkf45_integrand_F1
        use constants, only: pi
        use config, only: output_path
        use RKF45_mod, only: RKF45_1D_with_context
        use grid, only: rkf45_tol

        implicit none

        type(rkf45_integrand_context_t), intent(inout) :: context

        type(rkf45_config_t), intent(in) :: rkf45_conf
        real(dp), intent(out) :: result
        real(dp) :: norm_factor
        real(dp) :: rk45_res
        integer :: j,k

        result = 0.0d0
        norm_factor = (context%xlp1 - context%xlm1)* (context%xlpp1 - context%xlpm1) / 4.0d0 ! gauss integration normalizaton

        !$omp parallel do collapse(2) private(j, k, rk45_res) firstprivate(context) reduction(+:result)
        do j=1,rkf45_conf%Nxp ! xp 
            do k=1,rkf45_conf%Nx !x
                context%xp = 0.5d0 * ((context%xlpp1 - context%xlpm1) * rkf45_conf%x_xp(j) + context%xlpp1 + context%xlpm1)
                context%x = 0.5d0 * ((context%xlp1 - context%xlm1) * rkf45_conf%x_x(k) + context%xlp1 + context%xlm1)

                rk45_res = 0.0d0

                call RKF45_1D_with_context(rkf45_integrand_F1, f0, theta_0, theta_max, h0, rkf45_tol, rk45_res, context)

                result = result + rkf45_conf%w_xp(j) * rkf45_conf%w_x(k) * rk45_res
            end do
        end do
        !$omp end parallel do

        result = result * exp(- context%ks**2.0d0 * context%rhoT**2.0d0) * norm_factor

    end subroutine

    subroutine rkf45_integrate_F2(result, rkf45_conf, context)
    
        use constants, only: pi
        use RKF45_mod, only: RKF45_1D_with_context
        use electrostatic_integrands_rkf45_mod, only: rkf45_integrand_context_t, rkf45_integrand_F2
        use grid, only: rkf45_tol

        implicit none

        type(rkf45_integrand_context_t), intent(inout) :: context

        type(rkf45_config_t), intent(in) :: rkf45_conf
        real(dp), intent(out) :: result
        real(dp) :: norm_factor
        real(dp) :: rk45_res
        integer :: j,k

        result = 0.0d0
        norm_factor = (context%xlp1 - context%xlm1)* (context%xlpp1 - context%xlpm1) / 4.0d0 ! gauss integration normalizaton

        !$omp parallel do collapse(2) private(j, k, rk45_res) firstprivate(context) reduction(+:result)
        do j=1,rkf45_conf%Nxp ! xp 
            do k=1,rkf45_conf%Nx !x
                context%xp = 0.5d0 * ((context%xlpp1 - context%xlpm1) * rkf45_conf%x_xp(j) + context%xlpp1 + context%xlpm1)
                context%x = 0.5d0 * ((context%xlp1 - context%xlm1) * rkf45_conf%x_x(k) + context%xlp1 + context%xlm1)

                rk45_res = 0.0d0

                call RKF45_1D_with_context(rkf45_integrand_F2, f0, theta_0, theta_max, h0, rkf45_tol, rk45_res, context)

                result = result + rkf45_conf%w_xp(j) * rkf45_conf%w_x(k) * rk45_res
            end do
        end do
        !$omp end parallel do

        result = result * exp(- context%ks**2.0d0 * context%rhoT**2.0d0) * norm_factor &
                * (-pi) / (4.0d0 * context%rhoT**4.0d0)

    end subroutine

    subroutine rkf45_integrate_F3(result, rkf45_conf, context)
    
        use constants, only: pi
        use RKF45_mod, only: RKF45_1D_with_context
        use electrostatic_integrands_rkf45_mod, only: rkf45_integrand_context_t, rkf45_integrand_F3
        use grid, only: rkf45_tol

        implicit none

        type(rkf45_integrand_context_t), intent(inout) :: context
        type(rkf45_config_t), intent(in) :: rkf45_conf
        real(dp), intent(out) :: result
        real(dp) :: norm_factor
        real(dp) :: rk45_res
        integer :: j,k

        result = 0.0d0
        norm_factor = (context%xlp1 - context%xlm1)* (context%xlpp1 - context%xlpm1) / 4.0d0 ! gauss integration normalizaton

        !$omp parallel do collapse(2) private(j, k, rk45_res) firstprivate(context) reduction(+:result)
        do j=1,rkf45_conf%Nxp ! xp 
            do k=1,rkf45_conf%Nx !x
                context%xp = 0.5d0 * ((context%xlpp1 - context%xlpm1) * rkf45_conf%x_xp(j) + context%xlpp1 + context%xlpm1)
                context%x = 0.5d0 * ((context%xlp1 - context%xlm1) * rkf45_conf%x_x(k) + context%xlp1 + context%xlm1)

                rk45_res = 0.0d0

                call RKF45_1D_with_context(rkf45_integrand_F3, f0, theta_0, theta_max, h0, rkf45_tol, rk45_res, context)

                result = result + rkf45_conf%w_xp(j) * rkf45_conf%w_x(k) * rk45_res
            end do
        end do
        !$omp end parallel do

        result = result * exp(- context%ks**2.0d0 * context%rhoT**2.0d0) * norm_factor &
                * (-pi) / (2.0d0 * context%rhoT**4.0d0)

    end subroutine

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

end module electrostatic_integrals_rkf45_mod
