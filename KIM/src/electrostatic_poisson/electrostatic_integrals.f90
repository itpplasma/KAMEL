module electrostatic_integrals_gauss_mod
    ! module for integrals of the electrostatic problem
    ! uses gauss legendre quadrature for integration

    use KIM_kinds_m, only: dp

    implicit none

    type :: gauss_config_t
        integer :: Nx, Nxp, Ntheta ! Number of nodes in each dimension
        real(dp), allocatable, dimension(:) :: x_x, w_x
        real(dp), allocatable, dimension(:) :: x_theta, w_theta
        real(dp), allocatable, dimension(:) :: x_xp, w_xp
    end type gauss_config_t

    contains

    subroutine init_gauss_int(gauss_conf)

        implicit none

        type(gauss_config_t), intent(inout) :: gauss_conf

        allocate(gauss_conf%x_x(gauss_conf%Nx), gauss_conf%w_x(gauss_conf%Nx))
        call compute_nodes_weights(gauss_conf%Nx, gauss_conf%x_x, gauss_conf%w_x)

        allocate(gauss_conf%x_xp(gauss_conf%Nxp), gauss_conf%w_xp(gauss_conf%Nxp))
        call compute_nodes_weights(gauss_conf%Nxp, gauss_conf%x_xp, gauss_conf%w_xp)

        allocate(gauss_conf%x_theta(gauss_conf%Ntheta), gauss_conf%w_theta(gauss_conf%Ntheta))
        call compute_nodes_weights(gauss_conf%Ntheta, gauss_conf%x_theta, gauss_conf%w_theta)

    end subroutine

    subroutine gauss_integrate_F0(int_F, a, b, result, gauss_conf)

        use electrostatic_integrands_gauss_mod, only: gauss_int_F0_rho_phi_t

        implicit none

        real(dp), intent(in) :: a, b
        type(gauss_config_t), intent(in) :: gauss_conf
        real(dp), intent(out) :: result
        type(gauss_int_F0_rho_phi_t), intent(in) :: int_F

        integer :: i
        real(dp) :: xm, xr

        result = 0.0d0
        xm = 0.5d0 * (b + a)
        xr = 0.5d0 * (b - a)

        do i = 1, gauss_conf%Nx
            result = result + gauss_conf%w_x(i) * int_F%f(xr * gauss_conf%x_x(i) + xm)
        end do
        result = result * xr

    end subroutine

    subroutine gauss_integrate_F1(int_F1, result, gauss_conf)

        use electrostatic_integrands_gauss_mod, only: gauss_int_F1_rho_phi_t, calc_b_coef
        use constants_m, only: pi
        use config_m, only: output_path

        implicit none

        class(gauss_int_F1_rho_phi_t), intent(inout) :: int_F1

        type(gauss_config_t), intent(in) :: gauss_conf
        real(dp), intent(out) :: result
        real(dp) :: x_mapped, xp_mapped, theta_mapped
        real(dp) :: norm_factor
        integer :: i,j,k
        integer :: iunit
        logical :: first_call = .true.
        save :: first_call

        result = 0.0d0

        norm_factor = pi * (int_F1%int_point%xlp1 - int_F1%int_point%xlm1) & ! normalization due to integral range shift
                        * (int_F1%int_point%xlpp1 - int_F1%int_point%xlpm1) / 8.0d0
        
        ! Write integration nodes to file (only on first call)
        if (first_call) then
            first_call = .false.
            iunit = 1234
            open(unit=iunit, file=trim(output_path)//'grid/F1_integration_nodes.txt', status='replace')
            write(iunit, '(A)') '# F1 Integration Nodes - First call information'
            write(iunit, '(A,I4)') '# Number of theta nodes: ', gauss_conf%Ntheta
            write(iunit, '(A,I4)') '# Number of x nodes: ', gauss_conf%Nx
            write(iunit, '(A,I4)') '# Number of xp nodes: ', gauss_conf%Nxp
            write(iunit, '(A)') '# j index, xl boundaries: xlm1, xl, xlp1, xlpm1, xlp, xlpp1, rhoT'
            write(iunit, '(I4, 7F12.6)') int_F1%int_point%j, int_F1%int_point%xlm1, int_F1%int_point%xl, &
                                        int_F1%int_point%xlp1, int_F1%int_point%xlpm1, int_F1%int_point%xlp, &
                                        int_F1%int_point%xlpp1, int_F1%int_point%rhoT
            write(iunit, '(A)') '#'
            write(iunit, '(A)') '# Theta nodes (Gauss-Legendre on [-1,1] and mapped values):'
            write(iunit, '(A)') '# i, GL_node, theta_mapped, weight'
            do i = 1, gauss_conf%Ntheta
                theta_mapped = 0.5d0 * pi * (gauss_conf%x_theta(i) + 1.0d0)
                write(iunit, '(I4, 3F16.10)') i, gauss_conf%x_theta(i), theta_mapped, gauss_conf%w_theta(i)
            end do
            write(iunit, '(A)') '#'
            write(iunit, '(A)') '# X nodes (Gauss-Legendre on [-1,1] and mapped values for first iteration):'
            write(iunit, '(A)') '# k, GL_node, x_mapped, weight'
            do k = 1, gauss_conf%Nx
                x_mapped = 0.5d0 * ((int_F1%int_point%xlp1 - int_F1%int_point%xlm1) * gauss_conf%x_x(k) + &
                    int_F1%int_point%xlp1 + int_F1%int_point%xlm1)
                write(iunit, '(I4, 3F16.10)') k, gauss_conf%x_x(k), x_mapped, gauss_conf%w_x(k)
            end do
            write(iunit, '(A)') '#'
            write(iunit, '(A)') '# Xp nodes (Gauss-Legendre on [-1,1] and mapped values for first iteration):'
            write(iunit, '(A)') '# j, GL_node, xp_mapped, weight'
            do j = 1, gauss_conf%Nxp
                xp_mapped = 0.5d0 * ((int_F1%int_point%xlpp1 - int_F1%int_point%xlpm1) * gauss_conf%x_xp(j) + &
                    int_F1%int_point%xlpp1 + int_F1%int_point%xlpm1)
                write(iunit, '(I4, 3F16.10)') j, gauss_conf%x_xp(j), xp_mapped, gauss_conf%w_xp(j)
            end do
            write(iunit, '(A)') '#'
            write(iunit, '(A)') '# Integration domain information:'
            write(iunit, '(A,2F12.6)') '# x range: ', int_F1%int_point%xlm1, int_F1%int_point%xlp1
            write(iunit, '(A,2F12.6)') '# xp range: ', int_F1%int_point%xlpm1, int_F1%int_point%xlpp1
            write(iunit, '(A,F12.6)') '# Norm factor: ', norm_factor
            close(iunit)
            
            ! Also write a detailed sampling file for all integration points
            open(unit=iunit+1, file=trim(output_path)//'grid/F1_integration_sampling.txt', status='replace')
            write(iunit+1, '(A)') '# F1 Integration sampling points for first call'
            write(iunit+1, '(A)') '# i_theta, j_xp, k_x, theta, xp, x, integrand_value'
            close(iunit+1)
        end if
        
        do i=1,gauss_conf%Ntheta ! theta
            theta_mapped = 0.5d0 * pi * (gauss_conf%x_theta(i) + 1.0d0)
            int_F1%int_point%a_coef = sqrt(1.0d0 / (1.0d0 + cos(theta_mapped))) / abs(int_F1%int_point%rhoT)

            do j=1,gauss_conf%Nxp ! xp 
                xp_mapped = 0.5d0 * ((int_F1%int_point%xlpp1 - int_F1%int_point%xlpm1) * gauss_conf%x_xp(j) + &
                    int_F1%int_point%xlpp1 + int_F1%int_point%xlpm1)

                do k=1,gauss_conf%Nx !x
                    x_mapped = 0.5d0 * ((int_F1%int_point%xlp1 - int_F1%int_point%xlm1) * gauss_conf%x_x(k) + &
                        int_F1%int_point%xlp1 + int_F1%int_point%xlm1)

                    int_F1%int_point%b_coef = calc_b_coef(x_mapped, xp_mapped)

                    call int_F1%int_point%calc_Jrg1()

                    result = result + gauss_conf%w_theta(i) * gauss_conf%w_xp(j) * gauss_conf%w_x(k) &
                        * int_F1%f(x_mapped, xp_mapped, theta_mapped) &
                        * norm_factor
                end do
            end do
        end do

    end subroutine

    subroutine gauss_integrate_F2(int_F2, result, gauss_conf)
    
        use electrostatic_integrands_gauss_mod, only: gauss_int_F2_rho_phi_t, calc_b_coef
        use constants_m, only: pi

        implicit none

        class(gauss_int_F2_rho_phi_t), intent(inout) :: int_F2

        type(gauss_config_t), intent(in) :: gauss_conf
        real(dp), intent(out) :: result
        real(dp) :: x_mapped, xp_mapped, theta_mapped
        real(dp) :: norm_factor
        integer :: i,j,k

        result = 0.0d0

        norm_factor = pi * (int_F2%int_point%xlp1 - int_F2%int_point%xlm1) & ! normalization due to integral range shift
                        * (int_F2%int_point%xlpp1 - int_F2%int_point%xlpm1) / 8.0d0

        do i=1,gauss_conf%Ntheta ! theta
            theta_mapped = 0.5d0 * pi *(gauss_conf%x_theta(i) + 1.0d0)
            int_F2%int_point%a_coef = sqrt(1.0d0 / (1.0d0 + cos(theta_mapped))) / abs(int_F2%int_point%rhoT)

            do j=1,gauss_conf%Nxp ! xp 
                xp_mapped = 0.5d0 * ((int_F2%int_point%xlpp1 - int_F2%int_point%xlpm1) * gauss_conf%x_xp(j) + &
                    int_F2%int_point%xlpp1 + int_F2%int_point%xlpm1)

                do k=1,gauss_conf%Nx !x
                    x_mapped = 0.5d0 * ((int_F2%int_point%xlp1 - int_F2%int_point%xlm1) * gauss_conf%x_x(k) + &
                        int_F2%int_point%xlp1 + int_F2%int_point%xlm1)

                    int_F2%int_point%b_coef = calc_b_coef(x_mapped, xp_mapped)
                    int_F2%int_point%xl_mapped = x_mapped
                    int_F2%int_point%xlp_mapped = xp_mapped
        
                    call int_F2%int_point%calc_Jrg1()
                    call int_F2%int_point%calc_Jrg2()
                    call int_F2%int_point%calc_Jrg3()
                    call int_F2%int_point%calc_Jrg4()

                    result = result + gauss_conf%w_theta(i) * gauss_conf%w_xp(j) * gauss_conf%w_x(k) &
                        * int_F2%f(x_mapped, xp_mapped, theta_mapped) &
                        * norm_factor
                end do
            end do
        end do
    end subroutine

    subroutine gauss_integrate_F3(int_F3, result, gauss_conf)

        use electrostatic_integrands_gauss_mod, only: gauss_int_F3_rho_phi_t, calc_b_coef
        use constants_m, only: pi

        implicit none

        class(gauss_int_F3_rho_phi_t), intent(inout) :: int_F3

        type(gauss_config_t), intent(in) :: gauss_conf
        real(dp), intent(out) :: result
        real(dp) :: norm_factor
        real(dp) :: x_mapped, xp_mapped, theta_mapped
        integer :: i,j,k

        result = 0.0d0

        norm_factor = pi * (int_F3%int_point%xlp1 - int_F3%int_point%xlm1) & ! normalization due to integral range shift
                        * (int_F3%int_point%xlpp1 - int_F3%int_point%xlpm1) / 8.0d0

        do i=1,gauss_conf%Ntheta ! theta
            theta_mapped = 0.5d0 * pi * (gauss_conf%x_theta(i) + 1.0d0)
            int_F3%int_point%a_coef = sqrt(1.0d0 / (1.0d0 + cos(theta_mapped))) / abs(int_F3%int_point%rhoT)

            do j=1,gauss_conf%Nxp ! xp 
                xp_mapped = 0.5d0 * ((int_F3%int_point%xlpp1 - int_F3%int_point%xlpm1) * gauss_conf%x_xp(j) + &
                    int_F3%int_point%xlpp1 + int_F3%int_point%xlpm1)

                do k=1,gauss_conf%Nx !x
                    x_mapped = 0.5d0 * ((int_F3%int_point%xlp1 - int_F3%int_point%xlm1) * gauss_conf%x_x(k) + &
                        int_F3%int_point%xlp1 + int_F3%int_point%xlm1)

                    int_F3%int_point%b_coef = calc_b_coef(x_mapped, xp_mapped)
                    int_F3%int_point%xl_mapped = x_mapped
                    int_F3%int_point%xlp_mapped = xp_mapped

                    call int_F3%int_point%calc_Jrg1()
                    call int_F3%int_point%calc_Jrg2()
                    call int_F3%int_point%calc_Jrg3()
                    call int_F3%int_point%calc_Jrg4()

                    result = result + gauss_conf%w_theta(i) * gauss_conf%w_xp(j) * gauss_conf%w_x(k) &
                        * int_F3%f(x_mapped, xp_mapped, theta_mapped) &
                        * norm_factor
                end do
            end do
        end do 

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

end module electrostatic_integrals_gauss_mod
