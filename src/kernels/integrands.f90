module integrands
    ! Module containing the integration integrands. In particular, the kernel
    ! functions as well as the spline basis functions used for the basis
    ! transformation.

    use kernels, only: K_rho_phi_of_rg, kernel_rho_phi_of_kr_krp_rg, &
                       K_rho_B_of_rg, kernel_rho_B_of_kr_krp_rg
    use constants, only: com_unit
    use grid, only: varphi_lkr, xl_grid, rg_grid, kr_grid, krp_grid

    implicit none

    contains


    double complex function integrand_w_exp_facs_rho_phi(l, lp, i_kr, i_krp, i_rg)

        implicit none
        integer, intent(in) :: l, lp
        integer, intent(in) :: i_kr, i_krp
        integer, intent(in) :: i_rg

        integrand_w_exp_facs_rho_phi = 0.0d0

        integrand_w_exp_facs_rho_phi = varphi_lkr(i_kr, l) * conjg(varphi_lkr(i_krp, lp)) &
            * exp(com_unit * kr_grid%xb(i_kr) * (rg_grid%xb(i_rg) - xl_grid%xb(l)) &
            + com_unit * krp_grid%xb(i_krp) * (xl_grid%xb(lp) - rg_grid%xb(i_rg))) &
            * K_rho_phi_of_rg(i_krp, i_kr, i_rg)

    end function integrand_w_exp_facs_rho_phi

    double complex function integrand_w_exp_facs_rho_B(l, lp, i_kr, i_krp, i_rg)
        
        implicit none
        integer, intent(in) :: l, lp
        integer, intent(in) :: i_kr, i_krp
        integer, intent(in) :: i_rg

        integrand_w_exp_facs_rho_B = 0.0d0

        integrand_w_exp_facs_rho_B = varphi_lkr(i_kr, l) * conjg(varphi_lkr(i_krp, lp)) &
            * exp(com_unit * kr_grid%xb(i_kr) * (rg_grid%xb(i_rg) - xl_grid%xb(l)) &
            + com_unit * krp_grid%xb(i_krp) * (xl_grid%xb(lp) - rg_grid%xb(i_rg))) &
            * K_rho_B_of_rg(i_krp, i_kr, i_rg)

    end function integrand_w_exp_facs_rho_B


end module