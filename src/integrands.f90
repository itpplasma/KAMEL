module integrands
    ! Module containing the integration integrands. In particular, the kernel
    ! functions as well as the spline basis functions used for the basis
    ! transformation.

    implicit none

    contains


    double complex function integrand_w_exp_facs_rho_phi(l, lp, i_kr, i_krp, i_rg)

        use kernel_functions, only: kernel_rho_phi_of_kr_krp_rg
        use kernel, only: K_rho_phi_of_rg
        use constants, only: com_unit
        use grid, only: varphi_lkr, xl, kr, krp, rb

        implicit none
        integer, intent(in) :: l, lp
        integer, intent(in) :: i_kr, i_krp
        integer, intent(in) :: i_rg

        !double complex tilde_varphi_lkr

        !double complex :: kernel_rho_phi_of_kr_krp_rg
 
        integrand_w_exp_facs_rho_phi = 0.0d0

        ! return \tilde{\varphi}_l(k_r) \tilde{\varphi}_{l'}^*(k_r') K(k_r, k_r'; r_g) exp(i k_r(r_g - x_l)+i k_r'(x_l' - r_g)) 

        integrand_w_exp_facs_rho_phi = varphi_lkr(i_kr, l) &!rb(l), rb(l-1), rb(l+1), rb(l+2), kr(i_kr)) &
            * conjg(varphi_lkr(i_krp, lp)) &!, rb(lp-1), rb(lp+1), rb(lp+2), krp(i_krp))) &
            * exp(com_unit * kr(i_kr) * (rb(i_rg)-xl(l)) + com_unit * krp(i_krp) * (xl(lp)-rb(i_rg))) &
            * K_rho_phi_of_rg(i_krp, i_kr, i_rg)

    end function integrand_w_exp_facs_rho_phi

    double complex function integrand_w_exp_facs_rho_B(l, lp, i_kr, i_krp, i_rg)

        use kernel_functions, only: kernel_rho_B_of_kr_krp_rg
        use kernel, only: K_rho_B_of_rg
        use constants, only: com_unit
        use grid, only: varphi_lkr, xl, kr, krp, rb

        implicit none
        integer, intent(in) :: l, lp
        integer, intent(in) :: i_kr, i_krp
        integer, intent(in) :: i_rg

        !double complex tilde_varphi_lkr

        !double complex :: kernel_rho_phi_of_kr_krp_rg
 
        integrand_w_exp_facs_rho_B = 0.0d0

        ! return \tilde{\varphi}_l(k_r) \tilde{\varphi}_{l'}^*(k_r') K(k_r, k_r'; r_g) exp(i k_r(r_g - x_l)+i k_r'(x_l' - r_g)) 

        integrand_w_exp_facs_rho_B = varphi_lkr(i_kr, l) &!rb(l), rb(l-1), rb(l+1), rb(l+2), kr(i_kr)) &
            * conjg(varphi_lkr(i_krp, lp)) &!, rb(lp-1), rb(lp+1), rb(lp+2), krp(i_krp))) &
            * exp(com_unit * kr(i_kr) * (rb(i_rg)-xl(l)) + com_unit * krp(i_krp) * (xl(lp)-rb(i_rg))) &
            * K_rho_B_of_rg(i_krp, i_kr, i_rg)

    end function integrand_w_exp_facs_rho_B


end module