! integrands for both Krook and FP collision models
module integrands_gauss_m

    use KIM_kinds_m, only: dp

    implicit none

    type :: integration_point_t
        real(dp) :: rhoT
        real(dp) :: ks
        integer :: j, mphi
        real(dp) :: xlm1, xlp1, xl
        real(dp) :: xlpm1, xlpp1, xlp
        real(dp) :: xl_mapped, xlp_mapped
        real(dp) :: a_coef, xbar
        real(dp) :: Jrg1, Jrg2, Jrg23, Jrg3, Jrg4
        contains
            procedure :: calc_Jrg1
            procedure :: calc_Jrg23
            procedure :: calc_Jrg4
    end type

    type :: gauss_int_F0_rho_phi_t
        type(integration_point_t) :: int_point
        contains
            procedure :: f => gauss_integrand_F0_rho_phi
    end type

    type :: gauss_int_F1_rho_phi_t
        type(integration_point_t) :: int_point
        contains
            procedure :: f => gauss_integrand_F1_rho_phi
    end type

    type :: gauss_int_F2_rho_phi_t
        type(integration_point_t) :: int_point
        contains
            procedure :: f => gauss_integrand_F2_rho_phi
    end type

    type :: gauss_int_F3_rho_phi_t
        type(integration_point_t) :: int_point
        contains
            procedure :: f => gauss_integrand_F3_rho_phi
    end type

    contains

    function gauss_integrand_F0_rho_phi(this, x) result(val)

        use grid_m, only: rg_grid
        use functions_m, only: varphi_l
        use constants_m, only: pi
        use numerics_utils_m, only: erf_diff

        implicit none

        real(dp), intent(in) :: x
        real(dp) :: val
        class(gauss_int_F0_rho_phi_t), intent(in) :: this

        val = varphi_l(x, this%int_point%xlm1, this%int_point%xl, this%int_point%xlp1) &
            * varphi_l(x, this%int_point%xlpm1, this%int_point%xlp, this%int_point%xlpp1) &
            * (&
                  erf_diff((rg_grid%xb(this%int_point%j+1)-x)/(sqrt(2.0d0) * abs(this%int_point%rhoT)), &
                (rg_grid%xb(this%int_point%j) - x)/(sqrt(2.0d0) * abs(this%int_point%rhoT)))&
            ) &
            * pi**2.0d0 ! from benchmarking, error function difference gives 2
            

    end function

    function gauss_integrand_F1_rho_phi(this, x, xp, theta) result(val)

        use constants_m, only: pi
        use species_m, only: plasma
        use gsl_mod, only: erf => gsl_sf_erf
        use KIM_kinds_m, only: dp
        use functions_m, only: varphi_l
        use numerics_utils_m, only: erf_diff
        use grid_m, only: rg_grid

        implicit none

        class(gauss_int_F1_rho_phi_t), intent(inout) :: this
        real(dp), intent(in) :: x, xp, theta
        real(dp) :: val
        real(dp) :: ks_val

        ks_val = 0.5d0 * (plasma%ks(this%int_point%j) + plasma%ks(this%int_point%j+1))

        val = varphi_l(xp, this%int_point%xlpm1, this%int_point%xlp, this%int_point%xlpp1) &
            * varphi_l(x, this%int_point%xlm1, this%int_point%xl, this%int_point%xlp1) &
            * pi ** 1.5d0 * (1.0d0 + cos(theta)) * cos(this%int_point%mphi * theta) / (this%int_point%rhoT * sin(theta)) &
            * exp(- ks_val**2.0d0 * this%int_point%rhoT**2.0d0 &
                  - (x - xp)**2.0d0 / (4.0d0 * this%int_point%rhoT**2.0d0 * (1.0d0 - cos(theta)))) &
            * erf_diff( &
                ((x+xp)/ 2.0d0 - rg_grid%xb(this%int_point%j)) / (this%int_point%rhoT * sqrt(1.0d0 + cos(theta))), &
                ((x+xp)/ 2.0d0 - rg_grid%xb(this%int_point%j + 1)) / (this%int_point%rhoT * sqrt(1.0d0 + cos(theta))) &
                )& 
            / (2.0d0 * pi) ! is missing, found by benchmarking


    end function


    function gauss_integrand_F2_rho_phi(this, x, xp, theta) result(val)

        use constants_m, only: pi
        use species_m, only: plasma
        use gsl_mod, only: erf => gsl_sf_erf
        use KIM_kinds_m, only: dp
        use functions_m, only: varphi_l

        implicit none

        class(gauss_int_F2_rho_phi_t), intent(inout) :: this
        real(dp), intent(in) :: x, xp, theta
        real(dp) :: val
        real(dp) :: ks_val

        ks_val = 0.5d0 * (plasma%ks(this%int_point%j) + plasma%ks(this%int_point%j+1))

        val = varphi_l(xp, this%int_point%xlpm1, this%int_point%xlp, this%int_point%xlpp1) &
            * varphi_l(x, this%int_point%xlm1, this%int_point%xl, this%int_point%xlp1) &
            * exp(- ks_val**2.0d0 * this%int_point%rhoT**2.0d0 &
                  - (x - xp)**2.0d0 / (4.0d0 * this%int_point%rhoT**2.0d0 * (1.0d0 - cos(theta)))) & 
            * pi * cos(this%int_point%mphi * theta) / (this%int_point%rhoT * sin(theta)**5.0d0) &
            * ( &
                this%int_point%Jrg1 * (&
                    (ks_val**2.0d0 * this%int_point%rhoT**2.0d0 + 2.0d0) &
                    * sin(theta)**2.0d0 &
                ) &
                - (cos(theta)**2.0d0 + 1.0d0) * this%int_point%Jrg23 &
                + 4.0d0 * cos(theta) * this%int_point%Jrg4 &
            ) &
            / (2.0d0 * pi) ! is missing, found by benchmarking

    end function


    function gauss_integrand_F3_rho_phi(this, x, xp, theta) result(val)

        use constants_m, only: pi
        use species_m, only: plasma
        use gsl_mod, only: erf => gsl_sf_erf
        use KIM_kinds_m, only: dp
        use functions_m, only: varphi_l

        implicit none

        class(gauss_int_F3_rho_phi_t), intent(inout) :: this
        real(dp), intent(in) :: x, xp, theta
        real(dp) :: val
        real(dp) :: ks_val

        ks_val = 0.5d0 * (plasma%ks(this%int_point%j) + plasma%ks(this%int_point%j+1))

        val = varphi_l(xp, this%int_point%xlpm1, this%int_point%xlp, this%int_point%xlpp1) &
            * varphi_l(x, this%int_point%xlm1, this%int_point%xl, this%int_point%xlp1) &
            * exp(- ks_val**2.0d0 * this%int_point%rhoT**2.0d0 &
                  - (x - xp)**2.0d0 / (4.0d0 * this%int_point%rhoT**2.0d0 * (1.0d0 - cos(theta)))) & 
            * (2.0d0 * pi) * cos((this%int_point%mphi - 1) * theta) / (this%int_point%rhoT * sin(theta)**5.0d0) &
            * ( &
                cos(theta) * sin(theta)**2.0d0 * this%int_point%Jrg1 &
                - cos(theta) * this%int_point%Jrg23 &
                + (cos(theta)**2.0d0 + 1.0d0) * this%int_point%Jrg4 &
            ) &
            / (2.0d0 * pi) ! is missing, found by benchmarking


    end function

    subroutine calc_Jrg1(this, theta, x, xp)

        use constants_m, only: pi
        use grid_m, only: rg_grid
        use numerics_utils_m, only: erf_diff

        implicit none

        class(integration_point_t), intent(inout) :: this
        real(dp), intent(in) :: theta
        real(dp), intent(in) :: x, xp

        this%Jrg1 = sqrt(pi) / 2.0d0 * sqrt(1.0d0 + cos(theta)) &
            *(&
                erf_diff(((x+xp)/2.0d0 - rg_grid%xb(this%j)) / (this%rhoT * sqrt(1.0d0 + cos(theta))), &
                         ((x+xp)/2.0d0 - rg_grid%xb(this%j+1)) / (this%rhoT * sqrt(1.0d0 + cos(theta)))) &
            )

    end subroutine


    subroutine calc_Jrg23(this, theta, x, xp)

        use constants_m, only: pi
        use grid_m, only: rg_grid
        use numerics_utils_m, only: erf_diff

        implicit none

        class(integration_point_t), intent(inout) :: this
        real(dp), intent(in) :: theta
        real(dp), intent(in) :: x, xp

        ! this%Jrg23 = this%rhoT**3.0d0 * (1.0d0 + cos(theta))**1.5d0 * &
        this%Jrg23 = (1.0d0 + cos(theta))**1.5d0 * &
            ( &
                sqrt(pi) * &
                    ( &
                        (x-xp)**2.0d0 / (4.0d0 * this%rhoT**2.0d0 * (1.0d0 + cos(theta))) &
                        + 0.5d0 &
                    ) &
                    * erf_diff( &
                        ((x+xp)/ 2.0d0 - rg_grid%xb(this%j)) / (this%rhoT * sqrt(1.0d0 + cos(theta))), &
                        ((x+xp)/ 2.0d0 - rg_grid%xb(this%j + 1)) / (this%rhoT * sqrt(1.0d0 + cos(theta))) &
                        ) &
                + exp(- ((x+xp)/2 - rg_grid%xb(this%j))**2.0d0 / (this%rhoT**2.0d0 * (1.0d0 + cos(theta)))) &
                    * ( rg_grid%xb(this%j) - (x+xp)/2)/ (this%rhoT * sqrt(1.0d0 + cos(theta))) &
                - exp(- ((x+xp)/2 - rg_grid%xb(this%j+1))**2.0d0 / (this%rhoT**2.0d0 * (1.0d0 + cos(theta)))) &
                    * ( rg_grid%xb(this%j+1) - (x+xp)/2)/ (this%rhoT * sqrt(1.0d0 + cos(theta))) &
            )

    end subroutine

    subroutine calc_Jrg4(this, theta, x, xp)

        use constants_m, only: pi
        use grid_m, only: rg_grid
        use numerics_utils_m, only: erf_diff

        implicit none

        class(integration_point_t), intent(inout) :: this
        real(dp), intent(in) :: theta
        real(dp), intent(in) :: x, xp

        ! this%Jrg4 = 0.5d0 * this%rhoT**3.0d0 * (1.0d0 + cos(theta))**1.5d0 * &
        this%Jrg4 = 0.5d0 * (1.0d0 + cos(theta))**1.5d0 * &
            ( &
                sqrt(pi) * &
                    ( &
                        0.5d0 - &
                        (x-xp)**2.0d0 / (4.0d0 * this%rhoT**2.0d0 * (1.0d0 + cos(theta))) &
                    ) &
                    * erf_diff( &
                        ((x+xp)/ 2.0d0 - rg_grid%xb(this%j)) / (this%rhoT * sqrt(1.0d0 + cos(theta))), &
                        ((x+xp)/ 2.0d0 - rg_grid%xb(this%j + 1)) / (this%rhoT * sqrt(1.0d0 + cos(theta))) &
                        ) &
                + exp(- ((x+xp)/2 - rg_grid%xb(this%j))**2.0d0 / (this%rhoT**2.0d0 * (1.0d0 + cos(theta)))) &
                    * ( rg_grid%xb(this%j) - (x+xp)/2)/ (this%rhoT * sqrt(1.0d0 + cos(theta))) &
                - exp(- ((x+xp)/2 - rg_grid%xb(this%j+1))**2.0d0 / (this%rhoT**2.0d0 * (1.0d0 + cos(theta)))) &
                    * ( rg_grid%xb(this%j+1) - (x+xp)/2)/ (this%rhoT * sqrt(1.0d0 + cos(theta))) &
            )

    end subroutine

    function calc_xbar(x,xp) result(xbar)

        implicit none

        real(dp), intent(in) :: x, xp
        real(dp) :: xbar

        xbar = 0.5d0 * (xp + x)

    end function

end module integrands_gauss_m
