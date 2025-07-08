! integrands for both Krook and FP collision models
module electrostatic_integrands

    use KIM_kinds, only: dp

    implicit none

    type :: integration_point_t
        real(dp) :: rhoT
        integer :: j
        real(dp) :: xlm1, xlp1, xl
        real(dp) :: xlpm1, xlpp1, xlp
        real(dp) :: a_coef, b_coef
        real(dp) :: Jrg1, Jrg2, Jrg3, Jrg4
        contains
            procedure :: calc_Jrg1
            procedure :: calc_Jrg2
            procedure :: calc_Jrg3
            procedure :: calc_Jrg4
    end type

    type :: int_F0_rho_phi_t
        type(integration_point_t) :: int_point
        contains
            procedure :: f => integrand_F0_rho_phi
    end type

    type :: int_F1_rho_phi_t
        type(integration_point_t) :: int_point
        contains
            procedure :: f => integrand_F1_rho_phi
    end type

    type :: int_F2_rho_phi_t
        type(integration_point_t) :: int_point
        contains
            procedure :: f => integrand_F2_rho_phi
    end type

    type :: int_F3_rho_phi_t
        type(integration_point_t) :: int_point
        contains
            procedure :: f => integrand_F3_rho_phi
    end type

    contains

    function integrand_F0_rho_phi(this, x) result(val)

        use grid, only: rg_grid
        use functions, only: varphi_l
        use gsl_mod, only: erf => gsl_sf_erf
        use constants, only: pi

        implicit none

        real(dp), intent(in) :: x
        real(dp) :: val
        class(int_F0_rho_phi_t), intent(in) :: this

        val = varphi_l(x, this%int_point%xlm1, this%int_point%xl, this%int_point%xlp1) &
            * varphi_l(x, this%int_point%xlpm1, this%int_point%xlp, this%int_point%xlpp1) &
            * (&
                  erf((rg_grid%xb(this%int_point%j+1)-x)/(sqrt(2.0d0) * abs(this%int_point%rhoT))) &
                - erf((rg_grid%xb(this%int_point%j) - x)/(sqrt(2.0d0) * abs(this%int_point%rhoT)))&
            ) &
            * 2.0d0 * pi**2.0d0

    end function

    function integrand_F1_rho_phi(this, x, xp, theta) result(val)

        use constants, only: pi
        use species, only: plasma
        use gsl_mod, only: erf => gsl_sf_erf
        use KIM_kinds, only: dp
        use functions, only: varphi_l

        implicit none

        class(int_F1_rho_phi_t), intent(inout) :: this
        real(dp), intent(in) :: x, xp, theta
        real(dp) :: val
        real(dp) :: ks_val

        ks_val = 0.5d0 * (plasma%ks(this%int_point%j) + plasma%ks(this%int_point%j+1))

        val = varphi_l(xp, this%int_point%xlpm1, this%int_point%xlp, this%int_point%xlpp1) &
            * varphi_l(x, this%int_point%xlm1, this%int_point%xl, this%int_point%xlp1) &
            * 2.0d0 * pi / (this%int_point%rhoT**2.0d0 * sin(theta)) &
            * exp(- ks_val**2.0d0 * this%int_point%rhoT**2.0d0 &
                  - (x - xp)**2.0d0 / (4.0d0 * this%int_point%rhoT**2.0d0 * (1.0d0 - cos(theta)))) &
            * this%int_point%Jrg1

    end function


    function integrand_F2_rho_phi(this, x, xp, theta) result(val)

        use constants, only: pi
        use species, only: plasma
        use gsl_mod, only: erf => gsl_sf_erf
        use KIM_kinds, only: dp
        use functions, only: varphi_l

        implicit none

        class(int_F2_rho_phi_t), intent(inout) :: this
        real(dp), intent(in) :: x, xp, theta
        real(dp) :: val
        real(dp) :: ks_val

        ks_val = 0.5d0 * (plasma%ks(this%int_point%j) + plasma%ks(this%int_point%j+1))

        val = varphi_l(xp, this%int_point%xlpm1, this%int_point%xlp, this%int_point%xlpp1) &
            * varphi_l(x, this%int_point%xlm1, this%int_point%xl, this%int_point%xlp1) &
            * (-pi) / (4.0d0 * this%int_point%rhoT**4.0d0 * sin(theta)**5.0d0) &
            * exp(- ks_val**2.0d0 * this%int_point%rhoT**2.0d0 &
                  - (x - xp)**2.0d0 / (4.0d0 * this%int_point%rhoT**2.0d0 * (1.0d0 - cos(theta)))) & 
            * ( &
                this%int_point%Jrg1 * (&
                    4.0d0 * cos(2.0d0 * theta) * this%int_point%rhoT**2.0d0 *(ks_val**2.0d0 * this%int_point%rhoT**2.0d0 + 1.0d0) &
                    - ks_val**2.0d0 * this%int_point%rhoT**4.0d0 * cos(4.0d0 * theta)  &
                    - this%int_point%rhoT**2.0d0 * (3.0d0 * ks_val**2.0d0 * this%int_point%rhoT**2.0d0 + 4.0d0) &
                ) &
                + (2.0d0 * cos(2.0d0 * theta) + 6.0d0) * (this%int_point%Jrg2 + this%int_point%Jrg3) &
                - 16.0d0 * cos(theta) * this%int_point%Jrg4 &
            )

    end function


    function integrand_F3_rho_phi(this, x, xp, theta) result(val)

        use constants, only: pi
        use species, only: plasma
        use gsl_mod, only: erf => gsl_sf_erf
        use KIM_kinds, only: dp
        use functions, only: varphi_l

        implicit none

        class(int_F3_rho_phi_t), intent(inout) :: this
        real(dp), intent(in) :: x, xp, theta
        real(dp) :: val
        real(dp) :: ks_val

        ks_val = 0.5d0 * (plasma%ks(this%int_point%j) + plasma%ks(this%int_point%j+1))

        val = varphi_l(xp, this%int_point%xlpm1, this%int_point%xlp, this%int_point%xlpp1) &
            * varphi_l(x, this%int_point%xlm1, this%int_point%xl, this%int_point%xlp1) &
            * (-pi) * cos(theta) / ( 2.0d0 * this%int_point%rhoT**4.0d0 * sin(theta)**5.0d0) &
            * exp(- ks_val**2.0d0 * this%int_point%rhoT**2.0d0 &
                  - (x - xp)**2.0d0 / (4.0d0 * this%int_point%rhoT**2.0d0 * (1.0d0 - cos(theta)))) & 
            * ( &
                this%int_point%rhoT**2.0d0 * (cos(3.0d0 * theta) - cos(theta)) * this%int_point%Jrg1 &
                + 4.0d0 * cos(theta) * (this%int_point%Jrg2 + this%int_point%Jrg3) &
                - 2.0d0 * (cos(2.0d0 * theta) + 3.0d0) * this%int_point%Jrg4 &
            )

    end function

    subroutine calc_Jrg1(this)

        use constants, only: pi
        use grid, only: rg_grid
        use gsl_mod, only: erf => gsl_sf_erf

        implicit none

        class(integration_point_t), intent(inout) :: this

        this%Jrg1 = sqrt(pi) / (2.0d0 * this%a_coef) &
            *(erf(this%a_coef * (this%b_coef - rg_grid%xb(this%j))) &
            - erf(this%a_coef * (this%b_coef - rg_grid%xb(this%j+1))))

    end subroutine

    subroutine calc_Jrg2(this)

        use constants, only: pi
        use grid, only: rg_grid
        use gsl_mod, only: erf => gsl_sf_erf

        implicit none

        class(integration_point_t), intent(inout) :: this

        this%Jrg2 = 1.0d0 / (4.0d0 * this%a_coef**3.0d0) &
                    * ( &
                        sqrt(pi) * (2.0d0 * this%a_coef**2.0d0 * (this%b_coef - this%xl)**2.0d0 + 1.0d0) &
                            * (erf(this%a_coef * (this%b_coef - rg_grid%xb(this%j))) &
                                - erf(this%a_coef * (this%b_coef - rg_grid%xb(this%j+1)))) &
                        + 2.0d0 * this%a_coef * exp(-this%a_coef**2.0d0 * (this%b_coef - rg_grid%xb(this%j))**2.0d0) &
                            * (this%b_coef + rg_grid%xb(this%j) - 2.0d0 * this%xl) &
                        - 2.0d0 * this%a_coef * exp(-this%a_coef**2.0d0 * (this%b_coef - rg_grid%xb(this%j+1))**2.0d0) &
                            * (this%b_coef + rg_grid%xb(this%j+1) - 2.0d0 * this%xl) &
                    )

    end subroutine


    subroutine calc_Jrg3(this)

        use constants, only: pi
        use grid, only: rg_grid
        use gsl_mod, only: erf => gsl_sf_erf

        implicit none

        class(integration_point_t), intent(inout) :: this

        this%Jrg3 = 1.0d0 / (4.0d0 * this%a_coef**3.0d0) &
                    * ( &
                        sqrt(pi) * (2.0d0 * this%a_coef**2.0d0 * (this%b_coef - this%xlp)**2.0d0 + 1.0d0) &
                            * (erf(this%a_coef * (this%b_coef - rg_grid%xb(this%j))) &
                                - erf(this%a_coef * (this%b_coef - rg_grid%xb(this%j+1)))) &
                        + 2.0d0 * this%a_coef * exp(-this%a_coef**2.0d0 * (this%b_coef - rg_grid%xb(this%j))**2.0d0) &
                            * (this%b_coef + rg_grid%xb(this%j) - 2.0d0 * this%xlp) &
                        - 2.0d0 * this%a_coef * exp(-this%a_coef**2.0d0 * (this%b_coef - rg_grid%xb(this%j+1))**2.0d0) &
                            * (this%b_coef + rg_grid%xb(this%j+1) - 2.0d0 * this%xlp) & 
                    )

    end subroutine

    subroutine calc_Jrg4(this)

        use constants, only: pi
        use grid, only: rg_grid
        use gsl_mod, only: erf => gsl_sf_erf

        implicit none

        class(integration_point_t), intent(inout) :: this

        this%Jrg4 = 1.0d0 / (4.0d0 * this%a_coef**3.0d0) &
                    * ( &
                        (erf(this%a_coef * (this%b_coef - rg_grid%xb(this%j))) &
                            - erf(this%a_coef * (this%b_coef - rg_grid%xb(this%j+1)))) &
                            * sqrt(pi) * (2.0d0 * this%a_coef**2.0d0 * (this%b_coef - this%xl)*(this%b_coef - this%xlp)+1.0d0) &
                        + 2.0d0 * this%a_coef * exp(-this%a_coef**2.0d0 * (this%b_coef - rg_grid%xb(this%j))**2.0d0) &
                            * (this%b_coef + rg_grid%xb(this%j) - this%xl - this%xlp) &
                        + 2.0d0 * this%a_coef * exp(-this%a_coef**2.0d0 * (this%b_coef - rg_grid%xb(this%j+1))**2.0d0) &
                            * (-this%b_coef - rg_grid%xb(this%j+1) + this%xl + this%xlp) &
                    )

    end subroutine

    function calc_b_coef(x,xp) result(b_coef)

        implicit none

        real(dp), intent(in) :: x, xp
        real(dp) :: b_coef

        b_coef = 0.5d0 * (xp + x)

    end function

end module