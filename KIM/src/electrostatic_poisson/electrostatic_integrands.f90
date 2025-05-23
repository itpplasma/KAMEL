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


    type, extends(integration_point_t) :: int_F0_rho_phi_t
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
        real(dp) :: rhoT
        integer :: j
        real(dp) :: xlm1, xlp1, xl
        real(dp) :: xlpm1, xlpp1, xlp
        contains
            procedure :: f => integrand_F2_rho_phi
    end type

    type :: int_F3_rho_phi_t
        real(dp) :: rhoT
        integer :: j
        real(dp) :: xlm1, xlp1, xl
        real(dp) :: xlpm1, xlpp1, xlp
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
                erf((x - rg_grid%xb(this%int_point%j))/(sqrt(2.0d0) * this%int_point%rhoT)) &
                - erf((x - rg_grid%xb(this%int_point%j+1))/(sqrt(2.0d0) * this%int_point%rhoT))&
            ) &
            * 2.0d0 * pi**2.0d0

    end function

    function integrand_F1_rho_phi(this, x, xp, theta) result(val)

        use constants, only: pi
        use species, only: plasma
        use gsl_mod, only: erf => gsl_sf_erf
        use grid, only: rg_grid, xl_grid
        use KIM_kinds, only: dp
        use functions, only: varphi_l

        implicit none

        class(int_F1_rho_phi_t), intent(inout) :: this
        real(dp), intent(in) :: x, xp, theta
        real(dp) :: val
        real(dp) :: ks_val

        ks_val = 0.5d0 * (plasma%ks(this%int_point%j) + plasma%ks(this%int_point%j+1))
        this%int_point%a_coef = sqrt(1.0d0 / (1.0d0 + cos(theta))) / this%int_point%rhoT
        this%int_point%b_coef = 0.5d0 * (xp - x)
        call this%int_point%calc_Jrg1()

        val = varphi_l(xp, this%int_point%xlpm1, this%int_point%xlp, this%int_point%xlpp1) &
            * varphi_l(x, this%int_point%xlm1, this%int_point%xl, this%int_point%xlp1) &
            * 2.0d0 * pi / (this%int_point%rhoT**2.0d0 * sin(theta)) &
            * exp(- ks_val**2.0d0 * this%int_point%rhoT**2.0d0 &
                  - (x + xp)**2.0d0 / (4.0d0 * this%int_point%rhoT**2.0d0 * (1.0d0 - cos(theta)))) &
            * this%int_point%Jrg1

    end function


    function integrand_F2_rho_phi(this, x, xp, theta) result(val)

        use constants, only: pi
        use species, only: plasma
        use gsl_mod, only: erf => gsl_sf_erf
        use grid, only: rg_grid, xl_grid
        use KIM_kinds, only: dp
        use functions, only: varphi_l

        implicit none

        class(int_F2_rho_phi_t), intent(in) :: this
        real(dp), intent(in) :: x, xp, theta
        real(dp) :: val
        real(dp) :: ks_val

        ks_val = 0.5d0 * (plasma%ks(this%j) + plasma%ks(this%j+1))

        val = 0.0d0

    end function


    function integrand_F3_rho_phi(this, x, xp, theta) result(val)

        use constants, only: pi
        use species, only: plasma
        use gsl_mod, only: erf => gsl_sf_erf
        use grid, only: rg_grid, xl_grid
        use KIM_kinds, only: dp
        use functions, only: varphi_l

        implicit none

        class(int_F3_rho_phi_t), intent(in) :: this
        real(dp), intent(in) :: x, xp, theta
        real(dp) :: val
        real(dp) :: ks_val

        ks_val = 0.5d0 * (plasma%ks(this%j) + plasma%ks(this%j+1))

        val = 0.0d0

    end function

    subroutine calc_Jrg1(this)

        use constants, only: pi
        use grid, only: rg_grid

        implicit none

        class(integration_point_t), intent(inout) :: this

        this%Jrg1 = sqrt(pi) / (2.0d0 * this%a_coef) &
            *(erf(this%a_coef * (this%b_coef + rg_grid%xb(this%j+1))) &
            - erf(this%a_coef * (this%b_coef + rg_grid%xb(this%j))))

    end subroutine

    subroutine calc_Jrg2(this)

        use constants, only: pi
        use grid, only: rg_grid

        implicit none

        class(integration_point_t), intent(inout) :: this

        this%Jrg2 = - sqrt(pi) / (4.0d0 * this%a_coef**3.0d0) &
                    * ( &
                        (2.0d0 * this%a_coef**2.0d0 * (this%b_coef + this%xl)**2.0d0 + 1.0d0) &
                            * erf(this%a_coef * (this%b_coef + rg_grid%xb(this%j))) &
                        + sqrt(pi) * (2.0d0 * this%a_coef**2.0d0 * (this%b_coef + this%xl)**2.0d0 + 1.0d0) &
                            * erf(this%a_coef * (this%b_coef + rg_grid%xb(this%j+1))) &
                        - 2.0d0 * this%a_coef * exp(-this%a_coef**2.0d0 * (this%b_coef + rg_grid%xb(this%j))**2.0d0) &
                            * (this%b_coef - rg_grid%xb(this%j) + 2.0d0 * this%xl) &
                        + 2.0d0 * this%a_coef * exp(-this%a_coef**2.0d0 * (this%b_coef + rg_grid%xb(this%j+1))**2.0d0) &
                            * (this%b_coef - rg_grid%xb(this%j+1) + 2.0d0 * this%xl) &
                    )

    end subroutine


    subroutine calc_Jrg3(this)

        use constants, only: pi
        use grid, only: rg_grid

        implicit none

        class(integration_point_t), intent(inout) :: this

        this%Jrg3 = - sqrt(pi) / (4.0d0 * this%a_coef**3.0d0) &
                    * ( &
                        (2.0d0 * this%a_coef**2.0d0 * (this%b_coef - this%xlp)**2.0d0 + 1.0d0) &
                            * erf(this%a_coef * (this%b_coef + rg_grid%xb(this%j))) &
                        + sqrt(pi) * (2.0d0 * this%a_coef**2.0d0 * (this%b_coef - this%xlp)**2.0d0 + 1.0d0) &
                            * erf(this%a_coef * (this%b_coef + rg_grid%xb(this%j+1))) &
                        - 2.0d0 * this%a_coef * exp(-this%a_coef**2.0d0 * (this%b_coef + rg_grid%xb(this%j))**2.0d0) &
                            * (this%b_coef - rg_grid%xb(this%j) - 2.0d0 * this%xlp) &
                        + 2.0d0 * this%a_coef * exp(-this%a_coef**2.0d0 * (this%b_coef + rg_grid%xb(this%j+1))**2.0d0) &
                            * (this%b_coef - rg_grid%xb(this%j+1) - 2.0d0 * this%xlp) &
                    )

    end subroutine

    subroutine calc_Jrg4(this)

        use constants, only: pi
        use grid, only: rg_grid

        implicit none

        class(integration_point_t), intent(inout) :: this

        this%Jrg4 = - sqrt(pi) / (4.0d0 * this%a_coef**3.0d0) &
                    * ( &
                        erf(this%a_coef * (this%b_coef + rg_grid%xb(this%j))) &
                            * (-2.0d0 * this%a_coef**2.0d0 * (this%b_coef + this%xl)*(this%b_coef - this%xlp)-1.0d0) &
                        + sqrt(pi) * erf(this%a_coef * (this%b_coef + rg_grid%xb(this%j+1))) &
                            * (-2.0d0 * this%a_coef**2.0d0 * (this%b_coef + this%xl)*(this%b_coef - this%xlp)-1.0d0) &
                    )

    end subroutine

end module