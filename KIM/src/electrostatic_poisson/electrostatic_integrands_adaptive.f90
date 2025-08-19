! integrands for both Krook and FP collision models
module electrostatic_integrands_rkf45_mod

    use KIM_kinds, only: dp
    use electrostatic_integrands_gauss_mod, only: integration_point_t

    implicit none

    type :: rkf45_int_F1_rho_phi_t
        type(integration_point_t) :: int_point
        contains
            procedure :: f => rkf45_integrand_F1_rho_phi
    end type

    type :: rkf45_int_F2_rho_phi_t
        type(integration_point_t) :: int_point
        contains
            procedure :: f => rkf45_integrand_F2_rho_phi
    end type

    type :: rkf45_int_F3_rho_phi_t
        type(integration_point_t) :: int_point
        contains
            procedure :: f => rkf45_integrand_F3_rho_phi
    end type

    contains

    function rkf45_integrand_F1_rho_phi(this, x, xp, theta) result(val)

        use constants, only: pi
        use species, only: plasma
        use gsl_mod, only: erf => gsl_sf_erf
        use KIM_kinds, only: dp
        use functions, only: varphi_l

        implicit none

        class(rkf45_int_F1_rho_phi_t), intent(inout) :: this
        real(dp), intent(in) :: x, xp, theta
        real(dp) :: val
        real(dp) :: ks_val

        ks_val = 0.5d0 * (plasma%ks(this%int_point%j) + plasma%ks(this%int_point%j+1))
        this%int_point%a_coef = sqrt(1.0d0 / (1.0d0 + cos(theta))) / abs(this%int_point%rhoT)
        this%int_point%b_coef = calc_b_coef(x, xp)

        call this%int_point%calc_Jrg1()

        val = varphi_l(xp, this%int_point%xlpm1, this%int_point%xlp, this%int_point%xlpp1) &
            * varphi_l(x, this%int_point%xlm1, this%int_point%xl, this%int_point%xlp1) &
            * 2.0d0 * pi / (this%int_point%rhoT**2.0d0 * sin(theta)) &
            * exp(- ks_val**2.0d0 * this%int_point%rhoT**2.0d0 &
                  - (x - xp)**2.0d0 / (4.0d0 * this%int_point%rhoT**2.0d0 * (1.0d0 - cos(theta)))) &
            * this%int_point%Jrg1

    end function


    function rkf45_integrand_F2_rho_phi(this, x, xp, theta) result(val)

        use constants, only: pi
        use species, only: plasma
        use gsl_mod, only: erf => gsl_sf_erf
        use KIM_kinds, only: dp
        use functions, only: varphi_l

        implicit none

        class(rkf45_int_F2_rho_phi_t), intent(inout) :: this
        real(dp), intent(in) :: x, xp, theta
        real(dp) :: val
        real(dp) :: ks_val

        ks_val = 0.5d0 * (plasma%ks(this%int_point%j) + plasma%ks(this%int_point%j+1))
        this%int_point%a_coef = sqrt(1.0d0 / (1.0d0 + cos(theta))) / abs(this%int_point%rhoT)
        this%int_point%b_coef = calc_b_coef(x, xp)

        call this%int_point%calc_Jrg1()
        call this%int_point%calc_Jrg2()
        call this%int_point%calc_Jrg3()
        call this%int_point%calc_Jrg4()

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


    function rkf45_integrand_F3_rho_phi(this, x, xp, theta) result(val)

        use constants, only: pi
        use species, only: plasma
        use gsl_mod, only: erf => gsl_sf_erf
        use KIM_kinds, only: dp
        use functions, only: varphi_l

        implicit none

        class(rkf45_int_F3_rho_phi_t), intent(inout) :: this
        real(dp), intent(in) :: x, xp, theta
        real(dp) :: val
        real(dp) :: ks_val

        ks_val = 0.5d0 * (plasma%ks(this%int_point%j) + plasma%ks(this%int_point%j+1))
        this%int_point%a_coef = sqrt(1.0d0 / (1.0d0 + cos(theta))) / abs(this%int_point%rhoT)
        this%int_point%b_coef = calc_b_coef(x, xp)

        call this%int_point%calc_Jrg1()
        call this%int_point%calc_Jrg2()
        call this%int_point%calc_Jrg3()
        call this%int_point%calc_Jrg4()

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


    function calc_b_coef(x,xp) result(b_coef)

        implicit none

        real(dp), intent(in) :: x, xp
        real(dp) :: b_coef

        b_coef = 0.5d0 * (xp + x)

    end function

end module electrostatic_integrands_rkf45_mod