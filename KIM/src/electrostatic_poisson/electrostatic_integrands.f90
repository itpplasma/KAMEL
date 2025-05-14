module electrostatic_integrands

    use KIM_kinds, only: dp

    implicit none

    type :: int_F0_rho_phi_t
        real(dp) :: rhoT
        integer :: j
        real(dp) :: xlm1, xlp1, xl
        real(dp) :: xlpm1, xlpp1, xlp
        contains
            procedure :: f => integrand_F0_rho_phi
    end type

    type :: int_F1_rho_phi_t
        real(dp) :: rhoT
        integer :: j
        real(dp) :: xlm1, xlp1, xl
        real(dp) :: xlpm1, xlpp1, xlp
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

        use grid, only: xl_grid, rg_grid
        use functions, only: varphi_l
        use gsl_mod, only: erf => gsl_sf_erf
        use constants, only: pi

        implicit none

        real(dp), intent(in) :: x
        real(dp) :: val
        class(int_F0_rho_phi_t), intent(in) :: this

        val = varphi_l(x, this%xlm1, this%xl, this%xlp1) &
            * varphi_l(x, this%xlpm1, this%xlp, this%xlpp1) &
            * (&
                erf((x - rg_grid%xb(this%j))/(sqrt(2.0d0) * this%rhoT)) &
                - erf((x - rg_grid%xb(this%j+1))/(sqrt(2.0d0) * this%rhoT))&
            )

    end function

    function integrand_F1_rho_phi(this, x, xp, theta) result(val)

        use constants, only: pi
        use species, only: plasma
        use gsl_mod, only: erf => gsl_sf_erf
        use grid, only: rg_grid, xl_grid
        use KIM_kinds, only: dp
        use functions, only: varphi_l

        implicit none

        class(int_F1_rho_phi_t), intent(in) :: this
        real(dp), intent(in) :: x, xp, theta
        real(dp) :: val
        real(dp) :: ks_val

        ks_val = 0.5d0 * (plasma%ks(this%j) + plasma%ks(this%j+1))

        val = - sqrt(pi) / (sqrt(1.0d0 - cos(theta)) * this%rhoT) * exp(- ks_val**2.0d0 * this%rhoT**2.0d0) &
            * varphi_l(xp, this%xlpm1, this%xlp, this%xlpp1) &
            * varphi_l(x, this%xlm1, this%xl, this%xlp1) &
            * exp(- (x+xp)**2.0d0 / (4.0d0 * this%rhoT**2.0d0 * (1.0d0 - cos(theta)))) &
            * (&
                erf(0.5d0 * (xp - x) + rg_grid%xb(this%j+1))/(this%rhoT * sqrt(1.0d0 + cos(theta))) &
                -erf(0.5d0 * (xp - x) + rg_grid%xb(this%j))/(this%rhoT * sqrt(1.0d0 + cos(theta))) &
            ) 

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

end module