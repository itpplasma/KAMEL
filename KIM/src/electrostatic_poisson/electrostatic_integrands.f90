module electrostatic_integrands

    use KIM_kinds, only: dp

    implicit none

    type :: int_B0_rho_phi_t
        real(dp) :: rhoT
        integer :: j
        real(dp) :: xlm1, xlp1, xl
        real(dp) :: xlpm1, xlpp1, xlp
        contains
            procedure :: f => integrand_mathcal_B0_rho_phi_semi_analytic
    end type

    type :: int_B1_rho_phi_t
        real(dp) :: rhoT
        integer :: j
        real(dp) :: xlm1, xlp1, xl
        real(dp) :: xlpm1, xlpp1, xlp
        contains
            procedure :: f => integrand_mathcal_B1_rho_phi_semi_analytic
    end type

    type :: int_struct_t
        real(dp) :: xlm1, xlp1, xl
        real(dp) :: xlpm1, xlpp1, xlp
        real(dp) :: rgj, rgjp1
        real(dp) :: rhoT, ks
        integer :: sp
    end type

    contains

    function integrand_mathcal_B0_rho_phi_semi_analytic(this, x) result(val)

        use grid, only: xl_grid, rg_grid
        use functions, only: varphi_l
        use gsl_mod, only: erf => gsl_sf_erf
        use constants, only: pi

        implicit none

        real(dp), intent(in) :: x
        real(dp) :: val
        class(int_B0_rho_phi_t), intent(in) :: this

        val = varphi_l(x, this%xlm1, this%xl, this%xlp1) &
            * varphi_l(x, this%xlpm1, this%xlp, this%xlpp1) &
            * (&
                erf((x - rg_grid%xb(this%j))/(sqrt(2.0d0) * this%rhoT)) &
                - erf((x - rg_grid%xb(this%j+1))/(sqrt(2.0d0) * this%rhoT))&
            )

    end function

    function integrand_mathcal_B1_rho_phi_semi_analytic(this, x, xp, theta) result(val)

        use constants, only: pi
        use species, only: plasma
        use gsl_mod, only: erf => gsl_sf_erf
        use grid, only: rg_grid, xl_grid
        use KIM_kinds, only: dp
        use functions, only: varphi_l

        implicit none

        class(int_B1_rho_phi_t), intent(in) :: this
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

    function G1_rho_phi(j, spec) result(val)

        use KIM_kinds, only: dp
        use species, only: species_t, plasma

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: ks_val, lambda, kpar, A1, A2, z0, rhoT
        complex(dp) :: plasma_Z

        ks_val = 0.5d0 * (plasma%ks(j) + plasma%ks(j+1))
        lambda = 0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1))
        kpar = 0.5d0 * (plasma%kp(j) + plasma%kp(j+1))
        A1 = 0.5d0 * (spec%A1(j) + spec%A1(j+1))
        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        z0 = 0.5d0 * (spec%z0(j) + spec%z0(j+1))
        rhoT = 0.5d0 * (spec%rho_L(j) + spec%rho_L(j+1))

        val = ks_val * rhoT /(lambda**2.0d0 * abs(kpar) * sqrt(2.0d0)) &
            * (&
                A1 * plasma_Z(z0) + A2 * plasma_Z(z0) * (1.0d0 + z0**2.0d0) + z0 * A2 &
            )

    end function

    function G2_rho_phi(j, spec) result(val)

        use KIM_kinds, only: dp
        use species, only: species_t, plasma

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: ks_val, lambda, kpar, A1, A2, z0, rhoT
        complex(dp) :: plasma_Z

        ks_val = 0.5d0 * (plasma%ks(j) + plasma%ks(j+1))
        lambda = 0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1))
        kpar = 0.5d0 * (plasma%kp(j) + plasma%kp(j+1))
        A1 = 0.5d0 * (spec%A1(j) + spec%A1(j+1))
        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        z0 = 0.5d0 * (spec%z0(j) + spec%z0(j+1))
        rhoT = 0.5d0 * (spec%rho_L(j) + spec%rho_L(j+1))

        val = ks_val /(lambda**2.0d0 * abs(kpar) * sqrt(2.0d0)) * A2

    end function

    function G3_rho_phi(j, spec) result(val)

        use KIM_kinds, only: dp
        use species, only: species_t, plasma

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: ks_val, lambda, kpar, A1, A2, z0, rhoT
        complex(dp) :: plasma_Z ! plasma dispersion function

        ks_val = 0.5d0 * (plasma%ks(j) + plasma%ks(j+1))
        lambda = 0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1))
        kpar = 0.5d0 * (plasma%kp(j) + plasma%kp(j+1))
        A1 = 0.5d0 * (spec%A1(j) + spec%A1(j+1))
        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        z0 = 0.5d0 * (spec%z0(j) + spec%z0(j+1))
        rhoT = 0.5d0 * (spec%rho_L(j) + spec%rho_L(j+1))

        val = ks_val /(lambda**2.0d0 * abs(kpar) * sqrt(2.0d0)) * A2

    end function

    function G1_rho_B(j, spec) result(val)

        use KIM_kinds, only: dp
        use species, only: species_t, plasma

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: ks_val, lambda, kpar, A1, A2, z0, rhoT
        complex(dp) :: plasma_Z

        A1 = 0.5d0 * (spec%A1(j) + spec%A1(j+1))
        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        z0 = 0.5d0 * (spec%z0(j) + spec%z0(j+1))

        val = 0.5d0 * A1 * (z0 * plasma_Z(z0) + 1.0d0) + A2 * &
            (&
                0.5d0 + (z0 * plasma_Z(z0) + 1.0d0) * (1.0d0 + z0**2.0d0) &
            )

    end function

    function G2_rho_B(j, spec) result(val)

        use KIM_kinds, only: dp
        use species, only: species_t, plasma
    
        implicit none

        integer, intent(in) :: j
        type(species_t),intent(in) :: spec
        complex(dp) :: val
        real(dp) :: ks_val, lambda, kpar, A1, A2, z0, rhoT
        complex(dp) :: plasma_Z

        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        z0 = 0.5d0 * (spec%z0(j) + spec%z0(j+1))

        val = A2 * (z0 * plasma_Z(z0) + 1.0d0)

    end function


    function G0_rho_phi(j, spec) result(val)

        use species, only: plasma, species_t
        use constants, only: pi

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val
        real(dp) :: lambda

        lambda = 0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1))
        val = 2.0d0 * pi**2.0d0 / (-lambda**2.0d0)  !/ sqrt(2.0d0)

    end function


    ! completely numerical integration

    function int_B0_rho_phi(x, xp, rg, theta, rhoT, ks) result(val)

        use constants, only: pi
        use KIM_kinds, only: dp

        implicit none

        real(dp), intent(in) :: x, xp, rg, theta, rhoT, ks
        real(dp) :: val

        if (abs(x-xp) <= 1e-15) then
            val = (sqrt(2.0d0 * pi))**(1.5d0) / rhoT * exp(- (x - rg)**2.0d0 / (2.0d0 * rhoT**2.0d0))
        else
            val = 0.0d0
        end if

    end function

    function int_B1_rho_phi(x, xp, rg, theta, rhoT, ks) result(val)

        use KIM_kinds, only: dp

        implicit none

        real(dp), intent(in) :: x, xp, rg, theta, rhoT, ks
        real(dp) :: val

        val = 2.0d0 / (rhoT**2.0d0 * sin(theta)) * zeta(x, xp, rg, theta, rhoT, ks)

    end function

    function int_B2_rho_phi(x, xp, rg, theta, rhoT, ks) result(val)

        use KIM_kinds, only: dp

        implicit none

        real(dp), intent(in) :: x, xp, rg, theta, rhoT, ks
        real(dp) :: val

        val = -zeta(x, xp, rg, theta, rhoT, ks) /(4.0d0 * rhoT**4.0d0 * sin(theta)**5.0d0) &
            * (&
                2.0d0 * cos(2.0d0*theta) &
                * (&
                    (x - rg)**2.0d0 + (xp + rg)**2.0d0 + 2.0d0 * rhoT**2.0d0 * (ks**2.0d0 * rhoT**2.0d0 + 1.0d0) &
                )&
                - ks**2.0d0 * rhoT**4.0d0 * cos(4.0d0*theta) &
                - rhoT**2.0d0 * (3.0d0 * ks**2.0d0 * rhoT**2.0d0 + 4.0d0) &
                + 6.0d0 * ((x - rg)**2.0d0 + (xp + rg)**2.0d0) &
                + 16.0d0 * cos(theta) * (x - rg) * (xp + rg) &
            )

    end function


    function zeta(x, xp, rg, theta, rhoT, ks) result(val)

        use KIM_kinds, only: dp

        implicit none
        real(dp), intent(in) :: x, xp, rg, theta, rhoT, ks
        real(dp) :: val

        val = exp(-ks**2.0d0 * rhoT**2.0d0 &
            - (rg + 0.5d0 * (xp -x))**2.0d0 * (1.0d0 - cos(theta)) / (rhoT**2.0d0 * sin(theta)**2.0d0) &
            - (xp +x)**2.0d0 * (1.0d0 + cos(theta)) / (4.0d0 * rhoT**2.0d0 * sin(theta)**2.0d0))

    end function

end module