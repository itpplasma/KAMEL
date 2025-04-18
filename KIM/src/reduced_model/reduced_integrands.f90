module reduced_integrands

    use KIM_kinds, only: dp

    implicit none

    type :: int_B0_t
        real(dp) :: rhoT
        integer :: l, lp, j
        contains
            procedure :: f => integrand_mathcal_B0
    end type

    type :: int_B1_t
        real(dp) :: rhoT
        integer :: l, lp, j
        contains
            procedure :: f => integrand_mathcal_B1
    end type

    contains

    function integrand_mathcal_B0(this, r) result(val)

        use grid, only: xl_grid, rg_grid
        use functions, only: varphi_l
        use gsl_mod, only: erf => gsl_sf_erf

        implicit none

        real(dp), intent(in) :: r
        real(dp) :: val
        class(int_B0_t), intent(in) :: this

        ! You can access l, lp, j from the host scope here
        val = varphi_l(r, xl_grid%xb(this%l-1), xl_grid%xb(this%l), xl_grid%xb(this%l+1)) &
            * varphi_l(r, xl_grid%xb(this%lp-1), xl_grid%xb(this%lp), xl_grid%xb(this%lp+1)) &
            * (&
                erf((r - rg_grid%xb(this%j))/(sqrt(2.0d0) * this%rhoT)) - erf((r - rg_grid%xb(this%j+1))/(sqrt(2.0d0) * this%rhoT))&
            )

    end function integrand_mathcal_B0

    function integrand_mathcal_B1(this, x, xp, theta) result(val)

        use constants, only: pi
        use species, only: plasma
        use gsl_mod, only: erf => gsl_sf_erf
        use grid, only: rg_grid, xl_grid
        use KIM_kinds, only: dp
        use functions, only: varphi_l
        use, intrinsic :: ieee_arithmetic

        implicit none

        class(int_B1_t), intent(in) :: this
        real(dp), intent(in) :: x, xp, theta
        real(dp) :: val
        real(dp) :: ks_val

        ks_val = 0.5d0 * (plasma%ks(this%j) + plasma%ks(this%j+1))

        val = - sqrt(pi) / (this%rhoT * sqrt(1.0d0 - cos(theta))) * exp(- ks_val**2.0d0 * this%rhoT**2.0d0) &
            * varphi_l(xp, xl_grid%xb(this%lp-1), xl_grid%xb(this%lp), xl_grid%xb(this%lp+1)) &
            * varphi_l(x, xl_grid%xb(this%l-1), xl_grid%xb(this%l), xl_grid%xb(this%l+1)) &
            * exp(- (x+xp)**2.0d0 / (4.0d0 * this%rhoT**2.0d0 * (1.0d0 - cos(theta)))) &
            * (&
                erf((0.5d0 * (xp - x) + rg_grid%xb(this%j+1))/(this%rhoT * sin(theta)**2.0d0) * sqrt(1.0d0 - cos(theta))) &
                -erf((0.5d0 * (xp - x) + rg_grid%xb(this%j))/(this%rhoT * sin(theta)**2.0d0) * sqrt(1.0d0 - cos(theta))) &
            ) 

        if (isnan(val)) then
            print *, "integrand_mathcal_B1: NaN for j = ", this%j, " l = ", this%l, " lp = ", this%lp
            print *, "x = ", x
            print *, "xp = ", xp
            print *, "theta = ", theta
            print *, "rhoT = ", this%rhoT
            print *, "ks_val = ", ks_val
            print *, "cos(theta) = ", cos(theta)
            print *, "sin(theta) = ", sin(theta)
            print *, "exp(-ks_val**2.0d0 * this%rhoT**2.0d0) = ", exp(- ks_val**2.0d0 * this%rhoT**2.0d0)
            print *, "exp(- (x+xp)**2.0d0 / (4.0d0 * this%rhoT**2.0d0 * (cos(theta) - 1.0d0))) = ", &
                exp(- (x+xp)**2.0d0 / (4.0d0 * this%rhoT**2.0d0 * (cos(theta) - 1.0d0)))
            print *, "exp((x+xp)**2.0d0 * (cos(theta)+1.0d0)/ (4.0d0 * this%rhoT**2.0d0 * sin(theta)**2.0d0)) = ", &
                exp((x+xp)**2.0d0 * (cos(theta)+1.0d0)/ (4.0d0 * this%rhoT**2.0d0 * sin(theta)**2.0d0))
            stop
        end if
        if (val >= 1.0d-3) then
            print *, "val = ", val
        end if

    end function

    function mathcal_A0(j, spec) result(val)

        use grid, only: xl_grid, rg_grid
        use gsl_mod, only: erf => gsl_sf_erf
        use back_quants, only: lambda_De, lambda_Di
        use species, only: plasma, species_t
        use constants, only: pi

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val
        real(dp) :: lambda

        lambda = 0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1))
        val = 8.0d0 * pi**3.0d0 / (-lambda**2.0d0)  / (8.0d0 * pi)

    end function mathcal_A0

    function mathcal_A1(j, spec) result(val)

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

        val = ks_val * rhoT /(lambda**2.0d0 * abs(kpar) * sqrt(2.0d0)) &
            * (&
                A1 * plasma_Z(z0) + A2 * plasma_Z(z0) * (1.0d0 + z0**2.0d0) + z0 * A2 &
            )

    end function

end module