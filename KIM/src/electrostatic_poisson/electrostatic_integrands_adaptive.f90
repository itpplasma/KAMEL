! integrands for both Krook and FP collision models
module electrostatic_integrands_rkf45_mod

    use KIM_kinds_m, only: dp

    implicit none

    type :: rkf45_integrand_context_t
        integer :: j
        real(dp) :: rhoT
        real(dp) :: x, xp, xlpm1, xlp, xlpp1, xlm1, xlp1, xl
        real(dp) :: ks
    end type rkf45_integrand_context_t

    contains

    function rkf45_integrand_F0(context) result(val)

        use grid_m, only: rg_grid
        use functions_m, only: varphi_l
        use gsl_mod, only: erf => gsl_sf_erf
        use constants_m, only: pi

        implicit none

        type(rkf45_integrand_context_t), intent(in) :: context
        real(dp) :: val

        ! delta(x-x') results in same argument in both varphi_l functions
        val = varphi_l(context%x, context%xlm1, context%xl, context%xlp1) &
            * varphi_l(context%x, context%xlpm1, context%xlp, context%xlpp1) &
            * (&
                  erf((rg_grid%xb(context%j+1)-context%x)/(sqrt(2.0d0) * abs(context%rhoT))) &
                - erf((rg_grid%xb(context%j) - context%x)/(sqrt(2.0d0) * abs(context%rhoT)))&
            ) &
            * 2.0d0 * pi**2.0d0

    end function

    function rkf45_integrand_F1(theta, context) result(val)

        use constants_m, only: pi
        use species_m, only: plasma
        use gsl_mod, only: erf => gsl_sf_erf
        use KIM_kinds_m, only: dp
        use functions_m, only: varphi_l

        implicit none

        !type(rkf45_integrand_context_t), intent(in) :: context
        class(*), intent(in) :: context
        real(dp), intent(in) :: theta
        real(dp) :: val
        real(dp) :: a, b
        real(dp) :: sin_t, cos_t, denom_1mcos, denom_1pcos

        select type(context)
        type is(rkf45_integrand_context_t)
    
            ! Numerically stable half-angle forms
            sin_t = sin(theta)
            cos_t = cos(theta)
            denom_1pcos = 2.0d0 * max(cos(0.5d0*theta)**2.0d0, 1.0d-300)
            a = sqrt(1.0d0 / denom_1pcos) / max(abs(context%rhoT), 1.0d-300)
            b = calc_b_coef(context%x, context%xp)

            val = varphi_l(context%xp, context%xlpm1, context%xlp, context%xlpp1) &
                * varphi_l(context%x, context%xlm1, context%xl, context%xlp1) &
                * 2.0d0 * pi / (context%rhoT**2.0d0 * max(sin_t, 1.0d-300)) &
                * exp(- (context%x - context%xp)**2.0d0 / (4.0d0 * max(context%rhoT**2.0d0, 1.0d-300) * &
                           max(2.0d0 * sin(0.5d0*theta)**2.0d0, 1.0d-300))) &
                * Jrg1(a, b, context%j)

        class default
            error stop "Unsupported context type"
        end select


    end function


    function rkf45_integrand_F2(theta, context) result(val)

        use constants_m, only: pi
        use species_m, only: plasma
        use KIM_kinds_m, only: dp
        use functions_m, only: varphi_l

        implicit none

        !type(rkf45_integrand_context_t), intent(in) :: context
        class(*), intent(in) :: context
        real(dp), intent(in) :: theta
        real(dp) :: val
        real(dp) :: a, b
        real(dp) :: sin_t, cos_t

        select type(context)
        type is(rkf45_integrand_context_t)

            sin_t = sin(theta)
            cos_t = cos(theta)
            a = 1.0d0 / (max(abs(context%rhoT), 1.0d-300) * sqrt(max(2.0d0 * cos(0.5d0*theta)**2.0d0, 1.0d-300)))
            b = calc_b_coef(context%x, context%xp)

            val = varphi_l(context%xp, context%xlpm1, context%xlp, context%xlpp1) &
                * varphi_l(context%x, context%xlm1, context%xl, context%xlp1) &
                * 1.0d0 / max(sin_t, 1.0d-300)**5.0d0 &
                * exp(- (context%x - context%xp)**2.0d0 / (4.0d0 * max(context%rhoT**2.0d0, 1.0d-300) * &
                           max(2.0d0 * sin(0.5d0*theta)**2.0d0, 1.0d-300))) & 
                * ( &
                    Jrg1(a, b, context%j) * (&
                        4.0d0 * cos(2.0d0 * theta) * context%rhoT**2.0d0 *(context%ks**2.0d0 * context%rhoT**2.0d0 + 1.0d0) &
                        - context%ks**2.0d0 * context%rhoT**4.0d0 * cos(4.0d0 * theta)  &
                        - context%rhoT**2.0d0 * (3.0d0 * context%ks**2.0d0 * context%rhoT**2.0d0 + 4.0d0) &
                    ) &
                    + (2.0d0 * cos(2.0d0 * theta) + 6.0d0) &
                    * (Jrg2(a, b, context%j, context%x) + Jrg3(a, b, context%j, context%xp)) &
                    - 16.0d0 * cos(theta) * Jrg4(a, b, context%j, context%x, context%xp) &
                )

        class default
            error stop "Unsupported context type"
        end select

    end function


    function rkf45_integrand_F3(theta, context) result(val)

        use constants_m, only: pi
        use species_m, only: plasma
        use KIM_kinds_m, only: dp
        use functions_m, only: varphi_l

        implicit none

        !type(rkf45_integrand_context_t), intent(in) :: context
        class(*), intent(in) :: context
        real(dp), intent(in) :: theta
        real(dp) :: val
        real(dp) :: a, b
        real(dp) :: sin_t, cos_t

        select type(context)
        type is(rkf45_integrand_context_t)

            sin_t = sin(theta)
            cos_t = cos(theta)
            a = 1.0d0 / (max(abs(context%rhoT), 1.0d-300) * sqrt(max(2.0d0 * cos(0.5d0*theta)**2.0d0, 1.0d-300)))
            b = calc_b_coef(context%x, context%xp)

            val = varphi_l(context%xp, context%xlpm1, context%xlp, context%xlpp1) &
                * varphi_l(context%x, context%xlm1, context%xl, context%xlp1) &
                * cos_t / max(sin_t, 1.0d-300)**5.0d0 &
                * exp(- (context%x - context%xp)**2.0d0 / (4.0d0 * max(context%rhoT**2.0d0, 1.0d-300) * &
                           max(2.0d0 * sin(0.5d0*theta)**2.0d0, 1.0d-300))) & 
                * ( &
                    context%rhoT**2.0d0 * (cos(3.0d0 * theta) - cos_t) * Jrg1(a, b, context%j) &
                    + 4.0d0 * cos_t * (Jrg2(a, b, context%j, context%x) + Jrg3(a, b, context%j, context%xp)) &
                    - 2.0d0 * (cos(2.0d0 * theta) + 3.0d0) * Jrg4(a, b, context%j, context%x, context%xp) &
                )

        class default
            error stop "Unsupported context type"
        end select

    end function


    function calc_b_coef(x,xp) result(b_coef)

        implicit none

        real(dp), intent(in) :: x, xp
        real(dp) :: b_coef

        b_coef = 0.5d0 * (xp + x)

    end function

    function Jrg1(a, b, j)

        use constants_m, only: pi
        use grid_m, only: rg_grid
        use numerics_utils_m, only: erf_diff

        implicit none

        real(dp), intent(in) :: a, b
        integer, intent(in) :: j
        real(dp) :: Jrg1
        Jrg1 = sqrt(pi) / (2.0d0 * a) &
            * erf_diff(a * (b - rg_grid%xb(j)), &
                       a * (b - rg_grid%xb(j+1)))

    end function

    function Jrg2(a, b, j, xl)

        use constants_m, only: pi
        use grid_m, only: rg_grid
        use numerics_utils_m, only: erf_diff

        implicit none

        real(dp), intent(in) :: a, b, xl
        integer, intent(in) :: j
        real(dp) :: Jrg2

        Jrg2 = 1.0d0 / (4.0d0 * a**3.0d0) &
                    * ( &
                        sqrt(pi) * (2.0d0 * a**2.0d0 * (b - xl)**2.0d0 + 1.0d0) &
                            * erf_diff(a * (b - rg_grid%xb(j)), &
                                       a * (b - rg_grid%xb(j+1))) &
                        + 2.0d0 * a * exp(-a**2.0d0 * (b - rg_grid%xb(j))**2.0d0) &
                            * (b + rg_grid%xb(j) - 2.0d0 * xl) &
                        - 2.0d0 * a * exp(-a**2.0d0 * (b - rg_grid%xb(j+1))**2.0d0) &
                            * (b + rg_grid%xb(j+1) - 2.0d0 * xl) &
                    )

    end function


    function Jrg3(a, b, j, xlp)

        use constants_m, only: pi
        use grid_m, only: rg_grid
        use numerics_utils_m, only: erf_diff

        implicit none

        real(dp), intent(in) :: a, b, xlp
        integer, intent(in) :: j
        real(dp) :: Jrg3

        jrg3 = 1.0d0 / (4.0d0 * a**3.0d0) &
                    * ( &
                        sqrt(pi) * (2.0d0 * a**2.0d0 * (b - xlp)**2.0d0 + 1.0d0) &
                            * erf_diff(a * (b - rg_grid%xb(j)), &
                                       a * (b - rg_grid%xb(j+1))) &
                        + 2.0d0 * a * exp(-a**2.0d0 * (b - rg_grid%xb(j))**2.0d0) &
                            * (b + rg_grid%xb(j) - 2.0d0 * xlp) &
                        - 2.0d0 * a * exp(-a**2.0d0 * (b - rg_grid%xb(j+1))**2.0d0) &
                            * (b + rg_grid%xb(j+1) - 2.0d0 * xlp) & 
                    )

    end function

    function Jrg4(a, b, j, xl, xlp)

        use constants_m, only: pi
        use grid_m, only: rg_grid
        use numerics_utils_m, only: erf_diff

        implicit none

        real(dp), intent(in) :: a, b, xl, xlp
        integer, intent(in) :: j
        real(dp) :: Jrg4

        Jrg4 = 1.0d0 / (4.0d0 * a**3.0d0) &
                    * ( &
                        erf_diff(a * (b - rg_grid%xb(j)), &
                                 a * (b - rg_grid%xb(j+1))) &
                            * sqrt(pi) * (2.0d0 * a**2.0d0 * (b - xl) &
                            * (b - xlp)+1.0d0) &
                        + 2.0d0 * a * exp(-a**2.0d0 * (b - rg_grid%xb(j))**2.0d0) &
                            * (b + rg_grid%xb(j) - xl - xlp) &
                        + 2.0d0 * a * exp(-a**2.0d0 * (b - rg_grid%xb(j+1))**2.0d0) &
                            * (-b - rg_grid%xb(j+1) + xl + xlp) &
                    )

    end function


end module electrostatic_integrands_rkf45_mod
