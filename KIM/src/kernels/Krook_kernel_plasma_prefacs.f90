module Krook_kernel_plasma_prefacs_m

    contains

    function Krook_kappa_rho_phi(j, spec) result(val)

        use species_m, only: species_t
        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val
        real(dp) :: lambda

        lambda = 0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1))
        val = 1.0d0 / (lambda**2.0d0)  !/ sqrt(2.0d0)

    end function

    function Krook_kappa_rho_B(j, spec) result(val)

        use species_m, only: plasma, species_t
        use KIM_kinds_m, only: dp
        use constants_m, only: sol, com_unit

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val

        val = (0.5d0 * (spec%vT(j) + spec%vT(j+1)))**2.0d0 * com_unit &
            / (0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1)))**2.0d0 &
            / (0.5d0 * (spec%omega_c(j) + spec%omega_c(j+1))) &
            / (0.5d0 * (plasma%kp(j) + plasma%kp(j+1))) &
            / sol

    end function


    function Krook_G0_rho_phi(j, spec) result(val)

        use species_m, only: species_t
        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val

        ! Unused arguments j and spec - required for interface consistency
        val = -1.0d0 !/ sqrt(2.0d0)

    end function

    function Krook_G1_rho_phi(j, spec) result(val)

        use KIM_kinds_m, only: dp
        use species_m, only: species_t, plasma

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: ks_val, kpar, A1, A2, rhoT
        complex(dp) :: plasma_Z, z0

        ks_val = 0.5d0 * (plasma%ks(j) + plasma%ks(j+1))
        kpar = 0.5d0 * (plasma%kp(j) + plasma%kp(j+1))
        A1 = 0.5d0 * (spec%A1(j) + spec%A1(j+1))
        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        z0 = 0.5d0 * (spec%z0(j) + spec%z0(j+1))
        rhoT = 0.5d0 * (spec%rho_L(j) + spec%rho_L(j+1))

        val = ks_val * rhoT /(abs(kpar) * sqrt(2.0d0)) &
            * (&
                A1 * plasma_Z(z0) + A2 * plasma_Z(z0) * (1.0d0 + z0**2.0d0) + z0 * A2 &
            )

    end function

    function Krook_G2_rho_phi(j, spec) result(val)

        use KIM_kinds_m, only: dp
        use species_m, only: species_t, plasma

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: ks_val, kpar, A2, rhoT
        complex(dp) :: plasma_Z, z0

        ks_val = 0.5d0 * (plasma%ks(j) + plasma%ks(j+1))
        kpar = 0.5d0 * (plasma%kp(j) + plasma%kp(j+1))
        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        rhoT = 0.5d0 * (spec%rho_L(j) + spec%rho_L(j+1))
        z0 = 0.5d0 * (spec%z0(j) + spec%z0(j+1))

        val = - ks_val * rhoT * A2 * plasma_Z(z0) / (abs(kpar) * sqrt(2.0d0))

    end function

    function Krook_G3_rho_phi(j, spec) result(val)

        use KIM_kinds_m, only: dp
        use species_m, only: species_t, plasma

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: ks_val, kpar, A2, rhoT

        ks_val = 0.5d0 * (plasma%ks(j) + plasma%ks(j+1))
        kpar = 0.5d0 * (plasma%kp(j) + plasma%kp(j+1))
        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        rhoT = 0.5d0 * (spec%rho_L(j) + spec%rho_L(j+1))

        val = ks_val * rhoT /(abs(kpar) * sqrt(2.0d0)) * A2

    end function

    function Krook_G0_rho_B(j, spec) result(val)

        use species_m, only: species_t
        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val

        ! Unused arguments j and spec - required for interface consistency
        val = 0.0d0

    end function

    function Krook_G1_rho_B(j, spec) result(val)

        use KIM_kinds_m, only: dp
        use species_m, only: species_t

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: A1, A2
        complex(dp) :: plasma_Z, z0

        A1 = 0.5d0 * (spec%A1(j) + spec%A1(j+1))
        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        z0 = 0.5d0 * (spec%z0(j) + spec%z0(j+1))

        val = 0.5d0 * A1 * (z0 * plasma_Z(z0) + 1.0d0) + A2 * &
            (&
                0.5d0 + (z0 * plasma_Z(z0) + 1.0d0) * (1.0d0 + z0**2.0d0) &
            )

    end function

    function Krook_G2_rho_B(j, spec) result(val)

        use KIM_kinds_m, only: dp
        use species_m, only: species_t

        implicit none

        integer, intent(in) :: j
        type(species_t),intent(in) :: spec
        complex(dp) :: val
        real(dp) :: A2
        complex(dp) :: plasma_Z, z0

        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        z0 = 0.5d0 * (spec%z0(j) + spec%z0(j+1))

        val = - A2 * (z0 * plasma_Z(z0) + 1.0d0)

    end function

    function Krook_G3_rho_B(j, spec) result(val)

        use KIM_kinds_m, only: dp
        use species_m, only: species_t

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val

        val = 0.5d0 * (spec%A2(j) + spec%A2(j+1))

    end function

end module
