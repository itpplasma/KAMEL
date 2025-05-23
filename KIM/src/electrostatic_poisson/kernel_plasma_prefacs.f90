module kernel_plasma_prefacs


    contains


    function kappa_rho_phi(j, spec) result(val)

        use species, only: plasma, species_t
        use constants, only: pi
        use KIM_kinds, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val
        real(dp) :: lambda

        lambda = 0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1))
        val = 1.0d0 / (lambda**2.0d0)  !/ sqrt(2.0d0)

    end function

    function kappa_rho_B(j, spec) result(val)

        use species, only: plasma, species_t
        use KIM_kinds, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val

        val = (0.5d0 * (spec%vT(j) + spec%vT(j+1)))**2.0d0 &
            / (0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1)))**2.0d0 &
            / (0.5d0 * (spec%omega_c(j) + spec%omega_c(j+1))) &
            / abs(0.5d0 * (plasma%kp(j) + plasma%kp(j+1)))

    end function


    function G0_rho_phi(j, spec) result(val)

        use species, only: plasma, species_t
        use constants, only: pi
        use KIM_kinds, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val

        val = -1.0d0 

    end function

    function G1_rho_phi(j, spec) result(val)

        use KIM_kinds, only: dp
        use species, only: species_t, plasma

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: ks_val, kpar, A1, A2, z0, rhoT
        complex(dp) :: plasma_Z

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

    function G2_rho_phi(j, spec) result(val)

        use KIM_kinds, only: dp
        use species, only: species_t, plasma

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: ks_val, kpar, A1, A2, rhoT
        complex(dp) :: plasma_Z

        ks_val = 0.5d0 * (plasma%ks(j) + plasma%ks(j+1))
        kpar = 0.5d0 * (plasma%kp(j) + plasma%kp(j+1))
        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        rhoT = 0.5d0 * (spec%rho_L(j) + spec%rho_L(j+1))

        val = ks_val * rhoT /(abs(kpar) * sqrt(2.0d0)) * A2

    end function

    function G3_rho_phi(j, spec) result(val)

        use KIM_kinds, only: dp
        use species, only: species_t, plasma

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val

        val = G2_rho_phi(j, spec)

    end function

    function G1_rho_B(j, spec) result(val)

        use KIM_kinds, only: dp
        use species, only: species_t, plasma

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: A1, A2, z0
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
        real(dp) :: A2, z0
        complex(dp) :: plasma_Z

        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        z0 = 0.5d0 * (spec%z0(j) + spec%z0(j+1))

        val = A2 * (z0 * plasma_Z(z0) + 1.0d0)

    end function

    function G3_rho_B(j, spec) result(val)

        use KIM_kinds, only: dp
        use species, only: species_t, plasma

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val

        val = 0.5d0 * (spec%A2(j) + spec%A2(j+1))

    end function

end module